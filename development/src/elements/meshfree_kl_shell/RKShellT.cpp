/* RKShellT.cpp — meshfree RKPM Kirchhoff-Love shell element (epic #59)
 *
 * Tahoe-native element on ElementBaseT (see RKShellT.h for the architecture rationale).
 * Milestone A (this file): reparent + surface-meshfree setup -- resolve the field, read the
 * node cloud, build per-node neighbor lists, register the ragged neighbor connectivity /
 * equations with the framework, build the plane-stress tangent and the RKPM solver.
 * Milestone B (next): RHSDriver/LHSDriver assemble the validated KL-shell kernels per stencil.
 */
#include "RKShellT.h"

#include "ParameterT.h"
#include "LimitT.h"
#include "ExceptionT.h"
#include "ElementSupportT.h"
#include "ModelManagerT.h"
#include "FieldT.h"
#include "MLSSolverT.h"
#include "MeshFreeT.h"
#include "KLShellKernels.h"
#include "FlanaganTaylor.h"
#include "PlaneStressJ2.h"
#include "ElementMatrixT.h"
#include "eIntegratorT.h"
#include "OutputSetT.h"
#include "GeometryT.h"
#include "iArray2DT.h"

#include <cmath>
#include <vector>
#include <cstdlib>

using namespace Tahoe;
using namespace Tahoe::KLShell;

/* constructor */
RKShellT::RKShellT(const ElementSupportT& support):
	ElementBaseT(support),
	fThickness(0.0),
	fYoung(0.0),
	fPoisson(0.0),
	fSupportFac(3.0),
	fCompleteness(3),
	fNumNodes(0),
	fMLS(NULL)
{
	SetName("meshfree_kl_shell");
	fLoad[0] = fLoad[1] = fLoad[2] = 0.0;
	fDensity = 1.0;
	fYield = 0.0;
	fHardening = 0.0;
	fYieldSat = 0.0;
	fSatRate = 0.0;
	fFiniteStrain = 0;
	fThicknessUpdate = 0;
	fCorotational = 0;
	fMonitorNode = 0;
	fMonitorStride = 1;
	fMonitorCount = 1;
	fDamping = 0.0;
	fStabMode = 0;
	fStabMembrane = 1.0;
	fStabBending = 0.0;
	fStabNatural = 1.0;
	fOutputID = -1;
}

/* destructor */
RKShellT::~RKShellT(void) { delete fMLS; }

/* required ElementBaseT interface (minimal for a linear static shell) */
GlobalT::SystemTypeT RKShellT::TangentType(void) const { return GlobalT::kSymmetric; }
void RKShellT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
#pragma unused(field)
#pragma unused(node)
#pragma unused(force)
}
double RKShellT::InternalEnergy(void) { return 0.0; }
void RKShellT::SendOutput(int kincode) {
#pragma unused(kincode)
}

/* register the ragged neighbor connectivity / equations with the framework */
void RKShellT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
#pragma unused(eq_1)
	fEqnos.Configure(fNeighbors, NumDOF());
	Field().SetLocalEqnos(fNeighbors, fEqnos);
	eq_2.Append(&fEqnos);

	/* alias each element card's equations to its stencil eqnos (as MeshFreeElementSupportT does)
	 * so the framework -- including the explicit nodal integrator -- can resolve this element's
	 * active DOFs via CurrentElement().Equations() */
	for (int i = 0; i < fNumNodes; i++)
		fEqnos.RowAlias(i, fElementCards[i].Equations());
}

/* Register the displacement field for output on the background cell mesh (from the .geom), so
 * the deformed shell is viewable in ParaView (ExodusII / EnSight, format chosen in the deck).
 * Falls back to point output if no background cells exist. */
void RKShellT::RegisterOutput(void)
{
	ModelManagerT& model = ElementSupport().ModelManager();
	const ArrayT<StringT>& ids = model.ElementGroupIDs();

	/* add the plastic-strain field only when plasticity is active (keeps the elastic benchmark
	 * output unchanged so the Scordelis-Lo regression stays bit-exact) */
	int nf = (fYield > 0.0) ? ((fThicknessUpdate) ? 5 : 4) : 3;
	ArrayT<StringT> n_labels(nf);
	n_labels[0] = "D_X"; n_labels[1] = "D_Y"; n_labels[2] = "D_Z";
	if (nf >= 4) n_labels[3] = "EQ_PLASTIC_STRAIN";   /* max over thru-thickness J2 stations (surface) */
	if (nf >= 5) n_labels[4] = "THICKNESS";           /* current thickness (Algorithm 3 D33 accumulation) */

	if (ids.Length() > 0) {
		/* output on the background cells (a real surface mesh in ParaView) */
		ArrayT<StringT> block_ID(ids.Length());
		fOutputConn.Dimension(ids.Length());
		for (int b = 0; b < ids.Length(); b++) {
			block_ID[b] = ids[b];
			model.ReadConnectivity(ids[b]);            /* lazy-loaded: read before use */
			fOutputConn[b] = model.ElementGroupPointer(ids[b]);
		}
		/* a surface cell is a 2D manifold in 3D: pick the geometry from the node count
		 * (the model mis-guesses a 4-node-in-3D cell as a tetrahedron) */
		int nen = fOutputConn[0]->MinorDim();
		GeometryT::CodeT geo = (nen == 3) ? GeometryT::kTriangle
		                     : (nen == 4) ? GeometryT::kQuadrilateral
		                     :              model.ElementGroupGeometry(ids[0]);
		ArrayT<StringT> e_labels; /* none */
		OutputSetT output_set(geo, block_ID, fOutputConn, n_labels, e_labels, false);
		fOutputID = ElementSupport().RegisterOutput(output_set);
		fOutputNodesUsed = output_set.NodesUsed();
	} else {
		/* point cloud output */
		iArrayT pts(fNumNodes);
		for (int i = 0; i < fNumNodes; i++) pts[i] = fGlobalIDs[i];
		OutputSetT output_set(pts, n_labels, false);
		fOutputID = ElementSupport().RegisterOutput(output_set);
		fOutputNodesUsed = output_set.NodesUsed();
	}
}

void RKShellT::WriteOutput(void)
{
	const dArray2DT& disp = Field()[0];

	/* deflection summary (validation diagnostic); NaN-aware so an explicit blow-up is visible */
	double minuy = 0.0, maxmag = 0.0; bool blowup = false;
	for (int i = 0; i < fNumNodes; i++) {
		int g = fGlobalIDs[i];
		double uy = disp(g,1);
		double mag = std::sqrt(disp(g,0)*disp(g,0)+disp(g,1)*disp(g,1)+disp(g,2)*disp(g,2));
		if (uy != uy || mag != mag || mag > 1.0e30) { blowup = true; continue; }
		if (uy < minuy) minuy = uy;
		if (mag > maxmag) maxmag = mag;
	}
	if (blowup) { fprintf(stdout, "[RKShell] *** BLOW-UP (NaN/inf): unstable -- K not positive definite "
		"(hourglass) or dt too large ***\n"); fflush(stdout); return; }
	fprintf(stdout, "[RKShell] nodes=%d  min(u_y)=% .6e  max|u|=% .6e\n", fNumNodes, minuy, maxmag);

	/* reaction force at the monitored node (sum of every stencil's internal-force contribution
	 * there) vs its displacement -> the Fig 18 load-displacement curve */
	if (fMonitorNode > 0) {
		int mg0 = fMonitorNode - 1;
		int gtot = fGlobalToLocal.Length();
		/* TOTAL reaction summed over the load set: monitor_count nodes spaced monitor_stride apart
		 * (e.g. the whole driven generator line: start=first line node, stride=nt, count=line length).
		 * Default count=1 -> single node (the displacement x-axis uses the first monitor node). */
		double react[3] = {0,0,0};
		for (int l=0; l<fMonitorCount; l++){
			int mg = mg0 + l*fMonitorStride;
			if (mg < 0 || mg >= gtot) continue;
			for (int K=0; K<fNumNodes; K++){
				int nn = fNeighbors.MinorDim(K);
				if (nn < 6) continue;
				const int* gnb = fNeighbors(K);
				int pos = -1; for (int k=0;k<nn;k++) if (gnb[k]==mg){ pos=k; break; }
				if (pos < 0) continue;
				dArrayT ue(3*nn), f;
				for (int k=0;k<nn;k++) for (int d=0;d<3;d++) ue[k*3+d]=disp(gnb[k],d);
				if (fFiniteStrain) InternalForceFS(K, ue, f, false, false);  /* exclude bending penalty from reported reaction */
				else               InternalForce(K, ue, f, false);
				for (int d=0;d<3;d++) react[d] += f[pos*3+d];
			}
		}
		/* cleanly parseable: $2=node $3=ux $4=uy $5=uz $6=Rx $7=Ry $8=Rz (Rx = TOTAL over load set) */
		fprintf(stdout, "[RKShell-react] %d %.8e %.8e %.8e %.8e %.8e %.8e\n",
			fMonitorNode, disp(mg0,0),disp(mg0,1),disp(mg0,2), react[0],react[1],react[2]);
	}

	/* kinetic energy + peak nodal speed: quasi-static health. KE should be tiny vs the deformation
	 * work, and max|v| ~ the loading rate; large KE / max|v| >> loading rate = dynamic ringing. */
	if (Field().Order() >= 1) {
		const dArray2DT& vel = Field()[1];
		double KE=0.0, v2max=0.0;
		for (int i=0;i<fNumNodes;i++){
			int g=fGlobalIDs[i];
			double v2=vel(g,0)*vel(g,0)+vel(g,1)*vel(g,1)+vel(g,2)*vel(g,2);
			KE += 0.5*fLumpedMass[i]*v2;
			if (v2>v2max) v2max=v2;
		}
		fprintf(stdout, "[RKShell-energy] KE=%.6e  max|v|=%.6e\n", KE, std::sqrt(v2max));
	}
	fflush(stdout);

	/* write the displacement field (+ equivalent plastic strain when plasticity is active) */
	if (fOutputID < 0) return;
	int nf = (fYield > 0.0) ? ((fThicknessUpdate) ? 5 : 4) : 3;
	dArray2DT n_values(fOutputNodesUsed.Length(), nf);
	for (int k = 0; k < fOutputNodesUsed.Length(); k++) {
		int g = fOutputNodesUsed[k];
		for (int d = 0; d < 3; d++) n_values(k,d) = disp(g, d);
		int loc = (g < fGlobalToLocal.Length()) ? fGlobalToLocal[g] : -1;
		if (nf >= 4) {
			double ep = 0.0;
			if (loc >= 0 && loc < (int) fJ2ep.size())
				for (size_t q = 0; q < fJ2ep[loc].size(); q++) if (fJ2ep[loc][q] > ep) ep = fJ2ep[loc][q];
			n_values(k,3) = ep;
		}
		if (nf >= 5)
			n_values(k,4) = (loc >= 0 && loc < (int) fThicknessCur.size()) ? fThicknessCur[loc] : fThickness;
	}
	dArray2DT e_values; /* none */
	ElementSupport().WriteOutput(fOutputID, n_values, e_values);
}

void RKShellT::ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const
{
#pragma unused(connects_1)
	connects_2.Append(&fNeighbors);
}

void RKShellT::ConnectsX(AutoArrayT<const iArray2DT*>& connects) const
{
	/* expose the background cell mesh as geometry connectivity */
	for (int b = 0; b < fOutputConn.Length(); b++) connects.Append(fOutputConn[b]);
}

/* parameters */
void RKShellT::DefineParameters(ParameterListT& list) const
{
	ElementBaseT::DefineParameters(list);

	ParameterT thickness(fThickness, "shell_thickness");
	thickness.AddLimit(LimitT(0.0, LimitT::Lower));
	list.AddParameter(thickness);

	ParameterT young(fYoung, "Young_modulus");
	young.AddLimit(LimitT(0.0, LimitT::Lower));
	list.AddParameter(young);

	ParameterT poisson(fPoisson, "Poisson_ratio");
	list.AddParameter(poisson);

	ParameterT support(fSupportFac, "support_factor");
	support.SetDefault(3.0);
	list.AddParameter(support);

	ParameterT complete(fCompleteness, "completeness");
	complete.SetDefault(3);
	list.AddParameter(complete);

	/* uniform per-area applied load (gravity / surface pressure components) */
	ParameterT lx(fLoad[0], "load_x"); lx.SetDefault(0.0); list.AddParameter(lx);
	ParameterT ly(fLoad[1], "load_y"); ly.SetDefault(0.0); list.AddParameter(ly);
	ParameterT lz(fLoad[2], "load_z"); lz.SetDefault(0.0); list.AddParameter(lz);

	/* mass density (lumped mass for explicit dynamics; scale up for dynamic relaxation) */
	ParameterT density(fDensity, "density");
	density.AddLimit(LimitT(0.0, LimitT::Lower));
	density.SetDefault(1.0);
	list.AddParameter(density);

	/* stabilization (paper section 5): 0=default(membrane+bending), 3=alpha-scaled membrane
	 * (5.3), 2=pure bending (5.2) */
	ParameterT stab(fStabMode, "stabilization");
	stab.SetDefault(0);
	list.AddParameter(stab);
	ParameterT sm(fStabMembrane, "stab_membrane"); sm.SetDefault(1.0); list.AddParameter(sm);
	ParameterT sb(fStabBending,  "stab_bending");  sb.SetDefault(0.0); list.AddParameter(sb);
	ParameterT sn(fStabNatural,  "stab_natural"); sn.SetDefault(1.0); list.AddParameter(sn);

	/* plane-stress J2 plasticity (0 yield = elastic): Y(ep) = yield_stress + hardening*ep */
	ParameterT yld(fYield, "yield_stress"); yld.SetDefault(0.0); list.AddParameter(yld);
	ParameterT hrd(fHardening, "hardening_modulus"); hrd.SetDefault(0.0); list.AddParameter(hrd);
	ParameterT ysat(fYieldSat, "yield_saturation"); ysat.SetDefault(0.0); list.AddParameter(ysat);
	ParameterT srat(fSatRate, "saturation_rate"); srat.SetDefault(0.0); list.AddParameter(srat);

	/* finite-deformation kinematics (Green-Lagrange, current-config geometry); needed for the
	 * large-displacement elasto-plastic buckling (Fig 18). 0 = small-strain linear. */
	ParameterT fs(fFiniteStrain, "finite_strain"); fs.SetDefault(0); list.AddParameter(fs);
	ParameterT tu(fThicknessUpdate, "thickness_update"); tu.SetDefault(0); list.AddParameter(tu);
	ParameterT co(fCorotational, "corotational"); co.SetDefault(0); list.AddParameter(co);

	/* report the reaction force at this (1-based global) node each output step -> Fig 18 curve */
	ParameterT mn(fMonitorNode, "monitor_node"); mn.SetDefault(0); list.AddParameter(mn);
	/* sum the reaction over monitor_count nodes spaced monitor_stride apart (driven load line) */
	ParameterT mst(fMonitorStride, "monitor_stride"); mst.SetDefault(1); list.AddParameter(mst);
	ParameterT mct(fMonitorCount,  "monitor_count");  mct.SetDefault(1); list.AddParameter(mct);
	ParameterT dmp(fDamping, "damping"); dmp.SetDefault(0.0); list.AddParameter(dmp);
}

void RKShellT::TakeParameterList(const ParameterListT& list)
{
	/* base setup: resolves the field + integrator, collects block info, calls DefineElements */
	ElementBaseT::TakeParameterList(list);

	/* shell + meshfree parameters */
	fThickness    = list.GetParameter("shell_thickness");
	fYoung        = list.GetParameter("Young_modulus");
	fPoisson      = list.GetParameter("Poisson_ratio");
	fSupportFac   = list.GetParameter("support_factor");
	fCompleteness = list.GetParameter("completeness");
	fLoad[0] = list.GetParameter("load_x");
	fLoad[1] = list.GetParameter("load_y");
	fLoad[2] = list.GetParameter("load_z");
	fDensity = list.GetParameter("density");
	fStabMode     = list.GetParameter("stabilization");
	fStabMembrane = list.GetParameter("stab_membrane");
	fStabBending  = list.GetParameter("stab_bending");
	fStabNatural  = list.GetParameter("stab_natural");
	fYield        = list.GetParameter("yield_stress");
	fHardening    = list.GetParameter("hardening_modulus");
	fYieldSat     = list.GetParameter("yield_saturation");
	fSatRate      = list.GetParameter("saturation_rate");
	fFiniteStrain = list.GetParameter("finite_strain");
	fThicknessUpdate = list.GetParameter("thickness_update");
	fCorotational = list.GetParameter("corotational");
	fMonitorNode  = list.GetParameter("monitor_node");
	fMonitorStride = list.GetParameter("monitor_stride");
	fMonitorCount  = list.GetParameter("monitor_count");
	fDamping      = list.GetParameter("damping");

	/* plane-stress (sigma33=0) isotropic tangent; condensation acts in the LOCAL shell-normal
	 * frame at assembly time via ToVoigtLocal */
	double lam = fYoung*fPoisson/((1.0+fPoisson)*(1.0-2.0*fPoisson));
	double mu  = fYoung/(2.0*(1.0+fPoisson));
	for (int a=0;a<6;a++) for (int b=0;b<6;b++) fC[a][b] = 0.0;
	for (int a=0;a<3;a++) for (int b=0;b<3;b++) fC[a][b] = lam + (a==b ? 2.0*mu : 0.0);
	for (int a=3;a<6;a++) fC[a][a] = mu;
	double c22 = fC[2][2];
	double Cc[6][6];
	for (int a=0;a<6;a++) for (int b=0;b<6;b++) Cc[a][b] = fC[a][b] - fC[a][2]*fC[2][b]/c22;
	for (int a=0;a<6;a++) { Cc[a][2]=0.0; Cc[2][a]=0.0; }
	for (int a=0;a<6;a++) for (int b=0;b<6;b++) fC[a][b] = Cc[a][b];

	/* RKPM shape-function solver in the 2D local chart */
	dArrayT gwin(3);
	gwin[0] = (fCompleteness >= 3) ? 1.8 : 1.5;
	gwin[1] = 0.4;
	gwin[2] = 3.0;
	fMLS = new MLSSolverT(2, fCompleteness, false, MeshFreeT::kGaussian, gwin);
	fMLS->Initialize();

	/* precompute the per-node stencil stiffness (linear elastic, fixed geometry) */
	BuildElementStiffness();

	/* optional stabilization unit test (env KLSHELL_SELFTEST=1) */
	if (getenv("KLSHELL_SELFTEST")) RunStabSelfTest();

	/* lumped nodal mass (for the explicit central-difference solver) */
	BuildLumpedMass();

	/* background cell connectivity: exposed via ConnectsX so the framework's node-element graph
	 * (and the explicit nodal-update node set) sees this element's nodes as geometry */
	{
		ModelManagerT& model = ElementSupport().ModelManager();
		const ArrayT<StringT>& ids = model.ElementGroupIDs();
		fOutputConn.Dimension(ids.Length());
		for (int b = 0; b < ids.Length(); b++) {
			model.ReadConnectivity(ids[b]);
			fOutputConn[b] = model.ElementGroupPointer(ids[b]);
		}
	}
}

/* lumped nodal mass m_I = rho * A_I * h (diagonal) */
void RKShellT::BuildLumpedMass(void)
{
	fLumpedMass.Dimension(fNumNodes);
	for (int i = 0; i < fNumNodes; i++)
		fLumpedMass[i] = fDensity*fNodalArea[i]*fThickness;

	/* per-node current thickness (Algorithm 3 D33 accumulation); starts at the reference thickness */
	fThicknessCur.assign(fNumNodes, fThickness);

	/* co-rotational frames (Algorithm 2): R = reference tangent frame [E1 E2 N0] (columns), V = I */
	fFrameR.assign((size_t)fNumNodes*9, 0.0);
	fFrameV.assign((size_t)fNumNodes*9, 0.0);
	for (int i=0;i<fNumNodes;i++){
		double* R=&fFrameR[(size_t)i*9]; double* V=&fFrameV[(size_t)i*9];
		V[0]=V[4]=V[8]=1.0;
		R[0]=R[4]=R[8]=1.0;
		if ((int)fXref[i].size()>=15){
			const double* X=&fXref[i][0];
			double X1[3]={X[0],X[1],X[2]}, X2[3]={X[3],X[4],X[5]}, N0[3];
			Cross(X1,X2,N0); double nl=Norm(N0);
			if (nl>1.0e-30){ for(int d=0;d<3;d++) N0[d]/=nl;
				double E1[3],E2[3]; OrthoTangents(N0,E1,E2);
				for(int d=0;d<3;d++){ R[d*3+0]=E1[d]; R[d*3+1]=E2[d]; R[d*3+2]=N0[d]; }
			}
		}
	}
}

/* stabilization unit test: strain energy E = sum_K u_K^T fKe_K u_K of unit-norm modes */
void RKShellT::RunStabSelfTest(void)
{
	/* node spacing (for the checkerboard parity) */
	double hsp = 1.0e30;
	for (int j=1;j<fNumNodes;j++){ double d2=0; for(int k=0;k<3;k++){double dx=fCoords(j,k)-fCoords(0,k);d2+=dx*dx;}
		if (d2>1.0e-12 && d2<hsp) hsp=d2; }
	hsp = std::sqrt(hsp);

	const int NM=9;
	const char* names[NM]={"trans_x","trans_z","rot_z","stretch_x","shear_xy",
	                       "quad_mem_x","quad_bend_z","MEMBRANE_hg","BENDING_hg"};
	double E[NM];
	for (int m=0;m<NM;m++){
		std::vector<double> U(3*fNumNodes,0.0);
		double nrm=0.0;
		for (int i=0;i<fNumNodes;i++){
			double x=fCoords(i,0), y=fCoords(i,1);
			int px=(int)std::floor(x/hsp+0.5), py=(int)std::floor(y/hsp+0.5);
			double s=((px+py)&1)? -1.0 : 1.0;            /* checkerboard parity */
			double u[3]={0,0,0};
			switch(m){
				case 0: u[0]=1.0; break;                  /* rigid translation x   */
				case 1: u[2]=1.0; break;                  /* rigid translation z   */
				case 2: u[0]=-y; u[1]=x; break;           /* rigid rotation about z */
				case 3: u[0]=x; break;                    /* linear stretch (eps_xx)*/
				case 4: u[0]=y; break;                    /* linear shear           */
				case 5: u[0]=0.5*x*x; break;              /* quadratic membrane (resolved strain gradient) */
				case 6: u[2]=0.5*x*x; break;              /* quadratic bending (resolved CONSTANT curvature) */
				case 7: u[0]=s; break;                    /* in-plane MEMBRANE hourglass */
				case 8: u[2]=s; break;                    /* out-of-plane BENDING hourglass */
			}
			for(int d=0;d<3;d++){ U[3*i+d]=u[d]; nrm+=u[d]*u[d]; }
		}
		nrm=std::sqrt(nrm); if(nrm<1.0e-30) nrm=1.0;
		for (size_t k=0;k<U.size();k++) U[k]/=nrm;
		double En=0.0;
		for (int K=0;K<fNumNodes;K++){
			int nn=fNeighbors.MinorDim(K);
			if (nn<6) continue;
			const int* gnb=fNeighbors(K);
			std::vector<double> ue(3*nn);
			for(int k=0;k<nn;k++){ int loc=fGlobalToLocal[gnb[k]]; for(int d=0;d<3;d++) ue[3*k+d]=U[3*loc+d]; }
			const dMatrixT& Ke=fKe[K];
			for(int a=0;a<3*nn;a++){ double sa=0; for(int b=0;b<3*nn;b++) sa+=Ke(a,b)*ue[b]; En+=ue[a]*sa; }
		}
		E[m]=En;
	}
	fprintf(stdout,"\n=== STAB SELF-TEST  stab_membrane=%.4g stab_bending=%.4g  (unit-norm mode energies) ===\n",
		fStabMembrane,fStabBending);
	for (int m=0;m<NM;m++) fprintf(stdout,"   %-12s  E = % .6e\n",names[m],E[m]);
	fprintf(stdout,"   expect: rigid ~0 ; linear physical & stab-invariant ; hourglass ~0 w/o stab, >0 if caught\n\n");

	/* CURVATURE UNIT TEST (Algorithm 1 geometry): from the reference mid-surface derivatives in fXref,
	 * second fundamental form b_ab = x,ab . n, metric g_ab = x,a . x,b; principal curvatures = eig(g^-1 b).
	 * Cylinder of radius R -> max principal curvature = 1/R, min ~ 0 (flat -> both ~0). */
	{
		double ksum=0.0, kmin=1e30, kmax=0.0; int kn=0;
		for (int i=0;i<fNumNodes;i++){
			if ((int)fXref[i].size()<15) continue;
			const double* X=&fXref[i][0];
			double x1[3]={X[0],X[1],X[2]}, x2[3]={X[3],X[4],X[5]};
			double x11[3]={X[6],X[7],X[8]}, x22[3]={X[9],X[10],X[11]}, x12[3]={X[12],X[13],X[14]};
			double nv[3]; Cross(x1,x2,nv); double nl=Norm(nv); if (nl<1.0e-30) continue;
			for (int d=0;d<3;d++) nv[d]/=nl;
			double b11=Dot(x11,nv), b22=Dot(x22,nv), b12=Dot(x12,nv);
			double g11=Dot(x1,x1), g22=Dot(x2,x2), g12=Dot(x1,x2);
			double detg=g11*g22-g12*g12; if (std::fabs(detg)<1.0e-30) continue;
			double Hc=(b11*g22-2.0*b12*g12+b22*g11)/(2.0*detg);   /* mean curvature */
			double Kc=(b11*b22-b12*b12)/detg;                      /* Gaussian curvature */
			double disc=Hc*Hc-Kc; if (disc<0) disc=0;
			double k1=std::fabs(Hc+std::sqrt(disc)), k2=std::fabs(Hc-std::sqrt(disc));
			double km=(k1>k2)?k1:k2;
			ksum+=km; if(km<kmin)kmin=km; if(km>kmax)kmax=km; kn++;
			/* geometry-accuracy study (paper 4.1): per-node signed principal curvatures + normal */
			if (getenv("KLSHELL_VASE") && i<fCoords.MajorDim())
				fprintf(stdout,"[CURV-NODE] %.8e %.8e %.8e %.10e %.10e %.10e %.10e %.10e\n",
					fCoords(i,0),fCoords(i,1),fCoords(i,2), Hc+std::sqrt(disc), Hc-std::sqrt(disc), nv[0],nv[1],nv[2]);
		}
		if (kn>0) fprintf(stdout,"[CURVATURE-TEST] max-principal-curvature: mean=%.6e min=%.6e max=%.6e  (cylinder R -> 1/R)\n\n",
			ksum/kn, kmin, kmax);
	}
	fflush(stdout);
}

/* meshfree node set + neighbor lists (overrides the FE block element setup) */
void RKShellT::DefineElements(const ArrayT<StringT>& block_ID, const ArrayT<int>& mat_index)
{
#pragma unused(mat_index)
	ModelManagerT& model = ElementSupport().ModelManager();
	const dArray2DT& all_coords = model.Coordinates();
	int nsd = all_coords.MinorDim();
	if (nsd != 3)
		ExceptionT::GeneralFail("RKShellT::DefineElements",
			"KL shell requires 3D coordinates, got %d", nsd);

	/* gather the shell node set: the nodes referenced by the declared background block(s), or
	 * ALL model nodes when no block is given (single-shell deck) */
	int num_global = all_coords.MajorDim();
	iArrayT used(num_global);
	if (block_ID.Length() == 0)
		used = 1;
	else {
		used = 0;
		for (int b = 0; b < block_ID.Length(); b++) {
			const iArray2DT& conn = model.ElementGroup(block_ID[b]);
			const int* p = conn.Pointer();
			for (int k = 0; k < conn.Length(); k++) used[p[k]] = 1;
		}
	}

	/* global node id <-> local shell index */
	fGlobalToLocal.Dimension(num_global);
	fGlobalToLocal = -1;
	fNumNodes = 0;
	for (int n = 0; n < num_global; n++) if (used[n]) fGlobalToLocal[n] = fNumNodes++;

	fCoords.Dimension(fNumNodes, 3);
	fGlobalIDs.Dimension(fNumNodes);
	for (int n = 0; n < num_global; n++)
		if (fGlobalToLocal[n] >= 0) {
			int l = fGlobalToLocal[n];
			for (int d = 0; d < 3; d++) fCoords(l, d) = all_coords(n, d);
			fGlobalIDs[l] = n;
		}

	/* one element card per integration node (meshfree nodal integration); the card carries
	 * per-node state storage. (Mirrors SCNIMFT's one-card-per-node setup.) */
	fElementCards.Dimension(fNumNodes);
	for (int i = 0; i < fNumNodes; i++) fElementCards[i].SetMaterialNumber(0);

	/* build neighbor lists (local indices) + nodal areas */
	BuildNeighbors();

	/* convert neighbor lists to GLOBAL node ids for the framework connectivity */
	ArrayT<int> rowcounts(fNumNodes);
	for (int i = 0; i < fNumNodes; i++) rowcounts[i] = fNeighbors.MinorDim(i);
	RaggedArray2DT<int> global_neighbors;
	global_neighbors.Configure(rowcounts);
	for (int i = 0; i < fNumNodes; i++) {
		int m = fNeighbors.MinorDim(i);
		const int* loc = fNeighbors(i);
		int* glob = global_neighbors(i);
		for (int k = 0; k < m; k++) glob[k] = fGlobalIDs[loc[k]];
	}
	fNeighbors = global_neighbors;
}

/* per-node neighbor lists (3D distance, local indices) + nodal areas */
void RKShellT::BuildNeighbors(void)
{
	double hmin = 1.0e30;
	for (int j = 1; j < fNumNodes; j++) {
		double d2 = 0.0;
		for (int d = 0; d < 3; d++) { double dx = fCoords(j,d)-fCoords(0,d); d2 += dx*dx; }
		if (d2 < hmin) hmin = d2;
	}
	double spacing = std::sqrt(hmin);
	double support = fSupportFac*spacing;

	ArrayT<int> counts(fNumNodes);
	AutoArrayT<int> flat;
	for (int i = 0; i < fNumNodes; i++) {
		int c = 0;
		for (int j = 0; j < fNumNodes; j++) {
			double d2 = 0.0;
			for (int d = 0; d < 3; d++) { double dx = fCoords(j,d)-fCoords(i,d); d2 += dx*dx; }
			if (std::sqrt(d2) < 0.99*support) { flat.Append(j); c++; }
		}
		counts[i] = c;
	}
	fNeighbors.Configure(counts);
	int pos = 0;
	for (int i = 0; i < fNumNodes; i++) {
		int* row = fNeighbors(i);
		for (int k = 0; k < counts[i]; k++) row[k] = flat[pos++];
	}

	fNodalArea.Dimension(fNumNodes);
	fNodalArea = spacing*spacing;
}

/* precompute per-node stencil stiffness via the validated KL-shell kernels (linear elastic) */
void RKShellT::BuildElementStiffness(void)
{
	fKe.Dimension(fNumNodes);
	fIPB.assign(fNumNodes, std::vector<double>());
	fIPw.assign(fNumNodes, std::vector<double>());
	fIPstab.assign(fNumNodes, std::vector<char>());
	/* plane-stress J2 history: 3 through-thickness base points per node (over-allocated; the base
	 * loop appends exactly the BuildGeom-valid points, indexed in order by InternalForce) */
	fJ2sig.assign(fNumNodes, std::vector<double>(9, 0.0));
	fJ2ep.assign(fNumNodes, std::vector<double>(3, 0.0));
	fJ2eps.assign(fNumNodes, std::vector<double>(9, 0.0));
	fDphi.assign(fNumNodes, std::vector<double>());
	fXref.assign(fNumNodes, std::vector<double>());
	fBendR.assign(fNumNodes, std::vector<double>());
	fBendN.assign(fNumNodes, std::vector<double>(3,0.0));
	fBendCoeff.assign(fNumNodes, 0.0);
	double h = fThickness;

	double hmin = 1.0e30;
	for (int j = 1; j < fNumNodes; j++) {
		double d2 = 0.0;
		for (int d = 0; d < 3; d++) { double dx = fCoords(j,d)-fCoords(0,d); d2 += dx*dx; }
		if (d2 < hmin) hmin = d2;
	}
	double support = fSupportFac*std::sqrt(hmin);

	double xg[3] = {-std::sqrt(3.0/5.0), 0.0, std::sqrt(3.0/5.0)};
	double wg[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0};

	for (int i = 0; i < fNumNodes; i++) {
		int nn = fNeighbors.MinorDim(i);
		const int* gnb = fNeighbors(i);
		fKe[i].Dimension(3*nn);
		fKe[i] = 0.0;
		if (nn < 6) continue;

		ArrayT<int> loc(nn);
		std::vector<double> nbX(3*nn);
		for (int k = 0; k < nn; k++) {
			loc[k] = fGlobalToLocal[gnb[k]];
			for (int d = 0; d < 3; d++) nbX[3*k+d] = fCoords(loc[k], d);
		}
		double psi1[3], psi2[3], n0[3];
		PCAFrame(&nbX[0], nn, psi1, psi2, n0);

		dArray2DT lc(nn, 2);
		for (int k = 0; k < nn; k++) {
			double dv[3] = { fCoords(loc[k],0)-fCoords(i,0), fCoords(loc[k],1)-fCoords(i,1),
			                 fCoords(loc[k],2)-fCoords(i,2) };
			lc(k,0) = Dot(dv, psi1); lc(k,1) = Dot(dv, psi2);
		}
		dArray2DT np(nn,1); np = support;
		dArrayT vol(nn); for (int k=0;k<nn;k++) vol[k] = fNodalArea[loc[k]];
		dArrayT sample(2); sample[0]=0.0; sample[1]=0.0;
		if (!fMLS->SetField(lc, np, vol, sample, 3)) continue;
		const dArray2DT& Dp = fMLS->Dphi();
		const dArray2DT& DDp = fMLS->DDphi();
		const dArray2DT& DDDp = fMLS->DDDphi();

		double x1[3]={0,0,0},x2[3]={0,0,0},x11[3]={0,0,0},x22[3]={0,0,0},x12[3]={0,0,0};
		for (int I=0;I<nn;I++) {
			double Xq[3]={fCoords(loc[I],0),fCoords(loc[I],1),fCoords(loc[I],2)};
			for (int d=0;d<3;d++){
				x1[d]+=Dp(0,I)*Xq[d]; x2[d]+=Dp(1,I)*Xq[d];
				x11[d]+=DDp(0,I)*Xq[d]; x22[d]+=DDp(1,I)*Xq[d]; x12[d]+=DDp(2,I)*Xq[d];
			}
		}
		/* store stencil shape derivatives + reference mid-surface derivatives (finite-strain path) */
		fDphi[i].resize((size_t)nn*5);
		for (int I=0;I<nn;I++){ fDphi[i][I*5]=Dp(0,I); fDphi[i][I*5+1]=Dp(1,I);
			fDphi[i][I*5+2]=DDp(0,I); fDphi[i][I*5+3]=DDp(1,I); fDphi[i][I*5+4]=DDp(2,I); }
		fXref[i].resize(15);
		for (int d=0;d<3;d++){ fXref[i][d]=x1[d]; fXref[i][3+d]=x2[d];
			fXref[i][6+d]=x11[d]; fXref[i][9+d]=x22[d]; fXref[i][12+d]=x12[d]; }

		double A_K = fNodalArea[i];
		double V_K = A_K*h;
		double cell = std::sqrt(A_K);
		double Mmom = cell*cell/12.0;
		dMatrixT& Ke = fKe[i];

		/* in-plane characteristic length h_pl (mean neighbor distance) -> alpha = min(1, h/h_pl)
		 * for the section-5.3 limited-membrane-stabilization scaling */
		double h_pl = 0.0; int npc = 0;
		for (int k=0;k<nn;k++){ if(loc[k]==i) continue;
			double dd=0.0; for(int d=0;d<3;d++){double dx=fCoords(loc[k],d)-fCoords(i,d); dd+=dx*dx;}
			h_pl += std::sqrt(dd); npc++; }
		if (npc>0) h_pl /= npc;
		double alpha = (h_pl > 0.0 && h < h_pl) ? (h/h_pl) : 1.0;

		/* nodal integration (3-pt thru-thickness Gauss), local-frame plane stress */
		for (int g=0;g<3;g++) {
			ShellGeom G;
			if (!BuildGeom(x1,x2,x11,x22,x12,h,xg[g],G)) continue;
			double e1[3],e2[3]; OrthoTangents(G.n,e1,e2);
			double cw = wg[g]*(h/2.0)*A_K;
			std::vector<std::vector<double> > Bv(nn, std::vector<double>(18));
			for (int I=0;I<nn;I++){
				double B[3][3][3];
				BMatrix(G,Dp(0,I),Dp(1,I),DDp(0,I),DDp(2,I),DDp(1,I),B);
				double bv[6][3]; ToVoigtLocal(B,e1,e2,G.n,bv);
				for(int r=0;r<6;r++)for(int c=0;c<3;c++) Bv[I][r*3+c]=bv[r][c];
			}
			/* store this base (material) integration point for the stress-driven internal force */
			{ std::vector<double> Bf((size_t)6*3*nn);
			  for(int I=0;I<nn;I++)for(int r=0;r<6;r++)for(int c=0;c<3;c++) Bf[(size_t)r*3*nn+3*I+c]=Bv[I][r*3+c];
			  fIPB[i].insert(fIPB[i].end(),Bf.begin(),Bf.end());
			  fIPw[i].push_back(cw); fIPstab[i].push_back(0); }
			for (int I=0;I<nn;I++) for (int J=0;J<nn;J++){
				double kij[3][3]={{0,0,0},{0,0,0},{0,0,0}};
				for(int a=0;a<6;a++)for(int b=0;b<6;b++){double Cab=fC[a][b];if(Cab==0.0)continue;
					for(int ci=0;ci<3;ci++)for(int cj=0;cj<3;cj++) kij[ci][cj]+=Bv[I][a*3+ci]*Cab*Bv[J][b*3+cj];}
				for(int ci=0;ci<3;ci++)for(int cj=0;cj<3;cj++) Ke(3*I+ci,3*J+cj)+=cw*kij[ci][cj];
			}
		}

		/* NATURAL (Taylor-gradient) stabilization at xi3=0 -- the paper's Eq. 33 CONSISTENT stabilizer
		 * (replaces the curvature penalty). K_stab = sum_l (B_,xil)^T C (B_,xil) * V_K * Mmom * alpha,
		 * using only 2nd shape derivatives (the xi3 3rd-derivative term vanishes at xi3=0). It vanishes
		 * on smooth fields (Scordelis-Lo stays bit-exact) but fires on node-to-node hourglass content;
		 * scaled by Mmom~s^2 it steps out of the way of sub-grid physical folds -> no force inflation. */
		if (fStabNatural > 0.0) {
			ShellGeom G0;
			if (BuildGeom(x1,x2,x11,x22,x12,h,0.0,G0)) {
				double e1[3],e2[3]; OrthoTangents(G0.n,e1,e2);
				double Vmom = fStabNatural * (h*A_K) * Mmom * alpha;   /* V_K * (s^2/12) * alpha */
				std::vector<std::vector<double> > Bg(nn, std::vector<double>(36,0.0)); /* [l*18+r*3+c] */
				for (int I=0;I<nn;I++){
					double B[3][3][3];
					BMatrix(G0,Dp(0,I),Dp(1,I),DDp(0,I),DDp(2,I),DDp(1,I),B);
					double P1l[2]={DDp(0,I),DDp(2,I)};   /* Psi,1,xil = {Psi,11, Psi,12} */
					double P2l[2]={DDp(2,I),DDp(1,I)};   /* Psi,2,xil = {Psi,12, Psi,22} */
					double Bgt[3][3][3][2];
					BMatrixGradient(G0,Dp(0,I),Dp(1,I),P1l,P2l,B,Bgt);
					for (int l=0;l<2;l++){
						double Bgl[3][3][3];
						for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int c=0;c<3;c++) Bgl[a][b][c]=Bgt[a][b][c][l];
						double bv[6][3]; ToVoigtLocal(Bgl,e1,e2,G0.n,bv);
						for(int r=0;r<6;r++)for(int c=0;c<3;c++) Bg[I][l*18+r*3+c]=bv[r][c];
					}
				}
				/* store the two B_,xil operators as elastic stab points (fIPstab=1, weight Vmom) so the
				 * existing SCNI force loop in InternalForce/FS assembles the explicit stabilization force
				 * f = sum_l (B_,xil)^T C (B_,xil) u * Vmom automatically -- no separate force code needed */
				for (int l=0;l<2;l++){
					std::vector<double> Bf((size_t)6*3*nn);
					for(int I=0;I<nn;I++)for(int r=0;r<6;r++)for(int c=0;c<3;c++) Bf[(size_t)r*3*nn+3*I+c]=Bg[I][l*18+r*3+c];
					fIPB[i].insert(fIPB[i].end(),Bf.begin(),Bf.end());
					fIPw[i].push_back(Vmom); fIPstab[i].push_back(1);
				}
				for (int I=0;I<nn;I++) for (int J=0;J<nn;J++){
					double kij[3][3]={{0,0,0},{0,0,0},{0,0,0}};
					for (int l=0;l<2;l++)
						for(int a=0;a<6;a++)for(int b=0;b<6;b++){ double Cab=fC[a][b]; if(Cab==0.0) continue;
							for(int ci=0;ci<3;ci++)for(int cj=0;cj<3;cj++)
								kij[ci][cj]+=Bg[I][l*18+a*3+ci]*Cab*Bg[J][l*18+b*3+cj]; }
					for(int ci=0;ci<3;ci++)for(int cj=0;cj<3;cj++) Ke(3*I+ci,3*J+cj)+=Vmom*kij[ci][cj];
				}
			}
		}

		/* BENDING-hourglass control, FREQUENCY-SELECTIVE: penalize R_I = D_chord_I - D_LS_I.
		 * D_chord (4/L^2 chord Laplacian) catches ALL node-to-node curvature -- the spurious sawtooth
		 * AND the physical resolved curvature (dimple/ovalization). D_LS is a moment-matched 3x3
		 * least-squares Laplacian: it REPRODUCES a resolved quadratic curvature field exactly but
		 * SMOOTHS the node-to-node sawtooth to ~0. So on resolved/physical bending D_chord==D_LS ->
		 * R~0 (NO penalty, no force inflation); on the hourglass D_LS~0 while D_chord is large -> R
		 * large -> penalty preserved. This fixes the chord operator's over-penalization of the sharp
		 * physical load dimple (which was inflating the reaction ~10x) without losing hourglass control. */
		if (fStabBending > 0.0) {
			double nrm[3]; Cross(x1,x2,nrm); double Jn=Norm(nrm);
			if (Jn > 1.0e-300) {
				for (int d=0;d<3;d++) nrm[d]/=Jn;
				/* chord curvature operator (catches everything, incl. the sawtooth) */
				std::vector<double> Dop(nn,0.0); double W=0.0;
				for (int J=0;J<nn;J++){
					double L2=lc(J,0)*lc(J,0)+lc(J,1)*lc(J,1);
					if (L2 < 1.0e-12) continue;          /* self node (origin of the chart) */
					Dop[J]=4.0/L2; W+=1.0;
				}
				double Dself=0.0;
				if (W>0.0) for (int J=0;J<nn;J++){ Dop[J]/=W; Dself+=Dop[J]; }
				/* moment-matched LS Laplacian: kappa = M^-1 sum p_J (u_J-u_K).n with quadratic basis
				 * p_J=[1/2 x^2, x y, 1/2 y^2]; Laplacian coeff C_LS_J = (M^-1 p_J)[0] + (M^-1 p_J)[2] */
				double M[3][3]={{0,0,0},{0,0,0},{0,0,0}};
				std::vector<double> px(nn),py(nn),pz(nn);
				for (int J=0;J<nn;J++){ double x=lc(J,0),y=lc(J,1); px[J]=0.5*x*x; py[J]=x*y; pz[J]=0.5*y*y;
					M[0][0]+=px[J]*px[J]; M[0][1]+=px[J]*py[J]; M[0][2]+=px[J]*pz[J];
					M[1][1]+=py[J]*py[J]; M[1][2]+=py[J]*pz[J]; M[2][2]+=pz[J]*pz[J]; }
				M[1][0]=M[0][1]; M[2][0]=M[0][2]; M[2][1]=M[1][2];
				double det=M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
				          -M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
				          +M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]);
				std::vector<double> Cls(nn,0.0); double Clsself=0.0;
				if (std::fabs(det) > 1.0e-300) {        /* singular -> collinear stencil; fall back to chord (Cls=0) */
					double id=1.0/det;
					double Mi[3][3];
					Mi[0][0]=(M[1][1]*M[2][2]-M[1][2]*M[2][1])*id; Mi[0][1]=(M[0][2]*M[2][1]-M[0][1]*M[2][2])*id; Mi[0][2]=(M[0][1]*M[1][2]-M[0][2]*M[1][1])*id;
					Mi[2][0]=(M[1][0]*M[2][1]-M[1][1]*M[2][0])*id; Mi[2][1]=(M[0][1]*M[2][0]-M[0][0]*M[2][1])*id; Mi[2][2]=(M[0][0]*M[1][1]-M[0][1]*M[1][0])*id;
					for (int J=0;J<nn;J++){
						double q0=Mi[0][0]*px[J]+Mi[0][1]*py[J]+Mi[0][2]*pz[J];
						double q2=Mi[2][0]*px[J]+Mi[2][1]*py[J]+Mi[2][2]*pz[J];
						Cls[J]=q0+q2; Clsself+=Cls[J];
					}
				}
				std::vector<double> Rk(nn);
				for (int I=0;I<nn;I++){
					double L2=lc(I,0)*lc(I,0)+lc(I,1)*lc(I,1);
					double Dch=(L2 < 1.0e-12) ? -Dself   : Dop[I];
					double Dls=(L2 < 1.0e-12) ? -Clsself : Cls[I];
					Rk[I]=Dch-Dls;
				}
				double coeff=fStabBending*(fYoung*h*h*h/12.0)*A_K;
				for (int I=0;I<nn;I++) for (int Jp=0;Jp<nn;Jp++){
					double rr=coeff*Rk[I]*Rk[Jp];
					for (int a=0;a<3;a++) for (int b=0;b<3;b++) Ke(3*I+a,3*Jp+b)+=rr*nrm[a]*nrm[b];
				}
				/* store for the finite-strain force path (same rank-1 penalty on the total u) */
				fBendR[i].assign(Rk.begin(), Rk.end());
				for (int d=0;d<3;d++) fBendN[i][d]=nrm[d];
				fBendCoeff[i]=coeff;
			}
		}

		/* SCNI / NSNI cell-smoothed assumed-strain stabilization. For node K's square smoothing cell
		 * (side s = nodal spacing) in the PCA chart, the divergence-theorem cell-averaged gradient
		 * reduces to a CENTRAL DIFFERENCE of the shape functions sampled at the 4 edge midpoints:
		 *   d~Psi/dxi1 = [Psi(+e1) - Psi(-e1)]/s , and likewise the smoothed 2nd derivs from the
		 * boundary integral of the 1st derivs. The stabilization is the residual R = B_direct - B~
		 * (analytical point operator minus the cell-smoothed operator). For a smooth field the cell
		 * average matches the point value (divergence theorem) -> R ~ 0 -> membrane/bending energy
		 * UNPOLLUTED; the hourglass sawtooth has nonzero cell average but is invisible to the point
		 * sample -> R large -> a PSD penalty (R^T C R) that suppresses the mode. */
		{
			double s_cell = std::sqrt(hmin);              /* nodal spacing = smoothing-cell side */
			/* analytical shape derivs at K (copied before re-evaluating at the cell midpoints) */
			std::vector<double> P1K(nn),P2K(nn),P11K(nn),P22K(nn),P12K(nn);
			for (int I=0;I<nn;I++){ P1K[I]=Dp(0,I); P2K[I]=Dp(1,I);
				P11K[I]=DDp(0,I); P22K[I]=DDp(1,I); P12K[I]=DDp(2,I); }

			/* sample phi (values) and Dphi (1st derivs) at the 4 cell-edge midpoints (+/-e1, +/-e2) */
			double mids[4][2] = {{ s_cell/2,0},{-s_cell/2,0},{0, s_cell/2},{0,-s_cell/2}};
			std::vector<std::vector<double> > Phi(4, std::vector<double>(nn,0.0));
			std::vector<std::vector<double> > Dx(4, std::vector<double>(nn,0.0)), Dy(4, std::vector<double>(nn,0.0));
			bool okcell = true;
			for (int m=0;m<4 && okcell;m++){
				dArrayT sM(2); sM[0]=mids[m][0]; sM[1]=mids[m][1];
				if (!fMLS->SetField(lc, np, vol, sM, 3)) { okcell=false; break; }
				const dArrayT& ph=fMLS->phi(); const dArray2DT& dp=fMLS->Dphi();
				for (int I=0;I<nn;I++){ Phi[m][I]=ph[I]; Dx[m][I]=dp(0,I); Dy[m][I]=dp(1,I); }
			}

			if (okcell) {
				double inv_s = 1.0/s_cell;
				/* cell-smoothed 1st + 2nd shape derivatives (R=+e1, L=-e1, T=+e2, B=-e2 -> m=0,1,2,3) */
				std::vector<double> sm1(nn),sm2(nn),sm11(nn),sm22(nn),sm12(nn);
				for (int I=0;I<nn;I++){
					sm1[I]  = (Phi[0][I]-Phi[1][I])*inv_s;             /* d~Psi/dxi1 */
					sm2[I]  = (Phi[2][I]-Phi[3][I])*inv_s;             /* d~Psi/dxi2 */
					sm11[I] = (Dx[0][I]-Dx[1][I])*inv_s;              /* d~^2Psi/dxi1^2 */
					sm22[I] = (Dy[2][I]-Dy[3][I])*inv_s;              /* d~^2Psi/dxi2^2 */
					sm12[I] = 0.5*((Dx[2][I]-Dx[3][I]) + (Dy[0][I]-Dy[1][I]))*inv_s; /* symmetric mixed */
				}

				double afac = (fStabMode == 3) ? alpha : 1.0;
				for (int g=0; g<3; g++) {
					ShellGeom G;
					if (!BuildGeom(x1,x2,x11,x22,x12,h,xg[g],G)) continue;
					double e1[3],e2[3]; OrthoTangents(G.n,e1,e2);
					double cw = fStabMembrane*afac*wg[g]*(h/2.0)*A_K;
					if (cw == 0.0) continue;
					/* residual R_I = B_direct_I - B~_I (Voigt, in K's local frame) */
					std::vector<std::vector<double> > Rv(nn, std::vector<double>(18));
					for (int I=0;I<nn;I++){
						double Bd[3][3][3], Bs[3][3][3];
						BMatrix(G, P1K[I], P2K[I], P11K[I], P12K[I], P22K[I], Bd);
						BMatrix(G, sm1[I], sm2[I], sm11[I], sm12[I], sm22[I], Bs);
						double bvd[6][3], bvs[6][3];
						ToVoigtLocal(Bd,e1,e2,G.n,bvd); ToVoigtLocal(Bs,e1,e2,G.n,bvs);
						for(int r=0;r<6;r++)for(int c=0;c<3;c++) Rv[I][r*3+c]=bvd[r][c]-bvs[r][c];
					}
					/* store this stabilization (always-elastic) integration point */
					{ std::vector<double> Rf((size_t)6*3*nn);
					  for(int I=0;I<nn;I++)for(int r=0;r<6;r++)for(int c=0;c<3;c++) Rf[(size_t)r*3*nn+3*I+c]=Rv[I][r*3+c];
					  fIPB[i].insert(fIPB[i].end(),Rf.begin(),Rf.end());
					  fIPw[i].push_back(cw); fIPstab[i].push_back(1); }
					for (int I=0;I<nn;I++) for (int Jp=0;Jp<nn;Jp++){
						double kij[3][3]={{0,0,0},{0,0,0},{0,0,0}};
						for(int a=0;a<6;a++)for(int b=0;b<6;b++){double Cab=fC[a][b];if(Cab==0.0)continue;
							for(int ci=0;ci<3;ci++)for(int cj=0;cj<3;cj++) kij[ci][cj]+=Rv[I][a*3+ci]*Cab*Rv[Jp][b*3+cj];}
						for(int ci=0;ci<3;ci++)for(int cj=0;cj<3;cj++) Ke(3*I+ci,3*Jp+cj)+=cw*kij[ci][cj];
					}
				}
			}
		}
	}
}

/* stress-driven internal force f = sum_pt B^T sigma(B*ue) w. Elastic (sigma = fC*eps) -> reduces
 * exactly to fKe*ue; the base material points (fIPstab==0) are where the plane-stress J2 stress
 * update plugs in for the Fig 18 elasto-plastic track. */
void RKShellT::InternalForce(int i, const dArrayT& ue, dArrayT& fout, bool commit)
{
	int nn = fNeighbors.MinorDim(i);
	int ndof = 3*nn;
	fout.Dimension(ndof); fout = 0.0;
	int npt = (int) fIPw[i].size();
	if (npt == 0) return;
	const double* Ball = &fIPB[i][0];
	bool plastic = (fYield > 0.0);
	int nbase = (int) fJ2ep[i].size();
	int bp = 0;                                   /* through-thickness base-point counter */
	for (int p=0;p<npt;p++){
		const double* B = Ball + (size_t)p*6*ndof;
		/* strain eps = B*ue (local-frame Voigt; in-plane = [0]e11 [1]e22 [5]g12) */
		double eps[6];
		for (int r=0;r<6;r++){ const double* Br=B+(size_t)r*ndof; double s=0.0;
			for(int c=0;c<ndof;c++) s+=Br[c]*ue[c]; eps[r]=s; }

		double sig[6] = {0,0,0,0,0,0};
		if (plastic && fIPstab[i][p]==0 && bp<nbase) {
			/* base material point: plane-stress J2 on the in-plane components, per-point state.
			 * (sigma33 = sigma13 = sigma23 = 0: plane stress + Kirchhoff no-transverse-shear.) */
			double sip[3] = { fJ2sig[i][3*bp], fJ2sig[i][3*bp+1], fJ2sig[i][3*bp+2] };
			double ep = fJ2ep[i][bp];
			double deps[3] = { eps[0]-fJ2eps[i][3*bp], eps[1]-fJ2eps[i][3*bp+1], eps[5]-fJ2eps[i][3*bp+2] };
			PlaneStressJ2Return(sip, deps, ep, fYoung, fPoisson, fYield, fHardening, fYieldSat, fSatRate);
			sig[0]=sip[0]; sig[1]=sip[1]; sig[5]=sip[2];
			if (commit) {
				fJ2sig[i][3*bp]=sip[0]; fJ2sig[i][3*bp+1]=sip[1]; fJ2sig[i][3*bp+2]=sip[2];
				fJ2ep[i][bp]=ep;
				fJ2eps[i][3*bp]=eps[0]; fJ2eps[i][3*bp+1]=eps[1]; fJ2eps[i][3*bp+2]=eps[5];
			}
			bp++;
		} else {
			/* elastic point (stabilization residual, or elastic material): sigma = fC*eps */
			if (fIPstab[i][p]==0) bp++;           /* keep the base counter aligned in the elastic case */
			for (int r=0;r<6;r++){ double s=0.0; for(int cc=0;cc<6;cc++) s+=fC[r][cc]*eps[cc]; sig[r]=s; }
		}
		/* f += B^T sig * w */
		double w = fIPw[i][p];
		for (int c=0;c<ndof;c++){ double s=0.0;
			for(int r=0;r<6;r++) s+=B[(size_t)r*ndof+c]*sig[r]; fout[c]+=s*w; }
	}
}

/* finite-deformation internal force: objective Green-Lagrange strain E = 1/2(g - G) from the
 * CURRENT-config metric/curvature (g_ab = x,a . x,b is rigid-rotation invariant -> objective with
 * no co-rotational machinery), current-config B = dE/du via BMatrix on the deformed geometry, the
 * plane-stress stress (elastic or per-station J2 on the strain increment), plus the (reference-
 * config) SCNI stabilization for hourglass/explicit control. */
void RKShellT::InternalForceFS(int i, const dArrayT& ue, dArrayT& fout, bool commit, bool include_bend)
{
	int nn = fNeighbors.MinorDim(i);
	int ndof = 3*nn;
	fout.Dimension(ndof); fout = 0.0;
	if ((int) fDphi[i].size() != nn*5) return;
	double h = fThickness;
	const double* Dp = &fDphi[i][0];
	const double* Xr = &fXref[i][0];

	/* current mid-surface derivatives x,a = X,a(ref) + u,a (u,a = sum Dphi_I,a u_I) */
	double u1[3]={0,0,0},u2[3]={0,0,0},u11[3]={0,0,0},u22[3]={0,0,0},u12[3]={0,0,0};
	for (int I=0;I<nn;I++){
		double p1=Dp[I*5],p2=Dp[I*5+1],p11=Dp[I*5+2],p22=Dp[I*5+3],p12=Dp[I*5+4];
		for(int d=0;d<3;d++){ double u=ue[I*3+d];
			u1[d]+=p1*u; u2[d]+=p2*u; u11[d]+=p11*u; u22[d]+=p22*u; u12[d]+=p12*u; }
	}
	double X1[3],X2[3],X11[3],X22[3],X12[3], x1[3],x2[3],x11[3],x22[3],x12[3];
	for(int d=0;d<3;d++){
		X1[d]=Xr[d]; X2[d]=Xr[3+d]; X11[d]=Xr[6+d]; X22[d]=Xr[9+d]; X12[d]=Xr[12+d];
		x1[d]=X1[d]+u1[d]; x2[d]=X2[d]+u2[d]; x11[d]=X11[d]+u11[d]; x22[d]=X22[d]+u22[d]; x12[d]=X12[d]+u12[d];
	}

	/* THICKNESS UPDATE (Algorithm 3): use the per-node current thickness accumulated from the actual
	 * through-thickness strain D33 of the sigma33=0 update (set at the end of this routine on commit).
	 * The thinned section loses bending/membrane stiffness -> drives necking localization. */
	if (fThicknessUpdate && i < (int)fThicknessCur.size() && fThicknessCur[i] > 0.0)
		h = fThicknessCur[i];

	/* RATE form: strain INCREMENT deps = B_cur . du (du = u - u_prev) reuses the existing BMatrix
	 * strain measure at the current config -> consistent with the linear element + objective for
	 * small per-step rotations (slow mass-scaled loading). Stress is accumulated per station. */
	const int* gnb = fNeighbors(i);
	dArrayT due(ndof);
	for (int k=0;k<nn;k++){ int gg=gnb[k]; for(int d=0;d<3;d++) due[k*3+d]=ue[k*3+d]-fUprev(gg,d); }

	double A_K = fNodalArea[i];
	double xg[3] = {-std::sqrt(3.0/5.0), 0.0, std::sqrt(3.0/5.0)};
	double wg[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0};
	bool plastic = (fYield > 0.0);
	double c_ps = fYoung/(1.0 - fPoisson*fPoisson);
	double de33_mid = 0.0;   /* mid-surface through-thickness strain increment (Algorithm 3 thickness) */

	/* ALGORITHM 2 (co-rotational): advance the material frame R via Flanagan-Taylor from the
	 * velocity-gradient increment dL = du,a (x) g^a, re-align its normal to the geometric mid-normal,
	 * and rotate the stored in-plane stress into the new frame (objective). cr_e1/cr_e2 = tracked
	 * tangent directions used (instead of the arbitrary OrthoTangents) so the stress co-rotates. */
	double cr_e1[3]={0,0,0}, cr_e2[3]={0,0,0};
	bool use_corot = (fCorotational && (size_t)(i+1)*9 <= fFrameR.size());
	if (use_corot) {
		double* Rf=&fFrameR[(size_t)i*9]; double* Vf=&fFrameV[(size_t)i*9];
		double e1o[3]={Rf[0],Rf[3],Rf[6]}, e2o[3]={Rf[1],Rf[4],Rf[7]};   /* old frame tangents */
		double du1[3]={0,0,0}, du2[3]={0,0,0};
		for (int I=0;I<nn;I++){ double p1=Dp[I*5],p2=Dp[I*5+1];
			for(int d=0;d<3;d++){ double du=due[I*3+d]; du1[d]+=p1*du; du2[d]+=p2*du; } }
		double g11=Dot(x1,x1),g22=Dot(x2,x2),g12=Dot(x1,x2), dgm=g11*g22-g12*g12;
		double gu1[3]={0,0,0},gu2[3]={0,0,0};
		if (std::fabs(dgm)>1.0e-30){ double gi11=g22/dgm,gi22=g11/dgm,gi12=-g12/dgm;
			for(int d=0;d<3;d++){ gu1[d]=gi11*x1[d]+gi12*x2[d]; gu2[d]=gi12*x1[d]+gi22*x2[d]; } }
		double dL[3][3];
		for(int a=0;a<3;a++)for(int b=0;b<3;b++) dL[a][b]=du1[a]*gu1[b]+du2[a]*gu2[b];
		double R[3][3],V[3][3];
		for(int a=0;a<3;a++)for(int b=0;b<3;b++){ R[a][b]=Rf[a*3+b]; V[a][b]=Vf[a*3+b]; }
		if (commit) KLShell::FlanaganTaylorStep(dL,R,V);
		double ng[3]; Cross(x1,x2,ng); double nl=Norm(ng);   /* re-align R[:,2] to geometric normal */
		if (nl>1.0e-30){ for(int d=0;d<3;d++) ng[d]/=nl;
			double e1n[3]={R[0][0],R[1][0],R[2][0]};
			double dp=e1n[0]*ng[0]+e1n[1]*ng[1]+e1n[2]*ng[2];
			for(int d=0;d<3;d++) e1n[d]-=dp*ng[d];
			double el=Norm(e1n);
			if (el>1.0e-30){ for(int d=0;d<3;d++) e1n[d]/=el; double e2n[3]; Cross(ng,e1n,e2n);
				for(int d=0;d<3;d++){ R[d][0]=e1n[d]; R[d][1]=e2n[d]; R[d][2]=ng[d]; } }
		}
		cr_e1[0]=R[0][0];cr_e1[1]=R[1][0];cr_e1[2]=R[2][0];
		cr_e2[0]=R[0][1];cr_e2[1]=R[1][1];cr_e2[2]=R[2][1];
		if (commit && plastic && i<(int)fJ2sig.size()){   /* co-rotate stored in-plane stress old->new */
			for (int g=0; 3*g+2 < (int)fJ2sig[i].size(); g++){
				double s11=fJ2sig[i][3*g],s22=fJ2sig[i][3*g+1],s12=fJ2sig[i][3*g+2];
				double sg[3][3];
				for(int a=0;a<3;a++)for(int b=0;b<3;b++)
					sg[a][b]=s11*e1o[a]*e1o[b]+s22*e2o[a]*e2o[b]+s12*(e1o[a]*e2o[b]+e2o[a]*e1o[b]);
				double n11=0,n22=0,n12=0;
				for(int a=0;a<3;a++)for(int b=0;b<3;b++){
					n11+=cr_e1[a]*sg[a][b]*cr_e1[b]; n22+=cr_e2[a]*sg[a][b]*cr_e2[b]; n12+=cr_e1[a]*sg[a][b]*cr_e2[b]; }
				fJ2sig[i][3*g]=n11; fJ2sig[i][3*g+1]=n22; fJ2sig[i][3*g+2]=n12;
			}
		}
		if (commit) for(int a=0;a<3;a++)for(int b=0;b<3;b++){ Rf[a*3+b]=R[a][b]; Vf[a*3+b]=V[a][b]; }
	}

	for (int g=0; g<3; g++){
		ShellGeom G;
		if (!BuildGeom(x1,x2,x11,x22,x12,h,xg[g],G)) continue;
		double e1[3],e2[3];
		if (use_corot){   /* tracked tangents, re-orthogonalized to this station's normal */
			double dp=cr_e1[0]*G.n[0]+cr_e1[1]*G.n[1]+cr_e1[2]*G.n[2];
			for(int d=0;d<3;d++) e1[d]=cr_e1[d]-dp*G.n[d];
			double el=Norm(e1);
			if (el>1.0e-30){ for(int d=0;d<3;d++) e1[d]/=el; Cross(G.n,e1,e2); }
			else OrthoTangents(G.n,e1,e2);
		} else OrthoTangents(G.n,e1,e2);
		std::vector<std::vector<double> > Bv(nn, std::vector<double>(18));
		for (int I=0;I<nn;I++){
			double B[3][3][3];
			BMatrix(G, Dp[I*5], Dp[I*5+1], Dp[I*5+2], Dp[I*5+4], Dp[I*5+3], B);
			double bv[6][3]; ToVoigtLocal(B,e1,e2,G.n,bv);
			for(int r=0;r<6;r++)for(int cc=0;cc<3;cc++) Bv[I][r*3+cc]=bv[r][cc];
		}
		double deps[6];
		for (int r=0;r<6;r++){ double s=0.0;
			for(int I=0;I<nn;I++)for(int cc=0;cc<3;cc++) s+=Bv[I][r*3+cc]*due[I*3+cc]; deps[r]=s; }

		/* accumulated in-plane stress: plane-stress J2 increment, or elastic predictor */
		double sip[3]={fJ2sig[i][3*g],fJ2sig[i][3*g+1],fJ2sig[i][3*g+2]};
		double dip[3]={deps[0],deps[1],deps[5]};
		if (plastic && g < (int)fJ2ep[i].size()){
			double ep=fJ2ep[i][g];
			if (fThicknessUpdate) {   /* Algorithm 3: sigma33=0 via secant on D33, returns the thickness strain */
				double de33 = PlaneStressJ2_D33(sip,dip,ep,fYoung,fPoisson,fYield,fHardening,fYieldSat,fSatRate);
				if (g==1) de33_mid = de33;   /* mid-surface (xi3=0) -> the membrane thickness change */
			} else {
				PlaneStressJ2Return(sip,dip,ep,fYoung,fPoisson,fYield,fHardening,fYieldSat,fSatRate);
			}
			if (commit) fJ2ep[i][g]=ep;
		} else {
			sip[0]+=c_ps*(dip[0]+fPoisson*dip[1]);
			sip[1]+=c_ps*(fPoisson*dip[0]+dip[1]);
			sip[2]+=c_ps*(1.0-fPoisson)/2.0*dip[2];
		}
		if (commit){ fJ2sig[i][3*g]=sip[0]; fJ2sig[i][3*g+1]=sip[1]; fJ2sig[i][3*g+2]=sip[2]; }
		double sig[6]={sip[0],sip[1],0,0,0,sip[2]};

		double w = wg[g]*(h/2.0)*A_K;
		for (int I=0;I<nn;I++)
			for(int cc=0;cc<3;cc++){ double s=0.0; for(int r=0;r<6;r++) s+=Bv[I][r*3+cc]*sig[r]; fout[I*3+cc]+=s*w; }
	}

	/* Algorithm 3 thickness accumulation: t_{n+1} = t_n * exp(D33*dt) (committed steps only) */
	if (commit && fThicknessUpdate && i < (int)fThicknessCur.size()) {
		double tn = fThicknessCur[i]*std::exp(de33_mid);
		if (tn > 1.0e-6*fThickness) fThicknessCur[i] = tn;   /* guard against collapse */
	}

	/* SCNI stabilization (reference-config, linear) from the stored stabilization points */
	int npt = (int) fIPw[i].size();
	if (npt > 0){ const double* Ball=&fIPB[i][0];
		for (int p=0;p<npt;p++){ if (fIPstab[i][p]==0) continue;
			const double* B=Ball+(size_t)p*6*ndof;
			double eps[6]; for(int r=0;r<6;r++){const double*Br=B+(size_t)r*ndof;double s=0.0;for(int c=0;c<ndof;c++)s+=Br[c]*ue[c];eps[r]=s;}
			double sig[6]; for(int r=0;r<6;r++){double s=0.0;for(int cc=0;cc<6;cc++)s+=fC[r][cc]*eps[cc];sig[r]=s;}
			double w=fIPw[i][p];
			for(int c=0;c<ndof;c++){double s=0.0;for(int r=0;r<6;r++)s+=B[(size_t)r*ndof+c]*sig[r];fout[c]+=s*w;}
		}
	}

	/* bending-hourglass control (rank-1 penalty along the CURRENT-config normal so the out-of-plane
	 * penalty stays orthogonal to the deformed tangent plane; a static reference normal would inject
	 * a spurious in-plane/membrane component at large crush rotations and over-stiffen the hinges) */
	if (include_bend && fBendCoeff[i] != 0.0 && (int)fBendR[i].size()==nn) {
		const double* Rb=&fBendR[i][0];
		double nv[3]; Cross(x1,x2,nv); double nL=Norm(nv);
		if (nL>1.0e-300){ for(int d=0;d<3;d++) nv[d]/=nL;
			double dp=nv[0]*fBendN[i][0]+nv[1]*fBendN[i][1]+nv[2]*fBendN[i][2];
			if (dp<0.0) for(int d=0;d<3;d++) nv[d]=-nv[d]; }   /* keep orientation consistent with ref */
		else { for(int d=0;d<3;d++) nv[d]=fBendN[i][d]; }
		double kru=0.0;
		for (int I=0;I<nn;I++){ double nu=nv[0]*ue[I*3]+nv[1]*ue[I*3+1]+nv[2]*ue[I*3+2]; kru+=Rb[I]*nu; }
		double c=fBendCoeff[i]*kru;
		for (int I=0;I<nn;I++){ double cr=c*Rb[I]; for(int d=0;d<3;d++) fout[I*3+d]+=cr*nv[d]; }
	}
}

/* LHS: stiffness (implicit) and/or lumped mass (explicit). Uses the inherited Top/NextElement +
 * no-arg AssembleLHS() flow (CurrentElement().Equations()), as the standard meshfree element. */
void RKShellT::LHSDriver(GlobalT::SystemTypeT sys_type)
{
#pragma unused(sys_type)
	double constM = 0.0, constK = 0.0;
	int formM = fIntegrator->FormM(constM);
	int formK = fIntegrator->FormK(constK);
	if (!formM && !formK) return;

	Top();
	while (NextElement()) {
		int i = fElementCards.Position();
		int nn = fNeighbors.MinorDim(i);
		if (nn < 6) continue;

		if (formM) { /* lumped mass on node i's own dof, within its stencil block */
			const int* gnb = fNeighbors(i);
			int ki = -1;
			for (int k = 0; k < nn; k++) if (gnb[k] == fGlobalIDs[i]) { ki = k; break; }
			fLHS.Dimension(3*nn);
			fLHS.SetFormat(ElementMatrixT::kSymmetric);
			fLHS = 0.0;
			if (ki >= 0) for (int d = 0; d < 3; d++) fLHS(3*ki+d, 3*ki+d) = constM*fLumpedMass[i];
			AssembleLHS();
		}
		if (formK) {
			fLHS.Dimension(3*nn);
			fLHS.SetFormat(ElementMatrixT::kSymmetric); /* Ke = B^T C B is symmetric */
			dMatrixT& m = fLHS; m = fKe[i];
			if (constK != 1.0) m *= constK;
			AssembleLHS();
		}
	}
}

/* residual: external load minus internal force (R = f_ext - f_int), per stencil via AssembleRHS() */
void RKShellT::RHSDriver(void)
{
	double constKd = 0.0;
	int formKd = fIntegrator->FormKd(constKd);
	if (!formKd) return;
	const dArray2DT& disp = Field()[0];               /* current nodal displacement */

	/* rate-form finite strain: ensure the previous-step displacement buffer is sized */
	if (fFiniteStrain && (fUprev.MajorDim()!=disp.MajorDim() || fUprev.MinorDim()!=disp.MinorDim())) {
		fUprev.Dimension(disp.MajorDim(), disp.MinorDim()); fUprev = 0.0;
	}

	Top();
	while (NextElement()) {
		int i = fElementCards.Position();
		int nn = fNeighbors.MinorDim(i);
		if (nn < 6) continue;
		const int* gnb = fNeighbors(i);

		/* stress-driven internal force f_int = sum_pt B^T sigma(B*u) w (== fKe*u for elasticity) */
		dArrayT ue(3*nn);
		for (int k=0;k<nn;k++) for (int d=0;d<3;d++) ue[k*3+d] = disp(gnb[k], d);
		if (fFiniteStrain) InternalForceFS(i, ue, fRHS, true);  /* finite-deformation (Fig 18) */
		else               InternalForce(i, ue, fRHS, true);    /* small-strain (commit J2 once/step) */
		fRHS *= -constKd;                              /* residual gets -f_int */

		/* external per-area load on node i (its own dof within the stencil) */
		int ki = -1;
		for (int k = 0; k < nn; k++) if (gnb[k] == fGlobalIDs[i]) { ki = k; break; }
		if (ki >= 0) {
			double A_K = fNodalArea[i];
			for (int d=0;d<3;d++) fRHS[3*ki+d] += constKd*fLoad[d]*A_K;
			/* mass-proportional damping (dynamic relaxation -> quasi-static): -alpha * m_i * v_i */
			if (fDamping > 0.0 && Field().Order() >= 1) {
				const dArray2DT& vel = Field()[1];
				int gi = fGlobalIDs[i];
				for (int d=0;d<3;d++) fRHS[3*ki+d] -= constKd*fDamping*fLumpedMass[i]*vel(gi, d);
			}
		}
		AssembleRHS();
	}

	/* advance the rate-form reference: u_prev <- u (once per step, after all stencils) */
	if (fFiniteStrain) fUprev = disp;
}
