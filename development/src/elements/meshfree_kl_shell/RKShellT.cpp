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
#include "ElementMatrixT.h"
#include "eIntegratorT.h"
#include "OutputSetT.h"
#include "GeometryT.h"
#include "iArray2DT.h"

#include <cmath>
#include <vector>

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
	fStabMode = 0;
	fStabMembrane = 1.0;
	fStabBending = 1.0;
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

	ArrayT<StringT> n_labels(3);
	n_labels[0] = "D_X"; n_labels[1] = "D_Y"; n_labels[2] = "D_Z";

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
	fflush(stdout);

	/* write the displacement field for visualization */
	if (fOutputID < 0) return;
	dArray2DT n_values(fOutputNodesUsed.Length(), 3);
	for (int k = 0; k < fOutputNodesUsed.Length(); k++)
		for (int d = 0; d < 3; d++) n_values(k,d) = disp(fOutputNodesUsed[k], d);
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
	ParameterT sb(fStabBending,  "stab_bending");  sb.SetDefault(1.0); list.AddParameter(sb);
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
			for (int I=0;I<nn;I++) for (int J=0;J<nn;J++){
				double kij[3][3]={{0,0,0},{0,0,0},{0,0,0}};
				for(int a=0;a<6;a++)for(int b=0;b<6;b++){double Cab=fC[a][b];if(Cab==0.0)continue;
					for(int ci=0;ci<3;ci++)for(int cj=0;cj<3;cj++) kij[ci][cj]+=Bv[I][a*3+ci]*Cab*Bv[J][b*3+cj];}
				for(int ci=0;ci<3;ci++)for(int cj=0;cj<3;cj++) Ke(3*I+ci,3*J+cj)+=cw*kij[ci][cj];
			}
		}

		/* membrane + bending (curvature-gradient) stabilization at xi3=0 */
		ShellGeom G0;
		if (BuildGeom(x1,x2,x11,x22,x12,h,0.0,G0)) {
			double e1[3],e2[3]; OrthoTangents(G0.n,e1,e2);
			std::vector<std::vector<double> > Bg(nn, std::vector<double>(36)), Bk(nn, std::vector<double>(36));
			for (int I=0;I<nn;I++){
				double Bz[3][3][3];
				BMatrix(G0,Dp(0,I),Dp(1,I),DDp(0,I),DDp(2,I),DDp(1,I),Bz);
				double P1l[2]={DDp(0,I),DDp(2,I)}, P2l[2]={DDp(2,I),DDp(1,I)};
				double Bgr[3][3][3][2];
				BMatrixGradient(G0,Dp(0,I),Dp(1,I),P1l,P2l,Bz,Bgr);
				double Bkk[3][3][3][2];
				BMatrixCurvatureGradient(G0,DDp(0,I),DDp(2,I),DDp(1,I),DDDp(0,I),DDDp(2,I),DDDp(1,I),DDDp(3,I),Bkk);
				for (int l=0;l<2;l++){
					double Bm[3][3][3], Bc[3][3][3];
					for(int ii=0;ii<3;ii++)for(int jj=0;jj<3;jj++)for(int kk=0;kk<3;kk++){
						Bm[ii][jj][kk]=Bgr[ii][jj][kk][l]; Bc[ii][jj][kk]=Bkk[ii][jj][kk][l]; }
					double bvm[6][3], bvc[6][3];
					ToVoigtLocal(Bm,e1,e2,G0.n,bvm); ToVoigtLocal(Bc,e1,e2,G0.n,bvc);
					for(int r=0;r<6;r++)for(int c=0;c<3;c++){Bg[I][l*18+r*3+c]=bvm[r][c]; Bk[I][l*18+r*3+c]=bvc[r][c];}
				}
			}
			double Mw, Mw_bend;
			if (fStabMode == 3) {        /* section 5.3: alpha-scaled membrane, no bending */
				Mw = fStabMembrane*alpha*V_K*Mmom;
				Mw_bend = 0.0;
			} else if (fStabMode == 2) { /* section 5.2: pure bending, no membrane */
				Mw = 0.0;
				Mw_bend = fStabBending*(A_K*Mmom)*(h*h*h/12.0);
			} else {                     /* default = paper Eq 33: MEMBRANE stab only (bending energy
			                              * comes from the through-thickness Gauss base; the paper
			                              * avoids 3rd-deriv bending stab) */
				Mw = fStabMembrane*V_K*Mmom;
				Mw_bend = 0.0;
			}
			for (int l=0;l<2;l++) for (int I=0;I<nn;I++) for (int J=0;J<nn;J++){
				double km[3][3]={{0,0,0},{0,0,0},{0,0,0}}, kb[3][3]={{0,0,0},{0,0,0},{0,0,0}};
				for(int a=0;a<6;a++)for(int b=0;b<6;b++){double Cab=fC[a][b];if(Cab==0.0)continue;
					for(int ci=0;ci<3;ci++)for(int cj=0;cj<3;cj++){
						km[ci][cj]+=Bg[I][l*18+a*3+ci]*Cab*Bg[J][l*18+b*3+cj];
						kb[ci][cj]+=Bk[I][l*18+a*3+ci]*Cab*Bk[J][l*18+b*3+cj]; }}
				for(int ci=0;ci<3;ci++)for(int cj=0;cj<3;cj++) Ke(3*I+ci,3*J+cj)+=Mw*km[ci][cj]+Mw_bend*kb[ci][cj];
			}
		}
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

	Top();
	while (NextElement()) {
		int i = fElementCards.Position();
		int nn = fNeighbors.MinorDim(i);
		if (nn < 6) continue;
		const int* gnb = fNeighbors(i);

		/* internal force f_int = Ke * u_stencil */
		dArrayT ue(3*nn);
		for (int k=0;k<nn;k++) for (int d=0;d<3;d++) ue[k*3+d] = disp(gnb[k], d);
		fRHS.Dimension(3*nn);
		fKe[i].Multx(ue, fRHS);
		fRHS *= -constKd;                              /* residual gets -f_int */

		/* external per-area load on node i (its own dof within the stencil) */
		int ki = -1;
		for (int k = 0; k < nn; k++) if (gnb[k] == fGlobalIDs[i]) { ki = k; break; }
		if (ki >= 0) {
			double A_K = fNodalArea[i];
			for (int d=0;d<3;d++) fRHS[3*ki+d] += constKd*fLoad[d]*A_K;
		}
		AssembleRHS();
	}
}
