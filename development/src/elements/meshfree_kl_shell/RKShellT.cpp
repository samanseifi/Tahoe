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
}

/* meshfree element: no FE-block output (override ElementBaseT's block-based output).
 * WriteOutput prints a deflection summary (temporary validation diagnostic). */
void RKShellT::RegisterOutput(void) {}
void RKShellT::WriteOutput(void)
{
	const dArray2DT& disp = Field()[0];
	double minuy = 0.0, maxmag = 0.0;
	for (int i = 0; i < fNumNodes; i++) {
		int g = fGlobalIDs[i];
		double uy = disp(g,1);
		if (uy < minuy) minuy = uy;
		double mag = std::sqrt(disp(g,0)*disp(g,0)+disp(g,1)*disp(g,1)+disp(g,2)*disp(g,2));
		if (mag > maxmag) maxmag = mag;
	}
	fprintf(stdout, "[RKShell] nodes=%d  min(u_y)=% .6e  max|u|=% .6e\n", fNumNodes, minuy, maxmag);
	fflush(stdout);
}

void RKShellT::ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const
{
#pragma unused(connects_1)
	connects_2.Append(&fNeighbors);
}

void RKShellT::ConnectsX(AutoArrayT<const iArray2DT*>& connects) const
{
#pragma unused(connects)
	/* no fixed-size geometry connectivity for a meshfree element */
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
			double Mw = V_K*Mmom;
			double Mw_bend = (A_K*Mmom)*(h*h*h/12.0);
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

/* tangent stiffness: assemble each per-node stencil stiffness */
void RKShellT::LHSDriver(GlobalT::SystemTypeT sys_type)
{
#pragma unused(sys_type)
	double constK = 0.0;
	int formK = fIntegrator->FormK(constK);
	if (!formK) return;
	int group = Group();
	for (int i = 0; i < fNumNodes; i++) {
		int nn = fNeighbors.MinorDim(i);
		if (nn < 6) continue;
		fLHS.Dimension(3*nn);
		fLHS.SetFormat(ElementMatrixT::kSymmetric); /* Ke = B^T C B is symmetric */
		dMatrixT& m = fLHS;
		m = fKe[i];
		if (constK != 1.0) m *= constK;
		iArrayT eq; eq.Alias(3*nn, fEqnos(i));
		ElementSupport().AssembleLHS(group, fLHS, eq);
	}
}

/* residual: external load minus internal force (R = f_ext - f_int) */
void RKShellT::RHSDriver(void)
{
	double constKd = 0.0;
	int formKd = fIntegrator->FormKd(constKd);
	if (!formKd) return;
	int group = Group();
	const dArray2DT& disp = Field()[0];               /* current nodal displacement */
	const iArray2DT& field_eqnos = Field().Equations();

	/* internal force: per stencil, -constKd * Ke * u_stencil */
	for (int i = 0; i < fNumNodes; i++) {
		int nn = fNeighbors.MinorDim(i);
		if (nn < 6) continue;
		const int* gnb = fNeighbors(i);
		dArrayT ue(3*nn);
		for (int k=0;k<nn;k++) for (int d=0;d<3;d++) ue[k*3+d] = disp(gnb[k], d);
		dArrayT fint(3*nn);
		fKe[i].Multx(ue, fint);
		fint *= -constKd;
		iArrayT eq; eq.Alias(3*nn, fEqnos(i));
		ElementSupport().AssembleRHS(group, fint, eq);
	}

	/* external uniform per-area load: + constKd * load * A_i at node i */
	for (int i = 0; i < fNumNodes; i++) {
		double A_K = fNodalArea[i];
		dArrayT fext(3);
		for (int d=0;d<3;d++) fext[d] = constKd*fLoad[d]*A_K;
		iArrayT eq; eq.Alias(3, field_eqnos(fGlobalIDs[i]));
		ElementSupport().AssembleRHS(group, fext, eq);
	}
}
