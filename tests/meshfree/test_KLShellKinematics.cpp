/* test_KLShellKinematics.cpp — meshfree KL-shell kinematics & strain-displacement tests.
 *
 * Validates the auxiliary-tensor kinematic chain and the per-node strain-displacement
 * matrix B_Iijk (development/src/elements/meshfree_kl_shell, KLShellKernels.h) by driving
 * analytic deformation modes through the B-matrix and checking the resulting spatial
 * velocity gradient / rate-of-deformation:
 *
 *   - rigid rotation about the surface normal  -> D = 0 (objectivity);
 *   - in-plane stretch eps along psi1          -> psi1 . D . psi1 = eps;
 *   - cylindrical bending w = (1/2) c xi1^2     -> psi1 . D . psi1 = -(h/2) xi3 c
 *     (linear & antisymmetric through-thickness: the KL bending signature);
 *   - parametric gradient B_Iijkl (xi3=0) matches a finite difference of grad(v3D) on a
 *     curved (cylindrical) surface (exercises the curvature-derivative terms).
 */

#include "gtest/gtest.h"

#include <cmath>
#include <vector>

#include "D2OrthoMLS2DT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "KLShellKernels.h"

using namespace Tahoe;
using namespace Tahoe::KLShell;

namespace {

const double kThickness = 0.05;

/* a node cloud + its center node, PCA frame, neighbor list and RK derivatives */
struct Patch {
	std::vector<double> X, Y, Z;        /* node coordinates */
	std::vector<int> nb;                /* neighbors of the center node */
	double psi1[3], psi2[3], n0[3];     /* PCA tangent frame at the center */
	std::vector<double> phi, d1, d2, dd1, dd2, dd12; /* RK values + derivs at center */
};

/* build a flat n x n patch on [0,L]^2 and the RK data at the center node */
bool BuildFlatPatch(int n, double L, Patch& p)
{
	double h = L / (n - 1);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) {
			p.X.push_back(i * h);
			p.Y.push_back(j * h);
			p.Z.push_back(0.0);
		}
	int N = int(p.X.size());

	int P = -1;
	double best = 1e30;
	for (int i = 0; i < N; i++) {
		double d = std::fabs(p.X[i] - 0.5 * L) + std::fabs(p.Y[i] - 0.5 * L);
		if (d < best) { best = d; P = i; }
	}

	double R = 3.5 * h;
	for (int Q = 0; Q < N; Q++) {
		double dx = p.X[Q] - p.X[P], dy = p.Y[Q] - p.Y[P];
		if (std::sqrt(dx * dx + dy * dy) < 0.99 * R) p.nb.push_back(Q);
	}
	int nn = int(p.nb.size());

	std::vector<double> nbX(3 * nn);
	for (int k = 0; k < nn; k++) {
		nbX[3 * k] = p.X[p.nb[k]];
		nbX[3 * k + 1] = p.Y[p.nb[k]];
		nbX[3 * k + 2] = 0.0;
	}
	pcaFrame(&nbX[0], nn, p.psi1, p.psi2, p.n0);
	if (p.n0[2] < 0.0)
		for (int d = 0; d < 3; d++) { p.n0[d] = -p.n0[d]; p.psi2[d] = -p.psi2[d]; }

	dArray2DT lc(nn, 2);
	for (int k = 0; k < nn; k++) {
		double dxv[3] = {p.X[p.nb[k]] - p.X[P], p.Y[p.nb[k]] - p.Y[P], 0.0};
		lc(k, 0) = dot3(dxv, p.psi1);
		lc(k, 1) = dot3(dxv, p.psi2);
	}

	D2OrthoMLS2DT efg(2);
	efg.Initialize();
	dArrayT dmax(nn);
	dmax = R;
	dArrayT sample(2);
	sample[0] = 0.0;
	sample[1] = 0.0;
	if (!efg.SetField(lc, dmax, sample)) return false;

	const dArrayT& phi = efg.phi();
	const dArray2DT& Dp = efg.Dphi();
	const dArray2DT& DDp = efg.DDphi();
	p.phi.resize(nn); p.d1.resize(nn); p.d2.resize(nn);
	p.dd1.resize(nn); p.dd2.resize(nn); p.dd12.resize(nn);
	for (int I = 0; I < nn; I++) {
		p.phi[I] = phi[I];
		p.d1[I] = Dp(0, I); p.d2[I] = Dp(1, I);
		p.dd1[I] = DDp(0, I); p.dd2[I] = DDp(1, I); p.dd12[I] = DDp(2, I);
	}
	return true;
}

/* reference position parametric derivatives at the center node */
void RefDerivs(const Patch& p, double x1[3], double x2[3], double x11[3], double x22[3], double x12[3])
{
	for (int d = 0; d < 3; d++) { x1[d]=x2[d]=x11[d]=x22[d]=x12[d]=0.0; }
	for (int I = 0; I < int(p.nb.size()); I++) {
		double Xq[3] = {p.X[p.nb[I]], p.Y[p.nb[I]], p.Z[p.nb[I]]};
		for (int d = 0; d < 3; d++) {
			x1[d] += p.d1[I] * Xq[d];   x2[d] += p.d2[I] * Xq[d];
			x11[d] += p.dd1[I] * Xq[d]; x22[d] += p.dd2[I] * Xq[d]; x12[d] += p.dd12[I] * Xq[d];
		}
	}
}

/* spatial velocity gradient grad(v3D) = sum_I B_Iijk v_Ik at the center node, station xi3,
 * for a velocity field sampled at the neighbor nodes (nodalV: 3*nn). */
void VelocityGradient(const Patch& p, double xi3, const std::vector<double>& nodalV, double L[3][3])
{
	double x1[3], x2[3], x11[3], x22[3], x12[3];
	RefDerivs(p, x1, x2, x11, x22, x12);
	ShellGeom g;
	buildGeom(x1, x2, x11, x22, x12, kThickness, xi3, g);

	for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) L[i][j] = 0.0;
	for (int I = 0; I < int(p.nb.size()); I++) {
		double B[3][3][3];
		Bmatrix(g, p.d1[I], p.d2[I], p.dd1[I], p.dd12[I], p.dd2[I], B);
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				for (int k = 0; k < 3; k++)
					L[i][j] += B[i][j][k] * nodalV[3 * I + k];
	}
}

double Project(const double M[3][3], const double a[3], const double b[3])
{
	double s = 0.0;
	for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) s += a[i] * M[i][j] * b[j];
	return s;
}

} /* anonymous namespace */

TEST(KLShellKinematics, RigidRotationGivesZeroStrain)
{
	Patch p;
	ASSERT_TRUE(BuildFlatPatch(15, 1.0, p));
	int nn = int(p.nb.size());

	/* v(Q) = omega * (n x (X_Q - X_P)) */
	double XP[3] = {p.X[p.nb[0]], p.Y[p.nb[0]], 0.0};
	for (int i = 0; i < nn; i++) if (std::fabs(p.X[p.nb[i]] - 0.5) < 1e-9 && std::fabs(p.Y[p.nb[i]] - 0.5) < 1e-9)
		{ XP[0] = p.X[p.nb[i]]; XP[1] = p.Y[p.nb[i]]; }
	std::vector<double> v(3 * nn);
	for (int I = 0; I < nn; I++) {
		double r[3] = {p.X[p.nb[I]] - XP[0], p.Y[p.nb[I]] - XP[1], 0.0};
		double vv[3];
		cross3(p.n0, r, vv);
		for (int d = 0; d < 3; d++) v[3 * I + d] = 0.7 * vv[d];
	}

	double Lm[3][3];
	VelocityGradient(p, 0.0, v, Lm);
	double D = 0.0;
	for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) D = std::max(D, std::fabs(0.5 * (Lm[i][j] + Lm[j][i])));
	EXPECT_NEAR(D, 0.0, 1e-9);
}

TEST(KLShellKinematics, InPlaneStretch)
{
	Patch p;
	ASSERT_TRUE(BuildFlatPatch(15, 1.0, p));
	int nn = int(p.nb.size());
	double XP[3] = {0.5, 0.5, 0.0};
	const double eps = 0.3;

	std::vector<double> v(3 * nn);
	for (int I = 0; I < nn; I++) {
		double r[3] = {p.X[p.nb[I]] - XP[0], p.Y[p.nb[I]] - XP[1], 0.0};
		double s = eps * dot3(r, p.psi1);
		for (int d = 0; d < 3; d++) v[3 * I + d] = s * p.psi1[d];
	}

	double Lm[3][3], D[3][3];
	VelocityGradient(p, 0.0, v, Lm);
	for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) D[i][j] = 0.5 * (Lm[i][j] + Lm[j][i]);
	EXPECT_NEAR(Project(D, p.psi1, p.psi1), eps, 1e-6);
	EXPECT_NEAR(Project(D, p.psi2, p.psi2), 0.0, 1e-6);
}

TEST(KLShellKinematics, ThroughThicknessBending)
{
	Patch p;
	ASSERT_TRUE(BuildFlatPatch(15, 1.0, p));
	int nn = int(p.nb.size());
	double XP[3] = {0.5, 0.5, 0.0};
	const double c = 1.0; /* curvature parameter */

	/* transverse field w = 1/2 c xi1^2 along the normal */
	std::vector<double> v(3 * nn);
	for (int I = 0; I < nn; I++) {
		double r[3] = {p.X[p.nb[I]] - XP[0], p.Y[p.nb[I]] - XP[1], 0.0};
		double xi1 = dot3(r, p.psi1);
		double w = 0.5 * c * xi1 * xi1;
		for (int d = 0; d < 3; d++) v[3 * I + d] = w * p.n0[d];
	}

	double Lp[3][3], Lz[3][3], Lm[3][3], Dp[3][3], Dz[3][3], Dm[3][3];
	VelocityGradient(p, +1.0, v, Lp);
	VelocityGradient(p,  0.0, v, Lz);
	VelocityGradient(p, -1.0, v, Lm);
	for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
		Dp[i][j] = 0.5 * (Lp[i][j] + Lp[j][i]);
		Dz[i][j] = 0.5 * (Lz[i][j] + Lz[j][i]);
		Dm[i][j] = 0.5 * (Lm[i][j] + Lm[j][i]);
	}
	double sp = Project(Dp, p.psi1, p.psi1);
	double sz = Project(Dz, p.psi1, p.psi1);
	double sm = Project(Dm, p.psi1, p.psi1);

	EXPECT_NEAR(sp, -(kThickness / 2.0) * (+1.0) * c, 1e-4); /* KL bending strain */
	EXPECT_NEAR(sm, -(kThickness / 2.0) * (-1.0) * c, 1e-4);
	EXPECT_NEAR(sz, 0.0, 1e-4);                              /* zero membrane part */
	EXPECT_NEAR(sp + sm, 0.0, 1e-6);                         /* antisymmetric */
}

/* B_Iijkl (xi3=0) vs finite difference of grad(v3D) on a curved (cylinder) surface */
TEST(KLShellKinematics, ParametricGradientOnCylinder)
{
	const double Rc = 1.0, Zlen = 2.0;
	const int nt = 64, nz = 21;
	double dth = 2.0 * M_PI / nt, dz = Zlen / (nz - 1);
	double spacing = (Rc * dth > dz) ? Rc * dth : dz;
	double R = 3.5 * spacing;

	std::vector<double> X, Y, Z;
	for (int i = 0; i < nt; i++)
		for (int j = 0; j < nz; j++) {
			double th = i * dth;
			X.push_back(Rc * std::cos(th));
			Y.push_back(Rc * std::sin(th));
			Z.push_back(j * dz);
		}
	int N = int(X.size());

	int P = -1;
	double best = 1e30;
	for (int i = 0; i < N; i++) {
		double d = std::fabs(Z[i] - Zlen / 2) + std::fabs(X[i] - Rc);
		if (d < best) { best = d; P = i; }
	}
	std::vector<int> nb;
	for (int Q = 0; Q < N; Q++) {
		double d0 = X[Q] - X[P], d1 = Y[Q] - Y[P], d2 = Z[Q] - Z[P];
		if (std::sqrt(d0 * d0 + d1 * d1 + d2 * d2) < 0.99 * R) nb.push_back(Q);
	}
	int nn = int(nb.size());
	ASSERT_GE(nn, 6);

	std::vector<double> nbX(3 * nn);
	for (int k = 0; k < nn; k++) { nbX[3*k]=X[nb[k]]; nbX[3*k+1]=Y[nb[k]]; nbX[3*k+2]=Z[nb[k]]; }
	double psi1[3], psi2[3], n0[3];
	pcaFrame(&nbX[0], nn, psi1, psi2, n0);

	dArray2DT lc(nn, 2);
	for (int k = 0; k < nn; k++) {
		double dxv[3] = {X[nb[k]] - X[P], Y[nb[k]] - Y[P], Z[nb[k]] - Z[P]};
		lc(k, 0) = dot3(dxv, psi1);
		lc(k, 1) = dot3(dxv, psi2);
	}

	/* arbitrary smooth velocity field at the nodes */
	std::vector<double> V(3 * nn);
	for (int k = 0; k < nn; k++) {
		V[3*k]   = 0.3 * Y[nb[k]] * Z[nb[k]];
		V[3*k+1] = 0.2 * X[nb[k]];
		V[3*k+2] = 0.1 * X[nb[k]] * X[nb[k]];
	}

	D2OrthoMLS2DT efg(2);
	efg.Initialize();
	dArrayT dmax(nn);
	dmax = R;

	/* reconstruct geometry + velocity derivs + grad(v3D) at a shifted sample (xi3=0) */
	struct Recon {
		D2OrthoMLS2DT* efg; const dArray2DT* lc; const dArrayT* dmax;
		const std::vector<double>* X; const std::vector<double>* Y; const std::vector<double>* Z;
		const std::vector<int>* nb; const std::vector<double>* V;
		void gradV(double s0, double s1, double L[3][3]) {
			dArrayT sample(2); sample[0]=s0; sample[1]=s1;
			efg->SetField(*lc, *dmax, sample);
			const dArray2DT& Dp = efg->Dphi();
			int nn = lc->MajorDim();
			double x1[3]={0,0,0},x2[3]={0,0,0},x11[3]={0,0,0},x22[3]={0,0,0},x12[3]={0,0,0};
			double v1[3]={0,0,0},v2[3]={0,0,0};
			const dArray2DT& DDp = efg->DDphi();
			for (int I=0;I<nn;I++){
				double Xq[3]={(*X)[(*nb)[I]],(*Y)[(*nb)[I]],(*Z)[(*nb)[I]]};
				double Vq[3]={(*V)[3*I],(*V)[3*I+1],(*V)[3*I+2]};
				for(int d=0;d<3;d++){ x1[d]+=Dp(0,I)*Xq[d];x2[d]+=Dp(1,I)*Xq[d];
					x11[d]+=DDp(0,I)*Xq[d];x22[d]+=DDp(1,I)*Xq[d];x12[d]+=DDp(2,I)*Xq[d];
					v1[d]+=Dp(0,I)*Vq[d];v2[d]+=Dp(1,I)*Vq[d]; }
			}
			ShellGeom g; buildGeom(x1,x2,x11,x22,x12,kThickness,0.0,g);
			/* grad(v3D)_ij = sum_I B_Iijk v_Ik (rebuild B per neighbor at this sample) */
			for(int i=0;i<3;i++)for(int j=0;j<3;j++)L[i][j]=0.0;
			for(int I=0;I<nn;I++){ double B[3][3][3];
				Bmatrix(g,Dp(0,I),Dp(1,I),DDp(0,I),DDp(2,I),DDp(1,I),B);
				double Vq[3]={(*V)[3*I],(*V)[3*I+1],(*V)[3*I+2]};
				for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int k=0;k<3;k++)L[i][j]+=B[i][j][k]*Vq[k]; }
		}
	} rec;
	rec.efg=&efg; rec.lc=&lc; rec.dmax=&dmax; rec.X=&X; rec.Y=&Y; rec.Z=&Z; rec.nb=&nb; rec.V=&V;

	/* analytic sum_I B_Iijkl v_Ik via BmatrixGrad at xi3=0 */
	dArrayT sample(2); sample[0]=0; sample[1]=0;
	efg.SetField(lc, dmax, sample);
	const dArray2DT& Dp = efg.Dphi();
	const dArray2DT& DDp = efg.DDphi();
	double x1[3]={0,0,0},x2[3]={0,0,0},x11[3]={0,0,0},x22[3]={0,0,0},x12[3]={0,0,0};
	for(int I=0;I<nn;I++){ double Xq[3]={X[nb[I]],Y[nb[I]],Z[nb[I]]};
		for(int d=0;d<3;d++){x1[d]+=Dp(0,I)*Xq[d];x2[d]+=Dp(1,I)*Xq[d];x11[d]+=DDp(0,I)*Xq[d];x22[d]+=DDp(1,I)*Xq[d];x12[d]+=DDp(2,I)*Xq[d];}}
	ShellGeom g0; buildGeom(x1,x2,x11,x22,x12,kThickness,0.0,g0);
	double Bg_sum[3][3][2]={{{0,0},{0,0},{0,0}},{{0,0},{0,0},{0,0}},{{0,0},{0,0},{0,0}}};
	for(int I=0;I<nn;I++){
		double Bz[3][3][3]; Bmatrix(g0,Dp(0,I),Dp(1,I),DDp(0,I),DDp(2,I),DDp(1,I),Bz);
		double P1l[2]={DDp(0,I),DDp(2,I)}, P2l[2]={DDp(2,I),DDp(1,I)};
		double Bg[3][3][3][2]; BmatrixGrad(g0,Dp(0,I),Dp(1,I),DDp(0,I),DDp(2,I),DDp(1,I),P1l,P2l,Bz,Bg);
		double Vq[3]={V[3*I],V[3*I+1],V[3*I+2]};
		for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int l=0;l<2;l++)for(int k=0;k<3;k++)Bg_sum[i][j][l]+=Bg[i][j][k][l]*Vq[k];
	}

	/* finite difference of grad(v3D) in xi1, xi2 */
	double delta = 1e-4;
	double maxErr = 0.0, scale = 0.0;
	for (int l = 0; l < 2; l++) {
		double sp[2] = {0, 0}, sm[2] = {0, 0};
		sp[l] = delta; sm[l] = -delta;
		double Lp[3][3], Lmn[3][3];
		rec.gradV(sp[0], sp[1], Lp);
		rec.gradV(sm[0], sm[1], Lmn);
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++) {
				double fd = (Lp[i][j] - Lmn[i][j]) / (2.0 * delta);
				maxErr = std::max(maxErr, std::fabs(Bg_sum[i][j][l] - fd));
				scale = std::max(scale, std::fabs(fd));
			}
	}
	EXPECT_LT(maxErr / (scale + 1e-12), 1e-5);
}
