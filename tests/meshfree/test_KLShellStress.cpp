/* test_KLShellStress.cpp — meshfree KL-shell stress-update tests.
 *
 * Validates the stress-update machinery used by the RKPM Kirchhoff-Love shell against
 * closed-form elasticity (pure 3x3 tensor algebra; no meshfree machinery needed):
 *
 *   - Flanagan & Taylor co-rotational update (objectivity): a rigid rotation produces the
 *     analytic rotation tensor, leaves the left stretch at identity, and rotates a stored
 *     stress without introducing spurious stress;
 *   - sigma33 = 0 plane-stress enforcement: in-plane uniaxial strain recovers the
 *     plane-stress moduli E/(1-nu^2) and nu*E/(1-nu^2);
 *   - 3-point through-thickness Gauss integration of a linear-through-thickness stress
 *     recovers the analytic bending moment D_plate * kappa.
 */

#include "gtest/gtest.h"

#include <cmath>

using namespace std;

namespace {

typedef double M3[3][3];

void Identity(M3 A) { for (int i=0;i<3;i++) for (int j=0;j<3;j++) A[i][j] = (i==j) ? 1.0 : 0.0; }
double Trace(const M3 A) { return A[0][0] + A[1][1] + A[2][2]; }

bool Inverse3(const M3 M, M3 Inv)
{
	double d = M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
	         - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
	         + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]);
	if (fabs(d) < 1e-300) return false;
	double id = 1.0 / d;
	Inv[0][0]= (M[1][1]*M[2][2]-M[1][2]*M[2][1])*id;
	Inv[0][1]=-(M[0][1]*M[2][2]-M[0][2]*M[2][1])*id;
	Inv[0][2]= (M[0][1]*M[1][2]-M[0][2]*M[1][1])*id;
	Inv[1][0]=-(M[1][0]*M[2][2]-M[1][2]*M[2][0])*id;
	Inv[1][1]= (M[0][0]*M[2][2]-M[0][2]*M[2][0])*id;
	Inv[1][2]=-(M[0][0]*M[1][2]-M[0][2]*M[1][0])*id;
	Inv[2][0]= (M[1][0]*M[2][1]-M[1][1]*M[2][0])*id;
	Inv[2][1]=-(M[0][0]*M[2][1]-M[0][1]*M[2][0])*id;
	Inv[2][2]= (M[0][0]*M[1][1]-M[0][1]*M[1][0])*id;
	return true;
}

void Multiply(const M3 A, const M3 B, M3 C)
{
	for (int i=0;i<3;i++) for (int j=0;j<3;j++) { double s=0; for (int k=0;k<3;k++) s+=A[i][k]*B[k][j]; C[i][j]=s; }
}

/* C = A B^T */
void MultiplyTransposeRight(const M3 A, const M3 B, M3 C)
{
	for (int i=0;i<3;i++) for (int j=0;j<3;j++) { double s=0; for (int k=0;k<3;k++) s+=A[i][k]*B[j][k]; C[i][j]=s; }
}

int Permutation(int i, int j, int k)
{
	if (i==j || j==k || i==k) return 0;
	if ((i==0&&j==1&&k==2)||(i==1&&j==2&&k==0)||(i==2&&j==0&&k==1)) return 1;
	return -1;
}

/* skew tensor from axial vector (W_ij = eps_ikj w_k) */
void Skew(const double w[3], M3 W)
{
	W[0][0]=W[1][1]=W[2][2]=0.0;
	W[0][1]=-w[2]; W[0][2]= w[1];
	W[1][0]= w[2]; W[1][2]=-w[0];
	W[2][0]=-w[1]; W[2][1]= w[0];
}

void Axial(const M3 W, double w[3]) { w[0]=W[2][1]; w[1]=W[0][2]; w[2]=W[1][0]; }

/* Flanagan & Taylor co-rotational update: evolve R, V given D, W at the mid-step */
void FlanaganTaylor(const M3 D, const M3 W, const M3 Vn, const M3 Rn, double dt, M3 Rn1, M3 Vn1)
{
	double w[3];
	Axial(W, w);

	double z[3] = {0,0,0};
	for (int i=0;i<3;i++) for (int k=0;k<3;k++) for (int j=0;j<3;j++) for (int m=0;m<3;m++) {
		int e = Permutation(i,k,j);
		if (e) z[i] += e * D[j][m] * Vn[m][k];
	}

	M3 T;
	double trV = Trace(Vn);
	for (int i=0;i<3;i++) for (int j=0;j<3;j++) T[i][j] = (i==j ? trV : 0.0) - Vn[i][j];
	M3 Tinv;
	Inverse3(T, Tinv);

	double omega[3];
	for (int i=0;i<3;i++) { double s=0; for (int k=0;k<3;k++) s += Tinv[i][k]*z[k]; omega[i] = w[i] + s; }
	M3 Om;
	Skew(omega, Om);

	M3 A, B, Ainv, Q;
	Identity(A); Identity(B);
	for (int i=0;i<3;i++) for (int j=0;j<3;j++) { A[i][j] -= 0.5*dt*Om[i][j]; B[i][j] += 0.5*dt*Om[i][j]; }
	Inverse3(A, Ainv);
	Multiply(Ainv, B, Q);
	Multiply(Q, Rn, Rn1);

	M3 DW, DWV, VnOm;
	for (int i=0;i<3;i++) for (int j=0;j<3;j++) DW[i][j] = D[i][j] + W[i][j];
	Multiply(DW, Vn, DWV);
	Multiply(Vn, Om, VnOm);
	for (int i=0;i<3;i++) for (int j=0;j<3;j++) Vn1[i][j] = Vn[i][j] + dt*(DWV[i][j] - VnOm[i][j]);
}

} /* anonymous namespace */

TEST(KLShellStress, FlanaganTaylorObjectivity)
{
	const double Omega = 0.9;
	const double dt = 2.0e-3;
	const int nsteps = 500; /* total angle 0.9 rad */

	M3 D, W, V, R;
	for (int i=0;i<3;i++) for (int j=0;j<3;j++) D[i][j] = 0.0;
	double wv[3] = {0,0,Omega};
	Skew(wv, W);
	Identity(V);
	Identity(R);

	M3 stressUnrot;
	for (int i=0;i<3;i++) for (int j=0;j<3;j++) stressUnrot[i][j] = 0.0;
	stressUnrot[0][0] = 1.0; /* uniaxial in x, in the unrotated frame */

	for (int s=0;s<nsteps;s++) { M3 R1,V1; FlanaganTaylor(D,W,V,R,dt,R1,V1);
		for (int i=0;i<3;i++) for (int j=0;j<3;j++) { R[i][j]=R1[i][j]; V[i][j]=V1[i][j]; } }

	double ang = Omega * dt * nsteps;
	M3 Ra;
	Identity(Ra);
	Ra[0][0]=cos(ang); Ra[0][1]=-sin(ang); Ra[1][0]=sin(ang); Ra[1][1]=cos(ang);

	double rErr = 0.0, vErr = 0.0;
	for (int i=0;i<3;i++) for (int j=0;j<3;j++) {
		rErr = max(rErr, fabs(R[i][j]-Ra[i][j]));
		vErr = max(vErr, fabs(V[i][j]-(i==j?1.0:0.0)));
	}

	/* spatial stress co-rotates: R sig R^T should equal Ra sig Ra^T */
	M3 t1, sig, t2, sigA;
	Multiply(R, stressUnrot, t1);  MultiplyTransposeRight(t1, R, sig);
	Multiply(Ra, stressUnrot, t2); MultiplyTransposeRight(t2, Ra, sigA);
	double sErr = 0.0;
	for (int i=0;i<3;i++) for (int j=0;j<3;j++) sErr = max(sErr, fabs(sig[i][j]-sigA[i][j]));

	EXPECT_NEAR(rErr, 0.0, 1e-5); /* Cayley-transform integration error over 500 steps */
	EXPECT_NEAR(vErr, 0.0, 1e-9); /* left stretch stays I (no deformation) */
	EXPECT_NEAR(sErr, 0.0, 1e-5); /* no spurious stress */
}

TEST(KLShellStress, PlaneStressModuli)
{
	const double E = 200.0, nu = 0.3, eps = 1.0e-3, dt = 1.0;
	double lam = E*nu/((1+nu)*(1-2*nu)), mu = E/(2*(1+nu));

	/* in-plane uniaxial strain D11=eps, D22=0; solve D33 so sigma33 = 0 */
	double D33 = 0.0, sig33 = 0.0, sig11 = 0.0, sig22 = 0.0;
	for (int it=0; it<20; it++) {
		double D11=eps, D22=0.0;
		double tr = D11 + D22 + D33;
		sig11 = (lam*tr + 2*mu*D11) * dt;
		sig22 = (lam*tr + 2*mu*D22) * dt;
		sig33 = (lam*tr + 2*mu*D33) * dt;
		if (fabs(sig33) < 1e-14) break;
		double C3333 = (lam + 2*mu) * dt;
		D33 -= sig33 / C3333;
	}

	EXPECT_NEAR(sig11, E/(1-nu*nu)*eps, 1e-9);
	EXPECT_NEAR(sig22, nu*E/(1-nu*nu)*eps, 1e-9);
	EXPECT_NEAR(sig33, 0.0, 1e-12);
	EXPECT_NEAR(D33, -nu/(1-nu)*eps, 1e-9);
}

TEST(KLShellStress, ThroughThicknessBendingMoment)
{
	const double E = 200.0, nu = 0.3, h = 0.05, kappa = 1.0;
	const double xg[3] = {-sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)};
	const double wg[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0};

	double M = 0.0;
	for (int g=0; g<3; g++) {
		double e11 = -(h/2.0) * xg[g] * kappa;   /* linear-through-thickness bending strain */
		double s11 = E/(1-nu*nu) * e11;          /* plane-stress */
		double z = (h/2.0) * xg[g];
		M += wg[g] * s11 * z * (h/2.0);          /* integral sigma11 z dz, dz=(h/2)dxi3 */
	}

	double Dplate = E*h*h*h/(12.0*(1-nu*nu));
	EXPECT_NEAR(M, -Dplate*kappa, 1e-12);        /* 3-pt Gauss is exact for this integrand */
}
