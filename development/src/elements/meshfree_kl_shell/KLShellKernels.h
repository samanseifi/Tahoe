/* KLShellKernels.h — meshfree Kirchhoff-Love shell math kernels (epic #59)
 *
 * Header-only, dependency-free (plain C++) implementations of the Wang & Bazilevs
 * (2025) meshfree Kirchhoff-Love shell kernels, used by RKShellT and validated by the
 * gtest suite tests/meshfree/test_KLShell* (ctest -R KLShell):
 *
 *   - PCA local tangent-plane parameterization (symmetric 3x3 Jacobi eigensolver);
 *   - shell geometry tensors A, B1, B2 and their parametric derivatives (BuildGeom);
 *   - per-node strain-displacement matrix B_Iijk (BMatrix, Eq. 51) + Voigt form;
 *   - parametric gradient B_Iijkl at xi3=0 (BMatrixGradient, Eq. 54) for stabilization;
 *   - compressible Neo-Hookean Cauchy stress with sigma33=0 plane stress.
 *
 * Conventions: 3-vectors are double[3] = (x,y,z); second-derivative rows are ordered
 * [xi1xi1, xi2xi2, xi1xi2]; Voigt strain rows are [11,22,33,23,13,12] (engineering shear).
 */
#ifndef _KL_SHELL_KERNELS_H_
#define _KL_SHELL_KERNELS_H_

#include <cmath>

namespace Tahoe {
namespace KLShell {

/* ---- small vector / tensor helpers ---------------------------------------- */

inline double Dot(const double a[3], const double b[3])
{
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline void Cross(const double a[3], const double b[3], double out[3])
{
	out[0] = a[1]*b[2] - a[2]*b[1];
	out[1] = a[2]*b[0] - a[0]*b[2];
	out[2] = a[0]*b[1] - a[1]*b[0];
}

inline double Norm(const double a[3]) { return std::sqrt(Dot(a, a)); }

/* Levi-Civita permutation symbol */
inline int Permutation(int i, int j, int k)
{
	if (i == j || j == k || i == k) return 0;
	if ((i == 0 && j == 1 && k == 2) ||
	    (i == 1 && j == 2 && k == 0) ||
	    (i == 2 && j == 0 && k == 1)) return 1;
	return -1;
}

/* inverse of a 3x3 matrix; returns false if singular */
inline bool Inverse(const double M[3][3], double Inv[3][3])
{
	double det = M[0][0]*(M[1][1]*M[2][2] - M[1][2]*M[2][1])
	           - M[0][1]*(M[1][0]*M[2][2] - M[1][2]*M[2][0])
	           + M[0][2]*(M[1][0]*M[2][1] - M[1][1]*M[2][0]);
	if (std::fabs(det) < 1e-300) return false;
	double idet = 1.0/det;

	Inv[0][0] =  (M[1][1]*M[2][2] - M[1][2]*M[2][1])*idet;
	Inv[0][1] = -(M[0][1]*M[2][2] - M[0][2]*M[2][1])*idet;
	Inv[0][2] =  (M[0][1]*M[1][2] - M[0][2]*M[1][1])*idet;
	Inv[1][0] = -(M[1][0]*M[2][2] - M[1][2]*M[2][0])*idet;
	Inv[1][1] =  (M[0][0]*M[2][2] - M[0][2]*M[2][0])*idet;
	Inv[1][2] = -(M[0][0]*M[1][2] - M[0][2]*M[1][0])*idet;
	Inv[2][0] =  (M[1][0]*M[2][1] - M[1][1]*M[2][0])*idet;
	Inv[2][1] = -(M[0][0]*M[2][1] - M[0][1]*M[2][0])*idet;
	Inv[2][2] =  (M[0][0]*M[1][1] - M[0][1]*M[1][0])*idet;
	return true;
}

/* symmetric 3x3 Jacobi eigensolver. eval[i] = i-th eigenvalue;
 * evec[i] = i-th eigenvector (stored as a row). */
inline void JacobiEigen3(const double Ain[3][3], double eval[3], double evec[3][3])
{
	double a[3][3];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			a[i][j] = Ain[i][j];

	double v[3][3] = {{1,0,0}, {0,1,0}, {0,0,1}};

	for (int sweep = 0; sweep < 100; sweep++) {

		double off = std::fabs(a[0][1]) + std::fabs(a[0][2]) + std::fabs(a[1][2]);
		if (off < 1e-18) break;

		for (int p = 0; p < 3; p++)
			for (int q = p+1; q < 3; q++) {

				if (std::fabs(a[p][q]) < 1e-300) continue;

				double theta = (a[q][q] - a[p][p])/(2.0*a[p][q]);
				double t = (theta >= 0 ? 1.0 : -1.0)/(std::fabs(theta) + std::sqrt(theta*theta + 1.0));
				double c = 1.0/std::sqrt(t*t + 1.0);
				double sn = t*c;

				for (int k = 0; k < 3; k++) {
					double x = a[k][p], y = a[k][q];
					a[k][p] = c*x - sn*y;
					a[k][q] = sn*x + c*y;
				}
				for (int k = 0; k < 3; k++) {
					double x = a[p][k], y = a[q][k];
					a[p][k] = c*x - sn*y;
					a[q][k] = sn*x + c*y;
				}
				for (int k = 0; k < 3; k++) {
					double x = v[k][p], y = v[k][q];
					v[k][p] = c*x - sn*y;
					v[k][q] = sn*x + c*y;
				}
			}
	}

	for (int i = 0; i < 3; i++) {
		eval[i] = a[i][i];
		for (int k = 0; k < 3; k++)
			evec[i][k] = v[k][i];
	}
}

/* PCA tangent frame from neighbor reference positions (X is nn x 3, row-major).
 * Returns the orthonormal tangent basis psi1, psi2 (the largest-variance plane)
 * and the surface normal n0 = psi1 x psi2. */
inline void PCAFrame(const double* X, int nn, double psi1[3], double psi2[3], double n0[3])
{
	double mean[3] = {0,0,0};
	for (int k = 0; k < nn; k++)
		for (int d = 0; d < 3; d++)
			mean[d] += X[3*k+d];
	for (int d = 0; d < 3; d++)
		mean[d] /= nn;

	double C[3][3] = {{0,0,0}, {0,0,0}, {0,0,0}};
	for (int k = 0; k < nn; k++) {
		double dd[3] = {X[3*k] - mean[0], X[3*k+1] - mean[1], X[3*k+2] - mean[2]};
		for (int a = 0; a < 3; a++)
			for (int b = 0; b < 3; b++)
				C[a][b] += dd[a]*dd[b];
	}

	double ev[3], evec[3][3];
	JacobiEigen3(C, ev, evec);

	/* order eigenvectors by decreasing eigenvalue */
	int o[3] = {0,1,2};
	for (int a = 0; a < 3; a++)
		for (int b = a+1; b < 3; b++)
			if (ev[o[b]] > ev[o[a]]) { int t = o[a]; o[a] = o[b]; o[b] = t; }

	for (int d = 0; d < 3; d++) {
		psi1[d] = evec[o[0]][d];
		psi2[d] = evec[o[1]][d];
	}
	Cross(psi1, psi2, n0);
	double nm = Norm(n0);
	for (int d = 0; d < 3; d++)
		n0[d] /= nm;
}

/* shell geometry tensors at a sample point for a through-thickness station xi3 */
struct ShellGeom
{
	double n[3];            /**< unit surface normal */
	double B1[3][3];        /**< auxiliary tensor B1 (Eq. 40) */
	double B2[3][3];        /**< auxiliary tensor B2 (Eq. 40) */
	double B1m[2][3][3];    /**< parametric derivatives B1,xi_m (m = 0,1) */
	double B2m[2][3][3];    /**< parametric derivatives B2,xi_m (m = 0,1) */
	double F[3][3];         /**< 3D Jacobian F3D = d x3D / d xi (Eq. 18) */
	double Finv[3][3];      /**< inverse of F */
	double h, xi3;          /**< shell thickness, through-thickness station */
	double x3D_ab[2][2][3]; /**< x3D,xi_a,xi_b (reference, at xi3) */
	double n_l[2][3];       /**< parametric derivatives of the normal n,xi_l */
};

/* build the shell geometry tensors from the parametric position derivatives
 * x1 = x,xi1, x2 = x,xi2, x11 = x,xi1xi1, x22 = x,xi2xi2, x12 = x,xi1xi2. */
inline bool BuildGeom(const double x1[3], const double x2[3], const double x11[3],
                      const double x22[3], const double x12[3], double h, double xi3,
                      ShellGeom& g)
{
	g.h = h;
	g.xi3 = xi3;

	double cr[3];
	Cross(x1, x2, cr);
	double J = Norm(cr);
	if (J < 1e-300) return false;
	for (int d = 0; d < 3; d++)
		g.n[d] = cr[d]/J;

	/* A = (I - n (x) n) / J  (Eq. 41) */
	double A[3][3];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			A[i][j] = ((i == j ? 1.0 : 0.0) - g.n[i]*g.n[j])/J;

	/* parametric derivatives of (x,xi1 x x,xi2) */
	double dcr[2][3];
	{
		double t1[3], t2[3];
		Cross(x11, x2, t1); Cross(x1, x12, t2);
		for (int d = 0; d < 3; d++) dcr[0][d] = t1[d] + t2[d];
		Cross(x12, x2, t1); Cross(x1, x22, t2);
		for (int d = 0; d < 3; d++) dcr[1][d] = t1[d] + t2[d];
	}

	/* n,xi_m = (I - n (x) n)/J . dcr_m ;  J,xi_m = n . dcr_m */
	double Jm[2];
	for (int m = 0; m < 2; m++) {
		Jm[m] = Dot(g.n, dcr[m]);
		for (int i = 0; i < 3; i++) {
			double s = 0.0;
			for (int k = 0; k < 3; k++)
				s += ((i == k ? 1.0 : 0.0) - g.n[i]*g.n[k])*dcr[m][k];
			g.n_l[m][i] = s/J;
		}
	}

	/* A,xi_m (Eq. 45) */
	double Am[2][3][3];
	for (int m = 0; m < 2; m++)
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				Am[m][i][j] = -(g.n_l[m][i]*g.n[j] + g.n[i]*g.n_l[m][j])/J - A[i][j]*Jm[m]/J;

	/* B1_ik = A_ip eps_pkl x2_l ;  B2_ik = A_ip eps_pkl x1_k  (Eq. 40) */
	for (int i = 0; i < 3; i++)
		for (int k = 0; k < 3; k++) {
			double s1 = 0.0, s2 = 0.0;
			for (int p = 0; p < 3; p++)
				for (int l = 0; l < 3; l++) {
					int e = Permutation(p, k, l);
					if (e) s1 += A[i][p]*e*x2[l];
				}
			for (int p = 0; p < 3; p++)
				for (int kk = 0; kk < 3; kk++) {
					int e = Permutation(p, kk, k);
					if (e) s2 += A[i][p]*e*x1[kk];
				}
			g.B1[i][k] = s1;
			g.B2[i][k] = s2;
		}

	/* parametric derivatives B1,xi_m and B2,xi_m (Eq. 44) */
	const double* x2m[2] = {x12, x22};
	const double* x1m[2] = {x11, x12};
	for (int m = 0; m < 2; m++)
		for (int i = 0; i < 3; i++)
			for (int k = 0; k < 3; k++) {
				double s1 = 0.0, s2 = 0.0;
				for (int p = 0; p < 3; p++)
					for (int l = 0; l < 3; l++) {
						int e = Permutation(p, k, l);
						if (e) s1 += e*(Am[m][i][p]*x2[l] + A[i][p]*x2m[m][l]);
					}
				for (int p = 0; p < 3; p++)
					for (int kk = 0; kk < 3; kk++) {
						int e = Permutation(p, kk, k);
						if (e) s2 += e*(Am[m][i][p]*x1[kk] + A[i][p]*x1m[m][kk]);
					}
				g.B1m[m][i][k] = s1;
				g.B2m[m][i][k] = s2;
			}

	/* in-plane columns of F3D (Eq. 15) and the through-thickness column (h/2) n */
	const double* xmm1[2] = {x11, x12};
	const double* xmm2[2] = {x12, x22};
	const double* x_a[2]  = {x1, x2};
	for (int a = 0; a < 2; a++)
		for (int i = 0; i < 3; i++) {
			double t = 0.0;
			for (int k = 0; k < 3; k++)
				t += g.B1[i][k]*xmm1[a][k] + g.B2[i][k]*xmm2[a][k];
			g.F[i][a] = x_a[a][i] + (h/2.0)*xi3*t;
		}
	for (int i = 0; i < 3; i++)
		g.F[i][2] = (h/2.0)*g.n[i];
	if (!Inverse(g.F, g.Finv)) return false;

	const double* xab[2][2] = {{x11, x12}, {x12, x22}};
	for (int a = 0; a < 2; a++)
		for (int b = 0; b < 2; b++)
			for (int d = 0; d < 3; d++)
				g.x3D_ab[a][b][d] = xab[a][b][d];
	return true;
}

/* per-node strain-displacement operator B_Iijk (Eq. 51): grad(v3D)_ij = sum_I B[i][j][k] v_Ik.
 * The P* are the node's RK shape-function derivatives:
 *   P1 = Psi,xi1, P2 = Psi,xi2, P11 = Psi,xi1xi1, P12 = Psi,xi1xi2, P22 = Psi,xi2xi2. */
inline void BMatrix(const ShellGeom& g, double P1, double P2, double P11, double P12,
                    double P22, double B[3][3][3])
{
	double hx = (g.h/2.0)*g.xi3;
	for (int i = 0; i < 3; i++)
		for (int k = 0; k < 3; k++) {
			double d_ik = (i == k ? 1.0 : 0.0);
			double br0 = d_ik*P1
			           + hx*(g.B1m[0][i][k]*P1 + g.B1[i][k]*P11 + g.B2m[0][i][k]*P2 + g.B2[i][k]*P12);
			double br1 = d_ik*P2
			           + hx*(g.B1m[1][i][k]*P1 + g.B1[i][k]*P12 + g.B2m[1][i][k]*P2 + g.B2[i][k]*P22);
			double br2 = (g.h/2.0)*(g.B1[i][k]*P1 + g.B2[i][k]*P2);
			for (int j = 0; j < 3; j++)
				B[i][j][k] = br0*g.Finv[0][j] + br1*g.Finv[1][j] + br2*g.Finv[2][j];
		}
}

/* parametric gradient of the strain-displacement operator at xi3=0 (Eq. 54). The
 * xi3*S third-derivative term vanishes at xi3=0 (the natural stabilization uses one-point
 * through-thickness quadrature there), so no third derivatives are needed. Bz must be
 * BMatrix evaluated with g0 at xi3=0; P1l[l]=Psi,xi1xi_l, P2l[l]=Psi,xi2xi_l.
 * Output Bg[i][j][k][l] for l in {0,1}. */
inline void BMatrixGradient(const ShellGeom& g0, double P1, double P2,
                            double P1l[2], double P2l[2], const double Bz[3][3][3],
                            double Bg[3][3][3][2])
{
	for (int l = 0; l < 2; l++) {

		double Pa_l[2] = {P1l[l], P2l[l]};

		for (int i = 0; i < 3; i++)
			for (int k = 0; k < 3; k++) {

				/* in-plane (a = xi1, xi2) contributions */
				double br[3];
				for (int a = 0; a < 2; a++) {
					double Bdx = 0.0;
					for (int m = 0; m < 3; m++)
						Bdx += Bz[i][m][k]*g0.x3D_ab[a][l][m];
					br[a] = (i == k ? 1.0 : 0.0)*Pa_l[a] - Bdx;
				}

				/* through-thickness (xi3) contribution */
				double Bdx3 = 0.0;
				for (int m = 0; m < 3; m++)
					Bdx3 += Bz[i][m][k]*(g0.h/2.0)*g0.n_l[l][m];
				double t3 = (g0.h/2.0)*(g0.B1[i][k]*P1l[l] + g0.B1m[l][i][k]*P1
				                      + g0.B2[i][k]*P2l[l] + g0.B2m[l][i][k]*P2) - Bdx3;

				for (int j = 0; j < 3; j++)
					Bg[i][j][k][l] = br[0]*g0.Finv[0][j] + br[1]*g0.Finv[1][j] + t3*g0.Finv[2][j];
			}
	}
}

/* Voigt 6x3 form of B[3][3][3]; rows [11,22,33,23,13,12] with engineering shear */
inline void ToVoigt(const double B[3][3][3], double Bv[6][3])
{
	for (int k = 0; k < 3; k++) {
		Bv[0][k] = B[0][0][k];
		Bv[1][k] = B[1][1][k];
		Bv[2][k] = B[2][2][k];
		Bv[3][k] = B[1][2][k] + B[2][1][k];
		Bv[4][k] = B[0][2][k] + B[2][0][k];
		Bv[5][k] = B[0][1][k] + B[1][0][k];
	}
}

/* compressible Neo-Hookean Cauchy stress from the deformation gradient F (3x3):
 *   sigma = (mu/J)(b - I) + (lambda ln J / J) I,   b = F F^T,   J = det F. */
inline void NeoHookeCauchy(const double F[3][3], double lambda, double mu, double sig[3][3])
{
	double J = F[0][0]*(F[1][1]*F[2][2] - F[1][2]*F[2][1])
	         - F[0][1]*(F[1][0]*F[2][2] - F[1][2]*F[2][0])
	         + F[0][2]*(F[1][0]*F[2][1] - F[1][1]*F[2][0]);

	double b[3][3];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++) {
			double s = 0.0;
			for (int k = 0; k < 3; k++)
				s += F[i][k]*F[j][k];
			b[i][j] = s;
		}

	double c1 = mu/J;
	double c2 = lambda*std::log(J)/J;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			sig[i][j] = c1*(b[i][j] - (i == j ? 1.0 : 0.0)) + (i == j ? c2 : 0.0);
}

/* enforce sigma33 = 0 (plane stress) by a local Newton on the normal stretch F[2][2].
 * F must be expressed in a frame whose 3-direction is the shell normal. On return F[2][2]
 * is adjusted and sig holds the plane-stress Cauchy stress. */
inline void NeoHookePlaneStress(double F[3][3], double lambda, double mu, double sig[3][3])
{
	for (int it = 0; it < 20; it++) {

		NeoHookeCauchy(F, lambda, mu, sig);
		double r = sig[2][2];
		if (std::fabs(r) < 1e-12) break;

		/* d sigma33 / d F33 by forward difference */
		double F33 = F[2][2];
		double dF = 1e-7*(std::fabs(F33) + 1e-3);
		F[2][2] = F33 + dF;
		double sp[3][3];
		NeoHookeCauchy(F, lambda, mu, sp);
		F[2][2] = F33;

		double dsdF = (sp[2][2] - r)/dF;
		if (std::fabs(dsdF) < 1e-300) break;
		F[2][2] = F33 - r/dsdF;
	}
}

} /* namespace KLShell */
} /* namespace Tahoe */

#endif /* _KL_SHELL_KERNELS_H_ */
