/* ExplNeoHookeanT.cpp — batch Neo-Hookean Cauchy stress.
 *
 * Compressible Neo-Hookean:
 *   sigma = (mu/J)(b - I) + (kappa*(J-1)/J)*I
 *   b = F*F^T  (left Cauchy-Green tensor)
 *   J = det(F)
 *
 * The inner loops are written for compiler SIMD auto-vectorization.
 */
#include "ExplNeoHookeanT.h"
#include <cmath>

using namespace Tahoe;

ExplNeoHookeanT::ExplNeoHookeanT(double mu, double kappa, double density)
	: fMu(mu), fKappa(kappa), fDensity(density)
{
}

/*----------------------------------------------------------------------
 * 2D plane strain
 * F layout: F[0]=F11, F[1]=F12, F[2]=F21, F[3]=F22
 * Plane strain: F33=1, F13=F31=F23=F32=0
 *----------------------------------------------------------------------*/
void ExplNeoHookeanT::ComputeStress2D(
	int nel,
	const double F[][MVSIZ],
	double sig11[], double sig22[], double sig12[],
	double* history) const
{
	double mu = fMu;
	double kappa = fKappa;

	/* SIMD-vectorizable inner loop */
	for (int i = 0; i < nel; i++) {
		double f11 = F[0][i], f12 = F[1][i];
		double f21 = F[2][i], f22 = F[3][i];

		/* J = det(F) — plane strain: J = F11*F22 - F12*F21 (F33=1) */
		double J = f11*f22 - f12*f21;

		/* b = F*F^T (2D part) */
		double b11 = f11*f11 + f12*f12;
		double b22 = f21*f21 + f22*f22;
		double b12 = f11*f21 + f12*f22;
		/* b33 = 1 for plane strain (F33=1) */

		/* sigma = (mu/J)(b - I) + (kappa*(J-1)/J)*I */
		double muJ = mu / J;
		double pres = kappa * (J - 1.0) / J;

		sig11[i] = muJ * (b11 - 1.0) + pres;
		sig22[i] = muJ * (b22 - 1.0) + pres;
		sig12[i] = muJ * b12;
		/* sig33 = muJ*(b33 - 1) + pres = muJ*(1-1) + pres = pres
		 * (plane strain out-of-plane stress, not assembled) */
	}
}

/*----------------------------------------------------------------------
 * 3D
 * F layout: F[0..8] = F11,F12,F13,F21,F22,F23,F31,F32,F33 (row-major)
 * sig layout: sig[0..5] = s11,s22,s33,s23,s13,s12 (Voigt)
 *----------------------------------------------------------------------*/
void ExplNeoHookeanT::ComputeStress3D(
	int nel,
	const double F[][MVSIZ],
	double sig[][MVSIZ],
	double* history) const
{
	double mu = fMu;
	double kappa = fKappa;

	/* SIMD-vectorizable inner loop */
	for (int i = 0; i < nel; i++) {
		double f11 = F[0][i], f12 = F[1][i], f13 = F[2][i];
		double f21 = F[3][i], f22 = F[4][i], f23 = F[5][i];
		double f31 = F[6][i], f32 = F[7][i], f33 = F[8][i];

		/* J = det(F) */
		double J = f11*(f22*f33 - f23*f32)
		         - f12*(f21*f33 - f23*f31)
		         + f13*(f21*f32 - f22*f31);

		/* b = F*F^T */
		double b11 = f11*f11 + f12*f12 + f13*f13;
		double b22 = f21*f21 + f22*f22 + f23*f23;
		double b33 = f31*f31 + f32*f32 + f33*f33;
		double b12 = f11*f21 + f12*f22 + f13*f23;
		double b13 = f11*f31 + f12*f32 + f13*f33;
		double b23 = f21*f31 + f22*f32 + f23*f33;

		/* sigma = (mu/J)(b - I) + (kappa*(J-1)/J)*I */
		double muJ = mu / J;
		double pres = kappa * (J - 1.0) / J;

		sig[0][i] = muJ * (b11 - 1.0) + pres; /* s11 */
		sig[1][i] = muJ * (b22 - 1.0) + pres; /* s22 */
		sig[2][i] = muJ * (b33 - 1.0) + pres; /* s33 */
		sig[3][i] = muJ * b23;                 /* s23 */
		sig[4][i] = muJ * b13;                 /* s13 */
		sig[5][i] = muJ * b12;                 /* s12 */
	}
}
