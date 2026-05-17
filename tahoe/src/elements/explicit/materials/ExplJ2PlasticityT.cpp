/* $Id: ExplJ2PlasticityT.cpp,v 2.0 2026/05/08 samanseifi Exp $ */
/* ExplJ2PlasticityT.cpp — batch finite-strain J2 plasticity.
 *
 * Hughes-Winget incremental-objective update:
 *   f_rel = F * F_n^{-1}
 *   de = sym(f_rel - I)      (strain increment)
 *   dw = skew(f_rel - I)     (spin increment)
 *   Q  = (I - dw/2)^{-1} (I + dw/2)   (objective rotation)
 *   sig_rot = Q sig_n Q^T
 *   sig_tr  = sig_rot + 2 mu de + lambda tr(de) I
 *   (radial return on deviator if phi > 0)
 */
#include "ExplJ2PlasticityT.h"
#include <cmath>
#include <cstring>

using namespace Tahoe;

static const double SQRT_2_3 = 0.81649658092772603273;   /* sqrt(2/3) */
static const double SQRT_3_2 = 1.22474487139158904909;   /* sqrt(3/2) */

ExplJ2PlasticityT::ExplJ2PlasticityT(double mu, double kappa,
                                     double sigma_Y, double H, double density)
	: fMu(mu), fKappa(kappa),
	  fLambda(kappa - 2.0*mu/3.0),
	  fSigmaY(sigma_Y), fH(H), fDensity(density)
{
}

double ExplJ2PlasticityT::WaveSpeed(void) const
{
	/* elastic P-wave modulus (hardening does not raise wave speed) */
	double M = fKappa + 4.0*fMu/3.0;
	return sqrt(M / fDensity);
}

/*----------------------------------------------------------------------
 * InitializeHistory — set F_n = I for all IP-element pairs.
 * Layout: history[ip * nhist * total + var * total + elem]
 * Vars 0,4,8 are F11,F22,F33 (diagonals of F_n = I).
 *----------------------------------------------------------------------*/
void ExplJ2PlasticityT::InitializeHistory(int nip, int total, double* history) const
{
	const int nhist = NumHistoryVars();
	for (int ip = 0; ip < nip; ip++) {
		double* base = history + ip * nhist * total;
		/* F11 = F22 = F33 = 1 */
		for (int e = 0; e < total; e++) base[0 * total + e] = 1.0;
		for (int e = 0; e < total; e++) base[4 * total + e] = 1.0;
		for (int e = 0; e < total; e++) base[8 * total + e] = 1.0;
		/* everything else (sigma_n, eps_p, off-diag F) already zero */
	}
}

/*----------------------------------------------------------------------
 * 2D plane strain — NOT IMPLEMENTED (fall through to elastic trial only)
 *----------------------------------------------------------------------*/
void ExplJ2PlasticityT::ComputeStress2D(
	int nel,
	const double F[][MVSIZ],
	double sig11[], double sig22[], double sig12[],
	double* history) const
{
	/* placeholder: treat as small-strain linear elastic for now.
	 * Full 2D plane strain J2 is a follow-up. */
	double mu = fMu;
	double lam = fLambda;
	for (int i = 0; i < nel; i++) {
		double f11 = F[0][i], f12 = F[1][i];
		double f21 = F[2][i], f22 = F[3][i];
		double e11 = f11 - 1.0;
		double e22 = f22 - 1.0;
		double e12 = 0.5*(f12 + f21);
		double trE = e11 + e22;
		sig11[i] = lam*trE + 2.0*mu*e11;
		sig22[i] = lam*trE + 2.0*mu*e22;
		sig12[i] = 2.0*mu*e12;
	}
}

/*----------------------------------------------------------------------
 * 3D — Hughes-Winget finite-strain J2 with isotropic hardening
 *
 * History layout (per IP, SoA: history[var * MVSIZ + i]):
 *   [0-8]   F_n (row-major): F11,F12,F13,F21,F22,F23,F31,F32,F33
 *   [9-14]  sigma_n (Voigt): s11,s22,s33,s23,s13,s12
 *   [15]    eps_p
 *----------------------------------------------------------------------*/
void ExplJ2PlasticityT::ComputeStress3D(
	int nel,
	const double F[][MVSIZ],
	double sig[][MVSIZ],
	double* history) const
{
	const double mu = fMu;
	const double lam = fLambda;
	const double sigY0 = fSigmaY;
	const double H = fH;
	const double three_mu_H = 3.0*mu + H;

	/* SIMD-friendly inner loop */
	for (int i = 0; i < nel; i++) {

		/* ----- load F_n from history ----- */
		double Fn[9];
		for (int k = 0; k < 9; k++)
			Fn[k] = history[k * MVSIZ + i];

		/* ----- invert F_n analytically (3x3) ----- */
		double det = Fn[0]*(Fn[4]*Fn[8] - Fn[5]*Fn[7])
		           - Fn[1]*(Fn[3]*Fn[8] - Fn[5]*Fn[6])
		           + Fn[2]*(Fn[3]*Fn[7] - Fn[4]*Fn[6]);
		double invdet = 1.0/det;
		double Fi[9];
		Fi[0] =  (Fn[4]*Fn[8] - Fn[5]*Fn[7]) * invdet;
		Fi[1] = -(Fn[1]*Fn[8] - Fn[2]*Fn[7]) * invdet;
		Fi[2] =  (Fn[1]*Fn[5] - Fn[2]*Fn[4]) * invdet;
		Fi[3] = -(Fn[3]*Fn[8] - Fn[5]*Fn[6]) * invdet;
		Fi[4] =  (Fn[0]*Fn[8] - Fn[2]*Fn[6]) * invdet;
		Fi[5] = -(Fn[0]*Fn[5] - Fn[2]*Fn[3]) * invdet;
		Fi[6] =  (Fn[3]*Fn[7] - Fn[4]*Fn[6]) * invdet;
		Fi[7] = -(Fn[0]*Fn[7] - Fn[1]*Fn[6]) * invdet;
		Fi[8] =  (Fn[0]*Fn[4] - Fn[1]*Fn[3]) * invdet;

		/* ----- relative deformation: f_rel = F * F_n^{-1} ----- */
		double f11 = F[0][i]*Fi[0] + F[1][i]*Fi[3] + F[2][i]*Fi[6];
		double f12 = F[0][i]*Fi[1] + F[1][i]*Fi[4] + F[2][i]*Fi[7];
		double f13 = F[0][i]*Fi[2] + F[1][i]*Fi[5] + F[2][i]*Fi[8];
		double f21 = F[3][i]*Fi[0] + F[4][i]*Fi[3] + F[5][i]*Fi[6];
		double f22 = F[3][i]*Fi[1] + F[4][i]*Fi[4] + F[5][i]*Fi[7];
		double f23 = F[3][i]*Fi[2] + F[4][i]*Fi[5] + F[5][i]*Fi[8];
		double f31 = F[6][i]*Fi[0] + F[7][i]*Fi[3] + F[8][i]*Fi[6];
		double f32 = F[6][i]*Fi[1] + F[7][i]*Fi[4] + F[8][i]*Fi[7];
		double f33 = F[6][i]*Fi[2] + F[7][i]*Fi[5] + F[8][i]*Fi[8];

		/* ----- strain increment de = sym(f_rel - I) ----- */
		double de11 = f11 - 1.0;
		double de22 = f22 - 1.0;
		double de33 = f33 - 1.0;
		double de23 = 0.5*(f23 + f32);
		double de13 = 0.5*(f13 + f31);
		double de12 = 0.5*(f12 + f21);

		/* ----- spin half-increment dw/2 = 0.5 * skew(f_rel - I) ----- */
		double w23 = 0.25*(f23 - f32);  /* already divided by 2 */
		double w13 = 0.25*(f13 - f31);
		double w12 = 0.25*(f12 - f21);

		/* ----- Hughes-Winget rotation: Q = (I - dw/2)^{-1} (I + dw/2)
		 * For skew-symmetric W = [[0,-w12,w13],[w12,0,-w23],[-w13,w23,0]],
		 * Q = (I - W)^{-1}(I + W) is the Cayley transform.
		 * We compute Q directly by solving two 3x3 systems but for a small
		 * skew matrix, expand to second order:
		 *   Q ≈ I + 2*W*(I + W + W^2 + ...) / (det stuff)
		 * For efficiency and accuracy we compute Q via explicit inverse.
		 *
		 * Let A = I + dw/2, B = I - dw/2.  Then Q = B^{-1} A.
		 * B has form:
		 *   [[ 1,   w12, -w13],
		 *    [-w12,  1,   w23],
		 *    [ w13, -w23,  1 ]]
		 * ----- */
		double B[9] = { 1.0,  w12, -w13,
		               -w12,  1.0,  w23,
		                w13, -w23,  1.0};

		double Bdet = B[0]*(B[4]*B[8] - B[5]*B[7])
		            - B[1]*(B[3]*B[8] - B[5]*B[6])
		            + B[2]*(B[3]*B[7] - B[4]*B[6]);
		double Bi = 1.0/Bdet;
		double Binv[9];
		Binv[0] =  (B[4]*B[8] - B[5]*B[7]) * Bi;
		Binv[1] = -(B[1]*B[8] - B[2]*B[7]) * Bi;
		Binv[2] =  (B[1]*B[5] - B[2]*B[4]) * Bi;
		Binv[3] = -(B[3]*B[8] - B[5]*B[6]) * Bi;
		Binv[4] =  (B[0]*B[8] - B[2]*B[6]) * Bi;
		Binv[5] = -(B[0]*B[5] - B[2]*B[3]) * Bi;
		Binv[6] =  (B[3]*B[7] - B[4]*B[6]) * Bi;
		Binv[7] = -(B[0]*B[7] - B[1]*B[6]) * Bi;
		Binv[8] =  (B[0]*B[4] - B[1]*B[3]) * Bi;

		/* A = I + dw/2 */
		double A[9] = { 1.0, -w12,  w13,
		                w12,  1.0, -w23,
		               -w13,  w23,  1.0};

		/* Q = Binv * A */
		double Q[9];
		Q[0] = Binv[0]*A[0] + Binv[1]*A[3] + Binv[2]*A[6];
		Q[1] = Binv[0]*A[1] + Binv[1]*A[4] + Binv[2]*A[7];
		Q[2] = Binv[0]*A[2] + Binv[1]*A[5] + Binv[2]*A[8];
		Q[3] = Binv[3]*A[0] + Binv[4]*A[3] + Binv[5]*A[6];
		Q[4] = Binv[3]*A[1] + Binv[4]*A[4] + Binv[5]*A[7];
		Q[5] = Binv[3]*A[2] + Binv[4]*A[5] + Binv[5]*A[8];
		Q[6] = Binv[6]*A[0] + Binv[7]*A[3] + Binv[8]*A[6];
		Q[7] = Binv[6]*A[1] + Binv[7]*A[4] + Binv[8]*A[7];
		Q[8] = Binv[6]*A[2] + Binv[7]*A[5] + Binv[8]*A[8];

		/* ----- rotate sigma_n: sigma_rot = Q sigma_n Q^T -----
		 * sigma_n is symmetric; Voigt order [s11,s22,s33,s23,s13,s12] */
		double sn11 = history[ 9*MVSIZ + i];
		double sn22 = history[10*MVSIZ + i];
		double sn33 = history[11*MVSIZ + i];
		double sn23 = history[12*MVSIZ + i];
		double sn13 = history[13*MVSIZ + i];
		double sn12 = history[14*MVSIZ + i];

		/* T = Q * sigma_n (not symmetric) */
		double T11 = Q[0]*sn11 + Q[1]*sn12 + Q[2]*sn13;
		double T12 = Q[0]*sn12 + Q[1]*sn22 + Q[2]*sn23;
		double T13 = Q[0]*sn13 + Q[1]*sn23 + Q[2]*sn33;
		double T21 = Q[3]*sn11 + Q[4]*sn12 + Q[5]*sn13;
		double T22 = Q[3]*sn12 + Q[4]*sn22 + Q[5]*sn23;
		double T23 = Q[3]*sn13 + Q[4]*sn23 + Q[5]*sn33;
		double T31 = Q[6]*sn11 + Q[7]*sn12 + Q[8]*sn13;
		double T32 = Q[6]*sn12 + Q[7]*sn22 + Q[8]*sn23;
		double T33 = Q[6]*sn13 + Q[7]*sn23 + Q[8]*sn33;

		/* sigma_rot = T * Q^T (symmetric, six independent components) */
		double sr11 = T11*Q[0] + T12*Q[1] + T13*Q[2];
		double sr22 = T21*Q[3] + T22*Q[4] + T23*Q[5];
		double sr33 = T31*Q[6] + T32*Q[7] + T33*Q[8];
		double sr23 = T21*Q[6] + T22*Q[7] + T23*Q[8];
		double sr13 = T11*Q[6] + T12*Q[7] + T13*Q[8];
		double sr12 = T11*Q[3] + T12*Q[4] + T13*Q[5];

		/* ----- elastic trial: sigma_tr = sigma_rot + C:de ----- */
		double trDe = de11 + de22 + de33;
		double lamTr = lam * trDe;
		double st11 = sr11 + lamTr + 2.0*mu*de11;
		double st22 = sr22 + lamTr + 2.0*mu*de22;
		double st33 = sr33 + lamTr + 2.0*mu*de33;
		double st23 = sr23 + 2.0*mu*de23;
		double st13 = sr13 + 2.0*mu*de13;
		double st12 = sr12 + 2.0*mu*de12;

		/* ----- deviator / pressure split ----- */
		double p = (st11 + st22 + st33) / 3.0;
		double s11 = st11 - p;
		double s22 = st22 - p;
		double s33 = st33 - p;
		/* off-diag are deviator in Cauchy form */

		/* ----- yield check: q = sqrt(3/2 s:s), phi = q - sigma_Y(eps_p) ----- */
		double eps_p = history[15*MVSIZ + i];
		double sig_Y = sigY0 + H * eps_p;

		double s2 = s11*s11 + s22*s22 + s33*s33
		          + 2.0*(st23*st23 + st13*st13 + st12*st12);
		double q = sqrt(1.5 * s2);
		double phi = q - sig_Y;

		double new_s11, new_s22, new_s33, new_s23, new_s13, new_s12;
		double new_eps_p = eps_p;

		if (phi > 0.0) {
			/* radial return: n_hat = dev_tr / ||dev_tr||_F = sqrt(2/3)*3*s/(2q)
			 * Classical form: sigma_new = sigma_rot + C:de_elastic
			 * Equivalent deviatoric update with H-W:
			 *   dlam = phi / (3 mu + H)
			 *   s_new = s_trial * (1 - 3 mu dlam / q) */
			double dlam = phi / three_mu_H;
			double factor = 1.0 - 3.0*mu*dlam / q;
			new_s11 = s11 * factor;
			new_s22 = s22 * factor;
			new_s33 = s33 * factor;
			new_s23 = st23 * factor;
			new_s13 = st13 * factor;
			new_s12 = st12 * factor;
			new_eps_p = eps_p + dlam;
		} else {
			new_s11 = s11;
			new_s22 = s22;
			new_s33 = s33;
			new_s23 = st23;
			new_s13 = st13;
			new_s12 = st12;
		}

		/* ----- recombine: sigma = s_new + p I ----- */
		double new_sig11 = new_s11 + p;
		double new_sig22 = new_s22 + p;
		double new_sig33 = new_s33 + p;

		/* output */
		sig[0][i] = new_sig11;
		sig[1][i] = new_sig22;
		sig[2][i] = new_sig33;
		sig[3][i] = new_s23;
		sig[4][i] = new_s13;
		sig[5][i] = new_s12;

		/* ----- update history ----- */
		history[ 0*MVSIZ + i] = F[0][i];
		history[ 1*MVSIZ + i] = F[1][i];
		history[ 2*MVSIZ + i] = F[2][i];
		history[ 3*MVSIZ + i] = F[3][i];
		history[ 4*MVSIZ + i] = F[4][i];
		history[ 5*MVSIZ + i] = F[5][i];
		history[ 6*MVSIZ + i] = F[6][i];
		history[ 7*MVSIZ + i] = F[7][i];
		history[ 8*MVSIZ + i] = F[8][i];
		history[ 9*MVSIZ + i] = new_sig11;
		history[10*MVSIZ + i] = new_sig22;
		history[11*MVSIZ + i] = new_sig33;
		history[12*MVSIZ + i] = new_s23;
		history[13*MVSIZ + i] = new_s13;
		history[14*MVSIZ + i] = new_s12;
		history[15*MVSIZ + i] = new_eps_p;
	}
}
