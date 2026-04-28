/* test_ExplJ2Plasticity.cpp — unit tests for batch J2 plasticity.
 *
 * Tests the constitutive law in isolation (no FEM):
 *   - elastic limit (below yield, trial == elastic)
 *   - uniaxial tension: yield + linear hardening curve
 *   - rotation invariance: rotated F gives rotated sigma
 *   - InitializeHistory: F_n starts as identity
 *   - batch-vs-scalar: same result for single element at different batch slots
 */

#include "gtest/gtest.h"
#include "ExplJ2PlasticityT.h"
#include "ExplicitKernelT.h"  /* for MVSIZ */
#include <cmath>
#include <vector>

using namespace Tahoe;

static const int MVSIZ = ExplicitKernelT::MVSIZ;

/* helper: set F as identity for element i in batch */
static void SetIdentity(double F[9][MVSIZ], int i) {
	for (int k = 0; k < 9; k++) F[k][i] = 0.0;
	F[0][i] = 1.0; F[4][i] = 1.0; F[8][i] = 1.0;
}

/* helper: set F[9][MVSIZ] at slot i from a flat 3x3 row-major array */
static void SetF(double F[9][MVSIZ], int i, const double f[9]) {
	for (int k = 0; k < 9; k++) F[k][i] = f[k];
}

/* helper: material with typical steel parameters */
static ExplJ2PlasticityT MakeSteel() {
	double E = 200e3, nu = 0.3;
	double mu = E / (2.0*(1.0+nu));
	double kappa = E / (3.0*(1.0-2.0*nu));
	return ExplJ2PlasticityT(mu, kappa, /*sigma_Y*/ 250.0, /*H*/ 1000.0, /*rho*/ 7.85e-9);
}

/*=====================================================================
 * TEST 1: InitializeHistory sets F_n = I
 *===================================================================*/
TEST(ExplJ2, InitializeHistorySetsFnToIdentity) {
	ExplJ2PlasticityT mat = MakeSteel();
	int nip = 8, total = 10, nhist = mat.NumHistoryVars();
	std::vector<double> hist(nip * nhist * total, 0.0);
	mat.InitializeHistory(nip, total, hist.data());

	for (int ip = 0; ip < nip; ip++) {
		double* base = hist.data() + ip * nhist * total;
		for (int e = 0; e < total; e++) {
			EXPECT_DOUBLE_EQ(base[0*total + e], 1.0) << "F11 ip=" << ip << " e=" << e;
			EXPECT_DOUBLE_EQ(base[4*total + e], 1.0) << "F22 ip=" << ip << " e=" << e;
			EXPECT_DOUBLE_EQ(base[8*total + e], 1.0) << "F33 ip=" << ip << " e=" << e;
			EXPECT_DOUBLE_EQ(base[9*total + e],  0.0) << "sigma11";
			EXPECT_DOUBLE_EQ(base[15*total + e], 0.0) << "eps_p";
		}
	}
}

/*=====================================================================
 * TEST 2: Elastic limit — strain below yield gives purely elastic stress
 *
 * Apply small uniaxial stretch lambda = 1.001 so strain << yield.
 * Expected: sigma11 ~ E * eps, sigma22 = sigma33 = 0 for uniaxial stress.
 * But in this test we apply stretch with Poisson contraction — simplest:
 * apply infinitesimal volumetric strain and check elastic response.
 *
 * Since H-W with de = sym(f_rel - I):
 *   F = diag(1 + eps, 1 + eps, 1 + eps) → de = eps * I, dw = 0
 *   sigma_trial = (lam*3*eps + 2*mu*eps) * I = (kappa * 3 * eps) * I (pressure)
 *===================================================================*/
TEST(ExplJ2, ElasticBulkBelowYield) {
	ExplJ2PlasticityT mat = MakeSteel();
	int nhist = mat.NumHistoryVars();
	std::vector<double> hist(nhist * MVSIZ, 0.0);
	for (int i = 0; i < 1; i++) {
		hist[0*MVSIZ + i] = 1.0; hist[4*MVSIZ + i] = 1.0; hist[8*MVSIZ + i] = 1.0;
	}

	double F[9][MVSIZ];
	SetIdentity(F, 0);
	double eps = 1e-4;
	F[0][0] = 1.0 + eps; F[4][0] = 1.0 + eps; F[8][0] = 1.0 + eps;

	double sig[6][MVSIZ];
	mat.ComputeStress3D(/*nel*/1, F, sig, hist.data());

	double E = 200e3, nu = 0.3;
	double kappa = E / (3.0*(1.0-2.0*nu));
	double expected_pressure = 3.0 * kappa * eps;

	EXPECT_NEAR(sig[0][0], expected_pressure, 1e-6 * expected_pressure) << "sigma11";
	EXPECT_NEAR(sig[1][0], expected_pressure, 1e-6 * expected_pressure) << "sigma22";
	EXPECT_NEAR(sig[2][0], expected_pressure, 1e-6 * expected_pressure) << "sigma33";
	EXPECT_NEAR(sig[3][0], 0.0, 1e-10) << "sigma23";
	EXPECT_NEAR(sig[4][0], 0.0, 1e-10) << "sigma13";
	EXPECT_NEAR(sig[5][0], 0.0, 1e-10) << "sigma12";

	/* eps_p should remain 0 */
	EXPECT_NEAR(hist[15*MVSIZ + 0], 0.0, 1e-12);
}

/*=====================================================================
 * TEST 3: Uniaxial shear yielding
 *
 * Apply a simple shear increment F = I + gamma * (e1 x e2).
 * For small gamma: de12 = 0.5*gamma, others zero. Then:
 *   trial sigma has only sigma12 = 2*mu*de12 = mu*gamma (Cauchy shear)
 *   von Mises q = sqrt(3) * |sigma12|
 * Yield when q = sigma_Y → gamma_yield = sigma_Y / (mu*sqrt(3))
 *===================================================================*/
TEST(ExplJ2, PureShearYieldsAtMisesCriterion) {
	ExplJ2PlasticityT mat = MakeSteel();
	int nhist = mat.NumHistoryVars();

	double E = 200e3, nu = 0.3;
	double mu = E / (2.0*(1.0+nu));
	double sigY = 250.0;

	std::vector<double> hist(nhist * MVSIZ, 0.0);
	hist[0*MVSIZ + 0] = 1.0; hist[4*MVSIZ + 0] = 1.0; hist[8*MVSIZ + 0] = 1.0;

	/* apply shear just past yield */
	double gamma_yield = sigY / (mu * sqrt(3.0));
	double gamma = gamma_yield * 1.5;  /* 50% past yield */

	double F[9][MVSIZ];
	SetIdentity(F, 0);
	F[1][0] = gamma;   /* F12 = gamma */

	double sig[6][MVSIZ];
	mat.ComputeStress3D(1, F, sig, hist.data());

	/* check von Mises of returned stress does not exceed yield + H*eps_p */
	double s11 = sig[0][0], s22 = sig[1][0], s33 = sig[2][0];
	double s12 = sig[5][0];
	double p = (s11+s22+s33)/3.0;
	double d11 = s11-p, d22 = s22-p, d33 = s33-p;
	double q = sqrt(1.5 * (d11*d11 + d22*d22 + d33*d33 + 2*s12*s12));
	double eps_p = hist[15*MVSIZ + 0];
	double sig_Y_current = sigY + 1000.0 * eps_p;

	/* after radial return, von Mises should equal (sigma_Y + H*eps_p) */
	EXPECT_NEAR(q, sig_Y_current, 1e-4 * sig_Y_current)
		<< "Post-return q=" << q << " sig_Y_current=" << sig_Y_current;

	/* plastic strain must be positive */
	EXPECT_GT(eps_p, 0.0);
}

/*=====================================================================
 * TEST 4: Batch consistency — identical F at different slots gives identical sigma
 *===================================================================*/
TEST(ExplJ2, BatchConsistency) {
	ExplJ2PlasticityT mat = MakeSteel();
	int nhist = mat.NumHistoryVars();
	int nel = 64;
	std::vector<double> hist(nhist * MVSIZ, 0.0);
	for (int i = 0; i < nel; i++) {
		hist[0*MVSIZ + i] = 1.0; hist[4*MVSIZ + i] = 1.0; hist[8*MVSIZ + i] = 1.0;
	}

	double F[9][MVSIZ];
	double f_shear[9] = {1.0, 0.02, 0.0,
	                     0.0, 1.0,  0.0,
	                     0.0, 0.0,  1.0};
	for (int i = 0; i < nel; i++) SetF(F, i, f_shear);

	double sig[6][MVSIZ];
	mat.ComputeStress3D(nel, F, sig, hist.data());

	for (int i = 1; i < nel; i++) {
		for (int k = 0; k < 6; k++)
			EXPECT_NEAR(sig[k][i], sig[k][0], 1e-10)
				<< "slot " << i << " component " << k;
	}
}

/*=====================================================================
 * TEST 5: Rigid rotation objectivity
 *
 * Apply F = R (pure rotation, no stretch). After building up some plastic
 * state with shear, then rotating, the stress magnitude should rotate with
 * the body (Mises invariance) but not grow spuriously.
 *
 * Simplified check: start at zero stress, apply rigid rotation only.
 * The stress should remain zero through any rotation history.
 *===================================================================*/
TEST(ExplJ2, RigidRotationZeroStress) {
	ExplJ2PlasticityT mat = MakeSteel();
	int nhist = mat.NumHistoryVars();
	std::vector<double> hist(nhist * MVSIZ, 0.0);
	hist[0*MVSIZ + 0] = 1.0; hist[4*MVSIZ + 0] = 1.0; hist[8*MVSIZ + 0] = 1.0;

	double F[9][MVSIZ];
	SetIdentity(F, 0);
	double sig[6][MVSIZ];

	/* Apply SMALL rotation increments.  Hughes-Winget has O(d_theta^2) error
	 * per step from truncating sym(f_rel - I) at first order.  For d_theta =
	 * 0.001 rad (~0.057 deg), per-step de error ~ d_theta^2/2 = 5e-7, giving
	 * pressure drift ~ 3*kappa*5e-7 ~ 0.2 MPa per step.  Over 100 steps of
	 * 0.001 rad each (5.73 deg total), total drift < 20 MPa. */
	const int n_steps = 100;
	double d_theta = 0.001;  /* radians per step */
	double cumulative_theta = 0.0;
	for (int step = 0; step < n_steps; step++) {
		cumulative_theta += d_theta;
		double c = cos(cumulative_theta), s = sin(cumulative_theta);
		F[0][0] =  c; F[1][0] = -s; F[2][0] = 0.0;
		F[3][0] =  s; F[4][0] =  c; F[5][0] = 0.0;
		F[6][0] = 0.0; F[7][0] = 0.0; F[8][0] = 1.0;

		mat.ComputeStress3D(1, F, sig, hist.data());
	}

	/* Final stress at 5.73 deg rotation: drift < 30 MPa tolerates the
	 * known quadratic error.  This is a SANITY check for no catastrophic
	 * failure; strict rotation invariance requires exact Hencky or polar
	 * decomposition (tracked in #24 follow-up). */
	double max_sig = 0.0;
	for (int k = 0; k < 6; k++) max_sig = std::max(max_sig, fabs(sig[k][0]));
	EXPECT_LT(max_sig, 30.0)
		<< "after " << cumulative_theta*180/M_PI << " deg rotation: max sigma="
		<< max_sig << " MPa (H-W quadratic drift, documented).";

	/* eps_p should never grow from rigid rotation alone */
	EXPECT_NEAR(hist[15*MVSIZ + 0], 0.0, 1e-10);
}

/*=====================================================================
 * TEST 6: Uniaxial stretch past yield, linear hardening curve
 *
 * Apply axial stretch (incompressible flow approximation via Poisson) and
 * check sigma vs eps_p follows sigma_Y + H*eps_p.
 *===================================================================*/
TEST(ExplJ2, UniaxialStretchHardeningCurve) {
	double E = 117e3, nu = 0.35;
	double mu = E / (2.0*(1.0+nu));
	double kappa = E / (3.0*(1.0-2.0*nu));
	double sigY = 90.0, H = 150.0;
	ExplJ2PlasticityT mat(mu, kappa, sigY, H, 8.96e-9);
	int nhist = mat.NumHistoryVars();
	std::vector<double> hist(nhist * MVSIZ, 0.0);
	hist[0*MVSIZ + 0] = 1.0; hist[4*MVSIZ + 0] = 1.0; hist[8*MVSIZ + 0] = 1.0;

	double F[9][MVSIZ];
	SetIdentity(F, 0);
	double sig[6][MVSIZ];

	/* ramp axial stretch 1.0 -> 1.10 in 200 small steps
	 * For plastic incompressibility, lateral stretch ~ 1/sqrt(lambda_z). */
	const int n_steps = 200;
	double lam_final = 1.10;

	double last_q = 0.0;
	double last_eps_p = 0.0;
	for (int step = 1; step <= n_steps; step++) {
		double lam_z = 1.0 + (lam_final - 1.0) * step / n_steps;
		double lam_xy = 1.0/sqrt(lam_z);  /* isochoric */
		F[0][0] = lam_xy; F[4][0] = lam_xy; F[8][0] = lam_z;
		mat.ComputeStress3D(1, F, sig, hist.data());

		/* von Mises */
		double s11=sig[0][0], s22=sig[1][0], s33=sig[2][0];
		double p = (s11+s22+s33)/3.0;
		double d11 = s11-p, d22 = s22-p, d33 = s33-p;
		double s23=sig[3][0], s13=sig[4][0], s12=sig[5][0];
		double q = sqrt(1.5*(d11*d11+d22*d22+d33*d33 + 2*(s23*s23+s13*s13+s12*s12)));
		double eps_p = hist[15*MVSIZ + 0];
		last_q = q; last_eps_p = eps_p;
	}

	/* after full ramp, stress should be on yield surface sigma_Y + H*eps_p */
	double expected_q = sigY + H*last_eps_p;
	EXPECT_NEAR(last_q, expected_q, 0.05*expected_q)  /* 5% tolerance */
		<< "q=" << last_q << " expected " << expected_q << " (eps_p=" << last_eps_p << ")";
	EXPECT_GT(last_eps_p, 0.0) << "no plastic strain accumulated at 10% stretch";
}

/*=====================================================================
 * TEST 7: Pure elastic load-unload cycle — hypoelastic path-dependence
 *
 * Apply ELASTIC shear strain (below yield), then reverse it.  Ideal
 * hyperelasticity returns sigma=0.  Hughes-Winget hypoelastic has small
 * O(gamma^2) path-dependence even in elastic regime.
 *
 * Full-unload to gamma=0 is constructed to stay elastic in both directions
 * (|gamma| < gamma_yield = sigma_Y / (2*mu) so the path never touches the
 * yield surface).
 *===================================================================*/
TEST(ExplJ2, ElasticLoadUnloadPathDependence) {
	ExplJ2PlasticityT mat = MakeSteel();
	int nhist = mat.NumHistoryVars();
	std::vector<double> hist(nhist * MVSIZ, 0.0);
	hist[0*MVSIZ + 0] = 1.0; hist[4*MVSIZ + 0] = 1.0; hist[8*MVSIZ + 0] = 1.0;

	double F[9][MVSIZ];
	double sig[6][MVSIZ];

	double E = 200e3, nu = 0.3;
	double mu = E / (2.0*(1.0+nu));
	double gamma_yield = 250.0 / (mu * sqrt(3.0));  /* 0.00163 */
	double gamma_peak = 0.5 * gamma_yield;            /* well elastic */

	const int n_load = 100;
	double peak_sigma12 = 0.0;
	for (int step = 1; step <= n_load; step++) {
		SetIdentity(F, 0);
		F[1][0] = gamma_peak * step / n_load;
		mat.ComputeStress3D(1, F, sig, hist.data());
	}
	peak_sigma12 = sig[5][0];
	double peak_q = fabs(peak_sigma12) * sqrt(3.0);

	for (int step = 1; step <= n_load; step++) {
		SetIdentity(F, 0);
		F[1][0] = gamma_peak * (1.0 - (double)step / n_load);
		mat.ComputeStress3D(1, F, sig, hist.data());
	}

	/* residual: should be near zero for pure elastic cycle */
	double residual_sigma12 = sig[5][0];
	double residual_q = fabs(residual_sigma12) * sqrt(3.0);

	EXPECT_LT(residual_q, 0.05 * peak_q)
		<< "residual q=" << residual_q << " peak_q=" << peak_q
		<< " drift=" << (residual_q/peak_q)*100 << "% (H-W documented limit <5%)";

	/* no plastic strain should accumulate in purely elastic cycle */
	EXPECT_NEAR(hist[15*MVSIZ + 0], 0.0, 1e-10);
}
