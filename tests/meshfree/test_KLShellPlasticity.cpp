/* test_KLShellPlasticity.cpp — plane-stress J2 radial-return tests for the meshfree KL shell.
 *
 * Validates PlaneStressJ2Return (development/src/elements/meshfree_kl_shell/PlaneStressJ2.h),
 * the per-Gauss-point plane-stress (sigma33=0) von Mises stress update used on the explicit
 * elasto-plastic track (Fig 18). Pure material-point algebra; no meshfree machinery needed.
 *
 *   - elastic step below yield reproduces the plane-stress moduli;
 *   - past yield the stress returns onto the von Mises surface (seq == Y) and plastic strain grows;
 *   - linear isotropic hardening: the surface follows Y0 + H*ep;
 *   - uniaxial and pure-shear yield match the analytic von Mises limits;
 *   - the shear block stays decoupled from the in-plane normal block.
 */
#include "gtest/gtest.h"
#include <cmath>
#include "PlaneStressJ2.h"

using namespace Tahoe::KLShell;

namespace {

const double E = 70.0e3, nu = 0.3;   /* aluminium-like (MPa) */
const double Y0 = 250.0;

/* von Mises effective stress sqrt(s^T P s) */
double Seq(const double s[3]) { return std::sqrt(J2pq(s)); }

TEST(KLShellPlasticity, ElasticBelowYield)
{
	double sig[3] = {0,0,0}, ep = 0.0;
	/* small uniaxial strain, stress stays well under yield */
	double deps[3] = {1.0e-4, 0.0, 0.0};
	PlaneStressJ2Return(sig, deps, ep, E, nu, Y0, 0.0);
	double c = E/(1.0-nu*nu);
	EXPECT_NEAR(sig[0], c*deps[0], 1e-8);          /* s11 = E/(1-nu^2) * e11 */
	EXPECT_NEAR(sig[1], c*nu*deps[0], 1e-8);       /* s22 = nu*E/(1-nu^2) * e11 */
	EXPECT_NEAR(sig[2], 0.0, 1e-12);
	EXPECT_EQ(ep, 0.0);                            /* no plastic flow */
}

TEST(KLShellPlasticity, UniaxialStrainPerfectPlasticityReturnsToSurface)
{
	double sig[3] = {0,0,0}, ep = 0.0;
	/* uniaxial STRAIN far past yield (lateral constrained: e22=0), perfect plasticity (H=0) */
	double deps[3] = {0.05, 0.0, 0.0};
	PlaneStressJ2Return(sig, deps, ep, E, nu, Y0, 0.0);
	EXPECT_NEAR(Seq(sig), Y0, 1e-6);               /* lands on the von Mises surface */
	EXPECT_GT(ep, 0.0);                            /* plastic strain accumulated */
	/* constrained lateral strain => Poisson-induced s22>0; surface => s11 > seq = Y0 */
	EXPECT_GT(sig[1], 0.0);
	EXPECT_GT(sig[0], sig[1]);
	EXPECT_GT(sig[0], Y0);
	EXPECT_NEAR(sig[2], 0.0, 1e-6);                /* no shear */
}

TEST(KLShellPlasticity, ConsistencyOnYieldSurface)
{
	double sig[3] = {0,0,0}, ep = 0.0;
	double H = 500.0;
	double deps[3] = {0.02, -0.005, 0.003};        /* general multiaxial increment past yield */
	PlaneStressJ2Return(sig, deps, ep, E, nu, Y0, H);
	EXPECT_NEAR(Seq(sig), Y0 + H*ep, 1e-4);        /* stress sits on Y(ep) exactly */
	EXPECT_GT(ep, 0.0);
}

TEST(KLShellPlasticity, LinearHardeningGrowsTheSurface)
{
	double H = 1000.0;
	double sig[3] = {0,0,0}, ep = 0.0;
	double d1[3] = {0.01,0,0};
	PlaneStressJ2Return(sig, d1, ep, E, nu, Y0, H);
	double seq1 = Seq(sig), ep1 = ep;
	double d2[3] = {0.01,0,0};
	PlaneStressJ2Return(sig, d2, ep, E, nu, Y0, H);
	double seq2 = Seq(sig);
	EXPECT_GT(ep, ep1);                            /* more plastic strain */
	EXPECT_GT(seq2, seq1);                         /* hardened: surface grew */
	EXPECT_NEAR(seq2, Y0 + H*ep, 1e-4);
}

TEST(KLShellPlasticity, PureShearYield)
{
	double sig[3] = {0,0,0}, ep = 0.0;
	double deps[3] = {0.0, 0.0, 0.02};             /* engineering shear past yield */
	PlaneStressJ2Return(sig, deps, ep, E, nu, Y0, 0.0);
	EXPECT_NEAR(Seq(sig), Y0, 1e-6);
	/* pure shear von Mises: seq = sqrt(3)*|s12| -> |s12| = Y0/sqrt(3) */
	EXPECT_NEAR(std::fabs(sig[2]), Y0/std::sqrt(3.0), 1.0);
	EXPECT_NEAR(sig[0], 0.0, 1e-6);
	EXPECT_NEAR(sig[1], 0.0, 1e-6);
}

} /* namespace */
