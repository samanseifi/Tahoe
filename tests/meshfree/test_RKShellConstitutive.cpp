#include <gtest/gtest.h>

#include <cmath>

#include "FlanaganTaylor.h"
#include "PlaneStressJ2.h"

using namespace Tahoe::KLShell;

TEST(RKShellThickness, PadeUpdateHasExpectedProperties)
{
	EXPECT_DOUBLE_EQ(ThicknessStretchPade(0.0), 1.0);
	EXPECT_NEAR(ThicknessStretchPade(0.2)*ThicknessStretchPade(-0.2), 1.0, 1.0e-14);
	EXPECT_NEAR(ThicknessStretchPade(1.0e-3), std::exp(1.0e-3), 1.0e-9);
	EXPECT_DOUBLE_EQ(ThicknessStretchPade(2.0), 0.0);
}

TEST(RKShellPlaneStress, ElasticIncrementMatchesCondensedLaw)
{
	double sig[3] = {0.0, 0.0, 0.0};
	const double deps[3] = {1.0e-4, 0.0, 0.0};
	double ep = 0.0;
	const double E = 210000.0;
	const double nu = 0.3;

	PlaneStressJ2Return(sig, deps, ep, E, nu, 1.0e9, 0.0);

	const double c = E/(1.0 - nu*nu);
	EXPECT_NEAR(sig[0], c*deps[0], 1.0e-10);
	EXPECT_NEAR(sig[1], c*nu*deps[0], 1.0e-10);
	EXPECT_DOUBLE_EQ(sig[2], 0.0);
	EXPECT_DOUBLE_EQ(ep, 0.0);
}

TEST(RKShellCorotation, RigidQuarterTurnRotatesStressObjectively)
{
	const int increments = 200;
	const double dtheta = (0.5*std::acos(-1.0))/increments;
	double R[3][3] = {{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
	double V[3][3] = {{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
	const double dL[3][3] = {{0.0,-dtheta,0.0},{dtheta,0.0,0.0},{0.0,0.0,0.0}};
	for (int i = 0; i < increments; ++i) FlanaganTaylorStep(dL, R, V);

	double body_stress[3][3] = {{100.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
	double work[3][3];
	double spatial_stress[3][3];
	M3_mul(R, body_stress, work);
	M3_mulT(work, R, spatial_stress);

	EXPECT_NEAR(spatial_stress[0][0], 0.0, 1.0e-2);
	EXPECT_NEAR(spatial_stress[1][1], 100.0, 1.0e-2);
	EXPECT_NEAR(spatial_stress[0][1], 0.0, 1.0e-2);
}
