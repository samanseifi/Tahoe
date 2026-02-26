#include <gtest/gtest.h>
#include "NeoHookean.h"
#include "dArrayT.h"
#include "dSymMatrixT.h"

using namespace Tahoe;

// MeanStress(J) = 0.5*kappa*(J^2 - 1)  (from PotentialT base)
// At identity: lambda_bar = {1,1,1}, J=1 -> energy=0, stress=0

TEST(NeoHookean, EnergyAtIdentity) {
    NeoHookean mat;
    mat.SetKappaMu(100.0, 80.0);
    dArrayT lambda_bar(3);
    lambda_bar = 1.0;
    EXPECT_NEAR(mat.Energy(lambda_bar, 1.0), 0.0, 1.0e-12);
}

TEST(NeoHookean, DevStressZeroAtIdentity) {
    NeoHookean mat;
    mat.SetKappaMu(100.0, 80.0);
    dArrayT lambda_bar(3), tau(3);
    lambda_bar = 1.0;
    mat.DevStress(lambda_bar, tau);
    for (int i = 0; i < 3; i++)
        EXPECT_NEAR(tau[i], 0.0, 1.0e-12);
}

// Deviatoric Kirchhoff stress eigenvalues must sum to zero
TEST(NeoHookean, DevStressIsDeviatoric) {
    NeoHookean mat;
    mat.SetKappaMu(100.0, 80.0);
    dArrayT lambda_bar(3), tau(3);
    // volume-preserving uniaxial stretch: lambda_bar = {2, 1/sqrt(2), 1/sqrt(2)}
    lambda_bar[0] = 2.0; lambda_bar[1] = 1.0/sqrt(2.0); lambda_bar[2] = 1.0/sqrt(2.0);
    mat.DevStress(lambda_bar, tau);
    EXPECT_NEAR(tau[0] + tau[1] + tau[2], 0.0, 1.0e-10);
}

// MeanStress = 0.5*kappa*(J^2 - 1)
TEST(NeoHookean, MeanStressAtJ2) {
    NeoHookean mat;
    double kappa = 100.0;
    mat.SetKappaMu(kappa, 80.0);
    double J = 2.0;
    EXPECT_NEAR(mat.MeanStress(J), 0.5 * kappa * (J*J - 1.0), 1.0e-10);
}

// MeanStress is zero at J=1 (no volumetric strain)
TEST(NeoHookean, MeanStressZeroAtJ1) {
    NeoHookean mat;
    mat.SetKappaMu(100.0, 80.0);
    EXPECT_NEAR(mat.MeanStress(1.0), 0.0, 1.0e-12);
}

// DevMod diagonal entries should be positive for positive mu
TEST(NeoHookean, DevModDiagonalPositive) {
    NeoHookean mat;
    mat.SetKappaMu(100.0, 80.0);
    dArrayT lambda_bar(3);
    lambda_bar[0] = 1.5; lambda_bar[1] = 0.9; lambda_bar[2] = 0.7;
    dSymMatrixT eigenmod(dSymMatrixT::k3D);
    mat.DevMod(lambda_bar, eigenmod);
    EXPECT_GT(eigenmod[0], 0.0);
    EXPECT_GT(eigenmod[1], 0.0);
    EXPECT_GT(eigenmod[2], 0.0);
}
