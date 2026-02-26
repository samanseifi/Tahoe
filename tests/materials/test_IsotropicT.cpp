#include <gtest/gtest.h>
#include "IsotropicT.h"
#include "ExceptionT.h"

using namespace Tahoe;

TEST(IsotropicT, DefaultConstructorYieldsZeroModuli) {
    IsotropicT mat;
    EXPECT_DOUBLE_EQ(mat.Young(), 0.0);
    EXPECT_DOUBLE_EQ(mat.Mu(),    0.0);
}

TEST(IsotropicT, ShearModulusFromENu) {
    IsotropicT mat;
    mat.Set_E_nu(200.0e9, 0.3);
    double expected = 200.0e9 / (2.0 * (1.0 + 0.3));
    EXPECT_NEAR(mat.Mu(), expected, 1.0);
}

TEST(IsotropicT, LameModulusFromENu) {
    IsotropicT mat;
    mat.Set_E_nu(200.0e9, 0.3);
    double mu  = 200.0e9 / (2.0 * 1.3);
    double lam = 2.0 * mu * 0.3 / (1.0 - 0.6);
    EXPECT_NEAR(mat.Lambda(), lam, 1.0);
}

TEST(IsotropicT, YoungFromMuKappa) {
    IsotropicT mat;
    double mu = 80.0e9, kappa = 200.0e9;
    mat.Set_mu_kappa(mu, kappa);
    // E = 9*kappa*mu / (3*kappa + mu)
    double E_expected = 9.0 * kappa * mu / (3.0 * kappa + mu);
    EXPECT_NEAR(mat.Young(), E_expected, 1.0e3);
}

TEST(IsotropicT, NegativeYoungThrows) {
    IsotropicT mat;
    EXPECT_THROW(mat.Set_E_nu(-1.0, 0.3), ExceptionT::CodeT);
}

TEST(IsotropicT, PoissonOutOfRangeThrows) {
    IsotropicT mat;
    EXPECT_THROW(mat.Set_E_nu(100.0, 0.6), ExceptionT::CodeT);
}
