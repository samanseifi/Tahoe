#include <gtest/gtest.h>
#include "dSymMatrixT.h"
#include "dArrayT.h"
#include <algorithm>
#include <functional>

using namespace Tahoe;

TEST(dSymMatrixT, Dimension3DStorageSize) {
    dSymMatrixT S(dSymMatrixT::k3D);
    EXPECT_EQ(S.Length(), 6);
}

TEST(dSymMatrixT, Dimension2DStorageSize) {
    dSymMatrixT S(dSymMatrixT::k2D);
    EXPECT_EQ(S.Length(), 3);
}

TEST(dSymMatrixT, Trace3D) {
    dSymMatrixT S(dSymMatrixT::k3D);
    S = 0.0;
    S(0,0) = 1.0; S(1,1) = 2.0; S(2,2) = 3.0;
    EXPECT_DOUBLE_EQ(S.Trace(), 6.0);
}

TEST(dSymMatrixT, Det3DDiagonal) {
    dSymMatrixT S(dSymMatrixT::k3D);
    S = 0.0;
    S(0,0) = 2.0; S(1,1) = 3.0; S(2,2) = 4.0;
    EXPECT_NEAR(S.Det(), 24.0, 1.0e-12);
}

TEST(dSymMatrixT, DeviatoricTraceZero) {
    dSymMatrixT S(dSymMatrixT::k3D), dev(dSymMatrixT::k3D);
    S = 0.0;
    S(0,0) = 4.0; S(1,1) = 2.0; S(2,2) = 0.0;
    dev.Deviatoric(S);
    EXPECT_NEAR(dev.Trace(), 0.0, 1.0e-12);
}

// Note: dSymMatrixT::Eigenvalues ignores sort_descending for 3D Cardano path
// (#pragma unused). We sort ourselves before comparing.
TEST(dSymMatrixT, Eigenvalues3DDiagonal) {
    dSymMatrixT S(dSymMatrixT::k3D);
    S = 0.0;
    S(0,0) = 5.0; S(1,1) = 3.0; S(2,2) = 1.0;
    dArrayT evals(3);
    S.Eigenvalues(evals, true);
    double e[3] = {evals[0], evals[1], evals[2]};
    std::sort(e, e + 3, std::greater<double>());
    EXPECT_NEAR(e[0], 5.0, 1.0e-10);
    EXPECT_NEAR(e[1], 3.0, 1.0e-10);
    EXPECT_NEAR(e[2], 1.0, 1.0e-10);
}
