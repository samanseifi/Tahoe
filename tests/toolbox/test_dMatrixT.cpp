#include <gtest/gtest.h>
#include "dMatrixT.h"

using namespace Tahoe;

TEST(dMatrixT, Det2x2) {
    dMatrixT A(2);
    A(0,0) = 4.0; A(0,1) = 3.0;
    A(1,0) = 2.0; A(1,1) = 1.0;
    // det = 4*1 - 3*2 = -2
    EXPECT_DOUBLE_EQ(A.Det(), -2.0);
}

TEST(dMatrixT, Det3x3Identity) {
    dMatrixT I(3);
    I = 0.0;
    I(0,0) = 1.0; I(1,1) = 1.0; I(2,2) = 1.0;
    EXPECT_DOUBLE_EQ(I.Det(), 1.0);
}

TEST(dMatrixT, Trace) {
    dMatrixT A(3);
    A = 0.0;
    A(0,0) = 2.0; A(1,1) = 3.0; A(2,2) = 5.0;
    EXPECT_DOUBLE_EQ(A.Trace(), 10.0);
}

TEST(dMatrixT, Inverse2x2) {
    dMatrixT A(2), Ainv(2);
    A(0,0) = 2.0; A(0,1) = 1.0;
    A(1,0) = 5.0; A(1,1) = 3.0;
    Ainv.Inverse(A);
    // Verify A * Ainv = I
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++) {
            double sum = 0.0;
            for (int k = 0; k < 2; k++) sum += A(i,k) * Ainv(k,j);
            EXPECT_NEAR(sum, (i == j) ? 1.0 : 0.0, 1.0e-12);
        }
}

TEST(dMatrixT, Symmetrize) {
    dMatrixT A(2);
    A(0,0) = 1.0; A(0,1) = 2.0;
    A(1,0) = 4.0; A(1,1) = 3.0;
    A.Symmetrize();
    EXPECT_DOUBLE_EQ(A(0,1), A(1,0));
    EXPECT_DOUBLE_EQ(A(0,1), 3.0);  // (2+4)/2
}

TEST(dMatrixT, ScalarAssignment) {
    dMatrixT A(3);
    A = 0.0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            EXPECT_DOUBLE_EQ(A(i,j), 0.0);
}
