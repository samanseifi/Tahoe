#include <gtest/gtest.h>
#include "QuadT.h"
#include "dArrayT.h"

using namespace Tahoe;

// Node coordinates (from QuadT.cpp): ra={-1,1,1,-1}, sa={-1,-1,1,1}

TEST(QuadT, PartitionOfUnityAtCenter) {
    QuadT quad(4);
    dArrayT xi(2), Na(4);
    xi[0] = 0.0; xi[1] = 0.0;
    quad.EvaluateShapeFunctions(xi, Na);
    double sum = 0.0;
    for (int i = 0; i < 4; i++) sum += Na[i];
    EXPECT_NEAR(sum, 1.0, 1.0e-14);
}

TEST(QuadT, PartitionOfUnityAtGaussPoint) {
    QuadT quad(4);
    dArrayT xi(2), Na(4);
    const double g = 1.0 / sqrt(3.0);
    xi[0] = g; xi[1] = g;
    quad.EvaluateShapeFunctions(xi, Na);
    double sum = 0.0;
    for (int i = 0; i < 4; i++) sum += Na[i];
    EXPECT_NEAR(sum, 1.0, 1.0e-14);
}

// Node 0 is at (-1,-1): N_0 = 1, others = 0
TEST(QuadT, NodalRecoveryNode0) {
    QuadT quad(4);
    dArrayT xi(2), Na(4);
    xi[0] = -1.0; xi[1] = -1.0;
    quad.EvaluateShapeFunctions(xi, Na);
    EXPECT_NEAR(Na[0], 1.0, 1.0e-14);
    EXPECT_NEAR(Na[1], 0.0, 1.0e-14);
    EXPECT_NEAR(Na[2], 0.0, 1.0e-14);
    EXPECT_NEAR(Na[3], 0.0, 1.0e-14);
}

// Node 1 is at (1,-1): N_1 = 1, others = 0
TEST(QuadT, NodalRecoveryNode1) {
    QuadT quad(4);
    dArrayT xi(2), Na(4);
    xi[0] = 1.0; xi[1] = -1.0;
    quad.EvaluateShapeFunctions(xi, Na);
    EXPECT_NEAR(Na[0], 0.0, 1.0e-14);
    EXPECT_NEAR(Na[1], 1.0, 1.0e-14);
    EXPECT_NEAR(Na[2], 0.0, 1.0e-14);
    EXPECT_NEAR(Na[3], 0.0, 1.0e-14);
}

// All shape functions should be non-negative inside the element
TEST(QuadT, ShapeFunctionsNonNegativeInterior) {
    QuadT quad(4);
    dArrayT xi(2), Na(4);
    xi[0] = 0.3; xi[1] = -0.5;
    quad.EvaluateShapeFunctions(xi, Na);
    for (int i = 0; i < 4; i++)
        EXPECT_GE(Na[i], 0.0);
}
