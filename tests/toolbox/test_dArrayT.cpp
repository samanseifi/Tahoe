#include <gtest/gtest.h>
#include "dArrayT.h"

using namespace Tahoe;

TEST(dArrayT, DefaultConstructorIsEmpty) {
    dArrayT a;
    EXPECT_EQ(a.Length(), 0);
}

TEST(dArrayT, ExplicitDimension) {
    dArrayT a(5);
    EXPECT_EQ(a.Length(), 5);
}

TEST(dArrayT, ScalarAssignment) {
    dArrayT a(3);
    a = 3.14;
    EXPECT_DOUBLE_EQ(a[0], 3.14);
    EXPECT_DOUBLE_EQ(a[1], 3.14);
    EXPECT_DOUBLE_EQ(a[2], 3.14);
}

TEST(dArrayT, Magnitude) {
    dArrayT a(3);
    a[0] = 3.0; a[1] = 0.0; a[2] = 4.0;
    EXPECT_DOUBLE_EQ(a.Magnitude(), 5.0);
}

TEST(dArrayT, OperatorPlusEquals) {
    dArrayT a(2), b(2);
    a[0] = 1.0; a[1] = 2.0;
    b[0] = 0.5; b[1] = -1.0;
    a += b;
    EXPECT_DOUBLE_EQ(a[0], 1.5);
    EXPECT_DOUBLE_EQ(a[1], 1.0);
}

TEST(dArrayT, ScalarMultiply) {
    dArrayT a(2);
    a[0] = 2.0; a[1] = -3.0;
    a *= 2.0;
    EXPECT_DOUBLE_EQ(a[0],  4.0);
    EXPECT_DOUBLE_EQ(a[1], -6.0);
}

TEST(dArrayT, CopyAssignment) {
    dArrayT a(3), b;
    a[0] = 1.0; a[1] = 2.0; a[2] = 3.0;
    b = a;
    EXPECT_EQ(b.Length(), 3);
    EXPECT_DOUBLE_EQ(b[0], 1.0);
    EXPECT_DOUBLE_EQ(b[2], 3.0);
}
