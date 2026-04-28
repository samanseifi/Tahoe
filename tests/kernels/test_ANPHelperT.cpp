/* test_ANPHelperT.cpp — unit tests for the Bonet-Burton ANP helper.
 *
 * Verifies the J-gather / nodal-average / J-bar-scatter pipeline on
 * minimal hand-checkable meshes.
 */

#include "gtest/gtest.h"
#include "ANPHelperT.h"
#include <cmath>
#include <vector>

using namespace Tahoe;

/*=====================================================================
 * Single tet, J=2 → J_bar = 2 (degenerate case).
 *===================================================================*/
TEST(ANP, SingleElementUniform) {
	int nelem = 1, nnod = 4, nen = 4;
	int conn[4] = {0, 1, 2, 3};
	double V_ref[1] = {1.0/6.0};
	ANPHelperT anp;
	anp.Init(nelem, nnod, nen, conn, V_ref);

	double J_e[1] = {2.0};
	double J_bar[1];
	anp.ComputeJBar(J_e, J_bar);

	EXPECT_NEAR(J_bar[0], 2.0, 1e-12);
}

/*=====================================================================
 * Two tets sharing a face, both J=1.5 → J_bar should be 1.5 everywhere.
 * (Volume-weighted average of equal values is just the value.)
 *===================================================================*/
TEST(ANP, TwoElementsUniformJ) {
	int nelem = 2, nnod = 5, nen = 4;
	int conn[8] = {0, 1, 2, 3,    /* tet 0: nodes 0..3 */
	               1, 4, 2, 3};   /* tet 1: shares face (1,2,3) */
	double V_ref[2] = {1.0/6.0, 1.0/6.0};
	ANPHelperT anp;
	anp.Init(nelem, nnod, nen, conn, V_ref);

	double J_e[2] = {1.5, 1.5};
	double J_bar[2];
	anp.ComputeJBar(J_e, J_bar);

	EXPECT_NEAR(J_bar[0], 1.5, 1e-12);
	EXPECT_NEAR(J_bar[1], 1.5, 1e-12);
}

/*=====================================================================
 * Two equal-volume tets sharing a 3-node face, J_e = (1.0, 2.0).
 *
 * Volume-weighted node J:
 *   node 0 (only in tet 0): J_n = (V*1.0)/V = 1.0
 *   node 4 (only in tet 1): J_n = (V*2.0)/V = 2.0
 *   nodes 1,2,3 (shared):    J_n = (V*1.0 + V*2.0)/(V+V) = 1.5
 *
 * J_bar per element:
 *   tet 0 = (1.0 + 1.5 + 1.5 + 1.5)/4 = 1.375
 *   tet 1 = (1.5 + 2.0 + 1.5 + 1.5)/4 = 1.625
 *===================================================================*/
TEST(ANP, TwoElementsAveraging) {
	int nelem = 2, nnod = 5, nen = 4;
	int conn[8] = {0, 1, 2, 3,
	               1, 4, 2, 3};
	double V_ref[2] = {1.0, 1.0};   /* equal weights */
	ANPHelperT anp;
	anp.Init(nelem, nnod, nen, conn, V_ref);

	double J_e[2] = {1.0, 2.0};
	double J_bar[2];
	anp.ComputeJBar(J_e, J_bar);

	EXPECT_NEAR(J_bar[0], 1.375, 1e-12);
	EXPECT_NEAR(J_bar[1], 1.625, 1e-12);
}

/*=====================================================================
 * Volume weighting: same shared face, but tet 0 has 9× the reference
 * volume of tet 1.
 *
 *   shared nodes (1,2,3) get J_n = (9*1.0 + 1*2.0)/(9+1) = 1.1
 *   node 0 only in tet 0: J_n = 1.0
 *   node 4 only in tet 1: J_n = 2.0
 *
 *   J_bar tet 0 = (1.0 + 1.1 + 1.1 + 1.1)/4 = 1.075
 *   J_bar tet 1 = (1.1 + 2.0 + 1.1 + 1.1)/4 = 1.325
 *===================================================================*/
TEST(ANP, TwoElementsVolumeWeighted) {
	int nelem = 2, nnod = 5, nen = 4;
	int conn[8] = {0, 1, 2, 3,
	               1, 4, 2, 3};
	double V_ref[2] = {9.0, 1.0};
	ANPHelperT anp;
	anp.Init(nelem, nnod, nen, conn, V_ref);

	double J_e[2] = {1.0, 2.0};
	double J_bar[2];
	anp.ComputeJBar(J_e, J_bar);

	EXPECT_NEAR(J_bar[0], 1.075, 1e-12);
	EXPECT_NEAR(J_bar[1], 1.325, 1e-12);
}

/*=====================================================================
 * F-bar transformation: F_bar = (J_bar/J)^(1/3) * F.
 * For F = diag(2, 1, 1) → J = 2, with J_bar = 1, scale = (1/2)^(1/3).
 *===================================================================*/
TEST(ANP, ApplyFBar3D) {
	const int stride = 8;  /* small SoA-style stride */
	int nelem = 2;

	double F[9 * stride];
	for (int k = 0; k < 9; k++)
		for (int i = 0; i < stride; i++)
			F[k*stride + i] = 0.0;
	/* element 0: F = diag(2,1,1) */
	F[0*stride + 0] = 2.0;  /* F11 */
	F[4*stride + 0] = 1.0;  /* F22 */
	F[8*stride + 0] = 1.0;  /* F33 */
	/* element 1: F = diag(1,1,3) */
	F[0*stride + 1] = 1.0;
	F[4*stride + 1] = 1.0;
	F[8*stride + 1] = 3.0;

	double J_e[2]    = {2.0, 3.0};
	double J_bar[2]  = {1.0, 1.0};   /* force isochoric */
	double F_bar[9 * stride];

	ANPHelperT::ApplyFBar3D(nelem, stride, F, J_e, J_bar, F_bar);

	double s0 = std::cbrt(0.5), s1 = std::cbrt(1.0/3.0);
	EXPECT_NEAR(F_bar[0*stride + 0], 2.0 * s0, 1e-12);
	EXPECT_NEAR(F_bar[4*stride + 0], 1.0 * s0, 1e-12);
	EXPECT_NEAR(F_bar[8*stride + 1], 3.0 * s1, 1e-12);
	/* det(F_bar) = J_bar (by construction) */
	double det0 = (F_bar[0*stride+0]) * (F_bar[4*stride+0]) * (F_bar[8*stride+0]);
	double det1 = (F_bar[0*stride+1]) * (F_bar[4*stride+1]) * (F_bar[8*stride+1]);
	EXPECT_NEAR(det0, 1.0, 1e-12);
	EXPECT_NEAR(det1, 1.0, 1e-12);
}
