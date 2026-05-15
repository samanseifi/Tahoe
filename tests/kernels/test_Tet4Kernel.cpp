/* $Id: test_Tet4Kernel.cpp,v 2.0 2026/05/08 samanseifi Exp $ */
/* test_Tet4Kernel.cpp — unit tests for Tet4KernelT.
 *
 * Verifies geometry-only properties of the 4-node constant-strain tet:
 *   - Reference unit tet (corners at standard positions): volume = 1/6
 *   - Spatial derivatives sum to zero (partition of unity)
 *   - Linear field reproduction: a linear function f(x,y,z) = a*x + b*y + c*z
 *     interpolated through nodal values has dN/dx . f_n = a, etc.
 *   - Batch consistency: identical coords at different MVSIZ slots give identical output
 */

#include "gtest/gtest.h"
#include "Tet4KernelT.h"
#include <cmath>

using namespace Tahoe;

static const int MVSIZ = ExplicitKernelT::MVSIZ;
static const int MAX_NEN = ExplicitKernelT::MAX_NEN;

/* helper: Tahoe Tet4 reference unit tet (see Tet4KernelT.h header) */
static void SetReferenceTet(double xc[MAX_NEN][MVSIZ],
                            double yc[MAX_NEN][MVSIZ],
                            double zc[MAX_NEN][MVSIZ], int slot)
{
	/* node 0 at (1,0,0), node 1 at (0,1,0), node 2 at (0,0,0), node 3 at (0,0,1) */
	xc[0][slot] = 1.0;  yc[0][slot] = 0.0;  zc[0][slot] = 0.0;
	xc[1][slot] = 0.0;  yc[1][slot] = 1.0;  zc[1][slot] = 0.0;
	xc[2][slot] = 0.0;  yc[2][slot] = 0.0;  zc[2][slot] = 0.0;
	xc[3][slot] = 0.0;  yc[3][slot] = 0.0;  zc[3][slot] = 1.0;
}

TEST(Tet4Kernel, ReferenceVolume) {
	Tet4KernelT k;
	double xc[MAX_NEN][MVSIZ], yc[MAX_NEN][MVSIZ], zc[MAX_NEN][MVSIZ];
	double dNdx[MAX_NEN][MVSIZ], dNdy[MAX_NEN][MVSIZ], dNdz[MAX_NEN][MVSIZ];
	double detJ[MVSIZ];
	double w;

	SetReferenceTet(xc, yc, zc, 0);
	k.ComputeIPData(0, 1, xc, yc, zc, dNdx, dNdy, dNdz, detJ, w);

	/* For Tet4: V_element = w * detJ, w = 1/6, detJ = 1 for the reference unit tet
	 *   (Jacobian columns are the position vectors of nodes 0,1,3 from node 2,
	 *    here = e_x, e_y, e_z, so |det J| = 1)
	 * Total volume of unit tet = 1/6.  */
	EXPECT_NEAR(w * detJ[0], 1.0/6.0, 1e-12)
		<< "volume mismatch: w*detJ=" << w*detJ[0];
}

TEST(Tet4Kernel, PartitionOfUnityDerivatives) {
	/* Derivatives of shape functions sum to zero in each spatial direction
	 *   (because sum_n N_n = 1 always, so derivative is identically 0).  */
	Tet4KernelT k;
	double xc[MAX_NEN][MVSIZ], yc[MAX_NEN][MVSIZ], zc[MAX_NEN][MVSIZ];
	double dNdx[MAX_NEN][MVSIZ], dNdy[MAX_NEN][MVSIZ], dNdz[MAX_NEN][MVSIZ];
	double detJ[MVSIZ];
	double w;

	SetReferenceTet(xc, yc, zc, 0);
	k.ComputeIPData(0, 1, xc, yc, zc, dNdx, dNdy, dNdz, detJ, w);

	double sx = 0, sy = 0, sz = 0;
	for (int n = 0; n < 4; n++) {
		sx += dNdx[n][0];
		sy += dNdy[n][0];
		sz += dNdz[n][0];
	}
	EXPECT_NEAR(sx, 0.0, 1e-12);
	EXPECT_NEAR(sy, 0.0, 1e-12);
	EXPECT_NEAR(sz, 0.0, 1e-12);
}

TEST(Tet4Kernel, LinearFieldReproduction) {
	/* Linear field f(x) = a + b*x + c*y + d*z, sampled at the 4 nodes:
	 *   sum_n f_n * dN_n/dx = b      (the gradient)
	 * Tests that dN/dx is a correct gradient operator.
	 */
	Tet4KernelT k;
	double xc[MAX_NEN][MVSIZ], yc[MAX_NEN][MVSIZ], zc[MAX_NEN][MVSIZ];
	double dNdx[MAX_NEN][MVSIZ], dNdy[MAX_NEN][MVSIZ], dNdz[MAX_NEN][MVSIZ];
	double detJ[MVSIZ];
	double w;

	SetReferenceTet(xc, yc, zc, 0);
	k.ComputeIPData(0, 1, xc, yc, zc, dNdx, dNdy, dNdz, detJ, w);

	/* set field values: f(x,y,z) = 7 + 3*x - 2*y + 5*z */
	double a = 7.0, b = 3.0, c = -2.0, d = 5.0;
	double f[4];
	for (int n = 0; n < 4; n++)
		f[n] = a + b*xc[n][0] + c*yc[n][0] + d*zc[n][0];

	double gx = 0, gy = 0, gz = 0;
	for (int n = 0; n < 4; n++) {
		gx += f[n] * dNdx[n][0];
		gy += f[n] * dNdy[n][0];
		gz += f[n] * dNdz[n][0];
	}
	EXPECT_NEAR(gx, b, 1e-12);
	EXPECT_NEAR(gy, c, 1e-12);
	EXPECT_NEAR(gz, d, 1e-12);
}

TEST(Tet4Kernel, ScaledTetVolume) {
	/* Scale reference tet by alpha → volume scales by alpha^3 */
	Tet4KernelT k;
	double xc[MAX_NEN][MVSIZ], yc[MAX_NEN][MVSIZ], zc[MAX_NEN][MVSIZ];
	double dNdx[MAX_NEN][MVSIZ], dNdy[MAX_NEN][MVSIZ], dNdz[MAX_NEN][MVSIZ];
	double detJ[MVSIZ];
	double w;

	double alpha = 2.5;
	xc[0][0] = alpha;  yc[0][0] = 0.0;    zc[0][0] = 0.0;
	xc[1][0] = 0.0;    yc[1][0] = alpha;  zc[1][0] = 0.0;
	xc[2][0] = 0.0;    yc[2][0] = 0.0;    zc[2][0] = 0.0;
	xc[3][0] = 0.0;    yc[3][0] = 0.0;    zc[3][0] = alpha;

	k.ComputeIPData(0, 1, xc, yc, zc, dNdx, dNdy, dNdz, detJ, w);

	double V = w * detJ[0];
	double V_expected = alpha*alpha*alpha / 6.0;
	EXPECT_NEAR(V, V_expected, 1e-10)
		<< "scaled tet V=" << V << ", expected " << V_expected;
}

TEST(Tet4Kernel, BatchConsistency) {
	/* Identical coords at different MVSIZ slots produce identical output. */
	Tet4KernelT k;
	double xc[MAX_NEN][MVSIZ], yc[MAX_NEN][MVSIZ], zc[MAX_NEN][MVSIZ];
	double dNdx[MAX_NEN][MVSIZ], dNdy[MAX_NEN][MVSIZ], dNdz[MAX_NEN][MVSIZ];
	double detJ[MVSIZ];
	double w;

	int nel = 32;
	for (int i = 0; i < nel; i++) SetReferenceTet(xc, yc, zc, i);
	k.ComputeIPData(0, nel, xc, yc, zc, dNdx, dNdy, dNdz, detJ, w);

	for (int i = 1; i < nel; i++) {
		EXPECT_NEAR(detJ[i], detJ[0], 1e-12) << "slot " << i;
		for (int n = 0; n < 4; n++) {
			EXPECT_NEAR(dNdx[n][i], dNdx[n][0], 1e-12);
			EXPECT_NEAR(dNdy[n][i], dNdy[n][0], 1e-12);
			EXPECT_NEAR(dNdz[n][i], dNdz[n][0], 1e-12);
		}
	}
}
