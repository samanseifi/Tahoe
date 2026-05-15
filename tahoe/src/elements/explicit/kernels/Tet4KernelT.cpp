/* $Id: Tet4KernelT.cpp,v 2.0 2026/05/08 samanseifi Exp $ */
/* Tet4KernelT.cpp — 4-node constant-strain tet, 1-point integration.
 *
 * Parametric coords (r,s,t) on the reference unit tet:
 *   node 0:  (1, 0, 0)
 *   node 1:  (0, 1, 0)
 *   node 2:  (0, 0, 0)
 *   node 3:  (0, 0, 1)
 *
 * Shape functions: N_0=r, N_1=s, N_2=1-r-s-t, N_3=t
 * Parametric derivatives are constants (no Gauss loop):
 *
 *   dN/dr  = [ +1,  0, -1,  0 ]
 *   dN/ds  = [  0, +1, -1,  0 ]
 *   dN/dt  = [  0,  0, -1, +1 ]
 *
 * Single IP at centroid; weight = 1/6 (volume of unit reference tet).
 * The Jacobian is constant within each element.
 */
#include "Tet4KernelT.h"

using namespace Tahoe;

/* parametric derivatives (Tahoe node order — see header) */
static const double dNdr[4] = { +1.0,  0.0, -1.0,  0.0 };
static const double dNds[4] = {  0.0, +1.0, -1.0,  0.0 };
static const double dNdt[4] = {  0.0,  0.0, -1.0, +1.0 };

void Tet4KernelT::ComputeIPData(
	int ip, int nel,
	const double xc[][MVSIZ],
	const double yc[][MVSIZ],
	const double zc[][MVSIZ],
	double dNdx[][MVSIZ],
	double dNdy[][MVSIZ],
	double dNdz[][MVSIZ],
	double detJ[],
	double& weight) const
{
	weight = 1.0 / 6.0;   /* unit reference tet volume */

	/* Jacobian, inverse, spatial derivatives — vectorized over nel.
	 * Tet4 has constant parametric derivatives, so the Jacobian
	 * computation simplifies: J_ij = sum_n (x_n)_i * dN_n/dxi_j */
	for (int i = 0; i < nel; i++) {
		double J[3][3];
		J[0][0] = xc[0][i]*dNdr[0] + xc[1][i]*dNdr[1] + xc[2][i]*dNdr[2] + xc[3][i]*dNdr[3];
		J[0][1] = xc[0][i]*dNds[0] + xc[1][i]*dNds[1] + xc[2][i]*dNds[2] + xc[3][i]*dNds[3];
		J[0][2] = xc[0][i]*dNdt[0] + xc[1][i]*dNdt[1] + xc[2][i]*dNdt[2] + xc[3][i]*dNdt[3];
		J[1][0] = yc[0][i]*dNdr[0] + yc[1][i]*dNdr[1] + yc[2][i]*dNdr[2] + yc[3][i]*dNdr[3];
		J[1][1] = yc[0][i]*dNds[0] + yc[1][i]*dNds[1] + yc[2][i]*dNds[2] + yc[3][i]*dNds[3];
		J[1][2] = yc[0][i]*dNdt[0] + yc[1][i]*dNdt[1] + yc[2][i]*dNdt[2] + yc[3][i]*dNdt[3];
		J[2][0] = zc[0][i]*dNdr[0] + zc[1][i]*dNdr[1] + zc[2][i]*dNdr[2] + zc[3][i]*dNdr[3];
		J[2][1] = zc[0][i]*dNds[0] + zc[1][i]*dNds[1] + zc[2][i]*dNds[2] + zc[3][i]*dNds[3];
		J[2][2] = zc[0][i]*dNdt[0] + zc[1][i]*dNdt[1] + zc[2][i]*dNdt[2] + zc[3][i]*dNdt[3];

		double det = J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1])
		           - J[0][1]*(J[1][0]*J[2][2] - J[1][2]*J[2][0])
		           + J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
		detJ[i] = det;
		double inv = 1.0 / det;

		double Ji[3][3];
		Ji[0][0] =  (J[1][1]*J[2][2] - J[1][2]*J[2][1])*inv;
		Ji[0][1] = -(J[0][1]*J[2][2] - J[0][2]*J[2][1])*inv;
		Ji[0][2] =  (J[0][1]*J[1][2] - J[0][2]*J[1][1])*inv;
		Ji[1][0] = -(J[1][0]*J[2][2] - J[1][2]*J[2][0])*inv;
		Ji[1][1] =  (J[0][0]*J[2][2] - J[0][2]*J[2][0])*inv;
		Ji[1][2] = -(J[0][0]*J[1][2] - J[0][2]*J[1][0])*inv;
		Ji[2][0] =  (J[1][0]*J[2][1] - J[1][1]*J[2][0])*inv;
		Ji[2][1] = -(J[0][0]*J[2][1] - J[0][1]*J[2][0])*inv;
		Ji[2][2] =  (J[0][0]*J[1][1] - J[0][1]*J[1][0])*inv;

		/* dN_n/dX_i = sum_k (J^{-1})[k,i] * dN_n/dξ_k    (chain rule)
		 * Note: uses COLUMN i of Ji, NOT row i. */
		for (int n = 0; n < 4; n++) {
			dNdx[n][i] = Ji[0][0]*dNdr[n] + Ji[1][0]*dNds[n] + Ji[2][0]*dNdt[n];
			dNdy[n][i] = Ji[0][1]*dNdr[n] + Ji[1][1]*dNds[n] + Ji[2][1]*dNdt[n];
			dNdz[n][i] = Ji[0][2]*dNdr[n] + Ji[1][2]*dNds[n] + Ji[2][2]*dNdt[n];
		}
	}
}
