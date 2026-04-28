/* Hex8KernelT.cpp — 8-node hexahedral, 2x2x2 Gauss integration. */
#include "Hex8KernelT.h"

using namespace Tahoe;

static const double kGP = 0.5773502691896258; /* 1/sqrt(3) */

static const double gp_xi[8]  = { -kGP,+kGP,+kGP,-kGP, -kGP,+kGP,+kGP,-kGP };
static const double gp_eta[8] = { -kGP,-kGP,+kGP,+kGP, -kGP,-kGP,+kGP,+kGP };
static const double gp_mu[8]  = { -kGP,-kGP,-kGP,-kGP, +kGP,+kGP,+kGP,+kGP };

void Hex8KernelT::ComputeIPData(
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
	double xi  = gp_xi[ip];
	double eta = gp_eta[ip];
	double mu  = gp_mu[ip];
	weight = 1.0; /* 2x2x2 Gauss: all weights = 1 */

	/* parent domain shape derivatives */
	double ep = 1.0+eta, em = 1.0-eta;
	double mp = 1.0+mu,  mm = 1.0-mu;
	double xp = 1.0+xi,  xm = 1.0-xi;

	double dNdxi[8], dNdeta[8], dNdmu[8];
	dNdxi[0] = -0.125*em*mm;  dNdeta[0] = -0.125*xm*mm;  dNdmu[0] = -0.125*xm*em;
	dNdxi[1] = +0.125*em*mm;  dNdeta[1] = -0.125*xp*mm;  dNdmu[1] = -0.125*xp*em;
	dNdxi[2] = +0.125*ep*mm;  dNdeta[2] = +0.125*xp*mm;  dNdmu[2] = -0.125*xp*ep;
	dNdxi[3] = -0.125*ep*mm;  dNdeta[3] = +0.125*xm*mm;  dNdmu[3] = -0.125*xm*ep;
	dNdxi[4] = -0.125*em*mp;  dNdeta[4] = -0.125*xm*mp;  dNdmu[4] = +0.125*xm*em;
	dNdxi[5] = +0.125*em*mp;  dNdeta[5] = -0.125*xp*mp;  dNdmu[5] = +0.125*xp*em;
	dNdxi[6] = +0.125*ep*mp;  dNdeta[6] = +0.125*xp*mp;  dNdmu[6] = +0.125*xp*ep;
	dNdxi[7] = -0.125*ep*mp;  dNdeta[7] = +0.125*xm*mp;  dNdmu[7] = +0.125*xm*ep;

	/* Jacobian, inverse, spatial derivatives — vectorized over nel */
	for (int i = 0; i < nel; i++) {
		double J[3][3] = {{0.0}};
		for (int n = 0; n < 8; n++) {
			J[0][0] += xc[n][i]*dNdxi[n];  J[0][1] += xc[n][i]*dNdeta[n];  J[0][2] += xc[n][i]*dNdmu[n];
			J[1][0] += yc[n][i]*dNdxi[n];  J[1][1] += yc[n][i]*dNdeta[n];  J[1][2] += yc[n][i]*dNdmu[n];
			J[2][0] += zc[n][i]*dNdxi[n];  J[2][1] += zc[n][i]*dNdeta[n];  J[2][2] += zc[n][i]*dNdmu[n];
		}

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
		 * Note: uses COLUMN i of Ji (= J^{-1}), NOT row i.
		 * Bug fix Mar 2026: this was wrong (transposed) for non-axis-aligned
		 * Jacobians; symmetric Ji on aligned meshes hid the issue. */
		for (int n = 0; n < 8; n++) {
			dNdx[n][i] = Ji[0][0]*dNdxi[n] + Ji[1][0]*dNdeta[n] + Ji[2][0]*dNdmu[n];
			dNdy[n][i] = Ji[0][1]*dNdxi[n] + Ji[1][1]*dNdeta[n] + Ji[2][1]*dNdmu[n];
			dNdz[n][i] = Ji[0][2]*dNdxi[n] + Ji[1][2]*dNdeta[n] + Ji[2][2]*dNdmu[n];
		}
	}
}
