/* Q4KernelT.cpp — 4-node quadrilateral, 2x2 Gauss integration. */
#include "Q4KernelT.h"

using namespace Tahoe;

static const double kGP = 0.5773502691896258; /* 1/sqrt(3) */

static const double gp_xi[4]  = { -kGP, +kGP, +kGP, -kGP };
static const double gp_eta[4] = { -kGP, -kGP, +kGP, +kGP };

void Q4KernelT::ComputeIPData(
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
	weight = 1.0; /* 2x2 Gauss: all weights = 1 */

	/* parent domain shape derivatives */
	double dNdxi[4], dNdeta[4];
	dNdxi[0] = -0.25*(1.0-eta);  dNdeta[0] = -0.25*(1.0-xi);
	dNdxi[1] = +0.25*(1.0-eta);  dNdeta[1] = -0.25*(1.0+xi);
	dNdxi[2] = +0.25*(1.0+eta);  dNdeta[2] = +0.25*(1.0+xi);
	dNdxi[3] = -0.25*(1.0+eta);  dNdeta[3] = +0.25*(1.0-xi);

	/* Jacobian, inverse, spatial derivatives — vectorized over nel */
	for (int i = 0; i < nel; i++) {
		double J11 = 0.0, J12 = 0.0, J21 = 0.0, J22 = 0.0;
		for (int n = 0; n < 4; n++) {
			J11 += xc[n][i]*dNdxi[n];   J12 += xc[n][i]*dNdeta[n];
			J21 += yc[n][i]*dNdxi[n];   J22 += yc[n][i]*dNdeta[n];
		}
		double det = J11*J22 - J12*J21;
		detJ[i] = det;
		double inv = 1.0 / det;
		double Ji11 =  J22*inv, Ji12 = -J12*inv;
		double Ji21 = -J21*inv, Ji22 =  J11*inv;
		for (int n = 0; n < 4; n++) {
			dNdx[n][i] = Ji11*dNdxi[n] + Ji12*dNdeta[n];
			dNdy[n][i] = Ji21*dNdxi[n] + Ji22*dNdeta[n];
		}
	}
}
