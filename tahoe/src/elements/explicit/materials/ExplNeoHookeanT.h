/* $Id: ExplNeoHookeanT.h,v 2.0 2026/05/08 samanseifi Exp $ */
/* ExplNeoHookeanT.h — batch Neo-Hookean material for explicit dynamics.
 *
 * Compressible Neo-Hookean in current configuration (updated Lagrangian).
 * Cauchy stress: sigma = (mu/J)(b - I) + (kappa*(J-1)/J)*I
 * where b = F*F^T is the left Cauchy-Green tensor, J = det(F).
 *
 * Hyperelastic, no history variables.
 */

#ifndef _EXPL_NEO_HOOKEAN_T_H_
#define _EXPL_NEO_HOOKEAN_T_H_

#include "ExplicitMaterialT.h"

namespace Tahoe {

class ExplNeoHookeanT : public ExplicitMaterialT
{
public:
	ExplNeoHookeanT(double mu, double kappa, double density);

	double Density(void) const { return fDensity; }
	int NumHistoryVars(void) const { return 0; }
	double WaveSpeed(void) const;

	void ComputeStress2D(
		int nel,
		const double F[][MVSIZ],
		double sig11[], double sig22[], double sig12[],
		double* history
	) const;

	void ComputeStress3D(
		int nel,
		const double F[][MVSIZ],
		double sig[][MVSIZ],
		double* history
	) const;

private:
	double fMu;
	double fKappa;
	double fDensity;
};

} /* namespace Tahoe */
#endif /* _EXPL_NEO_HOOKEAN_T_H_ */
