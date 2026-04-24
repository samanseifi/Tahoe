/* ExplJ2PlasticityT.h — batch finite-strain J2 plasticity for explicit dynamics.
 *
 * Hughes-Winget incremental-objective algorithm with isotropic linear hardening.
 * Yield criterion: phi = sqrt(3*J2(dev sigma)) - (sigma_Y + H * eps_p)
 * Radial return in deviatoric space; pressure updated elastically.
 *
 * History per IP (3D, 16 doubles):
 *   [ 0- 8] F_n      previous deformation gradient (row-major 3x3)
 *   [ 9-14] sigma_n  previous Cauchy stress (Voigt: s11,s22,s33,s23,s13,s12)
 *   [   15] eps_p    equivalent plastic strain
 *
 * History per IP (2D plane strain, 9 doubles):
 *   [ 0- 3] F_n      previous in-plane F (F11,F12,F21,F22)
 *   [ 4- 7] sigma_n  previous Cauchy stress (s11,s22,s33,s12)
 *   [   8] eps_p
 */

#ifndef _EXPL_J2_PLASTICITY_T_H_
#define _EXPL_J2_PLASTICITY_T_H_

#include "ExplicitMaterialT.h"

namespace Tahoe {

class ExplJ2PlasticityT : public ExplicitMaterialT
{
public:
	/** Construct with bulk/shear moduli, yield stress, linear hardening modulus.
	 *  \param mu       shear modulus
	 *  \param kappa    bulk modulus
	 *  \param sigma_Y  initial yield stress
	 *  \param H        isotropic hardening modulus (linear)
	 *  \param density  mass density */
	ExplJ2PlasticityT(double mu, double kappa,
	                  double sigma_Y, double H, double density);

	double Density(void) const { return fDensity; }
	double WaveSpeed(void) const;

	/** 3D: 16 history vars per IP. 2D not yet implemented (returns 9). */
	int NumHistoryVars(void) const { return 16; }

	/** Set F_n = I at t=0 for every IP-element. */
	void InitializeHistory(int nip, int total, double* history) const;

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
	double fLambda;  /* Lame lambda = kappa - 2*mu/3 */
	double fSigmaY;
	double fH;
	double fDensity;
};

} /* namespace Tahoe */
#endif /* _EXPL_J2_PLASTICITY_T_H_ */
