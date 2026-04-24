/* ExplicitMaterialT.h — abstract batch material interface for explicit dynamics.
 *
 * A batch material computes Cauchy stress from the deformation gradient F
 * for MVSIZ elements simultaneously. The element driver is responsible for
 * computing F; the material only implements the constitutive law.
 *
 * All arrays use SoA layout for SIMD auto-vectorization. The virtual dispatch
 * happens once per IP per element group (not per element), so the v-table cost
 * is negligible — the inner loop inside each implementation is what gets
 * vectorized.
 *
 * Materials with history variables (plasticity) store them in a flat array
 * with nhist doubles per element per IP, managed by the element driver.
 */

#ifndef _EXPLICIT_MATERIAL_T_H_
#define _EXPLICIT_MATERIAL_T_H_

#include "ExplicitKernelT.h" /* for MVSIZ */

namespace Tahoe {

class ExplicitMaterialT
{
public:

	static const int MVSIZ = ExplicitKernelT::MVSIZ;

	virtual ~ExplicitMaterialT(void) {}

	/** material density */
	virtual double Density(void) const = 0;

	/** number of history variables per element per IP (0 for elastic) */
	virtual int NumHistoryVars(void) const = 0;

	/** P-wave speed: c = sqrt((kappa + 4mu/3) / rho).
	 *  Used for CFL time step and viscous hourglass control. */
	virtual double WaveSpeed(void) const = 0;

	/** Initialize history variables at t=0. Default: leave zero (correct for
	 *  most rate-independent hyperelastic and small-strain plasticity).
	 *  Override for rate plasticity where F_n must start as I.
	 *
	 *  \param nip      number of integration points per element
	 *  \param total    total number of elements
	 *  \param history  flat buffer of size (nip*NumHistoryVars()*total), zero-initialized
	 *                  Layout: history[ip * nhist * total + var * total + elem] */
	virtual void InitializeHistory(int nip, int total, double* history) const {}

	/** Compute Cauchy stress from deformation gradient F for 2D plane strain.
	 *
	 * \param nel     number of elements (<= MVSIZ)
	 * \param F       deformation gradient [4][MVSIZ]: F11,F12,F21,F22
	 * \param sig11   output Cauchy stress [MVSIZ]
	 * \param sig22   output Cauchy stress [MVSIZ]
	 * \param sig12   output Cauchy stress [MVSIZ]
	 * \param history per-element history at this IP (NULL if stateless)
	 *                layout: history[var * MVSIZ + elem] for SoA vectorization
	 */
	virtual void ComputeStress2D(
		int nel,
		const double F[][MVSIZ],
		double sig11[], double sig22[], double sig12[],
		double* history
	) const = 0;

	/** Compute Cauchy stress from deformation gradient F for 3D.
	 *
	 * \param F    deformation gradient [9][MVSIZ]: F11,F12,F13,F21,...,F33
	 * \param sig  output Cauchy stress [6][MVSIZ]: s11,s22,s33,s23,s13,s12
	 */
	virtual void ComputeStress3D(
		int nel,
		const double F[][MVSIZ],
		double sig[][MVSIZ],
		double* history
	) const = 0;
};

} /* namespace Tahoe */
#endif /* _EXPLICIT_MATERIAL_T_H_ */
