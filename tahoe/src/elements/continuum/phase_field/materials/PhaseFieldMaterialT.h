/* Phase-field fracture material */
#ifndef _PHASE_FIELD_MATERIAL_T_H_
#define _PHASE_FIELD_MATERIAL_T_H_

/* base class */
#include "ContinuumMaterialT.h"

/* direct members */
#include "dMatrixT.h"
#include "dArrayT.h"

namespace Tahoe {

/* forward declarations */
class PhaseFieldMatSupportT;

/** Interface for phase-field fracture materials.
 *  Stores the fracture toughness Gc, length scale ell, and residual
 *  stiffness parameter k_small. The phase-field equation (AT2 model) is:
 *
 *    -Gc*ell*Laplacian(d) + (Gc/ell)*d = 2*(1-d)*H
 *
 *  where H = max{psi_total} is the crack driving force (history variable). */
class PhaseFieldMaterialT: public ContinuumMaterialT
{
public:

	/** constructor */
	PhaseFieldMaterialT(void);

	/** set support */
	void SetPhaseFieldMatSupport(const PhaseFieldMatSupportT* support);

	/** \name material parameters */
	/*@{*/
	/** fracture toughness Gc [force/length] */
	double FractureToughness(void) const { return fGc; }

	/** length scale parameter ell [length] */
	double LengthScale(void) const { return fEll; }

	/** residual stiffness to prevent singularity */
	double ResidualStiffness(void) const { return fKSmall; }
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	virtual void DefineParameters(ParameterListT& list) const;
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** support for phase-field materials */
	const PhaseFieldMatSupportT* fPhaseFieldMatSupport;

	/** \name parameters */
	/*@{*/
	double fGc;     /**< fracture toughness */
	double fEll;    /**< length scale */
	double fKSmall; /**< residual stiffness (~1e-6) */
	/*@}*/
};

} // namespace Tahoe

#endif /* _PHASE_FIELD_MATERIAL_T_H_ */
