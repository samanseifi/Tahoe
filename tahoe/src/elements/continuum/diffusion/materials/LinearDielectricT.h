/* created for electromechanical coupling in DiffusionElementT */
#ifndef _LINEAR_DIELECTRIC_T_H_
#define _LINEAR_DIELECTRIC_T_H_

/* base class */
#include "DiffusionMaterialT.h"

namespace Tahoe {

/** Isotropic linear dielectric material for use as the voltage field material
 *  in electromechanical coupling. The only physical parameter is the electric
 *  permittivity epsilon, which is passed to both DiffusionElementT (for the
 *  reference-frame electric displacement residual) and to SimoQ1P0 (for the
 *  Maxwell stress and electrical tangent). */
class LinearDielectricT: public DiffusionMaterialT
{
public:

	/** constructor */
	LinearDielectricT(void);

	/** electric permittivity */
	double Permittivity(void) const { return fEpsilon; }

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:

	double fEpsilon; /**< electric permittivity */
};

} // namespace Tahoe
#endif /* _LINEAR_DIELECTRIC_T_H_ */
