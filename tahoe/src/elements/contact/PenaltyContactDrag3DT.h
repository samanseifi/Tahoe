/* $Id: PenaltyContactDrag3DT.h,v 1.3 2005/07/20 06:54:46 paklein Exp $ */
#ifndef _PENALTY_CONTACT_DRAG_3D_T_H_
#define _PENALTY_CONTACT_DRAG_3D_T_H_

/* base classes */
#include "PenaltyContact3DT.h"

/* direct members */
#include "InverseMapT.h"

namespace Tahoe {

/** penalty contact formulation with constant drag force */
class PenaltyContactDrag3DT: public PenaltyContact3DT
{
public:

	/** constructor */
	PenaltyContactDrag3DT(const ElementSupportT& support);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** construct the residual force vector */
	virtual void RHSDriver(void);

	/** construct the effective mass matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);

protected:

	/** magnitude of the drag traction */
	double fDrag;

	/** gap tolerance to enable adhesion */
	double fGapTolerance;

	/** slip tolerance to enable adhesion */
	double fSlipTolerance;
};

} /* namespace Tahoe */

#endif /* _PENALTY_CONTACT_DRAG_3D_T_H_ */
