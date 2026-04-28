/* $Id: PenaltyContact3DT.h,v 1.7 2004/07/15 08:26:08 paklein Exp $ */
/* created: paklein (02/09/2000) */
#ifndef _PENALTY_CONTACT3D_T_H_
#define _PENALTY_CONTACT3D_T_H_

/* base classes */
#include "Contact3DT.h"

namespace Tahoe {

/** 3D penalty contact element */
class PenaltyContact3DT: public Contact3DT
{
public:

	/** constructor */
	PenaltyContact3DT(const ElementSupportT& support);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** construct the effective mass matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);

	/** construct the residual force vector */
	virtual void RHSDriver(void);

protected:

	/** penalty "stiffness" */
	double fK;

	/** Coulomb friction coefficient (0 = frictionless, default).
	 *  When >0, applies regularised kinetic Coulomb friction:
	 *    f_t = -mu * |f_n| * v_t / sqrt(|v_t|^2 + eps^2)
	 *  with v_t the relative tangential velocity at the contact point.
	 *  See #26. */
	double fMu;

	/** Velocity regularisation for Coulomb friction (mm/us in our units).
	 *  Prevents the discontinuity at zero tangential velocity. */
	double fFrictionEps;

	/** \name element coords and displacements */
	/*@{*/
	dArray2DT fElCoord;
	dArray2DT fElRefCoord;
	dArray2DT fElDisp;
	/*@}*/
	
	/** \name work space */
	/*@{*/
	dMatrixT fdc_du;
	dMatrixT fdn_du;
	dMatrixT fM1;
	dMatrixT fM2;
	dArrayT  fV1;
	/*@}*/
};

} // namespace Tahoe

#endif /* _PENALTY_CONTACT3D_T_H_ */
