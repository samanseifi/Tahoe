/* $Id: TotalLagrangianT.h,v 1.8 2004/07/15 08:26:27 paklein Exp $ */
/* created: paklein (07/03/1996) */
#ifndef _TOTAL_LAGRANGRIAN_T_H_
#define _TOTAL_LAGRANGRIAN_T_H_

/* base class */
#include "FiniteStrainT.h"

/* direct members */
#include "dMatrixT.h"

namespace Tahoe {

/** total Lagrangian, finite strain element */
class TotalLagrangianT: public FiniteStrainT
{
public:

	/** constructors */
	TotalLagrangianT(const ElementSupportT& support);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/
		
protected:

	/** form the element stiffness matrix */
	virtual void FormStiffness(double constK);

	/** calculate the internal force contribution ("-k*d") */
	virtual void FormKd(double constK);

protected:

	/** \name workspace */
	/*@{*/
	dMatrixT fStressMat;   /**< space for a stress tensor */
	dMatrixT fStressStiff; /**< compact stress stiffness contribution */
	dMatrixT fGradNa;      /**< shape function gradients matrix */
	
	dArrayT   fTemp2;
	dMatrixT  fTempMat1, fTempMat2;
	dArray2DT fDNa_x;
	/*@}*/
};

} // namespace Tahoe 
#endif /* _TOTAL_LAGRANGRIAN_T_H_ */
