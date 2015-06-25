/* $Id: J2_C0Hardening1DT.h,v 1.1 2004/07/20 23:38:26 rdorgan Exp $ */
#ifndef _J2_C0_HARD_1D_T_H_
#define _J2_C0_HARD_1D_T_H_

/* base class */
#include "ParameterInterfaceT.h"

/* direct members */
#include "C1FunctionT.h"

namespace Tahoe {

/* forward declarations */
class dSymMatrixT;

/** parameters for J2 plasticity with general yield function */
class J2_C0Hardening1DT: virtual public ParameterInterfaceT
{
public:

	/** constructor */
	J2_C0Hardening1DT(void);

	/* destructor */
	virtual ~J2_C0Hardening1DT(void);

	/** returns the value value of the yield function */
	double YieldCondition(const dSymMatrixT& relstress, double alpha) const;

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** \name hardening functions and their 1st derivatives */
	/*@{*/
	double   H(double) const { return 0.0; }; // no kinematic hardening yet
	double  dH(double) const { return 0.0; };
	double  K(double a) const;
	double dK(double a) const;
	/*@}*/

protected:

	/** C1 isotropic hardening function */
	C1FunctionT* fK;	

	/** true if the hardening function is linear */
	bool fIsLinear;
};

/* hardening functions and their 1st derivatives */
inline double J2_C0Hardening1DT::K(double a) const { return fK->Function(a); }
inline double J2_C0Hardening1DT::dK(double a) const { return fK->DFunction(a); }

} /* namespace Tahoe */

#endif /* _J2_C0_HARD_1D_T_H_ */
