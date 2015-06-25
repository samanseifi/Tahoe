/* $Id: SSSolidMatList1DT.h,v 1.3 2004/07/20 23:21:30 rdorgan Exp $ */
#ifndef _SS_MATLIST_1D_T_H_
#define _SS_MATLIST_1D_T_H_

/* base classes */
#include "SolidMatListT.h"
#include "SolidT.h"

namespace Tahoe {

/* forward declaration */
class SSSolidMatT;

/** materials list for 1D structural analysis */
class SSSolidMatList1DT: public SolidMatListT, public SolidT
{
public:

	/** constructor */
	SSSolidMatList1DT(int length, const SSMatSupportT& support);
	SSSolidMatList1DT(void);

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

	/** construct the specified material or NULL if the request cannot be completed */
	SSSolidMatT* NewSSSolidMat(const StringT& name) const;

private:

	/** support for small strain materials */
	const SSMatSupportT* fSSMatSupport;
};

} /* namespace Tahoe */

#endif /* _SS_MATLIST_1D_T_H_ */
