/* $Id: ElementListT.h,v 1.10 2004/07/15 08:25:44 paklein Exp $ */
/* created: paklein (04/20/1998) */
#ifndef _ELEMENTLIST_T_H_
#define _ELEMENTLIST_T_H_

/* base classes */
#include "pArrayT.h"
#include "ParameterInterfaceT.h"

/* direct members */
#include "ElementSupportT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class ElementBaseT;
class eIntegratorT;
class StringT;
class FEManagerT;
class ElementSupportT;

/** list of elements. Constructs list of element objects and
 * provides some attributes. */
class ElementListT: public pArrayT<ElementBaseT*>, public ParameterInterfaceT
{
public:

	/** constructor */
	ElementListT(FEManagerT& fe);

	/** destructor */
	~ElementListT(void);

	/** returns true of ALL element groups have interpolant DOF's */
	bool InterpolantDOFs(void) const;

	/** returns true if contact group present */
	bool HasContact(void) const { return fHasContact; };

	/** change the active element groups.
	 * \param mask list with length of the \e total number of element
	 *        groups with true|false determining whether the element
	 *        group is active. */
	void SetActiveElementGroupMask(const ArrayT<bool>& mask);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list.
	 * \param list input parameter list, which should be validated using ParameterInterfaceT::ValidateParameterList
	 *        to ensure the list conforms to the description defined by the interface. */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:

	/** return a pointer to a new element group or NULL if the request cannot be completed */
	ElementBaseT* NewElement(const StringT& name) const;

private:

	/** data needed for element contructors */
	ElementSupportT fSupport;

	/** cached pointers to element groups */
	ArrayT<ElementBaseT*> fAllElementGroups;
	
	/** true if list contains contact elements */
	bool fHasContact;
};

} /* namespace Tahoe */

#endif /* _ELEMENTLIST_T_H_ */
