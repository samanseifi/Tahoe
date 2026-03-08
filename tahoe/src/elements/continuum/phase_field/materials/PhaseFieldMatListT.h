/* Phase-field material list */
#ifndef _PHASE_FIELD_MAT_LIST_T_H_
#define _PHASE_FIELD_MAT_LIST_T_H_

/* base class */
#include "MaterialListT.h"

namespace Tahoe {

/* forward declarations */
class PhaseFieldMatSupportT;
class PhaseFieldMaterialT;

/** list of materials for phase-field fracture analysis */
class PhaseFieldMatListT: public MaterialListT
{
public:

	/** constructors */
	PhaseFieldMatListT(int length, const PhaseFieldMatSupportT& support);
	PhaseFieldMatListT(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	virtual void DefineSubs(SubListT& sub_list) const;
	virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order,
		SubListT& sub_lists) const;
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

	/** construct the specified material or NULL */
	PhaseFieldMaterialT* NewPhaseFieldMaterial(const StringT& name) const;

private:

	/** support for phase-field materials */
	const PhaseFieldMatSupportT* fPhaseFieldMatSupport;
};

} // namespace Tahoe

#endif /* _PHASE_FIELD_MAT_LIST_T_H_ */
