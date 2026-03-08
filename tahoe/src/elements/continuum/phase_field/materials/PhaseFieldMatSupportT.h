/* Phase-field material support */
#ifndef _PHASE_FIELD_MAT_SUPPORT_T_H_
#define _PHASE_FIELD_MAT_SUPPORT_T_H_

/* base class */
#include "MaterialSupportT.h"

/* direct members */
#include "dArrayT.h"

namespace Tahoe {

/* forward declarations */
class PhaseFieldElementT;

/** support for phase-field materials */
class PhaseFieldMatSupportT: public MaterialSupportT
{
public:

	/** constructor */
	PhaseFieldMatSupportT(int ndof, int nip);

	/** \name field gradients */
	/*@{*/
	/** field gradient at the current integration point */
	const dArrayT& Gradient(void) const;

	/** field gradient at the specified integration point */
	const dArrayT& Gradient(int ip) const;

	/** set the source for the gradient information */
	void SetGradient(const ArrayT<dArrayT>* gradient_list);
	/*@}*/

	/** \name host code information */
	/*@{*/
	const PhaseFieldElementT* PhaseField(void) const { return fPhaseField; };
	virtual void SetContinuumElement(const ContinuumElementT* p);
	/*@}*/

private:

	/** field gradient: [nip] x [nsd] */
	const ArrayT<dArrayT>* fGradient_list;

	/** pointer to the phase-field element */
	const PhaseFieldElementT* fPhaseField;
};

/* inlines */
inline const dArrayT& PhaseFieldMatSupportT::Gradient(int ip) const
{
	return (*fGradient_list)[ip];
}

inline const dArrayT& PhaseFieldMatSupportT::Gradient(void) const
{
	return (*fGradient_list)[CurrIP()];
}

} /* namespace Tahoe */

#endif /* _PHASE_FIELD_MAT_SUPPORT_T_H_ */
