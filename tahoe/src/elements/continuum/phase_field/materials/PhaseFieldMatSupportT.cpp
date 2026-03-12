/* Phase-field material support */
#include "PhaseFieldMatSupportT.h"
#include "ElementsConfig.h"

#ifdef CONTINUUM_ELEMENT
#include "PhaseFieldElementT.h"
#endif

using namespace Tahoe;

/* constructor */
PhaseFieldMatSupportT::PhaseFieldMatSupportT(int ndof, int nip):
	MaterialSupportT(ndof, nip),
	fGradient_list(NULL),
	fPhaseField(NULL)
{

}

/* set the source for the gradient information */
void PhaseFieldMatSupportT::SetGradient(const ArrayT<dArrayT>* gradient_list)
{
	fGradient_list = gradient_list;
}

/* set the element group pointer */
void PhaseFieldMatSupportT::SetContinuumElement(const ContinuumElementT* p)
{
	/* inherited */
	MaterialSupportT::SetContinuumElement(p);

#ifdef CONTINUUM_ELEMENT
	fPhaseField = TB_DYNAMIC_CAST(const PhaseFieldElementT*, p);
#endif
}
