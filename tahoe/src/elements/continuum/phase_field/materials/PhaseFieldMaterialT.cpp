/* Phase-field fracture material */
#include "PhaseFieldMaterialT.h"
#include "PhaseFieldMatSupportT.h"

using namespace Tahoe;

/* array behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<PhaseFieldMaterialT>::fByteCopy = false;
DEFINE_TEMPLATE_STATIC const bool ArrayT<PhaseFieldMaterialT*>::fByteCopy = true;
} /* namespace Tahoe */

/* constructor */
PhaseFieldMaterialT::PhaseFieldMaterialT(void):
	ParameterInterfaceT("phase_field_fracture_material"),
	fPhaseFieldMatSupport(NULL),
	fGc(1.0),
	fEll(1.0),
	fKSmall(1.0e-6)
{

}

/* set support */
void PhaseFieldMaterialT::SetPhaseFieldMatSupport(const PhaseFieldMatSupportT* support)
{
	/* inherited */
	SetMaterialSupport(support);
	fPhaseFieldMatSupport = support;
}

/* describe the parameters needed by the interface */
void PhaseFieldMaterialT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ContinuumMaterialT::DefineParameters(list);

	list.AddParameter(fGc, "Gc");
	list.AddParameter(fEll, "length_scale");
	list.AddParameter(fKSmall, "residual_stiffness");
}

/* accept parameter list */
void PhaseFieldMaterialT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ContinuumMaterialT::TakeParameterList(list);

	fGc = list.GetParameter("Gc");
	fEll = list.GetParameter("length_scale");
	fKSmall = list.GetParameter("residual_stiffness");
}
