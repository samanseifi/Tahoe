/* Phase-field material list */
#include "PhaseFieldMatListT.h"
#include "PhaseFieldMatSupportT.h"

/* phase-field materials */
#include "PhaseFieldMaterialT.h"

using namespace Tahoe;

/* constructors */
PhaseFieldMatListT::PhaseFieldMatListT(int length, const PhaseFieldMatSupportT& support):
	MaterialListT(length),
	fPhaseFieldMatSupport(&support)
{
	SetName("phase_field_material");
}

PhaseFieldMatListT::PhaseFieldMatListT(void):
	fPhaseFieldMatSupport(NULL)
{
	SetName("phase_field_material");
}

/* information about subordinate parameter lists */
void PhaseFieldMatListT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	MaterialListT::DefineSubs(sub_list);

	/* an array of choices */
	sub_list.AddSub("phase_field_material_list", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void PhaseFieldMatListT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order,
	SubListT& sub_lists) const
{
	if (name == "phase_field_material_list")
	{
		order = ParameterListT::Choice;
		sub_lists.AddSub("phase_field_fracture_material");
	}
	else /* inherited */
		ParameterInterfaceT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* PhaseFieldMatListT::NewSub(const StringT& name) const
{
	PhaseFieldMaterialT* material = NewPhaseFieldMaterial(name);
	if (material)
		return material;
	else /* inherited */
		return MaterialListT::NewSub(name);
}

/* accept parameter list */
void PhaseFieldMatListT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	MaterialListT::TakeParameterList(list);

	/* construct materials */
	AutoArrayT<PhaseFieldMaterialT*> materials;
	const ArrayT<ParameterListT>& subs = list.Lists();
	for (int i = 0; i < subs.Length(); i++) {
		const ParameterListT& sub = subs[i];
		PhaseFieldMaterialT* mat = NewPhaseFieldMaterial(sub.Name());
		if (mat) {
			materials.Append(mat);
			mat->TakeParameterList(sub);
		}
	}

	/* transfer */
	Dimension(materials.Length());
	for (int i = 0; i < materials.Length(); i++)
		fArray[i] = materials[i];
}

/* construct the specified material or NULL */
PhaseFieldMaterialT* PhaseFieldMatListT::NewPhaseFieldMaterial(const StringT& name) const
{
	PhaseFieldMaterialT* mat = NULL;

	if (name == "phase_field_fracture_material")
		mat = new PhaseFieldMaterialT;

	/* set support */
	if (mat) mat->SetPhaseFieldMatSupport(fPhaseFieldMatSupport);

	return mat;
}
