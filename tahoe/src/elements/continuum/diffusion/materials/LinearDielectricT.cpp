/* created for electromechanical coupling in DiffusionElementT */
#include "LinearDielectricT.h"
#include "ParameterListT.h"

using namespace Tahoe;

/* constructor */
LinearDielectricT::LinearDielectricT(void):
	ParameterInterfaceT("linear_dielectric_material"),
	fEpsilon(1.0)
{

}

/* describe the parameters needed by the interface */
void LinearDielectricT::DefineParameters(ParameterListT& list) const
{
	/* inherited (registers density, specific_heat, conductivity) — skip those
	 * by calling the grandparent: ContinuumMaterialT::DefineParameters */
	ContinuumMaterialT::DefineParameters(list);

	list.AddParameter(fEpsilon, "epsilon");
}

/* accept parameter list */
void LinearDielectricT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ContinuumMaterialT::TakeParameterList(list);

	fEpsilon = list.GetParameter("epsilon");

	/* store epsilon as conductivity scalar so base-class k_ij() / q_i()
	 * stay consistent — fConductivity is dimensioned by SetDiffusionMatSupport */
	fConductivity.Identity(fEpsilon);
}
