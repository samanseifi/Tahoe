/* RKShellT.cpp — meshfree RKPM Kirchhoff–Love shell element (epic #59) */
#include "RKShellT.h"

#include "ParameterT.h"
#include "LimitT.h"
#include "ExceptionT.h"

using namespace Tahoe;

/* constructor */
RKShellT::RKShellT(const ElementSupportT& support):
	MeshFreeFSSolidT(support),
	fThickness(0.0)
{
	SetName("meshfree_kl_shell");
}

/* describe the parameters needed by the interface */
void RKShellT::DefineParameters(ParameterListT& list) const
{
	/* inherited (reuses the meshfree solid support + material structure) */
	MeshFreeFSSolidT::DefineParameters(list);

	/* shell mid-surface thickness (> 0) */
	ParameterT thickness(fThickness, "shell_thickness");
	LimitT positive(0.0, LimitT::Lower);
	thickness.AddLimit(positive);
	list.AddParameter(thickness);
}

/* accept parameter list */
void RKShellT::TakeParameterList(const ParameterListT& list)
{
	/* inherited setup (meshfree support, shapes, material list) */
	MeshFreeFSSolidT::TakeParameterList(list);

	/* shell thickness */
	fThickness = list.GetParameter("shell_thickness");
}

/* form the residual (internal force) */
void RKShellT::RHSDriver(void)
{
	ExceptionT::GeneralFail("RKShellT::RHSDriver",
		"KL-shell internal force is not yet ported into the element. The formulation "
		"kernels are validated standalone (contrib/kl_kinematics, kl_bmatrix, kl_stress, "
		"kl_assembly; run `ctest -R meshfree_`). Porting tracked in issue #66.");
}

/* form the tangent stiffness */
void RKShellT::LHSDriver(GlobalT::SystemTypeT sys_type)
{
#pragma unused(sys_type)
	ExceptionT::GeneralFail("RKShellT::LHSDriver",
		"KL-shell tangent stiffness is not yet ported into the element. See issue #66.");
}
