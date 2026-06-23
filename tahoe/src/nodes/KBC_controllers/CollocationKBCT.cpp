/* CollocationKBCT.cpp */
#include "CollocationKBCT.h"

#include "MeshFreeCollocationSupportT.h"
#include "ElementBaseT.h"
#include "FieldT.h"
#include "ScheduleT.h"
#include "ParameterUtils.h"   /* StringListT::Extract */
#include "RaggedArray2DT.h"
#include "dArray2DT.h"

using namespace Tahoe;

CollocationKBCT::CollocationKBCT(const BasicSupportT& support):
	KBC_ControllerT(support),
	fElementGroup(-1), fDOF(0), fScheduleNum(-1), fSchedule(NULL),
	fValue(0.0), fCode(KBC_CardT::kDsp), fProvider(NULL), fBuilt(false)
{
	SetName("collocation_KBC");
}

void CollocationKBCT::DefineParameters(ParameterListT& list) const
{
	KBC_ControllerT::DefineParameters(list);

	ParameterT field(ParameterT::Word, "field");
	field.SetDefault("displacement");
	list.AddParameter(field);

	ParameterT group(ParameterT::Integer, "element_group");   /* 1-based provider group; 0 = auto-discover */
	group.SetDefault(0);
	group.AddLimit(0, LimitT::LowerInclusive);
	list.AddParameter(group);

	ParameterT dof(ParameterT::Integer, "dof");               /* 1-based */
	dof.AddLimit(1, LimitT::LowerInclusive);
	list.AddParameter(dof);

	ParameterT sched(ParameterT::Integer, "schedule");
	sched.AddLimit(1, LimitT::LowerInclusive);
	list.AddParameter(sched);

	ParameterT value(ParameterT::Double, "value");            /* prescribed PHYSICAL value */
	value.SetDefault(0.0);
	list.AddParameter(value);

	ParameterT code(ParameterT::Word, "code");                /* "displacement" | "velocity" */
	code.SetDefault("displacement");
	list.AddParameter(code);
}

void CollocationKBCT::DefineSubs(SubListT& sub_list) const
{
	KBC_ControllerT::DefineSubs(sub_list);
	sub_list.AddSub("node_ID_list");   /* constrained node set(s) */
}

ParameterInterfaceT* CollocationKBCT::NewSub(const StringT& name) const
{
	return KBC_ControllerT::NewSub(name);   /* node_ID_list handled by the base */
}

void CollocationKBCT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "CollocationKBCT::TakeParameterList";
	KBC_ControllerT::TakeParameterList(list);

	fFieldName = list.GetParameter("field");
	fElementGroup = list.GetParameter("element_group"); fElementGroup--;
	fDOF = list.GetParameter("dof"); fDOF--;
	fScheduleNum = list.GetParameter("schedule"); fScheduleNum--;
	fValue = list.GetParameter("value");

	StringT code = list.GetParameter("code");
	if (code == "velocity") fCode = KBC_CardT::kVel;
	else                    fCode = KBC_CardT::kDsp;

	fSchedule = fSupport.Schedule(fScheduleNum);
	if (!fSchedule) ExceptionT::BadInputValue(caller, "could not resolve schedule %d", fScheduleNum+1);

	/* constrained nodes */
	StringListT::Extract(list.GetList("node_ID_list"), fID_List);
	GetNodes(fID_List, fNodes);
	if (fNodes.Length() == 0) ExceptionT::BadInputValue(caller, "empty constrained node set");

	/* one card per constrained node; values filled each InitStep. Use a unit dummy schedule so the
	 * integrator applies the solved coefficient directly (the real schedule is folded in by us). */
	fKBC_Cards.Dimension(fNodes.Length());
	for (int i = 0; i < fNodes.Length(); i++)
		fKBC_Cards[i].SetValues(fNodes[i], fDOF, fCode, NULL, 0.0);

	/* the collocation operator is built lazily on the first InitStep, once the element group that
	 * provides the meshfree shapes is guaranteed constructed (avoids setup-order coupling). */
	fBuilt = false;
}

void CollocationKBCT::Configure(const StringT& field, int dof, KBC_CardT::CodeT code,
	const ScheduleT* schedule, double value, const iArrayT& nodes)
{
	fFieldName = field;
	fDOF = dof;
	fCode = code;
	fSchedule = schedule;
	fValue = value;
	fElementGroup = -1;     /* auto-discover the meshfree provider */
	fNodes = nodes;

	fKBC_Cards.Dimension(fNodes.Length());
	for (int i = 0; i < fNodes.Length(); i++)
		fKBC_Cards[i].SetValues(fNodes[i], fDOF, fCode, NULL, 0.0);
	fBuilt = false;
}

void CollocationKBCT::BuildOperator(void)
{
	const char caller[] = "CollocationKBCT::BuildOperator";

	/* locate the meshfree collocation provider (cross-cast through RTTI) */
	if (fElementGroup >= 0) {
		ElementBaseT& group = fSupport.ElementGroup(fElementGroup);
		fProvider = dynamic_cast<const MeshFreeCollocationSupportT*>(&group);
		if (!fProvider)
			ExceptionT::GeneralFail(caller, "element_group %d does not provide meshfree collocation data",
				fElementGroup+1);
	} else {
		/* auto-discover: the single element group that implements the provider interface */
		for (int g = 0; g < fSupport.NumElementGroups() && !fProvider; g++)
			fProvider = dynamic_cast<const MeshFreeCollocationSupportT*>(&fSupport.ElementGroup(g));
		if (!fProvider)
			ExceptionT::GeneralFail(caller, "no element group provides meshfree collocation data "
				"(needed for collocation BC)");
	}

	/* shape rows phi_J(x_I) for the constrained nodes */
	RaggedArray2DT<int> support; RaggedArray2DT<double> phi;
	if (!fProvider->CollocationData(fNodes, support, phi))
		ExceptionT::GeneralFail(caller, "provider returned no collocation data for the constrained set");

	/* size the global->constrained map by the field node count */
	const FieldT* field = fSupport.Field(fFieldName);
	if (!field) ExceptionT::GeneralFail(caller, "could not resolve field \"%s\"", fFieldName.Pointer());
	int max_gid = (*field)[0].MajorDim();

	if (!fSolver.SetCollocation(fNodes, max_gid, support, phi))
		ExceptionT::GeneralFail(caller, "singular constrained collocation matrix (degenerate node set)");

	fBuilt = true;
}

void CollocationKBCT::InitStep(void)
{
	KBC_ControllerT::InitStep();
	if (!fBuilt) BuildOperator();

	/* current coefficients for this dof (free-node contributions to the RHS) */
	const FieldT* field = fSupport.Field(fFieldName);
	int order = (fCode == KBC_CardT::kVel) ? 1 : 0;
	const dArray2DT& d = (*field)[order];
	int max_gid = d.MajorDim();
	dArrayT coeff(max_gid);
	for (int g = 0; g < max_gid; g++) coeff[g] = d(g, fDOF);

	/* prescribed PHYSICAL value at every constrained node = value * schedule(t) */
	double ubar = fValue*(fSchedule ? fSchedule->Value() : 1.0);
	dArrayT prescribed(fNodes.Length());
	prescribed = ubar;

	/* solve Phi_c d_c = prescribed - Phi_free d_free, emit the coefficients as card values */
	dArrayT d_c;
	fSolver.Solve(prescribed, coeff, d_c);
	for (int a = 0; a < fNodes.Length(); a++)
		fKBC_Cards[a].SetValues(fNodes[a], fDOF, fCode, NULL, d_c[a]);
}
