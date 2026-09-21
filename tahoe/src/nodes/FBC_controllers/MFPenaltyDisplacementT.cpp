/* MFPenaltyDisplacementT.cpp */
#include "MFPenaltyDisplacementT.h"

#include "MeshFreeCollocationSupportT.h"
#include "ElementBaseT.h"
#include "FieldT.h"
#include "FieldSupportT.h"
#include "ModelManagerT.h"
#include "ScheduleT.h"
#include "eIntegratorT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"   /* StringListT::Extract */
#include "ofstreamT.h"

#include <cmath>

using namespace Tahoe;

/* constructor */
MFPenaltyDisplacementT::MFPenaltyDisplacementT(void):
	fDOF(0),
	fScheduleNum(-1),
	fSchedule(NULL),
	fValue(0.0),
	fPenalty(0.0),
	fDampingRatio(0.0),
	fElementGroup(-1),
	fProvider(NULL),
	fBuilt(false),
	fHaveGapStep(false),
	fCount(0),
	fStiffness(ElementMatrixT::kSymmetric)
{
	SetName("penalty_displacement_meshfree");
}

/* the constraint couples every coefficient in the support of a constrained node */
void MFPenaltyDisplacementT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
	/* inherited */
	FBC_ControllerT::Equations(eq_1, eq_2);

	/* collocation rows and equation numbers (element groups exist by now) */
	if (!fBuilt) BuildOperator();
	else {
		/* equation numbers may have been renumbered */
		const iArray2DT& eqnos = Field().Equations();
		for (int i = 0; i < fNodes.Length(); i++) {
			const int* J = fSupport(i);
			int* eq = fEqnos(i);
			for (int j = 0; j < fSupport.MinorDim(i); j++) eq[j] = eqnos(J[j], fDOF);
		}
	}

	eq_2.Append(&fEqnos);
}

void MFPenaltyDisplacementT::Connectivities(AutoArrayT<const iArray2DT*>& connects_1,
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2,
	AutoArrayT<const iArray2DT*>& equivalent_nodes) const
{
	/* inherited */
	FBC_ControllerT::Connectivities(connects_1, connects_2, equivalent_nodes);
	if (fBuilt) connects_2.Append(&fSupport);
}

void MFPenaltyDisplacementT::InitialCondition(void)
{
	if (!fBuilt) BuildOperator();
	fLambda = 0.0;
	fHaveGapStep = false;
	fUSum = 0.0; fLambdaSum = 0.0; fCount = 0;
}

/* tangent: k Phi Phi^T (+ c/dt-type term through the integrator's damping factor) */
void MFPenaltyDisplacementT::ApplyLHS(GlobalT::SystemTypeT sys_type)
{
#pragma unused(sys_type)
	double constK = 0.0, constC = 0.0;
	int formK = fIntegrator->FormK(constK);
	int formC = fIntegrator->FormC(constC);
	if (!formK && !formC) return;
	if (!fBuilt) BuildOperator();

	for (int i = 0; i < fNodes.Length(); i++)
	{
		int nn = fSupport.MinorDim(i);
		const double* phi = fPhi(i);
		double coeff = (formK ? constK*fPenalty : 0.0) + (formC ? constC*fDashpot[i] : 0.0);

		fStiffness.Dimension(nn);
		fStiffness.SetFormat(ElementMatrixT::kSymmetric);
		for (int a = 0; a < nn; a++)
			for (int b = 0; b < nn; b++)
				fStiffness(a,b) = coeff*phi[a]*phi[b];

		iArrayT eqnos(nn, fEqnos(i));
		FieldSupport().AssembleLHS(fGroup, fStiffness, eqnos);
	}
}

/* force: f_J = Phi_J lambda, lambda = -k g - c dg/dt */
void MFPenaltyDisplacementT::ApplyRHS(void)
{
	double constKd = 0.0;
	int formKd = fIntegrator->FormKd(constKd);
	if (!formKd) return;
	if (!fBuilt) BuildOperator();

	double ubar = fValue*(fSchedule ? fSchedule->Value() : 1.0);
	double dt = FieldSupport().TimeStep();

	for (int i = 0; i < fNodes.Length(); i++)
	{
		double g = Gap(i, ubar);
		double gdot = (fHaveGapStep && dt > 0.0) ? (g - fGapStep[i])/dt : 0.0;
		fLambda[i] = -fPenalty*g - fDashpot[i]*gdot;

		int nn = fSupport.MinorDim(i);
		const double* phi = fPhi(i);
		fForce.Dimension(nn);
		for (int j = 0; j < nn; j++) fForce[j] = constKd*phi[j]*fLambda[i];

		iArrayT eqnos(nn, fEqnos(i));
		FieldSupport().AssembleRHS(fGroup, fForce, eqnos);
	}
}

/* store the gap for the rate term and accumulate the output means */
void MFPenaltyDisplacementT::CloseStep(void)
{
	if (!fBuilt) return;
	double ubar = fValue*(fSchedule ? fSchedule->Value() : 1.0);
	for (int i = 0; i < fNodes.Length(); i++) {
		double g = Gap(i, ubar);
		fGapStep[i] = g;
		fUSum[i] += g + ubar;
		fLambdaSum[i] += fLambda[i];
	}
	fHaveGapStep = true;
	fCount++;
}

void MFPenaltyDisplacementT::WriteOutput(ostream& out) const
{
	if (!fBuilt) return;
	double ubar = fValue*(fSchedule ? fSchedule->Value() : 1.0);
	double n = (fCount > 0) ? double(fCount) : 1.0;
	int prec = out.precision();
	out.precision(8);
	for (int i = 0; i < fNodes.Length(); i++) {
		out << " penalty_displacement_meshfree: node " << fNodes[i] + 1
		    << " u = " << Gap(i, ubar) + ubar
		    << " lambda = " << fLambda[i]
		    << " u_mean = " << fUSum[i]/n
		    << " lambda_mean = " << fLambdaSum[i]/n << '\n';
		fUSum[i] = 0.0; fLambdaSum[i] = 0.0;
	}
	fCount = 0;
	out.precision(prec);
}

/* describe the parameters needed by the interface */
void MFPenaltyDisplacementT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	FBC_ControllerT::DefineParameters(list);

	ParameterT dof(ParameterT::Integer, "dof");                   /* 1-based */
	dof.AddLimit(1, LimitT::LowerInclusive);
	list.AddParameter(dof);

	ParameterT schedule(ParameterT::Integer, "schedule");
	schedule.AddLimit(1, LimitT::LowerInclusive);
	list.AddParameter(schedule);

	ParameterT value(ParameterT::Double, "value");               /* prescribed PHYSICAL value */
	value.SetDefault(0.0);
	list.AddParameter(value);

	ParameterT penalty(ParameterT::Double, "penalty");
	penalty.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(penalty);

	ParameterT damping(ParameterT::Double, "damping_ratio");     /* fraction of critical */
	damping.SetDefault(0.0);
	damping.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(damping);

	ParameterT group(ParameterT::Integer, "element_group");      /* 1-based; 0 = auto-discover */
	group.SetDefault(0);
	group.AddLimit(0, LimitT::LowerInclusive);
	list.AddParameter(group);
}

/* information about subordinate parameter lists */
void MFPenaltyDisplacementT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FBC_ControllerT::DefineSubs(sub_list);

	/* constrained node sets */
	sub_list.AddSub("node_ID_list");
}

/* accept parameter list */
void MFPenaltyDisplacementT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "MFPenaltyDisplacementT::TakeParameterList";

	/* inherited */
	FBC_ControllerT::TakeParameterList(list);

	fDOF = list.GetParameter("dof"); fDOF--;
	fScheduleNum = list.GetParameter("schedule"); fScheduleNum--;
	fValue = list.GetParameter("value");
	fPenalty = list.GetParameter("penalty");
	fDampingRatio = list.GetParameter("damping_ratio");
	fElementGroup = list.GetParameter("element_group"); fElementGroup--;

	if (fDOF >= Field().NumDOF())
		ExceptionT::BadInputValue(caller, "dof %d exceeds the field dimension %d", fDOF+1, Field().NumDOF());

	fSchedule = FieldSupport().Schedule(fScheduleNum);
	if (!fSchedule) ExceptionT::BadInputValue(caller, "could not resolve schedule %d", fScheduleNum+1);

	/* constrained nodes */
	StringListT::Extract(list.GetList("node_ID_list"), fID_List);
	ModelManagerT& model = FieldSupport().ModelManager();
	model.ManyNodeSets(fID_List, fNodes);
	if (fNodes.Length() == 0) ExceptionT::BadInputValue(caller, "empty constrained node set");

	/* the collocation operator is built once the element groups are guaranteed constructed */
	fBuilt = false;
}

/***********************************************************************
 * Private
 ***********************************************************************/

void MFPenaltyDisplacementT::BuildOperator(void)
{
	const char caller[] = "MFPenaltyDisplacementT::BuildOperator";

	/* locate the meshfree provider */
	if (fElementGroup >= 0) {
		ElementBaseT& group = FieldSupport().ElementGroup(fElementGroup);
		fProvider = dynamic_cast<const MeshFreeCollocationSupportT*>(&group);
		if (!fProvider)
			ExceptionT::GeneralFail(caller, "element_group %d does not provide meshfree collocation data",
				fElementGroup+1);
	} else {
		for (int g = 0; g < FieldSupport().NumElementGroups() && !fProvider; g++)
			fProvider = dynamic_cast<const MeshFreeCollocationSupportT*>(&FieldSupport().ElementGroup(g));
		if (!fProvider)
			ExceptionT::GeneralFail(caller, "no element group provides meshfree collocation data");
	}

	/* collocation rows Phi_J(x_A) */
	if (!fProvider->CollocationData(fNodes, fSupport, fPhi))
		ExceptionT::GeneralFail(caller, "provider returned no collocation data for the constrained set");

	/* equation numbers of the support coefficients */
	const iArray2DT& eqnos = Field().Equations();
	iArrayT counts(fNodes.Length());
	for (int i = 0; i < fNodes.Length(); i++) counts[i] = fSupport.MinorDim(i);
	fEqnos.Configure(counts);
	for (int i = 0; i < fNodes.Length(); i++) {
		const int* J = fSupport(i);
		int* eq = fEqnos(i);
		for (int j = 0; j < counts[i]; j++) eq[j] = eqnos(J[j], fDOF);
	}

	/* dashpots: c = ratio 2 sqrt(k m_eff), m_eff = 1/sum_J Phi_J^2/m_J */
	fDashpot.Dimension(fNodes.Length());
	fDashpot = 0.0;
	if (fDampingRatio > 0.0)
	{
		for (int i = 0; i < fNodes.Length(); i++) {
			int nn = fSupport.MinorDim(i);
			iArrayT support(nn, const_cast<int*>(fSupport(i)));
			dArrayT mass;
			if (!fProvider->NodalMass(support, mass))
				ExceptionT::GeneralFail(caller, "damping_ratio > 0 but the element group reports no nodal mass");
			const double* phi = fPhi(i);
			double inv = 0.0;
			for (int j = 0; j < nn; j++)
				if (mass[j] > 0.0) inv += phi[j]*phi[j]/mass[j];
			double m_eff = (inv > 0.0) ? 1.0/inv : 0.0;
			fDashpot[i] = fDampingRatio*2.0*sqrt(fPenalty*m_eff);
		}
	}

	/* state */
	fLambda.Dimension(fNodes.Length()); fLambda = 0.0;
	fGapStep.Dimension(fNodes.Length()); fGapStep = 0.0; fHaveGapStep = false;
	fUSum.Dimension(fNodes.Length()); fUSum = 0.0;
	fLambdaSum.Dimension(fNodes.Length()); fLambdaSum = 0.0;
	fCount = 0;

	/* report */
	ofstreamT& out = FieldSupport().Output();
	out << "\n M e s h f r e e   P e n a l t y   D i s p l a c e m e n t :\n\n";
	out << " Number of constrained nodes . . . . . . . . . = " << fNodes.Length() << '\n';
	out << " Constrained dof . . . . . . . . . . . . . . . = " << fDOF+1 << '\n';
	out << " Penalty stiffness . . . . . . . . . . . . . . = " << fPenalty << '\n';
	out << " Damping ratio . . . . . . . . . . . . . . . . = " << fDampingRatio << '\n';
	for (int i = 0; i < fNodes.Length(); i++)
		out << "    node " << fNodes[i]+1 << ": support " << fSupport.MinorDim(i)
		    << "  dashpot " << fDashpot[i] << '\n';

	fBuilt = true;
}

double MFPenaltyDisplacementT::Gap(int i, double ubar) const
{
	const dArray2DT& u = Field()[0];
	int nn = fSupport.MinorDim(i);
	const int* J = fSupport(i);
	const double* phi = fPhi(i);
	double g = -ubar;
	for (int j = 0; j < nn; j++) g += phi[j]*u(J[j], fDOF);
	return g;
}
