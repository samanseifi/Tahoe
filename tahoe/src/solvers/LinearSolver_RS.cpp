/* created: samansei (2008/05/26) on top of paklein (05/30/1996) */
#include "LinearSolver_RS.h"
#include "FEManagerT.h"

using namespace Tahoe;

/* constructors */
LinearSolver_RS::LinearSolver_RS(FEManagerT& fe_manager, int group):
	SolverT(fe_manager, group),
	fFormLHS(1)
{
	SetName("linear_solver_RS");
}

/* signal new reconfigured system */
void LinearSolver_RS::Initialize(int tot_num_eq, int loc_num_eq, int start_eq)
{
	/* inherited */
	SolverT::Initialize(tot_num_eq, loc_num_eq, start_eq);

	/* flag to reform LHS */
	fFormLHS = 1;
}

/* start solution step */
void LinearSolver_RS::InitStep(void)
{
	/* inherited */
	SolverT::InitStep();

	/* no iterations count */
	fNumIteration = 0;
}

/* solve the current step */
SolverT::SolutionStatusT LinearSolver_RS::Solve(int)
{

	try {
	/* initialize */
	fRHS = 0.0;

	/* form the residual force vector */
	fFEManager.FormRHS(Group());

	/* solve equation system */
	if (fFormLHS)
	{
		/* unlock */
		fLHS_lock = kOpen;

		/* initialize */
		fLHS->Clear();

		/* form the stiffness matrix */
		fFEManager.FormLHS(Group(), GlobalT::kNonSymmetric);

		/* flag not to reform */
		fFormLHS = 0;

		/* lock */
		fLHS_lock = kLocked;
	}

	//cout << "When? LinearSolver_RS::Solve" <<  endl;

		/* determine update vector */
	if (!fLHS->Solve(fRHS)) ExceptionT::BadJacobianDet("LinearSolver_RS::Solve");

	//cout << fRHS << endl;

	/* update displacements */
	fFEManager.Update(Group(), fRHS);

	/* relaxation */
	GlobalT::RelaxCodeT relaxcode = fFEManager.RelaxSystem(Group());

	/* relax for configuration change */
	if (relaxcode == GlobalT::kRelax) fFormLHS = 1;
			//NOTE: NLSolver calls "fFEManager.Reinitialize()". Should this happen
			//      here, too? For statics, should also reset the structure of
			//      global stiffness matrix, but since EFG only breaks connections
			//      and doesn't make new ones, this should be OK for now. PAK (03/04/99)

	/* trigger set of new equations */
	if (relaxcode == GlobalT::kReEQ ||
	    relaxcode == GlobalT::kReEQRelax)
		fFEManager.SetEquationSystem(Group());

	return kConverged;
	} /* end try */

	/* not OK */
	catch (ExceptionT::CodeT exc)
	{
		cout << "\n LinearSolver_RS::Solve: caught exception: "
		     << ExceptionT::ToString(exc) << endl;
		return kFailed;
	}
}

/* solve the current step, including ghost particle force */
#ifdef DEM_COUPLING_DEV
SolverT::SolutionStatusT LinearSolver_RS::Solve(int any, FEDEManagerT& fFEDEManager, ArrayT<FBC_CardT>& fGhostFBC)
{
	try {
	// initialize
	fRHS = 0.0;

	// form the residual force vector from ghost particles
	fFEDEManager.FormRHS(Group(), fGhostFBC);

	// form the residual force vector
	fFEManager.FormRHS(Group());

	// solve equation system
	if (fFormLHS)
	{
		// unlock
		fLHS_lock = kOpen;

		// initialize
		fLHS->Clear();

		// form the stiffness matrix
		fFEManager.FormLHS(Group(), GlobalT::kNonSymmetric);

		// flag not to reform
		fFormLHS = 0;

		// lock
		fLHS_lock = kLocked;
	}

	// determine update vector
	if (!fLHS->Solve(fRHS)) ExceptionT::BadJacobianDet("LinearSolver_RS::Solve");

	// update displacements
	fFEManager.Update(Group(), fRHS);

	// relaxation
	GlobalT::RelaxCodeT relaxcode = fFEManager.RelaxSystem(Group());

	// relax for configuration change
	if (relaxcode == GlobalT::kRelax) fFormLHS = 1;
			//NOTE: NLSolver calls "fFEManager.Reinitialize()". Should this happen
			//      here, too? For statics, should also reset the structure of
			//      global stiffness matrix, but since EFG only breaks connections
			//      and doesn't make new ones, this should be OK for now. PAK (03/04/99)

	// trigger set of new equations
	if (relaxcode == GlobalT::kReEQ ||
	    relaxcode == GlobalT::kReEQRelax)
		fFEManager.SetEquationSystem(Group());

	return kConverged;
	} // end try

	// not OK
	catch (ExceptionT::CodeT exc)
	{
		cout << "\n LinearSolver_RS::Solve: caught exception: "
		     << ExceptionT::ToString(exc) << endl;
		return kFailed;
	}
}
#endif

/* signal time step change */
void LinearSolver_RS::SetTimeStep(double dt)
{
	/* inherited */
	SolverT::SetTimeStep(dt);

	/* reform LHS matrix */
	fFormLHS = 1;
}
