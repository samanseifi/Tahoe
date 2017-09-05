/* $Id: LinearSolver_RS.h,v 1.11 2008/05/26 18:52:27 bcyansfn Exp $ */
/* created: paklein (05/30/1996) */
#ifndef _LINEAR_SOLVER_RS_H_
#define _LINEAR_SOLVER_RS_H_

/* base class */
#include "SolverT.h"

namespace Tahoe {

/** solver for linear problems */
class LinearSolver_RS: public SolverT
{
public:

	/** constructor */
	LinearSolver_RS(FEManagerT& fe_manager, int group);

	/** configure system */
	virtual void Initialize(int tot_num_eq, int loc_num_eq, int start_eq);

	/** start solution step */
	virtual void InitStep(void);

	/** solve the system over the current time increment.
	 * \param num_iterations maximum number of iterations to execute. Hitting this limit
	 *        does not signal a SolverT::kFailed status, unless solver's internal parameters
	 *        also indicate the solution procedure has failed.
	 * \return one of SolverT::IterationsStatusT */
	virtual SolutionStatusT Solve(int num_iterations);

#ifdef DEM_COUPLING_DEV
	virtual SolutionStatusT Solve(int num_iterations, FEDEManagerT& fFEDEManager, ArrayT<FBC_CardT>& fGhostFBC);
#endif

	/** signal time step change. Chance to clear cached values that may depend on the
	 * time increment. LinearSolver_RS::SetTimeStep triggers recalculation of the LHS
	 * matrix because some time integrators use an effective mass matrix that is a function
	 * of the time increment. */
	virtual void SetTimeStep(double dt);

private:

	/* flag to form RHS */
	int fFormLHS;
		// reform conditions:
		// (1) initially
		// (2) if equation system is reconfigured
		// (3) when the time step changes
};

} // namespace Tahoe
#endif /* _LINEAR_SOLVER_RS_H_ */
