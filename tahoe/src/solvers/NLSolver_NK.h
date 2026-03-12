/* created: (2026) */
#ifndef _NL_SOLVER_NK_H_
#define _NL_SOLVER_NK_H_

/* base class */
#include "NLSolver.h"

/* direct members */
#include "dArray2DT.h"
#include "ArrayT.h"

namespace Tahoe {

/** Newton-Krylov nonlinear solver.
 *
 * Each outer Newton iteration solves the linear system
 *
 *   K * delta_u = -F
 *
 * using restarted GMRES(m) with a Jacobi (diagonal) left preconditioner
 * instead of a direct factorization.  This avoids the cost of a full LU
 * or Cholesky factorization and can be memory-efficient for large problems.
 *
 * The matrix type should support GlobalMatrixT::Multx() and
 * GlobalMatrixT::CopyDiagonal().  All profile-solver (CCSMatrixT) and
 * MSRMatrixT-based (SPOOLES, MUMPS) matrices satisfy these requirements. */
class NLSolver_NK: public NLSolver
{
public:

	/** constructor */
	NLSolver_NK(FEManagerT& fe_manager, int group);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** do one iteration: GMRES inner solve followed by solution update */
	virtual void Iterate(void);

private:

	/** restarted GMRES(m) solver.
	 * Solves A*x = b using left Jacobi preconditioning.
	 * \param x on entry: initial guess (typically zero); on exit: solution
	 * \param b right-hand side (not modified)
	 * \return true if the specified tolerance was reached */
	bool GMRES(dArrayT& x, const dArrayT& b) const;

	/** apply a Givens rotation in place: [h1; h2] <- G * [h1; h2] */
	static void ApplyGivens(double c, double s, double& h1, double& h2);

private:

	/** \name GMRES parameters */
	/*@{*/
	int    fKrylovRestart;     /**< restart dimension m */
	double fLinearTolerance;   /**< relative tolerance for the inner solve */
	int    fMaxLinearIter;     /**< maximum total GMRES iterations (0 = unlimited) */
	/*@}*/

	/** \name GMRES workspace (allocated during Initialize) */
	/*@{*/
	mutable ArrayT<dArrayT> fQ;   /**< Krylov basis vectors q[0..m] */
	mutable dArray2DT       fH;   /**< upper Hessenberg matrix (m+1) x m */
	mutable dArrayT         fg;   /**< rotated RHS vector (m+1) */
	mutable dArrayT         fcs;  /**< cosines for Givens rotations */
	mutable dArrayT         fsn;  /**< sines  for Givens rotations */
	mutable dArrayT         fy;   /**< least-squares solution */
	mutable dArrayT         fz;   /**< preconditioned vector workspace */
	mutable dArrayT         fw;   /**< Arnoldi vector workspace */
	mutable dArrayT         fdiag;/**< diagonal preconditioner (1/K_ii) */
	/*@}*/
};

} // namespace Tahoe
#endif /* _NL_SOLVER_NK_H_ */
