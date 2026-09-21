/* CollocationSolverT.h — direct nodal-collocation transform for meshfree essential BCs.
 *
 * Meshfree RK/MLS shape functions are non-interpolatory: Phi_J(x_I) != delta_IJ. Hence setting
 * the bare nodal coefficient d_I does NOT impose the physical displacement at node I, which is
 *   u(x_I) = sum_J Phi_J(x_I) d_J.
 * "Direct nodal collocation" (Wang & Bazilevs 2024; Belytschko et al.) instead constrains the
 * physical value. For a set of constrained nodes {I} with prescribed physical values ubar_I, the
 * constrained coefficients d_c solve
 *   sum_{J in C} Phi_J(x_I) d_J  =  ubar_I  -  sum_{J not in C} Phi_J(x_I) d_J     (for all I in C)
 * i.e.  Phi_c d_c = ubar - Phi_free d_free, where Phi_c is the (Nc x Nc) restriction of the
 * collocation matrix to the constrained set. This class builds and LU-factors Phi_c once (the
 * reference-config shapes are fixed) and solves it per dof / per step.
 *
 * Reusable by ANY meshfree method: the caller supplies, per constrained node, the support node
 * ids J and the shape values Phi_J(x_I) (see MeshFreeCollocationSupportT). Depends only on
 * toolbox containers, so it is unit-testable in isolation.
 */
#ifndef _COLLOCATION_SOLVER_T_H_
#define _COLLOCATION_SOLVER_T_H_

#include "dArrayT.h"
#include "iArrayT.h"
#include "RaggedArray2DT.h"

namespace Tahoe {

class CollocationSolverT
{
public:

	CollocationSolverT(void);

	/** build + LU-factor the constrained collocation matrix Phi_c.
	 * \param constrained global ids of the constrained nodes (length Nc, no duplicates)
	 * \param max_global_id largest global node id + 1 (sizes the global->constrained map)
	 * \param support per-constrained-node support node ids (row a aligned with constrained[a])
	 * \param phi     per-constrained-node shape values Phi_J(x_{constrained[a]}), aligned with support
	 * Returns false if Phi_c is singular (degenerate constrained set). */
	bool SetCollocation(const iArrayT& constrained, int max_global_id,
		const RaggedArray2DT<int>& support, const RaggedArray2DT<double>& phi);

	/** number of constrained nodes */
	int NumConstrained(void) const { return fConstrained.Length(); }

	/** the constrained node ids (global), in solve order */
	const iArrayT& Constrained(void) const { return fConstrained; }

	/** solve Phi_c d_c = prescribed - Phi_free d_free for one dof.
	 * \param prescribed physical values ubar_I at the constrained nodes (length Nc, solve order)
	 * \param coeff_dof  CURRENT coefficient values for this dof at ALL nodes, indexed by global id
	 *                   (the free-node contribution to the RHS is read from here)
	 * \param d_c        returns the constrained coefficients (length Nc, solve order)
	 * The free contribution uses coeff_dof at non-constrained support nodes; constrained entries of
	 * coeff_dof are ignored (they are what we solve for). */
	void Solve(const dArrayT& prescribed, const dArrayT& coeff_dof, dArrayT& d_c) const;

	/** convenience: physical value at constrained node a given a full coefficient vector for one dof
	 * (for verification): u(x_a) = sum_J Phi_J(x_a) coeff_dof[J]. */
	double Reconstruct(int a, const dArrayT& coeff_dof) const;

private:

	iArrayT fConstrained;         /**< constrained global ids, solve order */
	iArrayT fGlobalToCon;         /**< global id -> constrained local index, -1 if free */
	RaggedArray2DT<int>    fSupport; /**< per-constrained support ids */
	RaggedArray2DT<double> fPhi;     /**< per-constrained shape values, aligned with fSupport */

	/* LU factorization of Phi_c (row-major Nc x Nc), partial pivoting */
	int             fNc;
	dArrayT         fLU;          /**< Nc*Nc, in-place LU */
	iArrayT         fPiv;         /**< pivot row indices */
};

} /* namespace Tahoe */

#endif /* _COLLOCATION_SOLVER_T_H_ */
