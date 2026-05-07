/* BonetTetT.h — average-nodal-pressure Tet4 (Bonet & Burton 1998).
 *
 * Classic-Tahoe (UpdatedLagrangianT) element with the Bonet-Burton F-bar
 * fix for volumetric locking.  Equivalent to LS-DYNA's ELFORM=13.
 *
 * Reuses the ANPHelperT helper (also used by ExplicitElementT) to perform
 * the J-gather → nodal-average → J-bar-scatter pipeline.  The helper is
 * integrator-agnostic; this element drives it from inside the Newton
 * iteration of the implicit / quasistatic solver.
 *
 * Algorithm per step (issue #29: lagged-J̄ + FD tangent):
 *   1. Step-start (InitStep): compute J_e = V_curr/V_ref from the current
 *      (start-of-step) coordinates, run the ANP helper to get J̄_e at
 *      every element.  Cache J̄_e in fJbarE.
 *   2. Inside Newton (RHSDriver / LHSDriver): use the cached J̄_e —
 *      do NOT update it during the inner iteration.  SetGlobalShape()
 *      computes F per element from current u; this override scales F to
 *      F̄ = (J̄_e / J)^(1/3) F using the cached J̄_e.
 *   3. FormStiffness builds the per-element 12×12 tangent by forward
 *      finite differences on the F̄ residual at frozen J̄ (RecomputeF +
 *      ComputeInternalForce per perturbation).  Both residual and tangent
 *      see the same fJbarE, so the local tangent IS the consistent
 *      ∂R(u; J̄_step-start)/∂u — quadratic Newton convergence on
 *      moderately near-incompressible cases.
 *
 * Trade-offs / limitations:
 *   - The converged solution satisfies R(u; J̄_step-start) = 0, not the
 *     fully-coupled R(u; J̄(u)) = 0.  For quasi-static with small load
 *     steps the lag is small; the residual at the next InitStep absorbs
 *     it as J̄ is refreshed.
 *   - The per-element FD tangent misses the cross-element ∂J̄/∂u
 *     coupling (J̄ is nodal-averaged, so a node perturbation changes J̄
 *     in every neighbouring element — that block can't fit in a 12×12
 *     element-local stiffness).  For SEVERELY near-incompressible
 *     problems (κ ≫ μ, e.g. compress_bonet_tet.xml with κ=1000, μ=1)
 *     this missing volumetric coupling makes the local tangent's bulk
 *     stiffness vanish and Newton diverges.  Bonet & Burton 1998 itself
 *     formulates only the residual — for explicit dynamics — and never
 *     derives an implicit tangent; a fully consistent tangent would
 *     follow Bonet/Marriott/Hassan 2001 or Caylak/Mahnken 2012, tracked
 *     as a future analytical-tangent extension.
 *
 * XML usage (drop-in replacement for <updated_lagrangian> when using tets):
 *   <bonet_tet field_name="displacement">
 *     <tetrahedron/>
 *     <large_strain_element_block>
 *       <large_strain_material_3D>
 *         <RG_split_general density="1.0">
 *           <rg_eq_potential><neo-hookean kappa="..." mu="..."/></rg_eq_potential>
 *         </RG_split_general>
 *       </large_strain_material_3D>
 *     </large_strain_element_block>
 *   </bonet_tet>
 */

#ifndef _BONET_TET_T_H_
#define _BONET_TET_T_H_

#include "UpdatedLagrangianT.h"

namespace Tahoe {

class ANPHelperT;

class BonetTetT: public UpdatedLagrangianT
{
public:
	BonetTetT(const ElementSupportT& support);
	~BonetTetT(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

	/** Override: at step start, compute J̄_e from current coordinates and
	 *  cache it for the inner Newton iteration (issue #29 lagged-J̄). */
	virtual void InitStep(void);

protected:

	/** Override: standard F + then F-bar correction (in-place on fF_List). */
	virtual void SetGlobalShape(void);

	/** Override: dispatch to UpdatedLagrangianT::RHSDriver which iterates
	 *  elements and calls our overridden SetGlobalShape (which uses the
	 *  step-start J̄_e cached by InitStep). */
	virtual void RHSDriver(void);

	/** Override: same as RHSDriver — uses the cached J̄_e from InitStep,
	 *  so the inherited F-based analytical tangent is consistent with
	 *  the residual ∂R/∂u at frozen J̄. */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);

	/** Override: numerical (forward-FD) tangent of the F-bar residual at
	 *  the cached (frozen) J̄.  Per-element 12×12 block, FD on 12 dofs.
	 *  With both residual and tangent reading the same fJbarE, the local
	 *  tangent is the consistent ∂R(u; J̄_step-start)/∂u — quadratic Newton
	 *  convergence on moderately incompressible cases.
	 *
	 *  Limitation: the cross-element ∂J̄/∂u coupling (a node perturbation
	 *  changes J̄ in every neighbouring element) cannot fit in a 12×12
	 *  block.  At severe near-incompressibility (κ ≫ μ) the missing
	 *  block kills the local volumetric stiffness and Newton diverges.
	 *  See header docstring for the future analytical-tangent path. */
	virtual void FormStiffness(double constK);

private:

	void BuildANPData(void);            /* allocates fVrefE, fFlatConn, ANP helper */
	void ComputeAllJe(void);            /* fills fJe in current configuration */
	void UpdateJBar(void);              /* ComputeAllJe + fANP->ComputeJBar */
	void ComputeInternalForce(double constK, dArrayT& force);
	void RecomputeF(void);              /* refresh F + F-bar from current
	                                     * fLocDisp / fLocCurrCoords WITHOUT
	                                     * reloading them from the global
	                                     * Field — for FD perturbation. */

	ANPHelperT* fANP;
	int*        fFlatConn;     /* [fNelem * 4] connectivity, 0-based */
	double*     fVrefE;        /* [fNelem] reference tet volume */
	double*     fJe;           /* [fNelem] J = V_cur/V_ref */
	double*     fJbarE;        /* [fNelem] nodal-averaged J — CACHED across
	                            * Newton iters of one step (#29) */
	int         fNelem;
};

} /* namespace Tahoe */
#endif /* _BONET_TET_T_H_ */
