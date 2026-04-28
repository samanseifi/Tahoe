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
 * Algorithm per Newton iteration:
 *   1. Pre-pass: loop over all elements, compute J_e = V_curr/V_ref
 *      (cheap closed-form for Tet4 from current and reference coords)
 *   2. ANP helper: J_n = sum V_e^ref J_e / sum V_e^ref, J_bar_e = avg of J_n
 *   3. Standard FormKd loop: SetGlobalShape() computes F per element;
 *      THIS class's override scales F to F_bar = (J_bar/J)^(1/3) F
 *      before the material is called.
 *
 * The numerical stiffness tangent (used by SimoQ1P0) gives quadratic Newton
 * convergence for this kind of F-bar element; a future analytical tangent
 * would be cheaper but not required for correctness.
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

protected:

	/** Override: standard F + then F-bar correction (in-place on fF_List). */
	virtual void SetGlobalShape(void);

	/** Override: pre-pass to compute J_e and J_bar_e for all elements,
	 *  then dispatch to UpdatedLagrangianT::RHSDriver which iterates and
	 *  calls our overridden SetGlobalShape per element. */
	virtual void RHSDriver(void);

	/** Override: same pre-pass, then UpdatedLagrangianT::LHSDriver. */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);

private:

	void BuildANPData(void);   /* allocates fVrefE, fFlatConn, ANP helper */
	void ComputeAllJe(void);   /* fills fJe in current configuration */

	ANPHelperT* fANP;
	int*        fFlatConn;     /* [fNelem * 4] connectivity, 0-based */
	double*     fVrefE;        /* [fNelem] reference tet volume */
	double*     fJe;           /* [fNelem] J = V_cur/V_ref each Newton iter */
	double*     fJbarE;        /* [fNelem] nodal-averaged J */
	int         fNelem;
};

} /* namespace Tahoe */
#endif /* _BONET_TET_T_H_ */
