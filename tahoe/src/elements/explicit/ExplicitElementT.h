/* ExplicitElementT.h — MVSIZ-batched explicit solid element.
 *
 * Derives from SolidElementT to reuse mass matrix assembly, equation
 * numbering, output, and time integration hooks. Overrides the internal
 * force computation (ElementRHSDriver) with a batched SoA loop using
 * ExplicitKernelT for geometry and ExplicitMaterialT for constitutive law.
 *
 * The batched loop processes elements in blocks of MVSIZ for SIMD
 * auto-vectorization by the compiler. The element driver computes F
 * (deformation gradient) for each batch and passes it to the material.
 *
 * XML usage:
 *   <explicit_solid field_name="displacement" mass_type="lumped_mass">
 *       <quadrilateral/>
 *       <explicit_material density="1.0">
 *           <neo_hookean mu="1.0" kappa="1000.0"/>
 *       </explicit_material>
 *   </explicit_solid>
 */

#ifndef _EXPLICIT_ELEMENT_T_H_
#define _EXPLICIT_ELEMENT_T_H_

#include "UpdatedLagrangianT.h"

namespace Tahoe {

/* forward declarations */
class ExplicitKernelT;
class ExplicitMaterialT;
class ANPHelperT;

class ExplicitElementT : public UpdatedLagrangianT
{
public:

	static const int MVSIZ = 128;

	/** constructor */
	ExplicitElementT(const ElementSupportT& support);

	/** destructor */
	~ExplicitElementT(void);

	/** \name ParameterInterfaceT */
	/*@{*/
	virtual void DefineParameters(ParameterListT& list) const;
	virtual void DefineSubs(SubListT& sub_list) const;
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** override RHSDriver (virtual) to use batched internal force */
	virtual void RHSDriver(void);

	/** override LHSDriver (virtual) to apply mass scaling to lumped mass */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);

private:

	/** Compute internal forces using MVSIZ batched SoA loop.
	 * \param constKd  force coefficient from integrator */
	void BatchedInternalForce(double constKd);

	/** build flat connectivity and equation number arrays at init */
	void BuildFlatArrays(void);

	/** Compute stable time step for all elements.
	 *  Returns the minimum dt across all elements. */
	double ComputeStableTimeStep(void) const;

	/** the element kernel (geometry/shape functions) */
	ExplicitKernelT* fKernel;

	/** the batch material */
	ExplicitMaterialT* fBatchMaterial;

	/** \name pre-computed flat arrays for fast gather/scatter */
	/*@{*/
	int fTotalElements;           /**< total elements across all blocks */
	int* fFlatConn;               /**< [fTotalElements * nen] node IDs */
	int* fFlatEqnos;              /**< [fTotalElements * nen * ndof] equation numbers */
	double* fGlobalRHS;           /**< cached pointer to global RHS vector */
	/*@}*/

	/** \name history variables (plasticity, etc.)
	 *  Layout: fHistory[ip * nhist * total + var * total + elem]
	 *  NULL / zero for stateless hyperelastic materials. */
	/*@{*/
	double* fHistory;
	int fNumHist;                 /**< history vars per IP (0 for elastic) */
	/*@}*/

	/** \name ANP / F-bar (Bonet-Burton 1998, LS-DYNA ELFORM=13).
	 *  When enabled, the dilatation J=det(F) is averaged at nodes and
	 *  fed back as F_bar = (J_bar/J)^(1/nsd) * F before the material call.
	 *  Removes Tet4 (or any reduced-integration element) volumetric locking
	 *  for nearly-incompressible materials.  See ANPHelperT. */
	/*@{*/
	bool fANPEnabled;
	ANPHelperT* fANP;             /**< owned; NULL if disabled */
	double* fVrefE;               /**< [fTotalElements] reference volumes */
	double* fJe;                  /**< [fTotalElements] J per element each step */
	double* fJbarE;               /**< [fTotalElements] nodal-averaged J */
	void BuildANPRefVolumes(void);
	void ComputeAllJe(void);
	/*@}*/

	/** \name hourglass control */
	/*@{*/
	enum HourglassTypeT { kNoHourglass = 0, kViscousHG, kStiffnessHG };
	HourglassTypeT fHourglassType;
	double fHourglassCoeff;       /**< hourglass coefficient (0.01-0.15 typical) */
	/*@}*/

	/** \name mass scaling */
	/*@{*/
	enum MassScalingTypeT { kNoMassScaling = 0, kFixedMassScaling, kAdaptiveMassScaling };
	MassScalingTypeT fMassScalingType;
	double fTargetDt;             /**< target time step for fixed mass scaling */
	double fDtScaleFactor;        /**< safety factor for adaptive (0.9 typical) */
	int fMassScaleInterval;       /**< update interval for adaptive (steps) */
	double* fMassScale;           /**< per-element mass scale factor [fTotalElements] */

	/** Apply mass scaling: increase element mass so that dt_elem >= target_dt.
	 *  Called at init for fixed, periodically for adaptive. */
	void ApplyMassScaling(void);
	/*@}*/
};

} /* namespace Tahoe */
#endif /* _EXPLICIT_ELEMENT_T_H_ */
