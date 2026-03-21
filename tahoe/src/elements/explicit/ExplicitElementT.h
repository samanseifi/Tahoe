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

	/** batched internal force computation — replaces per-element loop */
	virtual void ElementRHSDriver(void);

private:

	/** Compute internal forces using MVSIZ batched SoA loop.
	 * \param constKd  force coefficient from integrator */
	void BatchedInternalForce(double constKd);

	/** the element kernel (geometry/shape functions) */
	ExplicitKernelT* fKernel;

	/** the batch material */
	ExplicitMaterialT* fBatchMaterial;
};

} /* namespace Tahoe */
#endif /* _EXPLICIT_ELEMENT_T_H_ */
