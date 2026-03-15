/* Phase-field fracture element */
#ifndef _PHASE_FIELD_ELEMENT_T_H_
#define _PHASE_FIELD_ELEMENT_T_H_

/* base class */
#include "ContinuumElementT.h"

/* direct members */
#include "dArray2DT.h"
#include "LocalArrayT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

/* forward declarations */
class ShapeFunctionT;
class PhaseFieldMaterialT;
class PhaseFieldMatSupportT;
class StringT;

/** Phase-field fracture element (AT2 model).
 *
 *  Solves the phase-field evolution equation:
 *    -Gc*ell*Laplacian(d) + (Gc/ell)*d = 2*(1-d)*H
 *
 *  where d is the phase-field damage variable (0=intact, 1=broken),
 *  Gc is fracture toughness, ell is the length scale, and
 *  H = max{psi_total} is the crack driving force history variable.
 *
 *  The element stiffness and residual are:
 *    K_dd = integral[ Gc*ell*B^T*B + (Gc/ell + 2*H)*N^T*N ] dV
 *    f_d  = integral[ Gc*ell*B^T*grad(d) + (Gc/ell + 2*H)*N^T*d - 2*H*N^T ] dV
 *
 *  When coupled to mechanics, the element reads the displacement field
 *  to compute strain energy density at each integration point. */
class PhaseFieldElementT: public ContinuumElementT
{
public:

	/** list/index of nodal outputs */
	enum OutputCodeT {
		iNodalCoord = 0,  /**< (reference) nodal coordinates */
		iNodalDisp  = 1,  /**< nodal phase-field values */
		iMaterialData = 2 /**< material model output */
	};

	/** constructor */
	PhaseFieldElementT(const ElementSupportT& support);

	/** destructor */
	~PhaseFieldElementT(void);

	/** compute nodal force */
	virtual void AddNodalForce(const FieldT& field, int node, dArrayT& force);

	/** returns the stored energy */
	virtual double InternalEnergy(void);

	/** compute specified output parameter and send for smoothing */
	virtual void SendOutput(int kincode);

protected:

	/** initialization functions */
	/*@{*/
	virtual void SetLocalArrays(void);
	virtual void SetShape(void);
	/*@}*/
	virtual void SetGlobalShape(void);

	/** construct the effective mass matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);

	/** form the residual force vector */
	virtual void RHSDriver(void);

	/** increment current element */
	virtual bool NextElement(void);

	/** form the element stiffness matrix */
	virtual void FormStiffness(double constK);

	/** calculate the internal force contribution */
	virtual void FormKd(double constK);

	/** construct a new material support */
	virtual MaterialSupportT* NewMaterialSupport(MaterialSupportT* p = NULL) const;

	/** return a pointer to a new material list */
	virtual MaterialListT* NewMaterialList(const StringT& name, int size);

	/** driver for calculating output values */
	virtual void ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
	                           const iArrayT& e_codes, dArray2DT& e_values);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	virtual void DefineParameters(ParameterListT& list) const;
	virtual void DefineSubs(SubListT& sub_list) const;
	ParameterInterfaceT* NewSub(const StringT& name) const;
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

	/** extract the list of material parameters */
	virtual void CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const;

private:

	/** \name construct output labels array */
	/*@{*/
	virtual void SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const;
	virtual void SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const;
	virtual void GenerateOutputLabels(const iArrayT& n_counts,
		ArrayT<StringT>& n_labels, const iArrayT& e_counts, ArrayT<StringT>& e_labels) const;
	/*@}*/

	/** compute strain energy density at the current integration point
	 *  from the displacement field (if available). Returns 0 if no
	 *  mechanical coupling. */
	double ComputeStrainEnergyDensity(int ip) const;

protected:

	/** current material */
	PhaseFieldMaterialT* fCurrMaterial;

	/** \name mechanical coupling */
	/*@{*/
	/** displacement local array (from mechanical field) */
	LocalArrayT* fLocDisplacement;

	/** deformation gradient at integration points */
	ArrayT<dMatrixT> fF_List;
	dArrayT          fF_all;
	/*@}*/

	/** \name crack driving force history variable H = max{psi} at each IP */
	/*@{*/
	/** H values for the current step: [num_elements * nip] */
	dArrayT fH_current;

	/** H values at the last converged step */
	dArrayT fH_last;

	/** total number of integration points across all elements */
	int fTotalNumIP;
	/*@}*/

	/** \name work space */
	/*@{*/
	dMatrixT fD;  /**< constitutive matrix (identity scaled by Gc*ell) */
	dMatrixT fB;  /**< "strain-displacement" matrix (gradient of shape functions) */
	/*@}*/

	/** field gradients over the element */
	ArrayT<dArrayT> fGradient_list;

	/* parameters */
	static const int NumNodalOutputCodes;

	/** mechanical coupling parameters */
	bool fMechanicalCoupling;

	/** Neo-Hookean parameters for strain energy density computation.
	 *  When coupled to mechanics, these are read from the input file
	 *  so the phase-field element can independently compute psi(F). */
	double fMu;     /**< shear modulus */
	double fLambda;  /**< Lame's first parameter */

private:

	PhaseFieldMatSupportT* fPhaseFieldMatSupport;
};

} /* namespace Tahoe */

#endif /* _PHASE_FIELD_ELEMENT_T_H_ */
