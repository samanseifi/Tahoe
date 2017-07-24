/* $Id: SimoQ1P0.h,v 1.9 2004/07/15 08:26:27 paklein Exp $ */
#ifndef _SIMO_Q1_P0_H_
#define _SIMO_Q1_P0_H_

/* base classes */
#include "UpdatedLagrangianT.h"
#include "FSSolidMatT.h"

namespace Tahoe {

/** finite strain, mixed element formulation.
 * The formulation is due to Simo, Taylor, and Pister, CMAME \b 51,
 * 177-208, 1985. The basic idea behind the formulation is to
 * represent the pressure and dilatation \f$ \Theta \f$ as separate
 * fields from the displacement. For the continuous case, the determinant
 * of the deformation gradient
   \f[
       J = \det \mathbf{F}
         = \mathbf{1} + \frac{\partial \mathbf{u}}{\partial \mathbf{X}}
   \f]
 * is equal to the dilatation. However, when the displacement
 * field \f$ \mathbf{u} \f$ is restricted to a finite dimensional
 * representation, it may not contain enough degrees of freedom to
 * represent nearly incompressible deformations without making the
 * response overly stiff. Therefore, the dilatation and pressure are
 * represented as separate fields. The modified deformation gradient
 * is given by
   \f[
       \bar{\mathbf{F}} = \left( \frac{\Theta}{J} \right)^{1/3} \mathbf{F}.
   \f]
 * The remainder of the formulation results as a consequence. For Q1P0,
 * the (2D) bi- or (3D) trilinear displacement field (Q1) is combined
 * with piecewise constant (P0) pressure and dilatation fields. Since the
 * space for these fields is restricted to within element domains, these
 * degrees of freedom can be determined analytically at the element level
 * and substituted into the remaining element equations to results in
 * a purely displacement-based formulation.
 *
 * \note Several errors appear in the derivation
 * of the consistent tangent in the CMAME paper. Therefore,
 * the implementation of the tangent here does not match the
 * published formulation. */
class SimoQ1P0: public UpdatedLagrangianT
{
public:

	/** constructor */
	SimoQ1P0(const ElementSupportT& support);

	/** destructor */
	virtual ~SimoQ1P0();

	/** finalize current step - step is solved */
	virtual void CloseStep(void);

	/** restore last converged state */
	virtual GlobalT::RelaxCodeT ResetStep(void);

	/** read restart information from stream */
	virtual void ReadRestart(istream& in);

	/** write restart information from stream */
	virtual void WriteRestart(ostream& out) const;

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:
	virtual void SetShape();
	//
	virtual void SetLocalArrays();

	/** form shape functions and derivatives */
	virtual void SetGlobalShape(void);

	/** form the element stiffness matrix */
	virtual void FormStiffness(double constK);

	/** calculate the internal force contribution ("-k*d") */
	virtual void FormKd(double constK);

private:

	/** compute mean shape function gradient, H (reference volume), and
	 * current element volume, equation (2.20) */
	void SetMeanGradient(dArray2DT& mean_gradient, double& H, double& v) const;

	/** special mixed index term in the tangent. Needed to compute
	 * the term in the tangent resulting from
	 * \f$ \nabla \mathbf{u} \textrm{:} \left( \nabla \boldsymbol{\eta} \right)^T \f$.
	 */
	void bSp_bRq_to_KSqRp(const dMatrixT& b, dMatrixT& K) const;

	void MassMatrix();

	dSymMatrixT MaxwellStress(dArrayT E, const double epsilon);


protected:

	// Electric field stuff
	ArrayT<dArrayT> fE_List;
	dArrayT fE_all;

	dMatrixT fMaxwell;

	/** \name element volume */
	/*@{*/
	/** deformed element volume */
	dArrayT fElementVolume;

	/** deformed element volume from the last time step */
	dArrayT fElementVolume_last;
	/*@}*/

	/** element pressure. Calculated during SimoQ1P0::FormKd. */
	dArrayT fPressure;

	/** determinant of the deformation gradient for the current element */
	dArrayT fJacobian;

	/** \name work space */
	/*@{*/
	dArray2DT fMeanGradient; /**< mean gradient over element */
	dMatrixT fF_tmp; /**< F workspace */
	dMatrixT fNEEmat; /**< dimension of stiffness matrix */
	dMatrixT fdiff_b;
	dMatrixT fb_bar;
	dMatrixT fb_sig;

private:
	// For voltage field
	LocalArrayT fLocScalarPotential;

	dMatrixT fAmm_geo;
	dMatrixT fAmm_mat;
  	dMatrixT fMassMatrix;	// mass matrix for LHS

	dMatrixT fF_mech;
	dSymMatrixT D;
	// For case of electric field coupling for DE
  	const FieldT* fElectricScalarPotentialField;


	/*@}*/
};

} // namespace Tahoe
#endif /* _SIMO_Q1_P0_H_ */
