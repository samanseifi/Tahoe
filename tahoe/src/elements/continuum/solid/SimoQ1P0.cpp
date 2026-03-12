/* $Id: SimoQ1P0.cpp,v 1.15 2017/06/18 14:40:12 samanseifi Exp $ */
/* $Id: SimoQ1P0.cpp,v 1.14 2009/05/21 22:30:27 tdnguye Exp $ */
#include "SimoQ1P0.h"

#include "ShapeFunctionT.h"
#include "SolidMaterialT.h"
#include "SolidMatListT.h"

// #include "FSDEMatQ1P0T.h"
#include "FSSolidMatT.h"

#include "incQ1P0.h"

#include "TensorTransformT.h"

#include <iostream>

using namespace Tahoe;

/* constructor */
SimoQ1P0::SimoQ1P0(const ElementSupportT& support):
	UpdatedLagrangianT(support),
	fLocScalarPotential(LocalArrayT::kESP),
	fElectricScalarPotentialField(0),
	fElectricPermittivity(1.0),
	fPhaseFieldField(NULL),
	fLocPhaseField(LocalArrayT::kDisp),
	fPF_kSmall(1.0e-6)
{
	SetName("updated_lagrangian_Q1P0");
}

/* destructor */
SimoQ1P0::~SimoQ1P0()
{
		delete fCurrShapes;
		fCurrShapes = NULL;

}

/* finalize current step - step is solved */
void SimoQ1P0::CloseStep(void)
{
	/* inherited */
	UpdatedLagrangianT::CloseStep();

	/* store converged solution */
	fElementVolume_last = fElementVolume;
}

/* restore last converged state */
GlobalT::RelaxCodeT SimoQ1P0::ResetStep(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = UpdatedLagrangianT::ResetStep();

	/* store converged solution */
	fElementVolume = fElementVolume_last;

	return relax;
}

/* read restart information from stream */
void SimoQ1P0::ReadRestart(istream& in)
{
	/* inherited */
	UpdatedLagrangianT::ReadRestart(in);

	/* read restart data */
	in >> fElementVolume;

	/* reset last state */
	fElementVolume_last = fElementVolume;
}

/* write restart information from stream */
void SimoQ1P0::WriteRestart(ostream& out) const
{
	/* inherited */
	UpdatedLagrangianT::WriteRestart(out);

	/* read restart data */
	out << fElementVolume << '\n';
}

/* describe the parameters needed by the interface */
void SimoQ1P0::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	UpdatedLagrangianT::DefineParameters(list);

	/* electric permittivity for Maxwell stress and electrical tangent */
	ParameterT epsilon(fElectricPermittivity, "epsilon");
	epsilon.SetDefault(1.0);
	list.AddParameter(epsilon, ParameterListT::ZeroOrOnce);

	/* phase-field residual stiffness (optional) */
	ParameterT pf_k(fPF_kSmall, "pf_residual_stiffness");
	pf_k.SetDefault(1.0e-6);
	list.AddParameter(pf_k, ParameterListT::ZeroOrOnce);
}

/* accept parameter list */
void SimoQ1P0::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "SimoQ1P0::TakeParameterList";

	/* read epsilon before calling inherited (which may call SetGlobalShape) */
	const ParameterT* peps = list.Parameter("epsilon");
	fElectricPermittivity = peps ? double(*peps) : 1.0;

	/* read phase-field residual stiffness */
	const ParameterT* ppf_k = list.Parameter("pf_residual_stiffness");
	fPF_kSmall = ppf_k ? double(*ppf_k) : 1.0e-6;

	/* inherited */
	UpdatedLagrangianT::TakeParameterList(list);

	// Check if there is electric field coupling (For DE models)
	fElectricScalarPotentialField = ElementSupport().Field("electric_scalar_potential");
	if (!fElectricScalarPotentialField) {
	  std::cout << "There is no electric field coupling. Perhaps it's not DE model" << std::endl;
	}

	/* check for phase-field coupling */
	fPhaseFieldField = ElementSupport().Field("phase_field");
	if (fPhaseFieldField) {
		std::cout << "SimoQ1P0: phase-field fracture coupling detected." << std::endl;
	} else {
		std::cout << "SimoQ1P0: no phase-field coupling." << std::endl;
	}

	/* check geometry code and number of element nodes -> Q1 */
	if (GeometryCode() == GeometryT::kQuadrilateral) {
		if (NumElementNodes() != 4)
			ExceptionT::BadInputValue(caller, "expecting 4 node quad: %d", NumElementNodes());
	}
	else if (GeometryCode() == GeometryT::kHexahedron) {
		if (NumElementNodes() != 8)
			ExceptionT::BadInputValue(caller, "expecting 8 node hex: %d", NumElementNodes());
	}
	else
		ExceptionT::BadInputValue(caller, "expecting hex or quad geometry: %d", GeometryCode());

	/* need to store last deformed element volume */
	fElementVolume.Dimension(NumElements());
	fElementVolume = 0.0;
	fElementVolume_last.Dimension(NumElements());
	fElementVolume_last = 0.0;

	int nen = NumElementNodes();
	int nsd = NumSD();
	int nme = nen * nsd;	// # of mechanical DOFs per element

	/* getting ready for electric field calculations */
	const int nip = NumIP();
	fE_all.Dimension(nip*nsd);
	fE_all = 0.0;	// testing HSP
	fE_List.Dimension(nip);
	// Neccessary:
	for (int i = 0; i < nip; ++i) {
		fE_List[i].Alias(nsd, fE_all.Pointer(i * nsd));
	}
	fF_mech.Dimension(nsd);
	fF_mech.Identity();

	D.Dimension(2);

	fAmm_mat.Dimension(nme, nme);
	fAmm_geo.Dimension(nen, nen);	// dimensions changed for Q1P0!
	fMassMatrix.Dimension(nme, nme);

	/* electrical tangent in Voigt notation: 6x6 for 3D, 3x3 for 2D */
	int voigt_dim = (nsd == 3) ? 6 : 3;
	fTangentMechanical.Dimension(voigt_dim, voigt_dim);
	fTangentMechanical = 0.0;

	/* tensor transform for push-forward */
	fTransform.Dimension(nsd);

	/* phase-field IP storage */
	fPhaseField_ip.Dimension(nip);
	fPhaseField_ip = 0.0;
	fDegradation_ip.Dimension(nip);
	fDegradation_ip = 1.0;  /* no degradation by default */

	/* element pressure */
	fPressure.Dimension(NumElements());
	fPressure = 0.0;

	/* determinant of the deformation gradient */
	fJacobian.Dimension(NumIP());
	fJacobian = 1.0;

	/* dimension work space */
	fMeanGradient.Dimension(NumSD(), NumElementNodes());
	fNEEmat.Dimension(fLHS);
	fdiff_b.Dimension(fGradNa);
	fb_bar.Dimension(fGradNa);
	fb_sig.Dimension(fGradNa);
	fF_tmp.Dimension(NumSD());

	/* need to initialize previous volume */
	Top();
	while (NextElement())
	{
		/* inherited - computes gradients and standard
		 * deformation gradients */
		UpdatedLagrangianT::SetGlobalShape();

		/* compute mean of shape function gradients */
		double H; /* reference volume */
		double& v = fElementVolume_last[CurrElementNumber()];
		SetMeanGradient(fMeanGradient, H, v);
	}
}

/***********************************************************************
 * Protected
 ***********************************************************************/
/* form shape functions and derivatives */
void SimoQ1P0::SetGlobalShape(void)
{
	/* current element number */
	int elem = CurrElementNumber();

	/* inherited - computes gradients and standard
	 * deformation gradients */
	FiniteStrainT::SetGlobalShape();

	/* calculating electric field — only when the field is present */
	if (fElectricScalarPotentialField) {
		SetLocalU(fLocScalarPotential);
		for (int i = 0; i < NumIP(); i++) {
			dArrayT& E = fE_List[i];
			dMatrixT E1(1, NumSD());
			fShapes->GradU(fLocScalarPotential, E1, i);
			E1 *= -1.0;
			for (int j = 0; j < NumSD(); j++)
				E[j] = E1(0,j);
		}
	} else {
		/* no electric field — keep E = 0 at all integration points */
		fE_all = 0.0;
	}


	/* interpolate phase-field d at integration points */
	if (fPhaseFieldField) {
		SetLocalU(fLocPhaseField);
		dArrayT d_vec(1);
		for (int i = 0; i < NumIP(); i++) {
			fShapes->InterpolateU(fLocPhaseField, d_vec, i);
			double d = d_vec[0];
			/* clamp to [0, 1] */
			if (d < 0.0) d = 0.0;
			if (d > 1.0) d = 1.0;
			fPhaseField_ip[i] = d;
			fDegradation_ip[i] = (1.0 - d)*(1.0 - d) + fPF_kSmall;
		}
	} else {
		fDegradation_ip = 1.0;
	}

	/* shape function wrt current config */
	SetLocalX(fLocCurrCoords);
	fCurrShapes->SetDerivatives();

	/* compute mean of shape function gradients */
	double H; /* reference volume */
	double& v = fElementVolume[elem];
	SetMeanGradient(fMeanGradient, H, v);

	/* last deformed volume */
	double& v_last = fElementVolume_last[elem];

	/* what needs to get computed */
	int material_number = CurrentElement().MaterialNumber();
	bool needs_F = Needs_F(material_number);
	bool needs_F_last = Needs_F_last(material_number);

	/* B-bar exponent: 1/nsd ensures det(F_bar) = Theta (constant within element) */
	double bbar_exp = 1.0/NumSD();

	/* loop over integration points */
	for (int i = 0; i < NumIP(); i++)
	{
		/* deformation gradient */
		if (needs_F)
		{
			/* "replace" dilatation: F_bar = (Theta/J)^{1/nsd} * F */
			dMatrixT& F = fF_List[i];
			double J = F.Det();
			F *= pow(v/(H*J), bbar_exp);

			/* store Jacobian */
			fJacobian[i] = J;
		}

		/* "last" deformation gradient */
		if (needs_F_last)
		{
			/* "replace" dilatation */
			dMatrixT& F = fF_last_List[i];

			double J = F.Det();
			F *= pow(v_last/(H*J), bbar_exp);
		}
	}
}

/* initialization functions */
void SimoQ1P0::SetShape(void)
{
	/* inherited */
	FiniteStrainT::SetShape();

	/* linked shape functions */
	fCurrShapes = new ShapeFunctionT(*fShapes, fLocCurrCoords);
	if (!fCurrShapes) throw ExceptionT::kOutOfMemory ;

	fCurrShapes->Initialize();
}

/* bringing electrical model here */
void SimoQ1P0::SetLocalArrays()
{
	/* look for electric scalar potential field — optional */
	if (0 == fElectricScalarPotentialField)
		fElectricScalarPotentialField = ElementSupport().Field("electric_scalar_potential");

	/* inherited */
	FiniteStrainT::SetLocalArrays();

	/* allocate and set source */
	fLocCurrCoords.Dimension(NumElementNodes(), NumSD());
	ElementSupport().RegisterCoordinates(fLocCurrCoords);

	/* only register potential if the field exists */
	if (fElectricScalarPotentialField) {
		const int nen = NumElementNodes();
		fLocScalarPotential.Dimension(nen, 1);
		fElectricScalarPotentialField->RegisterLocal(fLocScalarPotential);
	}

	/* register phase-field local array */
	if (0 == fPhaseFieldField)
		fPhaseFieldField = ElementSupport().Field("phase_field");
	if (fPhaseFieldField) {
		const int nen = NumElementNodes();
		fLocPhaseField.Dimension(nen, 1);
		fPhaseFieldField->RegisterLocal(fLocPhaseField);
	}
}


/* recompute deformation gradients and B-bar state from current fLocDisp.
 * Call after perturbing displacement for numerical tangent. */
void SimoQ1P0::RecomputeDeformationState(void)
{
	const int nsd = NumSD();
	const int nip = NumIP();
	const int elem = CurrElementNumber();

	/* recompute current-config shape function derivatives.
	 * fCurrShapes is linked to fLocCurrCoords, which must be
	 * updated before this call. */
	fCurrShapes->SetDerivatives();

	/* recompute mean gradient and element volume */
	double H;
	SetMeanGradient(fMeanGradient, H, fElementVolume[elem]);

	/* recompute deformation gradients with B-bar correction */
	double bbar_exp = 1.0/nsd;
	for (int ip = 0; ip < nip; ip++)
	{
		dMatrixT& F = fF_List[ip];
		fShapes->GradU(fLocDisp, F, ip);
		F.PlusIdentity();

		double J = F.Det();
		fJacobian[ip] = J;
		F *= pow(fElementVolume[elem]/(H*J), bbar_exp);
	}
}

/* compute element internal force into provided vector.
 * Same physics as FormKd but writes to 'force' instead of fRHS. */
void SimoQ1P0::ComputeInternalForce(double constK, dArrayT& force)
{
	force = 0.0;

	const int elem = CurrElementNumber();
	double p_bar_tmp = 0.0;

	const double* Det    = fCurrShapes->IPDets();
	const double* Weight = fCurrShapes->IPWeights();

	/* mean B-bar gradient */
	fCurrShapes->GradNa(fMeanGradient, fb_bar);
	fShapes->TopIP();

	while (fShapes->NextIP())
	{
		/* B-bar strain-displacement matrix */
		Set_B_bar(fCurrShapes->Derivatives_U(), fMeanGradient, fB);

		/* Cauchy stress from base material */
		dSymMatrixT cauchy = fCurrMaterial->s_ij();

		/* Maxwell stress (if electrical coupling) */
		if (fElectricScalarPotentialField) {
			const dArrayT  E = fE_List[CurrIP()];
			const dMatrixT F = DeformationGradient();
			dSymMatrixT maxwell = s_electric_ij(E, F, fElectricPermittivity);
			cauchy += maxwell;
		}

		/* phase-field degradation */
		double g_d = fDegradation_ip[CurrIP()];
		cauchy *= g_d;

		/* B^T * sigma */
		fB.MultTx(cauchy, fNEEvec);

		/* J_bar / J correction */
		double J_bar = DeformationGradient().Det();
		double J_correction = J_bar/fJacobian[CurrIP()];

		/* integrate pressure */
		p_bar_tmp += (*Weight)*(*Det)*J_correction*cauchy.Trace()/NumSD();

		/* accumulate force */
		force.AddScaled(constK*(*Weight++)*(*Det++)*J_correction, fNEEvec);
	}

	/* volume-averaged pressure */
	fPressure[elem] = p_bar_tmp/fElementVolume[elem];
}

/* form the element stiffness matrix — numerical tangent via
 * forward finite differences of the internal force.
 * This gives the exact consistent tangent for the Q1P0
 * formulation, including all B-bar correction terms. */
void SimoQ1P0::FormStiffness(double constK)
{
	const int nsd = NumSD();
	const int nen = NumElementNodes();
	const int ndof_elem = nen*nsd;
	const int nip = NumIP();
	const int elem = CurrElementNumber();
	const double eps = 1.0e-7;

	/* ---- save state ---- */

	/* local displacements */
	dArrayT u_save(fLocDisp.Length());
	for (int k = 0; k < fLocDisp.Length(); k++)
		u_save[k] = fLocDisp[k];

	/* current coordinates */
	dArrayT x_save(fLocCurrCoords.Length());
	for (int k = 0; k < fLocCurrCoords.Length(); k++)
		x_save[k] = fLocCurrCoords[k];

	/* deformation gradients */
	ArrayT<dMatrixT> F_save(nip);
	for (int ip = 0; ip < nip; ip++) {
		F_save[ip].Dimension(nsd, nsd);
		F_save[ip] = fF_List[ip];
	}

	/* Jacobians, mean gradient, volume, pressure */
	dArrayT J_save(nip);
	J_save = fJacobian;
	dArray2DT mg_save(fMeanGradient);
	double v_save = fElementVolume[elem];
	double p_save = fPressure[elem];

	/* ---- reference internal force ---- */
	dArrayT f0(ndof_elem);
	ComputeInternalForce(constK, f0);

	/* ---- perturb each DOF and finite-difference ---- */
	dArrayT f_pert(ndof_elem);

	for (int a = 0; a < nen; a++) {
		for (int i = 0; i < nsd; i++) {
			/* column index in tangent matrix (equation ordering) */
			int col = a*nsd + i;

			/* perturb displacement and current coordinates */
			fLocDisp(a, i) += eps;
			fLocCurrCoords(a, i) += eps;

			/* recompute all deformation state */
			RecomputeDeformationState();

			/* compute perturbed internal force */
			ComputeInternalForce(constK, f_pert);

			/* forward finite difference → column of tangent */
			for (int j = 0; j < ndof_elem; j++)
				fLHS(j, col) += (f_pert[j] - f0[j])/eps;

			/* undo perturbation */
			fLocDisp(a, i) -= eps;
			fLocCurrCoords(a, i) -= eps;
		}
	}

	/* ---- restore full state ---- */
	for (int k = 0; k < fLocDisp.Length(); k++)
		fLocDisp[k] = u_save[k];
	for (int k = 0; k < fLocCurrCoords.Length(); k++)
		fLocCurrCoords[k] = x_save[k];
	for (int ip = 0; ip < nip; ip++)
		fF_List[ip] = F_save[ip];
	fJacobian = J_save;
	fMeanGradient = mg_save;
	fElementVolume[elem] = v_save;
	fPressure[elem] = p_save;
	fCurrShapes->SetDerivatives();
}

/* calculate the internal force contribution ("-k*d") */
void SimoQ1P0::FormKd(double constK)
{
	const double* Det    = fCurrShapes->IPDets();
	const double* Weight = fCurrShapes->IPWeights();

	/* collect incremental heat */
	bool need_heat = fElementHeat.Length() == fShapes->NumIP();

	/* current element number */
	int elem = CurrElementNumber();

	/* constant pressure */
	double& p_bar = fPressure[elem];
	p_bar = 0.0;

	fCurrShapes->TopIP();
	while ( fCurrShapes->NextIP() )
	{
		/* strain displacement matrix */
		Set_B_bar(fCurrShapes->Derivatives_U(), fMeanGradient, fB);

		/* B^T * Cauchy stress */
		dSymMatrixT cauchy = fCurrMaterial->s_ij();

		const dArrayT  E 		= fE_List[CurrIP()];
		const dMatrixT F 		= DeformationGradient();
		dSymMatrixT maxwell = s_electric_ij(E, F, fElectricPermittivity);
		cauchy += maxwell;

		/* apply phase-field degradation */
		double g_d = fDegradation_ip[CurrIP()];
		cauchy *= g_d;

		fB.MultTx(cauchy, fNEEvec);

		/* determinant of modified deformation gradient */
		double J_bar = DeformationGradient().Det();

		/* detF correction */
		double J_correction = J_bar/fJacobian[CurrIP()];

		/* integrate pressure */
		p_bar += (*Weight)*(*Det)*J_correction*cauchy.Trace()/NumSD();

		/* accumulate */
		fRHS.AddScaled(constK*(*Weight++)*(*Det++)*J_correction, fNEEvec);

		/* incremental heat generation */
		if (need_heat)
			fElementHeat[fShapes->CurrIP()] += fCurrMaterial->IncrementalHeat();
	}

	/* volume averaged */
	p_bar /= fElementVolume[CurrElementNumber()];
}

/* add Maxwell stress (with phase-field degradation) to output stress */
void SimoQ1P0::AddExtraStress(dSymMatrixT& cauchy) const
{
	if (!fElectricScalarPotentialField) return;

	/* Maxwell stress at current IP */
	const dArrayT  E = fE_List[CurrIP()];
	const dMatrixT F = DeformationGradient();

	/* NOTE: s_electric_ij is non-const because it uses workspace fStress.
	 * Use const_cast here since this is purely a query operation. */
	dSymMatrixT maxwell =
		const_cast<SimoQ1P0*>(this)->s_electric_ij(E, F, fElectricPermittivity);

	/* apply phase-field degradation */
	double g_d = fDegradation_ip[CurrIP()];
	maxwell *= g_d;

	cauchy += maxwell;
}


/***********************************************************************
 * Private
 ***********************************************************************/

void SimoQ1P0::MassMatrix()
{
	/* Calculate mass matrix for mechanical portion of LHS for FormStiffness */
	/* Implement lumped mass matrix only - better convergence for implicit dynamics
		as compared to consistent mass - see Hughes FEM book */
	int nen = NumElementNodes();
	int ndof = NumDOF();
	int nsd = NumSD();
    	int nme = nen * nsd;	// # of mechanical DOFs per element
	dArrayT NEEvec(nme);
	NEEvec = 0.0;
	double dsum = 0.0;
	double totmas = 0.0;
	fMassMatrix = 0.0;

	fShapes->TopIP();
	while (fShapes->NextIP() != 0) {

		/* integration factor - ignoring constM factor */
		double temp1 = fShapes->IPDet() * fShapes->IPWeight();
//		if (ip_weight) temp1 *= *ip_weight++;

		const double* Na = fShapes->IPShapeU();
		totmas += temp1;

		for (int lnd = 0; lnd < nen; lnd++) {
			double temp2 = temp1*Na[lnd]*Na[lnd];
			dsum += temp2;
			NEEvec[lnd] += temp2;
		}
	}

	/* scale diagonal to conserve total mass */
	double diagmass = totmas/dsum;

	/* lump mass onto diagonal */
	double* pmass = fMassMatrix.Pointer();
	int inc = fMassMatrix.Rows() + 1;
	for (int lnd = 0; lnd < nen; lnd++)
	{
		double temp = diagmass*NEEvec[lnd];
		for (int ed = 0; ed < ndof; ed++)
		{
			*pmass += temp;
			pmass += inc;
		}
	}

}

/* compute mean shape function gradient, Hughes (4.5.23) */
void SimoQ1P0::SetMeanGradient(dArray2DT& mean_gradient, double& H, double& v) const
{
	/* assume same integration rule defined for current and references
	 * shape functions */
	int nip = NumIP();
	const double*   det = fCurrShapes->IPDets();
	const double* det_0 = fShapes->IPDets();
	const double*     w = fShapes->IPWeights();

	/* H and current volume */
	H = 0.0;
	v = 0.0;

	for (int i = 0; i < nip; i++)
	{
		H += w[i]*det_0[i];
		v += w[i]*det[i];
	}

	/* initialize */
	mean_gradient = 0.0;

	/* integrate */
	for (int i = 0; i < nip; i++)
		mean_gradient.AddScaled(w[i]*det[i]/v, fCurrShapes->Derivatives_U(i));
}

void SimoQ1P0::bSp_bRq_to_KSqRp(const dMatrixT& b, dMatrixT& K) const
{
#if __option(extended_errorcheck)
	/* dimension check */
	if (b.Length() != K.Rows() ||
	    K.Rows() != K.Cols()) ExceptionT::SizeMismatch("SimoQ1P0::bSp_bRq_to_KSqRp");
#endif

	int dim = K.Rows();
	int sub_dim = b.Rows();
	int S = 0;
	int p = 0;
	for (int i = 0; i < dim; i++)
	{
		int R = 0;
		int q = 0;
		for (int j = 0; j < dim; j++)
		{
			K(i,j) = b(q,S)*b(p,R);

			q++;
			if (q == sub_dim) {
				R++;
				q = 0;
			}
		}
		p++;
		if (p == sub_dim) {
			S++;
			p = 0;
		}
	}
}


/* Calculating the total stress (elec+mech) in reference conf and then push it forward! */
dSymMatrixT SimoQ1P0::s_electric_ij(const dArrayT E, const dMatrixT F, const double epsilon)
{
	int nsd = NumSD();

	dSymMatrixT fStress(nsd);

	dMatrixT stress_temp(nsd);
	dMatrixT stress_temp2(nsd);

	if (nsd == 3){
		dMatrixT C(nsd); // Right Cauchy strain tensor

		C.MultATB(F, F);
		double J = F.Det();
		double I1 = C(0,0) + C(1,1) + C(2,2);

		stress_temp2 = 0.0;

		me_pk2_q1p0(epsilon, E.Pointer(), C.Pointer(), F.Pointer(), J, stress_temp2.Pointer());

		fStress.FromMatrix(stress_temp2);

	}	else if (nsd == 2) {

		dMatrixT F2D(nsd);	// Deformation Gradient
		dMatrixT C2D(nsd);	// Right Cauchy strain tensor

		F2D = F;
		C2D.MultATB(F, F);
		double J = F2D.Det();
		double det_C = C2D.Det();

		dMatrixT C3D(3), F3D(3), stress_temp(3), stress_temp2(3);
  		dArrayT E3D(3);

		stress_temp2 = 0.0;

		C3D[0] = C2D[0];
		C3D[1] = C2D[1];
		C3D[2] = 0.0;

		C3D[3] = C2D[2];
		C3D[4] = C2D[3];
		C3D[5] = 0.0;

		C3D[6] = 0.0;
		C3D[7] = 0.0;
		C3D[8] = 1.0;

		F3D[0] = F2D[0];
		F3D[1] = F2D[1];
		F3D[2] = 0.0;

		F3D[3] = F2D[2];
		F3D[4] = F2D[3];
		F3D[5] = 0.0;

		F3D[6] = 0.0;
		F3D[7] = 0.0;
		F3D[8] = 1.0;

		E3D[0] = E[0];
		E3D[1] = E[1];
		E3D[2] = 0.0;

		double I1 = C2D(0, 0) + C2D(1, 1) + 1.0/det_C;

		/* call C function for electrical part of PK2 stress */
		me_pk2_q1p0(epsilon, E3D.Pointer(), C3D.Pointer(), F3D.Pointer(), J, stress_temp2.Pointer());

		/* use stress_temp2 directly — stress_temp was uninitialized (bug fix) */
		fStress(0,0) = stress_temp2(0,0);
    	fStress(0,1) = stress_temp2(0,1);
		fStress(1,0) = stress_temp2(1,0);
    	fStress(1,1) = stress_temp2(1,1);
	}


	/* Getting ready to push forward the PK2 stress to Cauchy */
	dSymMatrixT S = fStress;
	double J = F.Det();

	fStress.MultQBQT(F, S);		// Push forward procedure
	fStress *= 1.0 / J;

	return fStress;

}

/* Correct 4th-order Voigt push-forward: c_spatial = (1/J) * P * C_mat * P^T
 *
 * Voigt ordering (0-based indices into the matrix):
 *   3D: (0,0)->0, (1,1)->1, (2,2)->2, (1,2)->3, (0,2)->4, (0,1)->5
 *   2D: (0,0)->0, (1,1)->1, (0,1)->2
 *
 * P[a,b] = F[vi[a], vi[b]] * F[vj[a], vj[b]]
 * where (vi[a], vj[a]) is the index pair for Voigt index a.
 */
dMatrixT SimoQ1P0::c_electrical_ijkl(const dArrayT E, const dMatrixT F, const double epsilon)
{
    int nsd = NumSD();
    double J = F.Det();

    dMatrixT C_mat_3D(6);
    C_mat_3D = 0.0;

    if (nsd == 3) {
        dMatrixT C(3);
        C.MultATB(F, F);
        me_tanmod_q1p0(epsilon, E.Pointer(), C.Pointer(), F.Pointer(), J, C_mat_3D.Pointer());

        /* 4th-order Voigt push-forward (3D):
         * Voigt pairs: 0=(0,0), 1=(1,1), 2=(2,2), 3=(1,2), 4=(0,2), 5=(0,1)
         * T[a,b] = F[vi[a],vi[b]]*F[vj[a],vj[b]]  (diagonal b)
         * T[a,b] += F[vi[a],vj[b]]*F[vj[a],vi[b]]  (off-diagonal b) */
        static const int vi3D[6] = {0, 1, 2, 1, 0, 0};
        static const int vj3D[6] = {0, 1, 2, 2, 2, 1};

        dMatrixT T(6);
        for (int a = 0; a < 6; a++)
            for (int b = 0; b < 6; b++) {
                T(a, b) = F(vi3D[a], vi3D[b]) * F(vj3D[a], vj3D[b]);
                if (vi3D[b] != vj3D[b])
                    T(a, b) += F(vi3D[a], vj3D[b]) * F(vj3D[a], vi3D[b]);
            }

		fTangentMechanical.SetToScaled(1.0 / J, fTransform.PushForward(F, C_mat_3D)); // finite difference c_ijkl

    }
    else if (nsd == 2) {
        dMatrixT F2D(nsd);
        dMatrixT C2D(nsd);
        F2D = F;
        C2D.MultATB(F, F);
        double J2D = F2D.Det();

        /* extend to 3D for the auto-generated C function */
        dMatrixT C3D(3), F3D(3);
        dArrayT E3D(3);

        C3D[0]=C2D[0]; C3D[1]=C2D[1]; C3D[2]=0.0;
        C3D[3]=C2D[2]; C3D[4]=C2D[3]; C3D[5]=0.0;
        C3D[6]=0.0;    C3D[7]=0.0;    C3D[8]=1.0;

        F3D[0]=F2D[0]; F3D[1]=F2D[1]; F3D[2]=0.0;
        F3D[3]=F2D[2]; F3D[4]=F2D[3]; F3D[5]=0.0;
        F3D[6]=0.0;    F3D[7]=0.0;    F3D[8]=1.0;

        E3D[0]=E[0]; E3D[1]=E[1]; E3D[2]=0.0;

        me_tanmod_q1p0(epsilon, E3D.Pointer(), C3D.Pointer(), F3D.Pointer(), J2D, C_mat_3D.Pointer());

        /* extract plane-strain 2D block (Voigt 3D ordering: 11->0, 22->1, 12->5) */
        dMatrixT C_mat_2D(3);
        C_mat_2D(0,0)=C_mat_3D(0,0); C_mat_2D(0,1)=C_mat_3D(0,1); C_mat_2D(0,2)=C_mat_3D(0,5);
        C_mat_2D(1,0)=C_mat_3D(1,0); C_mat_2D(1,1)=C_mat_3D(1,1); C_mat_2D(1,2)=C_mat_3D(1,5);
        C_mat_2D(2,0)=C_mat_3D(5,0); C_mat_2D(2,1)=C_mat_3D(5,1); C_mat_2D(2,2)=C_mat_3D(5,5);

        /* 4th-order Voigt push-forward (2D):
         * Voigt pairs: 0=(0,0), 1=(1,1), 2=(0,1)
         * T[a,b] includes cross-term for off-diagonal Voigt index b */
        static const int vi2D[3] = {0, 1, 0};
        static const int vj2D[3] = {0, 1, 1};

        dMatrixT T2D(3);
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++) {
                T2D(a, b) = F2D(vi2D[a], vi2D[b]) * F2D(vj2D[a], vj2D[b]);
                if (vi2D[b] != vj2D[b])
                    T2D(a, b) += F2D(vi2D[a], vj2D[b]) * F2D(vj2D[a], vi2D[b]);
            }


		fTangentMechanical.SetToScaled(1.0 / J2D, fTransform.PushForward(F2D, C_mat_2D)); // finite difference c_ijkl

    }

    return fTangentMechanical;
}
