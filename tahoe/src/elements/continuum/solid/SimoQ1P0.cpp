/* $Id: SimoQ1P0.cpp,v 1.14 2009/05/21 22:30:27 tdnguye Exp $ */
#include "SimoQ1P0.h"

#include "ShapeFunctionT.h"
#include "SolidMaterialT.h"
#include "SolidMatListT.h"

using namespace Tahoe;

/* constructor */
SimoQ1P0::SimoQ1P0(const ElementSupportT& support):
	UpdatedLagrangianT(support), fLocScalarPotential(LocalArrayT::kESP), fElectricScalarPotentialField(0)
{
	SetName("updated_lagrangian_Q1P0");
}

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

/* accept parameter list */
void SimoQ1P0::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "SimoQ1P0::TakeParameterList";


	/* inherited */
	UpdatedLagrangianT::TakeParameterList(list);

	// Check if there is electric field coupling (For DE models)
	fElectricScalarPotentialField = ElementSupport().Field("electric_scalar_potential");
	if (!fElectricScalarPotentialField) {
	  std::cout << "There is no electric field coupling. Perhaps it's not DE model" << std::endl;
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

    const int nip = NumIP();
    fE_all.Dimension(nip*nsd);
    fE_all = 0.0;	// testing HSP
    fE_List.Dimension(nip);

	// Neccessary
    for (int i = 0; i < nip; ++i) {
	    fE_List[i].Alias(nsd, fE_all.Pointer(i * nsd));
    }

    fAmm_mat.Dimension(nme, nme);
    fAmm_geo.Dimension(nen, nen);	// dimensions changed for Q1P0!

    fMassMatrix.Dimension(nme, nme);

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

	SetLocalU(fLocScalarPotential);


	/* inherited - computes gradients and standard
	 * deformation gradients */
	FiniteStrainT::SetGlobalShape();


	for (int i = 0; i < NumIP(); i++) {
        	// electric field
   	 	dArrayT& E = fE_List[i];
   	     	dMatrixT E1(1, NumSD());

   	     	fShapes->GradU(fLocScalarPotential, E1, i);
   	     	E1 *= -1.0;

   	     	for (int j = 0; j < NumSD(); j++)
   		     	E[j] = E1(0,j);
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

	/* loop over integration points */
	for (int i = 0; i < NumIP(); i++)
	{
		/* deformation gradient */
		if (needs_F)
		{
			/* "replace" dilatation */
			dMatrixT& F = fF_List[i];
			double J = F.Det();
			F *= pow(v/(H*J), 1.0/3.0);

			/* store Jacobian */
			fJacobian[i] = J;
		}

		/* "last" deformation gradient */
		if (needs_F_last)
		{
			/* "replace" dilatation */
			dMatrixT& F = fF_last_List[i];

			double J = F.Det();
			F *= pow(v_last/(H*J), 1.0/3.0);
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

void SimoQ1P0::SetLocalArrays()
{

	// look for an electric scalar potential field
	const FieldT* esp = 0;

	if (0 == fElectricScalarPotentialField) {
	  	esp = ElementSupport().Field("electric_scalar_potential");
	  	fElectricScalarPotentialField = esp;
	} else {
	  	esp = fElectricScalarPotentialField;
	}

	std::cout << "where?" << std::endl;

	// if (0 == esp) {
	//
	//   	std::cout << std::endl;
	//   	std::cout << "SimoQ1P0::SetLocalArrays: ";
	//   	std::cout << "Voltage field not found.";
	//   	std::cout << std::endl;
	//   	throw ExceptionT::kGeneralFail;
	// }

	/* inherited */
	FiniteStrainT::SetLocalArrays();

	const int nen = NumElementNodes();

	fLocScalarPotential.Dimension(nen, 1);

	// Register fields
	esp->RegisterLocal(fLocScalarPotential);

	/* allocate and set source */
	fLocCurrCoords.Dimension(NumElementNodes(), NumSD());
	ElementSupport().RegisterCoordinates(fLocCurrCoords);
}


/* form the element stiffness matrix */
void SimoQ1P0::FormStiffness(double constK)
{
//	int order = fIntegrator->Order();

	/* matrix format */
	dMatrixT::SymmetryFlagT format =
		(fLHS.Format() == ElementMatrixT::kNonSymmetric) ?
		dMatrixT::kWhole :
		dMatrixT::kUpperOnly;

	/* current element info */
	int el = CurrElementNumber();
	double v = fElementVolume[el];
	double p_bar = fPressure[el];

	fAmm_mat = 0.0;
	fAmm_geo = 0.0;

	/* integration */
	const double* Det    = fCurrShapes->IPDets();
	const double* Weight = fCurrShapes->IPWeights();

	/* initialize */
	fStressStiff = 0.0;

	fCurrShapes->GradNa(fMeanGradient, fb_bar);
	fShapes->TopIP();
	while ( fShapes->NextIP() )
	{
		/* double scale factor */
		double scale = constK*(*Det++)*(*Weight++);

	/* S T R E S S   S T I F F N E S S */
		/* compute Cauchy stress */
		const dSymMatrixT& cauchy = fCurrMaterial->s_ij();
		cauchy.ToMatrix(fCauchyStress);

		/* determinant of modified deformation gradient */
		double J_bar = DeformationGradient().Det();

		/* detF correction */
		double J_correction = J_bar/fJacobian[CurrIP()];
		double p = J_correction*cauchy.Trace()/3.0;

		/* get shape function gradients matrix */
		fCurrShapes->GradNa(fGradNa);
		fb_sig.MultAB(fCauchyStress, fGradNa);

		/* integration constants */
		fCauchyStress *= scale*J_correction;

		fAmm_geo.MultQTBQ(fGradNa, fCauchyStress, format, dMatrixT::kAccumulate);

		/* using the stress symmetry */
		fStressStiff.MultQTBQ(fGradNa, fCauchyStress,
			format, dMatrixT::kAccumulate);

	/* M A T E R I A L   S T I F F N E S S */
		/* strain displacement matrix */
		Set_B_bar(fCurrShapes->Derivatives_U(), fMeanGradient, fB);

		/* get D matrix */
		fD.SetToScaled(scale*J_correction, fCurrMaterial->c_ijkl());

		/* accumulate */
		fAmm_mat.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);

/*----------------- TEMPORARLY COMMENTED TO CHANGE THE IMPLEMENTATION ---------------- */
		/* $div div$ term */
		/* fNEEmat.Outer(fGradNa, fGradNa);
		fLHS.AddScaled(p_bar*scale, fNEEmat);

		fdiff_b.DiffOf(fGradNa, fb_bar);
		fNEEmat.Outer(fdiff_b, fdiff_b);
		fLHS.AddScaled(scale*2.0*p/3.0, fNEEmat);

		fNEEmat.Outer(fb_sig, fdiff_b);
		fNEEmat.Symmetrize();
		fLHS.AddScaled(-J_correction*scale*4.0/3.0, fNEEmat);

		bSp_bRq_to_KSqRp(fGradNa, fNEEmat);
		fLHS.AddScaled(scale*(p - p_bar), fNEEmat); */
	}
	//MassMatrix();
	//fLHS.AddBlock(0, 0, fMassMatrix);
	fAmm_mat.Expand(fAmm_geo, NumDOF(), dMatrixT::kAccumulate);
	fLHS.AddBlock(0, 0, fAmm_mat);


	// stress stiffness into fLHS
	/* fLHS.Expand(fStressStiff, NumDOF(), dMatrixT::kAccumulate);

	// $\bar{div}\bar{div}$ term
	fNEEmat.Outer(fb_bar, fb_bar);
	fLHS.AddScaled(-p_bar*v, fNEEmat); */
/* ------------------------------------------------------------------------------------- */
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
		const dSymMatrixT& cauchy = fCurrMaterial->s_ij();
		fB.MultTx(cauchy, fNEEvec);

		/* determinant of modified deformation gradient */
		double J_bar = DeformationGradient().Det();

		/* detF correction */
		double J_correction = J_bar/fJacobian[CurrIP()];

		/* integrate pressure */
		p_bar += (*Weight)*(*Det)*J_correction*cauchy.Trace()/3.0;

		/* accumulate */
		fRHS.AddScaled(constK*(*Weight++)*(*Det++)*J_correction, fNEEvec);

		/* incremental heat generation */
		if (need_heat)
			fElementHeat[fShapes->CurrIP()] += fCurrMaterial->IncrementalHeat();
	}

	/* volume averaged */
	p_bar /= fElementVolume[CurrElementNumber()];
}

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

/***********************************************************************
 * Private
 ***********************************************************************/

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
