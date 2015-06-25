/* $Id: CurrMixtureSpeciesT.cpp,v 1.11 2006/06/26 21:20:01 thao Exp $ */
#include "CurrMixtureSpeciesT.h"
#include "UpdatedLagMixtureT.h"
#include "Q1P0MixtureT.h"
#include "ShapeFunctionT.h"
#include "NLDiffusionMaterialT.h"
#include "MaterialListT.h"
#include "KBC_CardT.h"

//DEBUG
#include "ofstreamT.h"

using namespace Tahoe;

/* constructor */
CurrMixtureSpeciesT::CurrMixtureSpeciesT(const ElementSupportT& support):
	NLDiffusionElementT(support),
	fGradientOption(kGlobalProjection),
	fConcentration(kCurrent),
	fOutputMass(false),
	fUpdatedLagMixture(NULL),
	fQ1P0Mixture(NULL),
	fBackgroundSpecies(NULL),
	fIndex(-1),
	fLocCurrCoords(LocalArrayT::kCurrCoords)
{
	SetName("current_mixture_species");
}

void CurrMixtureSpeciesT::SendOutput(int kincode)
{
	/* output flags */
	iArrayT flags(fNodalOutputCodes.Length());

	/* set flags to get desired output */
	flags = IOBaseT::kAtNever;
	switch (kincode)
	{
		case iNodalDisp:
		    flags[iNodalDisp] = 1;
			break;
		case iMaterialData:  /*fluxes*/
		    flags[iMaterialData] = 1;
			break;
		default:
			cout << "\n DiffusionElementT::SendKinematic: invalid output code: ";
			cout << kincode << endl;
	}

	/* number of output values */
	iArrayT n_counts;
	SetNodalOutputCodes(IOBaseT::kAtInc, flags, n_counts);

	/* reset averaging workspace */
	ElementSupport().ResetAverage(n_counts.Sum());

	/* no element output */
	iArrayT e_counts(fElementOutputCodes.Length());
	e_counts = 0;

	/* generate nodal values */
	dArray2DT e_values, n_values;
	ComputeOutput(n_counts, n_values, e_counts, e_values);
}

/* construct output labels array */
void CurrMixtureSpeciesT::SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
	iArrayT& counts) const
{
	/* initialize */
	counts.Dimension(flags.Length());
	counts = 0;

	if (flags[iNodalCoord] == mode)
		counts[iNodalCoord] = NumSD();
	if (flags[iNodalDisp] == mode)
		counts[iNodalDisp] = NumDOF();
	if (flags[iMaterialData] == mode)
		counts[iMaterialData ] = NumDOF()*NumSD();
}

void CurrMixtureSpeciesT::SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
	iArrayT& counts) const
{
#pragma unused(mode)
#pragma unused(flags)
	if (counts.Sum() != 0)
		ExceptionT::BadInputValue("DiffusionElementT::SetElementOutputCodes", "not implemented");
}

/* driver for calculating output values */
void  CurrMixtureSpeciesT::ComputeOutput (const iArrayT& n_codes, dArray2DT& n_values,
	const iArrayT& e_codes, dArray2DT& e_values)
{
	/* number of output values */
	int n_out = n_codes.Sum();
	int e_out = e_codes.Sum();

	/* nothing to output */
	if (n_out == 0 && e_out == 0) return;

//TEMP
#pragma unused(e_values)
if (e_out > 0)
	ExceptionT::GeneralFail("CurrMixtureSpeciesT::ComputeOutput", "element output not supported");

	/* dimensions */
	int nen = NumElementNodes();
	int nsd = NumSD();

	/* reset averaging workspace */
	ElementSupport().ResetAverage(n_out);

	/* work arrays */
	dArray2DT nodal_space(nen, n_out);
	dArray2DT nodal_all(nen, n_out);
	dArray2DT coords, disp;
	dArray2DT flux_elem, matdat;

	/* ip values */
	dSymMatrixT cauchy(nsd);
	dArrayT flux_ip;

	/* set shallow copies */
	double* pall = nodal_space.Pointer();
	coords.Set(nen, n_codes[iNodalCoord], pall);
	pall += coords.Length();
	disp.Set(nen, n_codes[iNodalDisp], pall);
	pall += disp.Length();
	matdat.Set(nen, n_codes[iMaterialData], pall);

	Top();
	while (NextElement())
	{
		/* initialize */
	    nodal_space = 0.0;

		/* global shape function values */
		SetGlobalShape();
		SetLocalU(fLocDisp);
		
		/* coordinates and displacements all at once */
		if (n_codes[iNodalCoord]) fLocInitCoords.ReturnTranspose(coords);
		if (n_codes[iNodalDisp])  fLocDisp.ReturnTranspose(disp);

        /* material stuff */
        if (n_codes[iMaterialData])
        {
            flux_elem.Alias(NumIP(), NumSD(), fMassFlux(CurrElementNumber()));
            /* integrate */
            fShapes->TopIP();
            while (fShapes->NextIP())
            {
                flux_elem.RowAlias(fShapes->CurrIP(), flux_ip);
				fShapes->Extrapolate(flux_ip, matdat);
			}
		}

		/* copy in the cols (in sequence of output) */
		int colcount = 0;
		nodal_all.BlockColumnCopyAt(disp  , colcount); colcount += disp.MinorDim();
		nodal_all.BlockColumnCopyAt(coords, colcount); colcount += coords.MinorDim();
		nodal_all.BlockColumnCopyAt(matdat, colcount);

		/* accumulate - extrapolation done from ip's to corners => X nodes */
		ElementSupport().AssembleAverage(CurrentElement().NodesX(), nodal_all);

/*		const iArrayT nodes = CurrentElement().NodesX();
		for (int i = 0; i < nodes.Length(); i++)
			if (nodes[i] == 221) cout << "\nCompute Output 222: "<<nodal_all(i, 0);
*/
	}
	
	/* get nodally averaged values */
	ElementSupport().OutputUsedAverage(n_values);
/*
	cout << "\nCompute Output: n_values major dim: "<<n_values.MajorDim();
	cout << "\nCompute Output: n_values minor dim: "<<n_values.MinorDim();
	cout << "\nCompute Output: nvalues: ";
	for (int i =0; i<n_values.MinorDim(); i++)
		cout <<n_values(221,i)<<"\t";
*/
}

/* construct output labels array */
void CurrMixtureSpeciesT::GenerateOutputLabels(const iArrayT& n_codes,
	ArrayT<StringT>& n_labels, const iArrayT& e_codes, 
	ArrayT<StringT>& e_labels) const
{
//TEMP - no element labels for now
#pragma unused(e_labels)

	/* allocate node labels */
	n_labels.Dimension(n_codes.Sum());
	int count = 0;	 

	if (n_codes[iNodalDisp])
	{
		/* labels from the field */
		const ArrayT<StringT>& labels = Field().Labels();
		for (int i = 0; i < labels.Length(); i++)
			n_labels[count++] = labels[i];
	}

	if (n_codes[iNodalCoord])
	{
		const char* xlabels[] = {"x1", "x2", "x3"};
		for (int i = 0; i < NumSD(); i++)
			n_labels[count++] = xlabels[i];
	}

	/* material output labels */
	if (n_codes[iMaterialData])
	{
		const char* mlabels[] = {"m1", "m2", "m3"};

		const ArrayT<StringT>& field_labels = Field().Labels();
		for (int i = 0; i < field_labels.Length(); i++) {
            for (int j = 0; j < NumSD(); j++) {
                StringT label = mlabels[j];
                label.Append("_",field_labels[i]);
                n_labels[count++] = label;
            }
        }
	}
	
	if (e_codes.Sum() != 0)
		ExceptionT::GeneralFail("DiffusionElementT::GenerateOutputLabels", 
			"not expecting any element output codes");
}

/* write element output */
void CurrMixtureSpeciesT::WriteOutput(void)
{
	/* inherited */
	NLDiffusionElementT::WriteOutput();

/*	cout << "\nField Name: "<<Field().FieldName();
	const iArray2DT& eq = Field().Equations();
	cout << "\nNumber of equations: "<<eq.Length();
	cout << "\nNumber of nodes: "<<eq.MajorDim();
	cout << "\nNumber of dof: "<<eq.MinorDim();

	int cnt = 0;
	for (int i = 0; i<eq.MajorDim(); i++)
	{
		if (eq(i,0) < 0) {
			cout << "\n "<<i+1 << " "<<eq(i,0);
			cnt ++;
		}
	}
	cout << "\nnumber kcbs: "<<cnt;
*/	
	if (fOutputMass)
	{
		/* compute total species mass */
		double mass = 0.0;
		dArrayT ip_conc(1);
		Top();
		while (NextElement()) 
		{
			/* global shape function values */
			SetGlobalShape();
	
			/* collect nodal concentrations */
			SetLocalU(fLocDisp);
		
			/* loop over integration points */
			const double* j = fShapes->IPDets();
			const double* w = fShapes->IPWeights();
			fShapes->TopIP();
			while (fShapes->NextIP())
			{
				/* ip values */
				fShapes->InterpolateU(fLocDisp, ip_conc);		

				/* accumulate */
				mass += (*j++)*(*w++)*ip_conc[0];
			}
		}
	
		/* output */
		ofstreamT& out = ElementSupport().Output();
		int d_width = OutputWidth(out, &mass);
		out << '\n'
		    << setw(d_width) << ElementSupport().Time()
		    << setw(d_width) << mass
		    << ": time, mass of \"" << Field().FieldName() << "\"" << '\n';
	}
}

/* describe the parameters needed by the interface */
void CurrMixtureSpeciesT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	NLDiffusionElementT::DefineParameters(list);

	/* associated solid element group */
	list.AddParameter(ParameterT::Integer, "solid_element_group");

	/* species type */
	ParameterT species_opt(ParameterT::Enumeration, "species_type");
	species_opt.AddEnumeration("solid", kSolid);
	species_opt.AddEnumeration("fluid", kFluid);
	species_opt.AddEnumeration("solute", kSolute);
	species_opt.SetDefault(kSolid);
	list.AddParameter(species_opt);

	/* velocity of species is calculated wrt this reference frame */
	ParameterT frame(ParameterT::Word, "background_species");
	list.AddParameter(frame, ParameterListT::ZeroOrOnce);

	/* gradient option */
	ParameterT grad_opt(ParameterT::Enumeration, "stress_gradient_option");
	grad_opt.AddEnumeration("global_projection", kGlobalProjection);
	grad_opt.AddEnumeration("element_projection", kElementProjection);
	grad_opt.SetDefault(fGradientOption);
	list.AddParameter(grad_opt);

	/* output total species mass */
	ParameterT output_mass(fOutputMass, "output_mass");
	output_mass.SetDefault(fOutputMass);
	list.AddParameter(output_mass);
}

/* accept parameter list */
void CurrMixtureSpeciesT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "CurrMixtureSpeciesT::TakeParameterList";
 
	/* inherited */
	NLDiffusionElementT::TakeParameterList(list);

	/* resolve background solid element group */
	int solid_element_group = list.GetParameter("solid_element_group");
	solid_element_group--;
	ElementBaseT& element = ElementSupport().ElementGroup(solid_element_group);	
	fUpdatedLagMixture = TB_DYNAMIC_CAST(UpdatedLagMixtureT*, &element);
	fQ1P0Mixture       = TB_DYNAMIC_CAST(Q1P0MixtureT*, &element);
	if (!fUpdatedLagMixture && !fQ1P0Mixture)
		ExceptionT::GeneralFail(caller, "group %d \"%s\" is not a mixture", 
			solid_element_group+1, element.Name().Pointer());
	
	/* checks */
	if (fUpdatedLagMixture) {
		if (fUpdatedLagMixture->NumElements() != NumElements() ||
			fUpdatedLagMixture->NumElementNodes() != NumElementNodes() ||
			fUpdatedLagMixture->NumIP() != NumIP())
			ExceptionT::SizeMismatch(caller);
	} else {
		if (fQ1P0Mixture->NumElements() != NumElements() ||
			fQ1P0Mixture->NumElementNodes() != NumElementNodes() ||
			fQ1P0Mixture->NumIP() != NumIP())
			ExceptionT::SizeMismatch(caller);	
	}

	/* species type */
	int species_opt = list.GetParameter("species_type");
	if (species_opt == kSolid)
		fSpecies = kSolid;
	else if (species_opt == kFluid)
		fSpecies = kFluid;
	else if (species_opt == kSolute)
		fSpecies = kSolute;
	else
		ExceptionT::GeneralFail(caller, "unrecognized \"species type\" %d", species_opt);
	
	/* resolve background species */
	const ParameterT* bg_species = list.Parameter("background_species");
	if (bg_species)
	{
		StringT bg_species_name = *bg_species;
		if (bg_species_name == Field().FieldName())
			ExceptionT::GeneralFail(caller, "background_species must differ from this species \"%s\"",
				Field().FieldName().Pointer());
			
		int num_groups = ElementSupport().NumElementGroups();
		for (int i = 0; !fBackgroundSpecies && i < num_groups; i++) {
			ElementBaseT& element = ElementSupport().ElementGroup(i);
			if (element.Field().FieldName() == bg_species_name) {
				fBackgroundSpecies = TB_DYNAMIC_CAST(CurrMixtureSpeciesT*, &element);
			}
		}
		if (!fBackgroundSpecies)
			ExceptionT::GeneralFail(caller, "could not resolve background_species \"%s\"",
				bg_species_name.Pointer());
	}
	
	
	/* method used to compute stress gradient */
	int grad_opt = list.GetParameter("stress_gradient_option");
	if (grad_opt == kGlobalProjection)
		fGradientOption = kGlobalProjection;
	else if (grad_opt == kElementProjection)
		fGradientOption = kElementProjection;
	else
		ExceptionT::GeneralFail(caller, "unrecognized \"stress_gradient_option\" %d", grad_opt);

	/* output the total system mass */
	fOutputMass = list.GetParameter("output_mass");

    /* dimensions */
    int nip = NumIP();
    int nsd = NumSD();

    /* allocate/initialize */
    fcauchy_ip.Dimension(nip);
    fdcauchy_ip.Dimension(nip);
    for (int i = 0; i < nip; i++) {
        fcauchy_ip[i].Dimension(nsd);
        fdcauchy_ip[i].Dimension(nsd);
    }
    
	/* allocate work space for element-by-element stress projection */
	if (fGradientOption == kElementProjection) 
	{
		/* parent domain information */
		const ParentDomainT& parent_domain = fShapes->ParentDomain();
	
		/* allocate/initialize */
		fip_gradient.Dimension(nip);
		for (int i = 0; i < nip; i++) {
			fip_gradient[i].Dimension(nsd,nip);
			parent_domain.IPGradientTransform(i, fip_gradient[i]);
		}
	}

	/* resolve species index */
	fIndex = (fUpdatedLagMixture) ? 
		fUpdatedLagMixture->SpeciesIndex(Field().FieldName()) :
		fQ1P0Mixture->SpeciesIndex(Field().FieldName());
	if (fIndex < 0)
		ExceptionT::GeneralFail(caller, "could not resolve index of species \"%s\" in \"%s\"",
			Field().FieldName().Pointer(),
			((fUpdatedLagMixture) ? 
				fUpdatedLagMixture->Name().Pointer() : 
				fQ1P0Mixture->Name().Pointer()));

#if 0
	/* consistency check */
	double solid_density = fUpdatedLagMixture->Density(fIndex);
	for (int i = 0; i < fMaterialList->Length(); i++) {
		const ContinuumMaterialT* pcont_mat = (*fMaterialList)[i];
		const DiffusionMaterialT* pdiff_mat = TB_DYNAMIC_CAST(const DiffusionMaterialT*, pcont_mat);
		if (!pdiff_mat) ExceptionT::GeneralFail(caller, "error resolving density of material %d", i+1);
		if (fabs(pdiff_mat->Density() - solid_density) > kSmall)
			ExceptionT::GeneralFail(caller, 
				"density %g of the species at index %d of the solid mixture differs from %g",
					solid_density, fIndex+1, pdiff_mat->Density());
	}
#endif

	/* set concentration type */
	if (fUpdatedLagMixture)
		fUpdatedLagMixture->SetConcentration(fIndex, UpdatedLagMixtureT::kCurrent);
	else
		fQ1P0Mixture->SetConcentration(fIndex, Q1P0MixtureT::kCurrent);	

	/* dimension */
	fFluxVelocity.Dimension(NumElements(), NumIP()*NumSD());
	fMassFlux.Dimension(NumElements(), NumIP()*NumSD());
	fDivBackgroundVel.Dimension(NumElements(), NumIP());
	fNEEmat.Dimension(NumElementNodes());
	fNSDmat1.Dimension(NumSD());
	fNSDmat2.Dimension(NumSD());
	fNSDmat3.Dimension(NumSD());
	
	/* initialize */
	fFluxVelocity = 0.0;
    fMassFlux = 0.0;
	fDivBackgroundVel = 0.0;
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* allocate and initialize shape function objects */
void CurrMixtureSpeciesT::SetShape(void)
{
	/* clear existing */
	delete fShapes;

	/*set current configuration */
	const LocalArrayT& coords = fLocCurrCoords;

	/* construct shape functions in current coordinates*/
	fShapes = new ShapeFunctionT(GeometryCode(), NumIP(), coords);
	fShapes->Initialize();
}

/* allocate and initialize local arrays */
void CurrMixtureSpeciesT::SetLocalArrays(void)
{
	/* inherited */
	NLDiffusionElementT::SetLocalArrays();
	
	/* set up array of current nodal coordinates */
	fLocCurrCoords.Dimension(NumElementNodes(), NumSD());
	ElementSupport().RegisterCoordinates(fLocCurrCoords);
}

/* compute shape functions and derivatives */
void CurrMixtureSpeciesT::SetGlobalShape(void)
{
	/* collect current coordinates */
	SetLocalX(fLocCurrCoords);
	
	/* inherited */
	NLDiffusionElementT::SetGlobalShape();
	
	/* will need deformation gradient */
	if (fUpdatedLagMixture)	
		fUpdatedLagMixture->SetGlobalShape();
	else
		fQ1P0Mixture->SetGlobalShape();
}

/* reset loop */
void CurrMixtureSpeciesT::Top(void)
{
	/* inherited */
	NLDiffusionElementT::Top();

	/* synchronize solid element group */
	if (fUpdatedLagMixture)
		fUpdatedLagMixture->Top();
	else
		fQ1P0Mixture->Top();
}
	
/* advance to next element */ 
bool CurrMixtureSpeciesT::NextElement(void)
{
	/* inherited */
	bool next = NLDiffusionElementT::NextElement();

	/* synchronize solid element group */
	if (fUpdatedLagMixture)
		return fUpdatedLagMixture->NextElement() && next;
	else
		return fQ1P0Mixture->NextElement() && next;
}

/* form the residual force vector */
void CurrMixtureSpeciesT::RHSDriver(void)
{

   /* compute the flux velocities */
	ComputeMassFlux(false);
    
	/* inherited */
	NLDiffusionElementT::RHSDriver();
    
}

/* form group contribution to the stiffness matrix */
void CurrMixtureSpeciesT::LHSDriver(GlobalT::SystemTypeT sys_type)
{
#pragma unused(sys_type)

	/* compute the variation in flux velocities */
	ComputeMassFlux(true);

	/* inherited */
	NLDiffusionElementT::LHSDriver(sys_type);
}

/* calculate the internal force contribution ("-k*d") */
void CurrMixtureSpeciesT::FormKd(double constK)
{
	/* dimensions */
	int nsd = NumSD();
	int nip = NumIP();

	int e = CurrElementNumber();
	/* integration parameters fShapes is in current coordinates*/
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

    dArrayT Na;
	/* mass flux */
	dArray2DT m_e(nip, nsd, fMassFlux(e));
	dArrayT m, div_v;
	
	dMatrixT grad_conc;  /*gradient of concetration at integration point in current coords*/
	dArrayT conc;  /*concetration at integration point*/
	fShapes->TopIP();
	while (fShapes->NextIP())
	{
		int ip = fShapes->CurrIP();
	
		/* retrieve the mass flux */
		m_e.RowAlias(ip, m);

		/* set field gradient */
		grad_conc.Alias(1, nsd, fGradient_list[ip].Pointer());
		fShapes->GradU(fLocDisp, grad_conc);
		
		/* interpolate field */
		conc.Alias(1, fField_list.Pointer(ip));
		fShapes->InterpolateU(fLocDisp, conc);
        
        /*shapefunction*/
//		Na.Set(nen, (double*) fShapes->IPShapeU());

		/* get strain-displacement matrix */
		B(ip, fB);

		/* (div) flux contribution */
		fB.MultTx(m, fNEEvec);
		
		/* c div(v) contribution */
		/* deformation gradients */
		const dMatrixT& F = (fUpdatedLagMixture) ? 
		fUpdatedLagMixture->DeformationGradient(ip) : 
		fQ1P0Mixture->DeformationGradient(ip);
		const dMatrixT& F_last = (fUpdatedLagMixture) ?
		fUpdatedLagMixture->DeformationGradient_last(ip) :
		fQ1P0Mixture->DeformationGradient_last(ip);
			
		/* compute h and h^T h */
		fNSDmat1.DiffOf(F, F_last);           /* GRAD U (8.1.9) */
		fNSDmat2.Inverse(F);                  /* F^-1 */
		fNSDmat3.MultAB(fNSDmat1, fNSDmat2);  /* h (8.1.7) */
		fNSDmat1.MultATB(fNSDmat3, fNSDmat3); /* h^T h */

		/* compute velocity gradient (Simo: 8.1.22) and (Simo: 8.3.13) */
		double dt = ElementSupport().TimeStep();
		double by_dt = (dt > kSmall) ? 1.0/dt : 0.0;
		fNSDmat2.SetToCombination(by_dt, fNSDmat3, -0.5*by_dt, fNSDmat1);

		/* c div(v) */
		double c_div_v = conc[0]*fNSDmat2.Trace();
		fNEEvec.AddScaled(-c_div_v, fShapes->IPShapeU());

		if (fSpecies == kSolute) 
		{
			/*c div(v_f)*/
			dArray2DT v_e_bg;
			dArrayT v_bg(nsd);

			/*grad_x c dot v_f*/
			if (fBackgroundSpecies)
			{
				/* background flux velocity */
				const dArray2DT& bg_flux_velocity = fBackgroundSpecies->FluxVelocity();
				v_e_bg.Alias(nip, nsd, bg_flux_velocity(e));
			
                v_e_bg.RowAlias(ip,v_bg);
                dArrayT vec(1);
                grad_conc.Multx(v_bg,vec);
                fNEEvec.AddScaled(-vec[0], fShapes->IPShapeU());

                fNEEvec.AddScaled(-conc[0]*fDivBackgroundVel(e,ip), fShapes->IPShapeU());
            }
		}
		/* accumulate */
		fRHS.AddScaled(-constK*(*Weight++)*(*Det++), fNEEvec);
	}
}

/* form the element stiffness matrix */
void CurrMixtureSpeciesT::FormStiffness(double constK)
{
	/* must be nonsymmetric */
	if (fLHS.Format() != ElementMatrixT::kNonSymmetric)
		ExceptionT::GeneralFail("CurrMixtureSpeciesT::FormStiffness",
			"LHS matrix must be nonsymmetric");

	/* integration parameters */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	/* dimensions */
	int nsd = NumSD();
	int nip = NumIP();
	int nen = NumElementNodes();	

	/* integrate element stiffness */
	dMatrixT grad_conc;
	dArrayT conc;
	
	/*work space*/
	dArrayT Na;
  	/* (linearization) mass flux */
	dArray2DT dm_e(nip, nsd*nen, fDMassFlux(CurrElementNumber()));
	dMatrixT dm;
  	
 	fShapes->TopIP();
	while (fShapes->NextIP())
	{
		int ip = fShapes->CurrIP();

		double scale = -constK*(*Det++)*(*Weight++);

		/* set field gradient */
		grad_conc.Alias(1, nsd, fGradient_list[ip].Pointer());
		fShapes->GradU(fLocDisp, grad_conc);
		
		/* interpolate field */
		conc.Alias(1, fField_list.Pointer(ip));
		fShapes->InterpolateU(fLocDisp, conc);
        
		/* strain displacement matrix */
		B(ip, fB);

		/* shape function array */
//		Na.Alias(nen, fShapes->IPShapeU());			
		/* shape function array */
		Na.Set(nen, (double*) fShapes->IPShapeU());

		/*div(v)*/
		/* deformation gradients */
		const dMatrixT& F = (fUpdatedLagMixture) ? 
			fUpdatedLagMixture->DeformationGradient(ip) : 
			fQ1P0Mixture->DeformationGradient(ip);
		const dMatrixT& F_last = (fUpdatedLagMixture) ?
			fUpdatedLagMixture->DeformationGradient_last(ip) :
			fQ1P0Mixture->DeformationGradient_last(ip);
		
		/* compute h and h^T h */
		fNSDmat1.DiffOf(F, F_last);           /* GRAD U (8.1.9) */
		fNSDmat2.Inverse(F);                  /* F^-1 */
		fNSDmat3.MultAB(fNSDmat1, fNSDmat2);  /* h (8.1.7) */
		fNSDmat1.MultATB(fNSDmat3, fNSDmat3); /* h^T h */
		
		/* compute velocity gradient (Simo: 8.1.22) and (Simo: 8.3.13) */
		double dt = ElementSupport().TimeStep();
		double by_dt = (dt > kSmall) ? 1.0/dt : 0.0;
		fNSDmat2.SetToCombination(by_dt, fNSDmat3, -0.5*by_dt, fNSDmat1);
		
		double div_v = fNSDmat2.Trace();
		fLHS.Outer(Na, Na, scale*div_v, dMatrixT::kAccumulate);

		if(fSpecies == kFluid)  
		{
            dm.Alias(nsd, nen, dm_e(ip));
            
            fNEEmat.MultATB(fB,dm);
            fLHS.AddScaled(scale, fNEEmat);
        }
		else if (fSpecies == kSolute)
		{            
			/*c div(v_f)*/
			dArray2DT v_e_bg;
			dArrayT v_bg(nsd);
			/*grad_x c dot v_f*/
			if (fBackgroundSpecies)
			{
				/* background flux velocity */
				const dArray2DT& bg_flux_velocity = fBackgroundSpecies->FluxVelocity();
				v_e_bg.Alias(nip, nsd, bg_flux_velocity(CurrElementNumber()));
			
                v_e_bg.RowAlias(ip,v_bg);
                fB.MultTx(v_bg,fNEEvec);
                fLHS.Outer(Na, fNEEvec, -scale, dMatrixT::kAccumulate);

                fLHS.Outer(Na, Na, -scale*fDivBackgroundVel(CurrElementNumber(),ip), dMatrixT::kAccumulate);
            }
            
			/*Diffusion*/
			fD.SetToScaled(-scale, fCurrMaterial->k_ij());
			fLHS.MultQTBQ(fB, fD, dMatrixT::kWhole, dMatrixT::kAccumulate);
            
            /*get derivative of mass flux wrt to c!onc (dD*grad_c) */
            dm.Alias(nsd,nen,dm_e(ip));
            fNEEmat.MultATB(fB,dm);
            fLHS.AddScaled(scale, fNEEmat);
	    }
	}
}

/* compute the relative mass flux and velocities. Implemented now only for fluid species  for momentum driving forces only and for solute species for diffusion only */
void CurrMixtureSpeciesT::ComputeMassFlux(bool compute_dmass_flux)
{
	const char caller[] = "CurrMixtureSpeciesT::ComputeMassFlux";
	
	/* work space */
	int nen = NumElementNodes();
	int nsd = NumSD();
	int nip = NumIP();

	dArrayT ip_conc(1);
	dMatrixT grad_conc(1,nsd); /*current coords*/

    /*workspaces*/
	dArray2DT v_e, m_e,dm_e;
    dMatrixT dm;
	dArrayT v, m, Na;
    dArrayT vec(nsd);

    /*workspaces for calculating divergence*/
    dMatrixT ip_grad_x, ip_grad_cauchy;                   /*ip gradient operator*/
    ip_grad_x.Dimension(nsd,nip);
    dMatrixT jacobian(nsd);

	if(fSpecies == kFluid) 
	{
		/*work spaces*/
        /*workspaces for calculating divergence*/
		LocalArrayT cauchy(LocalArrayT::kUnspecified), dcauchy(LocalArrayT::kUnspecified); 
    
		dMatrixT ip_Grad_val, ip_Grad_val_j;
        dArrayT div_val(nsd);	
        dMatrixT d_divcauchy(nsd,nen), mat(nsd,nen);
		dMatrixT dcauchy_ip(nsd);
		dArray2DT dcauchy_avg;

		dArrayT body_force(nsd);
		LocalArrayT acc(LocalArrayT::kAcc, NumElementNodes(), nsd);
		dArrayT ip_acc(nsd);		
		dArrayT force(nsd); /*v=D*force, force = -phi*/;

		if (fGradientOption == kGlobalProjection)
		{	
			/* get concentration specific nodal Kirchhoff stresses (tau) */
			if (fUpdatedLagMixture) fUpdatedLagMixture->ProjectPartialCauchy(fIndex);
			else fQ1P0Mixture->ProjectPartialCauchy(fIndex);

			fcauchy_avg = ElementSupport().OutputAverage();
			cauchy.Dimension(NumElementNodes(), nsd*nsd);		
			cauchy.SetGlobal(fcauchy_avg);
			ip_Grad_val.Dimension(nsd*nsd, nsd);

			/* project variation in partial stresses to the nodes */
			if (compute_dmass_flux) {
				if (fUpdatedLagMixture) fUpdatedLagMixture->ProjectDPartialCauchy(fIndex);
				else fQ1P0Mixture->ProjectDPartialCauchy(fIndex);
				dcauchy_avg.Alias(ElementSupport().OutputAverage());
				dcauchy.Dimension(NumElementNodes(), nsd*nsd);		
				dcauchy.SetGlobal(dcauchy_avg);
			}

		}
        else {
            ip_grad_x.Dimension(nsd, nip);
        }

		/*calculate momentum driving force*/
		/* get the body force */
		if (fUpdatedLagMixture) fUpdatedLagMixture->BodyForce(body_force);
		else fQ1P0Mixture->BodyForce(body_force);
        
        /* initialize mass flux variation */
        if (compute_dmass_flux)
            fDMassFlux.Dimension(NumElements(), nip*nsd*nen);
        fDMassFlux = 0.0;

		Top();
		while (NextElement()) 
		{
			int e = CurrElementNumber();
			/* global shape function values */
			SetGlobalShape();
	
			/* collect nodal stresses */
			if (fGradientOption == kGlobalProjection) 
                SetLocalU(cauchy);

            /* collect integration point stresses - sets shapes functions over the element */
/*            if (fUpdatedLagMixture)
                fUpdatedLagMixture->IP_PartialStress(fIndex, &fcauchy_ip, (compute_dmass_flux) ? &fdcauchy_ip : NULL);
            else
				fQ1P0Mixture->IP_PartialStress(fIndex, &fcauchy_ip, (compute_dmass_flux) ? &fdcauchy_ip : NULL);
*/
            if (fUpdatedLagMixture)
                fUpdatedLagMixture->IP_PartialCauchy(fIndex, &fcauchy_ip, (compute_dmass_flux) ? &fdcauchy_ip : NULL);
            else
				fQ1P0Mixture->IP_PartialCauchy(fIndex, &fcauchy_ip, (compute_dmass_flux) ? &fdcauchy_ip : NULL);
				
            /* collect nodal accelerations */
            if (fUpdatedLagMixture)
                fUpdatedLagMixture->Acceleration(acc);
            else
                fQ1P0Mixture->Acceleration(acc);

			/* collect nodal concentrations */
			SetLocalU(fLocDisp);
		
			/* mass flux and velocity */
			v_e.Alias(nip, nsd, fFluxVelocity(e));
			m_e.Alias(nip, nsd, fMassFlux(e));
            if (compute_dmass_flux) dm_e.Alias(nip, nen*nsd, fDMassFlux(e));

			/* loop over integration points */
			fShapes->TopIP();
			while (fShapes->NextIP())
			{
				int ip = fShapes->CurrIP();

				/* ip values of current concentration and acceleration*/
				fShapes->InterpolateU(fLocDisp, ip_conc);
				fShapes->InterpolateU(acc, ip_acc);		

				/* inertial forces */
				force.SetToCombination(-ip_conc[0], ip_acc, ip_conc[0], body_force);
 
                /* compute stress divergence */
				if (fGradientOption == kGlobalProjection)
				{
					div_val = 0.0;

/*					if (fUpdatedLagMixture)
						fShapes->GradU(cauchy, ip_Grad_val);  
					else ExceptionT::GeneralFail(caller, "Not implemented for Q1P0");
*/
					fShapes->GradU(cauchy, ip_Grad_val);   /*gradient wrt to current configuration*/
										
					for (int j = 0; j < nsd; j++) {
						ip_Grad_val_j.Alias(nsd, nsd, ip_Grad_val(j));
						for (int i = 0; i < nsd; i++)
							div_val[i] += ip_Grad_val_j(i,j);
					}
				}
				else /* element-by-element gradient calculation */
				{
					/* transform gradient matrix to element (curr) coordinates */
					fShapes->ParentDomain().DomainJacobian(fLocCurrCoords, ip, jacobian);
					jacobian.Inverse();
					ip_grad_x.MultATB(jacobian, fip_gradient[ip]);

					/* compute divergence of tau */
                    /* dimensions */
                    int nsd = ip_grad_x.Rows();
                    int nip = ip_grad_x.Cols();

                    div_val = 0.0;
                    for (int k = 0; k < nip; k++) {
                        const dMatrixT& ip_stress = fcauchy_ip[k];
                        for (int i = 0 ; i < nsd; i++)
                            for (int j = 0; j < nsd; j++) /* div */
                                div_val[i] += ip_grad_x(j,k)*ip_stress(i,j);			
                    }

				}
				/* add divergence of stress term to force */
				force += div_val;
                    
				/* v = D div_sig

				v_e.RowAlias(ip, v);

				const dMatrixT& D = fCurrMaterial->k_ij();
				D.Multx(force, v); 

				m_e.RowAlias(ip, m);
				m.SetToScaled(ip_conc[0],v);
				*/

				/* m = D div_sig*/
				m_e.RowAlias(ip, m);
				v_e.RowAlias(ip, v);

				const dMatrixT& D = fCurrMaterial->k_ij();
				D.Multx(force, m); /* c*V */

				v.SetToScaled(1.0/ip_conc[0],m);
				
				
                /* mass flux variation */
                if (compute_dmass_flux)
                {
                    dm.Alias(nsd, nen, dm_e(ip));
    
                    /* shape function array */
                    Na.Alias(nen, fShapes->IPShapeU());			

					/* v = D div_sig*/
                    /* contribution from changing diffusivity -cf dB/dc phi*/
                    /*const dMatrixT& dD = fCurrMaterial->dk_ij();
                    dD.Multx(force, vec, ip_conc[0]); 
					*/
                    /*contribution form -B(2 cf dv/dt - 2 cf g - div sigma)*/
                    /*force.AddScaled(-ip_conc[0], ip_acc);
                    force.AddScaled(ip_conc[0], body_force);
                    const dMatrixT& D = fCurrMaterial->k_ij();
                    D.Multx(force, vec, 1.0, dMatrixT::kAccumulate);
                    */
					
					/*m = D div_sig*/
					/*contribution from changing diffusivity, -dB/dC phi*/
                    const dMatrixT& dD = fCurrMaterial->dk_ij();
                    dD.Multx(force, vec); 
					/*contribution form -B( dv/dt -  g )*/
					const dMatrixT& D = fCurrMaterial->k_ij();
					D.Multx(ip_acc, vec, -1.0, dMatrixT::kAccumulate);
					D.Multx(body_force, vec, dMatrixT::kAccumulate);
					
                    dm.Outer(vec, Na, 1.0, dMatrixT::kAccumulate);

                    /* compute divergence of dstress/dconc */
                    if (fGradientOption == kGlobalProjection)
                        ComputeDDivergence(fdcauchy_ip, d_divcauchy);
//                        ComputeDDivergence(dcauchy, d_divcauchy,dcauchy_ip);
                    else /* element-by-element gradient calculation */
                        ComputeDDivergence(ip_grad_x, fdcauchy_ip, d_divcauchy);

                    mat.MultAB(D, d_divcauchy);

					/* v = D div_sig*/
                    //dm.AddScaled(ip_conc[0], mat);

					/* m = D div_sig*/
                    dm.AddScaled(1.0, mat);
                }
            }
        }
    }
	else if (fSpecies == kSolute) 
	{		
		if (fGradientOption == kGlobalProjection){
			ProjectV();
			fv_bg_avg = ElementSupport().OutputAverage();
		}

        /* initialize mass flux variation */
        if (compute_dmass_flux)
            fDMassFlux.Dimension(NumElements(), nip*nsd*nen);
        fDMassFlux = 0.0;

		Top();
		while (NextElement()) 
		{
			int e = CurrElementNumber();

			/* global shape function values */
			SetGlobalShape();

			/* collect nodal concentrations */
			SetLocalU(fLocDisp);
			
			m_e.Alias(nip, nsd, fMassFlux(e));
			v_e.Alias(nip, nsd, fFluxVelocity(e));

			/* loop over integration points */
			fShapes->TopIP();
			while (fShapes->NextIP())
			{
                if (compute_dmass_flux) dm_e.Alias(nip, nen*nsd, fDMassFlux(e));

				int ip = fShapes->CurrIP();

				/* ip values of current concentration*/  //confirm that fLocDisp is current concentration*/
				fShapes->InterpolateU(fLocDisp, ip_conc);
				
				if (fGradientOption == kGlobalProjection){
					LocalArrayT v_bg(LocalArrayT::kUnspecified);
					dMatrixT ip_grad_v(nsd);

					v_bg.Dimension(NumElementNodes(), nsd);		
					v_bg.SetGlobal(fv_bg_avg);

					fShapes->GradU(v_bg, ip_grad_v);
                    fDivBackgroundVel(e,ip) = 0.0;
					for (int i = 0; i < nsd; i++)
						fDivBackgroundVel(e,ip) += ip_grad_v(i,i);
				}
				else{
					dArray2DT v_e_bg;
					dArrayT v_bg(nsd);
					/* retrieve background fluid velocity */
					if (fBackgroundSpecies)
					{
						/* background flux velocity */
						const dArray2DT& bg_flux_velocity = fBackgroundSpecies->FluxVelocity();
						v_e_bg.Alias(nip, nsd, bg_flux_velocity(e));
												 
                        dMatrixT jacobian(nsd);
                        fShapes->ParentDomain().DomainJacobian(fLocCurrCoords, ip, jacobian);
                        jacobian.Inverse();
                        ip_grad_x.MultATB(jacobian, fip_gradient[ip]);

                        fDivBackgroundVel(e,ip) = 0.0;
                        for (int k = 0; k < nip; k++) {
                            v_e_bg.RowAlias(k,v_bg);
                            for (int i = 0 ; i < nsd; i++)
                                fDivBackgroundVel(e,ip) += ip_grad_x(i,k)*v_bg[i];
                        }
                    }
//					else ExceptionT::GeneralFail(caller, "Background velocity not set");	/* work space */
				}

				/*ip values spatial gradient of current concentration*/
				fShapes->GradU(fLocDisp, grad_conc);   /*gradient wrt current coordinates*/
				
				/*calculate diffusive flux and velocities*/
				const dMatrixT& D = fCurrMaterial->k_ij();
				m_e.RowAlias(ip,m);
				D.Multx(grad_conc,m, -1.0);
				
				v_e.RowAlias(ip,v);
				v.SetToScaled(1.0/ip_conc[0],m);
                
                /* mass flux variation with respect to conc: only dD grad_c implemented for now */
                /* mass flux variation */
                if (compute_dmass_flux)
                {
                    dm.Alias(nsd, nen, dm_e(ip));
    
                    /* shape function array */
                    Na.Alias(nen, fShapes->IPShapeU());			

                    /* contribution from changing diffusivity -cf dB/dc phi*/
                    const dMatrixT& dD = fCurrMaterial->dk_ij();

                    dD.Multx(grad_conc, vec);
                    dm.Outer(vec, Na, 1.0, dMatrixT::kAccumulate);
                }
			}
		}
	}
}


/* project the given partial first Piola-Kirchoff stress to the nodes */
void CurrMixtureSpeciesT::ProjectV(void)
{
	const char caller[] = "CurrMixtureSpeciesT::ProjectV";

	/* dimensions */
	int nen = NumElementNodes();
	int nsd = NumSD();
	int nip = NumIP();
	
	dArray2DT nodal_v_bg(nen, nsd);
	dArray2DT v_e_bg;
	dArrayT v_bg;
	
	
	/* reset averaging workspace */
	ElementSupport().ResetAverage(nsd);
	
	Top();
	while (NextElement()) 
	{
		int e = CurrElementNumber();
		if (CurrentElement().Flag() != ElementCardT::kOFF)
		{
			/* global shape function values */
			SetGlobalShape();
			
			/* retrieve background fluid velocity */
			if (fBackgroundSpecies)
			{
				/* background flux velocity */
				const dArray2DT& bg_flux_velocity = fBackgroundSpecies->FluxVelocity();
				v_e_bg.Alias(nip, nsd, bg_flux_velocity(e));
										 
                /* extrapolate element background velocities */
                nodal_v_bg = 0.0;
                fShapes->TopIP();
                while (fShapes->NextIP())
                {
                    int ip = fShapes->CurrIP();
                    v_e_bg.RowAlias(ip, v_bg);
				
                    /* extrapolate to the nodes */
                    fShapes->Extrapolate(v_bg, nodal_v_bg);
                }
			}
//			else ExceptionT::GeneralFail(caller, "Background velocity not set");	/* work space */
			else nodal_v_bg = 0.0;	/* work space */
            
			/* accumulate - extrapolation done from ip's to corners => X nodes */
			ElementSupport().AssembleAverage(CurrentElement().NodesX(), nodal_v_bg);
		}
	}
}

/* compute variation of divergence with respect to the nodal values */
void CurrMixtureSpeciesT::ComputeDDivergence(const dMatrixT& ip_grad_transform, 
	const ArrayT<dMatrixT>& tensor_ip, dMatrixT& d_div) const
{
	/* dimensions */
	int nen = NumElementNodes();
	int nsd = ip_grad_transform.Rows();
	int nip = ip_grad_transform.Cols();	
	
	dArrayT Na;
	d_div = 0.0;
	for (int k = 0; k < nip; k++) 
	{
		/* shape function array */
		Na.Alias(nen, fShapes->IPShapeU(k));			
	
		const dMatrixT& A_k = tensor_ip[k];
		for (int i = 0 ; i < nsd; i++)
			for (int j = 0; j < nsd; j++) /* div */
				for (int n = 0; n < nen; n++)
					d_div(i,n) += ip_grad_transform(j,k)*A_k(i,j)*Na[n];
	}
}

void CurrMixtureSpeciesT::ComputeDDivergence(const ArrayT<dMatrixT>& tensor_ip, dMatrixT& d_div) const
{
	/* dimensions */
	int nen = NumElementNodes();
	int nip = fShapes->NumIP();
	int nsd = NumSD();
	
	/* shape function derivatives (at the current integration point) */
	const dArray2DT& DNa = fShapes->Derivatives_U();

	dArrayT Na;
	d_div = 0.0;
	for (int k = 0; k < nip; k++) 
	{
        /* extrapolation matrix */
        const dMatrixT& extrap = fShapes->Extrapolation();

		/* shape function array */
		Na.Alias(nen, fShapes->IPShapeU(k));			

		const dMatrixT& dval_ip = tensor_ip[k];

		for (int i = 0 ; i < nsd; i++)
			for (int j = 0; j < nsd; j++) /* div */
				for (int n = 0; n < nen; n++)
					for (int m = 0; m < nen; m++)
 						d_div(i,n) += DNa(j,m)*extrap(m,k)*dval_ip(i,j)*Na[n];
 	}
}

/*
void CurrMixtureSpeciesT::ComputeDDivergence(const LocalArrayT& nodal_dval, dMatrixT& d_div, dMatrixT& dval_ip) const
{
	int nen = NumElementNodes();
	int nip = fShapes->NumIP();
	int nsd = NumSD();
	
	const dMatrixT& extrap = fShapes->Extrapolation();
	const dArray2DT& DNa = fShapes->Derivatives_U();

	dArrayT Na;
	d_div = 0.0;
	for (int k = 0; k < nip; k++) 
	{
		Na.Alias(nen, fShapes->IPShapeU(k));			

		fShapes->InterpolateU(nodal_dval, dval_ip, k);

		for (int i = 0 ; i < nsd; i++)
			for (int j = 0; j < nsd; j++) 
				for (int n = 0; n < nen; n++)
					for (int m = 0; m < nen; m++)
 						d_div(i,n) += DNa(j,m)*extrap(m,k)*dval_ip(i,j)*Na[n];
 	}
}
*/
