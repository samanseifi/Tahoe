/* $Id: Q1P0MixtureT.cpp,v 1.2 2006/04/14 15:28:32 thao Exp $ */
#include "Q1P0MixtureT.h"
#include "ShapeFunctionT.h"
#include "FSSolidMixtureT.h"
#include "SolidMatListT.h"
#include "ScheduleT.h"
#include "eIntegratorT.h"

using namespace Tahoe;

/* constructor */
Q1P0MixtureT::Q1P0MixtureT(const ElementSupportT& support):
	SimoQ1P0(support)
{
	SetName("Q1P0_mixture");
}

/* resolve the species name into the index */
int Q1P0MixtureT::SpeciesIndex(const StringT& field_name) const {
	return FSSolidMixture().SpeciesIndex(field_name);
}

/* density of the given species */
double Q1P0MixtureT::Density(int i) {
	return FSSolidMixture().Density(i);
}

/* set concentration flag */
void Q1P0MixtureT::SetConcentration(int i, ConcentrationT conc)
{
	const char caller[] = "Q1P0MixtureT::SetConcentration";

	/* get material */
	FSSolidMixtureT& mixture = FSSolidMixture();

	/* set flag */
	if (conc == kReference)
		mixture.SetConcentration(i, FSSolidMixtureT::kReference);
	else if (conc == kCurrent)
		mixture.SetConcentration(i, FSSolidMixtureT::kCurrent);	
	else
		ExceptionT::GeneralFail(caller, "unrecognized flag %d", conc);
	
	/* (re-)set form of element stiffness matrix */
	GlobalT::SystemTypeT type = TangentType();
	if (type == GlobalT::kSymmetric)
		fLHS.SetFormat(ElementMatrixT::kSymmetricUpper);
	else if (type == GlobalT::kNonSymmetric)
		fLHS.SetFormat(ElementMatrixT::kNonSymmetric);
	else if (type == GlobalT::kDiagonal)
		fLHS.SetFormat(ElementMatrixT::kDiagonal);
}

/* project the given partial first Piola-Kirchoff stress to the nodes */
void Q1P0MixtureT::ProjectPartialStress(int i)
{
	/* dimensions */
	int nen = NumElementNodes();
	int nsd = NumSD();

	/* reset averaging workspace */
	ElementSupport().ResetAverage(nsd*nsd);

	/* work space */
	dMatrixT P(nsd);
	dArrayT P_1D;
	P_1D.Alias(P);

	/* loop over elements */
	dArray2DT nodal_P(nen, nsd*nsd);
	Top();
	while (NextElement())
		if (CurrentElement().Flag() != ElementCardT::kOFF)
		{
			/* get material */
			FSSolidMixtureT& mixture = FSSolidMixture();
		
			/* global shape function values */
			SetGlobalShape();
			
			/* collect concentration */
			mixture.UpdateConcentrations(i);

			/* extrapolate element stresses */
			nodal_P = 0.0;
			fShapes->TopIP();
			while (fShapes->NextIP())
			{
				/* Cauchy stress */
				const dSymMatrixT& cauchy = mixture.s_ij(i);
				
				/* Cauchy -> 1st PK stress */
				cauchy.ToMatrix(fStress);
				const dMatrixT& F = DeformationGradient();
				fF_inv.Inverse(F);
				P.MultABT(fStress, fF_inv);
				P *= F.Det();

				/* extrapolate to the nodes */
				fShapes->Extrapolate(P_1D, nodal_P);
			}

			/* accumulate - extrapolation done from ip's to corners => X nodes */
			ElementSupport().AssembleAverage(CurrentElement().NodesX(), nodal_P);
		}
}

/* project the given partial first Piola-Kirchoff stress to the nodes */
void Q1P0MixtureT::ProjectPartialCauchy(int i)
{
	/* dimensions */
	int nen = NumElementNodes();
	int nsd = NumSD();
	
	/* reset averaging workspace */
	ElementSupport().ResetAverage(nsd*nsd);
	
	/* work space */
	dMatrixT cauchy(nsd);
	dArrayT cauchy_1D;
	cauchy_1D.Alias(cauchy);

	/* loop over elements */
	dArray2DT nodal_cauchy(nen, nsd*nsd);

	Top();
	while (NextElement())
		if (CurrentElement().Flag() != ElementCardT::kOFF)
		{
			/* get material */
			FSSolidMixtureT& mixture = FSSolidMixture();
			
			/* global shape function values */
			SetGlobalShape();
			
			/* collect concentration */
			mixture.UpdateConcentrations(i);
			
			/* extrapolate element stresses */
			nodal_cauchy = 0.0;
			fCurrShapes->TopIP();
			while (fCurrShapes->NextIP())
			{
				/* Cauchy stress */
				const dSymMatrixT& stress = mixture.s_ij(i);
                stress.ToMatrix(cauchy);
                
                const dArrayT& conc = mixture.Get_IPConcentration();	

				/* extrapolate to the nodes */
				fCurrShapes->Extrapolate(cauchy_1D, nodal_cauchy);
			}
			
			/* accumulate - extrapolation done from ip's to corners => X nodes */
			ElementSupport().AssembleAverage(CurrentElement().NodesX(), nodal_cauchy);
		}
}

/* project the variation with concentration of the given partial first
 * Piola-Kirchoff stress to the nodes */
void Q1P0MixtureT::ProjectDPartialStress(int i)
{
	/* dimensions */
	int nen = NumElementNodes();
	int nsd = NumSD();

	/* reset averaging workspace */
	ElementSupport().ResetAverage(nsd*nsd);

	/* work space */
	dMatrixT P(nsd), F_inv(nsd), s(nsd);
	dArrayT P_1D;
	P_1D.Alias(P);

	/* loop over elements */
	dArray2DT nodal_P(nen, nsd*nsd);
	Top();
	while (NextElement())
		if (CurrentElement().Flag() != ElementCardT::kOFF)
		{
			/* get materials */
			FSSolidMixtureT& mixture = FSSolidMixture();

			/* global shape function values */
			SetGlobalShape();
			
			/* collect concentration */
			mixture.UpdateConcentrations(i);

			/* extrapolate element stresses */
			nodal_P = 0.0;
			fShapes->TopIP();
			while (fShapes->NextIP())
			{
				/* Cauchy stress */
				const dSymMatrixT& dcauchy = mixture.ds_ij_dc(i);
				
				/* Cauchy -> 1st PK stress */
				dcauchy.ToMatrix(fStress);
				const dMatrixT& F = DeformationGradient();
				fF_inv.Inverse(F);
				P.MultABT(fStress, fF_inv);
				P *= F.Det();

				/* extrapolate to the nodes */
				fShapes->Extrapolate(P_1D, nodal_P);
			}

			/* accumulate - extrapolation done from ip's to corners => X nodes */
			ElementSupport().AssembleAverage(CurrentElement().NodesX(), nodal_P);
		}
}

/* project the variation with concentration of the given partial first
 * Piola-Kirchoff stress to the nodes */
void Q1P0MixtureT::ProjectDPartialCauchy(int i)
{
	/* dimensions */
	int nen = NumElementNodes();
	int nsd = NumSD();

	/* reset averaging workspace */
	ElementSupport().ResetAverage(nsd*nsd);

	/* work space */
	dMatrixT cauchy(nsd), s(nsd);
	dArrayT cauchy_1D;
	cauchy_1D.Alias(fStress);

	/* loop over elements */
	dArray2DT nodal_cauchy(nen, nsd*nsd);
	Top();
	while (NextElement())
		if (CurrentElement().Flag() != ElementCardT::kOFF)
		{
			/* get materials */
			FSSolidMixtureT& mixture = FSSolidMixture();

			/* global shape function values */
			SetGlobalShape();
			
			/* collect concentration */
			mixture.UpdateConcentrations(i);

			/* extrapolate element stresses */
			nodal_cauchy = 0.0;
			fShapes->TopIP();
			while (fShapes->NextIP())
			{
				/* Cauchy stress */
				const dSymMatrixT& dcauchy = mixture.ds_ij_dc_exact(i);
				
				/* Cauchy -> 1st PK stress */
				dcauchy.ToMatrix(fStress);
                
				/* extrapolate to the nodes */
				fShapes->Extrapolate(cauchy_1D, nodal_cauchy);
			}

			/* accumulate - extrapolation done from ip's to corners => X nodes */
			ElementSupport().AssembleAverage(CurrentElement().NodesX(), nodal_cauchy);
		}
}

void Q1P0MixtureT::IP_PartialStress(int i, ArrayT<dMatrixT>* ip_stress, 
	ArrayT<dMatrixT>* ip_dstress)
{
	/* nothing wanted */
	if (!ip_stress && !ip_dstress)
		return;
	/* element is active */
	else if (CurrentElement().Flag() != ElementCardT::kOFF)
	{
		/* get materials */
		FSSolidMixtureT& mixture = FSSolidMixture();

		/* collect concentration */
		mixture.UpdateConcentrations(i);

		/* collect integration point element stresses */
		fShapes->TopIP();
		while (fShapes->NextIP())
		{
			/* destination */
			int ip = fShapes->CurrIP();
		
			/* deformation gradient */
			const dMatrixT& F = DeformationGradient();
			fF_inv.Inverse(F);

			/* stress */
			if (ip_stress)
			{
				/* Cauchy stress */
				const dSymMatrixT& cauchy = mixture.s_ij(i);
				cauchy.ToMatrix(fStress);
			
				/* Cauchy -> 1st PK stress */
				dMatrixT& P = (*ip_stress)[ip];
				P.MultABT(fStress, fF_inv);
				P *= F.Det();
			}

			/* stress variation */
			if (ip_dstress)
			{
				/* Cauchy stress */
				const dSymMatrixT& dcauchy = mixture.ds_ij_dc(i);
				dcauchy.ToMatrix(fStress);
			
				/* Cauchy -> 1st PK stress */
				dMatrixT& dP = (*ip_dstress)[ip];
				dP.MultABT(fStress, fF_inv);
				dP *= F.Det();
			}
		}
	}
	else /* zero them out */
	{
		if (ip_stress)
			for (int i = 0; i < ip_stress->Length(); i++)
				(*ip_stress)[i] = 0.0;

		if (ip_dstress)
			for (int i = 0; i < ip_dstress->Length(); i++)
				(*ip_dstress)[i] = 0.0;
	}
}

void Q1P0MixtureT::IP_PartialCauchy(int i, ArrayT<dMatrixT>* ip_stress, ArrayT<dMatrixT>* ip_dstress)
{
	/* nothing wanted */
	if (!ip_stress && !ip_dstress)
		return;
	/* element is active */
	else if (CurrentElement().Flag() != ElementCardT::kOFF)
	{
		/* get materials */
		FSSolidMixtureT& mixture = FSSolidMixture();
		
		/* collect concentration */
		mixture.UpdateConcentrations(i);
		
		/* collect integration point element stresses */
		fShapes->TopIP();
		while (fShapes->NextIP())
		{
			/* destination */
			int ip = fShapes->CurrIP();			

            /* stress */
			if (ip_stress)
			{
				/* Cauchy stress */
				const dSymMatrixT& cauchy = mixture.s_ij(i);
 				dMatrixT& mat_stress = (*ip_stress)[ip];
				cauchy.ToMatrix(mat_stress);
			}

			/* stress variation */
			if (ip_dstress)
			{
				/* Cauchy stress */
				const dSymMatrixT& dcauchy = mixture.ds_ij_dc_exact(i);
				dMatrixT& dstress = (*ip_dstress)[ip];
				dcauchy.ToMatrix(dstress);
			}
		}
	}
	else /* zero them out */
	{
		if (ip_stress)
			for (int i = 0; i < ip_stress->Length(); i++)
				(*ip_stress)[i] = 0.0;
		if (ip_dstress)
			for (int i = 0; i < ip_dstress->Length(); i++)
				(*ip_dstress)[i] = 0.0;
	}
}

/* return the nodal accelerations over the current element */
void Q1P0MixtureT::Acceleration(LocalArrayT& acc) const
{
	if (fIntegrator->Order() == 2) {
		acc.SetGlobal(fLocAcc.Global());
		SetLocalU(acc);
	}
	else
		acc = 0.0;
}

/* return the nodal velocities over the current element */
void Q1P0MixtureT::Velocity(LocalArrayT& vel) const
{
	/* should have been set during SetGlobalShape because FSSolidMixtureT
	 * "needs" velocity information */
	vel = fLocVel;
}

/* return the body force vector */
void Q1P0MixtureT::BodyForce(dArrayT& body_force) const
{
	if (fBodySchedule)
		body_force.SetToScaled(fBodySchedule->Value(), fBody);
	else
		body_force = 0.0;
}

/* accept parameter list */
void Q1P0MixtureT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SimoQ1P0::TakeParameterList(list);

	/* dimension work space */
	int nsd = NumSD();
	fF_inv.Dimension(nsd);
	fStress.Dimension(nsd);
}

/***********************************************************************
 * Protected
 ***********************************************************************/

const FSSolidMixtureT& Q1P0MixtureT::FSSolidMixture(void) const 
{
	const char caller[] = "Q1P0MixtureT::FSSolidMixture";

#if __option(extended_errorcheck)
	if (fMaterialList->Length() > 1)
		ExceptionT::GeneralFail(caller, "expecting only 1 material %d",
			fMaterialList->Length());
#endif

	const ContinuumMaterialT* pcont_mat = (*fMaterialList)[0]; /* just use first material */
	if (!pcont_mat) 
		ExceptionT::GeneralFail(caller, "material 0 is NULL");

	const FSSolidMixtureT* mixture = TB_DYNAMIC_CAST(const FSSolidMixtureT*, pcont_mat);
	if (!mixture)
		ExceptionT::GeneralFail(caller, "material \"%s\" is not a mixture",
			pcont_mat->Name().Pointer());

	return *mixture;
}
