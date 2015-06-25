/* $Id: EAMFCC3DMatT_edge.cpp,v 1.5 2011/12/01 20:38:15 beichuan Exp $ */
/* created: hspark (6/2/2009) */
#include "EAMFCC3DMatT_edge.h"

#include "EAMFCC3DSym_edge.h"
#include "dMatrixT.h"

#include <cmath>

using namespace Tahoe;

/* constructor */
EAMFCC3DMatT_edge::EAMFCC3DMatT_edge(void):
	ParameterInterfaceT("FCC_EAM_Surf"),
	fSurfaceThickness(-1),
	fEAM(NULL)
{

}

/* destructor */
EAMFCC3DMatT_edge::~EAMFCC3DMatT_edge(void) { delete fEAM; }

/* describe the parameters needed by the interface */
void EAMFCC3DMatT_edge::DefineParameters(ParameterListT& list) const
{	
	/* inherited */
	NL_E_MatT::DefineParameters(list);

	/* number of neighbor shells */
	ParameterT n_shells(ParameterT::Integer, "shells");
	n_shells.AddLimit(1, LimitT::LowerInclusive);
	list.AddParameter(n_shells);
	
	/* surface normal */
	ParameterT normal(ParameterT::Integer, "normal_code");
	normal.AddLimit(0, LimitT::LowerInclusive);
	normal.AddLimit(5, LimitT::UpperInclusive);
	list.AddParameter(normal);
}

/* describe the parameters needed by the interface */
void EAMFCC3DMatT_edge::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	NL_E_MatT::DefineSubs(sub_list);

	/* Cauchy-Born EAM parameters */
	sub_list.AddSub("FCC_EAM_Cauchy-Born");
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* EAMFCC3DMatT_edge::NewSub(const StringT& name) const
{
	if (name == "FCC_EAM_Cauchy-Born")
		return new EAMFCC3DSym_edge(0, 0);
	else /* inherited */
		return NL_E_MatT::NewSub(name);
}

/* accept parameter list */
void EAMFCC3DMatT_edge::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	NL_E_MatT::TakeParameterList(list);

	/* construct Cauchy-Born EAM solver */
	int shells = list.GetParameter("shells");
	int normal_code = list.GetParameter("normal_code");
	fEAM = new EAMFCC3DSym_edge(shells, normal_code);

	fEAM->TakeParameterList(list.GetList("FCC_EAM_Cauchy-Born"));

	/* TEMP - GET SURFACE THICKNESS FROM EAMFCC3D_SURF */
	fSurfaceThickness = fEAM->SurfaceThickness();
	
	/* reset density from the atomistic parameters */
	fDensity = fEAM->Density();
	
	/* Test Edge EAM Implementation - call EAM.cpp functions */
	dSymMatrixT strain(3), PK2(3);
	strain = 0.0;	// calculate energy density for 0 strain configuration
//	double value = ComputeEnergyDensity(strain);
	ComputePK2(strain,PK2);
	cout << "Edge PK2 = " << PK2 << endl;
}

/*************************************************************************
 * Private
 *************************************************************************/

void EAMFCC3DMatT_edge::ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli)
{
	/* Temporarily override for finite difference approximation */
	moduli = FSSolidMatT::C_IJKL();

	/* EAM solver */
	//fEAM->Moduli(moduli, E);
}

/* returns the stress corresponding to the given strain - the strain
* is the reduced index Green strain vector, while the stress is the
* reduced index 2nd Piola-Kirchhoff stress vector */
void EAMFCC3DMatT_edge::ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2)
{
	/* EAM solver */
	fEAM->SetStress(E, PK2);
}

/* returns the strain energy density for the specified strain */
double EAMFCC3DMatT_edge::ComputeEnergyDensity(const dSymMatrixT& E)
{
	/* EAM solver */
	return fEAM->EnergyDensity(E);
}
