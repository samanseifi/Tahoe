/* $Id: FCC3D.cpp,v 1.12 2011/12/01 21:11:37 bcyansfn Exp $ */
/* created: paklein (07/01/1996) */
#include "FCC3D.h"

#include "ElementsConfig.h"
#include "FCCLatticeT.h"
#include <cmath>

/* pair properties */
#ifdef PARTICLE_ELEMENT
#include "HarmonicPairT.h"
#include "LennardJonesPairT.h"
#else
#pragma message("FCC3D requires PARTICLE_ELEMENT")
#error "FCC3D requires PARTICLE_ELEMENT"
#endif

using namespace Tahoe;

/* constructor */
FCC3D::FCC3D(void):
	ParameterInterfaceT("FCC_3D"),
	fNearestNeighbor(-1),
	fFCCLattice(NULL),
	fPairProperty(NULL),
	fAtomicVolume(0),
	fBondTensor4(dSymMatrixT::NumValues(3)),
	fBondTensor2(dSymMatrixT::NumValues(3)),
	fFullDensityForStressOutput(true)	
{

}

/* destructor */
FCC3D::~FCC3D(void)
{
	delete fFCCLattice;
	delete fPairProperty;
}

/* describe the parameters needed by the interface */
void FCC3D::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	NL_E_MatT::DefineParameters(list);

	/* number of neighbor shells */
	ParameterT n_shells(ParameterT::Integer, "shells");
	n_shells.AddLimit(1, LimitT::LowerInclusive);
	list.AddParameter(n_shells);
}

/* information about subordinate parameter lists */
void FCC3D::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	NL_E_MatT::DefineSubs(sub_list);

	/* FCC lattice */
	sub_list.AddSub("CB_lattice_FCC");

	/* pair potential choice */
	sub_list.AddSub("FCC_3D_potential_choice", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void FCC3D::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	if (name == "FCC_3D_potential_choice")
	{
		order = ParameterListT::Choice;

		/* choice of potentials */
		sub_lists.AddSub("harmonic");
		sub_lists.AddSub("Lennard_Jones");
		sub_lists.AddSub("Paradyn_pair");
		sub_lists.AddSub("Matsui");
	}
	else /* inherited */
		NL_E_MatT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FCC3D::NewSub(const StringT& name) const
{
	/* try to construct pair property */
	PairPropertyT* pair_prop = PairPropertyT::New(name, fMaterialSupport);
	if (pair_prop)
		return pair_prop;
	else if (name == "CB_lattice_FCC")	
		return new FCCLatticeT(0);
	else /* inherited */
		return NL_E_MatT::NewSub(name);
}

/* accept parameter list */
void FCC3D::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "FCC3D::TakeParameterList";

	/* inherited */
	NL_E_MatT::TakeParameterList(list);

	/* number of shells */
	int nshells = list.GetParameter("shells");

	/* construct pair property */
	const ParameterListT& pair_prop = list.GetListChoice(*this, "FCC_3D_potential_choice");
	fPairProperty = PairPropertyT::New(pair_prop.Name(), &(MaterialSupport()));
	if (!fPairProperty) ExceptionT::GeneralFail(caller, "could not construct \"%s\"", pair_prop.Name().Pointer());
	fPairProperty->TakeParameterList(pair_prop);

	/* check */
	fNearestNeighbor = fPairProperty->NearestNeighbor();
	if (fNearestNeighbor < kSmall)
		ExceptionT::BadInputValue(caller, "nearest bond ! (%g > 0)", fNearestNeighbor);

	/* construct lattice */
	fFCCLattice = new FCCLatticeT(nshells);
	fFCCLattice->TakeParameterList(list.GetList("CB_lattice_FCC"));

	/* construct default bond density array */
	fFullDensity.Dimension(fFCCLattice->NumberOfBonds());
	fFullDensity = 1.0;
		
	/* compute the (approx) cell volume */
	double cube_edge = fNearestNeighbor*sqrt(2.0);
	fAtomicVolume = cube_edge*cube_edge*cube_edge/4.0;

	/* compute stress-free dilatation */
	double stretch = ZeroStressStretch();
	fNearestNeighbor *= stretch;
	cube_edge = fNearestNeighbor*sqrt(2.0);
	fAtomicVolume = cube_edge*cube_edge*cube_edge/4.0;

	/* reset the continuum density (4 atoms per unit cell) */
	fDensity = fPairProperty->Mass()/fAtomicVolume;
}

/* return a reference to the bond lattice */
const BondLatticeT& FCC3D::BondLattice(void) const {
	if (!fFCCLattice) ExceptionT::GeneralFail("FCC3D::BondLattice", "pointer not set");
	return *fFCCLattice;
}

/*************************************************************************
 * Protected
 *************************************************************************/

void FCC3D::ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli)
{	
	fFCCLattice->ComputeDeformedLengths(E);
	const dArrayT& bond_length = fFCCLattice->DeformedLengths();

	/* fetch function pointers */
	PairPropertyT::ForceFunction force = fPairProperty->getForceFunction();
	PairPropertyT::StiffnessFunction stiffness = fPairProperty->getStiffnessFunction();

	/* bond density */
	const double* density = fFullDensity.Pointer();
	int nb = fFCCLattice->NumberOfBonds();
	const ElementCardT* element = MaterialSupport().CurrentElement();
	if (element && element->IsAllocated()) {
		const dArrayT& d_array = element->DoubleData();
		density = d_array.Pointer(CurrIP()*nb);
	}
	
	/* sum over bonds */
	moduli = 0.0; 
	double R4byV = fNearestNeighbor*fNearestNeighbor*fNearestNeighbor*fNearestNeighbor/fAtomicVolume;
	for (int i = 0; i < nb; i++)
	{
		double ri = bond_length[i]*fNearestNeighbor;
		double coeff = (*density++)*(stiffness(ri, NULL, NULL) - force(ri, NULL, NULL)/ri)/ri/ri;
		fFCCLattice->BondComponentTensor4(i, fBondTensor4);
		moduli.AddScaled(R4byV*coeff, fBondTensor4);
	}
	
	/* symmetric */
	moduli.CopySymmetric();
}

/* 2nd Piola-Kirchhoff stress vector */
void FCC3D::ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2)
{
#if __option(extended_errorcheck)
	if (E.Rows() != PK2.Rows() ||
	   (E.Rows() != 2 && E.Rows() != 3))
	   ExceptionT::GeneralFail("FCC3D::ComputePK2");
#endif

	/* bond density */
	const double* density = fFullDensity.Pointer();
	int nb = fFCCLattice->NumberOfBonds();
	const ElementCardT* element = MaterialSupport().CurrentElement();
	bool keep_full_density = MaterialSupport().RunState() == GlobalT::kWriteOutput && fFullDensityForStressOutput;
	if (!keep_full_density && element && element->IsAllocated()) {
		const dArrayT& d_array = element->DoubleData();
		density = d_array.Pointer(CurrIP()*nb);
	}

	/* fetch function pointer */
	PairPropertyT::ForceFunction force = fPairProperty->getForceFunction();

	/* lattice properties */
	dArrayT& bond_length = fFCCLattice->DeformedLengths();
	const dArray2DT& bonds = fFCCLattice->Bonds();
	double R2byV = fNearestNeighbor*fNearestNeighbor/fAtomicVolume;

	double* pC = fFCCLattice->Stretch().Pointer();
	const double* pE = E.Pointer();
	double* pPK2 = PK2.Pointer();
	
	/* plane strain */
	if (E.Rows() == 2)
	{
		/* compute the stretch */
		pC[0] = 2.0*pE[0] + 1.0;
		pC[1] = 2.0*pE[1] + 1.0;
		pC[2] = 1.0;
		pC[3] = 0.0; /* 23 */
		pC[4] = 0.0; /* 13 */
		pC[5] = 2.0*pE[2]; /* 12 */	

		/* initialize */
		pPK2[0] = 0.0;
		pPK2[1] = 0.0;
		pPK2[2] = 0.0;

		/* sum over bonds */
		for (int i = 0; i < nb; i++)
		{
			/* bond vector */
			const double* R = bonds(i);
	
			/* deformed bond length */
			double& l = bond_length[i];
			l = sqrt(
				(pC[0]*R[0] + pC[5]*R[1] + pC[4]*R[2])*R[0] +
				(pC[5]*R[0] + pC[1]*R[1] + pC[3]*R[2])*R[1] +
				(pC[4]*R[0] + pC[3]*R[1] + pC[2]*R[2])*R[2]
			);

			/* accumulate */
			double ri = l*fNearestNeighbor;
			double coeff = R2byV*(*density++)*force(ri, NULL, NULL)/ri;
			pPK2[0] += coeff*R[0]*R[0];
			pPK2[1] += coeff*R[1]*R[1];
			pPK2[2] += coeff*R[0]*R[1];
		}
	}
	else /* 3D */ 
	{
		/* compute the stretch */
		pC[0] = 2.0*pE[0] + 1.0;
		pC[1] = 2.0*pE[1] + 1.0;
		pC[2] = 2.0*pE[2] + 1.0;
		pC[3] = 2.0*pE[3]; /* 23 */
		pC[4] = 2.0*pE[4]; /* 13 */
		pC[5] = 2.0*pE[5]; /* 12 */	

		/* initialize */
		pPK2[0] = 0.0;
		pPK2[1] = 0.0;
		pPK2[2] = 0.0;
		pPK2[3] = 0.0;
		pPK2[4] = 0.0;
		pPK2[5] = 0.0;

		/* sum over bonds */
		for (int i = 0; i < nb; i++)
		{
			/* bond vector */
			const double* R = bonds(i);
	
			/* deformed bond length */
			double& l = bond_length[i];
			l = sqrt(
				(pC[0]*R[0] + pC[5]*R[1] + pC[4]*R[2])*R[0] +
				(pC[5]*R[0] + pC[1]*R[1] + pC[3]*R[2])*R[1] +
				(pC[4]*R[0] + pC[3]*R[1] + pC[2]*R[2])*R[2]
			);

			/* accumulate */
			double ri = l*fNearestNeighbor;
			double coeff = R2byV*(*density++)*force(ri, NULL, NULL)/ri;
			pPK2[0] += coeff*R[0]*R[0];
			pPK2[1] += coeff*R[1]*R[1];
			pPK2[2] += coeff*R[2]*R[2];
			pPK2[3] += coeff*R[1]*R[2];
			pPK2[4] += coeff*R[0]*R[2];
			pPK2[5] += coeff*R[0]*R[1];
		}
	}
}

/* strain energy density */
double FCC3D::ComputeEnergyDensity(const dSymMatrixT& E)
{
	fFCCLattice->ComputeDeformedLengths(E);
	const dArrayT& bond_length = fFCCLattice->DeformedLengths();

	/* fetch function pointer */
	PairPropertyT::EnergyFunction energy = fPairProperty->getEnergyFunction();

	/* bond density */
	const double* density = fFullDensity.Pointer();
	int nb = fFCCLattice->NumberOfBonds();
	const ElementCardT* element = MaterialSupport().CurrentElement();
	if (element && element->IsAllocated()) {
		const dArrayT& d_array = element->DoubleData();
		density = d_array.Pointer(CurrIP()*nb);
	}

	/* sum over bonds */
	double tmpSum  = 0.;	
	for (int i = 0; i < nb; i++)
	{
		double r = bond_length[i]*fNearestNeighbor;
		tmpSum += (*density++)*energy(r, NULL, NULL);
	}
	tmpSum /= fAtomicVolume;
	
	return tmpSum;
}

/* return the equitriaxial stretch at which the stress is zero */
double FCC3D::ZeroStressStretch(void)
{
	const char caller[] = "FCC3D::ZeroStress";

	int nsd = 3;
	dSymMatrixT E(nsd), PK2(nsd);
	dMatrixT C(dSymMatrixT::NumValues(nsd));

	E = 0.0;
	ComputePK2(E, PK2);
	
	/* Newton iteration */
	int count = 0;
	double error, error0;
	error = error0 = fabs(PK2(0,0));
	while (count++ < 10 && error > kSmall && error/error0 > kSmall)
	{
		ComputeModuli(E, C);
		double dE = -PK2(0,0)/(C(0,0) + C(0,1) + C(0,2));
		E.PlusIdentity(dE);
		
		ComputePK2(E, PK2);
		error = fabs(PK2(0,0));
	}

	/* check convergence */
	if (error > kSmall && error/error0 > kSmall) {
		cout << "\n " << caller << ":\n";
		cout << " E =\n" << E << '\n';
		cout << " PK2 =\n" << PK2 << endl;
		ExceptionT::GeneralFail(caller, "failed to find stress-free state");
	}
	
	return sqrt(2.0*E(0,0) + 1.0);
}
