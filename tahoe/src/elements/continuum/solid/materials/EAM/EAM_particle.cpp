/* $Id: EAM_particle.cpp,v 1.6 2011/12/01 21:11:37 bcyansfn Exp $ */
/* created: hspark(02/25/2004) */
#include "EAM_particle.h"
#include <iostream> //TEMP
#include "CBLatticeT.h"

/* EAM property types */
#include "ParadynEAMT.h"

using namespace Tahoe;

/* constructor */
EAM_particle::EAM_particle(CBLatticeT& lattice, const StringT& param_file): 
	fLattice(lattice),
	fEAMProperty(NULL),
	fPairEnergy(NULL),
	fPairForce(NULL),
	fPairStiffness(NULL),
	fEmbedEnergy(NULL),
	fEmbedForce(NULL),
	fEmbedStiffness(NULL),
	fEDEnergy(NULL),
	fEDForce(NULL),
	fEDStiffness(NULL),
	fRepRho(3)
{
	/* construct EAM property - only the Paradyn EAM potentials are implemented */
	fEAMProperty = new ParadynEAMT(param_file);

	/* cache function pointers */
	fPairEnergy    = fEAMProperty->getPairEnergy();
	fPairForce     = fEAMProperty->getPairForce();
	fPairStiffness = fEAMProperty->getPairStiffness();
	fEDEnergy      = fEAMProperty->getElecDensEnergy();
	fEDForce       = fEAMProperty->getElecDensForce();
	fEDStiffness   = fEAMProperty->getElecDensStiffness();
	fEmbedEnergy   = fEAMProperty->getEmbedEnergy();
	fEmbedForce    = fEAMProperty->getEmbedForce();
	fEmbedStiffness = fEAMProperty->getEmbedStiffness();
	
	/* set lattice parameter */
	fLatticeParameter = fEAMProperty->GetLatticeParameter();
	
	/* set the atomic mass */
	fMass = fEAMProperty->Mass();
}

/* Destructor */
EAM_particle::~EAM_particle(void) {
	delete fEAMProperty;
}

/* set "glue" functions and dimension work space */
void EAM_particle::Initialize(int nsd, int numbonds)
{
	/* dimension work space */
	int nstrs = dSymMatrixT::NumValues(nsd);
	fBondTensor4.Dimension(nstrs);
	fAmn.Dimension(numbonds);
	fBondTensor2.Dimension(nstrs);
	fTensor2Table.Dimension(numbonds, nstrs);
	fBond1.Dimension(numbonds);
	fBond2.Dimension(numbonds);
	fBond3.Dimension(numbonds);
	
	/* Dimension surface cluster work space */
	const dArrayT& rb = fLattice.DeformedBulk();
	const dArrayT& rs1 = fLattice.DeformedSurf1();
	const dArrayT& rs2 = fLattice.DeformedSurf2();
	fBond4.Dimension(rb.Length());
	fBond5.Dimension(rs1.Length());
	fBond6.Dimension(rs2.Length());
	fBondTensor2b.Dimension(nstrs);
	//fIntType.Dimension(6,2);
	fIntType.Dimension(9,2);
}

/*
* Compute unit strain energy density:
*
*     unit strain energy = energy/atom
*/
double EAM_particle::ComputeUnitEnergy(void)
{
	const dArrayT& r = fLattice.DeformedLengths();
	const iArrayT& counts = fLattice.BondCounts();

	double rho = 0.0;
	double energy = 0.0;
	
	const int* pcount = counts.Pointer();
	int nb = r.Length();
	for (int i = 0; i < nb; i++)
	{
		int ci = *pcount++;

		double  ri = r[i];
		double a = fPairEnergy(ri, NULL, NULL);
		double phi = a*a/ri;
		double rho = fEDEnergy(ri, NULL, NULL);

		rho    += ci*rho;
		energy += ci*0.5*phi;
	}
	
	energy += fEmbedEnergy(rho, NULL, NULL);
	return energy;
}

/*
* Compute unit 2nd PK stress:
*
*     unit 2nd PK stress = SIJ*(volume per cell/atoms per cell)
*/
void EAM_particle::ComputeUnitStress(dSymMatrixT& stress)
{
	const dArrayT& r = fLattice.DeformedLengths();
	const iArrayT& counts = fLattice.BondCounts();

	/* total atomic density */
	double rho = TotalElectronDensity();
	double dFdrho = fEmbedForce(rho, NULL, NULL);	// checks out

	/* assemble stress */
	stress = 0.0;
	int nb = r.Length();	
	for (int i = 0; i < nb; i++)	// fNumBonds, ri, ci check out
	{
		double ri = r[i];
		int    ci = counts[i];	
		double a = fPairEnergy(ri, NULL, NULL);
		double Potential = a*a/ri;
		double DPotential = 2.0*fPairForce(ri, NULL, NULL)*a/ri-Potential/ri;
		double DDensity = fEDForce(ri, NULL, NULL);	
		double coeff = (1.0/ri)*ci*(0.5*DPotential + dFdrho*DDensity);
		fLattice.BondComponentTensor2(i,fBondTensor2);
		stress.AddScaled(coeff,fBondTensor2);
	}
}

/*
* Compute unit 2nd PK surface stress - used by surface CB solver
*
*     unit 2nd PK surface stress = SIJ*(area per cell/atoms per cell)
*/
void EAM_particle::ComputeUnitSurfaceStress(dSymMatrixT& stress)
{
	const dArrayT& r = fLattice.DeformedLengths();
	const iArrayT& counts = fLattice.BondCounts();
	const iArrayT& atom_type = fLattice.AtomTypes();	// added by HSP
	const dArray2DT& bonds = fLattice.Bonds();	// can access initial bonds
		
	/* calculate representative electron densities */
	ComputeElectronDensity();
	
	/* Create interaction table using calculated electron densities */
	fIntType(0,0) = (fRepRho[1]);	// surface1
	fIntType(0,1) = (fRepRho[1]);	// surface1
	fIntType(1,0) = (fRepRho[1]);	// surface1
	fIntType(1,1) = (fRepRho[2]);	// surface2
	fIntType(2,0) = (fRepRho[1]);	// surface1
	fIntType(2,1) = (fRepRho[0]);	// bulk
	fIntType(3,0) = (fRepRho[2]);	// surface2
	fIntType(3,1) = (fRepRho[1]);	// surface1
	fIntType(4,0) = (fRepRho[2]);	// surface2
	fIntType(4,1) = (fRepRho[2]);	// surface2
	fIntType(5,0) = (fRepRho[2]);	// surface2
	fIntType(5,1) = (fRepRho[0]);	// bulk
	fIntType(6,0) = (fRepRho[0]);	// bulk	- new additions for S3 and S4 layers
	fIntType(6,1) = (fRepRho[1]);	// surface1
	fIntType(7,0) = (fRepRho[0]);	// bulk
	fIntType(7,1) = (fRepRho[2]);	// surface2
	fIntType(8,0) = (fRepRho[0]);	// bulk
	fIntType(8,1) = (fRepRho[0]);	// bulk

	/* assemble stress */
	stress = 0.0;
	int nb = r.Length();
	dArrayT blah(2);

	for (int i = 0; i < nb; i++)
	{
		// define dFdrho_i and dFdrho_j here depending on rho_i and rho_j
		int type = atom_type[i];	// what is the interaction type?
		fIntType.RowCopy(type, blah);	// get correct electron densities for i-j interactions
		double dFdrho_i = fEmbedForce(blah[0], NULL, NULL);
		double dFdrho_j = fEmbedForce(blah[1], NULL, NULL);
		double ri = r[i];
		int    ci = counts[i];		
		double a = fPairEnergy(ri, NULL, NULL);
		double Potential = a*a/ri;
		double DPotential = 2.0*fPairForce(ri, NULL, NULL)*a/ri-Potential/ri;
		double DDensity = fEDForce(ri, NULL, NULL);
		double coeff = (1.0/ri)*ci*0.5*(DPotential + dFdrho_j*DDensity + dFdrho_i*DDensity);
		fLattice.BondComponentTensor2(i,fBondTensor2b);
		stress.AddScaled(coeff,fBondTensor2b);
	}
}
	   	    	
/*
* Compute unit material tangent moduli:
*
*   unit material tangent moduli = CIJKL*(volume per cell/atoms per cell)
*/
void EAM_particle::ComputeUnitModuli(dMatrixT& moduli)
{
	/* compute total electron density */
	double rho = TotalElectronDensity();

	/* initialize */
	moduli = 0.0;

	/* mixed bond contribution */
	FormMixedBondContribution(rho, moduli);
	
	/* single bond contribution */
	FormSingleBondContribution(rho, moduli);
}

/* calculate rho for bulk, surface1 and surface2 atoms */
void EAM_particle::ComputeElectronDensity(void)
{
	/* Get deformed lengths for bulk, surface1 and surface2 clusters */
	const dArrayT& rb = fLattice.DeformedBulk();
	const dArrayT& rs1 = fLattice.DeformedSurf1();
	const dArrayT& rs2 = fLattice.DeformedSurf2();
	
	/* Get bond counts for each cluster type */
	const iArrayT& countsb = fLattice.BulkCounts();
	const iArrayT& countss1 = fLattice.Surf1Counts();
	const iArrayT& countss2 = fLattice.Surf2Counts();
	
	/* values of electron density */
	double rhob = 0.0;
	double rhos1 = 0.0;
	double rhos2 = 0.0;

	const int* pcountb = countsb.Pointer();
	const int* pcounts1 = countss1.Pointer();
	const int* pcounts2 = countss2.Pointer();
	int nbb = rb.Length();
	int nbs1 = rs1.Length();
	int nbs2 = rs2.Length();

	/* bulk atoms */
	for (int i = 0; i < nbb; i++)
	{
		double ri = rb[i];
		double pedensb = fEDEnergy(ri, NULL, NULL);
		rhob += (*pcountb++)*(pedensb++);
	}
	
	/* surface1 atoms */
	for (int i = 0; i < nbs1; i++)
	{
		double rsurf1 = rs1[i];
		double pedenss1 = fEDEnergy(rsurf1, NULL, NULL);
		rhos1 += (*pcounts1++)*(pedenss1++);	
	}
	
	/* surface2 atoms */
	for (int i = 0; i < nbs2; i++)
	{
		double rsurf2 = rs2[i];
		double pedenss2 = fEDEnergy(rsurf2, NULL, NULL);
		rhos2 += (*pcounts2++)*(pedenss2++);
	}
	
	double blah[3] = {rhob,rhos1,rhos2};
	fRepRho = blah;
}

/*
* Compute the total electron density.
*/
double EAM_particle::TotalElectronDensity(void)
{
	const dArrayT& r = fLattice.DeformedLengths();
	const iArrayT& counts = fLattice.BondCounts();

	double rho = 0.0;
	const int* pcount = counts.Pointer();
	int nb = r.Length();
	for (int i = 0; i < nb; i++)
	{
		double ri = r[i];
		double pedensity = fEDEnergy(ri, NULL, NULL);
		rho += (*pcount++)*(pedensity);
	}
	return rho;
}

/* return the embedding force for a given electron density */
double EAM_particle::ReturnEmbeddingForce(double rho)
{
	return fEmbedForce(rho, NULL, NULL);
}

/**********************************************************************
 * Private
 **********************************************************************/

/*
* Form matrix of mixed pair potential and embedding
* energy derivatives.  NOTE: computes the UPPER triangle
* ONLY.
*/
void EAM_particle::FormMixedDerivatives(double rho)
{
	const dArrayT& r = fLattice.DeformedLengths();
	const iArrayT& counts = fLattice.BondCounts();

	double dFdrho = fEmbedForce(rho, NULL, NULL);
	double d2Fdrho2 = fEmbedStiffness(rho, NULL, NULL);

	/* form upper triangle only */
	int nb = r.Length();
	for (int j = 0; j < nb; j++)
	{
		double rj = r[j];
		int    cj = counts[j];

		double a = fPairEnergy(rj, NULL, NULL);
		double b = fPairForce(rj, NULL, NULL);
		double Potential = a*a/rj;
		double DPotential = 2.0*b*a/rj-Potential/rj;
		double DDPj = 2.0*fPairStiffness(rj, NULL, NULL)*a/rj;
		DDPj+=2.0*b*b/rj;
		DDPj-=(2.0*DPotential/rj);
		double DDDj = fEDStiffness(rj, NULL, NULL);
		double DDj = fEDForce(rj, NULL, NULL);
	
		for (int i = 0; i <= j; i++)
		{
			double Amn = 0.0;
			double ri = r[i];
			int    ci = counts[i];
					
			/* diagonal terms */
			if (i == j)
			{
				/* pair potential */
				Amn += 0.5*cj*DDPj;	
				
				/* embedding energy */
				Amn += cj*dFdrho*DDDj;
			}
			double DDensity = fEDForce(ri, NULL, NULL);
		
			/* mixed embedding energy term */
			Amn += ci*cj*d2Fdrho2*DDensity*DDj;
		
			fAmn(i,j) = Amn/(ri*rj);
		}
	}
	
	fAmn.CopySymmetric();
}	

/*
* Moduli tensor contributions.
*/
void EAM_particle::FormSingleBondContribution(double rho, dMatrixT& moduli)
{
	const dArrayT& r = fLattice.DeformedLengths();
	const iArrayT& counts = fLattice.BondCounts();

	/* Embedding energy derivative */
	double dFdrho = fEmbedForce(rho, NULL, NULL);
	int nb = r.Length();	
	for (int i = 0; i < nb; i++)
	{
		double ri = r[i];
		double a = fPairEnergy(ri, NULL, NULL);
		double Potential = a*a/ri;
		double DPotential = 2.0*fPairForce(ri, NULL, NULL)*a/ri-Potential/ri;
		double DDensity = fEDForce(ri, NULL, NULL);
		double coeff = -counts[i]*(0.5*DPotential + dFdrho*DDensity)/(ri*ri*ri);
		fLattice.BondComponentTensor4(i,fBondTensor4);		
		moduli.AddScaled(coeff,fBondTensor4);
	}
}

void EAM_particle::FormMixedBondContribution(double rho, dMatrixT& moduli)
{
	/* batch fetch */
	fLattice.BatchBondComponentTensor2(fTensor2Table);

	/* mixed bond derivatives */
	FormMixedDerivatives(rho);

/*	
	for (int n = 0; n < fNumBonds; n++)
	{
		for (int m = 0; m < fNumBonds; m++)
		{
			double Amn = fAmn(m,n);
			
			for (int j = 0; j < fModuliDim; j++)
			{
				double tempj = Amn*fTensor2Table(n,j);
			
				for (int i = 0; i <= j; i++)
					moduli(i,j) += tempj*fTensor2Table(m,i);
			}
		}
	}
*/
	
	/* reversing the loops */
	int dim = moduli.Rows();
	int nb = fLattice.NumberOfBonds();
	for (int j = 0; j < dim; j++)
		for (int i = 0; i <= j; i++)
		{
			register double cij = 0.0;
			double* pm = fTensor2Table(0) + i;
			double* pn = fTensor2Table(0) + j;
		
			for (int n = 0; n < nb; n++)
			{
				double* pAmn = fAmn(n);
				double* pm2  = pm;
				double  nj   = *pn;
			
				for (int m = 0; m < nb; m++)
				{
					cij += nj*(*pAmn)*(*pm2);

					pAmn++;
					pm2 += dim;
				}
				
				pn += dim;
			}
			
			moduli(i,j) += cij;
		}

	moduli.CopySymmetric();
}
