/* $Id: ZBSolverT.h,v 1.1 2007/11/08 22:42:23 hspark Exp $ */
/* created: hspark (11/8/2007) */
#ifndef _ZB_SOLVER_T_H_
#define _ZB_SOLVER_T_H_

/* base class */
#include "ParameterInterfaceT.h"

/* direct members */
#include "iArray2DT.h"
//#include "LengthsAndAnglesT.h"
//#include "SWDataT.h" //TEMP

namespace Tahoe {

/* forward declaration */
class ThermalDilatationT;

class ZBSolverT: public ParameterInterfaceT
{
public:

	/* potential type */
	enum PotentialTypeT {kSW = 0,
                       kPTHT = 1,
                    kZB = 2};

	/** constructor */
	ZBSolverT(const ThermalDilatationT* thermal);

	/** destructor */
	~ZBSolverT(void);

	/** \name constitutive properties */
	/*@{*/
	/** moduli - C_IJKL */
	void SetModuli(const dMatrixT& CIJ, dArrayT& Xsi, dMatrixT& moduli);

	/** stress - S_IJ (2nd PK) */
	void SetStress(const dMatrixT& CIJ, dArrayT& Xsi, dMatrixT& stress);

	/** strain energy density */
	double StrainEnergyDensity(const dMatrixT& CIJ, dArrayT& Xsi);
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;
 
  	/** information about subordinate parameter lists */
	//virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	//virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:

	/* Minimize the energy wrt Xsi using the initial value passed */
	void Equilibrate(const dMatrixT& CIJ, dArrayT& Xsi);

	/* set free dof - triggers recomputation */
	void SetdXsi(const dMatrixT& CIJ, const dArrayT& Xsi);

private:

	/* mod CB flag */
	bool fEquilibrate;

	/* pairs for 3-body potentials */
	const iArray2DT fPairs;

	/** lattice geometry */
	//LengthsAndAnglesT* fGeometry;
	
	/** thermal dilatation */
	const ThermalDilatationT* fThermal;

	/* derivatives wrt. Xsi */
	dArrayT		dXsi;
	dMatrixT	dXsidXsi;

	/* derivatives wrt. C */
	dMatrixT	dCdC_hat;

	/* mixed derivatives wrt. C and Xsi */
	dMatrixT	dCdXsi_hat;
	dMatrixT	fTempRank4;
	dMatrixT	fTempMixed;

#if 0	
	
	/* work space */
	ArrayT<nArrayT<double>*>  fMatrices; //linear combo
	dMatrixT	fMat1, fMat2;
	dMatrixT	fGradl_i;
//	dSymMatrixT	fSymMat1;
	dMatrixT	fGradl_C;
#endif

	dArrayT		fVec;
	dMatrixT    fMat1;	

	/** unit cell coordinates (reference). X,Y,Z coordinates of the unit
	 * cell atoms as column vectors */
	dMatrixT fUnitCellCoords;
	
	/** \name user-defined parameters */
	/*@{*/
	/** parameter vector to pass to C code */
	dArrayT fParams;
	
	/** lattice parameter for ZB crystals */
	double f_a0;
	
	/** Atomic volume */
	double f_omega0;
	
	/** Part of energy scaling for both attractive and repulsive terms */
	double f_D0;
	
	/** Part of energy scaling for both attractive and repulsive terms */
	double f_S0;
	
	/** Equilibrium dimer length */
	double f_r0;
	
	/** Mass */
	double fMass;
	
	/** Bond order exponent scaling term */
	double f_beta;
	
	/** Bond order parameter coefficient 1. */
	double f_gamma;
	
	/** Bond order parameter coefficient 2. */
	double f_c;
	
	/** Bond order parameter coefficient 3. */
	double f_d;
	
	/** Bond order parameter coefficient 4. */
	double f_h;
	
	/** Cutoff function length parameter 1. When atoms i & j are different species R = sqrt(R_i * R_j).*/
	double f_R;
	
	/** Cutoff function length parameter 2. When atoms i & j are different species D = sqrt(D_i * D_j).*/
	double f_cut;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _ZB_SOLVER_T_H_ */
