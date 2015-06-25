/* $Id: FSFiberMatViscT.h,v 1.6 2010/06/24 13:49:17 thao Exp $ */
/* created: paklein (06/09/1997) */
#ifndef _FD_FIBVISC_MAT_T_H_
#define _FD_FIBVISC_MAT_T_H_

/* base class */
#include "FSFiberMatT.h"

#include "C1FunctionT.h"
#include "PotentialT.h"

namespace Tahoe {

class FSFiberMatViscT: public FSFiberMatT
{
public:

	/** constructor */
	FSFiberMatViscT(void);
	~FSFiberMatViscT(void);

	/** required parameter flag. Indicates whether the constitutive model
	 * requires the deformation gradient from the previous time increment.
	 * \return false by default. */

	virtual GlobalT::SystemTypeT TangentType(void) const{return GlobalT::kNonSymmetric;};
	virtual bool Need_F_last(void) const { return false; };

	/** \name spatial description */
	/*@{*/
	/** spatial tangent modulus */
	virtual const dMatrixT& c_ijkl(void);

	/** Cauchy stress */
	virtual const dSymMatrixT& s_ij(void);

	/* material description */
	virtual const dMatrixT& C_IJKL(void); // material tangent moduli
	virtual const dSymMatrixT& S_IJ(void); // PK2 stress

	/** return true if the material has history variables */
	virtual bool HasHistory(void) const { return true; };
	
	/*Initialize history variable*/
	virtual bool NeedsPointInitialization(void) const {return true;}; 
	virtual void PointInitialize(void);              

	/* update/reset internal variables */
	virtual void UpdateHistory(void); // element at a time
	virtual void ResetHistory(void);  // element at a time

	void Load(ElementCardT& element, int ip);
	void Store(ElementCardT& element, int ip);
	
		/* Dimension internal state variables*/
	/*derived class must call RGViscoelaticity::SetStateVariables(fNumProcess)
	  to dimension internal state variable arrays if fNumProcess > 1 (default value)*/
	void SetStateVariables (const int numprocess);

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;
	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	
protected:

	/**************************************** matrix ***********************************************************/
	/*computes eq isotropic matrix stress*/
	virtual void ComputeMatrixStress (const dSymMatrixT& Stretch, dSymMatrixT& Stress);

	/*computes eq matrix moduli*/
	virtual void ComputeMatrixMod (const dSymMatrixT& Stretch, dSymMatrixT& Stress, dMatrixT& Mod);

	/*subsequent derived classes must define the following functions*/
	/*compute neq. values*/
	/*local newton loop for viscous stretch tensor of matrix*/ 
	virtual void ComputeMatrixCv(const dSymMatrixT& C, const dSymMatrixT& Cv_last, dSymMatrixT& Cv, const int process_index);

	/*computes isotropic matrix stress*/
	virtual void ComputeMatrixStress (const dSymMatrixT& Stretch, const dSymMatrixT& Stretch_v, 
				dSymMatrixT& Stress, const int process_index, const int fillmode = dSymMatrixT::kOverwrite);

	/*computes matrix moduli*/
	virtual void ComputeMatrixMod (const dSymMatrixT& Stretch, const dSymMatrixT& Stretch_v, dSymMatrixT& Stress,
				dMatrixT& Mod, const int process_index, const int fillmode = dSymMatrixT::kOverwrite);

	/********************************************** fiber ******************************************************/
	/*computes eq fiber stress in local frame*/
	virtual void ComputeFiberStress (const dSymMatrixT& Stretch, dSymMatrixT& Stress);
	
	/*computes eq fiber moduli in local frame*/
	virtual void ComputeFiberMod (const dSymMatrixT& Stretch, dSymMatrixT& dSymMatrixT, dMatrixT& Mod);

	/*defined pure virtual functions in FSFiberMatT*/
	/*computes fiber stress in local frame*/
	virtual void ComputeFiberStress (const dSymMatrixT& Stretch, const dSymMatrixT& Stretch_v, dSymMatrixT& Stress, 
				const int process_index) = 0;
	
	/*computes  moduli in local frame*/
	virtual void ComputeFiberMod (const dSymMatrixT& Stretch, const dSymMatrixT& Stretch_v, dSymMatrixT& Stress, 
				dMatrixT& Mod,  const int process_index) = 0;
				
	/*compute the algorithmic moduli dSNEQ/dCv deltaCv/deltadC in local fiber coord sys
	virtual void ComputeCalg (const dSymMatrixT& Stretch, const dSymMatrixT& Stretch_v,  
				dMatrixT& Calg, const int process_index) = 0; */
				
	/*local newton loop for viscous stretch tensor*/ 
	virtual void Compute_Cv(const dSymMatrixT& Stretch, const dSymMatrixT& Stretch_v_last, dSymMatrixT& Stretch_v,
				const int process_index) = 0;

protected:

	/* number of nonequilibrium processes*/
	/* must be set in derived classes before TakeParameterList is called*/
	/* default value is 1*/
	
	int fNumFibProcess;
	int fNumMatProcess;

	/*internal state variables. Dimension numprocess<nsd x nsd>*/
	ArrayT<dSymMatrixT> fC_v;
	ArrayT<dSymMatrixT> fC_vn;

  	/*constant viscosities for now*/
	ArrayT<C1FunctionT*> fVisc_m;
	
	/*moduli for Mooney Rivlin Potential, n_process x 2 (c1, c2)*/
	ArrayT<PotentialT*> fPot_m;

	/*algorithmic modulus*/
	dMatrixT fCalg;

	/*number of state variables*/
	int fnstatev;
	
	/* internal state variables array*/
	dArrayT fstatev;

	/*viscous stretch in plane of fiber*/
	dSymMatrixT fFiberStretch_v;
	dSymMatrixT fFiberStretch_vn;

 	/* work space needed to compute matrix viscoelastic behavior*/
	/* spectral operations for integration of matrix elastic stretch*/
	SpectralDecompT fSpectralDecompSpat;

	dSymMatrixT fInverse;
	dMatrixT fModMat;
	dMatrixT fMod3;

	dSymMatrixT fb;
	dArrayT fEigs;
	dArrayT	ftau;
	dSymMatrixT fdtau_dep;

  	dMatrixT fiK_m;
	dMatrixT fCalg_m;
};

/*equilibrium values*/
inline void FSFiberMatViscT::ComputeMatrixStress (const dSymMatrixT& Stretch, dSymMatrixT& Stress)
{ 
	ComputeMatrixStress(Stretch, Stretch, Stress, -1);
}

inline void FSFiberMatViscT::ComputeMatrixMod (const dSymMatrixT& Stretch, dSymMatrixT& Stress, dMatrixT& Mod)
{ 
	ComputeMatrixMod(Stretch, Stretch, Stress, Mod, -1);
}

inline void FSFiberMatViscT::ComputeFiberStress (const dSymMatrixT& Stretch, dSymMatrixT& Stress)
{	
	ComputeFiberStress(Stretch, Stretch, Stress, -1);
}

inline void FSFiberMatViscT::ComputeFiberMod (const dSymMatrixT& Stretch, dSymMatrixT& Stress, dMatrixT& Mod)
{ 
	ComputeFiberMod(Stretch, Stretch, Stress, Mod, -1);
}

} /* namespace Tahoe */

#endif /* _FD_STRUCT_MAT_T_H_ */
