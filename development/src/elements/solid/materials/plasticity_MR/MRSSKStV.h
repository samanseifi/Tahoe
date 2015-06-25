/* $Id: MRSSKStV.h,v 1.8 2010/07/21 19:58:20 regueiro Exp $ */
/* created: Majid T. Manzari (04/16/2003) */
#ifndef _MR_SS_KSTV_H_
#define _MR_SS_KSTV_H_

/* base classes */
#include "SSIsotropicMatT.h"
#include "HookeanMatT.h"

namespace Tahoe {

/* forward declarations */
class MRSSNLHardT;
class SSEnhLocMatSupportT;

class MRSSKStV: public SSIsotropicMatT, public HookeanMatT
{
  public:

	/* constructor */
	MRSSKStV(void);
	
	/* destructor */
	~MRSSKStV(void);

	/* form of tangent matrix (symmetric by default) */
	virtual GlobalT::SystemTypeT TangentType(void) const;
	
	/** access strains from previous time step */
	virtual bool Need_Strain_last(void) const { return true; };
	
	/** model has history variables */
	virtual bool HasHistory(void) const { return true; };

	/* update internal variables */
	virtual void UpdateHistory(void);

	/* reset internal variables to last converged solution */
	virtual void ResetHistory(void);
	
	/** returns elastic strain (3D) */
	virtual const dSymMatrixT& ElasticStrain(
                const dSymMatrixT& totalstrain, 
				const ElementCardT& element, int ip);

	/** \name spatial description */
	/*@{*/
	/** spatial tangent modulus */
	virtual const dMatrixT& c_ijkl(void);
	virtual const dMatrixT& ce_ijkl(void);
	virtual const dMatrixT& c_perfplas_ijkl(void);

	/** Cauchy stress */
	virtual const dSymMatrixT& s_ij(void);
	
	/** initial stress */
	void InitialStress(dSymMatrixT& Stress0);

	/** return the pressure associated with the last call to 
	 * SolidMaterialT::s_ij. See SolidMaterialT::Pressure
	 * for more information. */
	virtual double Pressure(void) const { return fStress.Trace()/3.0; };
	/*@}*/

	/* returns the strain energy density for the specified strain */
	virtual double StrainEnergyDensity(void);

	/* returns the number of variables computed for nodal extrapolation
	 * during for element output, ie. internal variables */
	virtual int  NumOutputVariables(void) const;
	virtual void OutputLabels(ArrayT<StringT>& labels) const;
	virtual void ComputeOutput(dArrayT& output);
	
	/*
	* Test for localization using "current" values for Cauchy
	* stress and the spatial tangent moduli. Returns true if the
	* determinant of the acoustic tensor is negative and returns
	* the normals and slipdirs. Returns false if the determinant is positive.
	*/
	bool IsLocalized(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipdirs, 
					AutoArrayT <double> &detAs, AutoArrayT <double> &dissipations_fact);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;
	
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;
	
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/
    
    /** set flag for writing iteration info */
	void GetIterationInfo(bool get_iters, int loc_iters);
	
protected:

	/* set modulus */
	virtual void SetModulus(dMatrixT& modulus); 
	int loccheck;
	
	// element localization flag
	int element_locflag;
	
	bool fGetItersInfo; // write global and local iteration info 
	
	/* 
	pointer to material support; does not actually return stresses from this class;
	stresses are calculated post-localization in SmallStrainEnhLocT
	*/
	const SSEnhLocMatSupportT* fSSEnhLocMatSupport;
 
  private:
  
  	/** pressure sensitive plasticity with nonlinear hardening and localization*/
	MRSSNLHardT* fMR;
  
  	/* return values */
  	dSymMatrixT	fStress;
  	dMatrixT	fModulus, fModulusCe;
    dMatrixT    fModulusPerfPlas;

};

} // namespace Tahoe 
#endif /* _MR_SS_KSTV_H_ */
