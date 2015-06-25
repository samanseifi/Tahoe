/* 3-invariant, single-surface dilation/compaction plasticity model
 * with isotropic and kinematic hardening
 * Implemented 8/02 Craig Foster
 */


#ifndef _FOSSUM_SS_ISOT_H_
#define _FOSSUM_SS_ISOT_H_

#include "SSIsotropicMatT.h"
#include "HookeanMatT.h"
#include "SSHookeanMatT.h"
#include "SpectralDecompT.h"

/* direct members */
#include "dMatrixT.h"
#include "dArrayT.h"
#include "LAdMatrixT.h"

namespace Tahoe {

/* forward declarations */
class SSEnhLocMatSupportT;

class FossumSSIsoT: public SSIsotropicMatT, public HookeanMatT //, public ParameterInterfaceT
{

public:

	/* constructor */
	FossumSSIsoT(void);

	/* destructor */
	virtual ~FossumSSIsoT(void);

	/* required parameter flags */
	virtual bool Need_Strain_last(void) const {return true;};

	/** material has history variables */
	virtual bool HasHistory(void) const { return true; };

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineParameters(ParameterListT& list) const;

	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;
	
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/*
	* Returns the value of the yield function given the
	* Cauchy stress vector and state variables, where alpha is
	* the deviatoric stress-like internal state variable
	*/
	double YieldCondition(const dSymMatrixT& stress, const double kappa, dSymMatrixT& alpha);
	double YieldFn(double I1, double J2, double J3, double kappa);

protected:
  
	double fA;         // Material Parameter for F_f part of yield fn
	double fB;
	double fC;
	double fTheta;
      
	double fR;         // Ratio of principal radii of ellipse F_c
	double fKappa0;    // Determines starting position of hardening cap

	double fW, fD1, fD2; // Isotropic cap hardening function parameters

	double fCalpha;    // determines rate of growth of back stress tensor

	double fPsi;       // Ratio of failure strength in tension to f. s. in compr.
	double fN;         // offset from initial yield to failure
	double fFluidity;   //fluidity parameter, relation time = fFluidity/(2*fmu)
	
	double fL;         // plastic potential parameter analogous to fB
	double fPhi;        // plastic potential parameter analogous to fTheta
	double fQ;         // plastic potential parameter analogous to fR
	bool fFossumDebug;

	/* form of tangent matrix (symmetric by default) */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/* update internal variables */
	virtual void UpdateHistory(void);

	/* reset internal variables to last converged solution */
	virtual void ResetHistory(void);

	/* moduli */
	virtual const dMatrixT& c_ijkl(void);
	virtual const dMatrixT& ce_ijkl(void);
	virtual const dMatrixT& con_ijkl(void);
	virtual const dMatrixT& c_perfplas_ijkl(void);
	virtual const dMatrixT& con_perfplas_ijkl(void);

	/* stress, including possibility of localized deformation */
	virtual const dSymMatrixT& s_ij(void);
	
	/* stress for return map */
	const dSymMatrixT& sigma_ij(void);

	/** return the pressure associated with the last call to 
	* SolidMaterialT::s_ij. See SolidMaterialT::Pressure
	* for more information. */
	virtual double Pressure(void) const { return fStress.Trace()/3.0; };

	/* returns the strain energy density for the specified strain */
	virtual double StrainEnergyDensity(void);

	/* returns the number of variables computed for nodal extrapolation
	* during for element output, ie. internal variables */
	virtual int  NumOutputVariables(void) const;
	virtual void OutputLabels(ArrayT<StringT>& labels) const;
	virtual void ComputeOutput(dArrayT& output);
	
	/* test for localization */
	bool IsLocalized(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipdirs, 	 
			AutoArrayT <double> &detAs, AutoArrayT <double> &dissipations_fact);
			
	/* test for discontinuous bifurcation (weak discontinuity) */
	bool IsLocalized_DB(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipdirs, 	 
			AutoArrayT <double> &detAs, AutoArrayT <double> &dissipations_fact);		

protected:

	/* set modulus */

	virtual void SetModulus(dMatrixT& modulus); 
	int loccheck, element_locflag;
 
private:
  
	/* return values */
	dSymMatrixT fStress; // stress choice between continuum model and localization model
	dSymMatrixT fSigma; // stress from itegration
	dSymMatrixT fStrain;
	dSymMatrixT fdGdSigma; // direction of plastic flow, not normalized
	dMatrixT fModulus;
	dMatrixT fModulusCe;
	dMatrixT fModulusPerfPlas;
	dMatrixT fModulusContinuum;
	dMatrixT fModulusContinuumPerfPlas;
	dArrayT fHardMod;
	
	// pointer to material support
	const SSEnhLocMatSupportT* fSSEnhLocMatSupport;


/*--------------------------------------------------------------------*/
// hardening functions

protected:

	/* status flags */
	enum LoadingStatusT {kIsPlastic = 0,
						kIsElastic = 1,
						kReset = 3}; // indicate not to repeat update

	/* returns elastic strain (3D) */
	virtual const dSymMatrixT& ElasticStrain(const dSymMatrixT& totalstrain, const ElementCardT& element, int ip);
                        
	/* return the correction to moduli due to plasticity (if any)
	*
	* Note: Return mapping occurs during the call to StressCorrection.
	*       The element passed in is already assumed to carry current
	*       internal variable values */
	const dMatrixT& ModuliCorrection(const ElementCardT& element, int ip); 

	/* Modulus for checking strong discontinuity bifurcation */
	const dMatrixT& ModuliCorrPerfPlas(const ElementCardT& element, int ip);
	
	/* plastic flow direction */
	const dSymMatrixT& PlasFlowDir(void);
	
	/* hardening modulus */
	const dArrayT& HardMod(void);

	/* return a pointer to a new plastic element object constructed with
	* the data from element */
	void AllocateElement(ElementCardT& element);

	enum InternalVariablesT {kkappa = 0,  // stress-like internal state variable
				 kdeltakappa = 1,  //increment of kappa 
				 kdgamma = 2}; // consistency parameter
							
	/* element level data */
	void Update(ElementCardT& element);
	void Reset(ElementCardT& element);

private:

	/* load element data for the specified integration point */
	void LoadData(const ElementCardT& element, int ip);

protected:

	/* element level internal state variables */
	dSymMatrixT fPlasticStrain; //total plastic strain (deviatoric and volumetric)
	dSymMatrixT fUnitNorm;      //unit normal to the yield surface
	dArrayT     fInternal;      //internal variables

private:

	/* number of integration points */
	int fNumIP;

	/* material parameters **/
	double fmu;
	double flambda;
	double fkappa;
	double fX_H;
	double fX;
	dSymMatrixT fBackStress;
	dSymMatrixT fDeltaAlpha;  

	/* spectral decomp parameters*/
	SpectralDecompT spectre;
	ArrayT<dSymMatrixT> m;
	dArrayT principalEqStress;

	/* return values */
	dSymMatrixT fElasticStrain;
	dMatrixT fModuliCorr;
	dMatrixT fModuliCorrPerfPlas;
                
	/* work space */
	dSymMatrixT One;  
	double fTimeFactor;

	int &fKappaCapped;
	int fKappaDummy;
        
private:

	/* auxiliary functions to yield function */
	double YieldFnGamma(double J2, double J3);
	double YieldFnFfMinusN(double I1);
	double YieldFnFc(double I1, const double kappa);
	int HeavisideFn(double arg);
	double Xfn(const double kappa);
	double YieldFnFf(double I1);
	
	double PlasticPotGfMinusN(double I1);
	double PlasticPotGc(double I1, const double kappa);
	double X_G(const double kappa);
	double PlasticPotGf(double I1);		

	/* auxiliaries to s_ij */
	bool StressPointIteration(double initialYieldCheck, dArrayT& iterationVars, dSymMatrixT workingBackStress, double workingKappa);
	bool ResidualIsConverged(dArrayT& residual, dArrayT& residual0);
	dArrayT CapKappa(const dArrayT &residual, const LAdMatrixT &dRdX, const double kappa);
	dArrayT CondenseAndSolve(const LAdMatrixT& dRdX, const dArrayT& residual);

	double ElasticConstant(int i, int j);	
	/* hardening functions and necessary derivatives */
	double Galpha(dSymMatrixT alpha);
	double KappaHardening(double I1, double kappa);
	double dfdDevStressA (double I1, double J2, double J3, double sigmaA);
	double dfdSigmaA(double I1, double J2, double J3, double sigmaA, double kappa);
	double dGdSigmaA(double I1, double J2, double J3, double sigmaA, double kappa);
	double dfdI1(double I1, double kappa);
	double dGdI1(double I1, double kappa);
	double dFfdI1(double I1);
	double dGfdI1(double I1);
	double dFcdI1(double I1, double kappa);
	double dGcdI1(double I1, double kappa);
	double dfdJ2(double J2, double J3);
	double dGammadJ2 (double J2, double J3);
	double dfdJ3(double J2, double J3);
	double dPlasticVolStraindX(double kappa);
	//double dXdKappa(double kappa);
	double dX_GdKappa(double kappa);


    /* Matrix for stress point NR iteration */
	LAdMatrixT FormdRdX(double I1, double J2, double J3, dArrayT principalEqStress, double workingKappa, dSymMatrixT workingStress, dSymMatrixT workingBackStress, double dGamma, ArrayT<dSymMatrixT> m);

    /* derivatives for FormdRdX */
	int KroneckerDelta (int A, int B);
	double d2GdSigmaBdSigmaC (double I1, double J2, double J3, double principalEqStressA, double principalEqStressB, int A, int B, double kappa);
	double d2GdDevStressdSigmaB (double I1, double J2, double J3, double principalEqStressA, double principalEqStressB, int A, int B);
	double d2GdI1dI1(double I1, double kappa);
	double d2GfdI1dI1(double I1);
	double d2GcdI1dI1(double I1, double kappa);
	double d2GdJ2dJ2 (double J2, double J3);
	double d2GammadJ2dJ2(double J2, double J3);
	double d2GdJ2dJ3 (double J2, double J3);
	double d2GammadJ2dJ3 (double J2);
	double dGammadJ3(double J2);
	double d2GdJ3dJ3 (double J2, double J3);
	double d2GdSigmaCdKappa (double I1, double kappa);
	double d2GcdI1dKappa(double I1, double kappa);
	double dGalphadAlphaB (dSymMatrixT alpha, dArrayT principalEqStress, int B, ArrayT<dSymMatrixT> m);
	double d2GdI1dKappa (double I1, double kappa);
	double dFcdKappa (double I1, double kappa);
	double dGcdKappa (double I1, double kappa);
	double d2X_GdKappadKappa( double kappa);
	double d2PlasticVolStraindXdX(double kappa);
	double dfdKappa(double I1, double kappa);
	double InnerProduct(dSymMatrixT A, dSymMatrixT B);

	/*tensor-valued derviatives for consistent tangent */
	dMatrixT D2GdSigmadSigma(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m);
	dSymMatrixT DfdSigma(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m);
	dSymMatrixT DGdSigma(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m);
	dSymMatrixT DfdAlpha(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m);
	dArrayT Hardening(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m, dSymMatrixT alpha);
};

const double sqrt3__ = sqrt(3.0);

/* Auxiliary Functions for yield function */
inline double FossumSSIsoT::YieldFnGamma(double J2, double J3)
{
	if (J2 <= 0.0)
		return 1.0;   //limit as s_ij -> 0

	double sin3Beta = -3.0*sqrt3__*J3/(2.0*J2*sqrt(J2));

	return .5*(1 + sin3Beta + 1/fPsi*(1.0 - sin3Beta)); 
}

inline double FossumSSIsoT::dGammadJ2 (double J2, double J3)
{
	return 9 * sqrt3__ * J3 * ( 1 - 1/fPsi) / (8 * J2*J2*sqrt(J2));
}

inline double FossumSSIsoT::dfdJ2(double J2, double J3)
{
	return YieldFnGamma (J2, J3) * YieldFnGamma (J2, J3) + 2 * J2 * YieldFnGamma (J2, J3)* dGammadJ2 (J2, J3);
}


} // namespace Tahoe 

#endif/*_FOSSUM_SS_ISOT_H_*/
