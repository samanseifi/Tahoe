/* $Id: DPSSLinHardLocT.h,v 1.6 2006/06/14 18:44:12 regueiro Exp $ */
/* created: myip (06/01/1999)                                      */

/*  
 * Interface for Drucker-Prager, nonassociative, small strain,
 * pressure-dependent plasticity model with linear isotropic hardening
 * and localization.
 *
 *	Note: all calculations are peformed in 3D.
 */

#ifndef _DP_SS_LIN_HARD_LOC_T_H_
#define _DP_SS_LIN_HARD_LOC_T_H_

/* base class */
#include "DPPrimitiveLocT.h"

/* direct members */
#include "dSymMatrixT.h"
#include "dMatrixT.h"
#include "dArrayT.h"

#include "DPSSKStVLoc.h"

namespace Tahoe {

/* forward declarations */
class ElementCardT;

class DPSSLinHardLocT: public DPPrimitiveLocT
{
public:

	/* constructor */
	DPSSLinHardLocT(int num_ip, double mu, double lambda);


	/* status flags */
	enum LoadingStatusT {kIsPlastic = 0,
						kIsElastic = 1,
						kIsLocalized = 2,
                        kReset = 3}; // indicate not to repeat update

	/* returns elastic strain (3D) */
	virtual const dSymMatrixT& ElasticStrain(const dSymMatrixT& totalstrain, 
		const ElementCardT& element, int ip);
			
	/* return correction to stress vector computed by mapping the
	 * stress back to the yield surface, if needed */
	const dSymMatrixT& StressCorrection(const dSymMatrixT& trialstrain, 
		ElementCardT& element, int ip, double dt); 

	/* return the correction to moduli due to plasticity (if any)
	 *
	 * Note: Return mapping occurs during the call to StressCorrection.
	 *       The element passed in is already assumed to carry current
	 *       internal variable values */
	const dMatrixT& ModuliCorrection(const ElementCardT& element, int ip, double dt); 
	const dMatrixT& ModuliCorrectionEP(const ElementCardT& element, int ip);
	/* Modulus for checking perfectly plastic bifurcation */
	const dMatrixT& ModuliCorrPerfPlas(const ElementCardT& element, int ip);

	/* return a pointer to a new plastic element object constructed with
	 * the data from element */
	void AllocateElement(ElementCardT& element);

	enum InternalVariablesT {kkappa = 0,  // stress-like internal state variable
							kstressnorm = 1,  // norm of stress
                            kdgamma = 2,  // consistency parameter
                            kftrial = 3, // yield function value
			    			kdgamma2 = 4, // 2nd consistency par. at vertex
			    			kgamma = 5}; // accumulated plastic strain
	
	/** internal variables */
	dArrayT& Internal(void) { return fInternal; };
	
	/* element level data */
	void Update(ElementCardT& element, double dt);
	void Reset(ElementCardT& element);

	/* returns 1 if the trial elastic strain state lies outside of the 
	 * yield surface */
	int PlasticLoading(const dSymMatrixT& trialstrain, ElementCardT& element, int ip);

	/* computes the deviatoric stress corresponding to the given element
	 * and elastic strain.  The function returns a reference to the
	 * stress in fDevStress */
	dSymMatrixT& DeviatoricStress(const dSymMatrixT& trialstrain, 
		const ElementCardT& element);

	/* computes the hydrostatic (mean) stress. */
	double MeanStress(const dSymMatrixT& trialstrain, const ElementCardT& element);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/
	
private:

	/* load element data for the specified integration point */
	void LoadData(const ElementCardT& element, int ip);

protected:

	/* element level internal state variables */
	dSymMatrixT fPlasticStrain; //total plastic strain (deviatoric and volumetric)
	dSymMatrixT fUnitNorm;      //unit normal to the yield surface
	dArrayT     fInternal;      //internal variables

private:

	/* pointer to calling DPSSKStV */
	//DPSSKStVLoc *fKStV;

	/* number of integration points */
	int fNumIP;

	/* material parameters **/
	double fmu;
	double flambda;
	double fK;
	double fX_H;
	double fX;
	double fMeanStress;
  
	/* return values */
	dSymMatrixT	fElasticStrain;
	dSymMatrixT	fStressCorr;
	dMatrixT	fModuliCorr;
	dMatrixT    fModuliCorrPerfPlas;
  		
	/* work space */
	dSymMatrixT fDevStress;
	dSymMatrixT fDevStrain; /* deviatoric part of the strain tensor */
	dSymMatrixT IdentityTensor2;  

	dMatrixT      fTensorTemp;
	dSymMatrixT   One;  
};

} // namespace Tahoe 
#endif /* _DP_SS_LIN_HARD_LOC_T_H_ */
