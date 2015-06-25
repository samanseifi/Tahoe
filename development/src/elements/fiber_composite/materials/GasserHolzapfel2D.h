 //$Id: GasserHolzapfel2D.h,v 1.1 2011/09/16 21:00:28 thao Exp $ 
// Model from Gasser TC, Ogden RW, Holzapfel GA. J R SocInterface 2006 3(6):15 – 35.

#ifndef _GASSER_HOLZAPFEL2D_H_
#define _GASSER_HOLZAPFEL2D_H_
 
/* base  class */
#include "SolidMaterialsConfig.h"
 
/* base class */
#include "FSFiberMatSplitT.h"

namespace Tahoe {

/** 2D Isotropic VIB solver using spectral decomposition formulation */
class GasserHolzapfel2D: public FSFiberMatSplitT
{
public:

	/* constructor */
	GasserHolzapfel2D(void);

	/* destructor */
	~GasserHolzapfel2D(void);
	
	/* strain energy density */
	virtual double StrainEnergyDensity(void);

	/*compute output variables*/
	virtual int NumOutputVariables() const;
	virtual void OutputLabels(ArrayT<StringT>& labels) const;
	virtual void ComputeOutput(dArrayT& output);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);

	/*@}*/

	/*calculates  matrix contribution to 2PK stress*/
	virtual void ComputeDevMatrixStress(const dSymMatrixT& Cbar,  dSymMatrixT& Stress);

	/*calculates matrix contribution to modulus*/
	virtual void ComputeDevMatrixMod(const dSymMatrixT& Cbar, dSymMatrixT& Stress, dMatrixT& Mod);
	
	/*calculates  matrix contribution to 2PK stress*/
	virtual double ComputeVolMatrixStress(const double I3);

	/*calculates matrix contribution to modulus*/
	virtual double ComputeVolMatrixMod(const double I3);

	/*computes integrated fiber stress in local frame*/
	virtual void ComputeFiberStress (const dSymMatrixT& Stretch, dSymMatrixT& Stress);
	
	/*computes integrated moduli in local frame*/
	virtual void ComputeFiberMod (const dSymMatrixT& Stretch, dSymMatrixT& Stress, dMatrixT& Mod);

protected:
	
	/*material parameters for matrix*/
	double fMu;  /*shear modulus*/
	double fBulkMod; /*bulk modulus*/
	double fb2;  /*bulk mod stiffening factor*/
	
	/*material parameters for fiber*/
	double fKappa; /*fiber dispersion parameters*/
	double fk1; /*fiber stiffness*/
	double fk2;/*fiber stiffening*/
	
};

}
#endif 
