/* $Id: SIMOD_2DT.h,v 1.3 2005/07/23 22:18:25 paklein Exp $ */
#ifndef _SIMOD_2D_T_H_
#define _SIMOD_2D_T_H_

/* enabled */
#ifdef __SIMOD__

/* base class */
#include "SurfacePotentialT.h"

/* forward declarations */
#include "SharedInterfaceModel.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;

/** 2D models in the SIMOD library */
class SIMOD_2DT: public SurfacePotentialT
{
public:

	/** constructor */
	SIMOD_2DT(void);
	
	/** destructor */
	virtual ~SIMOD_2DT(void);

	/** return the number of state variables needed by the model */
	int NumStateVariables(void) const;

	/** initialize the state variable array. By default, initialization
	 * involves only setting the array to zero. */
	virtual void InitStateVariables(ArrayT<double>& state);

	/** dissipated energy */
	virtual double FractureEnergy(const ArrayT<double>& state);

	/** potential energy */
	virtual double Potential(const dArrayT& jump_u, const ArrayT<double>& state);
	
	/** surface traction. Internal variables are integrated over the current
	 * time step. */	
	virtual const dArrayT& Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma, bool qIntegrate);

	/** tangent stiffness */
	virtual const dMatrixT& Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma);

	/** surface status */
	virtual StatusT Status(const dArrayT& jump_u, const ArrayT<double>& state);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;
	
	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:

	/** translate parameter types */
	static ValueT::TypeT SIMOD2Tahoe(parameter_item_spc::parameter_type simod_type);

private:

	/** SIMOD model */
	SharedInterfaceModel_spc::SharedInterfaceModel* fSIMOD;
};

} /* namespace Tahoe */

#endif /* __SIMOD__ */
#endif /* _SIMOD_2D_T_H_ */
