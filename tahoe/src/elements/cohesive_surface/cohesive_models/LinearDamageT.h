/* $Id: LinearDamageT.h,v 1.14 2006/06/03 16:25:14 tdnguye Exp $ */
/* created: paklein (08/26/2000) */
#ifndef _LINEAR_DAMAGE_T_H_
#define _LINEAR_DAMAGE_T_H_

/* base class */
#include "SurfacePotentialT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;

/** linear cohesive law. Tractions evolve with a linear 
 * damage-like decay with opening displacement. */
class LinearDamageT: public SurfacePotentialT
{
public:

	/** constructor.
	 * \param in input stream to read parameters
	 * \param init_traction location of traction on surface
	 *        at initialization */
	LinearDamageT(ifstreamT& in, const dArrayT& init_traction);

	/** return the number of state variables needed by the model */
	virtual int NumStateVariables(void) const;

	/** initialize the state variable array */
	virtual void InitStateVariables(ArrayT<double>& state);

	/** dissipated energy */
	virtual double FractureEnergy(const ArrayT<double>& state);

	/** potential energy */
	virtual double Potential(const dArrayT& jump, const ArrayT<double>& state);
	
	/** surface traction. Internal variables are integrated over the current
	 * time step. */	
	virtual const dArrayT& Traction(const dArrayT& jump, ArrayT<double>& state, const dArrayT& sigma, bool qIntegrate);

	/** tangent stiffness */
	virtual const dMatrixT& Stiffness(const dArrayT& jump, const ArrayT<double>& state, const dArrayT& sigma);

	/** surface status */
	virtual StatusT Status(const dArrayT& jump, const ArrayT<double>& state);
	
private:

	/** traction at initialization */
	const dArrayT& fInitTraction;

	/* traction potential parameters */
	double fd_c_n; /**< characteristic normal opening to failure */
	double fd_c_t; /**< characteristic tangential opening to failure */
	
	/* penetration stiffness */
	double fpenalty; /**< stiffening multiplier during interpenetration */
	double fK;       /**< calculated penetration stiffness */
	
};

} // namespace Tahoe 
#endif /* _LINEAR_DAMAGE_T_H_ */
