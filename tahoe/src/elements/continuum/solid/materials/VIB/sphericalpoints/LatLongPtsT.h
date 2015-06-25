/* $Id: LatLongPtsT.h,v 1.5 2009/05/21 22:30:27 tdnguye Exp $ */
/* created: paklein (10/31/1997) */
#ifndef _LATLONG_PTS_T_H_
#define _LATLONG_PTS_T_H_

/* base class */
#include "SpherePointsT.h"

namespace Tahoe {

/* forward declarations */
//class C1FunctionT;
class ifstreamT;

class LatLongPtsT: public SpherePointsT
{
public:

	/** constructor */
	LatLongPtsT(int n_phi, int n_theta);

	/** generate sphere points:
	 *
	 *   phi   = angle about z from x
	 *   theta = angle about x from z
	 *
	 * The final orientation is generated by applied the
	 * phi and theta rotations in succession about the local
	 * axes.
	 */
	virtual const dArray2DT& SpherePoints(double phi, double theta);

//	/** list of jacobian determinants */
//	virtual const dArrayT& Jacobians(const double theta, const double phi, const C1FunctionT* func_theta, const C1FunctionT* func_phi);
	
private:

	/* parameters */
	int	fNphi;
	int fNtheta;
};

} // namespace Tahoe 
#endif /* _LATLONG_PTS_T_H_ */
