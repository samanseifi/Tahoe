/* $Id: GradCrystalPlast2D.h,v 1.8 2005/01/21 16:51:21 paklein Exp $ */
#ifndef _GRAD_CRYSTAL_PLAST_2D_H_
#define _GRAD_CRYSTAL_PLAST_2D_H_

#include "GradCrystalPlast.h"

#include "ArrayT.h"
#include "dArray2DT.h"
#include "LocalArrayT.h"

namespace Tahoe {

class GradCrystalPlast2D: public GradCrystalPlast
{
 public:
  // constructor
  GradCrystalPlast2D(void);

  // crystal Cauchy stress
  virtual const dSymMatrixT& s_ij();

  // crystal modulus 
  virtual const dMatrixT& c_ijkl();

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

 protected: 

  // crystal Cauchy stress in 2D
  dSymMatrixT f2Ds_ij;
  
  // crystal tangent moduli in 2D
  dMatrixT f2Dc_ijkl;
};

} // namespace Tahoe 
#endif /* _GRAD_CRYSTAL_PLAST_2D_H_ */
