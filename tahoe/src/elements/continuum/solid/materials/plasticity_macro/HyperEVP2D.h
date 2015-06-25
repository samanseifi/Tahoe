/* $Id: HyperEVP2D.h,v 1.8 2011/12/01 21:11:38 bcyansfn Exp $ */
#ifndef _HYPER_EVP_2D_H_
#define _HYPER_EVP_2D_H_

#include "HyperEVP3D.h"

#include <iostream>
#include "dMatrixT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

class ifstreamT;
class SolidElementT;

class HyperEVP2D : public HyperEVP3D
{
 public:
  // constructor
  HyperEVP2D(ifstreamT& in, const FSMatSupportT& support);

  // Cauchy stress
  virtual const dSymMatrixT& s_ij();   

  // tangent modulus
  virtual const dMatrixT& c_ijkl();

 protected:
  // Cauchy stress in 2D
  dSymMatrixT f2Ds_ij;

  // tangent moduli in 2D
  dMatrixT f2Dc_ijkl; 
};

} // namespace Tahoe 
#endif /* _HYPER_EVP_2D_ */
