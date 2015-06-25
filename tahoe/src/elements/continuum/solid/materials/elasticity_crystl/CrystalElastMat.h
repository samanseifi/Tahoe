/* $Id: CrystalElastMat.h,v 1.5 2011/12/01 21:11:38 bcyansfn Exp $ */
/*
  File: CrystalElastMat.h
*/

#ifndef _CRYSTAL_ELAST_MAT_H_
#define _CRYSTAL_ELAST_MAT_H_

#include "ios_fwd_decl.h"

#include <fstream>
#include "ArrayT.h"
#include "LocalArrayT.h"


namespace Tahoe {

class CrystalElast;
class ifstreamT;
class dArrayT;
class dMatrixT;

class CrystalElastMat
{
 public:
  // constructor/virtual destructor
  CrystalElastMat(CrystalElast& poly);
  virtual ~CrystalElastMat();

  // compute elastic material constants
  virtual void ElasticityProps(dArrayT& matprop, double Temp_DegC, int elem, int intpt);

  // compute thermal material properties
  virtual void ThermalProps(dMatrixT& alpha, double Temp_DegC);

  // query for isotropic/anisotropic elasticity (default: false)
  virtual bool IsIsotropic() const;

 protected:
  // general stiffness coefficients
  double fC11;
  double fC12;
  double fC44;    // fC44=0.5*(fC11-fC12) for isotropic elasticity

  // temperature dependent moduli
  void CalculateModuli(double DegC);
  double TempDepModuli(double Temp, double const1, double const2, double const3);

  // temperature dependent thermal coefficients
  void CalculateAlpha(dMatrixT& alpha, double DegC);
};

} // namespace Tahoe 
#endif /* _CRYSTAL_ELAST_MAT_H_ */
