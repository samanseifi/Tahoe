// $Id: VoidGrowthModelImp.cpp,v 1.5 2011/12/01 21:11:38 bcyansfn Exp $
#include "VoidGrowthModelImp.h"

#include <iostream>


using namespace Tahoe;

// base class
VoidGrowthModelImp::VoidGrowthModelImp():
  fm (0.0)
{ }

VoidGrowthModelImp::~VoidGrowthModelImp() { }


// derived class: CocksVGModel
CocksVGModel::CocksVGModel(EVPFDBaseT& macro)
{ 
#pragma unused(macro)
}

CocksVGModel::~CocksVGModel() { }

void CocksVGModel::ACoefficients(const double& vvf, double& A1, double& A2) const
{
  if (vvf <= 1.5e-10) { A1 = 1.0; A2 = 0.0; return; }

  double mm = 2./(fm+1.);
  A1 = (1.+2./3.*vvf) * pow( (1.-vvf), -mm );
  A2 = 4.5/(1.+fm) * vvf/(1.+vvf) * pow( (1.-vvf), -mm );
}

void CocksVGModel::ADerivCoefficients(const double& vvf, double& dA1, double& dA2) const
{
  if (vvf <= 1.5e-10) { dA1 = 0.0; dA2 = 0.0; return; }

  double A1, A2;
  ACoefficients(vvf, A1, A2);
  dA1 = A1 * ( 1./(1.5+vvf) + 2./(1.+fm)/(1.-vvf) );
  dA2 = A2 * ( 1./vvf/(1.+vvf) + 2./(1.+fm)/(1.-vvf) );
}

//derived class: DuvaCrowVGModel
DuvaCrowVGModel::DuvaCrowVGModel(EVPFDBaseT& macro)
{ 
#pragma unused(macro)
}  

DuvaCrowVGModel::~DuvaCrowVGModel() { }

void DuvaCrowVGModel::ACoefficients(const double& vvf, double& A1, double& A2) const
{
  if (vvf <= 1.5e-10) { A1 = 1.0; A2 = 0.0; return; }

  double mm = 2./(fm+1.);
  A1 = (1.+2./3.*vvf) * pow( (1.-vvf), -mm );
  A2 = 2.25 * pow( 1./fm*(pow(vvf,-fm)-1.), -mm );
}

void DuvaCrowVGModel::ADerivCoefficients(const double& vvf, double& dA1, double& dA2) const
{
  if (vvf <= 1.5e-10) { dA1 = 0.0; dA2 = 0.0; return; }

  double A1, A2;
  ACoefficients(vvf, A1, A2);
  dA1 = A1 * ( 1./(1.5+vvf) + 2./(1.+fm)/(1.-vvf) );
  dA2 = A2 * 2.*fm/(1.+fm) / (vvf*(1.-pow(vvf,fm)));
}

//derived class: SofronisVGModel
SofronisVGModel::SofronisVGModel(EVPFDBaseT& macro)
{ 
#pragma unused(macro)
}  

SofronisVGModel::~SofronisVGModel() { }

void SofronisVGModel::ACoefficients(const double& vvf, double& A1, double& A2) const
{
  if (vvf <= 1.5e-10) { A1 = 1.0; A2 = 0.0; return; }

  double mm = 2./(fm+1.);
  A1 = pow( (1.-vvf)/(1.+vvf), -mm );
  A2 = 2.25 * pow( 1./fm*(pow(vvf,-fm)-1.), -mm );
}

void SofronisVGModel::ADerivCoefficients(const double& vvf, double& dA1, double& dA2) const
{
  if (vvf <= 1.5e-10) { dA1 = 0.0; dA2 = 0.0; return; }

  double A1, A2;
  ACoefficients(vvf, A1, A2);
  dA1 = A1 * 4./(1.+fm) / (1.-vvf*vvf);
  dA2 = A2 * 2.*fm/(1.+fm) / (vvf*(1.-pow(vvf,fm)));
}
