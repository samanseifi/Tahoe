/* $Id: MooneyRivlin.cpp,v 1.7 2011/12/01 21:11:37 bcyansfn Exp $ */
/* created:   TDN (5/31/2001) */
/* Phi(I1,J) = mu/2*(I1-3)+kappa/4*(J^2-1-2*ln(J)) */
/* I1 = trace(C); J=sqrt(det(C)) */
#include "MooneyRivlin.h"
#include "ExceptionT.h"

#include <cmath>
#include <iostream>


using namespace Tahoe;
const double third = 1.0/3.0;

MooneyRivlin::MooneyRivlin(void):
	fc1(0.0),
	fc2(0.0)
{
	SetName("mooney-rivlin");
}

/* set parameters */
void MooneyRivlin::SetKappaMu(double kappa, double mu)
{
	fMu = (fc1-fc2);
	SetKappa(kappa);
}

void MooneyRivlin::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	PotentialT::DefineParameters(list);

	list.AddParameter(fc1, "c1");
	list.AddParameter(fc2, "c2");
	
	/* set the description */
	list.SetDescription("Psi=0.5*c1(I1_bar-3) - 0.5*c2(I2_bar -3)");	
	/*Note that the constants c1 and c2 correspond to the constants of a 2 term Ogden model. * 
	 *Mooney Rivlin is conventionally written as Psi= c1b (I1_bar-3) + c2b (I2_bar-3),		 *
	 *where c1b=c1/2 and c2b = -c2/2.
	 I've verified the stress formulation several times.  It's correct.  TDN 6/29/2011 */
}

void MooneyRivlin::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	PotentialT::TakeParameterList(list);

	fc1 = list.GetParameter("c1");
	fc2 = list.GetParameter("c2");

	/* check */
	if (fc1-fc2 < kSmall) ExceptionT::BadInputValue("MooneyRivlin::TakeParameterList",
		"expecting a nonzero shear modulus: mu= c1-c2: %d", fc1);
}


double MooneyRivlin::Energy(const dArrayT& lambda_bar, const double& J,  double temperature)
{
  double I1 = lambda_bar[0]+lambda_bar[1]+lambda_bar[2];
  double I2 = lambda_bar[0]*lambda_bar[1]+lambda_bar[1]*lambda_bar[2]+lambda_bar[0]*lambda_bar[2];

  double phi =0.5*fc1*(I1-3.0)-0.5*fc2*(I2-3.0);
  phi += MeanEnergy(J);
  return(phi);
}
void MooneyRivlin::DevStress(const dArrayT& lambda_bar,dArrayT& tau,  double temperature)
{
  int nsd = tau.Length();
  
  const double& l0 = lambda_bar[0];
  const double& l1 = lambda_bar[1];
  const double& l2 = lambda_bar[2];
  
	/*c1(I1bar - 3.0) eigen kirchhoff stress*/
	tau[0] = fc1*third*(2.0*l0 - l1 - l2);
	tau[1] = fc1*third*(2.0*l1 - l2 - l0);

	/*c2(I2bar - 3.0) eigen kirchhoff stress*/
	tau[0] += fc2*third*(2.0/l0 - 1.0/l1 - 1.0/l2);
	tau[1] += fc2*third*(2.0/l1 - 1.0/l2 - 1.0/l0);

	if (nsd == 3)
	{
		tau[2] = fc1*third*(2.0*l2 - l0 - l1);
		tau[2] += fc2*third*(2.0/l2 - 1.0/l0 - 1.0/l1);
	}
}

void MooneyRivlin::DevMod(const dArrayT& lambda_bar, dSymMatrixT& eigenmodulus, double temperature)
{
  int nsd = eigenmodulus.Rows();
  double ninth = third*third;
  
  const double& l0 = lambda_bar[0];
  const double& l1 = lambda_bar[1];
  const double& l2 = lambda_bar[2];
  
	eigenmodulus[0] = 2.0*fc1*ninth*(4.0*l0 + l1 + l2);
	eigenmodulus[1] = 2.0*fc1*ninth*(4.0*l1 + l2 + l0);
	
	eigenmodulus[0] += -2.0*fc2*ninth*(4.0/l0 + 1.0/l1 + 1.0/l2);
	eigenmodulus[1] += -2.0*fc2*ninth*(4.0/l1 + 1.0/l2 + 1.0/l0);

	if (nsd == 2)
	{
		eigenmodulus[2] = 2.0*fc1*ninth*(-2.0*l0 - 2.0*l1 + l2);
		eigenmodulus[2] += -2.0*fc2*ninth*(-2.0/l0 - 2.0/l1 + 1.0/l2);
	}
	else 
	{
		eigenmodulus[2] = 2.0*fc1*ninth*(4.0*l2 + l0 + l1);
		eigenmodulus[3] = 2.0*fc1*ninth*(-2.0*l1 - 2.0*l2 + l0);
		eigenmodulus[4] = 2.0*fc1*ninth*(-2.0*l0 - 2.0*l2 + l1);
		eigenmodulus[5] = 2.0*fc1*ninth*(-2.0*l0 - 2.0*l1 + l2);

		eigenmodulus[2] += -2.0*fc2*ninth*(4.0/l2 + 1.0/l0 + 1.0/l1);
		eigenmodulus[3] += -2.0*fc2*ninth*(-2.0/l1 - 2.0/l2 + 1.0/l0);
		eigenmodulus[4] += -2.0*fc2*ninth*(-2.0/l0 - 2.0/l2 + 1.0/l1);
		eigenmodulus[5] += -2.0*fc2*ninth*(-2.0/l0 - 2.0/l1 + 1.0/l2);
	}
}
