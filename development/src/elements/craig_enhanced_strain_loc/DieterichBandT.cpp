
#include "DieterichBandT.h"

using namespace Tahoe;

/* constructor */
DieterichBandT::DieterichBandT(dArrayT normal, dArrayT slipDir, dArrayT
perpSlipDir, dArrayT coords,
				   double fH_delta_0, double
				   residCohesion, ArrayT<dSymMatrixT> stressList,
				   SSEnhLocDieterichT* element, double
			       theta_0):
BandT(normal, slipDir, perpSlipDir, coords,
				   fH_delta_0, residCohesion, stressList,
      element),
fTheta(theta_0),
fLastTheta(theta_0),
fSlipRate(0.0),
fLastSlipRate(0.0),
fLastThetaRate(1.0) //what should this be?
//fLastJumpIncrement(0.0),
//fDeltaTheta(0.0)
{
	//currentElement = dynamic_cast<SSEnhLocDieterichT*> (currentElement); 
}


double DieterichBandT::Theta()
{
  return fTheta;
}

/*
double DieterichBandT::DeltaTheta()
{
  return fDeltaTheta;
}
*/

void DieterichBandT::StoreTheta(double theta)
{
  fTheta = theta;
}

double DieterichBandT::SlipRateLast()
{
  return fLastSlipRate;
}

void DieterichBandT::CloseStep()
{
  /* no longer inherited - with new cohesion softening model, 
    only track initial residual cohesion*/
  //BandT::CloseStep();
    IncrementJump();

  /* update ISV theta */
  fLastTheta = fTheta;
  /* store jump increment for next time step */
  fLastSlipRate = fSlipRate;
  //fLastThetaRate = 1.0 -  fTheta * fSlipRate / 1.0; // fix this(currentElement -> fD_c);
}

void DieterichBandT::UpdateThetaRate(double d_c)
{
  fLastThetaRate = 1.0 -  fTheta * fSlipRate / d_c;
}
