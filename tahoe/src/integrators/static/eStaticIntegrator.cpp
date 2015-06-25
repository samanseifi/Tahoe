/* $Id: eStaticIntegrator.cpp,v 1.6 2006/11/14 03:28:56 paklein Exp $ */
/* created: paklein (10/14/1996) */
#include "eStaticIntegrator.h"
#include "Environment.h"
#include "ExceptionT.h"

using namespace Tahoe;

/* constructor */
eStaticIntegrator::eStaticIntegrator(void):
	fLHSMode(kNormal)
{ 

}

/* returns 1 if the algorithm requires M, C, or K and sets const equal
* to the coefficient for the linear combination of components in the
* element effective mass matrix */
int eStaticIntegrator::FormM(double& constM) const
{
	if (fLHSMode == kFormMOnly) {
		constM = 1;
		return 1;
	} else {
		return 0;
	}
}

int eStaticIntegrator::FormC(double& constC) const
{
#pragma unused(constC)
	return 0;
}

int eStaticIntegrator::FormK(double& constK) const
{
	if (fLHSMode == kFormMOnly) {
		return 0;
	} else {
		constK = 1.0;
		return 1;
	}
}

/* components of the internal force vector */
int eStaticIntegrator::FormMa(double& constMa) const
{
#pragma unused(constMa)
	return 0;
}

int eStaticIntegrator::FormCv(double& constCv) const
{
#pragma unused(constCv)
	return 0;
}

int eStaticIntegrator::FormKd(double& constKd) const
{
	constKd = 1.0;
	return 1;
}

/*************************************************************************
* Protected
*************************************************************************/

/* recalculate constants */
void eStaticIntegrator::eComputeParameters(void) { }
