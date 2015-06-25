/* $Id: ScaledCsch.cpp,v 1.4 2011/12/01 20:37:57 beichuan Exp $ */

#include "ScaledCsch.h"
#include <iostream>
#include <cmath>
#include "ExceptionT.h"
#include "dArrayT.h"

static const double Pi = 4.0*atan(1.0);

/* constructors */

using namespace Tahoe;

ScaledCsch::ScaledCsch(double A, double B):
	feta0(A),
	ftau0(B)
{ 
	SetName("scaled-csch");

}

ScaledCsch::ScaledCsch(void):feta0(0.0), ftau0(0.0) 
{ 
	SetName("scaled-csch");
}

/* I/O */
void ScaledCsch::Print(ostream& out) const
{

	/* parameters */
	out << " initial viscosity. . . . . . . . . . . . . . . . . . . . . = " << feta0 << '\n';
	out << " activation stress. . . . . . . . . . . . . . . . . . . . . . = " << ftau0 << '\n';
}

void ScaledCsch::PrintName(ostream& out) const
{
	out << "    ieta tau/tau0 Csch(tau/tau0) \n";
}

/*
* Returning values
*/
double ScaledCsch::Function(double tau) const
{
	double x = tau/ftau0;
	
	if (x > kSmall)
		return (feta0 * x/sinh(x));
	else return(feta0);
}

double ScaledCsch::DFunction(double tau) const
{
	double x = tau/ftau0;
	
	if (x > kSmall)
		return (feta0/ftau0 * (1.0 - x*cosh(x)/sinh(x))/sinh(x) );
	else return(0.0);
}

double ScaledCsch::DDFunction(double tau) const
{
	double x = tau/ftau0;
	
	if (x > kSmall)
	{
		double A = 3.0*x + x*cosh(2.0*x) - 2.0*sinh(2.0*x);
		double B = sinh(x);
		return (0.5*feta0/(ftau0*ftau0) * A/(B*B*B));
	}
	else return(-feta0/(3.0*ftau0*ftau0));
}

/* returning values in groups */
dArrayT& ScaledCsch::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pr = in.Pointer();
	double* pU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double tau = (*pr++);
		double x = tau/ftau0;
	
		if (x > kSmall)
			*pU++ = feta0 * x/sinh(x);
		else *pU++ = feta0;
	}
	return(out);
}

dArrayT& ScaledCsch::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pr = in.Pointer();
	double* pdU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double tau = (*pr++);
		double x = tau/ftau0;
	
		if (x > kSmall)
			*pdU++ = feta0/ftau0 * (1.0/sinh(x)-x/cosh(x));
		else *pdU++ = 0.0;
	}
	return(out);
}

dArrayT& ScaledCsch::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pr = in.Pointer();
	double* pddU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double tau = (*pr++);
		double x = tau/ftau0;

		if (x > kSmall)
		{
			double A = 3.0*x + x*cosh(2.0*x) - 2.0*sinh(2.0*x);
			double B = sinh(x);
			*pddU++ = 0.5*feta0/(ftau0*ftau0) * A/(B*B*B);
		}
		else *pddU++ = -feta0/(3.0*ftau0*ftau0);
	}
	return(out);
}

void ScaledCsch::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	C1FunctionT::DefineParameters(list);

	ParameterT eta0(ParameterT::Double, "viscosity");
	ParameterT tau0(ParameterT::Double, "activation_stress");
	
	eta0.AddLimit(0.0,LimitT::Lower);
	tau0.AddLimit(0.0,LimitT::Lower);
	list.AddParameter(eta0);
	list.AddParameter(tau0);

	/* set the description */
	list.SetDescription("f(tau) = eta0* tau/tau0 / Sinh(tau/tau0)");	
}

void ScaledCsch::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	C1FunctionT::TakeParameterList(list);

	/* bound */
	LimitT lower(0.0, LimitT::Lower);

	feta0 = list.GetParameter("viscosity");
	ftau0 = list.GetParameter("activation_stress");

}

