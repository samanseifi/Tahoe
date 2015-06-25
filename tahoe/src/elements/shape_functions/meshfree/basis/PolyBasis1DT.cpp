/* $Id: PolyBasis1DT.cpp,v 1.8 2005/04/22 00:33:39 paklein Exp $ */
/* created: paklein (12/11/1999) */
#include "PolyBasis1DT.h"

using namespace Tahoe;

/* constructor */
PolyBasis1DT::PolyBasis1DT(int complete):
	BasisT(complete, 1)
{
	if (fComplete < 0 || fComplete > 1)
		ExceptionT::OutOfRange("PolyBasis1DT::PolyBasis1DT",
			"completeness must be [0,1]: %d", complete);
}
	
/* return the number of basis functions */
int PolyBasis1DT::BasisDimension(void) const
{
	return fComplete + 1;
}

/* evaluate basis functions at coords */
void PolyBasis1DT::SetBasis(const dArray2DT& coords, int order)
{
#if __option(extended_errorcheck)
	/* dimension checking */
	const char caller[] = "PolyBasis1DT::SetBasis";
	if (coords.MinorDim() != fNumSD) ExceptionT::GeneralFail(caller);
	if (order > 3) ExceptionT::OutOfRange(caller); //kyonten (order increased to 3)
#endif

	/* dimensions */
	int nnd = coords.MajorDim();

	/* dimension work space */
	Dimension(nnd);

	switch (fComplete)
	{
		case 0: // constant basis
		{
			fP = 1.0;
			if (order > 0)
			{
				fDP[0] = 0.0;
				if (order > 1)
				{
					fDDP[0] = 0.0;
					if (order > 2) // kyonten (third derivative)
					{
						fDDDP[0] = 0.0;
					}
				}
			}
			break;
		}
		case 1: // linear basis
		{
			const double* px = coords.Pointer();
			double* pP0  = fP(0);
			double* pP1  = fP(1);
			double* pDP0 = (fDP[0])(0);
			double* pDP1 = (fDP[0])(1);
			for (int i = 0; i < nnd; i++)
			{
				*pP0++ = 1.0;
				*pP1++ = *px++;
				if (order > 0)
				{
					*pDP0++ = 0.0;
					*pDP1++ = 1.0;
				}
			}
			
			if (order > 1) fDDP[0] = 0.0;
			if (order > 2) fDDDP[0] = 0.0; // kyonten
			break;
		}
	}
}
