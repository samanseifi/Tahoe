/* $Id: PolyBasis2DT.cpp,v 1.8 2005/04/22 00:33:40 paklein Exp $ */
/* created: paklein (12/13/1999) */
#include "PolyBasis2DT.h"

using namespace Tahoe;

/* constructor */
PolyBasis2DT::PolyBasis2DT(int complete, bool cross_terms):
	BasisT(complete, 2),
	fCrossTerms(cross_terms)
{
	/* check */
	if (fComplete < 0 || fComplete > 3)
		ExceptionT::OutOfRange("PolyBasis2DT::PolyBasis2DT",
			"completeness must be [0,3]: %d", complete);
	if (fComplete >= 2 && fCrossTerms)
		ExceptionT::GeneralFail("PolyBasis2DT::PolyBasis2DT",
			"cross_terms not supported for completeness >= 2");
}
	
/* return the number of basis functions */
int PolyBasis2DT::BasisDimension(void) const
{
	switch (fComplete)
	{
		case 0:			
			return 1;
		case 1:
			if (!fCrossTerms)
				return 3;
			else
				return 4;
		case 2:
			if (!fCrossTerms)
				return 6;
			else
				return 9;
		case 3:
			if (!fCrossTerms)
				return 10;
			else
				return 16;
		default:
			ExceptionT::OutOfRange("PolyBasis2DT::BasisDimension");
	}
	
	return 0;
}

/* evaluate basis functions at coords */
void PolyBasis2DT::SetBasis(const dArray2DT& coords, int order)
{
#if __option(extended_errorcheck)
	/* dimension checking */
	const char caller[] = "PolyBasis2DT::SetBasis";
	if (coords.MinorDim() != fNumSD) ExceptionT::GeneralFail(caller);
	if (order > 3) ExceptionT::OutOfRange(caller); //kyonten (order increase to 3)
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
				fDP[1] = 0.0;
				if (order > 1)
				{
					fDDP[0] = 0.0;
					fDDP[1] = 0.0;
					fDDP[2] = 0.0;
					if (order > 2) // kyonten
					{
						fDDDP[0] = 0.0;
						fDDDP[1] = 0.0;
						fDDDP[2] = 0.0;
						fDDDP[3] = 0.0;
					}
				}
			}
			break;
		}
		case 1: // linear basis
		{
			const double* px = coords.Pointer();
			double*  pP0 = fP(0);
			double*  pP1 = fP(1);
			double*  pP2 = fP(2);
			double*  pP3 = NULL;
			
			double* pD0P0 = (fDP[0])(0);
			double* pD0P1 = (fDP[0])(1);
			double* pD0P2 = (fDP[0])(2);
			double* pD0P3 = NULL;

			double* pD1P0 = (fDP[1])(0);
			double* pD1P1 = (fDP[1])(1);
			double* pD1P2 = (fDP[1])(2);
			double* pD1P3 = NULL;

			if (fCrossTerms) {
				pP3   = (fCrossTerms) ?       fP(3) : NULL;
				pD0P3 = (fCrossTerms) ? (fDP[0])(3) : NULL;
				pD1P3 = (fCrossTerms) ? (fDP[1])(3) : NULL;
			}

			for (int i = 0; i < nnd; i++)
			{
				double x = *px++;
				double y = *px++;
			
				*pP0++ = 1.0;
				*pP1++ = x;
				*pP2++ = y;
				if (fCrossTerms) *pP3++ = x*y;
				
				if (order > 0)
				{
					*pD0P0++ = 0.0;
					*pD0P1++ = 1.0;
					*pD0P2++ = 0.0;

					*pD1P0++ = 0.0;
					*pD1P1++ = 0.0;
					*pD1P2++ = 1.0;
					
					if (fCrossTerms) {
						*pD0P3++ = y;
						*pD1P3++ = x;
					}

					if (order > 1)
					{
						fDDP[0] = 0.0;
						fDDP[1] = 0.0;
						fDDP[2] = 0.0;
						
						if (fCrossTerms)
							(fDDP[2])(3,i) = 1.0;
						
						if (order > 2)
						{
							fDDDP[0] = 0.0;
							fDDDP[1] = 0.0;
							fDDDP[2] = 0.0;
							fDDDP[3] = 0.0;
						}
					}
				}
			}
			break;
		}
		case 2: // quadratic basis: { 1, x, y, x^2, xy, y^2 } (no extra cross terms)
		{
			const double* px = coords.Pointer();
			for (int i = 0; i < nnd; i++)
			{
				double x = *px++;
				double y = *px++;

				fP(0,i) = 1.0;
				fP(1,i) = x;
				fP(2,i) = y;
				fP(3,i) = x*x;
				fP(4,i) = x*y;
				fP(5,i) = y*y;

				if (order > 0)
				{
					/* first derivatives: d/dx then d/dy */
					(fDP[0])(0,i) = 0.0; (fDP[0])(1,i) = 1.0; (fDP[0])(2,i) = 0.0;
					(fDP[0])(3,i) = 2.0*x; (fDP[0])(4,i) = y; (fDP[0])(5,i) = 0.0;

					(fDP[1])(0,i) = 0.0; (fDP[1])(1,i) = 0.0; (fDP[1])(2,i) = 1.0;
					(fDP[1])(3,i) = 0.0; (fDP[1])(4,i) = x; (fDP[1])(5,i) = 2.0*y;

					if (order > 1)
					{
						/* second derivatives: components [xx, yy, xy] */
						(fDDP[0])(0,i) = 0.0; (fDDP[0])(1,i) = 0.0; (fDDP[0])(2,i) = 0.0;
						(fDDP[0])(3,i) = 2.0; (fDDP[0])(4,i) = 0.0; (fDDP[0])(5,i) = 0.0;

						(fDDP[1])(0,i) = 0.0; (fDDP[1])(1,i) = 0.0; (fDDP[1])(2,i) = 0.0;
						(fDDP[1])(3,i) = 0.0; (fDDP[1])(4,i) = 0.0; (fDDP[1])(5,i) = 2.0;

						(fDDP[2])(0,i) = 0.0; (fDDP[2])(1,i) = 0.0; (fDDP[2])(2,i) = 0.0;
						(fDDP[2])(3,i) = 0.0; (fDDP[2])(4,i) = 1.0; (fDDP[2])(5,i) = 0.0;

						if (order > 2) /* quadratic -> all third derivatives vanish */
							for (int b = 0; b < 6; b++) {
								(fDDDP[0])(b,i) = 0.0; (fDDDP[1])(b,i) = 0.0;
								(fDDDP[2])(b,i) = 0.0; (fDDDP[3])(b,i) = 0.0;
							}
					}
				}
			}
			break;
		}
		case 3: // cubic basis: { 1, x, y, x^2, xy, y^2, x^3, x^2 y, x y^2, y^3 } (no extra cross)
		{
			const double* px = coords.Pointer();
			for (int i = 0; i < nnd; i++)
			{
				double x = *px++;
				double y = *px++;

				fP(0,i)=1.0; fP(1,i)=x; fP(2,i)=y;
				fP(3,i)=x*x; fP(4,i)=x*y; fP(5,i)=y*y;
				fP(6,i)=x*x*x; fP(7,i)=x*x*y; fP(8,i)=x*y*y; fP(9,i)=y*y*y;

				if (order > 0)
				{
					/* first derivatives */
					(fDP[0])(0,i)=0.0;(fDP[0])(1,i)=1.0;(fDP[0])(2,i)=0.0;
					(fDP[0])(3,i)=2.0*x;(fDP[0])(4,i)=y;(fDP[0])(5,i)=0.0;
					(fDP[0])(6,i)=3.0*x*x;(fDP[0])(7,i)=2.0*x*y;(fDP[0])(8,i)=y*y;(fDP[0])(9,i)=0.0;

					(fDP[1])(0,i)=0.0;(fDP[1])(1,i)=0.0;(fDP[1])(2,i)=1.0;
					(fDP[1])(3,i)=0.0;(fDP[1])(4,i)=x;(fDP[1])(5,i)=2.0*y;
					(fDP[1])(6,i)=0.0;(fDP[1])(7,i)=x*x;(fDP[1])(8,i)=2.0*x*y;(fDP[1])(9,i)=3.0*y*y;

					if (order > 1)
					{
						/* second derivatives: components [xx, yy, xy] */
						(fDDP[0])(0,i)=0.0;(fDDP[0])(1,i)=0.0;(fDDP[0])(2,i)=0.0;
						(fDDP[0])(3,i)=2.0;(fDDP[0])(4,i)=0.0;(fDDP[0])(5,i)=0.0;
						(fDDP[0])(6,i)=6.0*x;(fDDP[0])(7,i)=2.0*y;(fDDP[0])(8,i)=0.0;(fDDP[0])(9,i)=0.0;

						(fDDP[1])(0,i)=0.0;(fDDP[1])(1,i)=0.0;(fDDP[1])(2,i)=0.0;
						(fDDP[1])(3,i)=0.0;(fDDP[1])(4,i)=0.0;(fDDP[1])(5,i)=2.0;
						(fDDP[1])(6,i)=0.0;(fDDP[1])(7,i)=0.0;(fDDP[1])(8,i)=2.0*x;(fDDP[1])(9,i)=6.0*y;

						(fDDP[2])(0,i)=0.0;(fDDP[2])(1,i)=0.0;(fDDP[2])(2,i)=0.0;
						(fDDP[2])(3,i)=0.0;(fDDP[2])(4,i)=1.0;(fDDP[2])(5,i)=0.0;
						(fDDP[2])(6,i)=0.0;(fDDP[2])(7,i)=2.0*x;(fDDP[2])(8,i)=2.0*y;(fDDP[2])(9,i)=0.0;

						if (order > 2)
						{
							/* third derivatives (2D comps via ExpandIndex3): 0:xxx 1:xyy 2:xxy 3:yyy */
							for (int b = 0; b < 10; b++) {
								(fDDDP[0])(b,i)=0.0; (fDDDP[1])(b,i)=0.0;
								(fDDDP[2])(b,i)=0.0; (fDDDP[3])(b,i)=0.0;
							}
							(fDDDP[0])(6,i)=6.0; /* d3(x^3)/dx3   */
							(fDDDP[1])(8,i)=2.0; /* d3(x y^2)/dx dy2 */
							(fDDDP[2])(7,i)=2.0; /* d3(x^2 y)/dx2 dy */
							(fDDDP[3])(9,i)=6.0; /* d3(y^3)/dy3   */
						}
					}
				}
			}
			break;
		}
	}
}
