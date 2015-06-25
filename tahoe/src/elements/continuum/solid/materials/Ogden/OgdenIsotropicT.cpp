/* $Id: OgdenIsotropicT.cpp,v 1.15 2011/12/01 21:11:37 bcyansfn Exp $ */
/* created: paklein (10/01/2000) */
#include "OgdenIsotropicT.h"

#include <cmath>

using namespace Tahoe;

/* constructor */
OgdenIsotropicT::OgdenIsotropicT(void):
	ParameterInterfaceT("Ogden_isotropic"),
	fSpectralDecomp(NULL)
{

}

OgdenIsotropicT::~OgdenIsotropicT(void)
{
	delete fSpectralDecomp;
}

/* modulus */
const dMatrixT& OgdenIsotropicT::c_ijkl(void)
{
	/* deformation gradient */
	const dMatrixT& Fmat = F();
	
	/* transform */
	fModulus.SetToScaled(1.0/Fmat.Det(), PushForward(Fmat, OgdenIsotropicT::C_IJKL()));
	return fModulus;
} 
/**< \todo compute directly in spatial representation rather than transforming */
	
/* stresses */
const dSymMatrixT& OgdenIsotropicT::s_ij(void)
{
	/* deformation gradient */
	const dMatrixT& Fmat = F();
	
	/* transform */
	fStress.SetToScaled(1.0/Fmat.Det(), PushForward(Fmat, OgdenIsotropicT::S_IJ()));
	return fStress;
}
/**< \todo compute directly in spatial representation rather than transforming */

/* material description */
const dMatrixT& OgdenIsotropicT::C_IJKL(void)
{
	/* stretch */
	Compute_C(fC);

	/* spectral decomposition */
	fSpectralDecomp->SpectralDecomp_Jacobi(fC, false);

	/* principal values */
	const dArrayT& eigenstretch = fSpectralDecomp->Eigenvalues();
	ddWddE(eigenstretch, fdWdE, fddWddE);

	/* axial */
	fModulus = fSpectralDecomp->EigsToRank4(fddWddE);

	/* shear terms */
	const ArrayT<dArrayT>& eigenvectors = fSpectralDecomp->Eigenvectors();
	double diff_stretch, dtde;
	if (NumSD() == 2)
	{
		diff_stretch = eigenstretch[0] - eigenstretch[1];

		/* modulus coefficient */
		if (fabs(diff_stretch) > kSmall)
			dtde = (fdWdE[0] - fdWdE[1])/diff_stretch;
		else
			dtde = fddWddE(0,0) - fddWddE(0,1);
		MixedRank4_2D(eigenvectors[0], eigenvectors[1], fModMat);
		fModulus.AddScaled(2.0*dtde, fModMat);
	}
	else
	{
		/* 1,2 */
		diff_stretch = eigenstretch[0] - eigenstretch[1];
		if (fabs(diff_stretch) > kSmall)
			dtde = (fdWdE[0] - fdWdE[1])/diff_stretch;
		else
			dtde = fddWddE(0,0) - fddWddE(0,1);
		MixedRank4_3D(eigenvectors[0], eigenvectors[1], fModMat);
		fModulus.AddScaled(2.0*dtde, fModMat);
	
		/* 1,3 */
		diff_stretch = eigenstretch[0] - eigenstretch[2];
		if (fabs(diff_stretch) > kSmall)
			dtde = (fdWdE[0] - fdWdE[2])/diff_stretch;
		else
			dtde = fddWddE(0,0) - fddWddE(0,2);
		MixedRank4_3D(eigenvectors[0], eigenvectors[2], fModMat);
		fModulus.AddScaled(2.0*dtde, fModMat);
	
		/* 2,3 */
		diff_stretch = eigenstretch[1] - eigenstretch[2];
		if (fabs(diff_stretch) > kSmall)
			dtde = (fdWdE[1] - fdWdE[2])/diff_stretch;
		else
			dtde = fddWddE(1,1) - fddWddE(1,2);
		MixedRank4_3D(eigenvectors[1], eigenvectors[2], fModMat);
		fModulus.AddScaled(2.0*dtde, fModMat);
	}
	return fModulus;
}

const dSymMatrixT& OgdenIsotropicT::S_IJ(void)
{
	/* stretch */
	Compute_C(fC);

	/* spectral decomposition */
	fSpectralDecomp->SpectralDecomp_Jacobi(fC, false);
	//fSpectralDecomp->SpectralDecomp(C(), false); // closed-form decomposition
	const ArrayT<dArrayT>& eigenvectors = fSpectralDecomp->Eigenvectors();

	/* principal values */
	dWdE(fSpectralDecomp->Eigenvalues(), fdWdE);
	fStress = fSpectralDecomp->EigsToRank2(fdWdE);

	return (fStress);
}

/* return the pressure associated with the last call to stress */
double OgdenIsotropicT::Pressure(void) const
{
	return dArrayT::Dot(fSpectralDecomp->Eigenvalues(), fdWdE)/
			sqrt(fSpectralDecomp->Eigenvalues().Product());
}

/* accept parameter list */
void OgdenIsotropicT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	FSIsotropicMatT::TakeParameterList(list);

	/* construct spectral decomposition */
	int nsd = NumSD();
	fSpectralDecomp = new SpectralDecompT(nsd);

	/* dimension work space */
	fC.Dimension(nsd);
	fEigs.Dimension(nsd);
	fdWdE.Dimension(nsd);
	fddWddE.Dimension(nsd);
	fModMat.Dimension(dSymMatrixT::NumValues(nsd));
	fModulus.Dimension(dSymMatrixT::NumValues(nsd));
	fStress.Dimension(nsd);
}

/*************************************************************************
 * Private
 *************************************************************************/

/* construct symmetric rank-4 mixed-direction tensor (6.1.44) */
void OgdenIsotropicT::MixedRank4_2D(const dArrayT& a, const dArrayT& b,
	dMatrixT& rank4_ab) const
{
#if __option(extended_errorcheck)
	if (a.Length() != 2 ||
	    b.Length() != 2 ||
	    rank4_ab.Rows() != 3 ||
	    rank4_ab.Cols() != 3) throw ExceptionT::kSizeMismatch;
#endif

	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11;

//	z1 = A(1.);
//	z2 = A(2.);
//	z3 = B(1.);
//	z4 = B(2.);

	z1 = a[0];
	z2 = a[1];
	z3 = b[0];
	z4 = b[1];

	z5 = z1*z1;
	z6 = z2*z2;
	z7 = z3*z3;
	z8 = 2.*z1*z2*z3*z4;
	z9 = z4*z4;
	z3 = 2.*z3*z4;
	z4 = z3*z5;
	z3 = z3*z6;
	z10 = 2.*z1*z2*z7;
	z11 = 2.*z5*z7;
	z7 = z6*z7;
	z1 = 2.*z1*z2*z9;
	z2 = z5*z9;
	z5 = 2.*z6*z9;
	z4 = z10 + z4;
	z1 = z1 + z3;
	z2 = z2 + z7 + z8;
	z3 = 0.5*z4;
	z1 = 0.5*z1;
	z2 = 0.5*z2;

	//{{z11, z8, z3},
	// {z8, z5, z1},
	// {z3, z1, z2}}

	double* p = rank4_ab.Pointer();
	*p++ = z11;
	*p++ = z8;
	*p++ = z3;
	*p++ = z8;
	*p++ = z5;
	*p++ = z1;
	*p++ = z3;
	*p++ = z1;
	*p   = z2;
}

void OgdenIsotropicT::MixedRank4_3D(const dArrayT& a, const dArrayT& b,
	dMatrixT& rank4_ab) const
{
#if __option(extended_errorcheck)
	if (a.Length() != 3 ||
	    b.Length() != 3 ||
	    rank4_ab.Rows() != 6 ||
	    rank4_ab.Cols() != 6) throw ExceptionT::kSizeMismatch;
#endif

	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
	double z13, z14, z15, z16, z17, z18, z19, z20, z21, z22, z23, z24;
	double z25, z26, z27, z28, z29, z30, z31, z32, z33, z34, z35, z36;
	double z37, z38, z39, z40, z41;

//	z1 = A(1.);
//	z2 = A(2.);
//	z3 = A(3.);
//	z4 = B(1.);
//	z5 = B(2.);
//	z6 = B(3.);

	z1 = a[0];	
	z2 = a[1];
	z3 = a[2];
	z4 = b[0];
	z5 = b[1];
	z6 = b[2];
	z7 = z1*z1;
	z8 = z2*z2;
	z9 = z3*z3;
	z10 = 2.*z1*z4;
	z11 = z10*z2;
	z12 = z2*z4;
	z13 = z1*z12;
	z14 = z4*z4;
	z15 = z10*z5;
	z15 = z15*z3;
	z16 = z11*z5;
	z17 = z12*z5;
	z18 = 2.*z17;
	z17 = z17*z3;
	z18 = z18*z3;
	z19 = z1*z4*z5;
	z19 = z19*z3;
	z20 = z5*z5;
	z11 = z11*z6;
	z13 = z13*z6;
	z10 = z10*z3*z6;
	z12 = z12*z3*z6;
	z21 = 2.*z12;
	z22 = z1*z2*z5*z6;
	z23 = 2.*z22;
	z24 = z1*z3*z5*z6;
	z25 = 2.*z24;
	z26 = 2.*z2*z3*z5*z6;
	z27 = z6*z6;
	z28 = 2.*z1*z14;
	z29 = 2.*z1*z2;
	z30 = z14*z2;
	z11 = z11 + z15;
	z15 = z1*z20*z3;
	z31 = 2.*z2*z20*z3;
	z18 = z18 + z23;
	z21 = z21 + z25;
	z23 = z1*z2*z27;
	z1 = 2.*z1*z27*z3;
	z25 = 2.*z2*z27*z3;
	z2 = z2*z28;
	z28 = z28*z3;
	z29 = z20*z29;
	z3 = z3*z30;
	z30 = 2.*z14*z7;
	z32 = z20*z7;
	z33 = z27*z7;
	z34 = 2.*z4*z5*z7;
	z35 = 2.*z4*z6*z7;
	z7 = z5*z6*z7;
	z36 = z14*z8;
	z37 = 2.*z20*z8;
	z38 = z27*z8;
	z39 = 2.*z4*z5*z8;
	z40 = z4*z6*z8;
	z8 = 2.*z5*z6*z8;
	z14 = z14*z9;
	z20 = z20*z9;
	z27 = 2.*z27*z9;
	z41 = z4*z5*z9;
	z4 = 2.*z4*z6*z9;
	z5 = 2.*z5*z6*z9;
	z6 = 0.5*z11;
	z9 = 0.5*z18;
	z11 = 0.5*z21;
	z2 = z2 + z34;
	z18 = z28 + z35;
	z3 = z13 + z19 + z3 + z7;
	z7 = z16 + z32 + z36;
	z13 = z29 + z39;
	z15 = z15 + z17 + z22 + z40;
	z8 = z31 + z8;
	z14 = z10 + z14 + z33;
	z17 = z20 + z26 + z38;
	z12 = z12 + z23 + z24 + z41;
	z1 = z1 + z4;
	z4 = z25 + z5;
	z2 = 0.5*z2;
	z5 = 0.5*z18;
	z3 = 0.5*z3;
	z7 = 0.5*z7;
	z13 = 0.5*z13;
	z15 = 0.5*z15;
	z8 = 0.5*z8;
	z14 = 0.5*z14;
	z17 = 0.5*z17;
	z12 = 0.5*z12;
	z1 = 0.5*z1;
	z4 = 0.5*z4;
	
	//{{z30, z16, z10,  z6,  z5,  z2},
	// {z16, z37, z26,  z8,  z9, z13},
	// {z10, z26, z27,  z4,  z1, z11},
	// { z6,  z8,  z4, z17, z12, z15},
	// { z5,  z9,  z1, z12, z14,  z3},
	// { z2, z13, z11, z15,  z3,  z7}}
	
	double* p = rank4_ab.Pointer();
	*p++ = z30;
	*p++ = z16;
	*p++ = z10;
	*p++ = z6;
	*p++ = z5;
	*p++ = z2;
	*p++ = z16;
	*p++ = z37;
	*p++ = z26;
	*p++ = z8;
	*p++ = z9;
	*p++ = z13;
	*p++ = z10;
	*p++ = z26;
	*p++ = z27;
	*p++ = z4;
	*p++ = z1;
	*p++ = z11;
	*p++ = z6;
	*p++ = z8;
	*p++ = z4;
	*p++ = z17;
	*p++ = z12;
	*p++ = z15;
	*p++ = z5;
	*p++ = z9;
	*p++ = z1;
	*p++ = z12;
	*p++ = z14;
	*p++ = z3;
	*p++ = z2;
	*p++ = z13;
	*p++ = z11;
	*p++ = z15;
	*p++ = z3;
	*p  = z7;
}
