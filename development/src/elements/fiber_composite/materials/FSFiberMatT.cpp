/* $Id: FSFiberMatT.cpp,v 1.12 2013/02/01 19:16:06 tahoe.kziegler Exp $ */
/* created: paklein (06/09/1997) */
#include "FSFiberMatT.h"
#include "FSFiberMatSupportT.h"
#include "iArray2DT.h"

using namespace Tahoe;

/* constructor */
FSFiberMatT::FSFiberMatT(void):
	ParameterInterfaceT("fiber_composite_material"),
	fFSFiberMatSupport(NULL),
	fNumSD(3),
	fNumStress(3),
	fNumModuli(6)
{

	
}

/* set the material support or pass NULL to clear */
void FSFiberMatT::SetFSFiberMatSupport(const FSFiberMatSupportT* support)
{
	/* set inherited material support */
	FSSolidMatT::SetFSMatSupport(support);

	fFSFiberMatSupport = support;

}

/* modulus */
const dMatrixT& FSFiberMatT::C_IJKL(void)
{	
	/* stretch */
	Compute_C(fC);
	
	/*calculate matrix contribution*/
	ComputeMatrixMod(fC, fStress, fModulus);
	
	/*fiber contribution*/
	ComputeFiberStretch(fC, fFiberStretch);
	ComputeFiberMod(fFiberStretch, fFiberStress, fFiberMod);
				
	/* rotate stress to lab coordinates */
	AssembleFiberStress(fFiberStress, fStress);

	/* rotate modulus to lab coordinates */
	AssembleFiberModuli(fFiberMod, fModulus);

	return fModulus;
}
	
/* stress */
const dSymMatrixT& FSFiberMatT::S_IJ(void)
{
	/* stretch */
	Compute_C(fC);
	
	/*matrix contribution*/
	/*calculate matrix contribution*/
	ComputeMatrixStress(fC, fStress);

	/*fiber contribution*/
	ComputeFiberStretch(fC, fFiberStretch);
	ComputeFiberStress(fFiberStretch, fFiberStress);
	
	/* rotate stress to lab coordinates */
	AssembleFiberStress(fFiberStress, fStress);
	
	return(fStress);
}

/* material description */
const dMatrixT& FSFiberMatT::c_ijkl(void)
{
	/* deformation gradient */
	const dMatrixT& Fmat = F();
	
	/* transform */
	fModulus.SetToScaled(1.0/Fmat.Det(), PushForward(Fmat, C_IJKL()));
	return fModulus;
}
/**< \todo construct directly in material description */

const dSymMatrixT& FSFiberMatT::s_ij(void)
{
	/* deformation gradient */
	const dMatrixT& Fmat = F();
	
	/* transform */
	fStress.SetToScaled(1.0/Fmat.Det(), PushForward(Fmat, S_IJ()));
	return fStress;
}
	
/* describe the parameters needed by the interface */
void FSFiberMatT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	FSSolidMatT::DefineParameters(list);
	
	/* 2D option must be plain stress */
	ParameterT& constraint = list.GetParameter("constraint_2D");
	constraint.SetDefault(kPlaneStrain);
}

/* accept parameter list */
void FSFiberMatT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	FSSolidMatT::TakeParameterList(list);

	/* dimension work space */
	fC.Dimension(fNumSD);
	fModulus.Dimension(dSymMatrixT::NumValues(fNumSD));
	fStress.Dimension(fNumSD);

	/*assuming in-plane fibers*/
	fQ.Dimension(fNumSD);
	
	/*initialize*/
	fQ = 0.0;
	fQ[0] = 1.0;
	fQ[4] = 1.0;
	fQ[8] = 1.0;
	
	fNumFibStress = dSymMatrixT::NumValues(fNumSD);
	
	/* allocate memory */
	/*2D fiber stress and modulus*/
	fFiberStretch.Dimension(fNumSD);
	fFiberStress.Dimension(fNumSD);
	fFiberMod.Dimension(fNumFibStress);
	
}

/*********************************************************************************************
 *protected                                                                                  *
 *********************************************************************************************/
const dMatrixT& FSFiberMatT::GetRotation()
{
	const char caller[] = "FSFiberMatT::GetRotation";
	const dArray2DT& Fibers = FiberMatSupportT().Fiber_Vec();
		
	int num_fibers = Fibers.MajorDim();

	num_fibers--;
	
	if (num_fibers > 0)
	{
		int fiber_num = 0; /*build orthogonal basis from the first fiber vector and normal vector*/
		fQ(0,0) = Fibers(fiber_num,0);
		fQ(1,0) = Fibers(fiber_num,1);
		fQ(2,0) = Fibers(fiber_num,2);
			
		
		fQ(0,2) = Fibers(num_fibers,0);
		fQ(1,2) = Fibers(num_fibers,1);
		fQ(2,2) = Fibers(num_fibers,2);

		const double* A = fQ(0);
		const double* N = fQ(2);  /*normal vector*/

		fQ(0,1) = N[1]*A[2] - N[2]*A[1];
		fQ(1,1) = N[2]*A[0] - N[0]*A[2];
		fQ(2,1) = N[0]*A[1] - N[1]*A[0];
		
		const double* B = fQ(1);  /*right hand rule*/
        fQ(0,2) = A[1]*B[2] - A[2]*B[1];
        fQ(1,2) = A[2]*B[0] - A[0]*B[2];
        fQ(2,2) = A[0]*B[1] - A[1]*B[0];
		
	}
	else 
		ExceptionT::GeneralFail(caller, "empty side set: num sides = 0");
		
	return(fQ);
}

/* construct symmetric rank-4 mixed-direction tensor (6.1.44) */
void FSFiberMatT::MixedRank4_2D(const dArrayT& a, const dArrayT& b, 
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

void FSFiberMatT::MixedRank4_3D(const dArrayT& a, const dArrayT& b, 
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
	
	
//	cout << "a" << a[0] <<endl;
//	cout << "rank4_ab" << rank4_ab <<endl;
	
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


void FSFiberMatT::ComputeFiberStretch(const dSymMatrixT& C, dSymMatrixT& Cf)
{
	/*Rotate stretch to plane of fibril families
	  C* = QT C Q (Cf_ab = QIa CIJ QJb)
	*/
	
	const dMatrixT& Q = GetRotation();
	
	if (fNumFibStress == 6)
	{
		Cf.MultQTBQ(Q, C);
		
	}
	else if (fNumFibStress ==3)
	{
		const double* x = Q(0);
		const double* y = Q(1);

		Cf[0] = C[0]*x[0]*x[0] + C[1]*x[1]*x[1] + C[2]*x[2]*x[2] + 2.0*(C[3]*x[1]*x[2] + C[4]*x[0]*x[2] + C[5]*x[0]*x[1]);
		Cf[1] = C[0]*y[0]*y[0] + C[1]*y[1]*y[1] + C[2]*y[2]*y[2] + 2.0*(C[3]*y[1]*y[2] + C[4]*y[0]*y[2] + C[5]*y[0]*y[1]);
		Cf[2] = C[0]*x[0]*y[0] + C[1]*x[1]*y[1] + C[2]*x[2]*y[2] + C[3]*(x[1]*y[2] + y[1]*x[2]) + C[4]*(x[0]*y[2] + y[0]*x[2]) + C[5]*(x[0]*y[1] + y[0]*x[1]);
		
		//cout<<setprecision(16)<<"\nCompute_Fiber_Stretch"<<Cf;
		//cout<<setprecision(16)<<"\nStretch"<<C;

	}
}

void FSFiberMatT::AssembleFiberStress(const dSymMatrixT& sigf, dSymMatrixT& sig, const int fillmode)
{
	/* sig = Q.sf.QT */
	if (fillmode == dSymMatrixT::kOverwrite)
		sig = 0.0;
		
	const dMatrixT& Q = GetRotation();
	
	if (fNumFibStress == 6)
	{
		/* transform */
		sig[0] += sigf[0]*Q[0]*Q[0] + 2*sigf[5]*Q[0]*Q[3] + sigf[1]*Q[3]*Q[3] +
		          2*sigf[4]*Q[0]*Q[6] + 2*sigf[3]*Q[3]*Q[6] + sigf[2]*Q[6]*Q[6]; 
		sig[1] += sigf[0]*Q[1]*Q[1] + 2*sigf[5]*Q[1]*Q[4] + sigf[1]*Q[4]*Q[4] +
		          2*sigf[4]*Q[1]*Q[7] + 2*sigf[3]*Q[4]*Q[7] + sigf[2]*Q[7]*Q[7];
		sig[2] += sigf[0]*Q[2]*Q[2] + 2*sigf[5]*Q[2]*Q[5] + sigf[1]*Q[5]*Q[5] +
		          2*sigf[4]*Q[2]*Q[8] + 2*sigf[3]*Q[5]*Q[8] + sigf[2]*Q[8]*Q[8];
		
		sig[3] += Q[5]*(sigf[5]*Q[1] + sigf[1]*Q[4] + sigf[3]*Q[7]) +
		            Q[2]*(sigf[0]*Q[1] + sigf[5]*Q[4] + sigf[4]*Q[7]) +
		            Q[8]*(sigf[4]*Q[1] + sigf[3]*Q[4] + sigf[2]*Q[7]); 
		sig[4] += Q[5]*(sigf[5]*Q[0] + sigf[1]*Q[3] + sigf[3]*Q[6]) +
		            Q[2]*(sigf[0]*Q[0] + sigf[5]*Q[3] + sigf[4]*Q[6]) +
		            Q[8]*(sigf[4]*Q[0] + sigf[3]*Q[3] + sigf[2]*Q[6]); 
		sig[5] += Q[4]*(sigf[5]*Q[0] + sigf[1]*Q[3] + sigf[3]*Q[6]) +
		            Q[1]*(sigf[0]*Q[0] + sigf[5]*Q[3] + sigf[4]*Q[6]) +
		            Q[7]*(sigf[4]*Q[0] + sigf[3]*Q[3] + sigf[2]*Q[6]);

	}
	else if (fNumFibStress == 3)
	{
		const double& s1 = sigf[0]; /*sf11*/
		const double& s2 = sigf[1]; /*sf22*/
		const double& s3 = sigf[2]; /*sf12*/
	
		const double* x = Q(0);	
		const double* y = Q(1);

		
		/* sig = s11 a1 x a1 + s22 a2 x a2 + s12 (a1 x a2 + a2 x a1) */
		
		sig[0] += s1*x[0]*x[0] + s2*y[0]*y[0] + 2.0*s3*x[0]*y[0];
		sig[1] += s1*x[1]*x[1] + s2*y[1]*y[1] + 2.0*s3*x[1]*y[1];
		sig[2] += s1*x[2]*x[2] + s2*y[2]*y[2] + 2.0*s3*x[2]*y[2];
		sig[3] += s1*x[1]*x[2] + s2*y[1]*y[2] + s3*(x[1]*y[2] + x[2]*y[1]);
		sig[4] += s1*x[0]*x[2] + s2*y[0]*y[2] + s3*(x[0]*y[2] + x[2]*y[0]);
		sig[5] += s1*x[0]*x[1] + s2*y[0]*y[1] + s3*(x[0]*y[1] + x[1]*y[0]);


	}
	
	
}

void FSFiberMatT::AssembleFiberModuli(const dMatrixT& cf, dMatrixT& mod, const int fillmode)
{
	/*Rotate moduli from local frame (defined by fibrils) to global cartesian frame*/
	/*C_IJKL = QIa QJb QKc QLd Cf_abcd*/

	if (fillmode == dSymMatrixT::kOverwrite)
	mod = 0.0;

	const dMatrixT& Q = GetRotation();

	if (fNumFibStress == 6)
	{
		const double& c11 = cf(0,0);
		const double& c22 = cf(1,1);
		const double& c33 = cf(2,2);
		const double& c44 = cf(3,3);
		const double& c55 = cf(4,4);
		const double& c66 = cf(5,5);
	
		const double& c12 = cf(0,1);
		const double& c13 = cf(0,2);
		const double& c14 = cf(0,3);
		const double& c15 = cf(0,4);
		const double& c16 = cf(0,5);
		
		const double& c21 = cf(1,0);
		const double& c23 = cf(1,2);
		const double& c24 = cf(1,3);
		const double& c25 = cf(1,4);
		const double& c26 = cf(1,5);

		const double& c31 = cf(2,0);
		const double& c32 = cf(2,1);
		const double& c34 = cf(2,3);
		const double& c35 = cf(2,4);
		const double& c36 = cf(2,5);

		const double& c41 = cf(3,0);
		const double& c42 = cf(3,1);
		const double& c43 = cf(3,2);
		const double& c45 = cf(3,4);
		const double& c46 = cf(3,5);

		const double& c51 = cf(4,0);
		const double& c52 = cf(4,1);
		const double& c53 = cf(4,2);
		const double& c54 = cf(4,3);
		const double& c56 = cf(4,5);

		const double& c61 = cf(5,0);
		const double& c62 = cf(5,1);
		const double& c63 = cf(5,2);
		const double& c64 = cf(5,3);
		const double& c65 = cf(5,4);
	
		const double* x = Q(0);	
		const double* y = Q(1);
		const double* z = Q(2);

		mod(0,0) += c11*x[0]*x[0]*x[0]*x[0] + 2*c16*x[0]*x[0]*x[0]*y[0] + 2*c61*x[0]*x[0]*x[0]*y[0] + c12*x[0]*x[0]*y[0]*y[0] 
			+ c21*x[0]*x[0]*y[0]*y[0] + 4*c66*x[0]*x[0]*y[0]*y[0] + 2*c26*x[0]*y[0]*y[0]*y[0] 
			+ 2*c62*x[0]*y[0]*y[0]*y[0] + c22*y[0]*y[0]*y[0]*y[0] + 2*c15*x[0]*x[0]*x[0]*z[0] + 2*c51*x[0]*x[0]*x[0]*z[0] 
			+ 2*c14*x[0]*x[0]*y[0]*z[0] + 2*c41*x[0]*x[0]*y[0]*z[0] + 4*c56*x[0]*x[0]*y[0]*z[0] 
			+ 4*c65*x[0]*x[0]*y[0]*z[0] + 2*c25*x[0]*y[0]*y[0]*z[0] + 4*c46*x[0]*y[0]*y[0]*z[0] 
			+ 2*c52*x[0]*y[0]*y[0]*z[0] + 4*c64*x[0]*y[0]*y[0]*z[0] + 2*c24*y[0]*y[0]*y[0]*z[0] 
			+ 2*c42*y[0]*y[0]*y[0]*z[0] + c13*x[0]*x[0]*z[0]*z[0] + c31*x[0]*x[0]*z[0]*z[0] 
			+ 4*c55*x[0]*x[0]*z[0]*z[0] + 2*c36*x[0]*y[0]*z[0]*z[0] + 4*c45*x[0]*y[0]*z[0]*z[0] 
			+ 4*c54*x[0]*y[0]*z[0]*z[0] + 2*c63*x[0]*y[0]*z[0]*z[0] + c23*y[0]*y[0]*z[0]*z[0] 
			+ c32*y[0]*y[0]*z[0]*z[0] + 4*c44*y[0]*y[0]*z[0]*z[0] + 2*c35*x[0]*z[0]*z[0]*z[0] 
			+ 2*c53*x[0]*z[0]*z[0]*z[0] + 2*c34*y[0]*z[0]*z[0]*z[0] + 2*c43*y[0]*z[0]*z[0]*z[0] + c33*z[0]*z[0]*z[0]*z[0];
      
		mod(1,1) += c11*x[1]*x[1]*x[1]*x[1] + 2*c16*x[1]*x[1]*x[1]*y[1] + 2*c61*x[1]*x[1]*x[1]*y[1] + c12*x[1]*x[1]*y[1]*y[1] 
			+ c21*x[1]*x[1]*y[1]*y[1] + 4*c66*x[1]*x[1]*y[1]*y[1] + 2*c26*x[1]*y[1]*y[1]*y[1] 
			+ 2*c62*x[1]*y[1]*y[1]*y[1] + c22*y[1]*y[1]*y[1]*y[1] + 2*c15*x[1]*x[1]*x[1]*z[1] + 2*c51*x[1]*x[1]*x[1]*z[1] 
			+ 2*c14*x[1]*x[1]*y[1]*z[1] + 2*c41*x[1]*x[1]*y[1]*z[1] + 4*c56*x[1]*x[1]*y[1]*z[1] 
			+ 4*c65*x[1]*x[1]*y[1]*z[1] + 2*c25*x[1]*y[1]*y[1]*z[1] + 4*c46*x[1]*y[1]*y[1]*z[1] 
			+ 2*c52*x[1]*y[1]*y[1]*z[1] + 4*c64*x[1]*y[1]*y[1]*z[1] + 2*c24*y[1]*y[1]*y[1]*z[1] 
			+ 2*c42*y[1]*y[1]*y[1]*z[1] + c13*x[1]*x[1]*z[1]*z[1] + c31*x[1]*x[1]*z[1]*z[1] 
			+ 4*c55*x[1]*x[1]*z[1]*z[1] + 2*c36*x[1]*y[1]*z[1]*z[1] + 4*c45*x[1]*y[1]*z[1]*z[1] 
			+ 4*c54*x[1]*y[1]*z[1]*z[1] + 2*c63*x[1]*y[1]*z[1]*z[1] + c23*y[1]*y[1]*z[1]*z[1] 
			+ c32*y[1]*y[1]*z[1]*z[1] + 4*c44*y[1]*y[1]*z[1]*z[1] + 2*c35*x[1]*z[1]*z[1]*z[1] 
			+ 2*c53*x[1]*z[1]*z[1]*z[1] + 2*c34*y[1]*z[1]*z[1]*z[1] + 2*c43*y[1]*z[1]*z[1]*z[1] + c33*z[1]*z[1]*z[1]*z[1];
   
		mod(2,2) += c11*x[2]*x[2]*x[2]*x[2] + 2*c16*x[2]*x[2]*x[2]*y[2] + 2*c61*x[2]*x[2]*x[2]*y[2] + c12*x[2]*x[2]*y[2]*y[2] 
			+ c21*x[2]*x[2]*y[2]*y[2] + 4*c66*x[2]*x[2]*y[2]*y[2] + 2*c26*x[2]*y[2]*y[2]*y[2] 
			+ 2*c62*x[2]*y[2]*y[2]*y[2] + c22*y[2]*y[2]*y[2]*y[2] + 2*c15*x[2]*x[2]*x[2]*z[2] + 2*c51*x[2]*x[2]*x[2]*z[2] 
			+ 2*c14*x[2]*x[2]*y[2]*z[2] + 2*c41*x[2]*x[2]*y[2]*z[2] + 4*c56*x[2]*x[2]*y[2]*z[2] 
			+ 4*c65*x[2]*x[2]*y[2]*z[2] + 2*c25*x[2]*y[2]*y[2]*z[2] + 4*c46*x[2]*y[2]*y[2]*z[2] 
			+ 2*c52*x[2]*y[2]*y[2]*z[2] + 4*c64*x[2]*y[2]*y[2]*z[2] + 2*c24*y[2]*y[2]*y[2]*z[2] 
			+ 2*c42*y[2]*y[2]*y[2]*z[2] + c13*x[2]*x[2]*z[2]*z[2] + c31*x[2]*x[2]*z[2]*z[2] 
			+ 4*c55*x[2]*x[2]*z[2]*z[2] + 2*c36*x[2]*y[2]*z[2]*z[2] + 4*c45*x[2]*y[2]*z[2]*z[2] 
			+ 4*c54*x[2]*y[2]*z[2]*z[2] + 2*c63*x[2]*y[2]*z[2]*z[2] + c23*y[2]*y[2]*z[2]*z[2] 
			+ c32*y[2]*y[2]*z[2]*z[2] + 4*c44*y[2]*y[2]*z[2]*z[2] + 2*c35*x[2]*z[2]*z[2]*z[2] 
			+ 2*c53*x[2]*z[2]*z[2]*z[2] + 2*c34*y[2]*z[2]*z[2]*z[2] + 2*c43*y[2]*z[2]*z[2]*z[2] + c33*z[2]*z[2]*z[2]*z[2];
   
		mod(3,3) += c11*x[1]*x[1]*x[2]*x[2] + c16*x[1]*x[2]*x[2]*y[1] + c61*x[1]*x[2]*x[2]*y[1] 
		+ c66*x[2]*x[2]*y[1]*y[1] + c16*x[1]*x[1]*x[2]*y[2] + c61*x[1]*x[1]*x[2]*y[2] 
		+ c12*x[1]*x[2]*y[1]*y[2] + c21*x[1]*x[2]*y[1]*y[2] + 2*c66*x[1]*x[2]*y[1]*y[2] + c26*x[2]*y[1]*y[1]*y[2] 
		+ c62*x[2]*y[1]*y[1]*y[2] + c66*x[1]*x[1]*y[2]*y[2] + c26*x[1]*y[1]*y[2]*y[2] 
		+ c62*x[1]*y[1]*y[2]*y[2] + c22*y[1]*y[1]*y[2]*y[2] + c15*x[1]*x[2]*x[2]*z[1] 
		+ c51*x[1]*x[2]*x[2]*z[1] + c56*x[2]*x[2]*y[1]*z[1] + c65*x[2]*x[2]*y[1]*z[1] + c14*x[1]*x[2]*y[2]*z[1] 
		+ c41*x[1]*x[2]*y[2]*z[1] + c56*x[1]*x[2]*y[2]*z[1] + c65*x[1]*x[2]*y[2]*z[1] + c25*x[2]*y[1]*y[2]*z[1] + c46*x[2]*y[1]*y[2]*z[1] 
		+ c52*x[2]*y[1]*y[2]*z[1] + c64*x[2]*y[1]*y[2]*z[1] + c46*x[1]*y[2]*y[2]*z[1] + c64*x[1]*y[2]*y[2]*z[1] + c24*y[1]*y[2]*y[2]*z[1] 
		+ c42*y[1]*y[2]*y[2]*z[1] + c55*x[2]*x[2]*z[1]*z[1] + c45*x[2]*y[2]*z[1]*z[1] + c54*x[2]*y[2]*z[1]*z[1] + c44*y[2]*y[2]*z[1]*z[1] 
		+ c15*x[1]*x[1]*x[2]*z[2] + c51*x[1]*x[1]*x[2]*z[2] + c14*x[1]*x[2]*y[1]*z[2] + c41*x[1]*x[2]*y[1]*z[2] + c56*x[1]*x[2]*y[1]*z[2] 
		+ c65*x[1]*x[2]*y[1]*z[2] + c46*x[2]*y[1]*y[1]*z[2] + c64*x[2]*y[1]*y[1]*z[2] + c56*x[1]*x[1]*y[2]*z[2] + c65*x[1]*x[1]*y[2]*z[2] 
		+ c25*x[1]*y[1]*y[2]*z[2] + c46*x[1]*y[1]*y[2]*z[2] + c52*x[1]*y[1]*y[2]*z[2] + c64*x[1]*y[1]*y[2]*z[2] + c24*y[1]*y[1]*y[2]*z[2] 
		+ c42*y[1]*y[1]*y[2]*z[2] + c13*x[1]*x[2]*z[1]*z[2] + c31*x[1]*x[2]*z[1]*z[2] + 2*c55*x[1]*x[2]*z[1]*z[2] + c36*x[2]*y[1]*z[1]*z[2] 
		+ c45*x[2]*y[1]*z[1]*z[2] + c54*x[2]*y[1]*z[1]*z[2] + c63*x[2]*y[1]*z[1]*z[2] + c36*x[1]*y[2]*z[1]*z[2] + c45*x[1]*y[2]*z[1]*z[2] 
		+ c54*x[1]*y[2]*z[1]*z[2] + c63*x[1]*y[2]*z[1]*z[2] + c23*y[1]*y[2]*z[1]*z[2] + c32*y[1]*y[2]*z[1]*z[2] + 2*c44*y[1]*y[2]*z[1]*z[2] 
		+ c35*x[2]*z[1]*z[1]*z[2] + c53*x[2]*z[1]*z[1]*z[2] + c34*y[2]*z[1]*z[1]*z[2] + c43*y[2]*z[1]*z[1]*z[2] + c55*x[1]*x[1]*z[2]*z[2] 
		+ c45*x[1]*y[1]*z[2]*z[2] + c54*x[1]*y[1]*z[2]*z[2] + c44*y[1]*y[1]*z[2]*z[2] + c35*x[1]*z[1]*z[2]*z[2] + c53*x[1]*z[1]*z[2]*z[2] 
		+ c34*y[1]*z[1]*z[2]*z[2] + c43*y[1]*z[1]*z[2]*z[2] + c33*z[1]*z[1]*z[2]*z[2];
   
		mod(4,4) += c11*x[0]*x[0]*x[2]*x[2] + c16*x[0]*x[2]*x[2]*y[0] + c61*x[0]*x[2]*x[2]*y[0] + 
			c66*x[2]*x[2]*y[0]*y[0] + c16*x[0]*x[0]*x[2]*y[2] + c61*x[0]*x[0]*x[2]*y[2] + 
			c12*x[0]*x[2]*y[0]*y[2] + c21*x[0]*x[2]*y[0]*y[2] + 2*c66*x[0]*x[2]*y[0]*y[2] + c26*x[2]*y[0]*y[0]*y[2] + 
			c62*x[2]*y[0]*y[0]*y[2] + c66*x[0]*x[0]*y[2]*y[2] + c26*x[0]*y[0]*y[2]*y[2] + 
			c62*x[0]*y[0]*y[2]*y[2] + c22*y[0]*y[0]*y[2]*y[2] + c15*x[0]*x[2]*x[2]*z[0] + 
			c51*x[0]*x[2]*x[2]*z[0] + c56*x[2]*x[2]*y[0]*z[0] + c65*x[2]*x[2]*y[0]*z[0] + c14*x[0]*x[2]*y[2]*z[0] + 
			c41*x[0]*x[2]*y[2]*z[0] + c56*x[0]*x[2]*y[2]*z[0] + c65*x[0]*x[2]*y[2]*z[0] + c25*x[2]*y[0]*y[2]*z[0] + c46*x[2]*y[0]*y[2]*z[0] + 
			c52*x[2]*y[0]*y[2]*z[0] + c64*x[2]*y[0]*y[2]*z[0] + c46*x[0]*y[2]*y[2]*z[0] + c64*x[0]*y[2]*y[2]*z[0] + 
			c24*y[0]*y[2]*y[2]*z[0] + c42*y[0]*y[2]*y[2]*z[0] + c55*x[2]*x[2]*z[0]*z[0] + 
			c45*x[2]*y[2]*z[0]*z[0] + c54*x[2]*y[2]*z[0]*z[0] + c44*y[2]*y[2]*z[0]*z[0] + 
			c15*x[0]*x[0]*x[2]*z[2] + c51*x[0]*x[0]*x[2]*z[2] + c14*x[0]*x[2]*y[0]*z[2] + c41*x[0]*x[2]*y[0]*z[2] + 
			c56*x[0]*x[2]*y[0]*z[2] + c65*x[0]*x[2]*y[0]*z[2] + c46*x[2]*y[0]*y[0]*z[2] + c64*x[2]*y[0]*y[0]*z[2] + 
			c56*x[0]*x[0]*y[2]*z[2] + c65*x[0]*x[0]*y[2]*z[2] + c25*x[0]*y[0]*y[2]*z[2] + c46*x[0]*y[0]*y[2]*z[2] + 
			c52*x[0]*y[0]*y[2]*z[2] + c64*x[0]*y[0]*y[2]*z[2] + c24*y[0]*y[0]*y[2]*z[2] + c42*y[0]*y[0]*y[2]*z[2] + 
			c13*x[0]*x[2]*z[0]*z[2] + c31*x[0]*x[2]*z[0]*z[2] + 2*c55*x[0]*x[2]*z[0]*z[2] + c36*x[2]*y[0]*z[0]*z[2] + 
			c45*x[2]*y[0]*z[0]*z[2] + c54*x[2]*y[0]*z[0]*z[2] + c63*x[2]*y[0]*z[0]*z[2] + c36*x[0]*y[2]*z[0]*z[2] + c45*x[0]*y[2]*z[0]*z[2] + 
			c54*x[0]*y[2]*z[0]*z[2] + c63*x[0]*y[2]*z[0]*z[2] + c23*y[0]*y[2]*z[0]*z[2] + c32*y[0]*y[2]*z[0]*z[2] + 
			2*c44*y[0]*y[2]*z[0]*z[2] + c35*x[2]*z[0]*z[0]*z[2] + c53*x[2]*z[0]*z[0]*z[2] + 
			c34*y[2]*z[0]*z[0]*z[2] + c43*y[2]*z[0]*z[0]*z[2] + c55*x[0]*x[0]*z[2]*z[2] + 
			c45*x[0]*y[0]*z[2]*z[2] + c54*x[0]*y[0]*z[2]*z[2] + c44*y[0]*y[0]*z[2]*z[2] + 
			c35*x[0]*z[0]*z[2]*z[2] + c53*x[0]*z[0]*z[2]*z[2] + c34*y[0]*z[0]*z[2]*z[2] + 
			c43*y[0]*z[0]*z[2]*z[2] + c33*z[0]*z[0]*z[2]*z[2];
   
		mod(5,5) += c11*x[0]*x[0]*x[1]*x[1] + c16*x[0]*x[1]*x[1]*y[0] + c61*x[0]*x[1]*x[1]*y[0] + 
   c66*x[1]*x[1]*y[0]*y[0] + c16*x[0]*x[0]*x[1]*y[1] + c61*x[0]*x[0]*x[1]*y[1] + 
   c12*x[0]*x[1]*y[0]*y[1] + c21*x[0]*x[1]*y[0]*y[1] + 2*c66*x[0]*x[1]*y[0]*y[1] + c26*x[1]*y[0]*y[0]*y[1] + 
   c62*x[1]*y[0]*y[0]*y[1] + c66*x[0]*x[0]*y[1]*y[1] + c26*x[0]*y[0]*y[1]*y[1] + 
   c62*x[0]*y[0]*y[1]*y[1] + c22*y[0]*y[0]*y[1]*y[1] + c15*x[0]*x[1]*x[1]*z[0] + 
   c51*x[0]*x[1]*x[1]*z[0] + c56*x[1]*x[1]*y[0]*z[0] + c65*x[1]*x[1]*y[0]*z[0] + c14*x[0]*x[1]*y[1]*z[0] + 
   c41*x[0]*x[1]*y[1]*z[0] + c56*x[0]*x[1]*y[1]*z[0] + c65*x[0]*x[1]*y[1]*z[0] + c25*x[1]*y[0]*y[1]*z[0] + c46*x[1]*y[0]*y[1]*z[0] + 
   c52*x[1]*y[0]*y[1]*z[0] + c64*x[1]*y[0]*y[1]*z[0] + c46*x[0]*y[1]*y[1]*z[0] + c64*x[0]*y[1]*y[1]*z[0] + 
   c24*y[0]*y[1]*y[1]*z[0] + c42*y[0]*y[1]*y[1]*z[0] + c55*x[1]*x[1]*z[0]*z[0] + 
   c45*x[1]*y[1]*z[0]*z[0] + c54*x[1]*y[1]*z[0]*z[0] + c44*y[1]*y[1]*z[0]*z[0] + 
   c15*x[0]*x[0]*x[1]*z[1] + c51*x[0]*x[0]*x[1]*z[1] + c14*x[0]*x[1]*y[0]*z[1] + c41*x[0]*x[1]*y[0]*z[1] + 
   c56*x[0]*x[1]*y[0]*z[1] + c65*x[0]*x[1]*y[0]*z[1] + c46*x[1]*y[0]*y[0]*z[1] + c64*x[1]*y[0]*y[0]*z[1] + 
   c56*x[0]*x[0]*y[1]*z[1] + c65*x[0]*x[0]*y[1]*z[1] + c25*x[0]*y[0]*y[1]*z[1] + c46*x[0]*y[0]*y[1]*z[1] + 
   c52*x[0]*y[0]*y[1]*z[1] + c64*x[0]*y[0]*y[1]*z[1] + c24*y[0]*y[0]*y[1]*z[1] + c42*y[0]*y[0]*y[1]*z[1] + 
   c13*x[0]*x[1]*z[0]*z[1] + c31*x[0]*x[1]*z[0]*z[1] + 2*c55*x[0]*x[1]*z[0]*z[1] + c36*x[1]*y[0]*z[0]*z[1] + 
   c45*x[1]*y[0]*z[0]*z[1] + c54*x[1]*y[0]*z[0]*z[1] + c63*x[1]*y[0]*z[0]*z[1] + c36*x[0]*y[1]*z[0]*z[1] + c45*x[0]*y[1]*z[0]*z[1] + 
   c54*x[0]*y[1]*z[0]*z[1] + c63*x[0]*y[1]*z[0]*z[1] + c23*y[0]*y[1]*z[0]*z[1] + c32*y[0]*y[1]*z[0]*z[1] + 
   2*c44*y[0]*y[1]*z[0]*z[1] + c35*x[1]*z[0]*z[0]*z[1] + c53*x[1]*z[0]*z[0]*z[1] + 
   c34*y[1]*z[0]*z[0]*z[1] + c43*y[1]*z[0]*z[0]*z[1] + c55*x[0]*x[0]*z[1]*z[1] + 
   c45*x[0]*y[0]*z[1]*z[1] + c54*x[0]*y[0]*z[1]*z[1] + c44*y[0]*y[0]*z[1]*z[1] + 
   c35*x[0]*z[0]*z[1]*z[1] + c53*x[0]*z[0]*z[1]*z[1] + c34*y[0]*z[0]*z[1]*z[1] + 
   c43*y[0]*z[0]*z[1]*z[1] + c33*z[0]*z[0]*z[1]*z[1];
   
		mod(0,1) += c11*x[0]*x[0]*x[1]*x[1] + 2*c61*x[0]*x[1]*x[1]*y[0] + c21*x[1]*x[1]*y[0]*y[0] + 
   2*c16*x[0]*x[0]*x[1]*y[1] + 4*c66*x[0]*x[1]*y[0]*y[1] + 2*c26*x[1]*y[0]*y[0]*y[1] + 
   c12*x[0]*x[0]*y[1]*y[1] + 2*c62*x[0]*y[0]*y[1]*y[1] + c22*y[0]*y[0]*y[1]*y[1] + 
   2*c51*x[0]*x[1]*x[1]*z[0] + 2*c41*x[1]*x[1]*y[0]*z[0] + 4*c56*x[0]*x[1]*y[1]*z[0] + 
   4*c46*x[1]*y[0]*y[1]*z[0] + 2*c52*x[0]*y[1]*y[1]*z[0] + 2*c42*y[0]*y[1]*y[1]*z[0] + 
   c31*x[1]*x[1]*z[0]*z[0] + 2*c36*x[1]*y[1]*z[0]*z[0] + c32*y[1]*y[1]*z[0]*z[0] + 
   2*c15*x[0]*x[0]*x[1]*z[1] + 4*c65*x[0]*x[1]*y[0]*z[1] + 2*c25*x[1]*y[0]*y[0]*z[1] + 
   2*c14*x[0]*x[0]*y[1]*z[1] + 4*c64*x[0]*y[0]*y[1]*z[1] + 2*c24*y[0]*y[0]*y[1]*z[1] + 
   4*c55*x[0]*x[1]*z[0]*z[1] + 4*c45*x[1]*y[0]*z[0]*z[1] + 4*c54*x[0]*y[1]*z[0]*z[1] + 4*c44*y[0]*y[1]*z[0]*z[1] + 
   2*c35*x[1]*z[0]*z[0]*z[1] + 2*c34*y[1]*z[0]*z[0]*z[1] + c13*x[0]*x[0]*z[1]*z[1] + 
   2*c63*x[0]*y[0]*z[1]*z[1] + c23*y[0]*y[0]*z[1]*z[1] + 2*c53*x[0]*z[0]*z[1]*z[1] + 
   2*c43*y[0]*z[0]*z[1]*z[1] + c33*z[0]*z[0]*z[1]*z[1];
   
		mod(0,2) += c11*x[0]*x[0]*x[2]*x[2] + 2*c61*x[0]*x[2]*x[2]*y[0] + c21*x[2]*x[2]*y[0]*y[0] + 
   2*c16*x[0]*x[0]*x[2]*y[2] + 4*c66*x[0]*x[2]*y[0]*y[2] + 2*c26*x[2]*y[0]*y[0]*y[2] + 
   c12*x[0]*x[0]*y[2]*y[2] + 2*c62*x[0]*y[0]*y[2]*y[2] + c22*y[0]*y[0]*y[2]*y[2] + 
   2*c51*x[0]*x[2]*x[2]*z[0] + 2*c41*x[2]*x[2]*y[0]*z[0] + 4*c56*x[0]*x[2]*y[2]*z[0] + 
   4*c46*x[2]*y[0]*y[2]*z[0] + 2*c52*x[0]*y[2]*y[2]*z[0] + 2*c42*y[0]*y[2]*y[2]*z[0] + 
   c31*x[2]*x[2]*z[0]*z[0] + 2*c36*x[2]*y[2]*z[0]*z[0] + c32*y[2]*y[2]*z[0]*z[0] + 
   2*c15*x[0]*x[0]*x[2]*z[2] + 4*c65*x[0]*x[2]*y[0]*z[2] + 2*c25*x[2]*y[0]*y[0]*z[2] + 
   2*c14*x[0]*x[0]*y[2]*z[2] + 4*c64*x[0]*y[0]*y[2]*z[2] + 2*c24*y[0]*y[0]*y[2]*z[2] + 
   4*c55*x[0]*x[2]*z[0]*z[2] + 4*c45*x[2]*y[0]*z[0]*z[2] + 4*c54*x[0]*y[2]*z[0]*z[2] + 4*c44*y[0]*y[2]*z[0]*z[2] + 
   2*c35*x[2]*z[0]*z[0]*z[2] + 2*c34*y[2]*z[0]*z[0]*z[2] + c13*x[0]*x[0]*z[2]*z[2] + 
   2*c63*x[0]*y[0]*z[2]*z[2] + c23*y[0]*y[0]*z[2]*z[2] + 2*c53*x[0]*z[0]*z[2]*z[2] + 
   2*c43*y[0]*z[0]*z[2]*z[2] + c33*z[0]*z[0]*z[2]*z[2];
   
		mod(0,3) += c11*x[0]*x[0]*x[1]*x[2] + 2*c61*x[0]*x[1]*x[2]*y[0] + c21*x[1]*x[2]*y[0]*y[0] + c16*x[0]*x[0]*x[2]*y[1] + 
   2*c66*x[0]*x[2]*y[0]*y[1] + c26*x[2]*y[0]*y[0]*y[1] + c16*x[0]*x[0]*x[1]*y[2] + 2*c66*x[0]*x[1]*y[0]*y[2] + 
   c26*x[1]*y[0]*y[0]*y[2] + c12*x[0]*x[0]*y[1]*y[2] + 2*c62*x[0]*y[0]*y[1]*y[2] + 
   c22*y[0]*y[0]*y[1]*y[2] + 2*c51*x[0]*x[1]*x[2]*z[0] + 2*c41*x[1]*x[2]*y[0]*z[0] + 2*c56*x[0]*x[2]*y[1]*z[0] + 
   2*c46*x[2]*y[0]*y[1]*z[0] + 2*c56*x[0]*x[1]*y[2]*z[0] + 2*c46*x[1]*y[0]*y[2]*z[0] + 2*c52*x[0]*y[1]*y[2]*z[0] + 
   2*c42*y[0]*y[1]*y[2]*z[0] + c31*x[1]*x[2]*z[0]*z[0] + c36*x[2]*y[1]*z[0]*z[0] + 
   c36*x[1]*y[2]*z[0]*z[0] + c32*y[1]*y[2]*z[0]*z[0] + c15*x[0]*x[0]*x[2]*z[1] + 
   2*c65*x[0]*x[2]*y[0]*z[1] + c25*x[2]*y[0]*y[0]*z[1] + c14*x[0]*x[0]*y[2]*z[1] + 2*c64*x[0]*y[0]*y[2]*z[1] + 
   c24*y[0]*y[0]*y[2]*z[1] + 2*c55*x[0]*x[2]*z[0]*z[1] + 2*c45*x[2]*y[0]*z[0]*z[1] + 2*c54*x[0]*y[2]*z[0]*z[1] + 
   2*c44*y[0]*y[2]*z[0]*z[1] + c35*x[2]*z[0]*z[0]*z[1] + c34*y[2]*z[0]*z[0]*z[1] + 
   c15*x[0]*x[0]*x[1]*z[2] + 2*c65*x[0]*x[1]*y[0]*z[2] + c25*x[1]*y[0]*y[0]*z[2] + 
   c14*x[0]*x[0]*y[1]*z[2] + 2*c64*x[0]*y[0]*y[1]*z[2] + c24*y[0]*y[0]*y[1]*z[2] + 2*c55*x[0]*x[1]*z[0]*z[2] + 
   2*c45*x[1]*y[0]*z[0]*z[2] + 2*c54*x[0]*y[1]*z[0]*z[2] + 2*c44*y[0]*y[1]*z[0]*z[2] + c35*x[1]*z[0]*z[0]*z[2] + 
   c34*y[1]*z[0]*z[0]*z[2] + c13*x[0]*x[0]*z[1]*z[2] + 2*c63*x[0]*y[0]*z[1]*z[2] + 
   c23*y[0]*y[0]*z[1]*z[2] + 2*c53*x[0]*z[0]*z[1]*z[2] + 2*c43*y[0]*z[0]*z[1]*z[2] + c33*z[0]*z[0]*z[1]*z[2];
   
		mod(0,4) += c11*x[0]*x[0]*x[0]*x[2] + c16*x[0]*x[0]*x[2]*y[0] + 2*c61*x[0]*x[0]*x[2]*y[0] + 
   c21*x[0]*x[2]*y[0]*y[0] + 2*c66*x[0]*x[2]*y[0]*y[0] + c26*x[2]*y[0]*y[0]*y[0] + 
   c16*x[0]*x[0]*x[0]*y[2] + c12*x[0]*x[0]*y[0]*y[2] + 2*c66*x[0]*x[0]*y[0]*y[2] + 
   c26*x[0]*y[0]*y[0]*y[2] + 2*c62*x[0]*y[0]*y[0]*y[2] + c22*y[0]*y[0]*y[0]*y[2] + 
   c15*x[0]*x[0]*x[2]*z[0] + 2*c51*x[0]*x[0]*x[2]*z[0] + 2*c41*x[0]*x[2]*y[0]*z[0] + 2*c56*x[0]*x[2]*y[0]*z[0] + 
   2*c65*x[0]*x[2]*y[0]*z[0] + c25*x[2]*y[0]*y[0]*z[0] + 2*c46*x[2]*y[0]*y[0]*z[0] + 
   c14*x[0]*x[0]*y[2]*z[0] + 2*c56*x[0]*x[0]*y[2]*z[0] + 2*c46*x[0]*y[0]*y[2]*z[0] + 2*c52*x[0]*y[0]*y[2]*z[0] + 
   2*c64*x[0]*y[0]*y[2]*z[0] + c24*y[0]*y[0]*y[2]*z[0] + 2*c42*y[0]*y[0]*y[2]*z[0] + 
   c31*x[0]*x[2]*z[0]*z[0] + 2*c55*x[0]*x[2]*z[0]*z[0] + c36*x[2]*y[0]*z[0]*z[0] + 
   2*c45*x[2]*y[0]*z[0]*z[0] + c36*x[0]*y[2]*z[0]*z[0] + 2*c54*x[0]*y[2]*z[0]*z[0] + 
   c32*y[0]*y[2]*z[0]*z[0] + 2*c44*y[0]*y[2]*z[0]*z[0] + c35*x[2]*z[0]*z[0]*z[0] + 
   c34*y[2]*z[0]*z[0]*z[0] + c15*x[0]*x[0]*x[0]*z[2] + c14*x[0]*x[0]*y[0]*z[2] + 
   2*c65*x[0]*x[0]*y[0]*z[2] + c25*x[0]*y[0]*y[0]*z[2] + 2*c64*x[0]*y[0]*y[0]*z[2] + 
   c24*y[0]*y[0]*y[0]*z[2] + c13*x[0]*x[0]*z[0]*z[2] + 2*c55*x[0]*x[0]*z[0]*z[2] + 2*c45*x[0]*y[0]*z[0]*z[2] + 
   2*c54*x[0]*y[0]*z[0]*z[2] + 2*c63*x[0]*y[0]*z[0]*z[2] + c23*y[0]*y[0]*z[0]*z[2] + 2*c44*y[0]*y[0]*z[0]*z[2] + 
   c35*x[0]*z[0]*z[0]*z[2] + 2*c53*x[0]*z[0]*z[0]*z[2] + c34*y[0]*z[0]*z[0]*z[2] + 
   2*c43*y[0]*z[0]*z[0]*z[2] + c33*z[0]*z[0]*z[0]*z[2];
   
		mod(0,5) += c11*x[0]*x[0]*x[0]*x[1] + c16*x[0]*x[0]*x[1]*y[0] + 2*c61*x[0]*x[0]*x[1]*y[0] + 
   c21*x[0]*x[1]*y[0]*y[0] + 2*c66*x[0]*x[1]*y[0]*y[0] + c26*x[1]*y[0]*y[0]*y[0] + 
   c16*x[0]*x[0]*x[0]*y[1] + c12*x[0]*x[0]*y[0]*y[1] + 2*c66*x[0]*x[0]*y[0]*y[1] + 
   c26*x[0]*y[0]*y[0]*y[1] + 2*c62*x[0]*y[0]*y[0]*y[1] + c22*y[0]*y[0]*y[0]*y[1] + 
   c15*x[0]*x[0]*x[1]*z[0] + 2*c51*x[0]*x[0]*x[1]*z[0] + 2*c41*x[0]*x[1]*y[0]*z[0] + 2*c56*x[0]*x[1]*y[0]*z[0] + 
   2*c65*x[0]*x[1]*y[0]*z[0] + c25*x[1]*y[0]*y[0]*z[0] + 2*c46*x[1]*y[0]*y[0]*z[0] + 
   c14*x[0]*x[0]*y[1]*z[0] + 2*c56*x[0]*x[0]*y[1]*z[0] + 2*c46*x[0]*y[0]*y[1]*z[0] + 2*c52*x[0]*y[0]*y[1]*z[0] + 
   2*c64*x[0]*y[0]*y[1]*z[0] + c24*y[0]*y[0]*y[1]*z[0] + 2*c42*y[0]*y[0]*y[1]*z[0] + 
   c31*x[0]*x[1]*z[0]*z[0] + 2*c55*x[0]*x[1]*z[0]*z[0] + c36*x[1]*y[0]*z[0]*z[0] + 
   2*c45*x[1]*y[0]*z[0]*z[0] + c36*x[0]*y[1]*z[0]*z[0] + 2*c54*x[0]*y[1]*z[0]*z[0] + 
   c32*y[0]*y[1]*z[0]*z[0] + 2*c44*y[0]*y[1]*z[0]*z[0] + c35*x[1]*z[0]*z[0]*z[0] + 
   c34*y[1]*z[0]*z[0]*z[0] + c15*x[0]*x[0]*x[0]*z[1] + c14*x[0]*x[0]*y[0]*z[1] + 
   2*c65*x[0]*x[0]*y[0]*z[1] + c25*x[0]*y[0]*y[0]*z[1] + 2*c64*x[0]*y[0]*y[0]*z[1] + 
   c24*y[0]*y[0]*y[0]*z[1] + c13*x[0]*x[0]*z[0]*z[1] + 2*c55*x[0]*x[0]*z[0]*z[1] + 2*c45*x[0]*y[0]*z[0]*z[1] + 
   2*c54*x[0]*y[0]*z[0]*z[1] + 2*c63*x[0]*y[0]*z[0]*z[1] + c23*y[0]*y[0]*z[0]*z[1] + 2*c44*y[0]*y[0]*z[0]*z[1] + 
   c35*x[0]*z[0]*z[0]*z[1] + 2*c53*x[0]*z[0]*z[0]*z[1] + c34*y[0]*z[0]*z[0]*z[1] + 
   2*c43*y[0]*z[0]*z[0]*z[1] + c33*z[0]*z[0]*z[0]*z[1];
   
		mod(1,0) += c11*x[0]*x[0]*x[1]*x[1] + 2*c16*x[0]*x[1]*x[1]*y[0] + c12*x[1]*x[1]*y[0]*y[0] + 
   2*c61*x[0]*x[0]*x[1]*y[1] + 4*c66*x[0]*x[1]*y[0]*y[1] + 2*c62*x[1]*y[0]*y[0]*y[1] + 
   c21*x[0]*x[0]*y[1]*y[1] + 2*c26*x[0]*y[0]*y[1]*y[1] + c22*y[0]*y[0]*y[1]*y[1] + 
   2*c15*x[0]*x[1]*x[1]*z[0] + 2*c14*x[1]*x[1]*y[0]*z[0] + 4*c65*x[0]*x[1]*y[1]*z[0] + 
   4*c64*x[1]*y[0]*y[1]*z[0] + 2*c25*x[0]*y[1]*y[1]*z[0] + 2*c24*y[0]*y[1]*y[1]*z[0] + 
   c13*x[1]*x[1]*z[0]*z[0] + 2*c63*x[1]*y[1]*z[0]*z[0] + c23*y[1]*y[1]*z[0]*z[0] + 
   2*c51*x[0]*x[0]*x[1]*z[1] + 4*c56*x[0]*x[1]*y[0]*z[1] + 2*c52*x[1]*y[0]*y[0]*z[1] + 
   2*c41*x[0]*x[0]*y[1]*z[1] + 4*c46*x[0]*y[0]*y[1]*z[1] + 2*c42*y[0]*y[0]*y[1]*z[1] + 
   4*c55*x[0]*x[1]*z[0]*z[1] + 4*c54*x[1]*y[0]*z[0]*z[1] + 4*c45*x[0]*y[1]*z[0]*z[1] + 4*c44*y[0]*y[1]*z[0]*z[1] + 
   2*c53*x[1]*z[0]*z[0]*z[1] + 2*c43*y[1]*z[0]*z[0]*z[1] + c31*x[0]*x[0]*z[1]*z[1] + 
   2*c36*x[0]*y[0]*z[1]*z[1] + c32*y[0]*y[0]*z[1]*z[1] + 2*c35*x[0]*z[0]*z[1]*z[1] + 
   2*c34*y[0]*z[0]*z[1]*z[1] + c33*z[0]*z[0]*z[1]*z[1];
   
		mod(1,2) += c11*x[1]*x[1]*x[2]*x[2] + 2*c61*x[1]*x[2]*x[2]*y[1] + c21*x[2]*x[2]*y[1]*y[1] + 
   2*c16*x[1]*x[1]*x[2]*y[2] + 4*c66*x[1]*x[2]*y[1]*y[2] + 2*c26*x[2]*y[1]*y[1]*y[2] + 
   c12*x[1]*x[1]*y[2]*y[2] + 2*c62*x[1]*y[1]*y[2]*y[2] + c22*y[1]*y[1]*y[2]*y[2] + 
   2*c51*x[1]*x[2]*x[2]*z[1] + 2*c41*x[2]*x[2]*y[1]*z[1] + 4*c56*x[1]*x[2]*y[2]*z[1] + 
   4*c46*x[2]*y[1]*y[2]*z[1] + 2*c52*x[1]*y[2]*y[2]*z[1] + 2*c42*y[1]*y[2]*y[2]*z[1] + 
   c31*x[2]*x[2]*z[1]*z[1] + 2*c36*x[2]*y[2]*z[1]*z[1] + c32*y[2]*y[2]*z[1]*z[1] + 
   2*c15*x[1]*x[1]*x[2]*z[2] + 4*c65*x[1]*x[2]*y[1]*z[2] + 2*c25*x[2]*y[1]*y[1]*z[2] + 
   2*c14*x[1]*x[1]*y[2]*z[2] + 4*c64*x[1]*y[1]*y[2]*z[2] + 2*c24*y[1]*y[1]*y[2]*z[2] + 
   4*c55*x[1]*x[2]*z[1]*z[2] + 4*c45*x[2]*y[1]*z[1]*z[2] + 4*c54*x[1]*y[2]*z[1]*z[2] + 4*c44*y[1]*y[2]*z[1]*z[2] + 
   2*c35*x[2]*z[1]*z[1]*z[2] + 2*c34*y[2]*z[1]*z[1]*z[2] + c13*x[1]*x[1]*z[2]*z[2] + 
   2*c63*x[1]*y[1]*z[2]*z[2] + c23*y[1]*y[1]*z[2]*z[2] + 2*c53*x[1]*z[1]*z[2]*z[2] + 
   2*c43*y[1]*z[1]*z[2]*z[2] + c33*z[1]*z[1]*z[2]*z[2];
   
		mod(1,3) += c11*x[1]*x[1]*x[1]*x[2] + c16*x[1]*x[1]*x[2]*y[1] + 2*c61*x[1]*x[1]*x[2]*y[1] + 
   c21*x[1]*x[2]*y[1]*y[1] + 2*c66*x[1]*x[2]*y[1]*y[1] + c26*x[2]*y[1]*y[1]*y[1] + 
   c16*x[1]*x[1]*x[1]*y[2] + c12*x[1]*x[1]*y[1]*y[2] + 2*c66*x[1]*x[1]*y[1]*y[2] + 
   c26*x[1]*y[1]*y[1]*y[2] + 2*c62*x[1]*y[1]*y[1]*y[2] + c22*y[1]*y[1]*y[1]*y[2] + 
   c15*x[1]*x[1]*x[2]*z[1] + 2*c51*x[1]*x[1]*x[2]*z[1] + 2*c41*x[1]*x[2]*y[1]*z[1] + 2*c56*x[1]*x[2]*y[1]*z[1] + 
   2*c65*x[1]*x[2]*y[1]*z[1] + c25*x[2]*y[1]*y[1]*z[1] + 2*c46*x[2]*y[1]*y[1]*z[1] + 
   c14*x[1]*x[1]*y[2]*z[1] + 2*c56*x[1]*x[1]*y[2]*z[1] + 2*c46*x[1]*y[1]*y[2]*z[1] + 2*c52*x[1]*y[1]*y[2]*z[1] + 
   2*c64*x[1]*y[1]*y[2]*z[1] + c24*y[1]*y[1]*y[2]*z[1] + 2*c42*y[1]*y[1]*y[2]*z[1] + 
   c31*x[1]*x[2]*z[1]*z[1] + 2*c55*x[1]*x[2]*z[1]*z[1] + c36*x[2]*y[1]*z[1]*z[1] + 
   2*c45*x[2]*y[1]*z[1]*z[1] + c36*x[1]*y[2]*z[1]*z[1] + 2*c54*x[1]*y[2]*z[1]*z[1] + 
   c32*y[1]*y[2]*z[1]*z[1] + 2*c44*y[1]*y[2]*z[1]*z[1] + c35*x[2]*z[1]*z[1]*z[1] + 
   c34*y[2]*z[1]*z[1]*z[1] + c15*x[1]*x[1]*x[1]*z[2] + c14*x[1]*x[1]*y[1]*z[2] + 
   2*c65*x[1]*x[1]*y[1]*z[2] + c25*x[1]*y[1]*y[1]*z[2] + 2*c64*x[1]*y[1]*y[1]*z[2] + 
   c24*y[1]*y[1]*y[1]*z[2] + c13*x[1]*x[1]*z[1]*z[2] + 2*c55*x[1]*x[1]*z[1]*z[2] + 2*c45*x[1]*y[1]*z[1]*z[2] + 
   2*c54*x[1]*y[1]*z[1]*z[2] + 2*c63*x[1]*y[1]*z[1]*z[2] + c23*y[1]*y[1]*z[1]*z[2] + 2*c44*y[1]*y[1]*z[1]*z[2] + 
   c35*x[1]*z[1]*z[1]*z[2] + 2*c53*x[1]*z[1]*z[1]*z[2] + c34*y[1]*z[1]*z[1]*z[2] + 
   2*c43*y[1]*z[1]*z[1]*z[2] + c33*z[1]*z[1]*z[1]*z[2];
   
		mod(1,4) += c11*x[0]*x[1]*x[1]*x[2] + c16*x[1]*x[1]*x[2]*y[0] + 2*c61*x[0]*x[1]*x[2]*y[1] + 2*c66*x[1]*x[2]*y[0]*y[1] + 
   c21*x[0]*x[2]*y[1]*y[1] + c26*x[2]*y[0]*y[1]*y[1] + c16*x[0]*x[1]*x[1]*y[2] + 
   c12*x[1]*x[1]*y[0]*y[2] + 2*c66*x[0]*x[1]*y[1]*y[2] + 2*c62*x[1]*y[0]*y[1]*y[2] + c26*x[0]*y[1]*y[1]*y[2] + 
   c22*y[0]*y[1]*y[1]*y[2] + c15*x[1]*x[1]*x[2]*z[0] + 2*c65*x[1]*x[2]*y[1]*z[0] + 
   c25*x[2]*y[1]*y[1]*z[0] + c14*x[1]*x[1]*y[2]*z[0] + 2*c64*x[1]*y[1]*y[2]*z[0] + 
   c24*y[1]*y[1]*y[2]*z[0] + 2*c51*x[0]*x[1]*x[2]*z[1] + 2*c56*x[1]*x[2]*y[0]*z[1] + 2*c41*x[0]*x[2]*y[1]*z[1] + 
   2*c46*x[2]*y[0]*y[1]*z[1] + 2*c56*x[0]*x[1]*y[2]*z[1] + 2*c52*x[1]*y[0]*y[2]*z[1] + 2*c46*x[0]*y[1]*y[2]*z[1] + 
   2*c42*y[0]*y[1]*y[2]*z[1] + 2*c55*x[1]*x[2]*z[0]*z[1] + 2*c45*x[2]*y[1]*z[0]*z[1] + 2*c54*x[1]*y[2]*z[0]*z[1] + 
   2*c44*y[1]*y[2]*z[0]*z[1] + c31*x[0]*x[2]*z[1]*z[1] + c36*x[2]*y[0]*z[1]*z[1] + 
   c36*x[0]*y[2]*z[1]*z[1] + c32*y[0]*y[2]*z[1]*z[1] + c35*x[2]*z[0]*z[1]*z[1] + 
   c34*y[2]*z[0]*z[1]*z[1] + c15*x[0]*x[1]*x[1]*z[2] + c14*x[1]*x[1]*y[0]*z[2] + 
   2*c65*x[0]*x[1]*y[1]*z[2] + 2*c64*x[1]*y[0]*y[1]*z[2] + c25*x[0]*y[1]*y[1]*z[2] + c24*y[0]*y[1]*y[1]*z[2] + 
   c13*x[1]*x[1]*z[0]*z[2] + 2*c63*x[1]*y[1]*z[0]*z[2] + c23*y[1]*y[1]*z[0]*z[2] + 2*c55*x[0]*x[1]*z[1]*z[2] + 
   2*c54*x[1]*y[0]*z[1]*z[2] + 2*c45*x[0]*y[1]*z[1]*z[2] + 2*c44*y[0]*y[1]*z[1]*z[2] + 2*c53*x[1]*z[0]*z[1]*z[2] + 
   2*c43*y[1]*z[0]*z[1]*z[2] + c35*x[0]*z[1]*z[1]*z[2] + c34*y[0]*z[1]*z[1]*z[2] + c33*z[0]*z[1]*z[1]*z[2];
   
		mod(1,5) += c11*x[0]*x[1]*x[1]*x[1] + c16*x[1]*x[1]*x[1]*y[0] + c16*x[0]*x[1]*x[1]*y[1] + 2*c61*x[0]*x[1]*x[1]*y[1] + 
   c12*x[1]*x[1]*y[0]*y[1] + 2*c66*x[1]*x[1]*y[0]*y[1] + c21*x[0]*x[1]*y[1]*y[1] + 
   2*c66*x[0]*x[1]*y[1]*y[1] + c26*x[1]*y[0]*y[1]*y[1] + 2*c62*x[1]*y[0]*y[1]*y[1] + 
   c26*x[0]*y[1]*y[1]*y[1] + c22*y[0]*y[1]*y[1]*y[1] + c15*x[1]*x[1]*x[1]*z[0] + c14*x[1]*x[1]*y[1]*z[0] + 
   2*c65*x[1]*x[1]*y[1]*z[0] + c25*x[1]*y[1]*y[1]*z[0] + 2*c64*x[1]*y[1]*y[1]*z[0] + 
   c24*y[1]*y[1]*y[1]*z[0] + c15*x[0]*x[1]*x[1]*z[1] + 2*c51*x[0]*x[1]*x[1]*z[1] + 
   c14*x[1]*x[1]*y[0]*z[1] + 2*c56*x[1]*x[1]*y[0]*z[1] + 2*c41*x[0]*x[1]*y[1]*z[1] + 2*c56*x[0]*x[1]*y[1]*z[1] + 
   2*c65*x[0]*x[1]*y[1]*z[1] + 2*c46*x[1]*y[0]*y[1]*z[1] + 2*c52*x[1]*y[0]*y[1]*z[1] + 2*c64*x[1]*y[0]*y[1]*z[1] + 
   c25*x[0]*y[1]*y[1]*z[1] + 2*c46*x[0]*y[1]*y[1]*z[1] + c24*y[0]*y[1]*y[1]*z[1] + 
   2*c42*y[0]*y[1]*y[1]*z[1] + c13*x[1]*x[1]*z[0]*z[1] + 2*c55*x[1]*x[1]*z[0]*z[1] + 
   2*c45*x[1]*y[1]*z[0]*z[1] + 2*c54*x[1]*y[1]*z[0]*z[1] + 2*c63*x[1]*y[1]*z[0]*z[1] + c23*y[1]*y[1]*z[0]*z[1] + 
   2*c44*y[1]*y[1]*z[0]*z[1] + c31*x[0]*x[1]*z[1]*z[1] + 2*c55*x[0]*x[1]*z[1]*z[1] + 
   c36*x[1]*y[0]*z[1]*z[1] + 2*c54*x[1]*y[0]*z[1]*z[1] + c36*x[0]*y[1]*z[1]*z[1] + 
   2*c45*x[0]*y[1]*z[1]*z[1] + c32*y[0]*y[1]*z[1]*z[1] + 2*c44*y[0]*y[1]*z[1]*z[1] + 
   c35*x[1]*z[0]*z[1]*z[1] + 2*c53*x[1]*z[0]*z[1]*z[1] + c34*y[1]*z[0]*z[1]*z[1] + 
   2*c43*y[1]*z[0]*z[1]*z[1] + c35*x[0]*z[1]*z[1]*z[1] + c34*y[0]*z[1]*z[1]*z[1] + c33*z[0]*z[1]*z[1]*z[1];
   
		mod(2,0) += c11*x[0]*x[0]*x[2]*x[2] + 2*c16*x[0]*x[2]*x[2]*y[0] + c12*x[2]*x[2]*y[0]*y[0] + 
   2*c61*x[0]*x[0]*x[2]*y[2] + 4*c66*x[0]*x[2]*y[0]*y[2] + 2*c62*x[2]*y[0]*y[0]*y[2] + 
   c21*x[0]*x[0]*y[2]*y[2] + 2*c26*x[0]*y[0]*y[2]*y[2] + c22*y[0]*y[0]*y[2]*y[2] + 
   2*c15*x[0]*x[2]*x[2]*z[0] + 2*c14*x[2]*x[2]*y[0]*z[0] + 4*c65*x[0]*x[2]*y[2]*z[0] + 
   4*c64*x[2]*y[0]*y[2]*z[0] + 2*c25*x[0]*y[2]*y[2]*z[0] + 2*c24*y[0]*y[2]*y[2]*z[0] + 
   c13*x[2]*x[2]*z[0]*z[0] + 2*c63*x[2]*y[2]*z[0]*z[0] + c23*y[2]*y[2]*z[0]*z[0] + 
   2*c51*x[0]*x[0]*x[2]*z[2] + 4*c56*x[0]*x[2]*y[0]*z[2] + 2*c52*x[2]*y[0]*y[0]*z[2] + 
   2*c41*x[0]*x[0]*y[2]*z[2] + 4*c46*x[0]*y[0]*y[2]*z[2] + 2*c42*y[0]*y[0]*y[2]*z[2] + 
   4*c55*x[0]*x[2]*z[0]*z[2] + 4*c54*x[2]*y[0]*z[0]*z[2] + 4*c45*x[0]*y[2]*z[0]*z[2] + 4*c44*y[0]*y[2]*z[0]*z[2] + 
   2*c53*x[2]*z[0]*z[0]*z[2] + 2*c43*y[2]*z[0]*z[0]*z[2] + c31*x[0]*x[0]*z[2]*z[2] + 
   2*c36*x[0]*y[0]*z[2]*z[2] + c32*y[0]*y[0]*z[2]*z[2] + 2*c35*x[0]*z[0]*z[2]*z[2] + 
   2*c34*y[0]*z[0]*z[2]*z[2] + c33*z[0]*z[0]*z[2]*z[2];
   
		mod(2,1) += c11*x[1]*x[1]*x[2]*x[2] + 2*c16*x[1]*x[2]*x[2]*y[1] + c12*x[2]*x[2]*y[1]*y[1] + 
   2*c61*x[1]*x[1]*x[2]*y[2] + 4*c66*x[1]*x[2]*y[1]*y[2] + 2*c62*x[2]*y[1]*y[1]*y[2] + 
   c21*x[1]*x[1]*y[2]*y[2] + 2*c26*x[1]*y[1]*y[2]*y[2] + c22*y[1]*y[1]*y[2]*y[2] + 
   2*c15*x[1]*x[2]*x[2]*z[1] + 2*c14*x[2]*x[2]*y[1]*z[1] + 4*c65*x[1]*x[2]*y[2]*z[1] + 
   4*c64*x[2]*y[1]*y[2]*z[1] + 2*c25*x[1]*y[2]*y[2]*z[1] + 2*c24*y[1]*y[2]*y[2]*z[1] + 
   c13*x[2]*x[2]*z[1]*z[1] + 2*c63*x[2]*y[2]*z[1]*z[1] + c23*y[2]*y[2]*z[1]*z[1] + 
   2*c51*x[1]*x[1]*x[2]*z[2] + 4*c56*x[1]*x[2]*y[1]*z[2] + 2*c52*x[2]*y[1]*y[1]*z[2] + 
   2*c41*x[1]*x[1]*y[2]*z[2] + 4*c46*x[1]*y[1]*y[2]*z[2] + 2*c42*y[1]*y[1]*y[2]*z[2] + 
   4*c55*x[1]*x[2]*z[1]*z[2] + 4*c54*x[2]*y[1]*z[1]*z[2] + 4*c45*x[1]*y[2]*z[1]*z[2] + 4*c44*y[1]*y[2]*z[1]*z[2] + 
   2*c53*x[2]*z[1]*z[1]*z[2] + 2*c43*y[2]*z[1]*z[1]*z[2] + c31*x[1]*x[1]*z[2]*z[2] + 
   2*c36*x[1]*y[1]*z[2]*z[2] + c32*y[1]*y[1]*z[2]*z[2] + 2*c35*x[1]*z[1]*z[2]*z[2] + 
   2*c34*y[1]*z[1]*z[2]*z[2] + c33*z[1]*z[1]*z[2]*z[2];
   
		mod(2,3) += c11*x[1]*x[2]*x[2]*x[2] + c16*x[2]*x[2]*x[2]*y[1] + c16*x[1]*x[2]*x[2]*y[2] + 2*c61*x[1]*x[2]*x[2]*y[2] + 
   c12*x[2]*x[2]*y[1]*y[2] + 2*c66*x[2]*x[2]*y[1]*y[2] + c21*x[1]*x[2]*y[2]*y[2] + 
   2*c66*x[1]*x[2]*y[2]*y[2] + c26*x[2]*y[1]*y[2]*y[2] + 2*c62*x[2]*y[1]*y[2]*y[2] + 
   c26*x[1]*y[2]*y[2]*y[2] + c22*y[1]*y[2]*y[2]*y[2] + c15*x[2]*x[2]*x[2]*z[1] + c14*x[2]*x[2]*y[2]*z[1] + 
   2*c65*x[2]*x[2]*y[2]*z[1] + c25*x[2]*y[2]*y[2]*z[1] + 2*c64*x[2]*y[2]*y[2]*z[1] + 
   c24*y[2]*y[2]*y[2]*z[1] + c15*x[1]*x[2]*x[2]*z[2] + 2*c51*x[1]*x[2]*x[2]*z[2] + 
   c14*x[2]*x[2]*y[1]*z[2] + 2*c56*x[2]*x[2]*y[1]*z[2] + 2*c41*x[1]*x[2]*y[2]*z[2] + 2*c56*x[1]*x[2]*y[2]*z[2] + 
   2*c65*x[1]*x[2]*y[2]*z[2] + 2*c46*x[2]*y[1]*y[2]*z[2] + 2*c52*x[2]*y[1]*y[2]*z[2] + 2*c64*x[2]*y[1]*y[2]*z[2] + 
   c25*x[1]*y[2]*y[2]*z[2] + 2*c46*x[1]*y[2]*y[2]*z[2] + c24*y[1]*y[2]*y[2]*z[2] + 
   2*c42*y[1]*y[2]*y[2]*z[2] + c13*x[2]*x[2]*z[1]*z[2] + 2*c55*x[2]*x[2]*z[1]*z[2] + 
   2*c45*x[2]*y[2]*z[1]*z[2] + 2*c54*x[2]*y[2]*z[1]*z[2] + 2*c63*x[2]*y[2]*z[1]*z[2] + c23*y[2]*y[2]*z[1]*z[2] + 
   2*c44*y[2]*y[2]*z[1]*z[2] + c31*x[1]*x[2]*z[2]*z[2] + 2*c55*x[1]*x[2]*z[2]*z[2] + 
   c36*x[2]*y[1]*z[2]*z[2] + 2*c54*x[2]*y[1]*z[2]*z[2] + c36*x[1]*y[2]*z[2]*z[2] + 
   2*c45*x[1]*y[2]*z[2]*z[2] + c32*y[1]*y[2]*z[2]*z[2] + 2*c44*y[1]*y[2]*z[2]*z[2] + 
   c35*x[2]*z[1]*z[2]*z[2] + 2*c53*x[2]*z[1]*z[2]*z[2] + c34*y[2]*z[1]*z[2]*z[2] + 
   2*c43*y[2]*z[1]*z[2]*z[2] + c35*x[1]*z[2]*z[2]*z[2] + c34*y[1]*z[2]*z[2]*z[2] + c33*z[1]*z[2]*z[2]*z[2];
   
		mod(2,4) += c11*x[0]*x[2]*x[2]*x[2] + c16*x[2]*x[2]*x[2]*y[0] + c16*x[0]*x[2]*x[2]*y[2] + 2*c61*x[0]*x[2]*x[2]*y[2] + 
   c12*x[2]*x[2]*y[0]*y[2] + 2*c66*x[2]*x[2]*y[0]*y[2] + c21*x[0]*x[2]*y[2]*y[2] + 
   2*c66*x[0]*x[2]*y[2]*y[2] + c26*x[2]*y[0]*y[2]*y[2] + 2*c62*x[2]*y[0]*y[2]*y[2] + 
   c26*x[0]*y[2]*y[2]*y[2] + c22*y[0]*y[2]*y[2]*y[2] + c15*x[2]*x[2]*x[2]*z[0] + c14*x[2]*x[2]*y[2]*z[0] + 
   2*c65*x[2]*x[2]*y[2]*z[0] + c25*x[2]*y[2]*y[2]*z[0] + 2*c64*x[2]*y[2]*y[2]*z[0] + 
   c24*y[2]*y[2]*y[2]*z[0] + c15*x[0]*x[2]*x[2]*z[2] + 2*c51*x[0]*x[2]*x[2]*z[2] + 
   c14*x[2]*x[2]*y[0]*z[2] + 2*c56*x[2]*x[2]*y[0]*z[2] + 2*c41*x[0]*x[2]*y[2]*z[2] + 2*c56*x[0]*x[2]*y[2]*z[2] + 
   2*c65*x[0]*x[2]*y[2]*z[2] + 2*c46*x[2]*y[0]*y[2]*z[2] + 2*c52*x[2]*y[0]*y[2]*z[2] + 2*c64*x[2]*y[0]*y[2]*z[2] + 
   c25*x[0]*y[2]*y[2]*z[2] + 2*c46*x[0]*y[2]*y[2]*z[2] + c24*y[0]*y[2]*y[2]*z[2] + 
   2*c42*y[0]*y[2]*y[2]*z[2] + c13*x[2]*x[2]*z[0]*z[2] + 2*c55*x[2]*x[2]*z[0]*z[2] + 
   2*c45*x[2]*y[2]*z[0]*z[2] + 2*c54*x[2]*y[2]*z[0]*z[2] + 2*c63*x[2]*y[2]*z[0]*z[2] + c23*y[2]*y[2]*z[0]*z[2] + 
   2*c44*y[2]*y[2]*z[0]*z[2] + c31*x[0]*x[2]*z[2]*z[2] + 2*c55*x[0]*x[2]*z[2]*z[2] + 
   c36*x[2]*y[0]*z[2]*z[2] + 2*c54*x[2]*y[0]*z[2]*z[2] + c36*x[0]*y[2]*z[2]*z[2] + 
   2*c45*x[0]*y[2]*z[2]*z[2] + c32*y[0]*y[2]*z[2]*z[2] + 2*c44*y[0]*y[2]*z[2]*z[2] + 
   c35*x[2]*z[0]*z[2]*z[2] + 2*c53*x[2]*z[0]*z[2]*z[2] + c34*y[2]*z[0]*z[2]*z[2] + 
   2*c43*y[2]*z[0]*z[2]*z[2] + c35*x[0]*z[2]*z[2]*z[2] + c34*y[0]*z[2]*z[2]*z[2] + c33*z[0]*z[2]*z[2]*z[2];
   
		mod(2,5) += c11*x[0]*x[1]*x[2]*x[2] + c16*x[1]*x[2]*x[2]*y[0] + c16*x[0]*x[2]*x[2]*y[1] + 
   c12*x[2]*x[2]*y[0]*y[1] + 2*c61*x[0]*x[1]*x[2]*y[2] + 2*c66*x[1]*x[2]*y[0]*y[2] + 2*c66*x[0]*x[2]*y[1]*y[2] + 
   2*c62*x[2]*y[0]*y[1]*y[2] + c21*x[0]*x[1]*y[2]*y[2] + c26*x[1]*y[0]*y[2]*y[2] + 
   c26*x[0]*y[1]*y[2]*y[2] + c22*y[0]*y[1]*y[2]*y[2] + c15*x[1]*x[2]*x[2]*z[0] + 
   c14*x[2]*x[2]*y[1]*z[0] + 2*c65*x[1]*x[2]*y[2]*z[0] + 2*c64*x[2]*y[1]*y[2]*z[0] + c25*x[1]*y[2]*y[2]*z[0] + 
   c24*y[1]*y[2]*y[2]*z[0] + c15*x[0]*x[2]*x[2]*z[1] + c14*x[2]*x[2]*y[0]*z[1] + 
   2*c65*x[0]*x[2]*y[2]*z[1] + 2*c64*x[2]*y[0]*y[2]*z[1] + c25*x[0]*y[2]*y[2]*z[1] + c24*y[0]*y[2]*y[2]*z[1] + 
   c13*x[2]*x[2]*z[0]*z[1] + 2*c63*x[2]*y[2]*z[0]*z[1] + c23*y[2]*y[2]*z[0]*z[1] + 2*c51*x[0]*x[1]*x[2]*z[2] + 
   2*c56*x[1]*x[2]*y[0]*z[2] + 2*c56*x[0]*x[2]*y[1]*z[2] + 2*c52*x[2]*y[0]*y[1]*z[2] + 2*c41*x[0]*x[1]*y[2]*z[2] + 
   2*c46*x[1]*y[0]*y[2]*z[2] + 2*c46*x[0]*y[1]*y[2]*z[2] + 2*c42*y[0]*y[1]*y[2]*z[2] + 2*c55*x[1]*x[2]*z[0]*z[2] + 
   2*c54*x[2]*y[1]*z[0]*z[2] + 2*c45*x[1]*y[2]*z[0]*z[2] + 2*c44*y[1]*y[2]*z[0]*z[2] + 2*c55*x[0]*x[2]*z[1]*z[2] + 
   2*c54*x[2]*y[0]*z[1]*z[2] + 2*c45*x[0]*y[2]*z[1]*z[2] + 2*c44*y[0]*y[2]*z[1]*z[2] + 2*c53*x[2]*z[0]*z[1]*z[2] + 
   2*c43*y[2]*z[0]*z[1]*z[2] + c31*x[0]*x[1]*z[2]*z[2] + c36*x[1]*y[0]*z[2]*z[2] + 
   c36*x[0]*y[1]*z[2]*z[2] + c32*y[0]*y[1]*z[2]*z[2] + c35*x[1]*z[0]*z[2]*z[2] + 
   c34*y[1]*z[0]*z[2]*z[2] + c35*x[0]*z[1]*z[2]*z[2] + c34*y[0]*z[1]*z[2]*z[2] + 
   c33*z[0]*z[1]*z[2]*z[2];
   
		mod(3,0) += c11*x[0]*x[0]*x[1]*x[2] + 2*c16*x[0]*x[1]*x[2]*y[0] + c12*x[1]*x[2]*y[0]*y[0] + c61*x[0]*x[0]*x[2]*y[1] + 
   2*c66*x[0]*x[2]*y[0]*y[1] + c62*x[2]*y[0]*y[0]*y[1] + c61*x[0]*x[0]*x[1]*y[2] + 2*c66*x[0]*x[1]*y[0]*y[2] + 
   c62*x[1]*y[0]*y[0]*y[2] + c21*x[0]*x[0]*y[1]*y[2] + 2*c26*x[0]*y[0]*y[1]*y[2] + 
   c22*y[0]*y[0]*y[1]*y[2] + 2*c15*x[0]*x[1]*x[2]*z[0] + 2*c14*x[1]*x[2]*y[0]*z[0] + 2*c65*x[0]*x[2]*y[1]*z[0] + 
   2*c64*x[2]*y[0]*y[1]*z[0] + 2*c65*x[0]*x[1]*y[2]*z[0] + 2*c64*x[1]*y[0]*y[2]*z[0] + 2*c25*x[0]*y[1]*y[2]*z[0] + 
   2*c24*y[0]*y[1]*y[2]*z[0] + c13*x[1]*x[2]*z[0]*z[0] + c63*x[2]*y[1]*z[0]*z[0] + 
   c63*x[1]*y[2]*z[0]*z[0] + c23*y[1]*y[2]*z[0]*z[0] + c51*x[0]*x[0]*x[2]*z[1] + 
   2*c56*x[0]*x[2]*y[0]*z[1] + c52*x[2]*y[0]*y[0]*z[1] + c41*x[0]*x[0]*y[2]*z[1] + 2*c46*x[0]*y[0]*y[2]*z[1] + 
   c42*y[0]*y[0]*y[2]*z[1] + 2*c55*x[0]*x[2]*z[0]*z[1] + 2*c54*x[2]*y[0]*z[0]*z[1] + 2*c45*x[0]*y[2]*z[0]*z[1] + 
   2*c44*y[0]*y[2]*z[0]*z[1] + c53*x[2]*z[0]*z[0]*z[1] + c43*y[2]*z[0]*z[0]*z[1] + 
   c51*x[0]*x[0]*x[1]*z[2] + 2*c56*x[0]*x[1]*y[0]*z[2] + c52*x[1]*y[0]*y[0]*z[2] + 
   c41*x[0]*x[0]*y[1]*z[2] + 2*c46*x[0]*y[0]*y[1]*z[2] + c42*y[0]*y[0]*y[1]*z[2] + 2*c55*x[0]*x[1]*z[0]*z[2] + 
   2*c54*x[1]*y[0]*z[0]*z[2] + 2*c45*x[0]*y[1]*z[0]*z[2] + 2*c44*y[0]*y[1]*z[0]*z[2] + c53*x[1]*z[0]*z[0]*z[2] + 
   c43*y[1]*z[0]*z[0]*z[2] + c31*x[0]*x[0]*z[1]*z[2] + 2*c36*x[0]*y[0]*z[1]*z[2] + 
   c32*y[0]*y[0]*z[1]*z[2] + 2*c35*x[0]*z[0]*z[1]*z[2] + 2*c34*y[0]*z[0]*z[1]*z[2] + c33*z[0]*z[0]*z[1]*z[2];
   
		mod(3,1) += c11*x[1]*x[1]*x[1]*x[2] + 2*c16*x[1]*x[1]*x[2]*y[1] + c61*x[1]*x[1]*x[2]*y[1] + 
   c12*x[1]*x[2]*y[1]*y[1] + 2*c66*x[1]*x[2]*y[1]*y[1] + c62*x[2]*y[1]*y[1]*y[1] + 
   c61*x[1]*x[1]*x[1]*y[2] + c21*x[1]*x[1]*y[1]*y[2] + 2*c66*x[1]*x[1]*y[1]*y[2] + 
   2*c26*x[1]*y[1]*y[1]*y[2] + c62*x[1]*y[1]*y[1]*y[2] + c22*y[1]*y[1]*y[1]*y[2] + 
   2*c15*x[1]*x[1]*x[2]*z[1] + c51*x[1]*x[1]*x[2]*z[1] + 2*c14*x[1]*x[2]*y[1]*z[1] + 2*c56*x[1]*x[2]*y[1]*z[1] + 
   2*c65*x[1]*x[2]*y[1]*z[1] + c52*x[2]*y[1]*y[1]*z[1] + 2*c64*x[2]*y[1]*y[1]*z[1] + 
   c41*x[1]*x[1]*y[2]*z[1] + 2*c65*x[1]*x[1]*y[2]*z[1] + 2*c25*x[1]*y[1]*y[2]*z[1] + 2*c46*x[1]*y[1]*y[2]*z[1] + 
   2*c64*x[1]*y[1]*y[2]*z[1] + 2*c24*y[1]*y[1]*y[2]*z[1] + c42*y[1]*y[1]*y[2]*z[1] + 
   c13*x[1]*x[2]*z[1]*z[1] + 2*c55*x[1]*x[2]*z[1]*z[1] + 2*c54*x[2]*y[1]*z[1]*z[1] + 
   c63*x[2]*y[1]*z[1]*z[1] + 2*c45*x[1]*y[2]*z[1]*z[1] + c63*x[1]*y[2]*z[1]*z[1] + 
   c23*y[1]*y[2]*z[1]*z[1] + 2*c44*y[1]*y[2]*z[1]*z[1] + c53*x[2]*z[1]*z[1]*z[1] + 
   c43*y[2]*z[1]*z[1]*z[1] + c51*x[1]*x[1]*x[1]*z[2] + c41*x[1]*x[1]*y[1]*z[2] + 
   2*c56*x[1]*x[1]*y[1]*z[2] + 2*c46*x[1]*y[1]*y[1]*z[2] + c52*x[1]*y[1]*y[1]*z[2] + 
   c42*y[1]*y[1]*y[1]*z[2] + c31*x[1]*x[1]*z[1]*z[2] + 2*c55*x[1]*x[1]*z[1]*z[2] + 2*c36*x[1]*y[1]*z[1]*z[2] + 
   2*c45*x[1]*y[1]*z[1]*z[2] + 2*c54*x[1]*y[1]*z[1]*z[2] + c32*y[1]*y[1]*z[1]*z[2] + 2*c44*y[1]*y[1]*z[1]*z[2] + 
   2*c35*x[1]*z[1]*z[1]*z[2] + c53*x[1]*z[1]*z[1]*z[2] + 2*c34*y[1]*z[1]*z[1]*z[2] + 
   c43*y[1]*z[1]*z[1]*z[2] + c33*z[1]*z[1]*z[1]*z[2];
   
		mod(3,2) += c11*x[1]*x[2]*x[2]*x[2] + c61*x[2]*x[2]*x[2]*y[1] + 2*c16*x[1]*x[2]*x[2]*y[2] + c61*x[1]*x[2]*x[2]*y[2] + 
   c21*x[2]*x[2]*y[1]*y[2] + 2*c66*x[2]*x[2]*y[1]*y[2] + c12*x[1]*x[2]*y[2]*y[2] + 
   2*c66*x[1]*x[2]*y[2]*y[2] + 2*c26*x[2]*y[1]*y[2]*y[2] + c62*x[2]*y[1]*y[2]*y[2] + 
   c62*x[1]*y[2]*y[2]*y[2] + c22*y[1]*y[2]*y[2]*y[2] + c51*x[2]*x[2]*x[2]*z[1] + c41*x[2]*x[2]*y[2]*z[1] + 
   2*c56*x[2]*x[2]*y[2]*z[1] + 2*c46*x[2]*y[2]*y[2]*z[1] + c52*x[2]*y[2]*y[2]*z[1] + 
   c42*y[2]*y[2]*y[2]*z[1] + 2*c15*x[1]*x[2]*x[2]*z[2] + c51*x[1]*x[2]*x[2]*z[2] + 
   c41*x[2]*x[2]*y[1]*z[2] + 2*c65*x[2]*x[2]*y[1]*z[2] + 2*c14*x[1]*x[2]*y[2]*z[2] + 2*c56*x[1]*x[2]*y[2]*z[2] + 
   2*c65*x[1]*x[2]*y[2]*z[2] + 2*c25*x[2]*y[1]*y[2]*z[2] + 2*c46*x[2]*y[1]*y[2]*z[2] + 2*c64*x[2]*y[1]*y[2]*z[2] + 
   c52*x[1]*y[2]*y[2]*z[2] + 2*c64*x[1]*y[2]*y[2]*z[2] + 2*c24*y[1]*y[2]*y[2]*z[2] + 
   c42*y[1]*y[2]*y[2]*z[2] + c31*x[2]*x[2]*z[1]*z[2] + 2*c55*x[2]*x[2]*z[1]*z[2] + 
   2*c36*x[2]*y[2]*z[1]*z[2] + 2*c45*x[2]*y[2]*z[1]*z[2] + 2*c54*x[2]*y[2]*z[1]*z[2] + c32*y[2]*y[2]*z[1]*z[2] + 
   2*c44*y[2]*y[2]*z[1]*z[2] + c13*x[1]*x[2]*z[2]*z[2] + 2*c55*x[1]*x[2]*z[2]*z[2] + 
   2*c45*x[2]*y[1]*z[2]*z[2] + c63*x[2]*y[1]*z[2]*z[2] + 2*c54*x[1]*y[2]*z[2]*z[2] + 
   c63*x[1]*y[2]*z[2]*z[2] + c23*y[1]*y[2]*z[2]*z[2] + 2*c44*y[1]*y[2]*z[2]*z[2] + 
   2*c35*x[2]*z[1]*z[2]*z[2] + c53*x[2]*z[1]*z[2]*z[2] + 2*c34*y[2]*z[1]*z[2]*z[2] + 
   c43*y[2]*z[1]*z[2]*z[2] + c53*x[1]*z[2]*z[2]*z[2] + c43*y[1]*z[2]*z[2]*z[2] + c33*z[1]*z[2]*z[2]*z[2];
   
		mod(3,4) += c11*x[0]*x[1]*x[2]*x[2] + c16*x[1]*x[2]*x[2]*y[0] + c61*x[0]*x[2]*x[2]*y[1] + 
   c66*x[2]*x[2]*y[0]*y[1] + c16*x[0]*x[1]*x[2]*y[2] + c61*x[0]*x[1]*x[2]*y[2] + c12*x[1]*x[2]*y[0]*y[2] + 
   c66*x[1]*x[2]*y[0]*y[2] + c21*x[0]*x[2]*y[1]*y[2] + c66*x[0]*x[2]*y[1]*y[2] + c26*x[2]*y[0]*y[1]*y[2] + c62*x[2]*y[0]*y[1]*y[2] + 
   c66*x[0]*x[1]*y[2]*y[2] + c62*x[1]*y[0]*y[2]*y[2] + c26*x[0]*y[1]*y[2]*y[2] + 
   c22*y[0]*y[1]*y[2]*y[2] + c15*x[1]*x[2]*x[2]*z[0] + c65*x[2]*x[2]*y[1]*z[0] + c14*x[1]*x[2]*y[2]*z[0] + 
   c65*x[1]*x[2]*y[2]*z[0] + c25*x[2]*y[1]*y[2]*z[0] + c64*x[2]*y[1]*y[2]*z[0] + c64*x[1]*y[2]*y[2]*z[0] + 
   c24*y[1]*y[2]*y[2]*z[0] + c51*x[0]*x[2]*x[2]*z[1] + c56*x[2]*x[2]*y[0]*z[1] + c41*x[0]*x[2]*y[2]*z[1] + 
   c56*x[0]*x[2]*y[2]*z[1] + c46*x[2]*y[0]*y[2]*z[1] + c52*x[2]*y[0]*y[2]*z[1] + c46*x[0]*y[2]*y[2]*z[1] + 
   c42*y[0]*y[2]*y[2]*z[1] + c55*x[2]*x[2]*z[0]*z[1] + c45*x[2]*y[2]*z[0]*z[1] + c54*x[2]*y[2]*z[0]*z[1] + 
   c44*y[2]*y[2]*z[0]*z[1] + c15*x[0]*x[1]*x[2]*z[2] + c51*x[0]*x[1]*x[2]*z[2] + c14*x[1]*x[2]*y[0]*z[2] + 
   c56*x[1]*x[2]*y[0]*z[2] + c41*x[0]*x[2]*y[1]*z[2] + c65*x[0]*x[2]*y[1]*z[2] + c46*x[2]*y[0]*y[1]*z[2] + c64*x[2]*y[0]*y[1]*z[2] + 
   c56*x[0]*x[1]*y[2]*z[2] + c65*x[0]*x[1]*y[2]*z[2] + c52*x[1]*y[0]*y[2]*z[2] + c64*x[1]*y[0]*y[2]*z[2] + c25*x[0]*y[1]*y[2]*z[2] + 
   c46*x[0]*y[1]*y[2]*z[2] + c24*y[0]*y[1]*y[2]*z[2] + c42*y[0]*y[1]*y[2]*z[2] + c13*x[1]*x[2]*z[0]*z[2] + c55*x[1]*x[2]*z[0]*z[2] + 
   c45*x[2]*y[1]*z[0]*z[2] + c63*x[2]*y[1]*z[0]*z[2] + c54*x[1]*y[2]*z[0]*z[2] + c63*x[1]*y[2]*z[0]*z[2] + c23*y[1]*y[2]*z[0]*z[2] + 
   c44*y[1]*y[2]*z[0]*z[2] + c31*x[0]*x[2]*z[1]*z[2] + c55*x[0]*x[2]*z[1]*z[2] + c36*x[2]*y[0]*z[1]*z[2] + c54*x[2]*y[0]*z[1]*z[2] + 
   c36*x[0]*y[2]*z[1]*z[2] + c45*x[0]*y[2]*z[1]*z[2] + c32*y[0]*y[2]*z[1]*z[2] + c44*y[0]*y[2]*z[1]*z[2] + c35*x[2]*z[0]*z[1]*z[2] + 
   c53*x[2]*z[0]*z[1]*z[2] + c34*y[2]*z[0]*z[1]*z[2] + c43*y[2]*z[0]*z[1]*z[2] + c55*x[0]*x[1]*z[2]*z[2] + 
   c54*x[1]*y[0]*z[2]*z[2] + c45*x[0]*y[1]*z[2]*z[2] + c44*y[0]*y[1]*z[2]*z[2] + 
   c53*x[1]*z[0]*z[2]*z[2] + c43*y[1]*z[0]*z[2]*z[2] + c35*x[0]*z[1]*z[2]*z[2] + 
   c34*y[0]*z[1]*z[2]*z[2] + c33*z[0]*z[1]*z[2]*z[2];
   
		mod(3,5) += c11*x[0]*x[1]*x[1]*x[2] + c16*x[1]*x[1]*x[2]*y[0] + c16*x[0]*x[1]*x[2]*y[1] + c61*x[0]*x[1]*x[2]*y[1] + 
   c12*x[1]*x[2]*y[0]*y[1] + c66*x[1]*x[2]*y[0]*y[1] + c66*x[0]*x[2]*y[1]*y[1] + c62*x[2]*y[0]*y[1]*y[1] + 
   c61*x[0]*x[1]*x[1]*y[2] + c66*x[1]*x[1]*y[0]*y[2] + c21*x[0]*x[1]*y[1]*y[2] + c66*x[0]*x[1]*y[1]*y[2] + 
   c26*x[1]*y[0]*y[1]*y[2] + c62*x[1]*y[0]*y[1]*y[2] + c26*x[0]*y[1]*y[1]*y[2] + c22*y[0]*y[1]*y[1]*y[2] + 
   c15*x[1]*x[1]*x[2]*z[0] + c14*x[1]*x[2]*y[1]*z[0] + c65*x[1]*x[2]*y[1]*z[0] + c64*x[2]*y[1]*y[1]*z[0] + 
   c65*x[1]*x[1]*y[2]*z[0] + c25*x[1]*y[1]*y[2]*z[0] + c64*x[1]*y[1]*y[2]*z[0] + c24*y[1]*y[1]*y[2]*z[0] + 
   c15*x[0]*x[1]*x[2]*z[1] + c51*x[0]*x[1]*x[2]*z[1] + c14*x[1]*x[2]*y[0]*z[1] + c56*x[1]*x[2]*y[0]*z[1] + c56*x[0]*x[2]*y[1]*z[1] + 
   c65*x[0]*x[2]*y[1]*z[1] + c52*x[2]*y[0]*y[1]*z[1] + c64*x[2]*y[0]*y[1]*z[1] + c41*x[0]*x[1]*y[2]*z[1] + c65*x[0]*x[1]*y[2]*z[1] + 
   c46*x[1]*y[0]*y[2]*z[1] + c64*x[1]*y[0]*y[2]*z[1] + c25*x[0]*y[1]*y[2]*z[1] + c46*x[0]*y[1]*y[2]*z[1] + c24*y[0]*y[1]*y[2]*z[1] + 
   c42*y[0]*y[1]*y[2]*z[1] + c13*x[1]*x[2]*z[0]*z[1] + c55*x[1]*x[2]*z[0]*z[1] + c54*x[2]*y[1]*z[0]*z[1] + c63*x[2]*y[1]*z[0]*z[1] + 
   c45*x[1]*y[2]*z[0]*z[1] + c63*x[1]*y[2]*z[0]*z[1] + c23*y[1]*y[2]*z[0]*z[1] + c44*y[1]*y[2]*z[0]*z[1] + 
   c55*x[0]*x[2]*z[1]*z[1] + c54*x[2]*y[0]*z[1]*z[1] + c45*x[0]*y[2]*z[1]*z[1] + 
   c44*y[0]*y[2]*z[1]*z[1] + c53*x[2]*z[0]*z[1]*z[1] + c43*y[2]*z[0]*z[1]*z[1] + 
   c51*x[0]*x[1]*x[1]*z[2] + c56*x[1]*x[1]*y[0]*z[2] + c41*x[0]*x[1]*y[1]*z[2] + c56*x[0]*x[1]*y[1]*z[2] + 
   c46*x[1]*y[0]*y[1]*z[2] + c52*x[1]*y[0]*y[1]*z[2] + c46*x[0]*y[1]*y[1]*z[2] + c42*y[0]*y[1]*y[1]*z[2] + 
   c55*x[1]*x[1]*z[0]*z[2] + c45*x[1]*y[1]*z[0]*z[2] + c54*x[1]*y[1]*z[0]*z[2] + c44*y[1]*y[1]*z[0]*z[2] + 
   c31*x[0]*x[1]*z[1]*z[2] + c55*x[0]*x[1]*z[1]*z[2] + c36*x[1]*y[0]*z[1]*z[2] + c54*x[1]*y[0]*z[1]*z[2] + c36*x[0]*y[1]*z[1]*z[2] + 
   c45*x[0]*y[1]*z[1]*z[2] + c32*y[0]*y[1]*z[1]*z[2] + c44*y[0]*y[1]*z[1]*z[2] + c35*x[1]*z[0]*z[1]*z[2] + c53*x[1]*z[0]*z[1]*z[2] + 
   c34*y[1]*z[0]*z[1]*z[2] + c43*y[1]*z[0]*z[1]*z[2] + c35*x[0]*z[1]*z[1]*z[2] + c34*y[0]*z[1]*z[1]*z[2] + 
   c33*z[0]*z[1]*z[1]*z[2];
   
		mod(4,0) += c11*x[0]*x[0]*x[0]*x[2] + 2*c16*x[0]*x[0]*x[2]*y[0] + c61*x[0]*x[0]*x[2]*y[0] + 
   c12*x[0]*x[2]*y[0]*y[0] + 2*c66*x[0]*x[2]*y[0]*y[0] + c62*x[2]*y[0]*y[0]*y[0] + 
   c61*x[0]*x[0]*x[0]*y[2] + c21*x[0]*x[0]*y[0]*y[2] + 2*c66*x[0]*x[0]*y[0]*y[2] + 
   2*c26*x[0]*y[0]*y[0]*y[2] + c62*x[0]*y[0]*y[0]*y[2] + c22*y[0]*y[0]*y[0]*y[2] + 
   2*c15*x[0]*x[0]*x[2]*z[0] + c51*x[0]*x[0]*x[2]*z[0] + 2*c14*x[0]*x[2]*y[0]*z[0] + 2*c56*x[0]*x[2]*y[0]*z[0] + 
   2*c65*x[0]*x[2]*y[0]*z[0] + c52*x[2]*y[0]*y[0]*z[0] + 2*c64*x[2]*y[0]*y[0]*z[0] + 
   c41*x[0]*x[0]*y[2]*z[0] + 2*c65*x[0]*x[0]*y[2]*z[0] + 2*c25*x[0]*y[0]*y[2]*z[0] + 2*c46*x[0]*y[0]*y[2]*z[0] + 
   2*c64*x[0]*y[0]*y[2]*z[0] + 2*c24*y[0]*y[0]*y[2]*z[0] + c42*y[0]*y[0]*y[2]*z[0] + 
   c13*x[0]*x[2]*z[0]*z[0] + 2*c55*x[0]*x[2]*z[0]*z[0] + 2*c54*x[2]*y[0]*z[0]*z[0] + 
   c63*x[2]*y[0]*z[0]*z[0] + 2*c45*x[0]*y[2]*z[0]*z[0] + c63*x[0]*y[2]*z[0]*z[0] + 
   c23*y[0]*y[2]*z[0]*z[0] + 2*c44*y[0]*y[2]*z[0]*z[0] + c53*x[2]*z[0]*z[0]*z[0] + 
   c43*y[2]*z[0]*z[0]*z[0] + c51*x[0]*x[0]*x[0]*z[2] + c41*x[0]*x[0]*y[0]*z[2] + 
   2*c56*x[0]*x[0]*y[0]*z[2] + 2*c46*x[0]*y[0]*y[0]*z[2] + c52*x[0]*y[0]*y[0]*z[2] + 
   c42*y[0]*y[0]*y[0]*z[2] + c31*x[0]*x[0]*z[0]*z[2] + 2*c55*x[0]*x[0]*z[0]*z[2] + 2*c36*x[0]*y[0]*z[0]*z[2] + 
   2*c45*x[0]*y[0]*z[0]*z[2] + 2*c54*x[0]*y[0]*z[0]*z[2] + c32*y[0]*y[0]*z[0]*z[2] + 2*c44*y[0]*y[0]*z[0]*z[2] + 
   2*c35*x[0]*z[0]*z[0]*z[2] + c53*x[0]*z[0]*z[0]*z[2] + 2*c34*y[0]*z[0]*z[0]*z[2] + 
   c43*y[0]*z[0]*z[0]*z[2] + c33*z[0]*z[0]*z[0]*z[2];
   
		mod(4,1) += c11*x[0]*x[1]*x[1]*x[2] + c61*x[1]*x[1]*x[2]*y[0] + 2*c16*x[0]*x[1]*x[2]*y[1] + 2*c66*x[1]*x[2]*y[0]*y[1] + 
   c12*x[0]*x[2]*y[1]*y[1] + c62*x[2]*y[0]*y[1]*y[1] + c61*x[0]*x[1]*x[1]*y[2] + 
   c21*x[1]*x[1]*y[0]*y[2] + 2*c66*x[0]*x[1]*y[1]*y[2] + 2*c26*x[1]*y[0]*y[1]*y[2] + c62*x[0]*y[1]*y[1]*y[2] + 
   c22*y[0]*y[1]*y[1]*y[2] + c51*x[1]*x[1]*x[2]*z[0] + 2*c56*x[1]*x[2]*y[1]*z[0] + 
   c52*x[2]*y[1]*y[1]*z[0] + c41*x[1]*x[1]*y[2]*z[0] + 2*c46*x[1]*y[1]*y[2]*z[0] + 
   c42*y[1]*y[1]*y[2]*z[0] + 2*c15*x[0]*x[1]*x[2]*z[1] + 2*c65*x[1]*x[2]*y[0]*z[1] + 2*c14*x[0]*x[2]*y[1]*z[1] + 
   2*c64*x[2]*y[0]*y[1]*z[1] + 2*c65*x[0]*x[1]*y[2]*z[1] + 2*c25*x[1]*y[0]*y[2]*z[1] + 2*c64*x[0]*y[1]*y[2]*z[1] + 
   2*c24*y[0]*y[1]*y[2]*z[1] + 2*c55*x[1]*x[2]*z[0]*z[1] + 2*c54*x[2]*y[1]*z[0]*z[1] + 2*c45*x[1]*y[2]*z[0]*z[1] + 
   2*c44*y[1]*y[2]*z[0]*z[1] + c13*x[0]*x[2]*z[1]*z[1] + c63*x[2]*y[0]*z[1]*z[1] + 
   c63*x[0]*y[2]*z[1]*z[1] + c23*y[0]*y[2]*z[1]*z[1] + c53*x[2]*z[0]*z[1]*z[1] + 
   c43*y[2]*z[0]*z[1]*z[1] + c51*x[0]*x[1]*x[1]*z[2] + c41*x[1]*x[1]*y[0]*z[2] + 
   2*c56*x[0]*x[1]*y[1]*z[2] + 2*c46*x[1]*y[0]*y[1]*z[2] + c52*x[0]*y[1]*y[1]*z[2] + c42*y[0]*y[1]*y[1]*z[2] + 
   c31*x[1]*x[1]*z[0]*z[2] + 2*c36*x[1]*y[1]*z[0]*z[2] + c32*y[1]*y[1]*z[0]*z[2] + 2*c55*x[0]*x[1]*z[1]*z[2] + 
   2*c45*x[1]*y[0]*z[1]*z[2] + 2*c54*x[0]*y[1]*z[1]*z[2] + 2*c44*y[0]*y[1]*z[1]*z[2] + 2*c35*x[1]*z[0]*z[1]*z[2] + 
   2*c34*y[1]*z[0]*z[1]*z[2] + c53*x[0]*z[1]*z[1]*z[2] + c43*y[0]*z[1]*z[1]*z[2] + c33*z[0]*z[1]*z[1]*z[2];
   
		mod(4,2) += c11*x[0]*x[2]*x[2]*x[2] + c61*x[2]*x[2]*x[2]*y[0] + 2*c16*x[0]*x[2]*x[2]*y[2] + c61*x[0]*x[2]*x[2]*y[2] + 
   c21*x[2]*x[2]*y[0]*y[2] + 2*c66*x[2]*x[2]*y[0]*y[2] + c12*x[0]*x[2]*y[2]*y[2] + 
   2*c66*x[0]*x[2]*y[2]*y[2] + 2*c26*x[2]*y[0]*y[2]*y[2] + c62*x[2]*y[0]*y[2]*y[2] + 
   c62*x[0]*y[2]*y[2]*y[2] + c22*y[0]*y[2]*y[2]*y[2] + c51*x[2]*x[2]*x[2]*z[0] + c41*x[2]*x[2]*y[2]*z[0] + 
   2*c56*x[2]*x[2]*y[2]*z[0] + 2*c46*x[2]*y[2]*y[2]*z[0] + c52*x[2]*y[2]*y[2]*z[0] + 
   c42*y[2]*y[2]*y[2]*z[0] + 2*c15*x[0]*x[2]*x[2]*z[2] + c51*x[0]*x[2]*x[2]*z[2] + 
   c41*x[2]*x[2]*y[0]*z[2] + 2*c65*x[2]*x[2]*y[0]*z[2] + 2*c14*x[0]*x[2]*y[2]*z[2] + 2*c56*x[0]*x[2]*y[2]*z[2] + 
   2*c65*x[0]*x[2]*y[2]*z[2] + 2*c25*x[2]*y[0]*y[2]*z[2] + 2*c46*x[2]*y[0]*y[2]*z[2] + 2*c64*x[2]*y[0]*y[2]*z[2] + 
   c52*x[0]*y[2]*y[2]*z[2] + 2*c64*x[0]*y[2]*y[2]*z[2] + 2*c24*y[0]*y[2]*y[2]*z[2] + 
   c42*y[0]*y[2]*y[2]*z[2] + c31*x[2]*x[2]*z[0]*z[2] + 2*c55*x[2]*x[2]*z[0]*z[2] + 
   2*c36*x[2]*y[2]*z[0]*z[2] + 2*c45*x[2]*y[2]*z[0]*z[2] + 2*c54*x[2]*y[2]*z[0]*z[2] + c32*y[2]*y[2]*z[0]*z[2] + 
   2*c44*y[2]*y[2]*z[0]*z[2] + c13*x[0]*x[2]*z[2]*z[2] + 2*c55*x[0]*x[2]*z[2]*z[2] + 
   2*c45*x[2]*y[0]*z[2]*z[2] + c63*x[2]*y[0]*z[2]*z[2] + 2*c54*x[0]*y[2]*z[2]*z[2] + 
   c63*x[0]*y[2]*z[2]*z[2] + c23*y[0]*y[2]*z[2]*z[2] + 2*c44*y[0]*y[2]*z[2]*z[2] + 
   2*c35*x[2]*z[0]*z[2]*z[2] + c53*x[2]*z[0]*z[2]*z[2] + 2*c34*y[2]*z[0]*z[2]*z[2] + 
   c43*y[2]*z[0]*z[2]*z[2] + c53*x[0]*z[2]*z[2]*z[2] + c43*y[0]*z[2]*z[2]*z[2] + c33*z[0]*z[2]*z[2]*z[2];
   
		mod(4,3) += c11*x[0]*x[1]*x[2]*x[2] + c61*x[1]*x[2]*x[2]*y[0] + c16*x[0]*x[2]*x[2]*y[1] + 
   c66*x[2]*x[2]*y[0]*y[1] + c16*x[0]*x[1]*x[2]*y[2] + c61*x[0]*x[1]*x[2]*y[2] + c21*x[1]*x[2]*y[0]*y[2] + 
   c66*x[1]*x[2]*y[0]*y[2] + c12*x[0]*x[2]*y[1]*y[2] + c66*x[0]*x[2]*y[1]*y[2] + c26*x[2]*y[0]*y[1]*y[2] + c62*x[2]*y[0]*y[1]*y[2] + 
   c66*x[0]*x[1]*y[2]*y[2] + c26*x[1]*y[0]*y[2]*y[2] + c62*x[0]*y[1]*y[2]*y[2] + 
   c22*y[0]*y[1]*y[2]*y[2] + c51*x[1]*x[2]*x[2]*z[0] + c56*x[2]*x[2]*y[1]*z[0] + c41*x[1]*x[2]*y[2]*z[0] + 
   c56*x[1]*x[2]*y[2]*z[0] + c46*x[2]*y[1]*y[2]*z[0] + c52*x[2]*y[1]*y[2]*z[0] + c46*x[1]*y[2]*y[2]*z[0] + 
   c42*y[1]*y[2]*y[2]*z[0] + c15*x[0]*x[2]*x[2]*z[1] + c65*x[2]*x[2]*y[0]*z[1] + c14*x[0]*x[2]*y[2]*z[1] + 
   c65*x[0]*x[2]*y[2]*z[1] + c25*x[2]*y[0]*y[2]*z[1] + c64*x[2]*y[0]*y[2]*z[1] + c64*x[0]*y[2]*y[2]*z[1] + 
   c24*y[0]*y[2]*y[2]*z[1] + c55*x[2]*x[2]*z[0]*z[1] + c45*x[2]*y[2]*z[0]*z[1] + c54*x[2]*y[2]*z[0]*z[1] + 
   c44*y[2]*y[2]*z[0]*z[1] + c15*x[0]*x[1]*x[2]*z[2] + c51*x[0]*x[1]*x[2]*z[2] + c41*x[1]*x[2]*y[0]*z[2] + 
   c65*x[1]*x[2]*y[0]*z[2] + c14*x[0]*x[2]*y[1]*z[2] + c56*x[0]*x[2]*y[1]*z[2] + c46*x[2]*y[0]*y[1]*z[2] + c64*x[2]*y[0]*y[1]*z[2] + 
   c56*x[0]*x[1]*y[2]*z[2] + c65*x[0]*x[1]*y[2]*z[2] + c25*x[1]*y[0]*y[2]*z[2] + c46*x[1]*y[0]*y[2]*z[2] + c52*x[0]*y[1]*y[2]*z[2] + 
   c64*x[0]*y[1]*y[2]*z[2] + c24*y[0]*y[1]*y[2]*z[2] + c42*y[0]*y[1]*y[2]*z[2] + c31*x[1]*x[2]*z[0]*z[2] + c55*x[1]*x[2]*z[0]*z[2] + 
   c36*x[2]*y[1]*z[0]*z[2] + c54*x[2]*y[1]*z[0]*z[2] + c36*x[1]*y[2]*z[0]*z[2] + c45*x[1]*y[2]*z[0]*z[2] + c32*y[1]*y[2]*z[0]*z[2] + 
   c44*y[1]*y[2]*z[0]*z[2] + c13*x[0]*x[2]*z[1]*z[2] + c55*x[0]*x[2]*z[1]*z[2] + c45*x[2]*y[0]*z[1]*z[2] + c63*x[2]*y[0]*z[1]*z[2] + 
   c54*x[0]*y[2]*z[1]*z[2] + c63*x[0]*y[2]*z[1]*z[2] + c23*y[0]*y[2]*z[1]*z[2] + c44*y[0]*y[2]*z[1]*z[2] + c35*x[2]*z[0]*z[1]*z[2] + 
   c53*x[2]*z[0]*z[1]*z[2] + c34*y[2]*z[0]*z[1]*z[2] + c43*y[2]*z[0]*z[1]*z[2] + c55*x[0]*x[1]*z[2]*z[2] + 
   c45*x[1]*y[0]*z[2]*z[2] + c54*x[0]*y[1]*z[2]*z[2] + c44*y[0]*y[1]*z[2]*z[2] + 
   c35*x[1]*z[0]*z[2]*z[2] + c34*y[1]*z[0]*z[2]*z[2] + c53*x[0]*z[1]*z[2]*z[2] + 
   c43*y[0]*z[1]*z[2]*z[2] + c33*z[0]*z[1]*z[2]*z[2];
   
		mod(4,5) += c11*x[0]*x[0]*x[1]*x[2] + c16*x[0]*x[1]*x[2]*y[0] + c61*x[0]*x[1]*x[2]*y[0] + c66*x[1]*x[2]*y[0]*y[0] + 
   c16*x[0]*x[0]*x[2]*y[1] + c12*x[0]*x[2]*y[0]*y[1] + c66*x[0]*x[2]*y[0]*y[1] + c62*x[2]*y[0]*y[0]*y[1] + 
   c61*x[0]*x[0]*x[1]*y[2] + c21*x[0]*x[1]*y[0]*y[2] + c66*x[0]*x[1]*y[0]*y[2] + c26*x[1]*y[0]*y[0]*y[2] + 
   c66*x[0]*x[0]*y[1]*y[2] + c26*x[0]*y[0]*y[1]*y[2] + c62*x[0]*y[0]*y[1]*y[2] + c22*y[0]*y[0]*y[1]*y[2] + 
   c15*x[0]*x[1]*x[2]*z[0] + c51*x[0]*x[1]*x[2]*z[0] + c56*x[1]*x[2]*y[0]*z[0] + c65*x[1]*x[2]*y[0]*z[0] + c14*x[0]*x[2]*y[1]*z[0] + 
   c56*x[0]*x[2]*y[1]*z[0] + c52*x[2]*y[0]*y[1]*z[0] + c64*x[2]*y[0]*y[1]*z[0] + c41*x[0]*x[1]*y[2]*z[0] + c65*x[0]*x[1]*y[2]*z[0] + 
   c25*x[1]*y[0]*y[2]*z[0] + c46*x[1]*y[0]*y[2]*z[0] + c46*x[0]*y[1]*y[2]*z[0] + c64*x[0]*y[1]*y[2]*z[0] + c24*y[0]*y[1]*y[2]*z[0] + 
   c42*y[0]*y[1]*y[2]*z[0] + c55*x[1]*x[2]*z[0]*z[0] + c54*x[2]*y[1]*z[0]*z[0] + c45*x[1]*y[2]*z[0]*z[0] + 
   c44*y[1]*y[2]*z[0]*z[0] + c15*x[0]*x[0]*x[2]*z[1] + c14*x[0]*x[2]*y[0]*z[1] + c65*x[0]*x[2]*y[0]*z[1] + 
   c64*x[2]*y[0]*y[0]*z[1] + c65*x[0]*x[0]*y[2]*z[1] + c25*x[0]*y[0]*y[2]*z[1] + c64*x[0]*y[0]*y[2]*z[1] + 
   c24*y[0]*y[0]*y[2]*z[1] + c13*x[0]*x[2]*z[0]*z[1] + c55*x[0]*x[2]*z[0]*z[1] + c54*x[2]*y[0]*z[0]*z[1] + 
   c63*x[2]*y[0]*z[0]*z[1] + c45*x[0]*y[2]*z[0]*z[1] + c63*x[0]*y[2]*z[0]*z[1] + c23*y[0]*y[2]*z[0]*z[1] + c44*y[0]*y[2]*z[0]*z[1] + 
   c53*x[2]*z[0]*z[0]*z[1] + c43*y[2]*z[0]*z[0]*z[1] + c51*x[0]*x[0]*x[1]*z[2] + c41*x[0]*x[1]*y[0]*z[2] + 
   c56*x[0]*x[1]*y[0]*z[2] + c46*x[1]*y[0]*y[0]*z[2] + c56*x[0]*x[0]*y[1]*z[2] + c46*x[0]*y[0]*y[1]*z[2] + 
   c52*x[0]*y[0]*y[1]*z[2] + c42*y[0]*y[0]*y[1]*z[2] + c31*x[0]*x[1]*z[0]*z[2] + c55*x[0]*x[1]*z[0]*z[2] + 
   c36*x[1]*y[0]*z[0]*z[2] + c45*x[1]*y[0]*z[0]*z[2] + c36*x[0]*y[1]*z[0]*z[2] + c54*x[0]*y[1]*z[0]*z[2] + c32*y[0]*y[1]*z[0]*z[2] + 
   c44*y[0]*y[1]*z[0]*z[2] + c35*x[1]*z[0]*z[0]*z[2] + c34*y[1]*z[0]*z[0]*z[2] + c55*x[0]*x[0]*z[1]*z[2] + 
   c45*x[0]*y[0]*z[1]*z[2] + c54*x[0]*y[0]*z[1]*z[2] + c44*y[0]*y[0]*z[1]*z[2] + c35*x[0]*z[0]*z[1]*z[2] + 
   c53*x[0]*z[0]*z[1]*z[2] + c34*y[0]*z[0]*z[1]*z[2] + c43*y[0]*z[0]*z[1]*z[2] + c33*z[0]*z[0]*z[1]*z[2];
   
		mod(5,0) += c11*x[0]*x[0]*x[0]*x[1] + 2*c16*x[0]*x[0]*x[1]*y[0] + c61*x[0]*x[0]*x[1]*y[0] + 
   c12*x[0]*x[1]*y[0]*y[0] + 2*c66*x[0]*x[1]*y[0]*y[0] + c62*x[1]*y[0]*y[0]*y[0] + 
   c61*x[0]*x[0]*x[0]*y[1] + c21*x[0]*x[0]*y[0]*y[1] + 2*c66*x[0]*x[0]*y[0]*y[1] + 
   2*c26*x[0]*y[0]*y[0]*y[1] + c62*x[0]*y[0]*y[0]*y[1] + c22*y[0]*y[0]*y[0]*y[1] + 
   2*c15*x[0]*x[0]*x[1]*z[0] + c51*x[0]*x[0]*x[1]*z[0] + 2*c14*x[0]*x[1]*y[0]*z[0] + 2*c56*x[0]*x[1]*y[0]*z[0] + 
   2*c65*x[0]*x[1]*y[0]*z[0] + c52*x[1]*y[0]*y[0]*z[0] + 2*c64*x[1]*y[0]*y[0]*z[0] + 
   c41*x[0]*x[0]*y[1]*z[0] + 2*c65*x[0]*x[0]*y[1]*z[0] + 2*c25*x[0]*y[0]*y[1]*z[0] + 2*c46*x[0]*y[0]*y[1]*z[0] + 
   2*c64*x[0]*y[0]*y[1]*z[0] + 2*c24*y[0]*y[0]*y[1]*z[0] + c42*y[0]*y[0]*y[1]*z[0] + 
   c13*x[0]*x[1]*z[0]*z[0] + 2*c55*x[0]*x[1]*z[0]*z[0] + 2*c54*x[1]*y[0]*z[0]*z[0] + 
   c63*x[1]*y[0]*z[0]*z[0] + 2*c45*x[0]*y[1]*z[0]*z[0] + c63*x[0]*y[1]*z[0]*z[0] + 
   c23*y[0]*y[1]*z[0]*z[0] + 2*c44*y[0]*y[1]*z[0]*z[0] + c53*x[1]*z[0]*z[0]*z[0] + 
   c43*y[1]*z[0]*z[0]*z[0] + c51*x[0]*x[0]*x[0]*z[1] + c41*x[0]*x[0]*y[0]*z[1] + 
   2*c56*x[0]*x[0]*y[0]*z[1] + 2*c46*x[0]*y[0]*y[0]*z[1] + c52*x[0]*y[0]*y[0]*z[1] + 
   c42*y[0]*y[0]*y[0]*z[1] + c31*x[0]*x[0]*z[0]*z[1] + 2*c55*x[0]*x[0]*z[0]*z[1] + 2*c36*x[0]*y[0]*z[0]*z[1] + 
   2*c45*x[0]*y[0]*z[0]*z[1] + 2*c54*x[0]*y[0]*z[0]*z[1] + c32*y[0]*y[0]*z[0]*z[1] + 2*c44*y[0]*y[0]*z[0]*z[1] + 
   2*c35*x[0]*z[0]*z[0]*z[1] + c53*x[0]*z[0]*z[0]*z[1] + 2*c34*y[0]*z[0]*z[0]*z[1] + 
   c43*y[0]*z[0]*z[0]*z[1] + c33*z[0]*z[0]*z[0]*z[1];
   
		mod(5,1) += c11*x[0]*x[1]*x[1]*x[1] + c61*x[1]*x[1]*x[1]*y[0] + 2*c16*x[0]*x[1]*x[1]*y[1] + c61*x[0]*x[1]*x[1]*y[1] + 
   c21*x[1]*x[1]*y[0]*y[1] + 2*c66*x[1]*x[1]*y[0]*y[1] + c12*x[0]*x[1]*y[1]*y[1] + 
   2*c66*x[0]*x[1]*y[1]*y[1] + 2*c26*x[1]*y[0]*y[1]*y[1] + c62*x[1]*y[0]*y[1]*y[1] + 
   c62*x[0]*y[1]*y[1]*y[1] + c22*y[0]*y[1]*y[1]*y[1] + c51*x[1]*x[1]*x[1]*z[0] + c41*x[1]*x[1]*y[1]*z[0] + 
   2*c56*x[1]*x[1]*y[1]*z[0] + 2*c46*x[1]*y[1]*y[1]*z[0] + c52*x[1]*y[1]*y[1]*z[0] + 
   c42*y[1]*y[1]*y[1]*z[0] + 2*c15*x[0]*x[1]*x[1]*z[1] + c51*x[0]*x[1]*x[1]*z[1] + 
   c41*x[1]*x[1]*y[0]*z[1] + 2*c65*x[1]*x[1]*y[0]*z[1] + 2*c14*x[0]*x[1]*y[1]*z[1] + 2*c56*x[0]*x[1]*y[1]*z[1] + 
   2*c65*x[0]*x[1]*y[1]*z[1] + 2*c25*x[1]*y[0]*y[1]*z[1] + 2*c46*x[1]*y[0]*y[1]*z[1] + 2*c64*x[1]*y[0]*y[1]*z[1] + 
   c52*x[0]*y[1]*y[1]*z[1] + 2*c64*x[0]*y[1]*y[1]*z[1] + 2*c24*y[0]*y[1]*y[1]*z[1] + 
   c42*y[0]*y[1]*y[1]*z[1] + c31*x[1]*x[1]*z[0]*z[1] + 2*c55*x[1]*x[1]*z[0]*z[1] + 
   2*c36*x[1]*y[1]*z[0]*z[1] + 2*c45*x[1]*y[1]*z[0]*z[1] + 2*c54*x[1]*y[1]*z[0]*z[1] + c32*y[1]*y[1]*z[0]*z[1] + 
   2*c44*y[1]*y[1]*z[0]*z[1] + c13*x[0]*x[1]*z[1]*z[1] + 2*c55*x[0]*x[1]*z[1]*z[1] + 
   2*c45*x[1]*y[0]*z[1]*z[1] + c63*x[1]*y[0]*z[1]*z[1] + 2*c54*x[0]*y[1]*z[1]*z[1] + 
   c63*x[0]*y[1]*z[1]*z[1] + c23*y[0]*y[1]*z[1]*z[1] + 2*c44*y[0]*y[1]*z[1]*z[1] + 
   2*c35*x[1]*z[0]*z[1]*z[1] + c53*x[1]*z[0]*z[1]*z[1] + 2*c34*y[1]*z[0]*z[1]*z[1] + 
   c43*y[1]*z[0]*z[1]*z[1] + c53*x[0]*z[1]*z[1]*z[1] + c43*y[0]*z[1]*z[1]*z[1] + c33*z[0]*z[1]*z[1]*z[1];
   
		mod(5,2) += c11*x[0]*x[1]*x[2]*x[2] + c61*x[1]*x[2]*x[2]*y[0] + c61*x[0]*x[2]*x[2]*y[1] + 
   c21*x[2]*x[2]*y[0]*y[1] + 2*c16*x[0]*x[1]*x[2]*y[2] + 2*c66*x[1]*x[2]*y[0]*y[2] + 2*c66*x[0]*x[2]*y[1]*y[2] + 
   2*c26*x[2]*y[0]*y[1]*y[2] + c12*x[0]*x[1]*y[2]*y[2] + c62*x[1]*y[0]*y[2]*y[2] + 
   c62*x[0]*y[1]*y[2]*y[2] + c22*y[0]*y[1]*y[2]*y[2] + c51*x[1]*x[2]*x[2]*z[0] + 
   c41*x[2]*x[2]*y[1]*z[0] + 2*c56*x[1]*x[2]*y[2]*z[0] + 2*c46*x[2]*y[1]*y[2]*z[0] + c52*x[1]*y[2]*y[2]*z[0] + 
   c42*y[1]*y[2]*y[2]*z[0] + c51*x[0]*x[2]*x[2]*z[1] + c41*x[2]*x[2]*y[0]*z[1] + 
   2*c56*x[0]*x[2]*y[2]*z[1] + 2*c46*x[2]*y[0]*y[2]*z[1] + c52*x[0]*y[2]*y[2]*z[1] + c42*y[0]*y[2]*y[2]*z[1] + 
   c31*x[2]*x[2]*z[0]*z[1] + 2*c36*x[2]*y[2]*z[0]*z[1] + c32*y[2]*y[2]*z[0]*z[1] + 2*c15*x[0]*x[1]*x[2]*z[2] + 
   2*c65*x[1]*x[2]*y[0]*z[2] + 2*c65*x[0]*x[2]*y[1]*z[2] + 2*c25*x[2]*y[0]*y[1]*z[2] + 2*c14*x[0]*x[1]*y[2]*z[2] + 
   2*c64*x[1]*y[0]*y[2]*z[2] + 2*c64*x[0]*y[1]*y[2]*z[2] + 2*c24*y[0]*y[1]*y[2]*z[2] + 2*c55*x[1]*x[2]*z[0]*z[2] + 
   2*c45*x[2]*y[1]*z[0]*z[2] + 2*c54*x[1]*y[2]*z[0]*z[2] + 2*c44*y[1]*y[2]*z[0]*z[2] + 2*c55*x[0]*x[2]*z[1]*z[2] + 
   2*c45*x[2]*y[0]*z[1]*z[2] + 2*c54*x[0]*y[2]*z[1]*z[2] + 2*c44*y[0]*y[2]*z[1]*z[2] + 2*c35*x[2]*z[0]*z[1]*z[2] + 
   2*c34*y[2]*z[0]*z[1]*z[2] + c13*x[0]*x[1]*z[2]*z[2] + c63*x[1]*y[0]*z[2]*z[2] + 
   c63*x[0]*y[1]*z[2]*z[2] + c23*y[0]*y[1]*z[2]*z[2] + c53*x[1]*z[0]*z[2]*z[2] + 
   c43*y[1]*z[0]*z[2]*z[2] + c53*x[0]*z[1]*z[2]*z[2] + c43*y[0]*z[1]*z[2]*z[2] + 
   c33*z[0]*z[1]*z[2]*z[2];
   
		mod(5,3) += c11*x[0]*x[1]*x[1]*x[2] + c61*x[1]*x[1]*x[2]*y[0] + c16*x[0]*x[1]*x[2]*y[1] + c61*x[0]*x[1]*x[2]*y[1] + 
   c21*x[1]*x[2]*y[0]*y[1] + c66*x[1]*x[2]*y[0]*y[1] + c66*x[0]*x[2]*y[1]*y[1] + c26*x[2]*y[0]*y[1]*y[1] + 
   c16*x[0]*x[1]*x[1]*y[2] + c66*x[1]*x[1]*y[0]*y[2] + c12*x[0]*x[1]*y[1]*y[2] + c66*x[0]*x[1]*y[1]*y[2] + 
   c26*x[1]*y[0]*y[1]*y[2] + c62*x[1]*y[0]*y[1]*y[2] + c62*x[0]*y[1]*y[1]*y[2] + c22*y[0]*y[1]*y[1]*y[2] + 
   c51*x[1]*x[1]*x[2]*z[0] + c41*x[1]*x[2]*y[1]*z[0] + c56*x[1]*x[2]*y[1]*z[0] + c46*x[2]*y[1]*y[1]*z[0] + 
   c56*x[1]*x[1]*y[2]*z[0] + c46*x[1]*y[1]*y[2]*z[0] + c52*x[1]*y[1]*y[2]*z[0] + c42*y[1]*y[1]*y[2]*z[0] + 
   c15*x[0]*x[1]*x[2]*z[1] + c51*x[0]*x[1]*x[2]*z[1] + c41*x[1]*x[2]*y[0]*z[1] + c65*x[1]*x[2]*y[0]*z[1] + c56*x[0]*x[2]*y[1]*z[1] + 
   c65*x[0]*x[2]*y[1]*z[1] + c25*x[2]*y[0]*y[1]*z[1] + c46*x[2]*y[0]*y[1]*z[1] + c14*x[0]*x[1]*y[2]*z[1] + c56*x[0]*x[1]*y[2]*z[1] + 
   c46*x[1]*y[0]*y[2]*z[1] + c64*x[1]*y[0]*y[2]*z[1] + c52*x[0]*y[1]*y[2]*z[1] + c64*x[0]*y[1]*y[2]*z[1] + c24*y[0]*y[1]*y[2]*z[1] + 
   c42*y[0]*y[1]*y[2]*z[1] + c31*x[1]*x[2]*z[0]*z[1] + c55*x[1]*x[2]*z[0]*z[1] + c36*x[2]*y[1]*z[0]*z[1] + c45*x[2]*y[1]*z[0]*z[1] + 
   c36*x[1]*y[2]*z[0]*z[1] + c54*x[1]*y[2]*z[0]*z[1] + c32*y[1]*y[2]*z[0]*z[1] + c44*y[1]*y[2]*z[0]*z[1] + 
   c55*x[0]*x[2]*z[1]*z[1] + c45*x[2]*y[0]*z[1]*z[1] + c54*x[0]*y[2]*z[1]*z[1] + 
   c44*y[0]*y[2]*z[1]*z[1] + c35*x[2]*z[0]*z[1]*z[1] + c34*y[2]*z[0]*z[1]*z[1] + 
   c15*x[0]*x[1]*x[1]*z[2] + c65*x[1]*x[1]*y[0]*z[2] + c14*x[0]*x[1]*y[1]*z[2] + c65*x[0]*x[1]*y[1]*z[2] + 
   c25*x[1]*y[0]*y[1]*z[2] + c64*x[1]*y[0]*y[1]*z[2] + c64*x[0]*y[1]*y[1]*z[2] + c24*y[0]*y[1]*y[1]*z[2] + 
   c55*x[1]*x[1]*z[0]*z[2] + c45*x[1]*y[1]*z[0]*z[2] + c54*x[1]*y[1]*z[0]*z[2] + c44*y[1]*y[1]*z[0]*z[2] + 
   c13*x[0]*x[1]*z[1]*z[2] + c55*x[0]*x[1]*z[1]*z[2] + c45*x[1]*y[0]*z[1]*z[2] + c63*x[1]*y[0]*z[1]*z[2] + c54*x[0]*y[1]*z[1]*z[2] + 
   c63*x[0]*y[1]*z[1]*z[2] + c23*y[0]*y[1]*z[1]*z[2] + c44*y[0]*y[1]*z[1]*z[2] + c35*x[1]*z[0]*z[1]*z[2] + c53*x[1]*z[0]*z[1]*z[2] + 
   c34*y[1]*z[0]*z[1]*z[2] + c43*y[1]*z[0]*z[1]*z[2] + c53*x[0]*z[1]*z[1]*z[2] + c43*y[0]*z[1]*z[1]*z[2] + 
   c33*z[0]*z[1]*z[1]*z[2];
   
		mod(5,4) += c11*x[0]*x[0]*x[1]*x[2] + c16*x[0]*x[1]*x[2]*y[0] + c61*x[0]*x[1]*x[2]*y[0] + c66*x[1]*x[2]*y[0]*y[0] + 
   c61*x[0]*x[0]*x[2]*y[1] + c21*x[0]*x[2]*y[0]*y[1] + c66*x[0]*x[2]*y[0]*y[1] + c26*x[2]*y[0]*y[0]*y[1] + 
   c16*x[0]*x[0]*x[1]*y[2] + c12*x[0]*x[1]*y[0]*y[2] + c66*x[0]*x[1]*y[0]*y[2] + c62*x[1]*y[0]*y[0]*y[2] + 
   c66*x[0]*x[0]*y[1]*y[2] + c26*x[0]*y[0]*y[1]*y[2] + c62*x[0]*y[0]*y[1]*y[2] + c22*y[0]*y[0]*y[1]*y[2] + 
   c15*x[0]*x[1]*x[2]*z[0] + c51*x[0]*x[1]*x[2]*z[0] + c56*x[1]*x[2]*y[0]*z[0] + c65*x[1]*x[2]*y[0]*z[0] + c41*x[0]*x[2]*y[1]*z[0] + 
   c65*x[0]*x[2]*y[1]*z[0] + c25*x[2]*y[0]*y[1]*z[0] + c46*x[2]*y[0]*y[1]*z[0] + c14*x[0]*x[1]*y[2]*z[0] + c56*x[0]*x[1]*y[2]*z[0] + 
   c52*x[1]*y[0]*y[2]*z[0] + c64*x[1]*y[0]*y[2]*z[0] + c46*x[0]*y[1]*y[2]*z[0] + c64*x[0]*y[1]*y[2]*z[0] + c24*y[0]*y[1]*y[2]*z[0] + 
   c42*y[0]*y[1]*y[2]*z[0] + c55*x[1]*x[2]*z[0]*z[0] + c45*x[2]*y[1]*z[0]*z[0] + c54*x[1]*y[2]*z[0]*z[0] + 
   c44*y[1]*y[2]*z[0]*z[0] + c51*x[0]*x[0]*x[2]*z[1] + c41*x[0]*x[2]*y[0]*z[1] + c56*x[0]*x[2]*y[0]*z[1] + 
   c46*x[2]*y[0]*y[0]*z[1] + c56*x[0]*x[0]*y[2]*z[1] + c46*x[0]*y[0]*y[2]*z[1] + c52*x[0]*y[0]*y[2]*z[1] + 
   c42*y[0]*y[0]*y[2]*z[1] + c31*x[0]*x[2]*z[0]*z[1] + c55*x[0]*x[2]*z[0]*z[1] + c36*x[2]*y[0]*z[0]*z[1] + 
   c45*x[2]*y[0]*z[0]*z[1] + c36*x[0]*y[2]*z[0]*z[1] + c54*x[0]*y[2]*z[0]*z[1] + c32*y[0]*y[2]*z[0]*z[1] + c44*y[0]*y[2]*z[0]*z[1] + 
   c35*x[2]*z[0]*z[0]*z[1] + c34*y[2]*z[0]*z[0]*z[1] + c15*x[0]*x[0]*x[1]*z[2] + c14*x[0]*x[1]*y[0]*z[2] + 
   c65*x[0]*x[1]*y[0]*z[2] + c64*x[1]*y[0]*y[0]*z[2] + c65*x[0]*x[0]*y[1]*z[2] + c25*x[0]*y[0]*y[1]*z[2] + 
   c64*x[0]*y[0]*y[1]*z[2] + c24*y[0]*y[0]*y[1]*z[2] + c13*x[0]*x[1]*z[0]*z[2] + c55*x[0]*x[1]*z[0]*z[2] + 
   c54*x[1]*y[0]*z[0]*z[2] + c63*x[1]*y[0]*z[0]*z[2] + c45*x[0]*y[1]*z[0]*z[2] + c63*x[0]*y[1]*z[0]*z[2] + c23*y[0]*y[1]*z[0]*z[2] + 
   c44*y[0]*y[1]*z[0]*z[2] + c53*x[1]*z[0]*z[0]*z[2] + c43*y[1]*z[0]*z[0]*z[2] + c55*x[0]*x[0]*z[1]*z[2] + 
   c45*x[0]*y[0]*z[1]*z[2] + c54*x[0]*y[0]*z[1]*z[2] + c44*y[0]*y[0]*z[1]*z[2] + c35*x[0]*z[0]*z[1]*z[2] + 
   c53*x[0]*z[0]*z[1]*z[2] + c34*y[0]*z[0]*z[1]*z[2] + c43*y[0]*z[0]*z[1]*z[2] + c33*z[0]*z[0]*z[1]*z[2];
   
	}
	else if (fNumFibStress == 3)
	{
		const double& c11 = cf(0,0);
		const double& c22 = cf(1,1);
		const double& c33 = cf(2,2);
	
		const double& c23 = cf(1,2);
		const double& c32 = cf(2,1);
		const double& c13 = cf(0,2);
		const double& c31 = cf(2,0);
		const double& c12 = cf(0,1);
		const double& c21 = cf(1,0);
	
		const double* x = fQ(0);	
		const double* y = fQ(1);
		
		mod(0,0) += c11*x[0]*x[0]*x[0]*x[0] + 2*c13*x[0]*x[0]*x[0]*y[0] + 2*c31*x[0]*x[0]*x[0]*y[0] + c12*x[0]*x[0]*y[0]*y[0] + c21*x[0]*x[0]*y[0]*y[0] 
			+ 4*c33*x[0]*x[0]*y[0]*y[0] + 2*c23*x[0]*y[0]*y[0]*y[0] + 2*c32*x[0]*y[0]*y[0]*y[0] + c22*y[0]*y[0]*y[0]*y[0];

		mod(1,1) += c11*x[1]*x[1]*x[1]*x[1] + 2*c13*x[1]*x[1]*x[1]*y[1] + 2*c31*x[1]*x[1]*x[1]*y[1] + c12*x[1]*x[1]*y[1]*y[1] + c21*x[1]*x[1]*y[1]*y[1] 
			+ 4*c33*x[1]*x[1]*y[1]*y[1] + 2*c23*x[1]*y[1]*y[1]*y[1] + 2*c32*x[1]*y[1]*y[1]*y[1] + c22*y[1]*y[1]*y[1]*y[1];

		mod(2,2) += c11*x[2]*x[2]*x[2]*x[2] + 2*c13*x[2]*x[2]*x[2]*y[2] + 2*c31*x[2]*x[2]*x[2]*y[2] + c12*x[2]*x[2]*y[2]*y[2] + c21*x[2]*x[2]*y[2]*y[2] 
			+ 4*c33*x[2]*x[2]*y[2]*y[2] + 2*c23*x[2]*y[2]*y[2]*y[2] + 2*c32*x[2]*y[2]*y[2]*y[2] + c22*y[2]*y[2]*y[2]*y[2];

		mod(3,3) += c11*x[1]*x[1]*x[2]*x[2] + c13*x[1]*x[2]*x[2]*y[1] + c31*x[1]*x[2]*x[2]*y[1] + c33*x[2]*x[2]*y[1]*y[1] + c13*x[1]*x[1]*x[2]*y[2] 
			+ c31*x[1]*x[1]*x[2]*y[2] + c12*x[1]*x[2]*y[1]*y[2] + c21*x[1]*x[2]*y[1]*y[2] + 2*c33*x[1]*x[2]*y[1]*y[2] + c23*x[2]*y[1]*y[1]*y[2] + c32*x[2]*y[1]*y[1]*y[2] 
			+ c33*x[1]*x[1]*y[2]*y[2] + c23*x[1]*y[1]*y[2]*y[2] + c32*x[1]*y[1]*y[2]*y[2] + c22*y[1]*y[1]*y[2]*y[2];
		
		mod(4,4) += c11*x[0]*x[0]*x[2]*x[2] + c13*x[0]*x[2]*x[2]*y[0] + c31*x[0]*x[2]*x[2]*y[0] + c33*x[2]*x[2]*y[0]*y[0] + c13*x[0]*x[0]*x[2]*y[2] 
			+ c31*x[0]*x[0]*x[2]*y[2] + c12*x[0]*x[2]*y[0]*y[2] + c21*x[0]*x[2]*y[0]*y[2] + 2*c33*x[0]*x[2]*y[0]*y[2] + c23*x[2]*y[0]*y[0]*y[2] + c32*x[2]*y[0]*y[0]*y[2] 
			+ c33*x[0]*x[0]*y[2]*y[2] + c23*x[0]*y[0]*y[2]*y[2] + c32*x[0]*y[0]*y[2]*y[2] + c22*y[0]*y[0]*y[2]*y[2];

		mod(5,5) += c11*x[0]*x[0]*x[1]*x[1] + c13*x[0]*x[1]*x[1]*y[0] + c31*x[0]*x[1]*x[1]*y[0] + c33*x[1]*x[1]*y[0]*y[0] + c13*x[0]*x[0]*x[1]*y[1] 
			+ c31*x[0]*x[0]*x[1]*y[1] + c12*x[0]*x[1]*y[0]*y[1] + c21*x[0]*x[1]*y[0]*y[1] + 2*c33*x[0]*x[1]*y[0]*y[1] + c23*x[1]*y[0]*y[0]*y[1] + c32*x[1]*y[0]*y[0]*y[1] 
			+ c33*x[0]*x[0]*y[1]*y[1] + c23*x[0]*y[0]*y[1]*y[1] + c32*x[0]*y[0]*y[1]*y[1] + c22*y[0]*y[0]*y[1]*y[1];
		
		
		mod(0,1) += c11*x[0]*x[0]*x[1]*x[1] + 2*c31*x[0]*x[1]*x[1]*y[0] + c21*x[1]*x[1]*y[0]*y[0] + 2*c13*x[0]*x[0]*x[1]*y[1] + 4*c33*x[0]*x[1]*y[0]*y[1] 
			+ 2*c23*x[1]*y[0]*y[0]*y[1] + c12*x[0]*x[0]*y[1]*y[1] + 2*c32*x[0]*y[0]*y[1]*y[1] + c22*y[0]*y[0]*y[1]*y[1];

		mod(0,2) += c11*x[0]*x[0]*x[2]*x[2] + 2*c31*x[0]*x[2]*x[2]*y[0] + c21*x[2]*x[2]*y[0]*y[0] + 2*c13*x[0]*x[0]*x[2]*y[2] + 4*c33*x[0]*x[2]*y[0]*y[2] 
			+ 2*c23*x[2]*y[0]*y[0]*y[2] + c12*x[0]*x[0]*y[2]*y[2] + 2*c32*x[0]*y[0]*y[2]*y[2] + c22*y[0]*y[0]*y[2]*y[2];

		mod(0,3) += c11*x[0]*x[0]*x[1]*x[2] + 2*c31*x[0]*x[1]*x[2]*y[0] + c21*x[1]*x[2]*y[0]*y[0] + c13*x[0]*x[0]*x[2]*y[1] + 2*c33*x[0]*x[2]*y[0]*y[1] + c23*x[2]*y[0]*y[0]*y[1] 
			+ c13*x[0]*x[0]*x[1]*y[2] + 2*c33*x[0]*x[1]*y[0]*y[2] + c23*x[1]*y[0]*y[0]*y[2] + c12*x[0]*x[0]*y[1]*y[2] + 2*c32*x[0]*y[0]*y[1]*y[2] + c22*y[0]*y[0]*y[1]*y[2];

		mod(0,4) += c11*x[0]*x[0]*x[0]*x[2] + c13*x[0]*x[0]*x[2]*y[0] + 2*c31*x[0]*x[0]*x[2]*y[0] + c21*x[0]*x[2]*y[0]*y[0] + 2*c33*x[0]*x[2]*y[0]*y[0] + c23*x[2]*y[0]*y[0]*y[0] 
			+ c13*x[0]*x[0]*x[0]*y[2] + c12*x[0]*x[0]*y[0]*y[2] + 2*c33*x[0]*x[0]*y[0]*y[2] + c23*x[0]*y[0]*y[0]*y[2] + 2*c32*x[0]*y[0]*y[0]*y[2] + c22*y[0]*y[0]*y[0]*y[2];

		mod(0,5) += c11*x[0]*x[0]*x[0]*x[1] + c13*x[0]*x[0]*x[1]*y[0] + 2*c31*x[0]*x[0]*x[1]*y[0] + c21*x[0]*x[1]*y[0]*y[0] + 2*c33*x[0]*x[1]*y[0]*y[0] + c23*x[1]*y[0]*y[0]*y[0] 
			+ c13*x[0]*x[0]*x[0]*y[1] + c12*x[0]*x[0]*y[0]*y[1] + 2*c33*x[0]*x[0]*y[0]*y[1] + c23*x[0]*y[0]*y[0]*y[1] + 2*c32*x[0]*y[0]*y[0]*y[1] + c22*y[0]*y[0]*y[0]*y[1];


		mod(1,0) += c11*x[0]*x[0]*x[1]*x[1] + 2*c13*x[0]*x[1]*x[1]*y[0] + c12*x[1]*x[1]*y[0]*y[0] + 2*c31*x[0]*x[0]*x[1]*y[1] + 4*c33*x[0]*x[1]*y[0]*y[1] 
			+ 2*c32*x[1]*y[0]*y[0]*y[1] + c21*x[0]*x[0]*y[1]*y[1] + 2*c23*x[0]*y[0]*y[1]*y[1] + c22*y[0]*y[0]*y[1]*y[1];

		mod(1,2) += c11*x[1]*x[1]*x[2]*x[2] + 2*c31*x[1]*x[2]*x[2]*y[1] + c21*x[2]*x[2]*y[1]*y[1] + 2*c13*x[1]*x[1]*x[2]*y[2] + 4*c33*x[1]*x[2]*y[1]*y[2] 
			+ 2*c23*x[2]*y[1]*y[1]*y[2] + c12*x[1]*x[1]*y[2]*y[2] + 2*c32*x[1]*y[1]*y[2]*y[2] + c22*y[1]*y[1]*y[2]*y[2];

		mod(1,3) += c11*x[1]*x[1]*x[1]*x[2] + c13*x[1]*x[1]*x[2]*y[1] + 2*c31*x[1]*x[1]*x[2]*y[1] + c21*x[1]*x[2]*y[1]*y[1] + 2*c33*x[1]*x[2]*y[1]*y[1] + c23*x[2]*y[1]*y[1]*y[1] 
			+ c13*x[1]*x[1]*x[1]*y[2] + c12*x[1]*x[1]*y[1]*y[2] + 2*c33*x[1]*x[1]*y[1]*y[2] + c23*x[1]*y[1]*y[1]*y[2] + 2*c32*x[1]*y[1]*y[1]*y[2] + c22*y[1]*y[1]*y[1]*y[2];

		mod(1,4) += c11*x[0]*x[1]*x[1]*x[2] + c13*x[1]*x[1]*x[2]*y[0] + 2*c31*x[0]*x[1]*x[2]*y[1] + 2*c33*x[1]*x[2]*y[0]*y[1] + c21*x[0]*x[2]*y[1]*y[1] + c23*x[2]*y[0]*y[1]*y[1] 
			+ c13*x[0]*x[1]*x[1]*y[2] + c12*x[1]*x[1]*y[0]*y[2] + 2*c33*x[0]*x[1]*y[1]*y[2] + 2*c32*x[1]*y[0]*y[1]*y[2] + c23*x[0]*y[1]*y[1]*y[2] + c22*y[0]*y[1]*y[1]*y[2];
	
		mod(1,5) += c11*x[0]*x[1]*x[1]*x[1] + c13*x[1]*x[1]*x[1]*y[0] + c13*x[0]*x[1]*x[1]*y[1] + 2*c31*x[0]*x[1]*x[1]*y[1] + c12*x[1]*x[1]*y[0]*y[1] + 2*c33*x[1]*x[1]*y[0]*y[1] 
			+ c21*x[0]*x[1]*y[1]*y[1] + 2*c33*x[0]*x[1]*y[1]*y[1] + c23*x[1]*y[0]*y[1]*y[1] + 2*c32*x[1]*y[0]*y[1]*y[1] + c23*x[0]*y[1]*y[1]*y[1] + c22*y[0]*y[1]*y[1]*y[1];


		mod(2,0) += c11*x[0]*x[0]*x[2]*x[2] + 2*c13*x[0]*x[2]*x[2]*y[0] + c12*x[2]*x[2]*y[0]*y[0] + 2*c31*x[0]*x[0]*x[2]*y[2] + 4*c33*x[0]*x[2]*y[0]*y[2] 
			+ 2*c32*x[2]*y[0]*y[0]*y[2] + c21*x[0]*x[0]*y[2]*y[2] + 2*c23*x[0]*y[0]*y[2]*y[2] + c22*y[0]*y[0]*y[2]*y[2];

		mod(2,1) += c11*x[1]*x[1]*x[2]*x[2] + 2*c13*x[1]*x[2]*x[2]*y[1] + c12*x[2]*x[2]*y[1]*y[1] + 2*c31*x[1]*x[1]*x[2]*y[2] + 4*c33*x[1]*x[2]*y[1]*y[2] 
			+ 2*c32*x[2]*y[1]*y[1]*y[2] + c21*x[1]*x[1]*y[2]*y[2] + 2*c23*x[1]*y[1]*y[2]*y[2] + c22*y[1]*y[1]*y[2]*y[2];

		mod(2,3) += c11*x[1]*x[2]*x[2]*x[2] + c13*x[2]*x[2]*x[2]*y[1] + c13*x[1]*x[2]*x[2]*y[2] + 2*c31*x[1]*x[2]*x[2]*y[2] + c12*x[2]*x[2]*y[1]*y[2] + 2*c33*x[2]*x[2]*y[1]*y[2] 
			+ c21*x[1]*x[2]*y[2]*y[2] + 2*c33*x[1]*x[2]*y[2]*y[2] + c23*x[2]*y[1]*y[2]*y[2] + 2*c32*x[2]*y[1]*y[2]*y[2] + c23*x[1]*y[2]*y[2]*y[2] + c22*y[1]*y[2]*y[2]*y[2];

		mod(2,4) += c11*x[0]*x[2]*x[2]*x[2] + c13*x[2]*x[2]*x[2]*y[0] + c13*x[0]*x[2]*x[2]*y[2] + 2*c31*x[0]*x[2]*x[2]*y[2] + c12*x[2]*x[2]*y[0]*y[2] + 2*c33*x[2]*x[2]*y[0]*y[2] 
			+ c21*x[0]*x[2]*y[2]*y[2] + 2*c33*x[0]*x[2]*y[2]*y[2] + c23*x[2]*y[0]*y[2]*y[2] + 2*c32*x[2]*y[0]*y[2]*y[2] + c23*x[0]*y[2]*y[2]*y[2] + c22*y[0]*y[2]*y[2]*y[2];

		mod(2,5) += c11*x[0]*x[1]*x[2]*x[2] + c13*x[1]*x[2]*x[2]*y[0] + c13*x[0]*x[2]*x[2]*y[1] + c12*x[2]*x[2]*y[0]*y[1] + 2*c31*x[0]*x[1]*x[2]*y[2] + 2*c33*x[1]*x[2]*y[0]*y[2] 
			+ 2*c33*x[0]*x[2]*y[1]*y[2] + 2*c32*x[2]*y[0]*y[1]*y[2] + c21*x[0]*x[1]*y[2]*y[2] + c23*x[1]*y[0]*y[2]*y[2] + c23*x[0]*y[1]*y[2]*y[2] + c22*y[0]*y[1]*y[2]*y[2];


		mod(3,0) += c11*x[0]*x[0]*x[1]*x[2] + 2*c13*x[0]*x[1]*x[2]*y[0] + c12*x[1]*x[2]*y[0]*y[0] + c31*x[0]*x[0]*x[2]*y[1] + 2*c33*x[0]*x[2]*y[0]*y[1] + c32*x[2]*y[0]*y[0]*y[1] 
			+ c31*x[0]*x[0]*x[1]*y[2] + 2*c33*x[0]*x[1]*y[0]*y[2] + c32*x[1]*y[0]*y[0]*y[2] + c21*x[0]*x[0]*y[1]*y[2] + 2*c23*x[0]*y[0]*y[1]*y[2] + c22*y[0]*y[0]*y[1]*y[2];

		mod(3,1) += c11*x[1]*x[1]*x[1]*x[2] + 2*c13*x[1]*x[1]*x[2]*y[1] + c31*x[1]*x[1]*x[2]*y[1] + c12*x[1]*x[2]*y[1]*y[1] + 2*c33*x[1]*x[2]*y[1]*y[1] + c32*x[2]*y[1]*y[1]*y[1] 
			+ c31*x[1]*x[1]*x[1]*y[2] + c21*x[1]*x[1]*y[1]*y[2] + 2*c33*x[1]*x[1]*y[1]*y[2] + 2*c23*x[1]*y[1]*y[1]*y[2] + c32*x[1]*y[1]*y[1]*y[2] + c22*y[1]*y[1]*y[1]*y[2];

		mod(3,2) += c11*x[1]*x[2]*x[2]*x[2] + c31*x[2]*x[2]*x[2]*y[1] + 2*c13*x[1]*x[2]*x[2]*y[2] + c31*x[1]*x[2]*x[2]*y[2] + c21*x[2]*x[2]*y[1]*y[2] + 2*c33*x[2]*x[2]*y[1]*y[2] 
			+ c12*x[1]*x[2]*y[2]*y[2] + 2*c33*x[1]*x[2]*y[2]*y[2] + 2*c23*x[2]*y[1]*y[2]*y[2] + c32*x[2]*y[1]*y[2]*y[2] + c32*x[1]*y[2]*y[2]*y[2] + c22*y[1]*y[2]*y[2]*y[2];

		mod(3,4) += c11*x[0]*x[1]*x[2]*x[2] + c13*x[1]*x[2]*x[2]*y[0] + c31*x[0]*x[2]*x[2]*y[1] + c33*x[2]*x[2]*y[0]*y[1] + c13*x[0]*x[1]*x[2]*y[2] + c31*x[0]*x[1]*x[2]*y[2] 
			+ c12*x[1]*x[2]*y[0]*y[2] + c33*x[1]*x[2]*y[0]*y[2] + c21*x[0]*x[2]*y[1]*y[2] + c33*x[0]*x[2]*y[1]*y[2] + c23*x[2]*y[0]*y[1]*y[2] + c32*x[2]*y[0]*y[1]*y[2] + c33*x[0]*x[1]*y[2]*y[2] 
			+ c32*x[1]*y[0]*y[2]*y[2] + c23*x[0]*y[1]*y[2]*y[2] + c22*y[0]*y[1]*y[2]*y[2];

		mod(3,5) += c11*x[0]*x[1]*x[1]*x[2] + c13*x[1]*x[1]*x[2]*y[0] + c13*x[0]*x[1]*x[2]*y[1] + c31*x[0]*x[1]*x[2]*y[1] + c12*x[1]*x[2]*y[0]*y[1] + c33*x[1]*x[2]*y[0]*y[1] 
			+ c33*x[0]*x[2]*y[1]*y[1] + c32*x[2]*y[0]*y[1]*y[1] + c31*x[0]*x[1]*x[1]*y[2] + c33*x[1]*x[1]*y[0]*y[2] + c21*x[0]*x[1]*y[1]*y[2] + c33*x[0]*x[1]*y[1]*y[2] 
			+ c23*x[1]*y[0]*y[1]*y[2] + c32*x[1]*y[0]*y[1]*y[2] + c23*x[0]*y[1]*y[1]*y[2] + c22*y[0]*y[1]*y[1]*y[2];


		mod(4,0) += c11*x[0]*x[0]*x[0]*x[2] + 2*c13*x[0]*x[0]*x[2]*y[0] + c31*x[0]*x[0]*x[2]*y[0] + c12*x[0]*x[2]*y[0]*y[0] + 2*c33*x[0]*x[2]*y[0]*y[0] 
			+ c32*x[2]*y[0]*y[0]*y[0] + c31*x[0]*x[0]*x[0]*y[2] + c21*x[0]*x[0]*y[0]*y[2] + 2*c33*x[0]*x[0]*y[0]*y[2] + 2*c23*x[0]*y[0]*y[0]*y[2] + c32*x[0]*y[0]*y[0]*y[2] 
			+ c22*y[0]*y[0]*y[0]*y[2];

		mod(4,1) += c11*x[0]*x[1]*x[1]*x[2] + c31*x[1]*x[1]*x[2]*y[0] + 2*c13*x[0]*x[1]*x[2]*y[1] + 2*c33*x[1]*x[2]*y[0]*y[1] + c12*x[0]*x[2]*y[1]*y[1] + c32*x[2]*y[0]*y[1]*y[1] 
			+ c31*x[0]*x[1]*x[1]*y[2] + c21*x[1]*x[1]*y[0]*y[2] + 2*c33*x[0]*x[1]*y[1]*y[2] + 2*c23*x[1]*y[0]*y[1]*y[2] + c32*x[0]*y[1]*y[1]*y[2] + c22*y[0]*y[1]*y[1]*y[2];

		mod(4,2) += c11*x[0]*x[2]*x[2]*x[2] + c31*x[2]*x[2]*x[2]*y[0] + 2*c13*x[0]*x[2]*x[2]*y[2] + c31*x[0]*x[2]*x[2]*y[2] + c21*x[2]*x[2]*y[0]*y[2] + 2*c33*x[2]*x[2]*y[0]*y[2] 
			+ c12*x[0]*x[2]*y[2]*y[2] + 2*c33*x[0]*x[2]*y[2]*y[2] + 2*c23*x[2]*y[0]*y[2]*y[2] + c32*x[2]*y[0]*y[2]*y[2] + c32*x[0]*y[2]*y[2]*y[2] + c22*y[0]*y[2]*y[2]*y[2];

		mod(4,3) += c11*x[0]*x[1]*x[2]*x[2] + c31*x[1]*x[2]*x[2]*y[0] + c13*x[0]*x[2]*x[2]*y[1] + c33*x[2]*x[2]*y[0]*y[1] + c13*x[0]*x[1]*x[2]*y[2] + c31*x[0]*x[1]*x[2]*y[2] 
			+ c21*x[1]*x[2]*y[0]*y[2] + c33*x[1]*x[2]*y[0]*y[2] + c12*x[0]*x[2]*y[1]*y[2] + c33*x[0]*x[2]*y[1]*y[2] + c23*x[2]*y[0]*y[1]*y[2] + c32*x[2]*y[0]*y[1]*y[2] + c33*x[0]*x[1]*y[2]*y[2] 
			+ c23*x[1]*y[0]*y[2]*y[2] + c32*x[0]*y[1]*y[2]*y[2] + c22*y[0]*y[1]*y[2]*y[2];

		mod(4,5) += c11*x[0]*x[0]*x[1]*x[2] + c13*x[0]*x[1]*x[2]*y[0] + c31*x[0]*x[1]*x[2]*y[0] + c33*x[1]*x[2]*y[0]*y[0] + c13*x[0]*x[0]*x[2]*y[1] + c12*x[0]*x[2]*y[0]*y[1] 
			+ c33*x[0]*x[2]*y[0]*y[1] + c32*x[2]*y[0]*y[0]*y[1] + c31*x[0]*x[0]*x[1]*y[2] + c21*x[0]*x[1]*y[0]*y[2] + c33*x[0]*x[1]*y[0]*y[2] + c23*x[1]*y[0]*y[0]*y[2] 
			+ c33*x[0]*x[0]*y[1]*y[2] + c23*x[0]*y[0]*y[1]*y[2] + c32*x[0]*y[0]*y[1]*y[2] + c22*y[0]*y[0]*y[1]*y[2];
	
		mod(5,0) += c11*x[0]*x[0]*x[0]*x[1] + 2*c13*x[0]*x[0]*x[1]*y[0] + c31*x[0]*x[0]*x[1]*y[0] + c12*x[0]*x[1]*y[0]*y[0] + 2*c33*x[0]*x[1]*y[0]*y[0] + c32*x[1]*y[0]*y[0]*y[0] 
			+ c31*x[0]*x[0]*x[0]*y[1] + c21*x[0]*x[0]*y[0]*y[1] + 2*c33*x[0]*x[0]*y[0]*y[1] + 2*c23*x[0]*y[0]*y[0]*y[1] + c32*x[0]*y[0]*y[0]*y[1] + c22*y[0]*y[0]*y[0]*y[1];

		mod(5,1) += c11*x[0]*x[1]*x[1]*x[1] + c31*x[1]*x[1]*x[1]*y[0] + 2*c13*x[0]*x[1]*x[1]*y[1] + c31*x[0]*x[1]*x[1]*y[1] + c21*x[1]*x[1]*y[0]*y[1] + 2*c33*x[1]*x[1]*y[0]*y[1] 
			+ c12*x[0]*x[1]*y[1]*y[1] + 2*c33*x[0]*x[1]*y[1]*y[1] + 2*c23*x[1]*y[0]*y[1]*y[1] + c32*x[1]*y[0]*y[1]*y[1] + c32*x[0]*y[1]*y[1]*y[1] + c22*y[0]*y[1]*y[1]*y[1];

		mod(5,2) += c11*x[0]*x[1]*x[2]*x[2] + c31*x[1]*x[2]*x[2]*y[0] + c31*x[0]*x[2]*x[2]*y[1] + c21*x[2]*x[2]*y[0]*y[1] + 2*c13*x[0]*x[1]*x[2]*y[2] + 2*c33*x[1]*x[2]*y[0]*y[2] 
			+ 2*c33*x[0]*x[2]*y[1]*y[2] + 2*c23*x[2]*y[0]*y[1]*y[2] + c12*x[0]*x[1]*y[2]*y[2] + c32*x[1]*y[0]*y[2]*y[2] + c32*x[0]*y[1]*y[2]*y[2] + c22*y[0]*y[1]*y[2]*y[2];

		mod(5,3) += c11*x[0]*x[1]*x[1]*x[2] + c31*x[1]*x[1]*x[2]*y[0] + c13*x[0]*x[1]*x[2]*y[1] + c31*x[0]*x[1]*x[2]*y[1] + c21*x[1]*x[2]*y[0]*y[1] + c33*x[1]*x[2]*y[0]*y[1] + c33*x[0]*x[2]*y[1]*y[1] 
			+ c23*x[2]*y[0]*y[1]*y[1] + c13*x[0]*x[1]*x[1]*y[2] + c33*x[1]*x[1]*y[0]*y[2] + c12*x[0]*x[1]*y[1]*y[2] + c33*x[0]*x[1]*y[1]*y[2] + c23*x[1]*y[0]*y[1]*y[2] + c32*x[1]*y[0]*y[1]*y[2] 
			+ c32*x[0]*y[1]*y[1]*y[2] + c22*y[0]*y[1]*y[1]*y[2];

		mod(5,4) += c11*x[0]*x[0]*x[1]*x[2] + c13*x[0]*x[1]*x[2]*y[0] + c31*x[0]*x[1]*x[2]*y[0] + c33*x[1]*x[2]*y[0]*y[0] + c31*x[0]*x[0]*x[2]*y[1] + c21*x[0]*x[2]*y[0]*y[1] + c33*x[0]*x[2]*y[0]*y[1] 
			+ c23*x[2]*y[0]*y[0]*y[1] + c13*x[0]*x[0]*x[1]*y[2] + c12*x[0]*x[1]*y[0]*y[2] + c33*x[0]*x[1]*y[0]*y[2] + c32*x[1]*y[0]*y[0]*y[2] + c33*x[0]*x[0]*y[1]*y[2] 
			+ c23*x[0]*y[0]*y[1]*y[2] + c32*x[0]*y[0]*y[1]*y[2] + c22*y[0]*y[0]*y[1]*y[2];
	}
}
