#include "FSDE_incQ1P0.h"

#include <iostream>
#include <fstream>
using namespace std;

namespace Tahoe {

  inline FSDEMatQ1P0T::FSDEMatQ1P0T() :
    ParameterInterfaceT("Finite Strain Dielectric Elastomer Q1P0"),
        fFSDEMatSupportQ1P0(0)
  {
    SetName(FSDEMatQ1P0T::Name);
    Initialize();
  }

  // Set electrical permittivity
  inline void FSDEMatQ1P0T::SetElectricPermittivity(double epsilon)
  {	
    fElectricPermittivity = epsilon;
  }

  // Get electrical permittivity
  inline double FSDEMatQ1P0T::GetElectricPermittivity() const
  {
    return fElectricPermittivity;
  }

  //
  inline void FSDEMatQ1P0T::SetFSDEMatSupportQ1P0(
      const FSDEMatSupportQ1P0T* support)
  {
    fFSDEMatSupportQ1P0 = support;
  }

  //
  inline const dArrayT FSDEMatQ1P0T::ElectricField()
  {
    fElectricField = fFSDEMatSupportQ1P0->ElectricField();
	/*---------- Prescribing Electric Field ---------- */
	//fElectricField[0] = 0.0;
    //fElectricField[1] = -0.1;
    //fElectricField[2] = 0.0;
     /*-----------------------------------------------*/
    return fElectricField;
  }

  //
  inline const dArrayT FSDEMatQ1P0T::ElectricField(int ip)
  {
    fElectricField = fFSDEMatSupportQ1P0->ElectricField(ip);
	/*---------- Prescribing Electric Field ----------*/
 	//fElectricField[0] = 0.0;
    //fElectricField[1] = -0.1;
    //fElectricField[2] = 0.0;
    /*-----------------------------------------------*/
    return fElectricField;
  }

  /* ---------- Prescribed Deformation Tensor ---------- */
  inline const dMatrixT FSDEMatQ1P0T::DeformationMatrix()
  {
	  dMatrixT F = F_mechanical();
	  //dMatrixT F(3);
	  //F(0,0) = 1.05; F(1,0) = 0.0;  F(2,0) = 0.0;
	  //F(0,1) = 0.0;  F(1,1) = 0.95; F(2,1) = 0.0;
	  //F(0,2) = 0.0;  F(1,2) = 0.0;  F(2,2) = 1.0;
	  return F;
  }

  //
  inline const dMatrixT FSDEMatQ1P0T::RightCauchyGreenDeformation()
  {
    const dMatrixT& F = F_mechanical();
    dMatrixT FTF(3);
    FTF.MultATB(F, F);

    return FTF;
  }

  // material energy density
  inline double FSDEMatQ1P0T::StrainEnergyDensity()
  {
	  return 0.0;
  }

  // material mechanical tangent modulus
  inline const dMatrixT&
  FSDEMatQ1P0T::C_IJKL()
  {
    const dMatrixT& C = RightCauchyGreenDeformation();
    const dArrayT& E = ElectricField();
    const dMatrixT& F = F_mechanical();

	double I1 = C(0,0)+C(1,1)+C(2,2);
	double J = F.Det();

	/* call C function for mechanical part of tangent modulus */
 	mech_tanmod_q1p0(fParams.Pointer(), E.Pointer(), C.Pointer(), F.Pointer(), J, I1, fTangentMechanical.Pointer()); 
 	me_tanmod_q1p0(fParams.Pointer(), E.Pointer(), C.Pointer(), F.Pointer(), J, fTangentMechanicalElec.Pointer());
 	fTangentMechanical+=fTangentMechanicalElec;

 	/* -------------- Writing into a file -------------
 	ofstream myCIJKL;
 	myCIJKL.open("C_IJKL.txt");
 	for (int i = 0; i < fTangentMechanical.Rows(); i++) {
 		for (int j = 0; j < fTangentMechanical.Cols(); j++) {
 			myCIJKL << fTangentMechanical(i,j) << " ";
 		}
 		myCIJKL << endl;
 	}
 	myCIJKL.close();
 	 ------------------------------------------------- */


    return fTangentMechanical;
  }

  // Second Piola-Kirchhoff stress (mechanical)
  inline const dSymMatrixT&
  FSDEMatQ1P0T::S_IJ()
  {
    const dMatrixT& C = RightCauchyGreenDeformation();
	const dArrayT& E = ElectricField();
    const dMatrixT& F = F_mechanical();
	
	dMatrixT stress_temp(3);
	dMatrixT stress_temp2(3);
	double I1 = C(0,0)+C(1,1)+C(2,2);
	double J = C.Det();
	J = sqrt(J);
	
	/* call C function for mechanical part of PK2 stress */
 	mech_pk2_q1p0(fParams.Pointer(), E.Pointer(), C.Pointer(), F.Pointer(), J, I1, stress_temp.Pointer()); 
	me_pk2_q1p0(fParams.Pointer(), E.Pointer(), C.Pointer(), F.Pointer(), J, stress_temp2.Pointer());
	stress_temp+=stress_temp2;
	
	fStress.FromMatrix(stress_temp);

 	/* -------------- Writing into a file -------------
 	ofstream mySIJ;
 	mySIJ.open("S_IJ.txt");
 	for (int i = 0; i < fStress.Rows(); i++) {
 		for (int j = 0; j < fStress.Cols(); j++) {
 			mySIJ << fStress(i,j) << " ";
 		}
 		mySIJ << endl;
 	}
 	mySIJ.close();
 	 ------------------------------------------------- */

    return fStress;
  }

  // material electromechanical tangent modulus
  inline const dMatrixT&
  FSDEMatQ1P0T::E_IJK()
  {
    const dMatrixT& C = RightCauchyGreenDeformation();
	const dArrayT& E = ElectricField();
    const dMatrixT& F = F_mechanical();	
	double J = C.Det();
	J = sqrt(J);

	/* call C function for electromechanical tangent modulus */
 	me_mixedmodulus_q1p0(fParams.Pointer(), E.Pointer(), C.Pointer(), F.Pointer(), J, fTangentElectromechanical.Pointer());
 
 	/* -------------- Writing into a file -------------
 	ofstream myEIJK;
 	myEIJK.open("E_IJK.txt");
 	for (int i = 0; i < fTangentElectromechanical.Rows(); i++) {
 		for (int j = 0; j < fTangentElectromechanical.Cols(); j++) {
 			myEIJK << fTangentElectromechanical(i,j) << " ";
 		}
 		myEIJK << endl;
 	}
 	myEIJK.close();
 	 ------------------------------------------------- */


    return fTangentElectromechanical;
  }

  // spatial electromechanical tangent modulus
  inline const dMatrixT&
  FSDEMatQ1P0T::e_ijk()
  {
    const dMatrixT& C = RightCauchyGreenDeformation();
	const dArrayT& E = ElectricField();
	double J = C.Det();
	J = sqrt(J);  
    const dMatrixT& F = F_mechanical();
	
	/* call C function for (spatial) electromechanical tangent modulus */
 	me_mixedmodulus_q1p0spatial(fParams.Pointer(), E.Pointer(), C.Pointer(), F.Pointer(), J, fTangentElectromechanicalSpatial.Pointer());
 	
 	fTangentElectromechanicalSpatial /= J;

 	/* -------------- Writing into a file -------------
 	ofstream myeijk;
 	myeijk.open("e_ijk_2.txt");
 	for (int i = 0; i < fTangentElectromechanicalSpatial.Rows(); i++) {
 		for (int j = 0; j < fTangentElectromechanicalSpatial.Cols(); j++) {
 			myeijk << fTangentElectromechanicalSpatial(i,j) << " ";
 		}
 		myeijk << endl;
 	}
 	myeijk.close();
 	 ------------------------------------------------- */

    return fTangentElectromechanicalSpatial;
  }

  // material electric tangent modulus
  inline const dMatrixT&
  FSDEMatQ1P0T::B_IJ()
  {
    const dMatrixT& C = RightCauchyGreenDeformation();
	double J = C.Det();
	J = sqrt(J);

	dMatrixT Cinv(3);
	Cinv.Inverse(C);
	fTangentElectrical = Cinv;
	fTangentElectrical *= fElectricPermittivity;
	fTangentElectrical *= J;

 	/* -------------- Writing into a file -------------
 	ofstream myBIJ;
 	myBIJ.open("B_IJ.txt");
 	for (int i = 0; i < fTangentElectrical.Rows(); i++) {
 		for (int j = 0; j < fTangentElectrical.Cols(); j++) {
 			myBIJ << fTangentElectrical(i,j) << " ";
 		}
 		myBIJ << endl;
 	}
 	myBIJ.close();
 	 ------------------------------------------------- */

    return fTangentElectrical;
  }

  // spatial electric tangent modulus
  inline const dMatrixT&
  FSDEMatQ1P0T::b_ij()
  {
    const dMatrixT& F = F_mechanical();
    const dMatrixT& C = RightCauchyGreenDeformation();    
    const double J = F.Det();
	
	// repeat B_IJ first, then push forward
	dMatrixT Cinv(3), bij(3);
	Cinv.Inverse(C);
	bij = Cinv;
	bij *= fElectricPermittivity;
	bij *= J;
	
    // prevent aliasing
    const dMatrixT b = B_IJ();
    fTangentElectrical.MultABCT(F, b, F);
    fTangentElectrical /= J;

 	/* -------------- Writing into a file -------------
 	ofstream mybij;
 	mybij.open("b_ij_2.txt");
 	for (int i = 0; i < fTangentElectrical.Rows(); i++) {
 		for (int j = 0; j < fTangentElectrical.Cols(); j++) {
 			mybij << fTangentElectrical(i,j) << " ";
 		}
 		mybij << endl;
 	}
 	mybij.close();
 	 ------------------------------------------------- */

    return fTangentElectrical;
  }

  // Electric displacement 
  inline const dArrayT&
  FSDEMatQ1P0T::D_I()
  {
  	const dMatrixT& C = RightCauchyGreenDeformation();
  	const dArrayT& E = ElectricField();
    const dMatrixT& F = F_mechanical();
	/* -------------- Writing into a file -------------
	ofstream myE;
	myE.open("E.txt");
	for (int i = 0; i < E.Length(); i++) {
		myE << E[i] << endl;
	}
	myE.close();
 	 ------------------------------------------------- */
  	
  	double J = C.Det();
  	J = sqrt(J);
  
	/* call C function for electric stress (i.e. electric displacement D_{I}) */
 	elec_pk2_q1p0(fParams.Pointer(), E.Pointer(),  
 		C.Pointer(), F.Pointer(), J, fElectricDisplacement.Pointer()); 

 	/* -------------- Writing into a file -------------
 	ofstream myDI;
 	myDI.open("D_I.txt");
 	for (int i = 0; i < fElectricDisplacement.Length(); i++) {
 		myDI << fElectricDisplacement[i] << " ";
 	}
 	myDI.close();
 	 ------------------------------------------------- */

  	return fElectricDisplacement;
  }

  // spatial electric tangent modulus
  inline const dArrayT&
  FSDEMatQ1P0T::d_i()
  {
    const dMatrixT& F = F_mechanical();
    const double J = F.Det();
	
    // prevent aliasing
    const dArrayT D = D_I();
	F.Multx(D, fElectricDisplacement);
	fElectricDisplacement /= J;

 	/* -------------- Writing into a file -------------
 	ofstream mydi;
 	mydi.open("D_I_2.txt");
 	for (int i = 0; i < fElectricDisplacement.Length(); i++) {
 		mydi << fElectricDisplacement[i] << " ";
 	}
 	mydi.close();
 	 ------------------------------------------------- */

    return fElectricDisplacement;		// need to divide by J
  }

  // Electric field
  inline const dArrayT&
  FSDEMatQ1P0T::E_I()
  {
	/* -------------- Writing into a file -------------
	ofstream myEI;
	myEI.open("E_I.txt");
	for (int i = 0; i < fElectricField.Length(); i++) {
		myEI << fElectricField[i] << endl;
	}
	myEI.close();
 	 ------------------------------------------------- */

    return fElectricField;
  }

  // spatial tangent modulus
  inline const dMatrixT&
  FSDEMatQ1P0T::c_ijkl()
  {
//     const dMatrixT& F = F_mechanical();
//     const double J = F.Det();
// 
//     // prevent aliasing
//     const dMatrixT CIJKL = C_IJKL();
//     fTangentMechanical.SetToScaled(1.0 / J, PushForward(F, CIJKL));

	fTangentMechanical = FSSolidMatT::c_ijkl();

 	/* -------------- Writing into a file -------------
 	ofstream mycijkl;
 	mycijkl.open("c_ijkl_2.txt");
 	for (int i = 0; i < fTangentMechanical.Rows(); i++) {
 		for (int j = 0; j < fTangentMechanical.Cols(); j++) {
 			mycijkl << fTangentMechanical(i,j) << " ";
 		}
 		mycijkl << endl;
 	}
 	mycijkl.close();
 	 ------------------------------------------------- */

    return fTangentMechanical;

  }

  // Cauchy stress
  inline const dSymMatrixT&
  FSDEMatQ1P0T::s_ij()
  {
    const dMatrixT& F = F_mechanical();
    const double J = F.Det();
	
    // prevent aliasing
    const dSymMatrixT S = S_IJ();

    fStress.SetToScaled(1.0 / J, PushForward(F, S));

 	/* -------------- Writing into a file -------------
 	ofstream mysij;
 	mysij.open("s_ij_2.txt");
 	for (int i = 0; i < fStress.Rows(); i++) {
 		for (int j = 0; j < fStress.Cols(); j++) {
 			mysij << fStress(i,j) << " ";
 		}
 		mysij << endl;
 	}
 	mysij.close();
 	 ------------------------------------------------- */

    return fStress;
  }

  // pressure associated with the last computed stress
  inline double FSDEMatQ1P0T::Pressure() const
  {

    return 0.0;

  }

}	// namespace Tahoe
