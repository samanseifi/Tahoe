#include "ExceptionT.h"
#include "FSDEMat2DT.h"
#include "FSDE_inc2D.h"

namespace Tahoe {

  //
  //
  //
  static const char DE[] = "Dielectric_Elastomer_2D";
  const char* FSDEMat2DT::Name = DE;
  const int kNSD       = 2;
  const int kNumDOF    = 2;
  const int kStressDim =dSymMatrixT::NumValues(kNSD);
  //
  //
  //
  void FSDEMat2DT::Initialize()
  {
    fMu = 0.0;
    fElectricPermittivity = 0.0;
    fNrig = 0.0;
    fLambda = 0.0;
    fGamma = 0.0;
    t_0 = 0.0;
  }

  //
  //
  //
  void FSDEMat2DT::DefineParameters(ParameterListT& list) const
  {
	FSSolidMatT::DefineParameters(list);
	
	list.AddParameter(fMu, "mu");
	list.AddParameter(fElectricPermittivity, "epsilon");
 	list.AddParameter(fNrig, "Nrig");
 	list.AddParameter(fLambda, "lambda");
 	list.AddParameter(fGamma, "gamma");
 	list.AddParameter(t_0, "t_0");   // Ramping up the surface tension time scale
  }

  //
  //
  //
  void FSDEMat2DT::TakeParameterList(const ParameterListT& list)
  {
	FSSolidMatT::TakeParameterList(list);
	fMu = list.GetParameter("mu");
	fElectricPermittivity = list.GetParameter("epsilon");
 	fNrig = list.GetParameter("Nrig");
 	fLambda = list.GetParameter("lambda");
 	fGamma = list.GetParameter("gamma");
 	t_0 = list.GetParameter("t_0");

	/* write into vector to pass to C code for stress/modulus calculations */
	fParams.Dimension(4);
	fParams[0] = fMu;
	fParams[1] = fLambda;
 	fParams[2] = fElectricPermittivity;
 	fParams[3] = fNrig;
	
	/* dimension work space */
	fTangentMechanical.Dimension(kStressDim);
	fTangentMechanicalElec.Dimension(kStressDim);
	fStress.Dimension(kNumDOF);
	fTangentElectrical.Dimension(kNumDOF);
	fTangentElectromechanical.Dimension(kStressDim, kNumDOF);
	fTangentElectromechanicalSpatial.Dimension(kStressDim, kNumDOF);	
	fElectricDisplacement.Dimension(kNumDOF);
	fElectricField.Dimension(kNumDOF);
  }

  //
  // information about subordinate parameter lists
  //
  void FSDEMat2DT::DefineSubs(SubListT& sub_list) const
  {
	FSSolidMatT::DefineSubs(sub_list);
    return;
  }



} //namespace Tahoe
