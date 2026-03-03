/* DEDiffusionElementT.cpp */

#include "DEDiffusionElementT.h"
#include "ShapeFunctionT.h"
#include "ParameterListT.h"
#include "ParameterT.h"
#include "LimitT.h"
#include "ElementMatrixT.h"
#include "dMatrixT.h"
#include "dArrayT.h"

using namespace Tahoe;

/* constructor */
DEDiffusionElementT::DEDiffusionElementT(const ElementSupportT& support)
    : DiffusionElementT(support),
      fEpsilon(1.0)
{
    SetName("de_diffusion");
}

/* add epsilon to parameter list */
void DEDiffusionElementT::DefineParameters(ParameterListT& list) const
{
    /* inherited (field_name, geometry, output codes, etc.) */
    DiffusionElementT::DefineParameters(list);

    ParameterT eps(fEpsilon, "epsilon");
    eps.AddLimit(0.0, LimitT::LowerInclusive);
    eps.SetDefault(1.0);
    list.AddParameter(eps);
}

/* read epsilon */
void DEDiffusionElementT::TakeParameterList(const ParameterListT& list)
{
    /* inherited sets up fF_List, fLocDisplacement, mechanical_coupling, etc. */
    DiffusionElementT::TakeParameterList(list);

    fEpsilon = list.GetParameter("epsilon");
}

/* ---- tangent ----------------------------------------------------------- */
/* K_AB = integral( dNa/dxi * epsilon * J * Cinv_ij * dNb/dxj ) dV        */
void DEDiffusionElementT::FormStiffness(double constK)
{
    dMatrixT::SymmetryFlagT format =
        (fLHS.Format() == ElementMatrixT::kNonSymmetric)
        ? dMatrixT::kWhole : dMatrixT::kUpperOnly;

    const double* Det    = fShapes->IPDets();
    const double* Weight = fShapes->IPWeights();

    int nsd = NumSD();
    dMatrixT C(nsd), Cinv(nsd);

    fShapes->TopIP();
    while (fShapes->NextIP())
    {
        double scale = constK * (*Det++) * (*Weight++);

        /* strain-displacement (gradient) matrix */
        B(fShapes->CurrIP(), fB);

        /* permittivity tensor in deformed config: epsilon * J * C^{-1} */
        const dMatrixT& F = fF_List[CurrIP()];
        C.MultATB(F, F);
        double J = F.Det();
        Cinv.Inverse(C);

        fD.SetToScaled(fEpsilon * J * scale, Cinv);

        /* accumulate: K += B^T * D * B */
        fLHS.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);
    }
}

/* ---- residual ---------------------------------------------------------- */
/* R_A = integral( dNa/dxi * d_i ) dV                                      */
/* d_i = epsilon * J * C^{-1} * E,   E_i = -dphi/dxi                      */
void DEDiffusionElementT::FormKd(double constK)
{
    const double* Det    = fShapes->IPDets();
    const double* Weight = fShapes->IPWeights();

    int nsd = NumSD();
    dMatrixT grad, C(nsd), Cinv(nsd);
    dArrayT E(nsd), di(nsd);

    fShapes->TopIP();
    while (fShapes->NextIP())
    {
        /* electric field E = -grad(phi) at this IP */
        grad.Set(1, nsd, fGradient_list[CurrIP()].Pointer());
        IP_ComputeGradient(fLocDisp, grad);
        for (int i = 0; i < nsd; i++) E[i] = -grad(0, i);

        /* gradient matrix */
        B(fShapes->CurrIP(), fB);

        /* deformation gradient stored by parent's SetGlobalShape */
        const dMatrixT& F = fF_List[CurrIP()];
        C.MultATB(F, F);
        double J = F.Det();
        Cinv.Inverse(C);

        /* electric displacement: d_i = epsilon * J * C^{-1} * E */
        Cinv.Multx(E, di);
        di *= fEpsilon * J;

        /* assemble: fRHS -= constK * w * det * B^T * d_i */
        fB.MultTx(di, fNEEvec);
        fRHS.AddScaled(-constK * (*Weight++) * (*Det++), fNEEvec);
    }
}
