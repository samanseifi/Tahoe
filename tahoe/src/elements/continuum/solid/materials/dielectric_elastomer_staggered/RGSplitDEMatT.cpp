/* RGSplitDEMatT.cpp */

#include "RGSplitDEMatT.h"
#include "ParameterListT.h"
#include "ParameterT.h"
#include "LimitT.h"
#include "dMatrixT.h"
#include "ExceptionT.h"

using namespace Tahoe;

/* constructor */
RGSplitDEMatT::RGSplitDEMatT(void)
    : ParameterInterfaceT("RG_split_DE"),
      RGSplitT2(),
      fEpsilon(1.0)
{}

/* IElectricallyCouplable: store the current IP's electric field */
void RGSplitDEMatT::SetIPElectricField(const dArrayT& E)
{
    fE_current = E;
}

/* Cauchy stress: mechanical + Maxwell */
const dSymMatrixT& RGSplitDEMatT::s_ij(void)
{
    /* mechanical Cauchy stress from RGSplitT2 */
    fTotalStress = RGSplitT2::s_ij();

    /* add Maxwell stress only when |E| > 0 */
    int nsd_e = fE_current.Length();
    double Emag2 = 0.0;
    for (int i = 0; i < nsd_e; i++) Emag2 += fE_current[i] * fE_current[i];
    if (Emag2 > 0.0)
    {
        int nsd = NumSD();
        for (int i = 0; i < nsd; i++)
            for (int j = i; j < nsd; j++)
            {
                double mij = fEpsilon * (fE_current[i] * fE_current[j]
                             - 0.5 * Emag2 * (i == j ? 1.0 : 0.0));
                fTotalStress(i, j) += mij;
            }
    }

    return fTotalStress;
}

/* add epsilon parameter */
void RGSplitDEMatT::DefineParameters(ParameterListT& list) const
{
    /* inherited (includes density, constraint, etc.) */
    RGSplitT2::DefineParameters(list);

    /* dielectric permittivity */
    ParameterT eps(fEpsilon, "epsilon");
    eps.AddLimit(0.0, LimitT::LowerInclusive);
    eps.SetDefault(1.0);
    list.AddParameter(eps);
}

/* read epsilon, then dimension workspaces */
void RGSplitDEMatT::TakeParameterList(const ParameterListT& list)
{
    /* inherited sets up potentials, viscosities, etc. */
    RGSplitT2::TakeParameterList(list);

    fEpsilon = list.GetParameter("epsilon");

    /* dimension electric field and stress workspaces */
    int nsd = NumSD();
    fE_current.Dimension(nsd);
    fE_current = 0.0;
    fTotalStress.Dimension(nsd);
}
