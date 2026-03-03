/* RGSplitDEMatT.h
 *
 * Staggered dielectric elastomer mechanical material.
 *
 * Extends RGSplitT2 (the "RG_split_general" material) to include
 * Maxwell (electrostatic) stress when an electric field is present.
 * Designed for staggered coupling: the element provides the current
 * IP's electric field via SetIPElectricField() before calling s_ij().
 * c_ijkl() is left unchanged (pure mechanical tangent), which is
 * correct for a staggered scheme where E is frozen during the
 * mechanical Newton-Raphson solve.
 *
 * XML tag: "RG_split_DE"
 * XML parameters (in addition to RGSplitT2 params):
 *   epsilon  -- dielectric permittivity (default 1.0)
 *
 * The Maxwell stress added to the Cauchy stress is:
 *   sigma_ij^Maxwell = epsilon * (E_i E_j - 0.5 |E|^2 delta_ij)
 *
 * Pluggable potentials (NeoHookean, ArrudaBoyce, etc.) are selected
 * exactly as in RGSplitT2 via the <rg_eq_potential> sublist.
 */
#ifndef _RG_SPLIT_DE_MAT_T_H_
#define _RG_SPLIT_DE_MAT_T_H_

#include "RGSplitT2.h"
#include "IElectricallyCouplable.h"
#include "ParameterInterfaceT.h"
#include "dSymMatrixT.h"
#include "dArrayT.h"

namespace Tahoe {

class RGSplitDEMatT : public RGSplitT2, public IElectricallyCouplable
{
public:

    /** constructor */
    RGSplitDEMatT(void);

    /** destructor */
    virtual ~RGSplitDEMatT(void) {}

    /** IElectricallyCouplable: store E at current IP before s_ij() */
    virtual void SetIPElectricField(const dArrayT& E);

    /** Cauchy stress = mechanical (RGSplitT2) + Maxwell stress */
    virtual const dSymMatrixT& s_ij(void);

    /** \name ParameterInterfaceT */
    /*@{*/
    virtual void DefineParameters(ParameterListT& list) const;
    virtual void TakeParameterList(const ParameterListT& list);
    /*@}*/

private:

    /** dielectric permittivity */
    double fEpsilon;

    /** electric field at current integration point (set by element) */
    dArrayT fE_current;

    /** total stress workspace (mechanical + Maxwell) */
    dSymMatrixT fTotalStress;
};

} // namespace Tahoe
#endif /* _RG_SPLIT_DE_MAT_T_H_ */
