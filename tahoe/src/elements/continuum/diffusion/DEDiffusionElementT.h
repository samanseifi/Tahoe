/* DEDiffusionElementT.h
 *
 * Clean dielectric elastomer electric element for staggered coupling.
 *
 * Extends DiffusionElementT to implement the electric problem for a
 * dielectric elastomer in the deformed configuration.  All DE-specific
 * formulas are isolated here; DiffusionElementT remains a general solver.
 *
 * Physics:
 *   Electric field:        E_i = -dphi/dx_i   (phi = scalar potential)
 *   Electric displacement: d_i = epsilon * J * C^{-1}_{IJ} * E_J
 *                               = epsilon * b_ij * E_j
 *   Residual:              R_A = integral( dNa/dxi * d_i ) dV
 *   Tangent:               K_AB = integral( dNa/dxi * epsilon*J*C^{-1}_{ij} * dNb/dxj ) dV
 *
 * where J = det(F), C = F^T F, and F is the deformation gradient
 * obtained from the mechanical displacement field.
 *
 * XML tag: "de_diffusion"
 * XML attribute: epsilon (dielectric permittivity, default 1.0)
 *
 * The mechanical displacement field ("displacement") must exist in the
 * simulation; the element reads it automatically, exactly as the base
 * DiffusionElementT does.
 */
#ifndef _DE_DIFFUSION_ELEMENT_T_H_
#define _DE_DIFFUSION_ELEMENT_T_H_

#include "DiffusionElementT.h"

namespace Tahoe {

class DEDiffusionElementT : public DiffusionElementT
{
public:

    /** constructor */
    DEDiffusionElementT(const ElementSupportT& support);

    /** \name ParameterInterfaceT */
    /*@{*/
    virtual void DefineParameters(ParameterListT& list) const;
    virtual void TakeParameterList(const ParameterListT& list);
    /*@}*/

protected:

    /** tangent: K_AB = integral( B^T * (epsilon * J * C^{-1}) * B ) dV */
    virtual void FormStiffness(double constK);

    /** residual: R_A = integral( B^T * d_i ) dV
     *            d_i = epsilon * J * C^{-1} * E,  E = -grad(phi) */
    virtual void FormKd(double constK);

private:

    /** dielectric permittivity */
    double fEpsilon;
};

} // namespace Tahoe
#endif /* _DE_DIFFUSION_ELEMENT_T_H_ */
