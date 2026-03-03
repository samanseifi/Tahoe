/* IElectricallyCouplable.h
 *
 * Interface for finite-strain solid materials that can accept
 * an externally computed electric field at each integration point.
 * Used by SimoQ1P0 (and its derivatives) to push the current IP's
 * electric field into the material before calling s_ij(), enabling
 * a staggered dielectric-elastomer solve where Maxwell stress is
 * handled at the material level rather than the element level.
 */
#ifndef _I_ELECTRICALLY_COUPLABLE_H_
#define _I_ELECTRICALLY_COUPLABLE_H_

#include "dArrayT.h"

namespace Tahoe {

class IElectricallyCouplable
{
public:
    virtual ~IElectricallyCouplable() {}

    /** Set the electric field E at the current integration point.
     * Must be called before s_ij() so the material can include
     * Maxwell stress in the returned Cauchy stress. */
    virtual void SetIPElectricField(const dArrayT& E) = 0;
};

} // namespace Tahoe
#endif /* _I_ELECTRICALLY_COUPLABLE_H_ */
