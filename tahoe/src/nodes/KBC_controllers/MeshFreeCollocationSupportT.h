/* MeshFreeCollocationSupportT.h — interface a meshfree element implements so the collocation
 * BC controller (CollocationKBCT) can enforce essential BCs on the PHYSICAL displacement.
 *
 * Because RK/MLS shapes are non-interpolatory, the controller needs, for each constrained node I,
 * the support node ids J and the shape values Phi_J(x_I) evaluated at node I's own location. Any
 * element that can supply these (RKShellT, a future EFG/RKPM solid, ...) becomes usable with the
 * collocation BC without the BC module knowing the element's internals.
 */
#ifndef _MESHFREE_COLLOCATION_SUPPORT_T_H_
#define _MESHFREE_COLLOCATION_SUPPORT_T_H_

#include "iArrayT.h"
#include "RaggedArray2DT.h"

namespace Tahoe {

class MeshFreeCollocationSupportT
{
public:

	virtual ~MeshFreeCollocationSupportT(void) {}

	/** collocation rows for the requested global node ids.
	 * \param nodes   global node ids to provide rows for
	 * \param support returns, per requested node, the support node global ids J
	 * \param phi     returns, per requested node, Phi_J(x_I) aligned with support
	 * \return false if any requested node has no meshfree support (cannot collocate) */
	virtual bool CollocationData(const iArrayT& nodes,
		RaggedArray2DT<int>& support, RaggedArray2DT<double>& phi) const = 0;
};

} /* namespace Tahoe */

#endif /* _MESHFREE_COLLOCATION_SUPPORT_T_H_ */
