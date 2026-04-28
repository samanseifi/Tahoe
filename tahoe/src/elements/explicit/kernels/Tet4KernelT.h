/* Tet4KernelT.h — 4-node constant-strain tetrahedron kernel for 3D explicit
 * dynamics.  Single integration point at the centroid (parametric coords
 * r=s=t=1/4), weight = 1/6 (volume of the unit reference tet).
 *
 * Node ordering follows Tahoe convention (toolbox/.../TetrahedronT.cpp):
 *   node 0:  (r=1, s=0, t=0)
 *   node 1:  (r=0, s=1, t=0)
 *   node 2:  (r=0, s=0, t=0)   ← N_2 = 1 - r - s - t
 *   node 3:  (r=0, s=0, t=1)
 *
 * Limitations (intrinsic to constant-strain Tet4 — see #28 for ANP fix):
 *   - locks under near-incompressibility (J2 plasticity, rubber)
 *   - coarse approximation; converges slowly with mesh refinement
 *   - acceptable for elastic / mildly plastic problems and as a building block
 *     for the average-nodal-pressure tet (LS-DYNA ELFORM=13).
 */

#ifndef _TET4_KERNEL_T_H_
#define _TET4_KERNEL_T_H_

#include "ExplicitKernelT.h"

namespace Tahoe {

class Tet4KernelT : public ExplicitKernelT
{
public:
	int NodesPerElement(void) const { return 4; }
	int NumIP(void) const { return 1; }
	int NumSD(void) const { return 3; }

	void ComputeIPData(
		int ip, int nel,
		const double xc[][MVSIZ],
		const double yc[][MVSIZ],
		const double zc[][MVSIZ],
		double dNdx[][MVSIZ],
		double dNdy[][MVSIZ],
		double dNdz[][MVSIZ],
		double detJ[],
		double& weight
	) const;
};

} /* namespace Tahoe */
#endif /* _TET4_KERNEL_T_H_ */
