/* Hex8KernelT.h — 8-node hexahedral kernel for 3D explicit dynamics. */

#ifndef _HEX8_KERNEL_T_H_
#define _HEX8_KERNEL_T_H_

#include "ExplicitKernelT.h"

namespace Tahoe {

class Hex8KernelT : public ExplicitKernelT
{
public:
	int NodesPerElement(void) const { return 8; }
	int NumIP(void) const { return 8; }
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
#endif /* _HEX8_KERNEL_T_H_ */
