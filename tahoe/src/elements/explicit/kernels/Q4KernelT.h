/* Q4KernelT.h — 4-node quadrilateral kernel for 2D explicit dynamics. */

#ifndef _Q4_KERNEL_T_H_
#define _Q4_KERNEL_T_H_

#include "ExplicitKernelT.h"

namespace Tahoe {

class Q4KernelT : public ExplicitKernelT
{
public:
	int NodesPerElement(void) const { return 4; }
	int NumIP(void) const { return 4; }
	int NumSD(void) const { return 2; }

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
#endif /* _Q4_KERNEL_T_H_ */
