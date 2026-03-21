/* ExplicitKernelT.h — abstract element topology kernel for explicit dynamics.
 *
 * A kernel encapsulates the element geometry: shape functions, parent domain
 * integration points, Jacobian computation, and spatial derivative mapping.
 * It knows nothing about the material — that separation enables N kernels
 * x M materials = N*M combinations with only N+M implementations.
 *
 * All arrays use SoA layout with second dimension MVSIZ for SIMD vectorization.
 */

#ifndef _EXPLICIT_KERNEL_T_H_
#define _EXPLICIT_KERNEL_T_H_

namespace Tahoe {

class ExplicitKernelT
{
public:

	/** machine-tunable vector block size */
	static const int MVSIZ = 128;

	/** maximum nodes per element across all kernels */
	static const int MAX_NEN = 20;

	/** maximum spatial dimensions */
	static const int MAX_NSD = 3;

	virtual ~ExplicitKernelT(void) {}

	/** number of nodes per element */
	virtual int NodesPerElement(void) const = 0;

	/** number of integration points */
	virtual int NumIP(void) const = 0;

	/** number of spatial dimensions */
	virtual int NumSD(void) const = 0;

	/** Compute spatial derivatives and Jacobian determinant for nel elements
	 *  at integration point ip.
	 *
	 * \param ip     integration point index [0, NumIP())
	 * \param nel    number of elements in batch (<= MVSIZ)
	 * \param xc     current coordinates [MAX_NEN][MVSIZ], only [0..nen-1] used
	 * \param dNdx   output: dN/dx [MAX_NEN][MVSIZ]
	 * \param dNdy   output: dN/dy [MAX_NEN][MVSIZ]
	 * \param dNdz   output: dN/dz [MAX_NEN][MVSIZ] (unused for 2D)
	 * \param detJ   output: Jacobian determinant [MVSIZ]
	 * \param weight output: integration weight for this IP
	 */
	virtual void ComputeIPData(
		int ip, int nel,
		const double xc[][MVSIZ],
		const double yc[][MVSIZ],
		const double zc[][MVSIZ],
		double dNdx[][MVSIZ],
		double dNdy[][MVSIZ],
		double dNdz[][MVSIZ],
		double detJ[],
		double& weight
	) const = 0;
};

} /* namespace Tahoe */
#endif /* _EXPLICIT_KERNEL_T_H_ */
