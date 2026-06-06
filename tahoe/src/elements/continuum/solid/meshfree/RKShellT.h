/* RKShellT.h — meshfree RKPM Kirchhoff–Love shell element (epic #59)
 *
 * Implements the general-purpose meshfree KL shell of
 *   J. Wang, Y. Bazilevs, "A general-purpose meshfree Kirchhoff–Love shell
 *   formulation", Engineering with Computers 41:1379-1410 (2025).
 *
 * STATUS — scaffold. Registered as element type "meshfree_kl_shell". This class
 * subclasses MeshFreeFSSolidT to reuse the meshfree connectivity / equation-numbering /
 * assembly / material plumbing. The KL-shell physics — PCA local parameterization, the
 * auxiliary-tensor velocity-gradient kinematics, the naturally-stabilized nodal
 * integration, the co-rotational stress update and the sigma33=0 plane-stress loop —
 * have been implemented and validated as standalone kernels (contrib/kl_kinematics,
 * kl_bmatrix, kl_stress, kl_assembly; D2OrthoMLS2DT 2nd-derivative fix #69; CTest
 * "meshfree_*"). Porting those kernels into RHSDriver()/LHSDriver() is in progress (#66);
 * until then those methods fail loudly rather than silently running solid physics.
 */
#ifndef _RK_SHELL_T_H_
#define _RK_SHELL_T_H_

/* base class */
#include "MeshFreeFSSolidT.h"

namespace Tahoe {

class RKShellT: public MeshFreeFSSolidT
{
public:

	/** constructor */
	RKShellT(const ElementSupportT& support);

	/** \name ParameterInterfaceT */
	/*@{*/
	/** describe the parameters needed by the interface (adds shell_thickness) */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

	/** \name element-group drivers — KL-shell physics (porting in progress, #66) */
	/*@{*/
	/** form the residual (internal force) */
	virtual void RHSDriver(void);

	/** form the tangent stiffness */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);
	/*@}*/

protected:

	/** shell mid-surface thickness */
	double fThickness;
};

} /* namespace Tahoe */

#endif /* _RK_SHELL_T_H_ */
