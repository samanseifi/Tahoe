/* RKShellT.h — meshfree RKPM Kirchhoff-Love shell element (epic #59)
 *
 * Wang & Bazilevs (2025) general-purpose meshfree Kirchhoff-Love shell, as a Tahoe element
 * driven by XML input. Element type "meshfree_kl_shell".
 *
 * ARCHITECTURE: derives from ElementBaseT (NOT the 3D-solid MeshFreeFSSolidT). A KL shell is
 * a 2D manifold in 3D; the 3D-solid meshfree setup builds volume shape functions that reject a
 * surface node cloud ("geometry code does not match coordinates"). So this element builds its
 * own SURFACE meshfree support: per-node neighbor lists, per-node PCA tangent-plane chart, and
 * quadratic/cubic RKPM shape functions (MLSSolverT) in that local chart. The validated physics
 * lives in KLShellKernels.h (BuildGeom, BMatrix, ToVoigtLocal -- local-frame plane stress --,
 * the curvature-gradient stabilization). Nodal integration + 3-pt through-thickness Gauss.
 *
 * Assembly goes through the framework (Equations/ConnectsU register the ragged neighbor
 * connectivity; RHS/LHS assemble per-node stencils), so the global sparse matrix + solver
 * (profile/SPOOLES/MUMPS) are chosen in the input deck.
 */
#ifndef _RK_SHELL_T_H_
#define _RK_SHELL_T_H_

/* base class */
#include "ElementBaseT.h"

/* members */
#include "RaggedArray2DT.h"
#include "dArray2DT.h"
#include "dArrayT.h"
#include "iArrayT.h"

namespace Tahoe {

class MLSSolverT;

class RKShellT: public ElementBaseT
{
public:

	/** constructor */
	RKShellT(const ElementSupportT& support);

	/** destructor */
	virtual ~RKShellT(void);

	/** \name ElementBaseT required interface */
	/*@{*/
	virtual GlobalT::SystemTypeT TangentType(void) const;
	virtual void AddNodalForce(const FieldT& field, int node, dArrayT& force);
	virtual double InternalEnergy(void);
	virtual void SendOutput(int kincode);
	virtual void RegisterOutput(void);
	virtual void WriteOutput(void);
	/*@}*/

	/** \name connectivity / equations (register the ragged neighbor stencils) */
	/*@{*/
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);
	virtual void ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
		AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const;
	virtual void ConnectsX(AutoArrayT<const iArray2DT*>& connects) const;
	/*@}*/

	/** \name ParameterInterfaceT */
	/*@{*/
	virtual void DefineParameters(ParameterListT& list) const;
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** \name element-group drivers */
	/*@{*/
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);
	virtual void RHSDriver(void);
	/*@}*/

	/** override: this element defines its own meshfree node set, not FE blocks */
	virtual void DefineElements(const ArrayT<StringT>& block_ID, const ArrayT<int>& mat_index);

private:

	/** build per-node neighbor lists (3D distance within the support) */
	void BuildNeighbors(void);

private:

	/** \name shell + meshfree parameters */
	/*@{*/
	double fThickness;     /**< shell mid-surface thickness */
	double fYoung;         /**< Young's modulus */
	double fPoisson;       /**< Poisson ratio */
	double fSupportFac;    /**< support size in units of nodal spacing */
	int    fCompleteness;  /**< RKPM completeness (2 = quadratic, 3 = cubic) */
	/*@}*/

	/** \name surface meshfree data */
	/*@{*/
	int fNumNodes;                       /**< number of shell nodes */
	dArray2DT fCoords;                   /**< reference coordinates [fNumNodes] x 3 */
	dArrayT   fNodalArea;                /**< nodal integration weight per node */
	RaggedArray2DT<int> fNeighbors;      /**< [node] x [neighbor node ids] (0-based) */
	RaggedArray2DT<int> fEqnos;          /**< [node] x [neighbor dof equations] */
	MLSSolverT* fMLS;                    /**< RKPM shape-function solver (local chart) */
	/*@}*/

	/** plane-stress (sigma33=0, local normal frame) Voigt tangent */
	double fC[6][6];
};

} /* namespace Tahoe */

#endif /* _RK_SHELL_T_H_ */
