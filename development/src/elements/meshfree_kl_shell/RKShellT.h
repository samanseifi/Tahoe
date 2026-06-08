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
#include "dMatrixT.h"
#include "ArrayT.h"

#include <vector>

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

	/** RKShellT treats the nodal DOFs as the (quasi-interpolatory) displacements directly. */
	virtual int InterpolantDOFs(void) const { return 1; }
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

	/** precompute the per-node stencil stiffness K_e (linear elastic; geometry fixed) using the
	 * validated KL-shell kernels: surface PCA chart + RKPM shapes + nodal integration + membrane
	 * and curvature-gradient stabilization, with local-frame plane stress */
	void BuildElementStiffness(void);

	/** stress-driven internal force for node i's stencil: f_out = sum_pt B^T sigma(B*ue) * w.
	 * Base points use the material law (elastic now; plane-stress J2 on the Fig 18 track);
	 * stabilization points are always elastic. Reduces to fKe*ue for linear elasticity. */
	void InternalForce(int i, const dArrayT& ue, dArrayT& fout, bool commit);

	/** finite-deformation internal force for node i's stencil: recompute the current mid-surface
	 * geometry from the current positions, form the objective Green-Lagrange strain E = 1/2(g-G)
	 * (membrane) + (curvature change) from the metric, the current-config B = dE/du, the stress
	 * (elastic or per-point plane-stress J2 on the strain increment), and f = sum B^T S w. */
	void InternalForceFS(int i, const dArrayT& ue, dArrayT& fout, bool commit);

	/** lumped nodal mass m_I = rho * A_I * h (diagonal; for the explicit central-difference
	 * solver, with optional mass scaling via a large fDensity for quasi-static loading) */
	void BuildLumpedMass(void);

private:

	/** \name shell + meshfree parameters */
	/*@{*/
	double fThickness;     /**< shell mid-surface thickness */
	double fYoung;         /**< Young's modulus */
	double fPoisson;       /**< Poisson ratio */
	double fSupportFac;    /**< support size in units of nodal spacing */
	int    fCompleteness;  /**< RKPM completeness (2 = quadratic, 3 = cubic) */
	double fDensity;       /**< mass density (use a scaled value for explicit dynamic relaxation) */
	double fYield;         /**< J2 initial yield stress (0 = elastic, no plasticity) */
	double fHardening;     /**< J2 linear isotropic hardening modulus H: Y(ep) = Yield + H*ep */
	int    fFiniteStrain;  /**< 1 = finite-deformation (Green-Lagrange, current-config geometry) */
	int    fMonitorNode;   /**< 1-based global node to report reaction force at (0 = none); Fig 18 curve */
	/*@}*/

	/** \name stabilization — SCNI/NSNI cell-smoothed assumed-strain residual R = B_direct - B~tilde
	 * (B~tilde = divergence-theorem cell average over node K's quad cell). Consistent (vanishes on
	 * smooth fields -> membrane unpolluted) + PSD (explicit-stable). */
	/*@{*/
	int    fStabMode;      /**< 0 = full strength; 3 = alpha=min(1,h/h_pl) scaled (paper 5.3) */
	double fStabMembrane;  /**< stabilization-residual scale (coefficient on R^T C R) */
	double fStabBending;   /**< (reserved) */
	/*@}*/

	/** \name surface meshfree data */
	/*@{*/
	int fNumNodes;                       /**< number of shell nodes */
	dArray2DT fCoords;                   /**< reference coordinates [fNumNodes] x 3 (local order) */
	dArrayT   fNodalArea;                /**< nodal integration weight per node */
	dArrayT   fLumpedMass;               /**< lumped nodal mass m_I = rho * A_I * h (for explicit) */
	iArrayT   fGlobalToLocal;            /**< global node id -> local shell index (-1 if not a shell node) */
	iArrayT   fGlobalIDs;                /**< local shell index -> global node id */
	RaggedArray2DT<int> fNeighbors;      /**< [node] x [neighbor GLOBAL node ids] */
	RaggedArray2DT<int> fEqnos;          /**< [node] x [neighbor dof equations] */
	ArrayT<dMatrixT> fKe;                /**< per-node stencil stiffness (linear-elastic tangent / implicit LHS) */
	MLSSolverT* fMLS;                    /**< RKPM shape-function solver (local chart) */

	/** \name stress-driven internal force (f_int = sum B^T sigma): per-node integration points.
	 * Each entry stores a Voigt strain-displacement operator B [6 x 3nn] (row-major, flattened),
	 * its weight, and whether it is a stabilization point (always elastic) vs a base material
	 * point. Enables plugging the plane-stress J2 stress update into the base points (Fig 18). */
	/*@{*/
	std::vector<std::vector<double> > fIPB;   /**< [node] -> concatenated [npt * 6 * 3nn] */
	std::vector<std::vector<double> > fIPw;   /**< [node] -> [npt] integration weights */
	std::vector<std::vector<char> >   fIPstab;/**< [node] -> [npt] (1 = elastic stabilization point) */

	/** plane-stress J2 history per node per through-thickness BASE point (fIPstab==0): in-plane
	 * stress [s11,s22,s12], equivalent plastic strain, and previous in-plane strain [e11,e22,g12]
	 * (for the strain increment). Plastic localizes at the surfaces, so each base point is tracked
	 * independently -- never averaged to the mid-surface. */
	std::vector<std::vector<double> > fJ2sig;  /**< [node] -> [nbase * 3] */
	std::vector<std::vector<double> > fJ2ep;   /**< [node] -> [nbase] */
	std::vector<std::vector<double> > fJ2eps;  /**< [node] -> [nbase * 3] */

	/** finite-deformation data: stencil shape derivatives (recompute current geometry from current
	 * positions each step) + reference mid-surface derivatives (for the reference metric/curvature). */
	std::vector<std::vector<double> > fDphi;   /**< [node] -> [nn*5]: Dp0,Dp1,DDp0,DDp1,DDp2 per stencil node */
	std::vector<std::vector<double> > fXref;   /**< [node] -> [15]: ref x,1 x,2 x,11 x,22 x,12 (each 3) */
	dArray2DT fUprev;                          /**< previous-step nodal displacement (rate form: du = u - u_prev) */
	/*@}*/
	/*@}*/

	/** uniform per-area applied load (e.g. gravity); Scordelis-Lo: (0,-90,0) */
	double fLoad[3];

	/** \name output (displacement field on the background cell mesh, for ParaView) */
	/*@{*/
	int fOutputID;                       /**< registered output set id (-1 if none) */
	ArrayT<const iArray2DT*> fOutputConn;/**< background cell connectivity (referenced by the set) */
	iArrayT fOutputNodesUsed;            /**< nodes used by the output set (n_values ordering) */
	/*@}*/

	/** plane-stress (sigma33=0, local normal frame) Voigt tangent */
	double fC[6][6];
};

} /* namespace Tahoe */

#endif /* _RK_SHELL_T_H_ */
