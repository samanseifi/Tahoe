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
#include "MeshFreeCollocationSupportT.h"

/* members */
#include "RaggedArray2DT.h"
#include "dArray2DT.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "dMatrixT.h"
#include "ArrayT.h"
#include "StringT.h"

#include <vector>

namespace Tahoe {

class MLSSolverT;

class RKShellT: public ElementBaseT, public MeshFreeCollocationSupportT
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

	/** \name MeshFreeCollocationSupportT (direct nodal-collocation essential BCs, issue #70) */
	/*@{*/
	/** per-node support ids + Phi_J(x_I) at each requested node's own location (non-interpolatory
	 * RK shapes), so CollocationKBCT can impose the PHYSICAL displacement, not the bare coefficient. */
	virtual bool CollocationData(const iArrayT& nodes,
		RaggedArray2DT<int>& support, RaggedArray2DT<double>& phi) const;
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
	void InternalForceFS(int i, const dArrayT& ue, dArrayT& fout, bool commit, bool include_bend = true);

	/** lumped nodal mass m_I = rho * A_I * h (diagonal; for the explicit central-difference
	 * solver, with optional mass scaling via a large fDensity for quasi-static loading) */
	void BuildLumpedMass(void);

	/** total self-contact force on node i (paper sec 3.13, Eqs 83-86): repulsion from every node K
	 * that is NOT a meshfree neighbor of i and lies within the contact band [r_in, r_out], using the
	 * current node positions (reference + disp). Returns false (and f=0) if contact is disabled. */
	bool ComputeContactForce(int i, const dArray2DT& disp, double f[3]) const;

	/** pinball contact force-density psi(||r||) (Eq 85) with the c1/c2 smoothing of Eq 86 (p=2). */
	double ContactPsi(double r) const;

	/** resolve the rotational-continuity coupling interface node pairs (fCoupleIDa/b -> local pairs) */
	void BuildCoupling(void);

	/** add the sec 3.12 penalty-coupling force (preserve the kink angle) for both members of every
	 * interface pair into the per-node force fout indexed by the pair's stencil position. Returns the
	 * coupling force contribution to node `i` if `i` is one of the paired interface nodes. */
	void AddCouplingForce(const dArray2DT& disp);

	/** self-test for the new multi-body features (env KLSHELL_FEATURETEST): contact force law sign +
	 * monotonicity, and the coupling normal-rotation operator on a 90-degree fold. */
	void RunFeatureSelfTest(void);

	/** curved-surface patch test (env KLSHELL_PATCHTEST): apply the ANALYTIC uniform axial-stretch
	 * field u=(-nu*ea*x,-nu*ea*y, ea*z) to a cylinder mesh and measure the element's membrane strain
	 * error vs the exact (eps_axial=ea, eps_hoop=-nu*ea). Isolates the PCA flat-chart curvature error. */
	void RunCurvedPatchTest(void);

	/** stabilization unit test (triggered by env KLSHELL_SELFTEST): strain energy E = sum_K u_K^T
	 * fKe_K u_K for unit-norm rigid / linear / membrane-hourglass / bending-hourglass modes. Rigid
	 * ->~0; linear must be stab-invariant (consistency); hourglass ~0 without stab, >0 if caught. */
	void RunStabSelfTest(void);

private:

	/** \name shell + meshfree parameters */
	/*@{*/
	double fThickness;     /**< shell mid-surface thickness */
	double fYoung;         /**< Young's modulus */
	double fPoisson;       /**< Poisson ratio */
	double fSupportFac;    /**< support size in units of nodal spacing */
	int    fCompleteness;  /**< RKPM completeness (2 = quadratic, 3 = cubic) */
	int    fKernel;        /**< RK window: 0 = Gaussian (default), 1 = cubic B-spline C² (paper Eq 19) */
	double fDensity;       /**< mass density (use a scaled value for explicit dynamic relaxation) */
	double fYield;         /**< J2 initial yield stress (0 = elastic, no plasticity) */
	double fHardening;     /**< J2 linear isotropic hardening modulus H: Y(ep) = Yield + H*ep */
	double fYieldSat;      /**< J2 saturation yield Ysat (exponential hardening; 0 = linear only) */
	double fSatRate;       /**< J2 saturation rate delta: + (Ysat-Y0)(1-exp(-delta*ep)) */
	int    fFiniteStrain;  /**< 1 = finite-deformation (Green-Lagrange, current-config geometry) */
	int    fThicknessUpdate; /**< 1 = update thickness t=t0/J_area (plastic incompressibility) -> hinge thinning */
	int    fMonitorNode;   /**< 1-based global node to report reaction force at (0 = none); Fig 18 curve */
	int    fMonitorStride; /**< stride between summed reaction nodes (driven generator line = nt) */
	int    fMonitorCount;  /**< number of nodes to sum the reaction over (1 = single load node) */
	double fDamping;       /**< mass-proportional damping alpha (force -alpha*m*v); dynamic relaxation -> quasi-static */

	/** \name self-contact — pinball / volumetric-potential (paper sec 3.13, Eqs 83-86). Off by
	 * default (fContactStiffness=0). A node repels every NON-neighbor node within the contact band
	 * [r_in, r_out] so folding shell surfaces (tube accordion crush) cannot interpenetrate. */
	/*@{*/
	double fContactStiffness; /**< kc (force density scale, N/mm^4); 0 = self-contact OFF */
	double fContactRin;       /**< inner contact radius r_in (full repulsion below this) */
	double fContactRout;      /**< outer contact radius r_out (force tapers to 0 here) */
	/*@}*/

	/** \name rotational-continuity penalty coupling (paper sec 3.12, Eqs 80-82). Off by default
	 * (fPenaltyCoupling=0). Preserves the original angle between adjacent shell patches at a C0 kink
	 * by penalizing the relative normal-rotation jump of paired interface nodes (each side uses its own
	 * one-sided PCA chart). Interface pairs are read from fCoupleIDa / fCoupleIDb node sets. */
	/*@{*/
	double fPenaltyCoupling;  /**< dimensionless C in Cpen = C*E/hpl (Eq 82); 0 = coupling OFF */
	StringT fCoupleIDa;       /**< node_ID of the interface curve on patch (i) */
	StringT fCoupleIDb;       /**< node_ID of the matching interface curve on patch (j) (1:1 ordered) */
	iArrayT fCoupleA;         /**< local indices of the patch-(i) interface nodes */
	iArrayT fCoupleB;         /**< local indices of the patch-(j) interface nodes (paired with fCoupleA) */
	std::vector<double> fCoupleForce; /**< [3*fNumNodes] per-step penalty-coupling force, distributed
	                                   *   to each interface node's own dof in RHSDriver */
	/*@}*/
	/*@}*/

	/** \name stabilization — SCNI/NSNI cell-smoothed assumed-strain residual R = B_direct - B~tilde
	 * (B~tilde = divergence-theorem cell average over node K's quad cell). Consistent (vanishes on
	 * smooth fields -> membrane unpolluted) + PSD (explicit-stable). */
	/*@{*/
	double fStabNatural;   /**< Eq.33 natural Taylor-gradient stabilization scale (0=off, 1=full) -- MEMBRANE only */
	double fStabBending;   /**< curvature-gradient (bending) stabilization: bending analogue of Eq.33 via the in-surface parametric gradient of the curvature operator (BMatrixCurvatureGradient). Catches the bending hourglass Eq.33's one-point through-thickness (membrane-only) quadrature leaves for thin shells. 0=off, 1=consistent. */
	int    fStabSigGrad;   /**< 1 (default) = paper sec 3.10 sigma,xi stress-gradient stabilization ON (Eqs 67-72, reference-config B,xil); 0 = off. */
	double fStabSCNI;      /**< SCNI assumed-strain stabilization coefficient (PSD penalty R^T C R, R=B_direct-B_smoothed): variationally-consistent nodal integration that catches the sawtooth hourglass the Taylor stabilizer aliases over. 0 = off. */
	int    fSmoothedGrad;  /**< 1 = use the SCNI SMOOTHED (cell-averaged, divergence-theorem) shape gradients for the BASE strain integration: variationally consistent -> hourglass-free AND no spurious stiffness (unlike the penalty). 0 = direct nodal gradient. */
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
	std::vector<std::vector<double> > fSelfPhi; /**< [node] x [phi_J(x_I) at node I's own location],
	                                             *   aligned with fNeighbors -> collocation rows (#70) */
	RaggedArray2DT<int> fEqnos;          /**< [node] x [neighbor dof equations] */
	ArrayT<dMatrixT> fKe;                /**< per-node stencil stiffness (linear-elastic tangent / implicit LHS) */
	MLSSolverT* fMLS;                    /**< RKPM shape-function solver (local chart) */
	MLSSolverT* fMLSmem;                 /**< one-order-lower RK solver for the membrane B-bar (anti-locking) */

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
	std::vector<std::vector<double> > fDphiMem; /**< [node] -> [nn*2]: lower-order membrane 1st derivs d1m,d2m for the B-bar (= d1,d2 when B-bar off) */
	std::vector<std::vector<double> > fXref;   /**< [node] -> [15]: ref x,1 x,2 x,11 x,22 x,12 (each 3) */
	dArray2DT fUprev;                          /**< previous-step nodal displacement (rate form: du = u - u_prev) */
	std::vector<double> fThicknessCur;         /**< [node] current shell thickness (accumulates D33 from the
	                                            *   sigma33=0 update, Algorithm 3); init = fThickness */
	std::vector<double> fFrameR;               /**< [node] x 9: co-rotational rotation R (Flanagan-Taylor,
	                                            *   Algorithm 2); init = reference tangent frame [E1 E2 N0] */
	std::vector<double> fFrameV;               /**< [node] x 9: left-stretch V (Flanagan-Taylor); init = I */
	int    fCorotational;                      /**< 1 = Algorithm 2 co-rotational stress frame (objective at
	                                            *   large rotation); 0 = legacy OrthoTangents frame */

	/** accumulated mid-surface Cauchy STRESS GRADIENT history (paper sec 3.10, Eq 68): per node,
	 * [sigma,xi1: s11,s22,s12 ; sigma,xi2: s11,s22,s12] in the local shell frame. Advanced each step by
	 * sigma,xil += C~^P_mid . (B,xil . du) and co-rotated with the Flanagan-Taylor frame (Eq 72). Drives
	 * the paper-faithful membrane stabilization force f_stab = sum_l (B,xil)^T sigma,xil * V_K M_Kxil. */
	std::vector<std::vector<double> > fSigGrad; /**< [node] -> [6] */
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
