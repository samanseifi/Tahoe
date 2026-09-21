/* MFPenaltyDisplacementT.h -- penalty enforcement of a prescribed PHYSICAL displacement
 * component at meshfree nodes (issues #59, #70).
 *
 * For a non-interpolatory meshfree field the physical value at node A is
 * u(x_A) = sum_J Phi_J(x_A) d_J, so a prescribed displacement is the constraint
 * g_A = sum_J Phi_J(x_A) d_J - ubar = 0 over the whole support of x_A. This controller enforces
 * it variationally with a penalty: the generalized force f_J = Phi_J(x_A) lambda acts on EVERY
 * support coefficient, with lambda = -k g - c dg/dt. lambda is the physical point reaction
 * (work-conjugate to u(x_A)), which a coefficient KBC or CollocationKBCT does not provide.
 * The optional dashpot damps the ringing of the penalty spring in explicit dynamics; it is
 * sized as a fraction of critical, c = ratio * 2 sqrt(k m_eff), m_eff = 1 / sum_J Phi_J^2 / m_J,
 * using the lumped masses reported by the element through MeshFreeCollocationSupportT::NodalMass.
 *
 * Deck (under a field):
 *   <penalty_displacement_meshfree dof="1" schedule="1" value="-300.0" penalty="1.0e5"
 *                                  damping_ratio="1.0" element_group="0">
 *     <node_ID_list><String value="3"/></node_ID_list>
 *   </penalty_displacement_meshfree>
 *
 * Output (per output step, to the main output stream), one line per constrained node:
 *   penalty_displacement_meshfree: node N u = ... lambda = ... u_mean = ... lambda_mean = ...
 * where the means are over the steps since the previous output.
 */
#ifndef _MF_PENALTY_DISPLACEMENT_T_H_
#define _MF_PENALTY_DISPLACEMENT_T_H_

/* base class */
#include "FBC_ControllerT.h"

/* direct members */
#include "iArrayT.h"
#include "dArrayT.h"
#include "iArray2DT.h"
#include "RaggedArray2DT.h"
#include "ElementMatrixT.h"
#include "StringT.h"

namespace Tahoe {

/* forward declarations */
class ScheduleT;
class MeshFreeCollocationSupportT;

class MFPenaltyDisplacementT: public FBC_ControllerT
{
public:

	/** constructor */
	MFPenaltyDisplacementT(void);

	/* form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const { return GlobalT::kSymmetric; };

	/** the constraint couples every coefficient in the support of a constrained node */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);
	virtual void Connectivities(AutoArrayT<const iArray2DT*>& connects_1,
		AutoArrayT<const RaggedArray2DT<int>*>& connects_2,
		AutoArrayT<const iArray2DT*>& equivalent_nodes) const;

	/* initial condition/restart functions */
	virtual void InitialCondition(void);

	/** \name apply force and tangent contributions */
	/*@{*/
	virtual void ApplyLHS(GlobalT::SystemTypeT sys_type);
	virtual void ApplyRHS(void);
	/*@}*/

	/* initialize/finalize step */
	virtual void InitStep(void) {};
	virtual void CloseStep(void);

	/* reset to the last known solution */
	virtual void Reset(void) {};

	/** \name writing results */
	/*@{*/
	virtual void RegisterOutput(void) {};
	virtual void WriteOutput(ostream& out) const;
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	virtual void DefineParameters(ParameterListT& list) const;
	virtual void DefineSubs(SubListT& sub_list) const;
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:

	/** locate the meshfree provider, fetch the collocation rows, equation numbers and dashpots */
	void BuildOperator(void);

	/** constraint gap g = sum_J Phi_J d_J - ubar for constraint i at the current field values */
	double Gap(int i, double ubar) const;

private:

	/** \name input parameters */
	/*@{*/
	int    fDOF;             /**< constrained dof (0-based) */
	int    fScheduleNum;     /**< schedule number (0-based) */
	const ScheduleT* fSchedule;
	double fValue;           /**< prescribed physical value, scaled by the schedule */
	double fPenalty;         /**< penalty stiffness k */
	double fDampingRatio;    /**< dashpot as a fraction of critical damping (0 = none) */
	int    fElementGroup;    /**< provider element group (0-based); -1 = auto-discover */
	ArrayT<StringT> fID_List;
	/*@}*/

	/** \name collocation operator (reference configuration) */
	/*@{*/
	const MeshFreeCollocationSupportT* fProvider;
	bool fBuilt;
	iArrayT fNodes;                    /**< constrained global node ids */
	RaggedArray2DT<int>    fSupport;   /**< [constraint] support node ids J */
	RaggedArray2DT<double> fPhi;       /**< [constraint] Phi_J(x_A) */
	RaggedArray2DT<int>    fEqnos;     /**< [constraint] equation numbers of (J, fDOF) */
	dArrayT fDashpot;                  /**< [constraint] dashpot coefficient c */
	/*@}*/

	/** \name state */
	/*@{*/
	dArrayT fLambda;                   /**< [constraint] current multiplier (physical reaction) */
	dArrayT fGapStep;                  /**< [constraint] gap at the last closed step */
	bool    fHaveGapStep;
	mutable dArrayT fUSum, fLambdaSum; /**< running sums since the last output */
	mutable int fCount;
	/*@}*/

	/** \name work space */
	/*@{*/
	dArrayT fForce;
	ElementMatrixT fStiffness;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _MF_PENALTY_DISPLACEMENT_T_H_ */
