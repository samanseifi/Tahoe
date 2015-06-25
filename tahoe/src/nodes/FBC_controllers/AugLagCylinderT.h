/* $Id: AugLagCylinderT.h,v 1.2 2004/12/20 01:23:25 paklein Exp $ */
#ifndef _AUGLAG_CYLINDER_T_H_
#define _AUGLAG_CYLINDER_T_H_

/* base classes */
#include "PenaltyCylinderT.h"
#include "DOFElementT.h"

namespace Tahoe {

/* forward declarations */
class XDOF_ManagerT;

/** rigid cylinder enforced with a penalty formulation */
class AugLagCylinderT: public PenaltyCylinderT, public DOFElementT
{
public:

	/** constructor */
	AugLagCylinderT(void);

	/** initialize data */
	virtual void SetEquationNumbers(void);

	/* append element equations numbers to the list */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);
	virtual void Connectivities(AutoArrayT<const iArray2DT*>& connects_1,
		AutoArrayT<const RaggedArray2DT<int>*>& connects_2,
		AutoArrayT<const iArray2DT*>& equivalent_nodes) const;

	/* restarts */
	virtual void ReadRestart(istream& in);
	virtual void WriteRestart(ostream& out) const;

	/** initialize new time increment */
	virtual void InitStep(void);	

	/** finalize step */
	virtual void CloseStep(void);

	/** returns true if the internal force has been changed since
	 * the last time step. This is when the contact forces are
	 * recomputed for when solving using Uzawa. */
	virtual GlobalT::RelaxCodeT RelaxSystem(void);

	/** tangent
	 * \param sys_type "maximum" tangent type needed by the solver. The GlobalT::SystemTypeT
	 *        enum is ordered by generality. The solver should indicate the most general
	 *        system type that is actually needed. */
	virtual void ApplyLHS(GlobalT::SystemTypeT sys_type);

	/** \name implementation of the DOFElementT interface */
	/*@{*/
	/* returns the array for the DOF tags needed for the current config */
	virtual void SetDOFTags(void);
	virtual iArrayT& DOFTags(int tag_set);

	/* generate nodal connectivities - does nothing here */
	virtual void GenerateElementData(void);

	/* return the contact elements */
	virtual const iArray2DT& DOFConnects(int tag_set) const;

	/* restore the DOF values to the last converged solution */
	virtual void ResetDOF(dArray2DT& DOF, int tag_set) const;

	/* returns 1 if group needs to reconfigure DOF's, else 0 */
	virtual int Reconfigure(void);

	/** restore any state data to the previous converged state */
	virtual void ResetState(void) { };

	/** return the equation group to which the generate degrees of
	 * freedom belong. */
	virtual int Group(void) const;	
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;
	
	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/
	
protected:

	/** accumulate the contact force vector fContactForce */
	virtual void ComputeContactForce(double kforce);

protected:

	/* contact equation sets (shallow copy of contact node equs) */
	iArray2DT fContactEqnos2D;
	iArray2DT fContactTags;
	
	/** \name Augmented multiplier info */
	/*@{*/
	iArrayT fContactDOFtags; /**< contact DOF tags and DOF's */
	iArrayT fFloatingDOF;    /**< 1 if multiplier is attacted to node that has KBC's */
	dArrayT fLastDOF;        /**< multiplier history */
	/*@}*/

	/** \name parameters and data used with Uzawa method */
	/*@{*/
	/** do Uzawa iterations (1st order updates during AugLagWallT::RelaxSystem)
	 * otherwise solve concurrently */
	bool fUzawa;

	int fPrimalIterations;
	double fPenetrationTolerance;

	/** augmented multiplier (only used for Uzawa) */
	dArrayT fDOF;

	/** augmented multiplier from previous iteration. With line searching, the current
	 * value of the multplier should depend on the value from the previous iteration not
	 * on values calculated while performing the line search. Only used for Uzawa. */
	dArrayT fDOFi;
	
	/** iteration number associated with AugLagSphereT::fDOFi */
	int fIterationi;
	
	/** runtime flag */
	bool fRecomputeForce;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _AUGLAG_CYLINDER_T_H_ */
