/* $Id: PenaltyWallT.h,v 1.13 2005/11/20 23:08:29 regueiro Exp $ */
/* created: paklein (02/25/1997) */
#ifndef _PENATLY_WALL_T_H_
#define _PENATLY_WALL_T_H_

/* base class */
#include "PenaltyRegionT.h"

/* direct members */
#include "ElementMatrixT.h"

namespace Tahoe {

/** flat rigid barrier enforced with a penalized constraint */
class PenaltyWallT: public PenaltyRegionT
{
public:

	/** constructor */
	PenaltyWallT(void);

	/** tangent
	 * \param sys_type "maximum" tangent type needed by the solver. The GlobalT::SystemTypeT
	 *        enum is ordered by generality. The solver should indicate the most general
	 *        system type that is actually needed. */
	virtual void ApplyLHS(GlobalT::SystemTypeT sys_type);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:

	/* accumulate the contact force vector fContactForce */
	virtual void ComputeContactForce(double kforce);

protected:

	/* wall parameters */
	dArrayT		fnormal;
	dArrayT		ftangent;
	dArrayT		normal_displ;
	dArrayT		tangent_displ;
	
	/* wall normal and tangents */
	dArrayT		fntforce;
	dArrayT		fnforce;
	dArrayT		ftforce;
	dArrayT		fxyforce;
	dMatrixT	fQ;

	/* relative displacements */
	dArray2DT	fp_i; //relative displacement
	dArray2DT	fv_i; //relative velocity
	dArray2DT	fprev_i; // previous contact nodal coordinates
	dArray2DT	fcurr_i; // current contact nodal coordinates
	dArray2DT	fdispl_i; // current contact nodal displacements
	
	/* workspace */
	ElementMatrixT fLHS;  //tangent matrix
	dArrayT        fd_sh; //shallow
	iArrayT        fi_sh; //shallow
};

} // namespace Tahoe 
#endif /* _PENATLY_WALL_T_H_ */
