/* CollocationKBCT.h — kinematic BC controller that enforces an essential BC on the PHYSICAL
 * displacement (velocity) of a meshfree field by direct nodal collocation (issue #70).
 *
 * For interpolatory FE the deck'd use a plain `kinematic_BC` (sets the nodal coefficient d_I).
 * For non-interpolatory meshfree shapes that is wrong: the physical value is u(x_I)=sum_J Phi_J(x_I) d_J.
 * This controller, given a constrained node set + dof + a uniform prescribed value*schedule, solves
 * each step for the constrained COEFFICIENTS that produce the prescribed physical value (using the
 * CollocationSolverT + the element's MeshFreeCollocationSupportT shapes) and emits them as ordinary
 * kDsp/kVel cards. Reusable by any meshfree element that implements MeshFreeCollocationSupportT.
 *
 * Deck (under a field):
 *   <collocation_KBC field="displacement" element_group="1" dof="1"
 *                    schedule="1" value="-300.0" code="displacement">
 *     <node_ID_list><String value="3"/></node_ID_list>
 *   </collocation_KBC>
 */
#ifndef _COLLOCATION_KBC_T_H_
#define _COLLOCATION_KBC_T_H_

#include "KBC_ControllerT.h"
#include "CollocationSolverT.h"
#include "iArrayT.h"
#include "StringT.h"

namespace Tahoe {

class ScheduleT;
class MeshFreeCollocationSupportT;

class CollocationKBCT: public KBC_ControllerT
{
public:

	CollocationKBCT(const BasicSupportT& support);

	/** build the collocation operator (queries the element provider) */
	virtual void TakeParameterList(const ParameterListT& list);

	/** programmatic setup used by the `collocation="1"` switch on a plain kinematic_BC (the field
	 * routes flagged cards here). The meshfree provider is auto-discovered; no field/element_group
	 * need be named. \param field field name; \param dof 0-based; \param code kDsp/kVel;
	 * \param schedule resolved schedule (may be NULL); \param value prescribed physical value;
	 * \param nodes constrained global node ids. */
	void Configure(const StringT& field, int dof, KBC_CardT::CodeT code,
		const ScheduleT* schedule, double value, const iArrayT& nodes);

	/** recompute the constrained coefficients from the current free coefficients + prescribed value */
	virtual void InitStep(void);

	/** \name ParameterInterfaceT */
	/*@{*/
	virtual void DefineParameters(ParameterListT& list) const;
	virtual void DefineSubs(SubListT& sub_list) const;
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;
	/*@}*/

private:

	void BuildOperator(void);   /**< query provider + factor Phi_c (once, reference config) */

private:

	StringT  fFieldName;        /**< field whose coefficients are read/constrained */
	int      fElementGroup;     /**< 1-based element group providing collocation shapes */
	int      fDOF;              /**< constrained dof (0-based internally) */
	int      fScheduleNum;      /**< schedule (1-based in deck) */
	const ScheduleT* fSchedule;
	double   fValue;            /**< prescribed PHYSICAL value (scaled by the schedule) */
	KBC_CardT::CodeT fCode;     /**< kDsp or kVel */
	ArrayT<StringT> fID_List;   /**< constrained node set ids */

	const MeshFreeCollocationSupportT* fProvider;
	CollocationSolverT fSolver;
	iArrayT  fNodes;            /**< constrained global node ids (solve order) */
	bool     fBuilt;
};

} /* namespace Tahoe */

#endif /* _COLLOCATION_KBC_T_H_ */
