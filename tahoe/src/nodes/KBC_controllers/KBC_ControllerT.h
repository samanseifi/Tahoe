/* $Id: KBC_ControllerT.h,v 1.27 2009/05/26 13:03:26 tdnguye Exp $ */
/* created: paklein (09/05/2000) */
#ifndef _KBC_CONTROLLER_T_H_
#define _KBC_CONTROLLER_T_H_

/* base class */
#include "ParameterInterfaceT.h"

/* direct members */
#include "ArrayT.h"
#include "KBC_CardT.h"
#include "GlobalT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class ofstreamT;
class nIntegratorT;
class iArrayT;
class StringT;
class dArrayT;
class iArray2DT;
template <class TYPE> class AutoArrayT;
class BasicSupportT;

/** base class for all kinematic BC controllers. Classes that
 * implement more than simple boundary conditions */
class KBC_ControllerT: public ParameterInterfaceT
{
public:

	/** controller codes - derived classes */
	enum CodeT {   kNone =-1,
	            kK_Field = 0,
      kBimaterialK_Field = 1,
         kMappedPeriodic = 2,
              kTiedNodes = 3,
          kPeriodicNodes = 5,
             kPrescribed = 6,
    kScaledVelocityNodes = 7,
          kSetOfNodesKBC = 8,
                kTorsion = 9,
		       kConveyor = 10,
		       kConveyorSym = 11,
//	            kK_Field_3D = 12,
//				kAngledBC = 13
                };

	/** converts strings to KBC_ControllerT::CodeT */
	static CodeT Code(const char* name);

	/** constructor */
	KBC_ControllerT(const BasicSupportT& support);

	/** destructor */
	virtual ~KBC_ControllerT(void);

	/** boundary condition cards generated by the controller */
	const ArrayT<KBC_CardT>& KBC_Cards(void) const;

	/** non-const access boundary condition cards generated by the controller */
	ArrayT<KBC_CardT>& KBC_Cards(void);
	
	/** inform controller of external nodes */
	virtual void SetExternalNodes(const ArrayT<int>& ex_nodes) const;

	/** set to initial conditions */
	virtual void InitialCondition(void) {};

	/** \name restart functions */
	/*@{*/
	virtual void ReadRestart(ifstreamT& in);
	virtual void WriteRestart(ofstreamT& out) const;
	/*@}*/

	/** \name solution steps
	 * Methods signalling different stages of the solution process for
	 * a single time step. */
	/*@{*/
	/** initialize the current step */
	virtual void InitStep(void) { };

	/** computing residual force */
	virtual void FormRHS(void) { };

	/** computing tangent */
	virtual void FormLHS(GlobalT::SystemTypeT) { };
	
	/** apply the update to the solution. Does nothing by default. */
	virtual void Update(const dArrayT& update);

	/** signal that the solution has been found */
	virtual void CloseStep(void) { };

	/** solution for the current step failed. Restore system to its
	 * state at the beginning of the current time step. */
	virtual void Reset(void) { };
	/*@}*/

	/** returns true if the internal force has been changed since
	 * the last time step */
	virtual GlobalT::RelaxCodeT RelaxSystem(void);

	/** output current configuration */
	virtual void WriteOutput(ostream& out) const;

	/** append connectivities generated by the controller. By default, 
	 * does not append any connectivities to the list */
	virtual void Connectivities(AutoArrayT<const iArray2DT*>& connects, 
						AutoArrayT<const iArray2DT*>& equivalent_nodes) const;

	/** append equation number sets generated by the controller */
	virtual void Equations(AutoArrayT<const iArray2DT*>& equations) const;

	virtual void SetEquations(void);
	
	/** inform the Field if acting only as an IC_controller -- do not prescribe */
	virtual bool IsICController(void) { return false; }
	
protected:

	/** get nodes from the ModelManagerT
	 * \param in input stream listing the node ids
	 * \param id_list returns with the set id's of the nodes
	 * \param nodes returns with the nodes in the set id's */
	void GetNodes(const ArrayT<StringT>& id_list, iArrayT& nodes) const;

private:

	/** \name disallowed */
	/*@{*/
	KBC_ControllerT(KBC_ControllerT &);
	KBC_ControllerT& operator=(KBC_ControllerT &);
	/*@}*/
	
protected:

	/** host code support */
	const BasicSupportT& fSupport;

	/* boundary conditions cards - return value */
	ArrayT<KBC_CardT> fKBC_Cards;  	
};

/* boundary condition cards generated by the controller */
inline const ArrayT<KBC_CardT>& KBC_ControllerT::KBC_Cards(void) const
{
	return fKBC_Cards;
}

inline ArrayT<KBC_CardT>& KBC_ControllerT::KBC_Cards(void) { return fKBC_Cards; }

inline void KBC_ControllerT::Update(const dArrayT& update)
{
#pragma unused(update)
}

inline void KBC_ControllerT::Connectivities(AutoArrayT<const iArray2DT*>& connects, 
									AutoArrayT<const iArray2DT*>& equivalent_nodes) const
{
#pragma unused(connects)
#pragma unused(equivalent_nodes)
}

inline void KBC_ControllerT::Equations(AutoArrayT<const iArray2DT*>& equations) const
{
#pragma unused(equations)
}

inline void KBC_ControllerT::SetEquations(void)
{
}

inline void KBC_ControllerT::SetExternalNodes(const ArrayT<int>& ex_nodes) const
{
#pragma unused(ex_nodes)
}

} // namespace Tahoe 
#endif /* _KBC_CONTROLLER_T_H_ */
