/* $Id: TotalLagrangianCBSurfaceT.h,v 1.18 2009/06/04 22:45:03 hspark Exp $ */
#ifndef _TOTAL_LAGRANGRIAN_CB_SURFACE_T_H_
#define _TOTAL_LAGRANGRIAN_CB_SURFACE_T_H_

/* base class */
#include "TotalLagrangianT.h"

namespace Tahoe {

class FCC3D_Surf;
class EAMFCC3DMatT_surf;
class EAMFCC3DMatT_edge;
class CB_TersoffT_surf;
class CB_TersoffDimerT_surf;

/** total Lagrangian, finite strain element for working with Cauchy-Born approach
 * for modeling surface effects */
class TotalLagrangianCBSurfaceT: public TotalLagrangianT
{
public:

	/** constructor */
	TotalLagrangianCBSurfaceT(const ElementSupportT& support);

	/** destructor */
	virtual ~TotalLagrangianCBSurfaceT(void);

	/** \name writing output */
	/*@{*/
	/** register self for output */
	virtual void RegisterOutput(void);

	/** send output */
	virtual void WriteOutput(void);
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/* describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;
	
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

	/** accumulate the residual force on the specified node
	 * \param node test node
	 * \param force array into which to assemble to the residual force */
	virtual void AddNodalForce(const FieldT& field, int node, dArrayT& force);

protected:

	/** \name drivers called by ElementBaseT::FormRHS and ElementBaseT::FormLHS */
	/*@{*/
	/** form group contribution to the stiffness matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);

	/** form group contribution to the residual */
	virtual void RHSDriver(void);
	/*@}*/

	/** reduce the coordinates to a surface layer on the given face
	 *\param coords should enter with the coordinates of entire element and
	 *       returns with the coordinates defining the surface layer */
	void SurfaceLayer(LocalArrayT& coords, int face, double thickness) const;

protected:

	/** list of elements on the surface */
	iArrayT fSurfaceElements;

	/** elements neighbors */
	iArray2DT fSurfaceElementNeighbors;

	/** surface model number */
	iArray2DT fSurfaceElementFacesType;

	/** surface normals */
	ArrayT<dArrayT> fNormal;

	/** Tersoff CB model */
	ArrayT<CB_TersoffT_surf*> fTersoffSurfaceCB;

	/** Tersoff CB model for reconstructions */
	ArrayT<CB_TersoffDimerT_surf*> fTersoffDimerSurfaceCB;

    /** EAM Surface CB model */
	ArrayT<EAMFCC3DMatT_surf*> fEAMSurfaceCB;

    /** EAM Edge CB model */
	ArrayT<EAMFCC3DMatT_edge*> fEAMEdgeCB;

	/** surface Cauchy-Born models */
	ArrayT<FCC3D_Surf*> fSurfaceCB;

	/** support for the surface models */
	FSMatSupportT* fSurfaceCBSupport;

	/** deformation gradients at the surface integration points */
	ArrayT<dMatrixT> fF_Surf_List;

	/** indicator for EAM_CB or FCC_CB */
	StringT fIndicator;

	/** \name split integration */
	/*@{*/
	LocalArrayT fSplitInitCoords;
	ShapeFunctionT* fSplitShapes;
	/*@}*/
	
	/** \name surface output */
	/*@{*/
	/** ID obtained during ElementBaseT::RegisterOutput. Each surface type has its own output. */
	iArrayT fSurfaceOutputID;

	/** list of nodes on each surface type (by normal) */
	ArrayT<iArrayT> fSurfaceNodes;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _TOTAL_LAGRANGRIAN_CB_SURFACE_T_H_ */
