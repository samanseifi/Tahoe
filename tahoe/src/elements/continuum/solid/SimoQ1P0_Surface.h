
#if !defined(_SimoQ1P0_Surface_)
#define _SimoQ1P0_Surface_

#include <cassert>

/* base class */
#include "SimoQ1P0.h"

namespace Tahoe {

  // interface for finite deformation dielectric elastomers
  // based on 2008 JMPS paper of Suo et al.

  class SimoQ1P0_Surface: public SimoQ1P0 {

  public:

    // constructor
    SimoQ1P0_Surface(const ElementSupportT& support);

    // destructor
    virtual ~SimoQ1P0_Surface();

    // specify parameters needed by the interface
    virtual void DefineParameters(ParameterListT& list) const;

    // accept parameter list
    virtual void TakeParameterList(const ParameterListT& list);

    // element stiffness matrix
    virtual void FormStiffness(double constK);

    // internal force
    virtual void FormKd(double constK);

	/* TLCBSurfaceT stuff */
	virtual void DefineSubs(SubListT& sub_list) const;

  protected:

    // driver for calculating output values
    virtual void ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
        const iArrayT& e_codes, dArray2DT& e_values);

	// Return nodes for canonical element based on normal type
	iArrayT CanonicalNodes(const int node_index0, const int node_index1);

  private:

  protected:

	/** current coords with local ordering */
	LocalArrayT fLocCurrCoords;

	/** \name work space - from UpdatedLagrangianT.h */
	/*@{*/
	dMatrixT fGradNa;       /**< shape function gradients matrix: [nsd] x [nen] */
	/*@}*/

	/* TLCBSurfaceT Stuff */
	/** list of elements on the surface */
	iArrayT fSurfaceElements;

	/** elements neighbors */
	iArray2DT fSurfaceElementNeighbors;

	/** surface model number */
	iArray2DT fSurfaceElementFacesType;

	/** surface normals */
	ArrayT<dArrayT> fNormal;

	/** support for the surface models */
	FSMatSupportT* fSurfaceCBSupport;

	/** deformation gradients at the surface integration points */
	ArrayT<dMatrixT> fF_Surf_List;

	/** \name surface output */
	/*@{*/
	/** ID obtained during ElementBaseT::RegisterOutput. Each surface type has its own output. */
	iArrayT fSurfaceOutputID;

	/** list of nodes on each surface type (by normal) */
	ArrayT<iArrayT> fSurfaceNodes;
	/*@}*/

	dMatrixT fGrad_U;

  private:

    // Stiffness storage
    dMatrixT fAmm_mat2;	// mechanical material part of Hessian matrix
    dMatrixT fAmm_geo2;	// mechanical geometric part of Hessian matrix
    dMatrixT fB2, fD2, fLHS2;
    double fSurfTension, fNewSurfTension;
    double fT_0;

	dMatrixT tempstiff;
  };

} // namespace Tahoe

#endif // _SimoQ1P0_Surface_
