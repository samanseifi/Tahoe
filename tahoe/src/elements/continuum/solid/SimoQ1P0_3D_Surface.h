/* $Id: SimoQ1P0_3D_Surface.h,v 1.1 2026/03/02 00:00:00 samanseifi Exp $ */
#if !defined(_SimoQ1P0_3D_Surface_)
#define _SimoQ1P0_3D_Surface_

#include <cassert>

/* base class */
#include "SimoQ1P0.h"
#include "ParameterContainerT.h"

namespace Tahoe {

  /**
   * Updated Lagrangian Q1P0 element with surface tension for 3D hexahedral elements.
   *
   * Each free hex face (4-node bilinear quad) contributes a surface energy
   * W_s = gamma * int H dxi1 dxi2  (H = sqrt(EG - F^2))
   * integrated with 2x2 Gauss quadrature using the covariant metric tensor approach
   * (follows Henann & Bertoldi 2014 Abaqus UEL: uel_surften_3D_hex.for).
   *
   * Surface tension is specified per side set via <surface_tension> sub-lists in XML.
   */
  class SimoQ1P0_3D_Surface : public SimoQ1P0 {

  public:

    /* constructor */
    SimoQ1P0_3D_Surface(const ElementSupportT& support);

    /* destructor */
    virtual ~SimoQ1P0_3D_Surface();

    /* specify parameters needed by the interface */
    virtual void DefineParameters(ParameterListT& list) const;

    /* accept parameter list */
    virtual void TakeParameterList(const ParameterListT& list);

    /* element stiffness matrix */
    virtual void FormStiffness(double constK);

    /* internal force */
    virtual void FormKd(double constK);

    /* construct sub-list parameter interfaces */
    virtual ParameterInterfaceT* NewSub(const StringT& name) const;

    /* define subordinate parameter lists */
    virtual void DefineSubs(SubListT& sub_list) const;

  protected:

    /* driver for calculating output values */
    virtual void ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
        const iArrayT& e_codes, dArray2DT& e_values);

    /**
     * Fill dof_indices[12] with element-level DOF indices for the 4 face nodes.
     * With INTERLEAVED ordering (node-major): DOF for node n, dir d = 3*n + d.
     * face_nodes_index[a] = element-local node index for face node a (a=0..3).
     */
    void CanonicalNodes3D(const iArrayT& face_nodes_index, iArrayT& dof_indices) const;

  protected:

    /** current coords with local ordering */
    LocalArrayT fLocCurrCoords;

    /** shape function gradient workspace */
    dMatrixT fGradNa;

    /** list of elements on the surface */
    iArrayT fSurfaceElements;

    /** element neighbors (== -1 means free face) */
    iArray2DT fSurfaceElementNeighbors;

    /** surface model number per face */
    iArray2DT fSurfaceElementFacesType;

    /** surface normals (6 for hex) */
    ArrayT<dArrayT> fNormal;

    /** support for surface models */
    FSMatSupportT* fSurfaceCBSupport;

    /** deformation gradients at surface integration points */
    ArrayT<dMatrixT> fF_Surf_List;

    /** output ID per surface type */
    iArrayT fSurfaceOutputID;

    /** nodes on each surface type */
    ArrayT<iArrayT> fSurfaceNodes;

    dMatrixT fGrad_U;

  private:

    /** surface tension gamma per (surface element, face) — 0 means no contribution */
    dArray2DT fSurfaceGamma;

    /** ramp-up time t_0 per (surface element, face) */
    dArray2DT fSurfaceT0;

    /** mechanical Hessian storage: nme x nme = 24x24 for 8-node hex */
    dMatrixT fAmm_mat2;
  };

} // namespace Tahoe

#endif // _SimoQ1P0_3D_Surface_
