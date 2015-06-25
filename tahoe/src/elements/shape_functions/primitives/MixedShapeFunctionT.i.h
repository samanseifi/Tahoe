//
// $Id: MixedShapeFunctionT.i.h,v 1.1 2008/12/12 00:40:26 lxmota Exp $
//
// $Log: MixedShapeFunctionT.i.h,v $
// Revision 1.1  2008/12/12 00:40:26  lxmota
// Additional interpolation scheme for multi-field formulations. Initial sources.
//
//
namespace Tahoe {

  // compute jacobians of the nodal values at integration points
  inline void MixedShapeFunctionT::JacobianDets(dArrayT& dets)
  {
    ShapeFunctionT::JacobianDets(dets);
  }

} // namespace Tahoe
