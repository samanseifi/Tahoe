//
// $Id: FSPZMatSupportT.cpp,v 1.2 2008/12/12 18:58:15 amota Exp $
//
// $Log: FSPZMatSupportT.cpp,v $
// Revision 1.2  2008/12/12 18:58:15  amota
// Numerous changes.
//
// Revision 1.1  2008/06/16 18:21:41  lxmota
// Piezoelectric material support. Initial sources.
//
//

#include "FSPZMatSupportT.h"
#include "FSPiezoElectricSolidT.h"

namespace Tahoe{

  //
  //
  //
  void
  FSPZMatSupportT::SetContinuumElement(const ContinuumElementT* p)
  {

    //
    // inherited
    //
    FSMatSupportT::SetContinuumElement(p);

    //
    // cast to finite deformation piezoelectric pointer
    //
    fFSPiezoElectricSolid =
      dynamic_cast<const FSPiezoElectricSolidT*>(p);

  }

} //namespace Tahoe
