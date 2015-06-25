/* $Id: dArray2DT.cpp,v 1.11 2011/12/01 20:25:15 bcyansfn Exp $ */
/* created: paklein (07/16/1996) */

#include "dArray2DT.h"

#include <iostream>
#include <iomanip>

#include "toolboxConstants.h"
#include "iArrayT.h"
#include "dArrayT.h"

using namespace Tahoe;

/* array behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<dArray2DT*>::fByteCopy = true;
DEFINE_TEMPLATE_STATIC const bool ArrayT<const dArray2DT*>::fByteCopy = true;
DEFINE_TEMPLATE_STATIC const bool ArrayT<dArray2DT>::fByteCopy  = false;
} /* namespace Tahoe */

/* constructor */
dArray2DT::dArray2DT(void) { }
dArray2DT::dArray2DT(int majordim, int minordim):
	nArray2DT<double>(majordim, minordim) { }
dArray2DT::dArray2DT(int majordim, int minordim, const double* p):
	nArray2DT<double>(majordim, minordim, p) { }
dArray2DT::dArray2DT(const dArray2DT& source):
	nArray2DT<double>(source) { }
