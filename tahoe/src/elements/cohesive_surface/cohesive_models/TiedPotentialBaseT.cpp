/* $Id: TiedPotentialBaseT.cpp,v 1.3 2011/12/01 21:11:36 bcyansfn Exp $  */
/* created: cjkimme (10/23/2001) */
#include "TiedPotentialBaseT.h"

#include <iostream>
#include <cmath>

#include "ExceptionT.h"
#include "iArrayT.h"

using namespace Tahoe;

/* constants for state variable flags */
const double TiedPotentialBaseT::kTiedNode = -100.;
const double TiedPotentialBaseT::kReleaseNextStep = -10;
const double TiedPotentialBaseT::kFirstFreeStep = -1.;
const double TiedPotentialBaseT::kFreeNode = 1.;
const double TiedPotentialBaseT::kTieNextStep = 10.;

/* constructor */
TiedPotentialBaseT::TiedPotentialBaseT(void):
	iBulkGroups(0)
{
	// Do nothing
}

/* destructor */
TiedPotentialBaseT::~TiedPotentialBaseT(void) 
{
	// Nope 
}

const iArrayT& TiedPotentialBaseT::BulkGroups(void) const
{
  	return iBulkGroups;
}


