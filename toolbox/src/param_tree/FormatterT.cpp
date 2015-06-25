/* $Id: FormatterT.cpp,v 1.3 2011/12/01 20:25:17 bcyansfn Exp $ */
#include "FormatterT.h"

#include <cstring>

/* 1 less than length of fTabs */
const int kMaxTab = 10;

using namespace Tahoe;

/* constructor */
FormatterT::FormatterT(void):
	fTabCount(0)
{
	memset((void *) fTabs, '\0', (kMaxTab + 1)*sizeof(char));
}

/**********************************************************************
 * Private
 **********************************************************************/

const char* FormatterT::TabOut(void) const
{
	/* non-const this */
	FormatterT* non_const_this = (FormatterT*) this;

	/* grow length of tabs */
	if (fTabCount < (kMaxTab-1)) (non_const_this->fTabs)[(non_const_this->fTabCount)++] = '\t';
	return fTabs;
}

const char* FormatterT::TabIn(void) const
{
	/* non-const this */
	FormatterT* non_const_this = (FormatterT*) this;

	/* shorten length of tabs */
	if (fTabCount > 0) (non_const_this->fTabs)[--(non_const_this->fTabCount)] = '\0';
	return fTabs;
}

