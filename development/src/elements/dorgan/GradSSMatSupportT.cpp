/* $Id: GradSSMatSupportT.cpp,v 1.11 2004/09/02 18:25:04 rdorgan Exp $ */ 
#include "GradSSMatSupportT.h"
#include "ElementsConfig.h"

#include "GradSmallStrainT.h"

using namespace Tahoe;

/* constructor */
GradSSMatSupportT::GradSSMatSupportT(int ndof_disp, int ndof_pmultiplier, int nip_disp, int nip_pmultiplier):
	SSMatSupportT(ndof_disp, nip_disp),
	fPMultiplier_List(NULL),
	fPMultiplier_last_List(NULL),
	fGradPMultiplier_List(NULL),
	fGradPMultiplier_last_List(NULL),
	fLapPMultiplier_List(NULL),
	fLapPMultiplier_last_List(NULL),
	fGradSmallStrain(NULL)
{

}
 
/* destructor */
GradSSMatSupportT::~GradSSMatSupportT(void)
{

}

/* set source for the isotropic hadening */
void GradSSMatSupportT::SetLinearPMultiplier(const ArrayT<dMatrixT>* pmultiplier_List)
{
	/* keep pointer */
	fPMultiplier_List = pmultiplier_List;
}

/** set source for the isotropic hardening from the end of the previous time step */
void GradSSMatSupportT::SetLinearPMultiplier_last(const ArrayT<dMatrixT>* pmultiplier_last_List)
{
	/* keep pointer */
	fPMultiplier_last_List = pmultiplier_last_List;
}

/* set source for the gradient of isotropic hadening */
void GradSSMatSupportT::SetLinearGradPMultiplier(const ArrayT<dMatrixT>* gradpmultiplier_List)
{
	/* keep pointer */
	fGradPMultiplier_List = gradpmultiplier_List;
}

/** set source for the gradient of isotropic hardening from the end of the previous time step */
void GradSSMatSupportT::SetLinearGradPMultiplier_last(const ArrayT<dMatrixT>* gradpmultiplier_last_List)
{
	/* keep pointer */
	fGradPMultiplier_last_List = gradpmultiplier_last_List;
}

/* set source for the laplacian of isotropic hadening */
void GradSSMatSupportT::SetLinearLapPMultiplier(const ArrayT<dMatrixT>* lappmultiplier_List)
{
	/* keep pointer */
	fLapPMultiplier_List = lappmultiplier_List;
}

/** set source for the laplacian of isotropic hardening from the end of the previous time step */
void GradSSMatSupportT::SetLinearLapPMultiplier_last(const ArrayT<dMatrixT>* lappmultiplier_last_List)
{
	/* keep pointer */
	fLapPMultiplier_last_List = lappmultiplier_last_List;
}

/* set the element group pointer */
void GradSSMatSupportT::SetContinuumElement(const ContinuumElementT* p)
{
	/* inherited */
	SSMatSupportT::SetContinuumElement(p);

	/* cast to grad small strain pointer */
	fGradSmallStrain = TB_DYNAMIC_CAST(const GradSmallStrainT*, p);
}
