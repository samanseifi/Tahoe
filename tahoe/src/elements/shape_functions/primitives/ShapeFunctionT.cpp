/* $Id: ShapeFunctionT.cpp,v 1.22 2009/04/02 00:37:39 lxmota Exp $ */
/* created: paklein (06/26/1996) */
#include "ShapeFunctionT.h"
#include "ParentDomainT.h"
#include "LocalArrayT.h"
#include "dSymMatrixT.h"

using namespace Tahoe;

/* constructor */
ShapeFunctionT::ShapeFunctionT(GeometryT::CodeT geometry_code, int numIP,
    const LocalArrayT& coords, bool is_closed_set) :
      DomainIntegrationT(geometry_code, numIP, coords.NumberOfNodes(),
          is_closed_set), fCoords(coords), fGrad_x_temp(NULL), fStore(false),
      fCurrElementNumber(NULL)
{
	/* consistency */
	if (GeometryT::GeometryToNumSD(geometry_code) != fCoords.MinorDim())
		ExceptionT::GeneralFail("ShapeFunctionT::ShapeFunctionT",
			"geometry code %d does not match coordinates in %d dimensions",
			geometry_code, fCoords.MinorDim());

	/* configure workspace */
	Construct();
}

ShapeFunctionT::ShapeFunctionT(GeometryT::CodeT geometry_code, int numIP,
    const LocalArrayT& coords, int dummy_flag, bool is_closed_set) :
      DomainIntegrationT(geometry_code, numIP, coords.NumberOfNodes(),
          is_closed_set), fCoords(coords), fGrad_x_temp(NULL), fStore(false),
      fCurrElementNumber(NULL)
{
	/* consistency */
	if (GeometryT::GeometryToNumSD(geometry_code) != fCoords.MinorDim())
		ExceptionT::GeneralFail("ShapeFunctionT::ShapeFunctionT",
			"geometry code %d does not match coordinates in %d dimensions",
			geometry_code, fCoords.MinorDim());

	/* configure workspace */
	Construct_DN_DDN();
}

ShapeFunctionT::ShapeFunctionT(const ShapeFunctionT& link,
    const LocalArrayT& coords, bool is_closed_set) :
	DomainIntegrationT(link, is_closed_set),
	fCoords(coords),
	fGrad_x_temp(NULL),
	fStore(false),
	fCurrElementNumber(NULL)
{
	/* configure workspace */
	Construct();
}

/* compute local shape functions and derivatives */
void ShapeFunctionT::SetDerivatives(void)
{
	/* fetch from storage */
	if (fStore)
	{
		/* check */
		if (!fCurrElementNumber)
			ExceptionT::GeneralFail("ShapeFunctionT::SetDerivatives",
				"current element not set");

		/* get Jocobian information */
		fDet_store.RowCopy(*fCurrElementNumber, fDet);

		/* get shape function derivatives */
		for (int i = 0; i < fDNaX_store.Length(); i++)
			fDNaX_store[i].RowCopy(*fCurrElementNumber, fDNaX[i]);
	}
	else /* compute values */
		fDomain->ComputeDNa(fCoords, fDNaX, fDet);
}

void ShapeFunctionT::SetDerivatives_DN_DDN(void)
{
	/* fetch from storage */
	if (fStore)
	{
		/* check */
		if (!fCurrElementNumber)
			ExceptionT::GeneralFail("ShapeFunctionT::SetDerivatives_DN_DDN",
				"current element not set");

		/* get Jocobian information */
		fDet_store.RowCopy(*fCurrElementNumber, fDet);

		/* get shape function derivatives */
		for (int i = 0; i < fDNaX_store.Length(); i++)
			fDNaX_store[i].RowCopy(*fCurrElementNumber, fDNaX[i]);
	}
	else /* compute values */
		fDomain->ComputeDNa_DDNa(fCoords, fDNaX, fDDNaX, fDet);
}

/* field gradients at specific parent domain coordinates. */
void ShapeFunctionT::GradU(const LocalArrayT& nodal, dMatrixT& grad_U, const dArrayT& coord,
	dArrayT& Na, dArray2DT& DNa) const
{
	/* dimensions */
	int nnd = fCoords.NumberOfNodes();
	int nsd = fCoords.MinorDim();
	Na.Dimension(nnd);
	DNa.Dimension(nsd, nnd);

	/* compute the shape function derivatives */
	EvaluateShapeFunctions(coord, Na, DNa);

	/* compute the gradient */
	fDomain->Jacobian(nodal, DNa, grad_U);
}

/************************ for the current integration point *********************/
void ShapeFunctionT::InterpolateU(const LocalArrayT& nodal,
	ArrayT<double>& u) const
{
#if __option(extended_errorcheck)
	if (nodal.MinorDim() != u.Length() ||
	    nodal.NumberOfNodes() != pNaU->MinorDim())
	    ExceptionT::SizeMismatch("ShapeFunctionT::InterpolateU");
#endif

	int num_u = nodal.MinorDim();
	for (int i = 0; i < num_u; i++)
		u[i] = pNaU->DotRow(fCurrIP, nodal(i));
}

void ShapeFunctionT::InterpolateU(const LocalArrayT& nodal,
	ArrayT<double>& u, int ip) const
{
#if __option(extended_errorcheck)
	if (nodal.MinorDim() != u.Length() ||
	    nodal.NumberOfNodes() != pNaU->MinorDim())
	    ExceptionT::SizeMismatch("ShapeFunctionT::InterpolateU");
#endif

	int num_u = nodal.MinorDim();
	for (int i = 0; i < num_u; i++)
		u[i] = pNaU->DotRow(ip, nodal(i));
}

// interpolate field values to arbitrary point in parametric coordinates
// \param nodal array of nodal values: [nnd] x [nu]
// \param u interpolation of the nodal values
// \param ip integration point number
void ShapeFunctionT::InterpolateU(const LocalArrayT& nodal, ArrayT<double>& u,
    const dArrayT& coordinates) const
{
  int nnd = nodal.NumberOfNodes();
  int nnv = nodal.MinorDim();

#if __option(extended_errorcheck)
  if (nnv != u.Length() || nnd != pNaU->MinorDim()) {
    ExceptionT::SizeMismatch("ShapeFunctionT::InterpolateU");
  }
#endif

  dArrayT Na;
  Na.Dimension(nnd);

  EvaluateShapeFunctions(coordinates, Na);

  for (int i = 0; i < nnv; ++i) {
    u[i] = 0.0;
    for (int j = 0; j < nnd; ++j) {
      u[i] += nodal(j, i) * Na[j];
    }
  }
}

/* second derivative of shape functions */
void ShapeFunctionT::Grad_GradNa(const dArray2DT& DDNa, dMatrixT& grad_grad_Na) const
{
#if __option(extended_errorcheck)
	if (DDNa.MajorDim() != grad_grad_Na.Rows() ||
	    DDNa.MinorDim() != grad_grad_Na.Cols())
	    ExceptionT::SizeMismatch("ShapeFunctionT::Grad_GradNa");
#endif

	int numsd_x_2 = DDNa.MajorDim();
	int numnodes  = DDNa.MinorDim();
	double* p     = grad_grad_Na.Pointer();

	if (numsd_x_2 == 6)
	{
		const double* pNaxx = DDNa(0);
		const double* pNayy = DDNa(1);
		const double* pNazz = DDNa(2);
		const double* pNayz = DDNa(3);
		const double* pNaxz = DDNa(4);
		const double* pNaxy = DDNa(5);

		for (int i = 0; i < numnodes; i++)
		{
			*p++ = *pNaxx++;
			*p++ = *pNayy++;
			*p++ = *pNazz++;
			*p++ = *pNayz++;
			*p++ = *pNaxz++;
			*p++ = *pNaxy++;
		}
	}
}



/* shape function gradients matrix (Hughes,4.90) */
void ShapeFunctionT::GradNa(const dArray2DT& DNa, dMatrixT& grad_Na) const
{
#if __option(extended_errorcheck)
	if (DNa.MajorDim() != grad_Na.Rows() ||
	    DNa.MinorDim() != grad_Na.Cols())
	    ExceptionT::SizeMismatch("ShapeFunctionT::GradNa");
#endif

	int numsd    = DNa.MajorDim();
	int numnodes = DNa.MinorDim();
	double* p    = grad_Na.Pointer();

	if (numsd == 2)
	{
		const double* pNax = DNa(0);
		const double* pNay = DNa(1);

		for (int i = 0; i < numnodes; i++)
		{
			*p++ = *pNax++;
			*p++ = *pNay++;
		}
	}
	else if (numsd == 3)
	{
		const double* pNax = DNa(0);
		const double* pNay = DNa(1);
		const double* pNaz = DNa(2);

		for (int i = 0; i < numnodes; i++)
		{
			*p++ = *pNax++;
			*p++ = *pNay++;
			*p++ = *pNaz++;
		}
	}
	else /* 1D case defaults to this */
	{
		for (int i = 0; i < numsd; i++)
			for (int a = 0; a < numnodes; a++)
				grad_Na(i,a) = DNa(i,a);
	}
}


/* laplacian of field at current ip */
/*
void ShapeFunctionT::LaplaceU(const LocalArrayT& field, dArrayT& laplacian) const
{

}
*/

/* laplacian of strain at current ip */
/*
void ShapeFunctionT::LaplaceStrain(const dSymMatrixT& strain, dSymMatrixT& laplacian) const
{

}
*/


/********************************************************************************/

/* print the shape function values to the output stream */
void ShapeFunctionT::Print(ostream& out) const
{
	/* inherited */
	DomainIntegrationT::Print(out);

	out << "\n Domain shape function derivatives:\n";
	for (int i = 0; i < fDNaX.Length(); i++)
		fDNaX[i].WriteNumbered(out);

	out << "\n Field shape function derivatives:\n";
	if (fDNaX.Pointer() == pDNaU->Pointer())
	    out << " isoparametric \n";
	else
	    for (int i = 0; i < pDNaU->Length(); i++)
			(*pDNaU)[i].WriteNumbered(out);
}

void ShapeFunctionT::InitStore(int num_elements, const int* curr_element)
{
	const char caller[] = "ShapeFunctionT::InitStore";
	if (num_elements > 0 && !curr_element)
		ExceptionT::GeneralFail(caller);
	fCurrElementNumber = curr_element;

	/* allocate work space */
	fDet_store.Dimension(num_elements, fDet.Length());
	fDNaX_store.Dimension(fDNaX.Length());
	for (int i = 0; i < fDNaX_store.Length(); i++)
		fDNaX_store[i].Dimension(num_elements, fDNaX[i].Length());

	/* set flag */
	fStore = false;
}

void ShapeFunctionT::Store(void)
{
	/* check */
	if (!fCurrElementNumber)
		ExceptionT::GeneralFail("ShapeFunctionT::Store", "current element not set");

	/* store Jocobian information */
	fDet_store.SetRow(*fCurrElementNumber, fDet);

	/* store shape function derivatives */
	for (int i = 0; i < fDNaX_store.Length(); i++)
		fDNaX_store[i].SetRow(*fCurrElementNumber, fDNaX[i]);
}

void ShapeFunctionT::CloseStore(void)
{
	/* set flag */
	fStore = true;
}

void ShapeFunctionT::TransformDerivatives(const dMatrixT& changeofvar,
	const dArray2DT& original, dArray2DT& transformed) const
{
	int  numnodes = original.MinorDim();

	/* allocate memory */
	transformed.Dimension(original.MajorDim(),numnodes);

	/* not so const workspace */
	dArrayT& v1 = const_cast<dArrayT&>(fv1);
	dArrayT& v2 = const_cast<dArrayT&>(fv2);

	/* apply chain rule derivative */
	for (int i = 0; i < numnodes; i++)
	{
		/* fetch values */
		original.ColumnCopy(i,v1);

		/* transform */
		changeofvar.MultTx(v1,v2);

		/* write back */
		transformed.SetColumn(i,v2);
	}
}

/***********************************************************************
* Private
***********************************************************************/

/* configure work space arrays - initializes shape function to be
* isoparametric */
void ShapeFunctionT::Construct(void)
{
	/* check local array type (ambiguous for linear geometry) */
	if (fCoords.Type() != LocalArrayT::kInitCoords &&
	    fCoords.Type() != LocalArrayT::kCurrCoords)
	    ExceptionT::GeneralFail("ShapeFunctionT::Construct");

	/* dimensions */
	int numXnodes = fCoords.NumberOfNodes();
	int numUnodes = numXnodes; // assume isoparametric
	int numsd     = fCoords.MinorDim();

	/* parent domain jacobian */
	fDet.Dimension(fNumIP),

	/* memory for the derivatives */
	fDNaX.Dimension(fNumIP);
	for (int i = 0; i < fNumIP; i++)
		fDNaX[i].Dimension(numsd, numXnodes);

	/* initialize to isoparametric */
	pNaU  = &(fDomain->Na());
	pDNaU = &fDNaX;

	/* work space */
	fv1.Dimension(numsd);
	fv2.Dimension(numsd);
}

void ShapeFunctionT::Construct_DN_DDN(void)
{
	/* check local array type (ambiguous for linear geometry) */
	if (fCoords.Type() != LocalArrayT::kInitCoords &&
	    fCoords.Type() != LocalArrayT::kCurrCoords)
	    ExceptionT::GeneralFail("ShapeFunctionT::Construct");

	/* dimensions */
	int numXnodes = fCoords.NumberOfNodes();
	int numUnodes = numXnodes; // assume isoparametric
	int numsd     = fCoords.MinorDim();

	/* parent domain jacobian */
	fDet.Dimension(fNumIP),

	/* memory for the first derivatives */
	fDNaX.Dimension(fNumIP);
	for (int i = 0; i < fNumIP; i++)
		fDNaX[i].Dimension(numsd, numXnodes);

	/* memory for the second derivatives */
	fDDNaX.Dimension(fNumIP);
	for (int i = 0; i < fNumIP; i++)
		fDDNaX[i].Dimension(numsd*2, numXnodes);

	/* initialize to isoparametric */
	pNaU   = &(fDomain->Na());
	pDNaU  = &fDNaX;
	pDDNaU = &fDDNaX;

	/* work space */
	fv1.Dimension(numsd);
	fv2.Dimension(numsd);
}

//
// compute jacobians wrt parametric space at integration points
//
void ShapeFunctionT::JacobianDets(dArrayT& dets)
{

  const int nsd = NumSD();
  const int nip = NumIP();

  dets.Dimension(nip);

  fDomain->ComputeDNa(fCoords, fDNaX, dets);

}
