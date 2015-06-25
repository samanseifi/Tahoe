/* $Id: D3MeshFreeShapeFunctionT.h,v 1.10 2005/12/23 14:15:57 kyonten Exp $ */
/* created: paklein (10/23/1999) */
#ifndef _D3_MF_SHAPE_T_H_
#define _D3_MF_SHAPE_T_H_

/* base class */
#include "D2MeshFreeShapeFunctionT.h"

namespace Tahoe {

/* forward declarations */
class D3MeshFreeSupportT;

class D3MeshFreeShapeFunctionT: public D2MeshFreeShapeFunctionT
{
public:

	/* constructors */
	D3MeshFreeShapeFunctionT(GeometryT::CodeT geometry_code, int numIP,
		const LocalArrayT& coords, const dArray2DT& all_coords,
		const iArray2DT& connects, const iArrayT& nongridnodes,
		const int& currelement, const ParameterListT& mf_support_params);

	/** class-dependent initializations */
	virtual void Initialize(void);

	/* compute global shape derivatives */ 	
	virtual void SetDerivatives(void);
	int SetDerivativesAt(const dArrayT& x, AutoArrayT<int>& nodes);
		// returns 0 if MLS fit fails

	/* 3rd order shape function gradients matrix */
	void D3GradNa(dMatrixT& D3_grad_Na) const;

	/* 3rd order spatial gradients */
	void GradGradGradU(const LocalArrayT& nodal, dMatrixT& gradgradgrad_U) const;
	void GradGradGradU(const LocalArrayT& nodal, dMatrixT& gradgradgrad_U, int ip) const;

	/* 3rd derivatives of shape functions at IP */
	const dArray2DT& DDDerivatives_U(int ip) const { return fDDDNaU[ip]; };
	const dArray2DT& DDDerivatives_U(void) const { return fDDDNaU[fCurrIP]; };

	/* reconstruct displacement field and all derivatives */
	void NodalField(const dArray2DT& DOF, dArray2DT& field, dArray2DT& Dfield,
		dArray2DT& DDfield, dArray2DT& DDDfield, iArrayT& nodes);
		
	/* reference to the support */
	D3MeshFreeSupportT& D3MeshFreeSupport(void) const;

protected:

	/* meshfree database support */
	D3MeshFreeSupportT* fD3MFSupport;

	ArrayT<dArray2DT> fDDDNaU;
	
	/* work space for blended shape functions */
	ArrayT<dArray2DT> fDDDNa_tmp;
	
	/* number of derivatives to calculate */
	int fNumDeriv;
	
};

/* inlines */
inline D3MeshFreeSupportT& D3MeshFreeShapeFunctionT::D3MeshFreeSupport(void) const
{
	if (!fD3MFSupport) throw ExceptionT::kGeneralFail;
	return *fD3MFSupport;
}

/* spatial gradients */
inline void D3MeshFreeShapeFunctionT::GradGradGradU(const LocalArrayT& nodal,
	dMatrixT& gradgradgrad_U) const
{
#if __option(extended_errorcheck)
	if (gradgradgrad_U.Rows() != nodal.MinorDim() ||
		gradgradgrad_U.Cols() != fDDDNaU[fCurrIP].MajorDim())
    	throw ExceptionT::kSizeMismatch;
#endif

	fDomain->Jacobian(nodal, fDDDNaU[fCurrIP], gradgradgrad_U);	
}

inline void D3MeshFreeShapeFunctionT::GradGradGradU(const LocalArrayT& nodal,
	dMatrixT& gradgradgrad_U, int ip) const
{
#if __option(extended_errorcheck)
	if (gradgradgrad_U.Rows() != nodal.MinorDim() ||
		gradgradgrad_U.Cols() != fDDDNaU[ip].MajorDim())
    	throw ExceptionT::kSizeMismatch;
#endif
	
	fDomain->Jacobian(nodal, fDDDNaU[ip], gradgradgrad_U);	
}

} // namespace Tahoe 
#endif /* _D3_MF_SHAPE_T_H_ */
