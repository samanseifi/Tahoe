/* $Id: D2MeshFreeFSSolidT.h,v 1.7 2002/11/14 17:05:54 paklein Exp $ */
/* created: paklein (10/23/1999) */
#ifndef _D2_EFG_FDELASTIC_T_H_
#define _D2_EFG_FDELASTIC_T_H_

/* base class */
#include "MeshFreeFSSolidT.h"

namespace Tahoe {

/* forward declaration */
class D2MeshFreeShapeFunctionT;
class D2VIB2D;

class D2MeshFreeFSSolidT: public MeshFreeFSSolidT
{
public:

	/* constructor */
	D2MeshFreeFSSolidT(const ElementSupportT& support, const FieldT& field);

	/* accessors */
	const D2MeshFreeShapeFunctionT& D2MLSShapeFunction() const;

	/* check material's list */
	virtual void Initialize(void);

//DEV - no need to override
#if 0
	/* form the residual force vector */
	virtual void RHSDriver(void);
	void ElementRHSDriver(void);
#endif

protected:

	/* initialization functions */
	virtual void SetShape(void);

	/* increment current element */
	virtual bool NextElement(void);

	/* calculate the internal force contribution */
	virtual void FormKd(double constK);

private:

	/* 3rd rank tensor (double) contraction */
	void A_ijk_B_jkl(const dMatrixT& A, const dMatrixT& B, dMatrixT& C);

	/* write displacement field and gradients */
	 virtual void WriteField(void);

protected:

	/* shape functions */
	D2MeshFreeShapeFunctionT* fD2MFShapes;

	/* one shot wonder for now */
	D2VIB2D* pD2VIB2D;

	/* work space */
	dMatrixT fDW;  // like PK1
	dMatrixT fDDW; // gradient stress term
	dMatrixT fD2GradNa;
	nVariMatrixT<double> fD2GradNa_wrap;
	
	//debugging flag
	int DoPrint;
};

/* inline  */
inline const D2MeshFreeShapeFunctionT& D2MeshFreeFSSolidT::D2MLSShapeFunction() const
{
	if (!fD2MFShapes) throw ExceptionT::kGeneralFail;
	return *fD2MFShapes;
}

} // namespace Tahoe 
#endif /* _D2_EFG_FDELASTIC_T_H_ */
