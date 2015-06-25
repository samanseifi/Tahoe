/* $Id: EpetraCRSMatrixT.h,v 1.4 2006/11/25 22:06:11 paklein Exp $ */
#ifndef _EPETRA_CRS_MATRIX_T_H_
#define _EPETRA_CRS_MATRIX_T_H_

/* library support */
#ifdef __TRILINOS__

/* base class */
#include "GlobalMatrixT.h"
#include "Epetra_CrsMatrix.h"

/* direct members */
#include "AutoArrayT.h"

namespace Tahoe {

/* forward declarations */
class MSRBuilderT;

/** interface for solvers in Trilinos */
class EpetraCRSMatrixT: public GlobalMatrixT
{
public:

	/** constructor */
	EpetraCRSMatrixT(ostream& out, int check_code, const CommunicatorT& comm);

	/** copy constructor */
	EpetraCRSMatrixT(const EpetraCRSMatrixT& rhs);

	/** destructor */
	~EpetraCRSMatrixT(void);
	
	/** translate this to an Epetra_CrsMatrix. Requestor is responsible for
	 * freeing the returned object */
	Epetra_CrsMatrix* Translate(void) const;

	/** return a pointer to the processor map */
	Epetra_Map* Map(void) { return fepetra_map; }

	/** set the internal matrix structure */
	virtual void Initialize(int tot_num_eq, int loc_num_eq, int start_eq);

	/** write information to output stream after EpetraCRSMatrixT::Initialize
	/* set all matrix values to 0.0 */
	virtual void Clear(void);
	
	/** \name set matrix structure
	 * active equations are eq > 0 */
	/*@{*/
	void AddEquationSet(const iArray2DT& eqnos);
	void AddEquationSet(const RaggedArray2DT<int>& eqnos);
	/*@}*/

	/** \name assembly operators
	 * active equations are eq > 0 */	
	/*@{*/
	virtual void Assemble(const ElementMatrixT& elMat, const ArrayT<int>& eqnos);
	virtual void Assemble(const ElementMatrixT& elMat, const ArrayT<int>& row_eqnos,
		const ArrayT<int>& col_eqnos);
	virtual void Assemble(const nArrayT<double>& diagonal_elMat, const ArrayT<int>& eqnos);
	/*@}*/

	/* number scope and reordering */
	virtual EquationNumberScopeT EquationNumberScope(void) const;
	virtual bool RenumberEquations(void) const;

	/** return the form of the matrix */
	virtual GlobalT::SystemTypeT MatrixType(void) const { return GlobalT::kNonSymmetric; };

	/** assignment operator */
	EpetraCRSMatrixT& operator=(const EpetraCRSMatrixT& rhs);

	/** return a clone of self */
	virtual GlobalMatrixT* Clone(void) const;

protected:

	/** solution driver. Not implemented. Added only because BackSubstitute is a pure virtual
	 * function, and we would like to use this class to generate translate to Epetra_CrsMatrix  */
	virtual void BackSubstitute(dArrayT& result);

	/** \name check functions */
	/*@{*/
	virtual void PrintAllPivots(void) const;
	virtual void PrintZeroPivots(void) const;
	virtual void PrintLHS(bool force) const;
	/*@}*/

	/** element accessor. Returns a pointer to the given element in the matrix
	 * or NULL if the element is not in the set of non-zero elements. The search
	 * for columns is performed through recursive bisection of the indicies of
	 * columns with non-zero elements in the specified row. */
	double* operator()(int row, int col);

protected:

	/** matrix structure builder */
	MSRBuilderT* fBuilder;

	/** \name Epetra components */
	/*@{*/
	Epetra_Comm* fepetra_comm;
	Epetra_Map* fepetra_map;
	/*@}*/

	/** \name matrix and factors in SuperLU formats */
	/*@{*/
	AutoArrayT<int> factive;
	AutoArrayT<int> frowptr;
	AutoArrayT<int> fcolind;
	AutoArrayT<double> fnzval;
	/*@}*/

	/** \name factorization flags */
	/*@{*/
	/** true if matrix has been symbolically factorized */
	bool fIsSymFactorized;

	/** true if matrix has been numerically factorized */
	bool fIsNumFactorized;
	/*@}*/
};

/* element accessor */
inline double* EpetraCRSMatrixT::operator()(int row, int col)
{
	const char caller[] = "EpetraCRSMatrixT::operator()";
	int loc_row = row - fStartEQ + 1; /* fStartEQ is 1,... */

#if __option(extended_errorcheck)
	/* range checks */
	if (loc_row < 0 || loc_row >= fLocNumEQ) ExceptionT::OutOfRange(caller, "row_loc %d < 0 or >= %d", loc_row, fLocNumEQ);
	if (col < 0 || col >= fTotNumEQ) ExceptionT::OutOfRange(caller, "col %d < 0 || >= %d", col, fTotNumEQ);
#endif

	/* non-zero columns */
	int r_dex = frowptr[loc_row];
	int* pcolind = fcolind.Pointer(r_dex);

	/* range */
	int min = 0;
	int max = frowptr[loc_row+1] - r_dex - 1;

	/* is last value */
	if (pcolind[max] == col)
		return fnzval.Pointer(r_dex + max);

	/* bisection */
	int c_dex = (max + min)/2;
	int c_hit = pcolind[c_dex];
	while (c_hit != col && max != min+1) {
	
		/* shift bounds */
		if (c_hit > col)
			max = c_dex;
		else /* painds[cdex] < col */
			min = c_dex;
	
		/* bisect */
		c_dex = (max + min)/2;
		c_hit = pcolind[c_dex];
	}
	
	/* found */
	if (c_hit == col)
		return fnzval.Pointer(r_dex + c_dex);
	else /* not found */
		return NULL;
}

} /* namespace Tahoe */

#endif /* __TRILINOS__ */
#endif /* _EPETRA_CRS_MATRIX_T_H_ */
