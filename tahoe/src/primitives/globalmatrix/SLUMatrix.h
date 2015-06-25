/* $Id: SLUMatrix.h,v 1.14 2005/04/13 21:49:58 paklein Exp $ */
/* created: rbridson (06/30/2000) */
#ifndef _SLU_MATRIX_H_
#define _SLU_MATRIX_H_

/* project headers */
#include "Environment.h"

#undef __SUPERLU__
//#define __SUPERLU__

/* library support */
#ifdef __SUPERLU__

/* base class */
#include "GlobalMatrixT.h"

/* direct members */
#include "LinkedListT.h"

/* SuperLU 3.0 types and function prototypes */
#include "dsp_defs.h"

namespace Tahoe {

/** interface to SuperLU 3.0 serial linear solver */
class SLUMatrix: public GlobalMatrixT
{
public:

	/** constructor */
	SLUMatrix(ostream& out, int check_code, const CommunicatorT& comm);

	/** destructor */
	~SLUMatrix(void);

	/* set the internal matrix structure.
	 * NOTE: do not call Initialize() until equation topology has been set
	 * with AddEquationSet() for all equation sets */
	virtual void Initialize(int tot_num_eq, int loc_num_eq, int start_eq);

	/** write information to output stream after SLUMatrix::Initialize
	 * has been called */
	virtual void Info(ostream& out);
	
	/* set all matrix values to 0.0 */
	virtual void Clear(void);
	
	/* add element group equations to the overall topology.
	 * NOTE: assembly positions (equation numbers) = 1...fNumEQ
	 * equations can be of fixed size (iArray2DT) or
	 * variable length (RaggedArray2DT) */
	virtual void AddEquationSet(const iArray2DT& eqset);
	virtual void AddEquationSet(const RaggedArray2DT<int>& eqset);
	
	/* assemble the element contribution into the LHS matrix - assumes
	 * that elMat is square (n x n) and that eqnos is also length n.
	 * NOTE: assembly positions (equation numbers) = 1...fNumEQ */
	virtual void Assemble(const ElementMatrixT& elMat, const ArrayT<int>& eqnos);
	virtual void Assemble(const ElementMatrixT& elMat, const ArrayT<int>& row_eqnos,
		const ArrayT<int>& col_eqnos);
	virtual void Assemble(const nArrayT<double>& diagonal_elMat, const ArrayT<int>& eqnos);
	
	/* element accessor - READ ONLY */
	double Element(int row, int col) const;

	/* number scope and reordering */
	virtual EquationNumberScopeT EquationNumberScope(void) const;
	virtual bool RenumberEquations(void) const;

	/** return the form of the matrix */
	virtual GlobalT::SystemTypeT MatrixType(void) const { return GlobalT::kNonSymmetric; };

	/** assignment operator */
	SLUMatrix& operator=(const SLUMatrix&);
	
	/** return a clone of self */
	virtual GlobalMatrixT* Clone(void) const;

protected:

	/** \name tuning parameters */
	/*@{*/
	static const int kSupernodeRelax = 5;
	static const int kPanelSize = 10;
	/*@}*/

	/** output in sparse format */
	friend ostream& operator<<(ostream& out, const SLUMatrix& matrix) {
		NCformat *A = (NCformat*) (matrix.fA.Store);
		int i, j;
		for (i = 0; i < matrix.fA.ncol; ++i) {
			out << -1 << " " << 0.0 << "\n";
			for (j = A->colptr[i]; j < A->colptr[i+1]; ++j) {
				out << A->rowind[j]+1 << " ";
				out << ((double*)A->nzval)[j] << "\n";
			}
		}
		return out;
	};

	/** factorize the matrix. SLUMatrix::Factorize does nothing because an
	 * all-in-one driver provided with SuperLU 3.0 is called in 
	 * SLUMatrix::BackSubstitute */
	virtual void Factorize(void);
	
	/** solution driver. Calls all-in-one driver provided with SuperLU 3.0 which 
	 * can be called for solving multiple right-hand sides or just resolving
	 * a matrix with the same sparsity pattern as a previous solve. This driver
	 * routine is adapted from dlinsolx2.c provided in the SuperLU 3.0 examples. */
	virtual void BackSubstitute(dArrayT& result);

	/** \name check functions */
	/*@{*/
	virtual void PrintAllPivots(void) const;
	virtual void PrintZeroPivots(void) const;
	virtual void PrintLHS(void) const;
	/*@}*/

	/** element accessor - read and write, for assembly. Exception for */
	/* access to unstored zero. */
	double& operator()(int row, int col);

	/** (over)estimate number of nonzeros from equation sets */
	virtual void EstimateNNZ (int *colLength, int &nnz);

	/** insert all the element equations into A */
	virtual void InsertEquations (NCformat *A, int *colLength, int &nnz);

	/** insert the list nzlist of nonzeros (with length nzlen) into column c
	 * of matrix A. Use colLengths to keep track of column lengths, and
	 * keeping nnz up-to-date. The columns of A will have nonzeros in
	 * ascending order. */
	virtual void InsertNZ (NCformat *A, int *colLength, int &nnz, int c,
	   int nzlen, int *nzlist);

	/** compress columns in A */
	virtual void CompressColumns (NCformat *A, const int *colLength);

private:

	/** no copy constructor */
	SLUMatrix(const SLUMatrix&);

protected:

	/** \name matrix and its factors in SuperLU column format (Harwell-Boeing) */
	/*@{*/
	SuperMatrix fA;
	SuperMatrix fL;
	SuperMatrix fU;

	/** flag to see if SLUMatrix::fL and SLUMatrix::fU have been allocated */
	bool fLUallocated;
	/*@}*/

	/** \name factorization workspace */
	/*@{*/
	/** column permutations */
	int* fperm_c;

	/** row permutations */
	int* fperm_r;
	
	/** true if fperm_c has been set */
	bool fIsColOrdered;

	/** symbolic information used in factorization */
	int *fetree;
	/*@}*/

	/** SuperLU options */
	superlu_options_t foptions;

	/** \name equation sets */
	/*@{*/
	LinkedListT<const iArray2DT*> fEqnos;
	LinkedListT<const RaggedArray2DT<int>*> fRaggedEqnos;
	/*@}*/
};

/* element accessor - read and write, for assembly. Exception for */
/* access to unstored zero. row and col are */
inline double& SLUMatrix::operator()(int row, int col)
{
	const char caller[] = "SLUMatrix::operator()";

#if __option(extended_errorcheck)
	/* range checks */
	if (row < 0 || row >= fLocNumEQ) ExceptionT::GeneralFail(caller);
	if (col < 0 || col >= fLocNumEQ) ExceptionT::GeneralFail(caller);
#endif

	NCformat *A = (NCformat*) fA.Store;

	/* look through column col for the given row index */
	int start = A->colptr[col];
	int stop  = A->colptr[col+1];
	for (int i = start; i < stop; ++i)
	{
		if (row == A->rowind[i])
			return ((double*)A->nzval)[i];
	}

	/* otherwise this nonzero wasn't present */
	ExceptionT::OutOfRange(caller, "(%d,%d) not present", row, col);
	
	/* dummy */
	return ((double*)A->nzval)[0];
}

/* element accessor - READ ONLY */
inline double SLUMatrix::Element(int row, int col) const
{
	const char caller[] = "SLUMatrix::Element()";

#if __option(extended_errorcheck)
	/* range checks */
	if (row < 0 || row >= fLocNumEQ) ExceptionT::GeneralFail(caller);
	if (col < 0 || col >= fLocNumEQ) ExceptionT::GeneralFail(caller);
#endif

	NCformat *A = (NCformat*) fA.Store;

	/* look through column col for the given row index */
	int start = A->colptr[col];
	int stop  = A->colptr[col+1];
	for (int i = start; i < stop; ++i)
	{
		if (row == A->rowind[i])
			return ((double*)A->nzval)[i];
	}

	/* otherwise this nonzero wasn't present */
	return 0.0;
}

} /* namespace Tahoe */
#endif /* __SUPERLU__ */
#endif /* _SLU_MATRIX_H_ */
