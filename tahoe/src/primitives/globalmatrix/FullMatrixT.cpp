/* $Id: FullMatrixT.cpp,v 1.24 2011/12/01 21:11:40 bcyansfn Exp $ */
/* created: paklein (03/07/1998) */
#include "FullMatrixT.h"
#include <iostream>
#include <iomanip>
#include "toolboxConstants.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "ElementMatrixT.h"
#include "StringT.h"
#include "ofstreamT.h"
#include "CommunicatorT.h"

using namespace Tahoe;

/* constructor */
FullMatrixT::FullMatrixT(ostream& out,int check_code, const CommunicatorT& comm):
	GlobalMatrixT(out, check_code, comm),
	fIsFactorized(false)
{

}

/* copy constructor */
FullMatrixT::FullMatrixT(const FullMatrixT& source):
	GlobalMatrixT(source),
	fMatrix(source.fMatrix),
	fIsFactorized(source.fIsFactorized)
{

}

/* set the internal matrix structure.
* NOTE: do not call Initialize() equation topology has been set
* with AddEquationSet() for all equation sets */
void FullMatrixT::Initialize(int tot_num_eq, int loc_num_eq, int start_eq)
{
	/* inherited */
	GlobalMatrixT::Initialize(tot_num_eq, loc_num_eq, start_eq);

	/* check */
	if (tot_num_eq != loc_num_eq)
		ExceptionT::GeneralFail("FullMatrixT::Initialize",
			"total equations %d != local equations %d", tot_num_eq, loc_num_eq);

	/* allocate work space */
	fMatrix.Dimension(fLocNumEQ);
	fIsFactorized = false;
}

/* set all matrix values to 0.0 */
void FullMatrixT::Clear(void)
{
	/* inherited */
	GlobalMatrixT::Clear();
	
	/* clear values */
	fMatrix = 0.0;	

	/* set flag */
	fIsFactorized = false;
}

/* add element group equations to the overall topology.
* NOTE: assembly positions (equation numbers) = 1...fNumEQ
* equations can be of fixed size (iArray2DT) or
* variable length (RaggedArray2DT) */
void FullMatrixT::AddEquationSet(const iArray2DT& eqset)
{
#pragma unused(eqset)
// no equation data needed
}

void FullMatrixT::AddEquationSet(const RaggedArray2DT<int>& eqset)
{
#pragma unused(eqset)
// no equation data needed
}

/* assemble the element contribution into the LHS matrix - assumes
* that elMat is square (n x n) and that eqnos is also length n.
* NOTE: assembly positions (equation numbers) = 1...fNumEQ */
void FullMatrixT::Assemble(const ElementMatrixT& elMat, const ArrayT<int>& eqnos)
{
#if __option(extended_errorcheck)
	/* dimension checking */
	if (elMat.Rows() != eqnos.Length() ||
	    elMat.Cols() != eqnos.Length()) throw ExceptionT::kSizeMismatch;
#endif

	/* element matrix format */
	ElementMatrixT::FormatT format = elMat.Format();

	if (format == ElementMatrixT::kDiagonal)
	{
		/* from diagonal only */
		const double* pelMat = elMat.Pointer();
		int inc = elMat.Rows() + 1;
	
		const int* peq = eqnos.Pointer();		
		for (int i = 0; i < elMat.Length(); i++)
		{
			int eq = *peq++;
			
			/* active dof's only */
			if (eq-- > 0) fMatrix(eq,eq) += *pelMat;
			
			pelMat += inc;
		}
	}
	else
	{
		/* copy to full symmetric */
		if (elMat.Format() == ElementMatrixT::kSymmetricUpper) elMat.CopySymmetric();

		const int* peq = eqnos.Pointer();	
		for (int col = 0; col < elMat.Cols(); col++)
		{
			int eqc = eqnos[col];
			const double* pelMat = elMat(col);
		
			/* active dof's only */
			if (eqc-- > 0)
			{		
				const int* peqr = eqnos.Pointer();
				for (int row = 0; row < elMat.Rows(); row++)
				{
					int eqr = *peqr++;
				
					/* active dof's only */
					if (eqr-- > 0) fMatrix(eqr,eqc) += *pelMat;
					
					pelMat++;
				}
			}
		}
	}
}

void FullMatrixT::Assemble(const ElementMatrixT& elMat, const ArrayT<int>& row_eqnos,
	const ArrayT<int>& col_eqnos)
{
#if __option(extended_errorcheck)
	/* dimension check */
	if (elMat.Rows() != row_eqnos.Length() ||
	    elMat.Cols() != col_eqnos.Length()) throw ExceptionT::kSizeMismatch;
#endif

	/* element matrix format */
	ElementMatrixT::FormatT format = elMat.Format();

	if (format == ElementMatrixT::kDiagonal)
	{
		cout << "\n FullMatrixT::Assemble(m, r, c): cannot assemble diagonal matrix" << endl;
		throw ExceptionT::kGeneralFail;
	}
	else
	{
		/* copy to full symmetric */
		if (elMat.Format() == ElementMatrixT::kSymmetricUpper) elMat.CopySymmetric();

		/* assemble active degrees of freedom */
		int n_c = col_eqnos.Length();
		int n_r = row_eqnos.Length();
		for (int col = 0; col < n_c; col++)
		{
			int ceqno = col_eqnos[col] - 1;	
			if (ceqno > -1)	
				for (int row = 0; row < n_r; row++)
				{
					int reqno = row_eqnos[row] - 1;
					if (reqno > -1)
						fMatrix(reqno,ceqno) += elMat(row,col);
				}
		}
	}
}

void FullMatrixT::Assemble(const nArrayT<double>& diagonal_elMat, const ArrayT<int>& eqnos)
{
#if __option(extended_errorcheck)
	/* dimension check */
	if (diagonal_elMat.Length() != eqnos.Length()) throw ExceptionT::kSizeMismatch;
#endif

	for (int i = 0; i < eqnos.Length(); i++)
	{
		int eqno = eqnos[i] - 1;
		if (eqno > -1)
			fMatrix(eqno,eqno) += diagonal_elMat[i];
	}
}

/* strong manipulation functions */
void FullMatrixT::OverWrite(const ElementMatrixT& elMat, const nArrayT<int>& eqnos)
{
#if __option(extended_errorcheck)
	/* dimension checking */
	if (elMat.Rows() != eqnos.Length() ||
	    elMat.Cols() != eqnos.Length()) throw ExceptionT::kSizeMismatch;
#endif

	/* copy to full symmetric */
	if (elMat.Format() == ElementMatrixT::kSymmetricUpper) elMat.CopySymmetric();

	const int* peq = eqnos.Pointer();	
	for (int col = 0; col < elMat.Cols(); col++)
	{
		int eqc = eqnos[col];
		const double* pelMat = elMat(col);
		
		/* active dof's only */
		if (eqc-- > 0)
		{		
			const int* peqr = eqnos.Pointer();
			for (int row = 0; row < elMat.Rows(); row++)
			{
				int eqr = *peqr++;
			
				/* active dof's only */
				if (eqr-- > 0) fMatrix(eqr,eqc) += *pelMat;
				
				pelMat++;
			}
		}
	}
}

void FullMatrixT::Disassemble(dMatrixT& elMat, const nArrayT<int>& eqnos) const
{
#if __option(extended_errorcheck)
	/* dimension checking */
	if (elMat.Rows() != eqnos.Length() ||
	    elMat.Cols() != eqnos.Length()) throw ExceptionT::kSizeMismatch;
#endif

	const int* peq = eqnos.Pointer();
	for (int col = 0; col < elMat.Cols(); col++)
	{
		int        eqc = eqnos[col];
		double* pelMat = elMat(0);
		
		/* active dof's only */
		if (eqc-- > 0)
		{		
			const int* peqr = eqnos.Pointer();
			for (int row = 0; row < elMat.Rows(); row++)
			{
				int eqr = *peqr++;
			
				/* active dof's only */
				if (eqr-- > 0)
					*pelMat = fMatrix(eqr,eqc);
				else
					*pelMat = 0.0;
				
				pelMat++;
			}
		}
		else
			/* clear column */
			elMat.SetCol(col, 0.0);
	}
}

void FullMatrixT::DisassembleDiagonal(dArrayT& diagonals, const nArrayT<int>& eqnos) const
{
#if __option(extended_errorcheck)
	/* dimension checking */
	if (diagonals.Length() != eqnos.Length()) throw ExceptionT::kSizeMismatch;
#endif

	for (int i = 0; i < eqnos.Length(); i++)
	{
		/* ignore requests for inactive equations */	
		int eq = eqnos[i];	
		if (eq-- > 0)
			diagonals[i] = fMatrix(eq,eq);
		else
			diagonals[i] = 0.0;
	}
}

/* number scope and reordering */
GlobalMatrixT::EquationNumberScopeT FullMatrixT::EquationNumberScope(void) const
{
	return kLocal;
}

bool FullMatrixT::RenumberEquations(void) const { return false; }

/* assignment operator */
FullMatrixT& FullMatrixT::operator=(const FullMatrixT& rhs)
{
	/* no copies of self */
	if (&rhs == this) return *this;

	/* inherited */
	GlobalMatrixT::operator=(rhs);

	fMatrix = rhs.fMatrix;
	fIsFactorized = rhs.fIsFactorized;

	return *this;
}
	
/* return a clone of self. Caller is responsible for disposing of the matrix */
GlobalMatrixT* FullMatrixT::Clone(void) const
{
	FullMatrixT* new_mat = new FullMatrixT(*this);
	return new_mat;
}

void FullMatrixT::Multx(const dArrayT& x, dArrayT& b) const
{
	/* already factorized */
	if (fIsFactorized)
		ExceptionT::GeneralFail("FullMatrixT::Multx", "matrix is factorized");

	/* calculate product */
	fMatrix.Multx(x, b);
}

void FullMatrixT::MultTx(const dArrayT& x, dArrayT& b) const
{
	/* already factorized */
	if (fIsFactorized)
		ExceptionT::GeneralFail("FullMatrixT::MultTx", "matrix is factorized");

	/* calculate product */
	fMatrix.MultTx(x, b);
}

/* vector-matrix-vector product */
double FullMatrixT::MultmBn(const dArrayT& m, const dArrayT& n) const
{
	/* already factorized */
	if (fIsFactorized)
		ExceptionT::GeneralFail("FullMatrixT::MultmBn", "matrix is factorized");

	/* calculate product */
	return fMatrix.MultmBn(m, n);
}

/**************************************************************************
 * Protected
 **************************************************************************/

/* determine new search direction and put the results in result */
void FullMatrixT::BackSubstitute(dArrayT& result)
{
	if (fIsFactorized)
		ExceptionT::GeneralFail("FullMatrixT::BackSubstitute", "no multiple solves");

	fMatrix.LinearSolve(result);
	fIsFactorized = true;
}

/* rank check functions */
void FullMatrixT::PrintAllPivots(void) const
{
//TEMP: no full, nonsymmetric factorization implemented
}

void FullMatrixT::PrintZeroPivots(void) const
{
//TEMP: no full, nonsymmetric factorization implemented
}

void FullMatrixT::PrintLHS(bool force) const
{
	if (!force && fCheckCode != GlobalMatrixT::kPrintLHS) return;

	/* output stream */
	StringT file = fstreamT::Root();
	if (fPrintTag != "") 
		file.Append(fPrintTag, sOutputCount);
	else
		file.Append("FullMatrixT.LHS.", sOutputCount);
	if (fComm.Size() > 1) file.Append(".p", fComm.Rank());	
	ofstreamT out(file);
	out.precision(14);

	/* write non-zero values in RCV format */
	for (int r = 0; r < fMatrix.Rows(); r++)
		for (int c = 0; c < fMatrix.Cols(); c++)
			if (fMatrix(r,c) != 0.0)
				out << r+1 << " " << c+1 << " " << fMatrix(r,c) << '\n';
	
	/* increment count */
	sOutputCount++;
}
