/* MUMPSMatrixT.h — MUMPS sequential sparse direct solver wrapper for Tahoe */
#ifndef _MUMPS_MATRIX_T_H_
#define _MUMPS_MATRIX_T_H_

/* base class — same as SPOOLESMatrixT */
#include "MSRMatrixT.h"

#ifdef __MUMPS__

#include "iArrayT.h"
#include "dArrayT.h"
#include "dmumps_c.h"

namespace Tahoe {

/** Wrapper around the MUMPS sequential (no-MPI) sparse direct solver.
 *
 *  Inherits MSRMatrixT for sparse assembly, then converts to COO triplet
 *  format (irn/jcn/a) for MUMPS via MSRMatrixT::GenerateRCV().  The index
 *  shift from 0-based (Tahoe) to 1-based (MUMPS/Fortran) is done in
 *  Factorize().
 *
 *  Usage in XML input file:
 *    <MUMPS_matrix message_level="silent" always_symmetric="false" ordering="PORD"/>
 *
 *  This is the default linear solver when a deck names no matrix (#80).
 *  The symbolic analysis is reused while the sparsity pattern of the
 *  assembled matrix is unchanged, so repeated Newton iterations only
 *  refactorize numerically.
 */
class MUMPSMatrixT: public MSRMatrixT
{
public:

    /** MUMPS ICNTL(7) fill-reducing orderings */
    enum OrderingT { kAMD = 0, kAMF = 2, kSCOTCH = 3, kPORD = 4,
                     kMETIS = 5, kQAMD = 6, kAuto = 7 };

    /** name of an ordering, for output */
    static const char* OrderingName(int ordering);

    MUMPSMatrixT(ostream& out, int check_code, bool symmetric,
        int message_level, int ordering, const CommunicatorT& comm);

    /** copy constructor — not supported */
    MUMPSMatrixT(const MUMPSMatrixT& source);

    virtual ~MUMPSMatrixT(void);

    /** MUMPS solve does not overwrite the assembled matrix values */
    virtual bool SolvePreservesData(void) const { return true; }

    /** clear values for next assembly and reset the MUMPS instance */
    virtual void Clear(void);

    /** return the form of the matrix (always non-symmetric storage) */
    virtual GlobalT::SystemTypeT MatrixType(void) const { return GlobalT::kNonSymmetric; }

    /** return a clone of self */
    virtual GlobalMatrixT* Clone(void) const;

protected:

    /** LU factorize the assembled matrix: MUMPS job 4 (analysis and
     *  factorization) when the sparsity pattern changed, job 2 otherwise */
    virtual void Factorize(void);

    /** back-substitute: result is overwritten with the solution (MUMPS job 3) */
    virtual void BackSubstitute(dArrayT& result);

    virtual void PrintAllPivots(void) const {}
    virtual void PrintZeroPivots(void) const {}
    virtual void PrintLHS(bool force = false) const {}

private:

    void Initialize(void);   /* MUMPS job -1 */
    void Finalize(void);     /* MUMPS job -2 */

    DMUMPS_STRUC_C fId;      /**< MUMPS control structure */
    /** run one MUMPS job; returns false on error (infog[0] < 0) */
    bool RunJob(int job);

    int  fMessageLevel;
    int  fOrdering;      /**< requested ICNTL(7) */
    bool fSymmetric;
    bool fIsInitialized;
    bool fIsFactorized;
    bool fIsAnalyzed;    /**< fRowIdx/fColIdx hold the analyzed pattern */
    bool fInitedMPI;     /**< true if this object called MPI_Init */

    /** 1-based COO storage kept alive between Factorize and BackSubstitute */
    iArrayT fRowIdx;   /**< irn: row indices (1-based) */
    iArrayT fColIdx;   /**< jcn: column indices (1-based) */
    dArrayT fValues;   /**< a:   nonzero values */
};

} // namespace Tahoe

#endif /* __MUMPS__ */
#endif /* _MUMPS_MATRIX_T_H_ */
