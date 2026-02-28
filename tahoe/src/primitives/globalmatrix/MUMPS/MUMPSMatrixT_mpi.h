/* MUMPSMatrixT_mpi.h — MUMPS MPI distributed sparse direct solver for Tahoe */
#ifndef _MUMPS_MATRIX_T_MPI_H_
#define _MUMPS_MATRIX_T_MPI_H_

#include "MSRMatrixT.h"

#if defined(__MUMPS__) && defined(__TAHOE_MPI__)

#include "iArrayT.h"
#include "dArrayT.h"
#include "dmumps_c.h"
/* mpi.h is included transitively via CommunicatorT.h when __TAHOE_MPI__ is set */

namespace Tahoe {

/** Wrapper around the MUMPS parallel (MPI) sparse direct solver.
 *
 *  Uses MPI_COMM_WORLD so MUMPS distributes the factorization across all
 *  MPI ranks.  The matrix data is provided in centralized assembled mode
 *  (icntl[17]=0): rank 0 holds all irn/jcn/a entries and receives the
 *  solution; other ranks participate in computation only.
 *
 *  Requires TAHOE_MUMPS=ON and TAHOE_MPI=ON.
 *
 *  Usage in XML input file (MPI run only):
 *    <MUMPS_MPI_matrix message_level="silent" always_symmetric="false"/>
 */
class MUMPSMatrixT_mpi : public MSRMatrixT
{
public:

    MUMPSMatrixT_mpi(ostream& out, int check_code, bool symmetric,
        int message_level, const CommunicatorT& comm);

    MUMPSMatrixT_mpi(const MUMPSMatrixT_mpi& source);

    virtual ~MUMPSMatrixT_mpi(void);

    virtual bool SolvePreservesData(void) const { return true; }

    virtual void Clear(void);

    virtual GlobalT::SystemTypeT MatrixType(void) const { return GlobalT::kNonSymmetric; }

    virtual GlobalMatrixT* Clone(void) const;

protected:

    virtual void Factorize(void);
    virtual void BackSubstitute(dArrayT& result);

    virtual void PrintAllPivots(void) const {}
    virtual void PrintZeroPivots(void) const {}
    virtual void PrintLHS(bool force = false) const {}

private:

    void Initialize(void);
    void Finalize(void);

    DMUMPS_STRUC_C fId;
    int  fMessageLevel;
    bool fSymmetric;
    bool fIsInitialized;
    bool fIsFactorized;

    iArrayT fRowIdx;
    iArrayT fColIdx;
    dArrayT fValues;
};

} // namespace Tahoe

#endif /* __MUMPS__ && __TAHOE_MPI__ */
#endif /* _MUMPS_MATRIX_T_MPI_H_ */
