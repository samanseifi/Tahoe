/* MUMPSMatrixT_mpi.cpp — MUMPS MPI distributed sparse direct solver for Tahoe */
#include "MUMPSMatrixT_mpi.h"

#if defined(__MUMPS__) && defined(__TAHOE_MPI__)

#include "CommunicatorT.h"
#include "ExceptionT.h"
#include <cstring>

using namespace Tahoe;

/* -----------------------------------------------------------------------
 * Constructor / destructor
 * ----------------------------------------------------------------------- */

MUMPSMatrixT_mpi::MUMPSMatrixT_mpi(ostream& out, int check_code, bool symmetric,
    int message_level, const CommunicatorT& comm)
  : MSRMatrixT(out, check_code, symmetric, comm),
    fMessageLevel(message_level),
    fSymmetric(symmetric),
    fIsInitialized(false),
    fIsFactorized(false)
{
    memset(&fId, 0, sizeof(DMUMPS_STRUC_C));
}

MUMPSMatrixT_mpi::MUMPSMatrixT_mpi(const MUMPSMatrixT_mpi& source)
  : MSRMatrixT(source),
    fMessageLevel(source.fMessageLevel),
    fSymmetric(source.fSymmetric),
    fIsInitialized(false),
    fIsFactorized(false)
{
    ExceptionT::GeneralFail("MUMPSMatrixT_mpi::MUMPSMatrixT_mpi(copy)", "not implemented");
}

MUMPSMatrixT_mpi::~MUMPSMatrixT_mpi(void)
{
    if (fIsInitialized) Finalize();
}

/* -----------------------------------------------------------------------
 * Public
 * ----------------------------------------------------------------------- */

void MUMPSMatrixT_mpi::Clear(void)
{
    MSRMatrixT::Clear();
    if (!fIsInitialized) Initialize();
    fIsFactorized = false;
}

GlobalMatrixT* MUMPSMatrixT_mpi::Clone(void) const
{
    return new MUMPSMatrixT_mpi(*this);
}

/* -----------------------------------------------------------------------
 * Protected
 * ----------------------------------------------------------------------- */

void MUMPSMatrixT_mpi::Factorize(void)
{
    if (fIsFactorized) return;

    /* Rank 0 holds the full matrix; other ranks pass null pointers.
     * MUMPS centralized assembled input (icntl[17]=0) requires only rank 0
     * to provide irn/jcn/a. */
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        iArrayT r, c;
        dArrayT v;
        GenerateRCV(r, c, v, 1.0e-15);

        const int nnz = r.Length();
        fRowIdx.Dimension(nnz);
        fColIdx.Dimension(nnz);
        fValues.Dimension(nnz);
        for (int i = 0; i < nnz; i++) {
            fRowIdx[i] = r[i] + 1;
            fColIdx[i] = c[i] + 1;
            fValues[i] = v[i];
        }

        fId.n   = fTotNumEQ;
        fId.nz  = nnz;
        fId.irn = fRowIdx.Pointer();
        fId.jcn = fColIdx.Pointer();
        fId.a   = fValues.Pointer();
    }

    /* Job 4: Analysis + Numerical factorization (combined) */
    fId.job = 4;
    dmumps_c(&fId);
    if (fId.infog[0] != 0)
        ExceptionT::BadJacobianDet("MUMPSMatrixT_mpi::Factorize",
            "MUMPS analysis+factorization (job 4) failed: infog[0]=%d infog[1]=%d "
            "(negative value means singular or near-singular matrix)",
            fId.infog[0], fId.infog[1]);

    fIsFactorized = true;
}

void MUMPSMatrixT_mpi::BackSubstitute(dArrayT& result)
{
    if (!fIsFactorized)
        ExceptionT::GeneralFail("MUMPSMatrixT_mpi::BackSubstitute",
            "matrix has not been factorized");

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        fId.rhs  = result.Pointer();
        fId.nrhs = 1;
        fId.lrhs = fTotNumEQ;
    }

    fId.job = 3;
    dmumps_c(&fId);
    if (fId.infog[0] != 0)
        ExceptionT::BadJacobianDet("MUMPSMatrixT_mpi::BackSubstitute",
            "MUMPS solve (job 3) failed: infog[0]=%d infog[1]=%d",
            fId.infog[0], fId.infog[1]);
}

/* -----------------------------------------------------------------------
 * Private
 * ----------------------------------------------------------------------- */

void MUMPSMatrixT_mpi::Initialize(void)
{
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    memset(&fId, 0, sizeof(DMUMPS_STRUC_C));

    fId.job          = -1;
    fId.par          =  1;
    fId.sym          = fSymmetric ? 2 : 0;
    fId.comm_fortran = (MUMPS_INT) MPI_Comm_c2f(MPI_COMM_WORLD);
    dmumps_c(&fId);

    if (fId.infog[0] != 0)
        ExceptionT::GeneralFail("MUMPSMatrixT_mpi::Initialize",
            "MUMPS init (job -1) failed: infog[0]=%d infog[1]=%d",
            fId.infog[0], fId.infog[1]);

    /* Only rank 0 prints — suppress output on all worker ranks to avoid
     * N-fold duplication of diagnostic messages with N MPI processes. */
    const bool print = (rank == 0) && (fMessageLevel > 0);
    fId.icntl[0] = print ? 6 : -1;
    fId.icntl[1] = (rank == 0 && fMessageLevel > 1) ? 6 : -1;
    fId.icntl[2] = print ? 6 : -1;
    fId.icntl[3] = print ? fMessageLevel : 0;

    /* Centralized assembled input — rank 0 provides all data */
    fId.icntl[17] = 0;

    /* AMD fill-reducing ordering */
    fId.icntl[6] = 0;

    /* Allow 200% of estimated workspace (default 20% can be too tight) */
    fId.icntl[13] = 100;

    fIsInitialized = true;
}

void MUMPSMatrixT_mpi::Finalize(void)
{
    fId.job = -2;
    dmumps_c(&fId);
    fIsInitialized = false;
    fIsFactorized  = false;
}

#endif /* __MUMPS__ && __TAHOE_MPI__ */
