/* MUMPSMatrixT.cpp — MUMPS sequential sparse direct solver wrapper for Tahoe */
#include "MUMPSMatrixT.h"

#ifdef __MUMPS__

#include "ExceptionT.h"
#include "mumps_mpi_util.h"
#include "BLASRuntimeT.h"
#include <cstring>

using namespace Tahoe;

/* -----------------------------------------------------------------------
 * Constructor / destructor
 * ----------------------------------------------------------------------- */

const char* MUMPSMatrixT::OrderingName(int ordering)
{
    switch (ordering) {
        case kAMD:    return "AMD";
        case kAMF:    return "AMF";
        case kSCOTCH: return "SCOTCH";
        case kPORD:   return "PORD";
        case kMETIS:  return "METIS";
        case kQAMD:   return "QAMD";
        case kAuto:   return "automatic";
        default:      return "unknown";
    }
}

int MUMPSMatrixT::EffectiveOrdering(int ordering, int num_equations)
{
    return (ordering == kPORD && num_equations < kPORDMinEquations) ? kAMD : ordering;
}

MUMPSMatrixT::MUMPSMatrixT(ostream& out, int check_code, bool symmetric,
    int message_level, int ordering, const CommunicatorT& comm)
  : MSRMatrixT(out, check_code, symmetric, comm),
    fMessageLevel(message_level),
    fOrdering(ordering),
    fSymmetric(symmetric),
    fIsInitialized(false),
    fIsFactorized(false),
    fInitedMPI(false)
{
    memset(&fId, 0, sizeof(DMUMPS_STRUC_C));
}

MUMPSMatrixT::MUMPSMatrixT(const MUMPSMatrixT& source)
  : MSRMatrixT(source),
    fMessageLevel(source.fMessageLevel),
    fOrdering(source.fOrdering),
    fSymmetric(source.fSymmetric),
    fIsInitialized(false),
    fIsFactorized(false),
    fInitedMPI(false)
{
    ExceptionT::GeneralFail("MUMPSMatrixT::MUMPSMatrixT(copy)", "not implemented");
}

MUMPSMatrixT::~MUMPSMatrixT(void)
{
    if (fIsInitialized) Finalize();
}

/* -----------------------------------------------------------------------
 * Public
 * ----------------------------------------------------------------------- */

void MUMPSMatrixT::Clear(void)
{
    /* inherited — zeroes the MSR value array */
    MSRMatrixT::Clear();

    /* checks — must be serial (same as SPOOLESMatrixT) */
    if (fTotNumEQ != fLocNumEQ)
        ExceptionT::GeneralFail("MUMPSMatrixT::Clear",
            "total equations (%d) != local equations (%d)",
            fTotNumEQ, fLocNumEQ);
    if (fStartEQ != 1)
        ExceptionT::GeneralFail("MUMPSMatrixT::Clear",
            "expecting first equation number to be 1 not %d", fStartEQ);

    /* Initialize MUMPS once (MPI + job=-1); subsequent steps reuse the same
     * MUMPS instance — Factorize() reruns jobs 1+2 with the new matrix. */
    if (!fIsInitialized) Initialize();

    fIsFactorized = false;
}

GlobalMatrixT* MUMPSMatrixT::Clone(void) const
{
    return new MUMPSMatrixT(*this);
}

/* -----------------------------------------------------------------------
 * Protected
 * ----------------------------------------------------------------------- */

void MUMPSMatrixT::Factorize(void)
{
    if (fIsFactorized) return;

    /* --- Extract COO triplets from MSR storage ---
     * GenerateRCV returns 0-based indices; MUMPS needs 1-based. */
    iArrayT r, c;
    dArrayT v;
    GenerateRCV(r, c, v, 1.0e-15);

    const int nnz = r.Length();
    fRowIdx.Dimension(nnz);
    fColIdx.Dimension(nnz);
    fValues.Dimension(nnz);
    for (int i = 0; i < nnz; i++) {
        fRowIdx[i] = r[i] + 1;   /* 0→1 based */
        fColIdx[i] = c[i] + 1;
        fValues[i] = v[i];
    }

    /* --- Pass matrix data to MUMPS --- */
    fId.n   = fTotNumEQ;
    fId.nz  = nnz;
    fId.irn = fRowIdx.Pointer();
    fId.jcn = fColIdx.Pointer();
    fId.a   = fValues.Pointer();

    /* --- Analysis + numerical factorization (job 4) ---
     * The analysis is redone every time.  Reusing it (job 2) while the
     * pattern is unchanged saved about 5 % on the Brinell deck but is
     * unsafe: the analysis computes a column permutation and scaling from
     * the matrix values (ICNTL(6), ICNTL(8)), and with the values of a
     * later Newton iterate MUMPS 5.6 declared a regular matrix singular
     * (#80, level.1 WLC deck). */
    if (fId.icntl[6] != kAuto)
        fId.icntl[6] = EffectiveOrdering(fOrdering, fTotNumEQ);
    bool ok = RunJob(4);

    /* an ordering this MUMPS build lacks: fall back to MUMPS's own choice */
    if (!ok && fId.icntl[6] != kAuto) {
        const int requested = fId.icntl[6];
        fId.icntl[6] = kAuto;
        ok = RunJob(4);
        if (ok)
            fOut << "\n MUMPS: " << OrderingName(requested)
                 << " ordering failed; using MUMPS automatic ordering\n";
        else
            fId.icntl[6] = requested;
    }

    if (!ok)
        ExceptionT::BadJacobianDet("MUMPSMatrixT::Factorize",
            "MUMPS analysis+factorization (job 4) failed: infog[0]=%d infog[1]=%d "
            "(negative value means singular or near-singular matrix)",
            fId.infog[0], fId.infog[1]);

    fIsFactorized = true;
}

void MUMPSMatrixT::BackSubstitute(dArrayT& result)
{
    if (!fIsFactorized)
        ExceptionT::GeneralFail("MUMPSMatrixT::BackSubstitute",
            "matrix has not been factorized");

    /* MUMPS overwrites fId.rhs in place with the solution */
    fId.rhs  = result.Pointer();
    fId.nrhs = 1;
    fId.lrhs = fTotNumEQ;

    fId.job = 3;   /* solve */
    dmumps_c(&fId);
    if (fId.infog[0] < 0)
        ExceptionT::BadJacobianDet("MUMPSMatrixT::BackSubstitute",
            "MUMPS solve (job 3) failed: infog[0]=%d infog[1]=%d",
            fId.infog[0], fId.infog[1]);
}

/* -----------------------------------------------------------------------
 * Private
 * ----------------------------------------------------------------------- */

void MUMPSMatrixT::Initialize(void)
{
    /* The parallel libmumps requires MPI to be initialized even when using
     * the sequential dummy communicator (comm_fortran=-987654).  If MPI is
     * not yet running (e.g. serial build without -DTAHOE_MPI=ON), we call
     * the C helper which performs MPI_Init without pulling in OpenMPI's
     * C++ binding headers that conflict with CommunicatorT's stub macros. */
    fInitedMPI = tahoe_mumps_mpi_init();

    memset(&fId, 0, sizeof(DMUMPS_STRUC_C));

    fId.job          = -1;  /* initialize */
    fId.par          =  1;  /* host participates in factorization */
    fId.sym          = fSymmetric ? 2 : 0;
    /* MPI_COMM_SELF (single rank) — compatible with parallel libmumps-dev.
     * The -987654 sentinel is only valid for the sequential libmumps-seq. */
    fId.comm_fortran = (MUMPS_INT) tahoe_mumps_comm_self();
    dmumps_c(&fId);

    if (fId.infog[0] < 0)
        ExceptionT::GeneralFail("MUMPSMatrixT::Initialize",
            "MUMPS init (job -1) failed: infog[0]=%d infog[1]=%d",
            fId.infog[0], fId.infog[1]);

    /* Output verbosity:
     *   icntl[0..2] = output unit for errors/diagnostics/global info (-1 = suppress)
     *   icntl[3]    = print level (0=nothing, 1=errors, 2=diagnostics, 3=stats, 4=info) */
    fId.icntl[0] = (fMessageLevel > 0) ? 6 : -1;
    fId.icntl[1] = (fMessageLevel > 1) ? 6 : -1;
    fId.icntl[2] = (fMessageLevel > 0) ? 6 : -1;
    fId.icntl[3] = (fMessageLevel > 0) ? fMessageLevel : 0;

    /* icntl[17]=0 : centralized assembled input (host provides irn/jcn/a) */
    fId.icntl[17] = 0;

    /* icntl[6] : fill-reducing ordering, PORD by default (#80) */
    fId.icntl[6] = EffectiveOrdering(fOrdering, fTotNumEQ);

    /* icntl[13]=100 : allow 200% of estimated workspace (default 20% often too tight) */
    fId.icntl[13] = 100;

    fIsInitialized = true;

    /* BLAS threads and a record of the BLAS actually in use */
    BLASRuntimeT::SetDefaultThreads();
    fOut << "\n MUMPS sparse direct solver:\n"
         << "    ordering . . . . . . . . . . . . . . . . . . . = " << OrderingName(fId.icntl[6]);
    if (fId.icntl[6] != fOrdering)
        fOut << " (" << OrderingName(fOrdering) << " is used from "
             << int(kPORDMinEquations) << " equations up)";
    fOut << '\n'
         << "    BLAS . . . . . . . . . . . . . . . . . . . . . = " << BLASRuntimeT::Describe() << '\n';
}

bool MUMPSMatrixT::RunJob(int job)
{
    fId.job = job;
    dmumps_c(&fId);
    return fId.infog[0] >= 0;   /* > 0 are warnings */
}

void MUMPSMatrixT::Finalize(void)
{
    fId.job = -2;   /* finalize — free MUMPS internal memory */
    dmumps_c(&fId);
    fIsInitialized = false;
    fIsFactorized  = false;

    /* Finalize MPI only if we were the ones who initialized it */
    tahoe_mumps_mpi_finalize(fInitedMPI);
    fInitedMPI = false;
}

#endif /* __MUMPS__ */
