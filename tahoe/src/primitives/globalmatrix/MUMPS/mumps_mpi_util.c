/* mumps_mpi_util.c — MPI init/finalize helpers for the MUMPS sequential wrapper.
 *
 * This is a plain C file (not C++).  Including mpi.h in a C file avoids the
 * OpenMPI C++ binding headers that conflict with Tahoe's CommunicatorT stubs
 * in non-MPI builds.
 *
 * No Tahoe headers are included here — only the system MPI header.
 */
#ifdef __MUMPS__

/* Suppress C++ MPI binding inclusion (harmless in C, but good practice) */
#define OMPI_SKIP_MPICXX  1
#define MPICH_SKIP_MPICXX 1

#include "mpi.h"
#include "mumps_mpi_util.h"

int tahoe_mumps_mpi_init(void)
{
    int inited = 0;
    MPI_Initialized(&inited);
    if (!inited) {
        MPI_Init(NULL, NULL);
        return 1;
    }
    return 0;
}

void tahoe_mumps_mpi_finalize(int we_inited)
{
    /* Do NOT call MPI_Finalize here.  MPI can only be initialized and
     * finalized ONCE per process lifetime.  If multiple MUMPS objects are
     * created (e.g. initial-condition solve + main solve), the first one's
     * destructor would finalize MPI and the second one would then call MPI
     * functions on a dead environment.  Instead, leave MPI alive — the OS
     * cleans up at process exit, which is acceptable for a serial solver. */
    (void) we_inited;
}

int tahoe_mumps_comm_self(void)
{
    return (int) MPI_Comm_c2f(MPI_COMM_SELF);
}

#endif /* __MUMPS__ */
