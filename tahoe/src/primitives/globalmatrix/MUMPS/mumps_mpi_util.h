/* mumps_mpi_util.h — MPI init/finalize helpers for the MUMPS sequential wrapper.
 *
 * Isolates MPI header inclusion into a .c file to avoid conflicts between
 * OpenMPI's C++ bindings and Tahoe's CommunicatorT stub macros (MPI_Comm=long).
 */
#ifndef MUMPS_MPI_UTIL_H
#define MUMPS_MPI_UTIL_H

#ifdef __cplusplus
extern "C" {
#endif

/** Initialize MPI if it has not already been initialized.
 *  Returns 1 if this call performed MPI_Init, 0 if MPI was already running. */
int tahoe_mumps_mpi_init(void);

/** Finalize MPI only if we_inited is non-zero (i.e. tahoe_mumps_mpi_init
 *  returned 1) and MPI has not already been finalized. */
void tahoe_mumps_mpi_finalize(int we_inited);

/** Return the Fortran handle for MPI_COMM_SELF (single-rank communicator).
 *  Use this instead of the -987654 sentinel: that value is only valid for
 *  libmumps-seq; the parallel libmumps-dev requires a real communicator. */
int tahoe_mumps_comm_self(void);

#ifdef __cplusplus
}
#endif

#endif /* MUMPS_MPI_UTIL_H */
