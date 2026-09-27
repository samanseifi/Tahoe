/* BLASRuntimeT.h — report and configure the BLAS used by MUMPS at run time (#80) */
#ifndef _BLAS_RUNTIME_T_H_
#define _BLAS_RUNTIME_T_H_

#include <string>

namespace Tahoe {

/** Run-time view of the BLAS library the process actually uses.
 *
 *  MUMPS calls dgemm_ and friends through whatever BLAS the dynamic loader
 *  resolved, which depends on how tahoe was linked (see #80 in
 *  tahoe/CMakeLists.txt) and on the system's libblas.so.3.  This class asks
 *  the loader instead of trusting the build configuration, so the report is
 *  correct for any BLAS (OpenBLAS, MKL, BLIS, Accelerate, reference).
 *
 *  All lookups use dlsym/dladdr, so no BLAS header or vendor library is
 *  needed at compile time. */
class BLASRuntimeT
{
public:

	/** one-line description: vendor, version or configuration, the library
	 *  file dgemm_ resolved to, and the number of BLAS threads */
	static std::string Describe(void);

	/** use one BLAS thread unless the environment sets a count
	 *  (OPENBLAS_NUM_THREADS, GOTO_NUM_THREADS, MKL_NUM_THREADS or
	 *  OMP_NUM_THREADS).  A single thread was fastest for the sparse
	 *  factorizations measured in #80 and avoids oversubscription with
	 *  Tahoe's own OpenMP threads.  Applied once per process. */
	static void SetDefaultThreads(void);
};

} /* namespace Tahoe */

#endif /* _BLAS_RUNTIME_T_H_ */
