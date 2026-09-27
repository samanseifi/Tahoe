/* BLASRuntimeT.cpp — report and configure the BLAS used by MUMPS at run time (#80) */
#include "BLASRuntimeT.h"

#include <cstdlib>
#include <cstdio>
#include <dlfcn.h>

using namespace Tahoe;

namespace {

/* symbol lookup across everything the process has loaded */
void* Symbol(const char* name) { return dlsym(RTLD_DEFAULT, name); }

bool EnvSet(const char* name)
{
	const char* v = getenv(name);
	return v != NULL && v[0] != '\0';
}

/* the library file that provides dgemm_ */
std::string DgemmLibrary(void)
{
	void* f = Symbol("dgemm_");
	Dl_info info;
	if (f && dladdr(f, &info) && info.dli_fname) return info.dli_fname;
	return "unresolved";
}

int OpenBLASThreads(void)
{
	typedef int (*get_t)(void);
	get_t get = (get_t) Symbol("openblas_get_num_threads");
	return get ? get() : -1;
}

int MKLThreads(void)
{
	typedef int (*get_t)(void);
	get_t get = (get_t) Symbol("MKL_Get_Max_Threads");
	return get ? get() : -1;
}

} /* namespace */

std::string BLASRuntimeT::Describe(void)
{
	char buf[512];
	std::string lib = DgemmLibrary();

	typedef char* (*config_t)(void);
	config_t oconfig = (config_t) Symbol("openblas_get_config");
	if (oconfig)
	{
		snprintf(buf, sizeof(buf), "OpenBLAS (%s), %d thread(s), %s",
			oconfig(), OpenBLASThreads(), lib.c_str());
		return buf;
	}

	typedef void (*mklver_t)(char*, int);
	mklver_t mklver = (mklver_t) Symbol("MKL_Get_Version_String");
	if (mklver)
	{
		char ver[128];
		mklver(ver, sizeof(ver));
		ver[sizeof(ver) - 1] = '\0';
		snprintf(buf, sizeof(buf), "%s, %d thread(s), %s", ver, MKLThreads(), lib.c_str());
		return buf;
	}

	typedef const char* (*blisver_t)(void);
	blisver_t blisver = (blisver_t) Symbol("bli_info_get_version_str");
	if (blisver)
	{
		snprintf(buf, sizeof(buf), "BLIS %s, %s", blisver(), lib.c_str());
		return buf;
	}

	snprintf(buf, sizeof(buf), "%s (no optimized-BLAS entry points found: "
		"probably the reference BLAS; see TAHOE_BLAS in the build)", lib.c_str());
	return buf;
}

void BLASRuntimeT::SetDefaultThreads(void)
{
	static bool done = false;
	if (done) return;
	done = true;

	typedef void (*set_t)(int);

	if (!EnvSet("OPENBLAS_NUM_THREADS") && !EnvSet("GOTO_NUM_THREADS") && !EnvSet("OMP_NUM_THREADS"))
	{
		set_t set = (set_t) Symbol("openblas_set_num_threads");
		if (set) set(1);
	}

	if (!EnvSet("MKL_NUM_THREADS") && !EnvSet("OMP_NUM_THREADS"))
	{
		set_t set = (set_t) Symbol("MKL_Set_Num_Threads");
		if (set) set(1);
	}
}
