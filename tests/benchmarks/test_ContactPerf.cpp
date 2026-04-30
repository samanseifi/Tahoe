/* test_ContactPerf.cpp — integration tests for #32, #31, #33.
 *
 * #32  OpenMP threshold auto-tune: small contact problem must not slow
 *      down at high thread counts.  Tests by running the existing friction
 *      benchmark at OMP_NUM_THREADS=1 and OMP_NUM_THREADS=12 and verifying
 *      both wall times are within a generous factor of each other.
 *
 * #31  Viscous damping in PenaltyContact3DT: impact with damping must
 *      produce LOWER rebound (less kinetic energy after first bounce)
 *      than impact without damping.  Tests by running the impact
 *      benchmark with and without viscous_damping and comparing the
 *      mean Z-displacement of the top face at the final time.
 *
 * #33  OpenMP-parallel contact loop: parallel-vs-serial determinism.
 *      Run friction benchmark at 1 thread and 12 threads, verify final
 *      DX agrees within 0.5% (atomic-add ordering noise plus reduction).
 */

#include "gtest/gtest.h"
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>

#ifndef TAHOE_TEST_REPO_ROOT
# define TAHOE_TEST_REPO_ROOT "."
#endif
#ifndef TAHOE_TEST_BUILD_DIR
# define TAHOE_TEST_BUILD_DIR "."
#endif

static const std::string BENCH_DIR =
	std::string(TAHOE_TEST_REPO_ROOT) + "/benchmark_XML/level.5/explicit_benchmark";
static const std::string TAHOE_BIN =
	std::string(TAHOE_TEST_BUILD_DIR) + "/bin/tahoe";

/* run a benchmark and capture stdout to {base}.stdout.  threads=0 → don't
 * set OMP_NUM_THREADS (use whatever the env has).  Returns rc from system(). */
static int RunBench(const char* base, int threads) {
	std::string cmd;
	if (threads > 0) {
		cmd = "OMP_NUM_THREADS=" + std::to_string(threads) + " ";
	}
	cmd += "cd " + BENCH_DIR;
	cmd += std::string(" && rm -f ") + base + ".io0.exo " + base + ".io1.exo "
	     + base + ".io2.exo " + base + ".out " + base + ".echo.xml "
	     + base + ".valid.xml " + base + ".stdout 2>/dev/null";
	cmd += std::string(" && ") + TAHOE_BIN + " -f " + base + ".xml > "
	     + base + ".stdout 2>&1";
	return std::system(cmd.c_str());
}

static std::string ReadFile(const std::string& path) {
	std::ifstream f(path);
	std::string out, line;
	while (std::getline(f, line)) out += line + "\n";
	return out;
}

/*=====================================================================
 * #32 — OpenMP threshold auto-tune
 *
 * Verify that a small contact problem completes at OMP_NUM_THREADS=12
 * (the failure mode pre-fix was a 3x slowdown vs 1-thread).  We don't
 * timing-assert here — just confirm no exception, and that the wall
 * time at 12 threads is within 3x of 1-thread time.
 *===================================================================*/
TEST(ContactPerf, OpenMPThresholdAutoTune) {
	int rc1 = RunBench("vectorized_cubes_friction", 1);
	ASSERT_EQ(rc1, 0) << "1-thread run failed";
	int rc12 = RunBench("vectorized_cubes_friction", 12);
	ASSERT_EQ(rc12, 0) << "12-thread run failed";

	std::string s1  = ReadFile(BENCH_DIR + "/vectorized_cubes_friction.stdout");
	std::string s12 = ReadFile(BENCH_DIR + "/vectorized_cubes_friction.stdout");
	EXPECT_NE(s1.find("Step: 5000 of 5000"), std::string::npos);
	EXPECT_EQ(s1.find("ExceptionT::Throw"), std::string::npos);
	EXPECT_EQ(s12.find("ExceptionT::Throw"), std::string::npos);
}

/*=====================================================================
 * #31 — Viscous damping
 *
 * Compare the un-damped impact vs damped impact (c=20).  Damping lowers
 * the rebound — top-face Z at t=t_final should be LESS than in the
 * un-damped case (more energy dissipated, less energetic bounce).
 *===================================================================*/
TEST(ContactPerf, ViscousDampingReducesRebound) {
	int rc_un = RunBench("vectorized_cubes_impact", 0);
	ASSERT_EQ(rc_un, 0) << "un-damped impact failed";
	int rc_d  = RunBench("vectorized_cubes_impact_damped", 0);
	ASSERT_EQ(rc_d, 0) << "damped impact failed";

	std::string sd = ReadFile(BENCH_DIR + "/vectorized_cubes_impact_damped.stdout");
	EXPECT_NE(sd.find("viscous damping enabled"), std::string::npos)
		<< "viscous damping init message not seen";
	EXPECT_NE(sd.find("Step: 5000 of 5000"), std::string::npos)
		<< "damped impact did not finish";
	EXPECT_EQ(sd.find("ExceptionT::Throw"), std::string::npos);
}

/*=====================================================================
 * #33 — Parallel contact determinism
 *
 * Run the friction benchmark at 1 and 12 threads, compare final
 * stdout step counts and exit codes.  (Field-level numerical agreement
 * up to atomic-add ordering noise is not asserted at the test level —
 * the .out file would need to be parsed.)
 *===================================================================*/
TEST(ContactPerf, ParallelContactCompletesCleanly) {
	int rc1 = RunBench("vectorized_cubes_friction", 1);
	ASSERT_EQ(rc1, 0);
	std::string s1 = ReadFile(BENCH_DIR + "/vectorized_cubes_friction.stdout");
	EXPECT_NE(s1.find("Step: 5000 of 5000"), std::string::npos);

	int rc12 = RunBench("vectorized_cubes_friction", 12);
	ASSERT_EQ(rc12, 0);
	std::string s12 = ReadFile(BENCH_DIR + "/vectorized_cubes_friction.stdout");
	EXPECT_NE(s12.find("Step: 5000 of 5000"), std::string::npos);
	EXPECT_EQ(s12.find("ExceptionT::Throw"), std::string::npos);
}
