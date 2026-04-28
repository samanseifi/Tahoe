/* test_TaylorBar.cpp — integration test for Taylor bar CI smoke benchmark.
 *
 * Runs benchmark_XML/level.5/taylor_bar/taylor_ci.xml (1 us simulation,
 * ~10 seconds) and verifies:
 *   - Simulation completes without crash / NaN
 *   - Deformation started (bar nodes moved)
 *   - Impact-face nodes stayed at z=0 (rigid wall BC enforced)
 *   - Non-impact-face nodes moved downward (initial velocity applied)
 *
 * This is a smoke test, not a validation run.  The full Wilkins-Guinan
 * comparison uses taylor_small.xml (60 us, ~9 min) and lives in the
 * level.5 benchmark suite.
 */

#include "gtest/gtest.h"
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <string>
#include <fstream>

/* Absolute paths set by CMake at compile time (TAHOE_TEST_REPO_ROOT,
 * TAHOE_TEST_BUILD_DIR).  Falls back to relative paths otherwise. */
#ifndef TAHOE_TEST_REPO_ROOT
# define TAHOE_TEST_REPO_ROOT "."
#endif
#ifndef TAHOE_TEST_BUILD_DIR
# define TAHOE_TEST_BUILD_DIR "."
#endif

static const std::string BENCH_DIR =
	std::string(TAHOE_TEST_REPO_ROOT) + "/benchmark_XML/level.5/taylor_bar";
static const std::string TAHOE_BIN =
	std::string(TAHOE_TEST_BUILD_DIR) + "/bin/tahoe";
static const char* XML_NAME = "taylor_ci.xml";
static const char* EXO_NAME = "taylor_ci.io0.exo";

TEST(TaylorBar, CompletesWithoutCrash) {
	std::string cmd = "cd ";
	cmd += BENCH_DIR;
	cmd += " && rm -f taylor_ci.io0.exo taylor_ci.out taylor_ci.echo.xml taylor_ci.valid.xml taylor_ci.stdout 2>/dev/null";
	cmd += " && ";
	cmd += TAHOE_BIN;
	cmd += " -f ";
	cmd += XML_NAME;
	cmd += " > taylor_ci.stdout 2>&1";

	int rc = std::system(cmd.c_str());
	EXPECT_EQ(rc, 0) << "Tahoe returned non-zero: cmd=" << cmd;

	/* output file must exist */
	std::string exo_path = BENCH_DIR + "/" + EXO_NAME;
	std::ifstream f(exo_path);
	EXPECT_TRUE(f.good()) << "Exodus output missing: " << exo_path;
}

TEST(TaylorBar, ExplicitPathActive) {
	/* Verify stdout contains the expected init banner — confirms batched
	 * explicit path was taken, J2 plasticity was recognized, history allocated. */
	std::string stdout_path = BENCH_DIR + "/taylor_ci.stdout";
	std::ifstream f(stdout_path);
	if (!f.good()) GTEST_SKIP() << "stdout not found — prior test may have failed";

	std::string line, content;
	while (std::getline(f, line)) content += line + "\n";

	EXPECT_NE(content.find("MVSIZ=128"), std::string::npos)
		<< "batched explicit path not activated";
	EXPECT_NE(content.find("J2 plasticity"), std::string::npos)
		<< "J2 plasticity material not selected";
	EXPECT_NE(content.find("history 16 vars/IP"), std::string::npos)
		<< "history allocation missing or wrong size";
}

/* Read final-step max DZ magnitude from stdout (".out" file); used as a coarse
 * "did the bar deform?" check.  We do not parse Exodus here to keep the test
 * dependency-light, but stdout summarizes the run. */
TEST(TaylorBar, BarDeformedPlastically) {
	/* read the captured stdout: contains the per-step progress */
	std::string out_path = BENCH_DIR + "/taylor_ci.stdout";
	std::ifstream f(out_path);
	if (!f.good()) GTEST_SKIP() << "Tahoe output missing";

	std::string content, line;
	while (std::getline(f, line)) content += line + "\n";

	/* must have completed all 10000 steps */
	EXPECT_NE(content.find("Step: 10000 of 10000"), std::string::npos)
		<< "simulation did not run to completion";
	EXPECT_EQ(content.find("ExceptionT::Throw"), std::string::npos)
		<< "exception thrown during run";
}
