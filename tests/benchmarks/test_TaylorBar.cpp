/* test_TaylorBar.cpp — integration test for Taylor bar CI smoke benchmark.
 *
 * Runs benchmark_XML/level.5/taylor_bar/taylor_ci.xml (1 us simulation,
 * ~30 seconds on CI) and verifies stack health:
 *   - Tahoe runs end-to-end without crash / NaN
 *   - Exodus output is produced
 *   - Batched explicit path was activated
 *   - J2 plasticity material was recognized and history allocated
 *   - All 10000 steps completed
 *   - No ExceptionT::Throw in stdout
 *
 * Note: this is a single TEST() so ctest --parallel can't reorder it
 * relative to subordinate assertions.  Splitting into multiple tests
 * caused the assertion-only tests to race ahead of the run test on CI.
 *
 * Full Wilkins-Guinan validation uses taylor_small.xml (60 us, ~10 min)
 * and lives in the level.5 benchmark suite, not in CI.
 */

#include "gtest/gtest.h"
#include <cstdlib>
#include <fstream>
#include <string>

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

TEST(TaylorBar, RunsAndDeformsPlastically) {
	/* run tahoe end-to-end on the CI smoke XML */
	std::string cmd = "cd " + BENCH_DIR
		+ " && rm -f taylor_ci.io0.exo taylor_ci.out taylor_ci.echo.xml"
		  " taylor_ci.valid.xml taylor_ci.stdout 2>/dev/null"
		+ " && " + TAHOE_BIN + " -f taylor_ci.xml > taylor_ci.stdout 2>&1";

	int rc = std::system(cmd.c_str());
	ASSERT_EQ(rc, 0) << "Tahoe returned non-zero exit; cmd=" << cmd;

	/* Exodus file must exist */
	std::ifstream exo(BENCH_DIR + "/taylor_ci.io0.exo");
	ASSERT_TRUE(exo.good()) << "Exodus output missing";

	/* read the captured stdout */
	std::ifstream f(BENCH_DIR + "/taylor_ci.stdout");
	ASSERT_TRUE(f.good()) << "stdout missing";
	std::string content, line;
	while (std::getline(f, line)) content += line + "\n";

	/* batched explicit path activated and J2 history allocated */
	EXPECT_NE(content.find("MVSIZ=128"), std::string::npos)
		<< "batched explicit path not activated";
	EXPECT_NE(content.find("J2 plasticity"), std::string::npos)
		<< "J2 plasticity not selected";
	EXPECT_NE(content.find("history 16 vars/IP"), std::string::npos)
		<< "history allocation missing or wrong size";

	/* simulation ran to completion without exception */
	EXPECT_NE(content.find("Step: 10000 of 10000"), std::string::npos)
		<< "did not reach final step";
	EXPECT_EQ(content.find("ExceptionT::Throw"), std::string::npos)
		<< "exception thrown during run";
}
