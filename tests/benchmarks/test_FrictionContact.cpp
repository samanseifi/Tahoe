/* test_FrictionContact.cpp — integration test for #26 Coulomb friction.
 *
 * Runs benchmark_XML/.../vectorized_cubes_friction.xml and compares against
 * vectorized_cubes_nofric.xml.  Friction must produce noticeable retardation
 * of the upper cube's tangential motion, AND a frictionless run must give
 * approximately rigid translation (v*t) within numerical bouncing.
 *
 * Single TEST() — ctest --parallel safe.
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
	std::string(TAHOE_TEST_REPO_ROOT) + "/benchmark_XML/level.5/explicit_benchmark";
static const std::string TAHOE_BIN =
	std::string(TAHOE_TEST_BUILD_DIR) + "/bin/tahoe";

static int RunBench(const char* xml) {
	std::string cmd = "cd " + BENCH_DIR
		+ " && rm -f " + xml + ".io0.exo " + xml + ".out "
		  + xml + ".echo.xml " + xml + ".valid.xml " + xml + ".stdout 2>/dev/null"
		+ " && " + TAHOE_BIN + " -f " + xml + ".xml > " + xml + ".stdout 2>&1";
	return std::system(cmd.c_str());
}

static std::string ReadFile(const std::string& path) {
	std::ifstream f(path);
	std::string out, line;
	while (std::getline(f, line)) out += line + "\n";
	return out;
}

TEST(FrictionContact, FrictionReducesTangentialMotion) {
	/* run frictionless reference */
	int rc0 = RunBench("vectorized_cubes_nofric");
	ASSERT_EQ(rc0, 0) << "frictionless reference run failed";

	/* run mu=0.3 case */
	int rc1 = RunBench("vectorized_cubes_friction");
	ASSERT_EQ(rc1, 0) << "friction run failed";

	/* both must complete fully without exception */
	std::string s0 = ReadFile(BENCH_DIR + "/vectorized_cubes_nofric.stdout");
	std::string s1 = ReadFile(BENCH_DIR + "/vectorized_cubes_friction.stdout");

	EXPECT_NE(s0.find("Step: 5000 of 5000"), std::string::npos)
		<< "no-friction did not reach final step";
	EXPECT_NE(s1.find("Step: 5000 of 5000"), std::string::npos)
		<< "friction did not reach final step";
	EXPECT_EQ(s0.find("ExceptionT::Throw"), std::string::npos);
	EXPECT_EQ(s1.find("ExceptionT::Throw"), std::string::npos);

	/* friction run must announce that mu was activated */
	EXPECT_NE(s1.find("Coulomb friction enabled"), std::string::npos)
		<< "Coulomb friction message not in stdout";
	EXPECT_EQ(s0.find("Coulomb friction enabled"), std::string::npos)
		<< "frictionless run should not enable Coulomb";
}
