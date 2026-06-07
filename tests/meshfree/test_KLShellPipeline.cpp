/* test_KLShellPipeline.cpp — in-pipeline regression for the meshfree_kl_shell ELEMENT.
 *
 * Runs the Scordelis-Lo roof deck through the actual tahoe binary (XML -> FEManagerT ->
 * SolverT -> RKShellT::LHS/RHSDriver -> framework sparse solve) and checks the free-edge
 * deflection. This tests the real element code path end-to-end (not just the standalone
 * kernels), keeping TDD on the Tahoe-native element. Reference deflection 0.3006; the 15x15
 * meshfree solution is ~0.275 (converges toward 0.30 under refinement). */

#include "gtest/gtest.h"

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>

#ifndef TAHOE_TEST_REPO_ROOT
#  define TAHOE_TEST_REPO_ROOT "."
#endif
#ifndef TAHOE_TEST_BUILD_DIR
#  define TAHOE_TEST_BUILD_DIR "."
#endif

TEST(KLShellPipeline, ScordelisLoThroughTahoe)
{
	const std::string dir = std::string(TAHOE_TEST_REPO_ROOT) + "/applications/meshfree_kl_shell";
	const std::string bin = std::string(TAHOE_TEST_BUILD_DIR) + "/bin/tahoe";

	std::string cmd = "cd " + dir
		+ " && rm -f scordelis_lo.stdout scordelis_lo.echo.xml scordelis_lo.out 2>/dev/null"
		+ " && " + bin + " -f scordelis_lo.xml > scordelis_lo.stdout 2>&1";
	int rc = std::system(cmd.c_str());
	ASSERT_EQ(rc, 0) << "tahoe run failed";

	std::ifstream f((dir + "/scordelis_lo.stdout").c_str());
	ASSERT_TRUE(f.good()) << "no stdout captured";
	std::stringstream ss; ss << f.rdbuf();
	std::string out = ss.str();

	/* no exceptions during the run */
	EXPECT_EQ(out.find("ExceptionT::Throw"), std::string::npos) << "tahoe threw:\n" << out;

	/* parse the element's deflection summary: "[RKShell] ... min(u_y)= <value> ..." */
	std::string key = "min(u_y)=";
	size_t p = out.find(key);
	ASSERT_NE(p, std::string::npos) << "no RKShell deflection summary:\n" << out;
	double uy = std::atof(out.c_str() + p + key.size());

	/* free-edge deflection is downward; 15x15 meshfree -> ~ -0.275 (ref -0.3006). Accept the
	 * converging band, comfortably bracketing both, but tight enough to catch a broken element. */
	double mag = (uy < 0.0) ? -uy : uy;
	EXPECT_GT(mag, 0.24) << "deflection too small (over-stiff): " << uy;
	EXPECT_LT(mag, 0.32) << "deflection too large (unstable/soft): " << uy;
}
