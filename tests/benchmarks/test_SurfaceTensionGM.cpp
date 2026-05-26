/* test_SurfaceTensionGM.cpp — integration test for issue #54
 * Gurtin-Murdoch strain-dependent surface elasticity.
 *
 * Runs benchmark_XML/level.4/surface_tension/gm_strip/ and verifies that
 *   1. yl_baseline_implicit.xml (E_s = 0) still runs cleanly      — backward
 *      compatibility with the legacy Young-Laplace path,
 *   2. gm_strip_implicit.xml    (E_s = 10) runs to completion,
 *   3. both produce Exodus output and reach the final step without
 *      exceptions or NaNs in stdout.
 *
 * Deep math (residual = -dPsi/du, K = d2Psi/du2 against JAX autodiff) is
 * covered by benchmark_XML/level.4/surface_tension/verify_GM_surface.py;
 * end-to-end implicit/explicit and GM/YL cross-checks are in
 * benchmark_XML/level.4/surface_tension/gm_strip/verify_gm_strip.py.
 *
 * This test exists so CI catches regressions in the C++ integration
 * (parameter reading, force assembly, tangent for E_s != 0).
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
    std::string(TAHOE_TEST_REPO_ROOT)
    + "/benchmark_XML/level.4/surface_tension/gm_strip";
static const std::string TAHOE_BIN =
    std::string(TAHOE_TEST_BUILD_DIR) + "/bin/tahoe";


static void run_and_check(const std::string& xml_stem, int expected_steps) {
    /* clean prior outputs, then run tahoe */
    const std::string xml = xml_stem + ".xml";
    const std::string stdout_file = xml_stem + ".stdout";

    std::string cmd = "cd " + BENCH_DIR
        + " && rm -f " + xml_stem + ".io0.exo " + xml_stem + ".out "
        + xml_stem + ".echo.xml " + xml_stem + ".valid.xml "
        + xml_stem + ".log " + stdout_file + " 2>/dev/null"
        + " && " + TAHOE_BIN + " -f " + xml + " > " + stdout_file + " 2>&1";

    int rc = std::system(cmd.c_str());
    ASSERT_EQ(rc, 0) << "tahoe exited non-zero on " << xml << "; cmd=" << cmd;

    /* Exodus output must exist */
    std::ifstream exo(BENCH_DIR + "/" + xml_stem + ".io0.exo");
    ASSERT_TRUE(exo.good()) << "Exodus output missing for " << xml_stem;

    /* read stdout for stack-health markers */
    std::ifstream f(BENCH_DIR + "/" + stdout_file);
    ASSERT_TRUE(f.good()) << "stdout missing for " << xml_stem;
    std::string content, line;
    while (std::getline(f, line)) content += line + "\n";

    /* no exception during run */
    EXPECT_EQ(content.find("ExceptionT::Throw"), std::string::npos)
        << "exception thrown during " << xml_stem << " run";
    EXPECT_EQ(content.find("nan"), std::string::npos)
        << "NaN appeared in " << xml_stem << " stdout";
    EXPECT_EQ(content.find("NaN"), std::string::npos)
        << "NaN appeared in " << xml_stem << " stdout";

    /* reached the final step */
    const std::string done_marker =
        "Step: " + std::to_string(expected_steps)
        + " of " + std::to_string(expected_steps);
    EXPECT_NE(content.find(done_marker), std::string::npos)
        << "did not reach final step for " << xml_stem
        << " (expected '" << done_marker << "')";
}


TEST(SurfaceTensionGM, YoungLaplaceBaselineStillRuns) {
    /* E_s = 0:  must reproduce the pre-#54 Young-Laplace path exactly. */
    run_and_check("yl_baseline_implicit", 2000);
}


TEST(SurfaceTensionGM, GMSurfaceElasticityRuns) {
    /* E_s = 10: exercises the new strain-dependent surface stress
     * sigma_s = gamma_0 + E_s*eps_eng plus its tangent contribution. */
    run_and_check("gm_strip_implicit", 2000);
}
