/* test_J2SimoTangent.cpp — consistent-tangent check for Simo_J2 (#78).
 *
 * Runs benchmark_XML/level.0/matrix_check/j2_simo_tangent.xml in a scratch
 * directory under the build tree.  The deck solves one distorted hex8 in
 * tension plus shear with J2 plasticity and check_code="check_LHS", so at
 * every Newton iteration the solver writes the analytic stiffness
 * (FullMatrixT.LHS.<even>) and a forward-difference stiffness assembled from
 * the residual (FullMatrixT.LHS.<odd>).
 *
 * With a consistent tangent the two agree to the truncation error of the
 * difference, about 3e-8 relative.  Before #78 the beta_2 term of the
 * Simo & Hughes (1998) Box 9.2 tangent lacked the factor ||s_trial|| and the
 * plastic iterations differed by about 1e-3.
 */

#include "gtest/gtest.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <utility>

#ifndef TAHOE_TEST_REPO_ROOT
# define TAHOE_TEST_REPO_ROOT "."
#endif
#ifndef TAHOE_TEST_BUILD_DIR
# define TAHOE_TEST_BUILD_DIR "."
#endif

namespace {

const std::string DECK_DIR =
	std::string(TAHOE_TEST_REPO_ROOT) + "/benchmark_XML/level.0/matrix_check";
const std::string WORK_DIR =
	std::string(TAHOE_TEST_BUILD_DIR) + "/test_J2SimoTangent";
const std::string TAHOE_BIN =
	std::string(TAHOE_TEST_BUILD_DIR) + "/bin/tahoe";

typedef std::map<std::pair<int,int>, double> SparseT;

/* read a matrix written by FullMatrixT::PrintLHS (row col value) */
bool ReadRCV(const std::string& file, SparseT& K)
{
	std::ifstream in(file.c_str());
	if (!in.good()) return false;
	K.clear();
	int r, c;
	double v;
	while (in >> r >> c >> v) K[std::make_pair(r, c)] = v;
	return true;
}

/* ||A - B||_F / ||B||_F */
double RelativeDifference(const SparseT& A, const SparseT& B)
{
	SparseT D = B;
	for (SparseT::const_iterator i = A.begin(); i != A.end(); ++i)
		D[i->first] -= i->second;
	double num = 0.0, den = 0.0;
	for (SparseT::const_iterator i = D.begin(); i != D.end(); ++i)
		num += i->second*i->second;
	for (SparseT::const_iterator i = B.begin(); i != B.end(); ++i)
		den += i->second*i->second;
	return std::sqrt(num/den);
}

} /* namespace */

TEST(J2SimoTangent, AnalyticMatchesFiniteDifference)
{
	std::string cmd = "rm -rf " + WORK_DIR + " && mkdir -p " + WORK_DIR
		+ " && cp " + DECK_DIR + "/j2_simo_tangent.xml " + DECK_DIR + "/j2_simo_tangent.geom " + WORK_DIR
		+ " && cd " + WORK_DIR
		+ " && " + TAHOE_BIN + " -f j2_simo_tangent.xml > j2_simo_tangent.stdout 2>&1";
	int rc = std::system(cmd.c_str());
	ASSERT_EQ(rc, 0) << "tahoe returned non-zero exit; cmd=" << cmd;

	/* the run must finish without an exception */
	std::ifstream f((WORK_DIR + "/j2_simo_tangent.stdout").c_str());
	ASSERT_TRUE(f.good());
	std::stringstream log;
	log << f.rdbuf();
	ASSERT_NE(log.str().find("Step: 1 of 1"), std::string::npos) << "step not reached";
	ASSERT_EQ(log.str().find("exit on exception"), std::string::npos) << "run aborted";

	/* compare every analytic / finite-difference pair */
	int pairs = 0;
	for (int k = 0; ; k += 2)
	{
		std::ostringstream a, b;
		a << WORK_DIR << "/FullMatrixT.LHS." << k;
		b << WORK_DIR << "/FullMatrixT.LHS." << k + 1;
		SparseT K_an, K_fd;
		if (!ReadRCV(a.str(), K_an) || !ReadRCV(b.str(), K_fd)) break;
		ASSERT_FALSE(K_fd.empty());

		double rel = RelativeDifference(K_an, K_fd);
		EXPECT_LT(rel, 1.0e-6) << "iteration " << k/2
			<< ": analytic and finite-difference stiffness differ by " << rel;
		pairs++;
	}

	/* the elastic predictor plus at least three plastic iterations */
	EXPECT_GE(pairs, 4) << "too few Newton iterations were checked";
}
