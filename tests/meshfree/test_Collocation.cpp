/* test_Collocation.cpp — unit tests for CollocationSolverT (direct nodal-collocation BC transform).
 *
 * Verifies that solving Phi_c d_c = ubar - Phi_free d_free yields constrained coefficients whose
 * RECONSTRUCTED physical value sum_J Phi_J(x_I) d_J equals the prescribed ubar_I. Covers the
 * interpolatory limit (Phi = delta -> d_c = ubar), a fully non-interpolatory constrained block,
 * and a case with free-node contributions in the support. */
#include "gtest/gtest.h"
#include "CollocationSolverT.h"
#include "RaggedArray2DT.h"
#include "iArrayT.h"
#include "dArrayT.h"

using namespace Tahoe;

namespace {

/* build a RaggedArray2DT from a vector-of-rows */
template <class T>
void BuildRagged(RaggedArray2DT<T>& r, const std::vector<std::vector<T> >& rows)
{
	iArrayT counts(rows.size());
	for (size_t i = 0; i < rows.size(); i++) counts[i] = (int) rows[i].size();
	r.Configure(counts);
	for (size_t i = 0; i < rows.size(); i++) {
		T* p = r(i);
		for (size_t k = 0; k < rows[i].size(); k++) p[k] = rows[i][k];
	}
}

/* reconstruct the full physical field at every constrained node and compare to prescribed */
void ExpectDelivers(const CollocationSolverT& solver, const dArrayT& prescribed,
	dArrayT& coeff /* modified: constrained entries overwritten with the solution */)
{
	dArrayT d_c;
	solver.Solve(prescribed, coeff, d_c);
	const iArrayT& con = solver.Constrained();
	for (int a = 0; a < con.Length(); a++) coeff[con[a]] = d_c[a];
	for (int a = 0; a < con.Length(); a++)
		EXPECT_NEAR(solver.Reconstruct(a, coeff), prescribed[a], 1.0e-10);
}

} /* namespace */

/* Phi_J(x_I) = delta_IJ  ->  coefficients must equal the prescribed physical values */
TEST(Collocation, InterpolatoryLimit)
{
	iArrayT con(3); con[0]=0; con[1]=1; con[2]=2;
	std::vector<std::vector<int> > sup(3); sup[0]={0}; sup[1]={1}; sup[2]={2};
	std::vector<std::vector<double> > phi(3); phi[0]={1.0}; phi[1]={1.0}; phi[2]={1.0};
	RaggedArray2DT<int> S; RaggedArray2DT<double> P; BuildRagged(S,sup); BuildRagged(P,phi);

	CollocationSolverT solver;
	ASSERT_TRUE(solver.SetCollocation(con, 3, S, P));

	dArrayT ubar(3); ubar[0]=5.0; ubar[1]=7.0; ubar[2]=9.0;
	dArrayT coeff(3); coeff=0.0;
	dArrayT d_c; solver.Solve(ubar, coeff, d_c);
	EXPECT_NEAR(d_c[0],5.0,1.0e-12);
	EXPECT_NEAR(d_c[1],7.0,1.0e-12);
	EXPECT_NEAR(d_c[2],9.0,1.0e-12);
}

/* fully non-interpolatory constrained block, no free nodes:
 * Phi_c = [[0.6,0.4],[0.3,0.7]], solve so the physical value at each node = prescribed */
TEST(Collocation, NonInterpolatoryNoFree)
{
	iArrayT con(2); con[0]=0; con[1]=1;
	std::vector<std::vector<int> > sup(2); sup[0]={0,1}; sup[1]={0,1};
	std::vector<std::vector<double> > phi(2); phi[0]={0.6,0.4}; phi[1]={0.3,0.7};
	RaggedArray2DT<int> S; RaggedArray2DT<double> P; BuildRagged(S,sup); BuildRagged(P,phi);

	CollocationSolverT solver;
	ASSERT_TRUE(solver.SetCollocation(con, 2, S, P));

	dArrayT ubar(2); ubar[0]=1.0; ubar[1]=2.0;
	dArrayT coeff(2); coeff=0.0;
	ExpectDelivers(solver, ubar, coeff);

	/* a bare coefficient clamp (d=ubar) would NOT deliver -> confirm collocation differs */
	dArrayT bare(2); bare[0]=1.0; bare[1]=2.0;
	EXPECT_GT(std::fabs(solver.Reconstruct(0,bare) - 1.0), 0.05);
}

/* constrained nodes 0,1 with a FREE node 2 (coeff 10) inside their supports:
 * the free contribution must be subtracted into the RHS so reconstruction still hits prescribed */
TEST(Collocation, FreeNodeContribution)
{
	iArrayT con(2); con[0]=0; con[1]=1;
	std::vector<std::vector<int> > sup(2); sup[0]={0,1,2}; sup[1]={0,1,2};
	std::vector<std::vector<double> > phi(2); phi[0]={0.6,0.3,0.1}; phi[1]={0.2,0.7,0.1};
	RaggedArray2DT<int> S; RaggedArray2DT<double> P; BuildRagged(S,sup); BuildRagged(P,phi);

	CollocationSolverT solver;
	ASSERT_TRUE(solver.SetCollocation(con, 3, S, P));

	dArrayT ubar(2); ubar[0]=1.0; ubar[1]=2.0;
	dArrayT coeff(3); coeff=0.0; coeff[2]=10.0;   /* free node 2 holds coefficient 10 */
	ExpectDelivers(solver, ubar, coeff);          /* coeff[2] stays 10 inside ExpectDelivers */
	EXPECT_NEAR(coeff[2], 10.0, 1.0e-12);
}

/* singular constrained block (duplicate rows) -> SetCollocation reports failure, no crash */
TEST(Collocation, SingularReported)
{
	iArrayT con(2); con[0]=0; con[1]=1;
	std::vector<std::vector<int> > sup(2); sup[0]={0,1}; sup[1]={0,1};
	std::vector<std::vector<double> > phi(2); phi[0]={0.5,0.5}; phi[1]={0.5,0.5}; /* identical rows */
	RaggedArray2DT<int> S; RaggedArray2DT<double> P; BuildRagged(S,sup); BuildRagged(P,phi);
	CollocationSolverT solver;
	EXPECT_FALSE(solver.SetCollocation(con, 2, S, P));
}
