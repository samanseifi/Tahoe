/* test_KLShellShape.cpp — meshfree KL-shell shape-function & geometry tests.
 *
 * Validates the geometry/approximation kernels of the RKPM Kirchhoff-Love shell
 * (development/src/elements/meshfree_kl_shell):
 *
 *   - second-derivative reproduction of the orthogonal-MLS basis (D2OrthoMLS2DT):
 *     a quadratic field's Hessian must be reproduced exactly (the bending prerequisite,
 *     guarded by the #69 fix to the diagonal 2nd-derivative recursion);
 *   - PCA local tangent-plane parameterization + surface normal via RK first derivatives,
 *     converging on a curved (cylindrical) point cloud;
 *   - principal curvatures via the first and second fundamental forms, recovering the
 *     analytic {1, 0} of a unit cylinder.
 */

#include "gtest/gtest.h"

#include <cmath>
#include <vector>

#include "D2OrthoMLS2DT.h"
#include "MLSSolverT.h"
#include "MeshFreeT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "KLShellKernels.h"

using namespace Tahoe;
using namespace Tahoe::KLShell;

namespace {

/* monomial value at (x,y); m: 0:1 1:x 2:y 3:x^2 4:xy 5:y^2 */
double Monomial(int m, double x, double y)
{
	switch (m) {
		case 0: return 1.0;
		case 1: return x;
		case 2: return y;
		case 3: return x * x;
		case 4: return x * y;
		case 5: return y * y;
	}
	return 0.0;
}

/* analytic gradient (d/dx, d/dy) of monomial m evaluated at the local origin (0,0) */
void MonomialGrad(int m, double& gx, double& gy)
{
	gx = gy = 0.0;
	if (m == 1) gx = 1.0;       /* x   */
	else if (m == 2) gy = 1.0;  /* y   */
	/* x^2, xy, y^2 have zero gradient at the origin */
}

/* analytic Hessian (xx, yy, xy) of monomial m */
void MonomialHessian(int m, double& hxx, double& hyy, double& hxy)
{
	hxx = hyy = hxy = 0.0;
	if (m == 3) hxx = 2.0;
	else if (m == 5) hyy = 2.0;
	else if (m == 4) hxy = 1.0;
}

/* neighbors of (px,py) strictly inside radius R over a flat n x n grid on [0,L]^2 */
void FlatGridNeighbors(int n, double L, double px, double py, double R,
	dArray2DT& localCoords, std::vector<double>& gx, std::vector<double>& gy)
{
	double h = L / (n - 1);
	std::vector<double> lx, ly;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) {
			double X = i * h;
			double Y = j * h;
			if (std::sqrt((X - px) * (X - px) + (Y - py) * (Y - py)) < 0.99 * R) {
				lx.push_back(X - px);
				ly.push_back(Y - py);
				gx.push_back(X);
				gy.push_back(Y);
			}
		}
	localCoords.Dimension(int(lx.size()), 2);
	for (int k = 0; k < int(lx.size()); k++) {
		localCoords(k, 0) = lx[k];
		localCoords(k, 1) = ly[k];
	}
}

} /* anonymous namespace */

/* The orthogonal-MLS quadratic basis must reproduce the Hessian of a quadratic field
 * to machine precision. (Diagonal 2nd derivatives were broken before the #69 fix.) */
TEST(KLShellShape, QuadraticHessianReproduction)
{
	const int n = 13;
	const double L = 1.0;
	const double h = L / (n - 1);
	const double R = 3.5 * h;

	D2OrthoMLS2DT efg(2);
	efg.Initialize();

	const double pts[5][2] = {{0.5, 0.5}, {0.37, 0.62}, {0.5, 0.34}, {0.42, 0.5}, {0.6, 0.45}};
	double maxValueErr = 0.0;
	double maxHessianErr = 0.0;

	for (int ip = 0; ip < 5; ip++) {

		dArray2DT lc;
		std::vector<double> gx, gy;
		FlatGridNeighbors(n, L, pts[ip][0], pts[ip][1], R, lc, gx, gy);
		int nn = lc.MajorDim();

		dArrayT dmax(nn);
		dmax = R;
		dArrayT sample(2);
		sample[0] = 0.0;
		sample[1] = 0.0;
		ASSERT_TRUE(efg.SetField(lc, dmax, sample) != 0);

		const dArrayT& phi = efg.phi();
		const dArray2DT& DDphi = efg.DDphi(); /* rows: xx, yy, xy */

		for (int m = 0; m < 6; m++) {

			double sumValue = 0.0, sumXX = 0.0, sumYY = 0.0, sumXY = 0.0;
			for (int I = 0; I < nn; I++) {
				double p = Monomial(m, lc(I, 0), lc(I, 1));
				sumValue += phi[I] * p;
				sumXX += DDphi(0, I) * p;
				sumYY += DDphi(1, I) * p;
				sumXY += DDphi(2, I) * p;
			}

			double hxx, hyy, hxy;
			MonomialHessian(m, hxx, hyy, hxy);
			maxValueErr = std::max(maxValueErr, std::fabs(sumValue - Monomial(m, 0.0, 0.0)));
			maxHessianErr = std::max(maxHessianErr,
				std::fabs(sumXX - hxx) + std::fabs(sumYY - hyy) + std::fabs(sumXY - hxy));
		}
	}

	EXPECT_NEAR(maxValueErr, 0.0, 1e-9);
	EXPECT_NEAR(maxHessianErr, 0.0, 1e-6);
}

/* RKPM (MLSSolverT) with quadratic completeness must also reproduce a quadratic field's
 * Hessian — this exercises the extended PolyBasis2DT (completeness 2, #61) and the general
 * moment-matrix inverse, and is the paper's actual basis (RKPM rather than EFG). */
TEST(KLShellShape, RKPMQuadraticHessianReproduction)
{
	const int n = 13;
	const double L = 1.0;
	const double h = L / (n - 1);
	const double R = 3.5 * h;

	dArrayT window(1);
	window[0] = 1.2; /* cubic-spline support scaling */
	MLSSolverT rkpm(2, 2, false, MeshFreeT::kCubicSpline, window);
	rkpm.Initialize();

	const double pts[3][2] = {{0.5, 0.5}, {0.37, 0.62}, {0.42, 0.5}};
	double maxValueErr = 0.0, maxGradErr = 0.0, maxHessianErr = 0.0;

	for (int ip = 0; ip < 3; ip++) {

		dArray2DT lc;
		std::vector<double> gx, gy;
		FlatGridNeighbors(n, L, pts[ip][0], pts[ip][1], R, lc, gx, gy);
		int nn = lc.MajorDim();

		dArray2DT nodalParam(nn, 1);
		nodalParam = R;
		dArrayT volume(nn);
		volume = h * h;
		dArrayT sample(2);
		sample[0] = 0.0;
		sample[1] = 0.0;
		ASSERT_TRUE(rkpm.SetField(lc, nodalParam, volume, sample, 2) != 0);

		const dArrayT& phi = rkpm.phi();
		const dArray2DT& Dphi = rkpm.Dphi();
		const dArray2DT& DDphi = rkpm.DDphi();

		for (int m = 0; m < 6; m++) {
			double sv = 0, sgx = 0, sgy = 0, sxx = 0, syy = 0, sxy = 0;
			for (int I = 0; I < nn; I++) {
				double p = Monomial(m, lc(I, 0), lc(I, 1));
				sv += phi[I] * p;
				sgx += Dphi(0, I) * p;
				sgy += Dphi(1, I) * p;
				sxx += DDphi(0, I) * p;
				syy += DDphi(1, I) * p;
				sxy += DDphi(2, I) * p;
			}
			double gxv, gyv, hxx, hyy, hxy;
			MonomialGrad(m, gxv, gyv);
			MonomialHessian(m, hxx, hyy, hxy);
			maxValueErr = std::max(maxValueErr, std::fabs(sv - Monomial(m, 0.0, 0.0)));
			maxGradErr = std::max(maxGradErr, std::fabs(sgx - gxv) + std::fabs(sgy - gyv));
			maxHessianErr = std::max(maxHessianErr,
				std::fabs(sxx - hxx) + std::fabs(syy - hyy) + std::fabs(sxy - hxy));
		}
	}

	EXPECT_NEAR(maxValueErr, 0.0, 1e-9);
	EXPECT_NEAR(maxGradErr, 0.0, 1e-7);
	EXPECT_NEAR(maxHessianErr, 0.0, 1e-5);
}

/* PCA parameterization + RK first-derivative normal converge on a cylinder. */
TEST(KLShellShape, CylinderNormalAndCurvature)
{
	const double Rc = 1.0;
	const double Zlen = 2.0;

	double prevNormalErr = -1.0;
	double prevCurvErr = -1.0;
	int refinements = 0;

	const int levels[3][2] = {{24, 9}, {48, 17}, {96, 33}};
	for (int lv = 0; lv < 3; lv++) {

		int nt = levels[lv][0];
		int nz = levels[lv][1];
		double dth = 2.0 * M_PI / nt;
		double dz = Zlen / (nz - 1);
		double spacing = (Rc * dth > dz) ? Rc * dth : dz;
		double R = 3.5 * spacing;

		std::vector<double> X, Y, Z, TH;
		for (int i = 0; i < nt; i++)
			for (int j = 0; j < nz; j++) {
				double th = i * dth;
				X.push_back(Rc * std::cos(th));
				Y.push_back(Rc * std::sin(th));
				Z.push_back(j * dz);
				TH.push_back(th);
			}
		int N = int(X.size());

		D2OrthoMLS2DT efg(2);
		efg.Initialize();

		double sumNormalSq = 0.0;
		double sumCurvSq = 0.0;
		int count = 0;

		for (int P = 0; P < N; P++) {

			/* interior in the axial direction only */
			if (Z[P] < 0.5 * dz || Z[P] > Zlen - 0.5 * dz) continue;

			std::vector<int> nb;
			for (int Q = 0; Q < N; Q++) {
				double d0 = X[Q] - X[P], d1 = Y[Q] - Y[P], d2 = Z[Q] - Z[P];
				if (std::sqrt(d0 * d0 + d1 * d1 + d2 * d2) < 0.99 * R) nb.push_back(Q);
			}
			int nn = int(nb.size());
			if (nn < 6) continue;

			std::vector<double> nbX(3 * nn);
			for (int k = 0; k < nn; k++) {
				nbX[3 * k] = X[nb[k]];
				nbX[3 * k + 1] = Y[nb[k]];
				nbX[3 * k + 2] = Z[nb[k]];
			}
			double psi1[3], psi2[3], n0[3];
			PCAFrame(&nbX[0], nn, psi1, psi2, n0);

			dArray2DT lc(nn, 2);
			for (int k = 0; k < nn; k++) {
				double dxv[3] = {X[nb[k]] - X[P], Y[nb[k]] - Y[P], Z[nb[k]] - Z[P]};
				lc(k, 0) = Dot(dxv, psi1);
				lc(k, 1) = Dot(dxv, psi2);
			}

			dArrayT dmax(nn);
			dmax = R;
			dArrayT sample(2);
			sample[0] = 0.0;
			sample[1] = 0.0;
			if (!efg.SetField(lc, dmax, sample)) continue;
			const dArray2DT& Dphi = efg.Dphi();
			const dArray2DT& DDphi = efg.DDphi();

			double x1[3] = {0,0,0}, x2[3] = {0,0,0};
			double x11[3] = {0,0,0}, x22[3] = {0,0,0}, x12[3] = {0,0,0};
			for (int I = 0; I < nn; I++) {
				double Xq[3] = {X[nb[I]], Y[nb[I]], Z[nb[I]]};
				for (int d = 0; d < 3; d++) {
					x1[d] += Dphi(0, I) * Xq[d];
					x2[d] += Dphi(1, I) * Xq[d];
					x11[d] += DDphi(0, I) * Xq[d];
					x22[d] += DDphi(1, I) * Xq[d];
					x12[d] += DDphi(2, I) * Xq[d];
				}
			}

			double nvec[3];
			Cross(x1, x2, nvec);
			double nm = Norm(nvec);
			if (nm < 1e-14) continue;
			for (int d = 0; d < 3; d++) nvec[d] /= nm;

			/* normal error vs analytic (cos t, sin t, 0) */
			double na[3] = {std::cos(TH[P]), std::sin(TH[P]), 0.0};
			double sgn = (Dot(nvec, na) < 0.0) ? -1.0 : 1.0;
			double ne = 0.0;
			for (int d = 0; d < 3; d++) ne += (sgn * nvec[d] - na[d]) * (sgn * nvec[d] - na[d]);
			sumNormalSq += ne;

			/* principal curvatures via first/second fundamental forms */
			double E = Dot(x1, x1), F = Dot(x1, x2), G = Dot(x2, x2);
			double Lf = Dot(x11, nvec), Mf = Dot(x12, nvec), Nf = Dot(x22, nvec);
			double detI = E * G - F * F;
			if (std::fabs(detI) > 1e-14) {
				double a2 = detI;
				double a1 = -(E * Nf - 2.0 * F * Mf + G * Lf);
				double a0 = Lf * Nf - Mf * Mf;
				double disc = a1 * a1 - 4.0 * a2 * a0;
				if (disc < 0.0) disc = 0.0;
				double k1 = std::fabs((-a1 + std::sqrt(disc)) / (2.0 * a2));
				double k2 = std::fabs((-a1 - std::sqrt(disc)) / (2.0 * a2));
				double kmax = (k1 > k2) ? k1 : k2;
				double kmin = (k1 > k2) ? k2 : k1;
				double ce = std::fabs(kmax - 1.0) + std::fabs(kmin); /* cylinder: {1, 0} */
				sumCurvSq += ce * ce;
			}
			count++;
		}

		ASSERT_GT(count, 0);
		double normalErr = std::sqrt(sumNormalSq / count);
		double curvErr = std::sqrt(sumCurvSq / count);

		if (prevNormalErr > 0.0) EXPECT_LT(normalErr, prevNormalErr); /* converging */
		if (prevCurvErr > 0.0) EXPECT_LT(curvErr, prevCurvErr);
		prevNormalErr = normalErr;
		prevCurvErr = curvErr;
		refinements++;
	}

	ASSERT_EQ(refinements, 3);
	EXPECT_LT(prevNormalErr, 1e-5); /* finest level */
	EXPECT_LT(prevCurvErr, 5e-3);
}
