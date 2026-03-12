/* Unit tests for the GMRES linear solver underlying NLSolver_NK.
 *
 * These tests exercise GMRES on small, dense linear systems whose solutions
 * are known analytically so that correctness can be verified without
 * building the full FEM infrastructure. */
#include <gtest/gtest.h>
#include <cmath>
#include <vector>

/* -----------------------------------------------------------------------
 * Minimal, self-contained GMRES(m) implementation mirroring NLSolver_NK.
 * This lets us test the algorithm in isolation.
 * ----------------------------------------------------------------------- */

/** Apply Givens rotation: [h1; h2] <- [c s; -s c] [h1; h2] */
static void ApplyGivens(double c, double s, double& h1, double& h2)
{
    double t = c*h1 + s*h2;
    h2       = -s*h1 + c*h2;
    h1       = t;
}

/** Dense matrix-vector product y = A*x (row-major, n x n). */
static void matvec(const std::vector<std::vector<double>>& A,
                   const std::vector<double>& x,
                   std::vector<double>& y)
{
    int n = static_cast<int>(x.size());
    y.assign(n, 0.0);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            y[i] += A[i][j] * x[j];
}

/** dot product */
static double dot(const std::vector<double>& a, const std::vector<double>& b)
{
    double s = 0.0;
    for (size_t i = 0; i < a.size(); i++) s += a[i]*b[i];
    return s;
}

/** vector 2-norm */
static double norm2(const std::vector<double>& v)
{
    return std::sqrt(dot(v, v));
}

/**
 * GMRES(m) with diagonal (Jacobi) left preconditioner.
 *
 * \param A   n x n matrix (row-major nested vector)
 * \param x   initial guess on entry; solution on exit
 * \param b   right-hand side
 * \param m   restart dimension
 * \param tol relative residual tolerance
 * \return    achieved relative residual
 */
static double GMRES(const std::vector<std::vector<double>>& A,
                    std::vector<double>& x,
                    const std::vector<double>& b,
                    int m, double tol)
{
    int n = static_cast<int>(b.size());

    /* diagonal preconditioner: diag[i] = 1/A[i][i] */
    std::vector<double> diag(n, 1.0);
    for (int i = 0; i < n; i++)
        if (std::fabs(A[i][i]) > 1.0e-14) diag[i] = 1.0/A[i][i];

    /* compute initial residual norm for convergence check */
    std::vector<double> tmp(n);
    matvec(A, x, tmp);
    std::vector<double> r0(n);
    for (int i = 0; i < n; i++) r0[i] = (b[i] - tmp[i]) * diag[i];
    double beta0 = norm2(r0);
    if (beta0 < 1.0e-14) return 0.0;

    int max_iter = 10 * n;
    int total    = 0;

    /* storage */
    std::vector<std::vector<double>> Q(m+1, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> H(m+1, std::vector<double>(m, 0.0));
    std::vector<double> g(m+1), cs(m), sn(m), y(m);

    double res_ratio = 1.0;

    while (total < max_iter && res_ratio > tol)
    {
        /* r = M^{-1}(b - A*x) */
        matvec(A, x, tmp);
        for (int i = 0; i < n; i++) r0[i] = (b[i] - tmp[i]) * diag[i];
        double beta = norm2(r0);
        if (beta < 1.0e-14) break;

        res_ratio = beta / beta0;
        if (res_ratio <= tol) break;

        for (int i = 0; i <= m; i++) std::fill(Q[i].begin(), Q[i].end(), 0.0);
        for (int i = 0; i <= m; i++) std::fill(H[i].begin(), H[i].end(), 0.0);
        std::fill(g.begin(), g.end(), 0.0);
        g[0] = beta;
        for (int i = 0; i < n; i++) Q[0][i] = r0[i] / beta;

        int j_end = 0;
        bool inner_done = false;
        for (int j = 0; j < m && total < max_iter; j++, total++)
        {
            j_end = j;
            /* w = M^{-1} A q[j] */
            matvec(A, Q[j], tmp);
            for (int i = 0; i < n; i++) tmp[i] *= diag[i];

            /* modified Gram-Schmidt */
            for (int i = 0; i <= j; i++) {
                double h = dot(tmp, Q[i]);
                H[i][j]  = h;
                for (int k = 0; k < n; k++) tmp[k] -= h * Q[i][k];
            }
            double hj1 = norm2(tmp);
            H[j+1][j]  = hj1;
            if (hj1 > 0.0)
                for (int k = 0; k < n; k++) Q[j+1][k] = tmp[k] / hj1;

            /* apply previous Givens rotations */
            for (int i = 0; i < j; i++)
                ApplyGivens(cs[i], sn[i], H[i][j], H[i+1][j]);

            /* new rotation */
            double denom = std::sqrt(H[j][j]*H[j][j] + H[j+1][j]*H[j+1][j]);
            if (denom > 0.0) { cs[j]=H[j][j]/denom; sn[j]=H[j+1][j]/denom; }
            else              { cs[j]=1.0; sn[j]=0.0; }
            ApplyGivens(cs[j], sn[j], H[j][j], H[j+1][j]);
            ApplyGivens(cs[j], sn[j], g[j],    g[j+1]);

            res_ratio = std::fabs(g[j+1]) / beta0;
            if (res_ratio <= tol) { inner_done = true; break; }
        }

        /* back-substitution */
        int k = j_end;
        std::fill(y.begin(), y.end(), 0.0);
        for (int i = k; i >= 0; i--) {
            double s = g[i];
            for (int l = i+1; l <= k; l++) s -= H[i][l] * y[l];
            y[i] = (std::fabs(H[i][i]) > 1.0e-14) ? s / H[i][i] : 0.0;
        }
        for (int i = 0; i <= k; i++)
            for (int p = 0; p < n; p++)
                x[p] += y[i] * Q[i][p];

        if (inner_done) break;
    }

    return res_ratio;
}

/* -----------------------------------------------------------------------
 * Tests
 * ----------------------------------------------------------------------- */

/** 2x2 symmetric positive-definite system */
TEST(GMRES, SPD_2x2) {
    std::vector<std::vector<double>> A = {{4.0, 1.0}, {1.0, 3.0}};
    std::vector<double> b = {1.0, 2.0};
    std::vector<double> x(2, 0.0);
    double res = GMRES(A, x, b, 10, 1.0e-10);
    // exact solution: x = [1/11, 7/11]
    EXPECT_NEAR(x[0],  1.0/11.0, 1.0e-10);
    EXPECT_NEAR(x[1],  7.0/11.0, 1.0e-10);
    EXPECT_LT(res, 1.0e-10);
}

/** 3x3 diagonal system (solved in one iteration) */
TEST(GMRES, Diagonal_3x3) {
    std::vector<std::vector<double>> A = {
        {2.0, 0.0, 0.0},
        {0.0, 3.0, 0.0},
        {0.0, 0.0, 5.0}};
    std::vector<double> b = {4.0, 9.0, 5.0};
    std::vector<double> x(3, 0.0);
    double res = GMRES(A, x, b, 10, 1.0e-12);
    EXPECT_NEAR(x[0], 2.0, 1.0e-12);
    EXPECT_NEAR(x[1], 3.0, 1.0e-12);
    EXPECT_NEAR(x[2], 1.0, 1.0e-12);
    EXPECT_LT(res, 1.0e-12);
}

/** 5x5 non-symmetric system */
TEST(GMRES, NonSymmetric_5x5) {
    // Strictly diagonally dominant → unique solution
    std::vector<std::vector<double>> A = {
        {10.0, -1.0,  2.0,  0.0,  0.0},
        {-1.0,  11.0, -1.0,  3.0,  0.0},
        { 2.0, -1.0, 10.0, -1.0,  0.0},
        { 0.0,  3.0, -1.0,  8.0, -2.0},
        { 0.0,  0.0,  0.0, -2.0,  6.0}};
    std::vector<double> x_exact = {1.0, 2.0, 3.0, 4.0, 5.0};
    // compute b = A * x_exact
    std::vector<double> b(5, 0.0);
    matvec(A, x_exact, b);
    std::vector<double> x(5, 0.0);
    double res = GMRES(A, x, b, 10, 1.0e-10);
    for (int i = 0; i < 5; i++)
        EXPECT_NEAR(x[i], x_exact[i], 1.0e-8);
    EXPECT_LT(res, 1.0e-10);
}

/** Restart test: m=2, forcing multiple restart cycles for a 4x4 system */
TEST(GMRES, Restart_4x4) {
    std::vector<std::vector<double>> A = {
        {4.0, 1.0, 0.0, 0.0},
        {1.0, 4.0, 1.0, 0.0},
        {0.0, 1.0, 4.0, 1.0},
        {0.0, 0.0, 1.0, 4.0}};
    std::vector<double> x_exact = {1.0, -1.0, 2.0, -2.0};
    std::vector<double> b(4, 0.0);
    matvec(A, x_exact, b);
    std::vector<double> x(4, 0.0);
    /* use restart=2 to exercise the restart logic */
    double res = GMRES(A, x, b, 2, 1.0e-10);
    for (int i = 0; i < 4; i++)
        EXPECT_NEAR(x[i], x_exact[i], 1.0e-8);
    EXPECT_LT(res, 1.0e-10);
}

/** Zero right-hand side: solution is zero */
TEST(GMRES, ZeroRHS) {
    std::vector<std::vector<double>> A = {{3.0, 1.0}, {1.0, 2.0}};
    std::vector<double> b = {0.0, 0.0};
    std::vector<double> x = {0.0, 0.0};
    double res = GMRES(A, x, b, 5, 1.0e-12);
    EXPECT_NEAR(x[0], 0.0, 1.0e-14);
    EXPECT_NEAR(x[1], 0.0, 1.0e-14);
    EXPECT_LT(res, 1.0e-12);
}

/** Givens rotation correctness */
TEST(GMRES, GivensRotation) {
    double h1 = 3.0, h2 = 4.0;
    double denom = std::sqrt(h1*h1 + h2*h2); // 5
    double c = h1/denom, s = h2/denom;
    ApplyGivens(c, s, h1, h2);
    EXPECT_NEAR(h1,  5.0, 1.0e-14);
    EXPECT_NEAR(h2,  0.0, 1.0e-14);
}
