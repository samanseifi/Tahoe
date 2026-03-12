/* created: (2026) */
#include "NLSolver_NK.h"

#include <cmath>
#include <iostream>

#include "ExceptionT.h"
#include "FEManagerT.h"
#include "GlobalMatrixT.h"
#include "ParameterListT.h"
#include "LimitT.h"

using namespace Tahoe;

/* constructor */
NLSolver_NK::NLSolver_NK(FEManagerT& fe_manager, int group):
	NLSolver(fe_manager, group),
	fKrylovRestart(30),
	fLinearTolerance(1.0e-3),
	fMaxLinearIter(0)
{
	SetName("newton_krylov_solver");

	/* console variables */
	iAddVariable("gmres_restart",       fKrylovRestart);
	iAddVariable("linear_tolerance",    fLinearTolerance);
	iAddVariable("max_linear_iter",     fMaxLinearIter);
}

/* describe the parameters needed by the interface */
void NLSolver_NK::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	NLSolver::DefineParameters(list);

	/* GMRES restart dimension */
	ParameterT restart(fKrylovRestart, "gmres_restart");
	restart.AddLimit(1, LimitT::LowerInclusive);
	restart.SetDefault(fKrylovRestart);
	list.AddParameter(restart);

	/* inner linear solve tolerance */
	ParameterT lin_tol(fLinearTolerance, "linear_tolerance");
	lin_tol.AddLimit(0.0, LimitT::Lower);
	lin_tol.SetDefault(fLinearTolerance);
	list.AddParameter(lin_tol);

	/* maximum total GMRES iterations (0 = no limit) */
	ParameterT max_lin(fMaxLinearIter, "max_linear_iter");
	max_lin.AddLimit(0, LimitT::LowerInclusive);
	max_lin.SetDefault(fMaxLinearIter);
	list.AddParameter(max_lin);
}

/* accept parameter list */
void NLSolver_NK::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	NLSolver::TakeParameterList(list);

	fKrylovRestart    = list.GetParameter("gmres_restart");
	fLinearTolerance  = list.GetParameter("linear_tolerance");
	fMaxLinearIter    = list.GetParameter("max_linear_iter");
}

/*************************************************************************
 * Protected
 *************************************************************************/

/* one Newton-Krylov iteration: solve K*du = rhs using GMRES,
 * then apply the update */
void NLSolver_NK::Iterate(void)
{
	int n = fRHS.Length();

	/* allocate workspace if needed */
	int m = fKrylovRestart;
	if (fQ.Length() != m + 1) {
		fQ.Dimension(m + 1);
		for (int i = 0; i <= m; i++) fQ[i].Dimension(n);
		fH.Dimension(m + 1, m);
		fg.Dimension(m + 1);
		fcs.Dimension(m);
		fsn.Dimension(m);
		fy.Dimension(m);
		fz.Dimension(n);
		fw.Dimension(n);
		fdiag.Dimension(n);
	}

	/* copy RHS (GMRES solves in place; we need the original RHS intact) */
	dArrayT rhs(n);
	rhs = fRHS;

	/* initial guess = zero */
	dArrayT update(n);
	update = 0.0;

	/* solve K * update = rhs using GMRES */
	bool converged = GMRES(update, rhs);
	if (!converged)
		cout << "\n NLSolver_NK::Iterate: GMRES did not converge to "
		     << fLinearTolerance << " (continuing Newton iteration)\n";

	/* apply update to the system */
	Update(update, NULL);
}

/*************************************************************************
 * Private
 *************************************************************************/

/** restarted GMRES(m) with left Jacobi preconditioning.
 *
 * Solves  A * x = b  where A = fLHS, b = rhs.
 * The left preconditioner is M = diag(A); i.e., we actually solve
 *   (M^{-1} A) x = M^{-1} b.
 *
 * Algorithm: Saad & Schultz (1986), SIAM J. Sci. Stat. Comput. 7(3).
 */
bool NLSolver_NK::GMRES(dArrayT& x, const dArrayT& b) const
{
	int n  = x.Length();
	int m  = fKrylovRestart;

	/* -----------------------------------------------------------
	 * build Jacobi preconditioner: fdiag[i] = 1/A_{ii}
	 * fall back to identity when diagonal is zero or unavailable
	 * ----------------------------------------------------------- */
	fdiag = 1.0;
	dArrayT raw_diag(n);
	if (fLHS->CopyDiagonal(raw_diag)) {
		for (int i = 0; i < n; i++)
			fdiag[i] = (fabs(raw_diag[i]) > 1.0e-14) ? 1.0/raw_diag[i] : 1.0;
	}

	/* maximum total iterations */
	int max_iter = (fMaxLinearIter > 0) ? fMaxLinearIter : (3 * n);

	int total_iter = 0;
	bool done      = false;

	/* outer restart loop */
	while (!done && total_iter < max_iter)
	{
		/* r = b - A*x */
		fLHS->Multx(x, fw);                           /* fw = A*x     */
		fw.SetToCombination(1.0, b, -1.0, fw);        /* fw = b - A*x */

		/* apply preconditioner: fw <- M^{-1} fw */
		for (int i = 0; i < n; i++) fw[i] *= fdiag[i];

		double beta = std::sqrt(nArrayT<double>::Dot(fw, fw));
		if (beta == 0.0) { done = true; break; }

		/* check convergence on initial residual */
		if (beta < fLinearTolerance) { done = true; break; }

		/* q[0] = fw / beta */
		fQ[0].SetToScaled(1.0/beta, fw);

		/* initialize rotated RHS */
		fg    = 0.0;
		fg[0] = beta;

		int j_end = 0; /* last Arnoldi step taken */

		/* inner Arnoldi loop */
		for (int j = 0; j < m && total_iter < max_iter; j++, total_iter++)
		{
			j_end = j;

			/* w = M^{-1} A q[j] */
			fLHS->Multx(fQ[j], fw);                       /* fw = A*q[j]         */
			for (int i = 0; i < n; i++) fw[i] *= fdiag[i]; /* fw = M^{-1}*A*q[j] */

			/* modified Gram-Schmidt orthogonalization */
			for (int i = 0; i <= j; i++)
			{
				double hij = nArrayT<double>::Dot(fw, fQ[i]);
				fH(i, j)   = hij;
				/* fw -= hij * q[i] */
				for (int k = 0; k < n; k++) fw[k] -= hij * fQ[i][k];
			}
			double h_jp1_j = std::sqrt(nArrayT<double>::Dot(fw, fw));
			fH(j+1, j)     = h_jp1_j;

			/* new basis vector */
			if (h_jp1_j > 0.0)
				fQ[j+1].SetToScaled(1.0/h_jp1_j, fw);

			/* apply previous Givens rotations to column j */
			for (int i = 0; i < j; i++)
				ApplyGivens(fcs[i], fsn[i], fH(i,j), fH(i+1,j));

			/* new Givens rotation to eliminate H(j+1,j) */
			double denom = std::sqrt(fH(j,j)*fH(j,j) + fH(j+1,j)*fH(j+1,j));
			if (denom > 0.0) {
				fcs[j] =  fH(j,j)   / denom;
				fsn[j] =  fH(j+1,j) / denom;
			} else {
				fcs[j] = 1.0;
				fsn[j] = 0.0;
			}
			ApplyGivens(fcs[j], fsn[j], fH(j,j), fH(j+1,j));
			ApplyGivens(fcs[j], fsn[j], fg[j],   fg[j+1]);

			double res_norm = fabs(fg[j+1]);
			if (res_norm < fLinearTolerance) {
				j_end = j;
				done  = true;
				break;
			}
		}

		/* -----------------------------------------------------------
		 * back-substitute to find y: H(0..j_end,0..j_end) * y = g
		 * H is upper triangular after Givens rotations
		 * ----------------------------------------------------------- */
		int k = j_end;
		fy = 0.0;
		for (int i = k; i >= 0; i--)
		{
			double sum = fg[i];
			for (int l = i+1; l <= k; l++)
				sum -= fH(i,l) * fy[l];
			fy[i] = (fabs(fH(i,i)) > 1.0e-14) ? sum / fH(i,i) : 0.0;
		}

		/* x += Q[:,0..k] * y */
		for (int i = 0; i <= k; i++)
			for (int p = 0; p < n; p++)
				x[p] += fy[i] * fQ[i][p];
	}

	return done;
}

/* apply Givens rotation [h1; h2] <- [c s; -s c] [h1; h2] */
void NLSolver_NK::ApplyGivens(double c, double s, double& h1, double& h2)
{
	double temp = c*h1 + s*h2;
	h2          = -s*h1 + c*h2;
	h1          = temp;
}
