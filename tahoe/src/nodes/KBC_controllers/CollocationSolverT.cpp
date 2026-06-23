/* CollocationSolverT.cpp */
#include "CollocationSolverT.h"
#include <cmath>
#include "ExceptionT.h"

using namespace Tahoe;

CollocationSolverT::CollocationSolverT(void): fNc(0) {}

bool CollocationSolverT::SetCollocation(const iArrayT& constrained, int max_global_id,
	const RaggedArray2DT<int>& support, const RaggedArray2DT<double>& phi)
{
	fConstrained = constrained;
	fSupport = support;
	fPhi = phi;
	fNc = constrained.Length();

	/* global id -> constrained local index */
	fGlobalToCon.Dimension(max_global_id);
	fGlobalToCon = -1;
	for (int a = 0; a < fNc; a++) {
		int g = constrained[a];
		if (g < 0 || g >= max_global_id)
			ExceptionT::OutOfRange("CollocationSolverT::SetCollocation", "node %d out of range", g);
		fGlobalToCon[g] = a;
	}

	/* assemble Phi_c (Nc x Nc, row-major): Phi_c[a][b] = Phi_{constrained[b]}(x_{constrained[a]}) */
	fLU.Dimension(fNc*fNc);
	fLU = 0.0;
	for (int a = 0; a < fNc; a++) {
		int m = fSupport.MinorDim(a);
		const int* sj = fSupport(a);
		const double* pj = fPhi(a);
		for (int k = 0; k < m; k++) {
			int b = (sj[k] >= 0 && sj[k] < max_global_id) ? fGlobalToCon[sj[k]] : -1;
			if (b >= 0) fLU[a*fNc + b] += pj[k];
		}
	}

	/* LU factorization with partial pivoting (Doolittle), in place */
	fPiv.Dimension(fNc);
	for (int i = 0; i < fNc; i++) fPiv[i] = i;
	for (int col = 0; col < fNc; col++) {
		/* pivot: largest magnitude in this column at/below the diagonal */
		int piv = col; double best = std::fabs(fLU[col*fNc + col]);
		for (int r = col+1; r < fNc; r++) {
			double v = std::fabs(fLU[r*fNc + col]);
			if (v > best) { best = v; piv = r; }
		}
		if (best < 1.0e-300) return false;   /* singular constrained set */
		if (piv != col) {
			for (int c = 0; c < fNc; c++) {
				double t = fLU[col*fNc + c]; fLU[col*fNc + c] = fLU[piv*fNc + c]; fLU[piv*fNc + c] = t;
			}
			int ti = fPiv[col]; fPiv[col] = fPiv[piv]; fPiv[piv] = ti;
		}
		double diag = fLU[col*fNc + col];
		for (int r = col+1; r < fNc; r++) {
			double f = fLU[r*fNc + col] / diag;
			fLU[r*fNc + col] = f;
			for (int c = col+1; c < fNc; c++) fLU[r*fNc + c] -= f*fLU[col*fNc + c];
		}
	}
	return true;
}

void CollocationSolverT::Solve(const dArrayT& prescribed, const dArrayT& coeff_dof, dArrayT& d_c) const
{
	if (prescribed.Length() != fNc)
		ExceptionT::SizeMismatch("CollocationSolverT::Solve");
	d_c.Dimension(fNc);

	/* RHS = prescribed - Phi_free . d_free (only the NON-constrained support contributions) */
	dArrayT rhs(fNc);
	for (int a = 0; a < fNc; a++) {
		double free_sum = 0.0;
		int m = fSupport.MinorDim(a);
		const int* sj = fSupport(a);
		const double* pj = fPhi(a);
		for (int k = 0; k < m; k++) {
			int g = sj[k];
			if (g < 0 || g >= fGlobalToCon.Length()) continue;
			if (fGlobalToCon[g] < 0)                     /* free node -> contributes to RHS */
				free_sum += pj[k]*coeff_dof[g];
		}
		rhs[a] = prescribed[a] - free_sum;
	}

	/* apply row permutation, then forward/back substitution */
	dArrayT y(fNc);
	for (int i = 0; i < fNc; i++) y[i] = rhs[fPiv[i]];
	for (int i = 0; i < fNc; i++) {                  /* L y (unit diagonal) */
		double s = y[i];
		for (int j = 0; j < i; j++) s -= fLU[i*fNc + j]*y[j];
		y[i] = s;
	}
	for (int i = fNc-1; i >= 0; i--) {               /* U x = y */
		double s = y[i];
		for (int j = i+1; j < fNc; j++) s -= fLU[i*fNc + j]*d_c[j];
		d_c[i] = s / fLU[i*fNc + i];
	}
}

double CollocationSolverT::Reconstruct(int a, const dArrayT& coeff_dof) const
{
	double u = 0.0;
	int m = fSupport.MinorDim(a);
	const int* sj = fSupport(a);
	const double* pj = fPhi(a);
	for (int k = 0; k < m; k++)
		if (sj[k] >= 0 && sj[k] < coeff_dof.Length()) u += pj[k]*coeff_dof[sj[k]];
	return u;
}
