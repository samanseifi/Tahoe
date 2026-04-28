/* ANPHelperT.cpp — Bonet-Burton 1998 average nodal pressure helper.
 *
 * The two-pass loop below uses #pragma omp atomic on the nodal scatter so
 * the helper is safe inside the OpenMP-parallel batch loop in
 * ExplicitElementT.  For classic single-threaded use the atomics compile
 * down to plain stores.
 */

#include "ANPHelperT.h"
#include <cstring>
#include <cmath>

using namespace Tahoe;

ANPHelperT::ANPHelperT(void)
	: fNelem(0), fNnod(0), fNen(0),
	  fConn(NULL), fVrefE(NULL),
	  fSumVrefN(NULL), fJN(NULL)
{
}

ANPHelperT::~ANPHelperT(void)
{
	delete[] fSumVrefN;
	delete[] fJN;
}

void ANPHelperT::Init(int nelem, int nnod, int nen,
                      const int* conn, const double* V_ref_e)
{
	fNelem = nelem;
	fNnod  = nnod;
	fNen   = nen;
	fConn  = conn;
	fVrefE = V_ref_e;

	delete[] fSumVrefN;
	delete[] fJN;
	fSumVrefN = new double[nnod];
	fJN       = new double[nnod];

	/* precompute per-node denominator: sum_{e ∋ n} V_e^ref */
	std::memset(fSumVrefN, 0, sizeof(double) * nnod);
	for (int e = 0; e < nelem; e++) {
		double V = V_ref_e[e];
		const int* ec = conn + e * nen;
		for (int n = 0; n < nen; n++)
			fSumVrefN[ec[n]] += V;
	}
}

void ANPHelperT::ComputeJBar(const double* J_e, double* J_bar_e)
{
	/* PASS 1: scatter V_e^ref * J_e into nodal accumulator. */
	std::memset(fJN, 0, sizeof(double) * fNnod);
	#pragma omp parallel for if(fNelem > 1024)
	for (int e = 0; e < fNelem; e++) {
		double VJ = fVrefE[e] * J_e[e];
		const int* ec = fConn + e * fNen;
		for (int n = 0; n < fNen; n++) {
			#pragma omp atomic
			fJN[ec[n]] += VJ;
		}
	}

	/* normalise: J_n /= sum_{e ∋ n} V_e^ref */
	#pragma omp parallel for if(fNnod > 1024)
	for (int n = 0; n < fNnod; n++) {
		double denom = fSumVrefN[n];
		fJN[n] = (denom > 0.0) ? fJN[n] / denom : 1.0;
	}

	/* PASS 2: gather J_n back to J_bar_e (average over the element's nen nodes). */
	double inv_nen = 1.0 / (double)fNen;
	#pragma omp parallel for if(fNelem > 1024)
	for (int e = 0; e < fNelem; e++) {
		const int* ec = fConn + e * fNen;
		double sum = 0.0;
		for (int n = 0; n < fNen; n++)
			sum += fJN[ec[n]];
		J_bar_e[e] = sum * inv_nen;
	}
}

/* F_bar = (J_bar / J)^(1/3) * F  for 3D.
 * F is laid out as 9 contiguous rows of length stride. */
void ANPHelperT::ApplyFBar3D(int nelem, int stride,
	const double* F_e_ptr,
	const double* J_e,
	const double* J_bar_e,
	double* F_bar_e_ptr)
{
	for (int i = 0; i < nelem; i++) {
		double J = J_e[i];
		double Jb = J_bar_e[i];
		double scale = (J > 0.0) ? std::cbrt(Jb / J) : 1.0;
		for (int k = 0; k < 9; k++)
			F_bar_e_ptr[k * stride + i] = scale * F_e_ptr[k * stride + i];
	}
}
