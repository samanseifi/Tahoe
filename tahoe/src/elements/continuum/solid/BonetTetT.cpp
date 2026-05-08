/* BonetTetT.cpp */
#include "BonetTetT.h"

#include "ANPHelperT.h"
#include "ElementSupportT.h"
#include "ModelManagerT.h"
#include "FiniteStrainT.h"
#include "ElementCardT.h"
#include "ShapeFunctionT.h"
#include "iArrayT.h"
#include "dArray2DT.h"
#include "dMatrixT.h"

#include <iostream>
#include <cmath>

using namespace Tahoe;

BonetTetT::BonetTetT(const ElementSupportT& support)
	: UpdatedLagrangianT(support),
	  fANP(NULL),
	  fFlatConn(NULL),
	  fVrefE(NULL),
	  fJe(NULL),
	  fJbarE(NULL),
	  fNelem(0)
{
	SetName("bonet_tet");
}

BonetTetT::~BonetTetT(void)
{
	delete fANP;
	delete[] fFlatConn;
	delete[] fVrefE;
	delete[] fJe;
	delete[] fJbarE;
}

void BonetTetT::TakeParameterList(const ParameterListT& list)
{
	/* base class — sets up connectivity, shape funcs, materials, etc. */
	UpdatedLagrangianT::TakeParameterList(list);

	/* this element only makes sense for Tet4 */
	if (NumElementNodes() != 4 || NumSD() != 3) {
		std::cout << "BonetTetT: WARNING — only valid for 3D Tet4 (nen=4, nsd=3); "
		          << "got nen=" << NumElementNodes() << " nsd=" << NumSD()
		          << ".  ANP F-bar disabled." << std::endl;
		return;
	}

	BuildANPData();
	std::cout << "BonetTetT: ANP-Tet4 (LS-DYNA ELFORM=13) active for "
	          << fNelem << " tets" << std::endl;

	/* Populate J̄_e from the reference coordinates so the very first
	 * FormRHS (which can fire before the first InitStep depending on
	 * the time-stepping driver) sees a sane J̄ = 1.0 across the mesh
	 * — F̄ = F at undeformed, no spurious volumetric scaling. */
	UpdateJBar();
}

/* Build flat connectivity, reference volumes, and instantiate ANP helper. */
void BonetTetT::BuildANPData(void)
{
	/* count total elements */
	fNelem = 0;
	for (int b = 0; b < fBlockData.Length(); b++)
		fNelem += fBlockData[b].Dimension();

	const int nen = 4;
	delete[] fFlatConn; delete[] fVrefE; delete[] fJe; delete[] fJbarE;
	fFlatConn = new int[fNelem * nen];
	fVrefE    = new double[fNelem];
	fJe       = new double[fNelem]();
	fJbarE    = new double[fNelem]();

	const dArray2DT& ref = ElementSupport().InitialCoordinates();
	for (int e = 0; e < fNelem; e++) {
		const ElementCardT& card = ElementCard(e);
		const iArrayT& nodes = card.NodesX();
		for (int n = 0; n < nen; n++)
			fFlatConn[e * nen + n] = nodes[n];

		double x[4], y[4], z[4];
		for (int n = 0; n < nen; n++) {
			int gn = fFlatConn[e * nen + n];
			x[n] = ref(gn, 0); y[n] = ref(gn, 1); z[n] = ref(gn, 2);
		}
		double a1=x[0]-x[2], a2=y[0]-y[2], a3=z[0]-z[2];
		double b1=x[1]-x[2], b2=y[1]-y[2], b3=z[1]-z[2];
		double c1=x[3]-x[2], c2=y[3]-y[2], c3=z[3]-z[2];
		fVrefE[e] = std::fabs(a1*(b2*c3-b3*c2)
		                    - a2*(b1*c3-b3*c1)
		                    + a3*(b1*c2-b2*c1)) / 6.0;
	}

	int numnod = ref.MajorDim();
	delete fANP;
	fANP = new ANPHelperT();
	fANP->Init(fNelem, numnod, nen, fFlatConn, fVrefE);
}

/* Compute J_e = V_curr / V_ref for every element. */
void BonetTetT::ComputeAllJe(void)
{
	const int nen = 4;
	const dArray2DT& cur = ElementSupport().CurrentCoordinates();
	#pragma omp parallel for if(fNelem > 1024)
	for (int e = 0; e < fNelem; e++) {
		const int* ec = fFlatConn + e * nen;
		double x[4], y[4], z[4];
		for (int n = 0; n < nen; n++) {
			x[n] = cur(ec[n], 0); y[n] = cur(ec[n], 1); z[n] = cur(ec[n], 2);
		}
		double a1=x[0]-x[2], a2=y[0]-y[2], a3=z[0]-z[2];
		double b1=x[1]-x[2], b2=y[1]-y[2], b3=z[1]-z[2];
		double c1=x[3]-x[2], c2=y[3]-y[2], c3=z[3]-z[2];
		double V = std::fabs(a1*(b2*c3-b3*c2)
		                   - a2*(b1*c3-b3*c1)
		                   + a3*(b1*c2-b2*c1)) / 6.0;
		fJe[e] = (fVrefE[e] > 0.0) ? V / fVrefE[e] : 1.0;
	}
}

/* Override SetGlobalShape — call parent first, then apply F-bar in-place. */
void BonetTetT::SetGlobalShape(void)
{
	UpdatedLagrangianT::SetGlobalShape();
	if (!fANP || NumElementNodes() != 4) return;

	int e = CurrElementNumber();
	if (e < 0 || e >= fNelem) return;

	double Jb = fJbarE[e];
	int nip = fF_List.Length();
	for (int ip = 0; ip < nip; ip++) {
		dMatrixT& F = fF_List[ip];
		double J = F.Det();
		if (J <= 0.0) continue;
		double scale = std::cbrt(Jb / J);
		double* p = F.Pointer();
		for (int k = 0; k < 9; k++) p[k] *= scale;
	}
}

/* Wrapper used by InitStep — and by TakeParameterList for the very first
 * residual evaluation, before InitStep has been called by the time loop. */
void BonetTetT::UpdateJBar(void)
{
	if (!fANP) return;
	ComputeAllJe();
	fANP->ComputeJBar(fJe, fJbarE);
}

/* Step-start hook: compute J̄_e from the current (start-of-step)
 * coordinates and cache it in fJbarE.  RHSDriver / LHSDriver read this
 * frozen value during every Newton iteration of the step, making the
 * FD tangent (which also reads fJbarE) consistent with the residual
 * — quadratic convergence on moderately near-incompressible cases (#29). */
void BonetTetT::InitStep(void)
{
	UpdatedLagrangianT::InitStep();
	UpdateJBar();
}

/* RHS / LHS: the lagged-J̄ scheme (#29) — J̄ is computed once at step
 * start by InitStep and held frozen for every Newton iteration of the
 * step.  Both residual and FD tangent see the same J̄, so the local
 * tangent is the consistent ∂R(u; J̄_step)/∂u and Newton converges
 * quadratically.  The trade-off is that the converged solution
 * satisfies R(u; J̄_step-start) = 0, not R(u; J̄(u)) = 0; the lag is
 * absorbed across time steps. */
void BonetTetT::RHSDriver(void)
{
	UpdatedLagrangianT::RHSDriver();
}

void BonetTetT::LHSDriver(GlobalT::SystemTypeT sys_type)
{
	UpdatedLagrangianT::LHSDriver(sys_type);
}

/* Refresh F + F-bar from the perturbed fLocDisp / fLocCurrCoords without
 * reloading them from the global Field.  Mirrors SimoQ1P0::
 * RecomputeDeformationState — SolidElementT::SetGlobalShape would call
 * SetLocalU(fLocDisp), clobbering the FD perturbation. */
void BonetTetT::RecomputeF(void)
{
	const int nip = NumIP();
	const int elem = CurrElementNumber();

	/* current-config shape function derivatives (read fLocCurrCoords) */
	fCurrShapes->SetDerivatives();

	/* F = I + grad u, computed via reference-config shapes */
	for (int ip = 0; ip < nip; ip++) {
		dMatrixT& F = fF_List[ip];
		fShapes->GradU(fLocDisp, F, ip);
		F.PlusIdentity();
	}

	/* F-bar scaling with the FROZEN J̄_e cached by InitStep */
	if (fANP && elem >= 0 && elem < fNelem) {
		double Jb = fJbarE[elem];
		for (int ip = 0; ip < nip; ip++) {
			dMatrixT& F = fF_List[ip];
			double J = F.Det();
			if (J <= 0.0) continue;
			double scale = std::cbrt(Jb / J);
			double* p = F.Pointer();
			for (int k = 0; k < 9; k++) p[k] *= scale;
		}
	}
}

/* Element internal force into a caller-provided buffer (writes, does not
 * accumulate).  Same physics as the inherited UpdatedLagrangianT::FormKd:
 * integrates B^T · σ over IPs, with σ from the F-bar in fF_List. */
void BonetTetT::ComputeInternalForce(double constK, dArrayT& force)
{
	const int n = fRHS.Length();
	dArrayT fRHS_save(n);
	for (int i = 0; i < n; i++) fRHS_save[i] = fRHS[i];

	fRHS = 0.0;
	UpdatedLagrangianT::FormKd(constK);
	for (int i = 0; i < n; i++) force[i] = fRHS[i];

	for (int i = 0; i < n; i++) fRHS[i] = fRHS_save[i];
}

/* Numerical (forward-FD) tangent of the F-bar residual at frozen J̄ —
 * the lagged-J̄ scheme's consistent tangent.  See header for the
 * derivation rationale.  Per-element, no cross-element terms. */
void BonetTetT::FormStiffness(double constK)
{
	if (!fANP) {                           /* fall back if not Tet4 */
		UpdatedLagrangianT::FormStiffness(constK);
		return;
	}

	const int nsd = NumSD();
	const int nen = NumElementNodes();
	const int ndof_elem = nen * nsd;
	const int nip = NumIP();
	const double eps = 1.0e-6;

	/* save state (the inner FD perturbs fLocDisp / fLocCurrCoords / fF_List) */
	dArrayT u_save(fLocDisp.Length());
	for (int k = 0; k < fLocDisp.Length(); k++) u_save[k] = fLocDisp[k];
	dArrayT x_save(fLocCurrCoords.Length());
	for (int k = 0; k < fLocCurrCoords.Length(); k++) x_save[k] = fLocCurrCoords[k];
	ArrayT<dMatrixT> F_save(nip);
	for (int ip = 0; ip < nip; ip++) {
		F_save[ip].Dimension(nsd, nsd);
		F_save[ip] = fF_List[ip];
	}

	/* baseline residual at the unperturbed state (uses cached F-bar in fF_List) */
	dArrayT f0(ndof_elem);
	ComputeInternalForce(constK, f0);

	/* column-by-column FD */
	dArrayT f_pert(ndof_elem);
	fLHS = 0.0;
	for (int a = 0; a < nen; a++) {
		for (int i = 0; i < nsd; i++) {
			fLocDisp(a, i)       += eps;
			fLocCurrCoords(a, i) += eps;
			RecomputeF();
			ComputeInternalForce(constK, f_pert);

			const int col = a * nsd + i;
			for (int j = 0; j < ndof_elem; j++)
				fLHS(j, col) += (f_pert[j] - f0[j]) / eps;

			fLocDisp(a, i)       -= eps;
			fLocCurrCoords(a, i) -= eps;
		}
	}

	/* restore state */
	for (int k = 0; k < fLocDisp.Length(); k++) fLocDisp[k] = u_save[k];
	for (int k = 0; k < fLocCurrCoords.Length(); k++) fLocCurrCoords[k] = x_save[k];
	for (int ip = 0; ip < nip; ip++) fF_List[ip] = F_save[ip];
	fCurrShapes->SetDerivatives();
}
