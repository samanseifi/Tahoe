/* $Id: PenaltyContact3DT.cpp,v 1.17 2011/12/01 21:11:36 bcyansfn Exp $ */
/* $Id: PenaltyContact3DT.cpp,v 2.0 2026/05/08 samanseifi Exp $ */
/* created: paklein (02/09/2000) */
#include "PenaltyContact3DT.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "ifstreamT.h"
#include "eIntegratorT.h"

using namespace Tahoe;

/* vector functions */
inline static void CrossProduct(const double* A, const double* B, double* AxB)
{   AxB[0] = A[1]*B[2] - A[2]*B[1];
	AxB[1] = A[2]*B[0] - A[0]*B[2];
	AxB[2] = A[0]*B[1] - A[1]*B[0];
};

inline static double Dot(const double* A, const double* B)
{ return A[0]*B[0] + A[1]*B[1] + A[2]*B[2]; };

inline static void Vector(const double* start, const double* end, double* v)
{
	v[0] = end[0] - start[0];
	v[1] = end[1] - start[1];
	v[2] = end[2] - start[2];
};

/* constructor */
PenaltyContact3DT::PenaltyContact3DT(const ElementSupportT& support):
	Contact3DT(support),
	fK(0.0),
	fMu(0.0),
	fFrictionEps(1.0e-6),
	fViscousDamping(0.0),
	fImplicitFriction(false)
{
	SetName("contact_3D_penalty");
}

/* describe the parameters needed by the interface */
void PenaltyContact3DT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	Contact3DT::DefineParameters(list);

	/* penalty stiffness */
	ParameterT stiffness(ParameterT::Double, "penalty_stiffness");
	stiffness.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(stiffness);

	/* Coulomb friction coefficient (0 = frictionless, default).
	 * When >0, regularised kinetic Coulomb is applied:
	 *   f_t = -mu * |f_n| * v_t / sqrt(|v_t|^2 + eps^2)  */
	ParameterT mu(ParameterT::Double, "friction_coefficient");
	mu.SetDefault(0.0);
	mu.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(mu, ParameterListT::ZeroOrOnce);

	ParameterT eps(ParameterT::Double, "friction_epsilon_velocity");
	eps.SetDefault(1.0e-6);
	eps.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(eps, ParameterListT::ZeroOrOnce);

	/* Normal-direction viscous damping (see #31).  Dimensional: force per
	 * unit velocity per unit area.  Default 0 = no damping (existing
	 * behaviour).  Typical value ~ 2 * sqrt(rho * E) for critical damping
	 * of the contact spring — user-calibrated. */
	ParameterT visc(ParameterT::Double, "viscous_damping");
	visc.SetDefault(0.0);
	visc.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(visc, ParameterListT::ZeroOrOnce);
}

/* accept parameter list */
void PenaltyContact3DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	Contact3DT::TakeParameterList(list);

	/* contact stiffness */
	fK = list.GetParameter("penalty_stiffness");

	/* Coulomb friction (defaults: mu=0 frictionless, eps=1e-6 regularisation).
	 * In the explicit branch eps is a velocity scale; in the implicit branch
	 * (#40) it is a slip-increment scale.  Same XML name to keep input files
	 * portable between integrators. */
	const ParameterT* p_mu = list.Parameter("friction_coefficient");
	if (p_mu) fMu = *p_mu;
	const ParameterT* p_eps = list.Parameter("friction_epsilon_velocity");
	if (p_eps) fFrictionEps = *p_eps;

	/* Detect integrator family by inspecting the field's kinematic order.
	 *   Order() == 0  →  no velocity stored  →  implicit (Newton-Raphson)
	 *   Order() >= 1  →  velocity stored     →  explicit (central-difference, Verlet, …) */
	fImplicitFriction = (Field().Order() < 1);
	if (fMu > 0.0)
		std::cout << "PenaltyContact3DT: Coulomb friction enabled, mu="
		          << fMu << ", eps=" << fFrictionEps
		          << " (" << (fImplicitFriction ? "implicit/slip-based" : "explicit/velocity-based")
		          << ", #" << (fImplicitFriction ? 40 : 26) << ")" << std::endl;

	/* Allocate per-striker slip history (issue #40).  Only used in the
	 * implicit branch but cheap to keep — a few KB on typical meshes. */
	if (fImplicitFriction && fMu > 0.0) {
		int num_strikers = fStrikerTags.Length();
		fPrevStrikerPos.Dimension(num_strikers, NumSD());
		fPrevFacetCentroid.Dimension(num_strikers, NumSD());
		fHasHistory.Dimension(num_strikers);
		fPrevStrikerPos    = 0.0;
		fPrevFacetCentroid = 0.0;
		fHasHistory        = 0;
	}

	/* Viscous normal damping (default 0 = no damping). */
	const ParameterT* p_visc = list.Parameter("viscous_damping");
	if (p_visc) fViscousDamping = *p_visc;
	if (fViscousDamping > 0.0)
		std::cout << "PenaltyContact3DT: viscous damping enabled, c="
		          << fViscousDamping << " (force per unit velocity per unit area)"
		          << std::endl;

	/* dimension workspace */
	fElCoord.Dimension(fNumFacetNodes + 1, NumSD());
	fElRefCoord.Dimension(fNumFacetNodes + 1, NumSD());
	fElDisp.Dimension(fNumFacetNodes + 1, NumDOF());

	fdc_du.Dimension(NumSD(), fElDisp.Length());
	fdn_du.Dimension(NumSD(), fElDisp.Length());
	fM1.Dimension(NumSD());
	fM2.Dimension(NumSD(), fElDisp.Length());
	fV1.Dimension(fElDisp.Length());

	double third = 1.0/3.0;
	double* p = fdc_du.Pointer();
	*p++ =-third;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ =-third;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ =-third;
	*p++ =-third;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ =-third;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ =-third;
	*p++ =-third;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ =-third;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ =-third;
	*p++ = 1;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ = 1;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p   = 1;
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* called by FormRHS and FormLHS */
void PenaltyContact3DT::LHSDriver(GlobalT::SystemTypeT)
{
	/* time integrator parameters */
	double constK = 0.0;
	int formK = fIntegrator->FormK(constK);
	if (!formK) return;

	/* references to global nodal data */
	const dArray2DT& init_coords = ElementSupport().InitialCoordinates();
	const dArray2DT& disp = Field()[0]; /* displacements */

	/* work space */
	dArrayT c(3), n(3);
	double a[3], b[3];
	iArrayT eqnos;
	
	/* loop over active elements */
	const iArray2DT& connects = *(fConnectivities[0]);
	for (int i = 0; i < connects.MajorDim(); i++)
	{
		const int* pelem = connects(i);

		/* collect element configuration */
		fElRefCoord.RowCollect(pelem, init_coords);
		fElDisp.RowCollect(pelem, disp);

		/* current configuration using effective displacement */
		fElCoord.SumOf(fElRefCoord, fElDisp);
	
		/* get facet and striker coords */
		fElCoord.RowAlias(0, fx1);
		fElCoord.RowAlias(1, fx2);
		fElCoord.RowAlias(2, fx3);
		fElCoord.RowAlias(3, fStriker);

		/* facet normal (direction) = a x b */
		Vector(fx1.Pointer(), fx2.Pointer(), a);
		Vector(fx1.Pointer(), fx3.Pointer(), b);
		CrossProduct(a, b, n.Pointer());
		double mag = sqrt(Dot(n.Pointer(), n.Pointer()));
		n[0] /= mag;
		n[1] /= mag;
		n[2] /= mag;

		/* (store) distance to facet */
		c[0] = fStriker[0] - (fx1[0] + fx2[0] + fx3[0])/3.0;
		c[1] = fStriker[1] - (fx1[1] + fx2[1] + fx3[1])/3.0;
		c[2] = fStriker[2] - (fx1[2] + fx2[2] + fx3[2])/3.0;
		double h = Dot(n.Pointer(), c.Pointer());

		/* contact */
		if (h < 0.0)
		{
			/* initialize */
			fRHS = 0.0;
			fLHS = 0.0;
	
			/* second variation of gap */
			DDg_tri_facet(
				fElRefCoord(0), fElRefCoord(1), fElRefCoord(2), fElRefCoord(3),
				fElDisp(0), fElDisp(1), fElDisp(2), fElDisp(3),
				fLHS);
			int striker_index = fStrikerTags_map.Map(pelem[3]);
			fLHS *= fK*h*constK*fStrikerArea[striker_index];

			/* d_c */
			fdc_du.MultTx(n, fV1);
			fRHS.SetToScaled(1.0, fV1);

			/* d_normal */
			Set_dn_du(fElCoord, fdn_du);					
			fM1.Outer(n, n);
			fM1.PlusIdentity(-1.0);
			fM2.MultATB(fM1, fdn_du);
			fM2.MultTx(c, fV1);
			fRHS.AddScaled(-1.0/mag, fV1);

			/* add term g^T g */
			fLHS.Outer(fRHS, fRHS, fK*constK*fStrikerArea[striker_index], dMatrixT::kAccumulate);

			/* get equation numbers */
			fEqnos[0].RowAlias(i, eqnos);
			
			/* assemble */
			ElementSupport().AssembleLHS(Group(), fLHS, eqnos);
		}
	}
}

void PenaltyContact3DT::RHSDriver(void)
{
	/* time integration parameters */
	double constKd = 0.0;
	int     formKd = fIntegrator->FormKd(constKd);
	if (!formKd) return;

	/* references to global nodal data */
	const dArray2DT& init_coords = ElementSupport().InitialCoordinates();
	const dArray2DT& disp = Field()[0]; /* displacements */

	/* reset tracking data */
	int num_contact = 0;
	double h_max = 0.0;

	/* clear force */
	fStrikerForce2D = 0.0;

	/* OpenMP-parallel contact loop (#33).  Workspaces (elCoord, dn_du, M1,
	 * M2, V1, RHS) are declared thread-local inside the loop so each thread
	 * has its own.  Class member fdc_du is constant after init — shared
	 * read-only.  fStrikerForce2D writes and AssembleRHS are serialised in
	 * a critical section (small per-pair cost, exact same scatter as serial).
	 *
	 * Threshold matches the explicit-element auto-tune (#32): only enable
	 * threading when the loop has enough work per thread to amortise the
	 * sync cost. */
	const int* base = fConnectivities[0]->Pointer();
	int rowlength = fConnectivities[0]->MinorDim();
	int N = fConnectivities[0]->MajorDim();

	int max_threads = 1;
#ifdef _OPENMP
	max_threads = omp_get_max_threads();
#endif
	bool use_omp = (N >= 32 * max_threads);

	#pragma omp parallel for schedule(dynamic, 16) \
		reduction(+:num_contact) reduction(min:h_max) if(use_omp)
	for (int i = 0; i < N; i++)
	{
		const int* pelem = base + i * rowlength;

		/* thread-local workspaces (replicate the class members) */
		dArray2DT elCoord(fNumFacetNodes + 1, NumSD());
		dArray2DT elDisp (fNumFacetNodes + 1, NumDOF());
		dMatrixT  dn_du  (NumSD(), elDisp.Length());
		dMatrixT  M1     (NumSD());
		dMatrixT  M2     (NumSD(), elDisp.Length());
		dArrayT   V1     (elDisp.Length());
		dArrayT   RHS    (elDisp.Length());
		dArrayT fx1, fx2, fx3, fStriker;
		dArrayT c(3), n(3);
		iArrayT eqnos;
		double a[3], b[3];

		/* collect element configuration */
		elCoord.RowCollect(pelem, init_coords);
		elDisp.RowCollect(pelem, disp);

		/* current configuration using effective displacement */
		elCoord.AddScaled(constKd, elDisp);

		/* get facet and striker coords */
		elCoord.RowAlias(0, fx1);
		elCoord.RowAlias(1, fx2);
		elCoord.RowAlias(2, fx3);
		elCoord.RowAlias(3, fStriker);

		/* facet normal (direction) = a x b */
		Vector(fx1.Pointer(), fx2.Pointer(), a);
		Vector(fx1.Pointer(), fx3.Pointer(), b);
		CrossProduct(a, b, n.Pointer());
		double mag = sqrt(Dot(n.Pointer(), n.Pointer()));
		n[0] /= mag;
		n[1] /= mag;
		n[2] /= mag;

		/* (store) distance to facet */
		c[0] = fStriker[0] - (fx1[0] + fx2[0] + fx3[0])/3.0;
		c[1] = fStriker[1] - (fx1[1] + fx2[1] + fx3[1])/3.0;
		c[2] = fStriker[2] - (fx1[2] + fx2[2] + fx3[2])/3.0;
		double h = Dot(n.Pointer(), c.Pointer());

		/* contact */
		if (h < 0.0)
		{
			/* tracking data — reductions handle accumulation */
			num_contact++;
			h_max = (h < h_max) ? h : h_max;

			/* penetration force */
			int striker_index = fStrikerTags_map.Map(pelem[3]);
			double dphi = -fK*h*fStrikerArea[striker_index];

			/* d_c (fdc_du is const member, read-only here) */
			fdc_du.MultTx(n, V1);
			RHS.SetToScaled(dphi, V1);

			/* d_normal */
			Set_dn_du(elCoord, dn_du);
			M1.Outer(n, n);
			M1.PlusIdentity(-1.0);
			M2.MultATB(M1, dn_du);
			M2.MultTx(c, V1);
			RHS.AddScaled(-dphi/mag, V1);

			/* ----- Coulomb friction (regularised kinetic) -----
			 * Two branches share the same Newton-3rd-law force distribution:
			 *   striker:    -f_t   (resists striker tangential motion)
			 *   facet 1..3: +f_t/3 (equal nodal split, preserves momentum)
			 *
			 * Explicit (#26): velocity-based, f_t = -mu·|f_n|·v_t / sqrt(|v_t|² + eps²)
			 * Implicit (#40): slip-based,    f_t = -mu·|f_n|·Δu_t / sqrt(|Δu_t|² + eps²)
			 *
			 * The implicit branch needs per-striker history for Δu_t; see
			 * fPrevStrikerPos / fPrevFacetCentroid / fHasHistory.  History
			 * is updated at CloseStep so the slip increment is over the
			 * just-converged load step.
			 *
			 * IMPLICIT TANGENT (NOT YET IN LHSDriver — known limitation #40):
			 * the residual contribution below is correct, but the analytical
			 * tangent ∂f_t/∂u of the slip-friction force is NOT added to
			 * LHSDriver.  Newton therefore relies only on the inherited normal
			 * contact stiffness; it converges only when the per-step slip
			 * increment stays small relative to friction_epsilon_velocity
			 * (i.e. small load steps), or when the tangential drag is small
			 * compared with the normal contact spring.  Production use of
			 * implicit Coulomb at large slip per step needs an analytical
			 * tangent  ∂f_t/∂Δu_t = s·(I − Δu_t⊗Δu_t/smooth²)  and  ∂f_t/∂|f_n|
			 * coupling — tracked as follow-up to #40. */
			if (fMu > 0.0 && !fImplicitFriction) {
				const dArray2DT& vel = Field()[1];
				double v_str[3] = {vel(pelem[3], 0), vel(pelem[3], 1), vel(pelem[3], 2)};
				double v_fc[3]  = {
					(vel(pelem[0], 0) + vel(pelem[1], 0) + vel(pelem[2], 0)) / 3.0,
					(vel(pelem[0], 1) + vel(pelem[1], 1) + vel(pelem[2], 1)) / 3.0,
					(vel(pelem[0], 2) + vel(pelem[1], 2) + vel(pelem[2], 2)) / 3.0
				};
				double v_rel[3] = {v_str[0]-v_fc[0], v_str[1]-v_fc[1], v_str[2]-v_fc[2]};
				double vrn = v_rel[0]*n[0] + v_rel[1]*n[1] + v_rel[2]*n[2];
				double vt[3] = {v_rel[0]-vrn*n[0], v_rel[1]-vrn*n[1], v_rel[2]-vrn*n[2]};
				double vt_mag2 = vt[0]*vt[0] + vt[1]*vt[1] + vt[2]*vt[2];
				double inv_smooth = 1.0 / sqrt(vt_mag2 + fFrictionEps*fFrictionEps);
				double f_n_mag = fabs(dphi);
				double s = -fMu * f_n_mag * inv_smooth;  /* opposes v_t */
				double ft[3] = {s*vt[0], s*vt[1], s*vt[2]};

				double third = 1.0 / 3.0;
				RHS[0]  += -ft[0] * third;  RHS[1]  += -ft[1] * third;  RHS[2]  += -ft[2] * third;
				RHS[3]  += -ft[0] * third;  RHS[4]  += -ft[1] * third;  RHS[5]  += -ft[2] * third;
				RHS[6]  += -ft[0] * third;  RHS[7]  += -ft[1] * third;  RHS[8]  += -ft[2] * third;
				RHS[9]  += ft[0];           RHS[10] += ft[1];           RHS[11] += ft[2];
			}
			else if (fMu > 0.0 && fImplicitFriction && fHasHistory[striker_index]) {
				/* slip increment Δr = (current relative position) − (previous relative position) */
				double r_cur[3] = {fStriker[0] - (fx1[0]+fx2[0]+fx3[0])/3.0,
				                   fStriker[1] - (fx1[1]+fx2[1]+fx3[1])/3.0,
				                   fStriker[2] - (fx1[2]+fx2[2]+fx3[2])/3.0};
				double r_prev[3] = {fPrevStrikerPos(striker_index, 0) - fPrevFacetCentroid(striker_index, 0),
				                    fPrevStrikerPos(striker_index, 1) - fPrevFacetCentroid(striker_index, 1),
				                    fPrevStrikerPos(striker_index, 2) - fPrevFacetCentroid(striker_index, 2)};
				double dr[3] = {r_cur[0]-r_prev[0], r_cur[1]-r_prev[1], r_cur[2]-r_prev[2]};
				double drn = dr[0]*n[0] + dr[1]*n[1] + dr[2]*n[2];
				double dut[3] = {dr[0]-drn*n[0], dr[1]-drn*n[1], dr[2]-drn*n[2]};
				double dut_mag2 = dut[0]*dut[0] + dut[1]*dut[1] + dut[2]*dut[2];
				double inv_smooth = 1.0 / sqrt(dut_mag2 + fFrictionEps*fFrictionEps);
				double f_n_mag = fabs(dphi);
				double s = -fMu * f_n_mag * inv_smooth;
				double ft[3] = {s*dut[0], s*dut[1], s*dut[2]};

				double third = 1.0 / 3.0;
				RHS[0]  += -ft[0] * third;  RHS[1]  += -ft[1] * third;  RHS[2]  += -ft[2] * third;
				RHS[3]  += -ft[0] * third;  RHS[4]  += -ft[1] * third;  RHS[5]  += -ft[2] * third;
				RHS[6]  += -ft[0] * third;  RHS[7]  += -ft[1] * third;  RHS[8]  += -ft[2] * third;
				RHS[9]  += ft[0];           RHS[10] += ft[1];           RHS[11] += ft[2];
			}

			/* ----- Normal viscous damping (#31) -----
			 * f_visc_n = -c * v_n_relative * area
			 * Damps the spring-mass ringing of the penalty contact.
			 *
			 * Distribution: same Newton-3rd-law split as Coulomb friction —
			 *   striker:  +f_visc * n (away from facet, resists approach)
			 *   facet 1..3: -f_visc * n / 3
			 *
			 * Note: a striker that's APPROACHING (v_n_rel < 0) gets a
			 * REPULSIVE damp force (along +n), reinforcing the normal spring;
			 * a striker that's SEPARATING (v_n_rel > 0) gets an ATTRACTIVE
			 * damp force.  This mimics LS-DYNA VDC behaviour.
			 */
			if (fViscousDamping > 0.0) {
				const dArray2DT& vel = Field()[1];
				double vsr[3] = {vel(pelem[3], 0), vel(pelem[3], 1), vel(pelem[3], 2)};
				double vfc[3] = {
					(vel(pelem[0], 0) + vel(pelem[1], 0) + vel(pelem[2], 0)) / 3.0,
					(vel(pelem[0], 1) + vel(pelem[1], 1) + vel(pelem[2], 1)) / 3.0,
					(vel(pelem[0], 2) + vel(pelem[1], 2) + vel(pelem[2], 2)) / 3.0
				};
				double vrn = (vsr[0]-vfc[0])*n[0] + (vsr[1]-vfc[1])*n[1] + (vsr[2]-vfc[2])*n[2];
				double area = fStrikerArea[striker_index];
				double f_visc = -fViscousDamping * vrn * area;
				double third = 1.0 / 3.0;
				RHS[0]  += -f_visc * n[0] * third;
				RHS[1]  += -f_visc * n[1] * third;
				RHS[2]  += -f_visc * n[2] * third;
				RHS[3]  += -f_visc * n[0] * third;
				RHS[4]  += -f_visc * n[1] * third;
				RHS[5]  += -f_visc * n[2] * third;
				RHS[6]  += -f_visc * n[0] * third;
				RHS[7]  += -f_visc * n[1] * third;
				RHS[8]  += -f_visc * n[2] * third;
				RHS[9]  +=  f_visc * n[0];
				RHS[10] +=  f_visc * n[1];
				RHS[11] +=  f_visc * n[2];
			}

			/* get equation numbers (read-only — fEqnos is shared) */
			fEqnos[0].RowAlias(i, eqnos);

			/* serialise the global-RHS scatter and per-striker force store
			 * so multiple threads contributing to the same equation row or
			 * the same striker index don't race. */
			#pragma omp critical(PenaltyContact3DT_assemble)
			{
				ElementSupport().AssembleRHS(Group(), RHS, eqnos);
				fStrikerForce2D(striker_index, 0) = dphi*n[0];
				fStrikerForce2D(striker_index, 1) = dphi*n[1];
				fStrikerForce2D(striker_index, 2) = dphi*n[2];
			}
		}
	}

	/* set tracking */
	SetTrackingData(num_contact, h_max);
}

/* Step-end hook: snapshot the contact-pair geometry per striker so the next
 * step's RHSDriver can compute the tangential slip increment Δu_t for
 * implicit Coulomb friction (#40).  We loop over the same pair connectivity
 * the residual sees, in serial — cost is O(num_pairs) and only one pass per
 * accepted load step. */
void PenaltyContact3DT::CloseStep(void)
{
	Contact3DT::CloseStep();

	if (!fImplicitFriction || fMu <= 0.0) return;

	const dArray2DT& init_coords = ElementSupport().InitialCoordinates();
	const dArray2DT& disp        = Field()[0];

	const int* base = fConnectivities[0]->Pointer();
	int rowlength   = fConnectivities[0]->MinorDim();
	int N           = fConnectivities[0]->MajorDim();

	for (int i = 0; i < N; i++) {
		const int* pelem = base + i * rowlength;
		int s_idx = fStrikerTags_map.Map(pelem[3]);

		double s_pos[3] = {
			init_coords(pelem[3], 0) + disp(pelem[3], 0),
			init_coords(pelem[3], 1) + disp(pelem[3], 1),
			init_coords(pelem[3], 2) + disp(pelem[3], 2)
		};
		double f_cent[3] = {0.0, 0.0, 0.0};
		for (int k = 0; k < 3; k++)
			for (int j = 0; j < 3; j++)
				f_cent[j] += (init_coords(pelem[k], j) + disp(pelem[k], j)) / 3.0;

		fPrevStrikerPos(s_idx, 0)    = s_pos[0];
		fPrevStrikerPos(s_idx, 1)    = s_pos[1];
		fPrevStrikerPos(s_idx, 2)    = s_pos[2];
		fPrevFacetCentroid(s_idx, 0) = f_cent[0];
		fPrevFacetCentroid(s_idx, 1) = f_cent[1];
		fPrevFacetCentroid(s_idx, 2) = f_cent[2];
		fHasHistory[s_idx]           = 1;
	}
}
