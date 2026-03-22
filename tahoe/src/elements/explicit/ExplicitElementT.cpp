/* ExplicitElementT.cpp — MVSIZ-batched explicit solid element. */
#include "ExplicitElementT.h"

#include "UpdatedLagrangianT.h"
#include "ExplicitKernelT.h"
#include "Q4KernelT.h"
#include "Hex8KernelT.h"
#include "ExplicitMaterialT.h"
#include "ExplNeoHookeanT.h"

#include "ElementSupportT.h"
#include "ModelManagerT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"
#include "eIntegratorT.h"
#include "MaterialListT.h"
#include "ContinuumMaterialT.h"

#include <iostream>
#include <cstring>
#include <chrono>

using namespace Tahoe;

ExplicitElementT::ExplicitElementT(const ElementSupportT& support)
	: UpdatedLagrangianT(support),
	  fKernel(NULL),
	  fBatchMaterial(NULL),
	  fTotalElements(0),
	  fFlatConn(NULL),
	  fFlatEqnos(NULL),
	  fGlobalRHS(NULL),
	  fHourglassType(kNoHourglass),
	  fHourglassCoeff(0.1),
	  fMassScalingType(kNoMassScaling),
	  fTargetDt(0.0),
	  fDtScaleFactor(0.9),
	  fMassScaleInterval(100),
	  fMassScale(NULL)
{
	SetName("explicit_solid");
}

ExplicitElementT::~ExplicitElementT(void)
{
	delete fKernel;
	delete fBatchMaterial;
	delete[] fFlatConn;
	delete[] fFlatEqnos;
	delete[] fMassScale;
}

/*----------------------------------------------------------------------
 * ParameterInterfaceT
 *----------------------------------------------------------------------*/
void ExplicitElementT::DefineParameters(ParameterListT& list) const
{
	/* inherited — picks up field_name, mass_type, etc. */
	UpdatedLagrangianT::DefineParameters(list);
}

void ExplicitElementT::DefineSubs(SubListT& sub_list) const
{
	UpdatedLagrangianT::DefineSubs(sub_list);
	sub_list.AddSub("hourglass_control", ParameterListT::ZeroOrOnce);
	sub_list.AddSub("mass_scaling", ParameterListT::ZeroOrOnce);
}

ParameterInterfaceT* ExplicitElementT::NewSub(const StringT& name) const
{
	if (name == "hourglass_control") {
		ParameterContainerT* hg = new ParameterContainerT(name);

		ParameterT type(ParameterT::Enumeration, "type");
		type.AddEnumeration("viscous", kViscousHG);
		type.AddEnumeration("stiffness", kStiffnessHG);
		type.SetDefault(kViscousHG);
		hg->AddParameter(type);

		ParameterT coeff(ParameterT::Double, "coefficient");
		coeff.SetDefault(0.1);
		hg->AddParameter(coeff);

		return hg;
	}
	if (name == "mass_scaling") {
		ParameterContainerT* ms = new ParameterContainerT(name);

		ParameterT type(ParameterT::Enumeration, "type");
		type.AddEnumeration("fixed", kFixedMassScaling);
		type.AddEnumeration("adaptive", kAdaptiveMassScaling);
		type.SetDefault(kFixedMassScaling);
		ms->AddParameter(type);

		ParameterT target_dt(ParameterT::Double, "target_dt");
		target_dt.SetDefault(0.0);
		ms->AddParameter(target_dt);

		ParameterT scale_factor(ParameterT::Double, "scale_factor");
		scale_factor.SetDefault(0.9);
		ms->AddParameter(scale_factor);

		ParameterT interval(ParameterT::Integer, "update_interval");
		interval.SetDefault(100);
		ms->AddParameter(interval);

		return ms;
	}
	return UpdatedLagrangianT::NewSub(name);
}

void ExplicitElementT::TakeParameterList(const ParameterListT& list)
{
	/* inherited — sets up connectivity, shape functions, materials, etc.
	 * Uses the same XML format as updated_lagrangian. */
	UpdatedLagrangianT::TakeParameterList(list);

	/* create kernel based on geometry */
	int nsd = NumSD();
	int nen = NumElementNodes();
	if (nsd == 2 && nen == 4)
		fKernel = new Q4KernelT;
	else if (nsd == 3 && nen == 8)
		fKernel = new Hex8KernelT;
	else {
		std::cout << "ExplicitElementT: unsupported topology nsd="
		          << nsd << " nen=" << nen
		          << " — using legacy element loop" << std::endl;
		fKernel = NULL;
	}

	/* create batch material from the first block's XML parameters.
	 * Scan recursively for mu/kappa/density in the material sub-tree.
	 * Works with both direct neo-hookean and RG_split_general wrappers. */
	int num_blocks = list.NumLists("large_strain_element_block");
	if (num_blocks > 0) {
		const ParameterListT& block = list.GetList("large_strain_element_block");
		const ParameterListT& mat_choice = block.GetListChoice(*this, "large_strain_material_choice");

		/* recursive search for mu/kappa in the material parameter tree */
		double mu = 0.0, kappa = 0.0, density = 1.0;
		bool found_mu = false, found_kappa = false;

		/* search function: scan this list and all its sub-lists */
		struct ParamFinder {
			static void Find(const ParameterListT& p,
				double& mu, double& kappa, double& density,
				bool& found_mu, bool& found_kappa)
			{
				const ParameterT* pm = p.Parameter("mu");
				const ParameterT* pk = p.Parameter("kappa");
				const ParameterT* pd = p.Parameter("density");
				if (pm) { mu = *pm; found_mu = true; }
				if (pk) { kappa = *pk; found_kappa = true; }
				if (pd) { density = *pd; }
				/* recurse into sub-lists */
				const ArrayT<ParameterListT>& subs = p.Lists();
				for (int i = 0; i < subs.Length(); i++)
					Find(subs[i], mu, kappa, density, found_mu, found_kappa);
			}
		};

		ParamFinder::Find(mat_choice, mu, kappa, density, found_mu, found_kappa);

		if (found_mu && found_kappa)
			fBatchMaterial = new ExplNeoHookeanT(mu, kappa, density);
	}

	/* hourglass control */
	if (list.NumLists("hourglass_control") > 0) {
		const ParameterListT& hg = list.GetList("hourglass_control");
		int hg_type = hg.GetParameter("type");
		fHourglassType = (HourglassTypeT)hg_type;
		fHourglassCoeff = hg.GetParameter("coefficient");
	}

	/* mass scaling */
	if (list.NumLists("mass_scaling") > 0) {
		const ParameterListT& ms = list.GetList("mass_scaling");
		int ms_type = ms.GetParameter("type");
		fMassScalingType = (MassScalingTypeT)ms_type;
		fTargetDt = ms.GetParameter("target_dt");
		fDtScaleFactor = ms.GetParameter("scale_factor");
		fMassScaleInterval = ms.GetParameter("update_interval");
	}

	if (fKernel && fBatchMaterial) {
		BuildFlatArrays();
		double dt_cfl = ComputeStableTimeStep();

		/* apply initial mass scaling if requested */
		if (fMassScalingType != kNoMassScaling) {
			fMassScale = new double[fTotalElements];
			for (int e = 0; e < fTotalElements; e++)
				fMassScale[e] = 1.0;
			if (fTargetDt > 0.0)
				ApplyMassScaling();
		}
		std::cout << "ExplicitElementT: MVSIZ=" << MVSIZ
		          << " batched path active (nsd=" << nsd
		          << " nen=" << nen
		          << " nel=" << fTotalElements
		          << " CFL_dt=" << dt_cfl;
		if (fHourglassType != kNoHourglass)
			std::cout << " hourglass="
			          << (fHourglassType == kViscousHG ? "viscous" : "stiffness")
			          << " coeff=" << fHourglassCoeff;
		if (fMassScalingType != kNoMassScaling) {
			std::cout << " mass_scaling="
			          << (fMassScalingType == kFixedMassScaling ? "fixed" : "adaptive")
			          << " target_dt=" << fTargetDt;
			if (fMassScalingType == kAdaptiveMassScaling)
				std::cout << " interval=" << fMassScaleInterval;
			/* count how many elements were scaled */
			if (fMassScale) {
				int n_scaled = 0;
				double max_scale = 1.0;
				for (int e = 0; e < fTotalElements; e++) {
					if (fMassScale[e] > 1.001) n_scaled++;
					if (fMassScale[e] > max_scale) max_scale = fMassScale[e];
				}
				std::cout << " scaled=" << n_scaled << "/" << fTotalElements
				          << " max_factor=" << max_scale;
			}
		}
		std::cout << ")" << std::endl;
	} else
		std::cout << "ExplicitElementT: falling back to legacy element loop" << std::endl;
}

/*----------------------------------------------------------------------
 * BuildFlatArrays — pre-compute connectivity and equation numbers
 *----------------------------------------------------------------------*/
void ExplicitElementT::BuildFlatArrays(void)
{
	const int nen = fKernel->NodesPerElement();
	const int ndof = NumDOF();

	/* count total elements */
	fTotalElements = 0;
	for (int b = 0; b < fBlockData.Length(); b++)
		fTotalElements += fBlockData[b].Dimension();

	/* allocate flat arrays */
	fFlatConn = new int[fTotalElements * nen];
	fFlatEqnos = new int[fTotalElements * nen * ndof];

	/* fill connectivity from ElementCards */
	for (int e = 0; e < fTotalElements; e++) {
		const ElementCardT& card = ElementCard(e);
		const iArrayT& nodes = card.NodesX();
		for (int n = 0; n < nen; n++)
			fFlatConn[e * nen + n] = nodes[n];
	}

	/* equation numbers are set up after TakeParameterList when the
	 * field equations are initialized. We'll fill them lazily on
	 * first call to BatchedInternalForce. */
}

/*----------------------------------------------------------------------
 * ComputeStableTimeStep — CFL condition for all elements
 *
 * For each element: dt = h / c  where
 *   h = characteristic length = sqrt(area) for 2D, cbrt(vol) for 3D
 *   c = wave speed = sqrt((kappa + 4mu/3) / rho)
 * Returns the global minimum dt.
 *----------------------------------------------------------------------*/
double ExplicitElementT::ComputeStableTimeStep(void) const
{
	if (!fKernel || !fBatchMaterial) return 1.0e30;

	const int nsd = NumSD();
	const int nen = fKernel->NodesPerElement();
	const dArray2DT& ref_coords = ElementSupport().InitialCoordinates();
	double rho = fBatchMaterial->Density();

	/* approximate P-wave speed: need kappa and mu from the material.
	 * For now, compute from the reference configuration Jacobian. */
	/* TODO: get kappa and mu from the material interface */
	double dt_min = 1.0e30;

	for (int e = 0; e < fTotalElements; e++) {
		const int* ec = fFlatConn + e * nen;

		if (nsd == 2) {
			/* compute approximate element area from node coordinates */
			double x[4], y[4];
			for (int n = 0; n < nen; n++) {
				x[n] = ref_coords(ec[n], 0);
				y[n] = ref_coords(ec[n], 1);
			}
			/* shoelace formula for Q4 area */
			double area = 0.5 * fabs(
				(x[0]-x[2])*(y[1]-y[3]) - (x[1]-x[3])*(y[0]-y[2]));
			double h = sqrt(area);

			/* dt = h / c where c = sqrt((kappa + 4mu/3) / rho)
			 * For now estimate c from element size and density */
			/* We don't have direct kappa/mu access from material interface.
			 * Use a conservative estimate based on the stored material. */
			/* TODO: add wave_speed() method to ExplicitMaterialT */
			double c = fBatchMaterial->WaveSpeed();
			double dt_elem = h / c;
			if (dt_elem < dt_min) dt_min = dt_elem;
		}
		else {
			/* 3D: approximate volume from hex coords */
			double x[8], y[8], z[8];
			for (int n = 0; n < nen; n++) {
				x[n] = ref_coords(ec[n], 0);
				y[n] = ref_coords(ec[n], 1);
				z[n] = ref_coords(ec[n], 2);
			}
			/* approximate volume: 1/6 * |diag1 x diag2 . diag3| */
			double dx1 = x[6]-x[0], dy1 = y[6]-y[0], dz1 = z[6]-z[0];
			double dx2 = x[7]-x[1], dy2 = y[7]-y[1], dz2 = z[7]-z[1];
			double dx3 = x[5]-x[3], dy3 = y[5]-y[3], dz3 = z[5]-z[3];
			double vol = fabs(dx1*(dy2*dz3-dz2*dy3)
			                - dy1*(dx2*dz3-dz2*dx3)
			                + dz1*(dx2*dy3-dy2*dx3)) / 6.0;
			double h = cbrt(vol);
			double c = 1.0; /* placeholder */
			double dt_elem = h / c;
			if (dt_elem < dt_min) dt_min = dt_elem;
		}
	}

	return dt_min;
}

/*----------------------------------------------------------------------
 * ApplyMassScaling — increase element mass to meet target time step
 *
 * For each element: dt_elem = h / c.
 * If dt_elem < target_dt, scale mass by (target_dt / dt_elem)^2.
 * Since dt ~ sqrt(M/K), scaling M by alpha^2 scales dt by alpha.
 *
 * The mass scale factors are stored per element. The actual mass matrix
 * is modified by the parent class mass assembly — we modify the density
 * seen by FormMass. For now, we report the scaling but the actual mass
 * modification requires hooking into the LHS assembly (future work).
 *----------------------------------------------------------------------*/
void ExplicitElementT::ApplyMassScaling(void)
{
	if (!fKernel || !fBatchMaterial || !fMassScale) return;

	const int nsd = NumSD();
	const int nen = fKernel->NodesPerElement();
	const dArray2DT& ref_coords = ElementSupport().InitialCoordinates();
	double c = fBatchMaterial->WaveSpeed();
	double target = fTargetDt * fDtScaleFactor; /* apply safety factor */

	int n_scaled = 0;
	double max_scale = 1.0;

	for (int e = 0; e < fTotalElements; e++) {
		const int* ec = fFlatConn + e * nen;
		double h;

		if (nsd == 2) {
			double x[4], y[4];
			for (int n = 0; n < nen; n++) {
				x[n] = ref_coords(ec[n], 0);
				y[n] = ref_coords(ec[n], 1);
			}
			double area = 0.5 * fabs(
				(x[0]-x[2])*(y[1]-y[3]) - (x[1]-x[3])*(y[0]-y[2]));
			h = sqrt(area);
		} else {
			double x[8], y[8], z[8];
			for (int n = 0; n < nen; n++) {
				x[n] = ref_coords(ec[n], 0);
				y[n] = ref_coords(ec[n], 1);
				z[n] = ref_coords(ec[n], 2);
			}
			double dx1=x[6]-x[0], dy1=y[6]-y[0], dz1=z[6]-z[0];
			double dx2=x[7]-x[1], dy2=y[7]-y[1], dz2=z[7]-z[1];
			double dx3=x[5]-x[3], dy3=y[5]-y[3], dz3=z[5]-z[3];
			double vol = fabs(dx1*(dy2*dz3-dz2*dy3)
			                - dy1*(dx2*dz3-dz2*dx3)
			                + dz1*(dx2*dy3-dy2*dx3)) / 6.0;
			h = cbrt(vol);
		}

		double dt_elem = h / c;
		if (dt_elem < target) {
			double alpha = target / dt_elem;
			fMassScale[e] = alpha * alpha; /* M_new = alpha^2 * M_old */
			n_scaled++;
			if (fMassScale[e] > max_scale) max_scale = fMassScale[e];
		} else {
			fMassScale[e] = 1.0;
		}
	}

	if (n_scaled > 0)
		std::cout << "ExplicitElementT: mass scaling applied to "
		          << n_scaled << "/" << fTotalElements
		          << " elements (max factor=" << max_scale << ")" << std::endl;
}

/*----------------------------------------------------------------------
 * RHSDriver — override the virtual dispatch point
 *----------------------------------------------------------------------*/
void ExplicitElementT::RHSDriver(void)
{
	/* inherited: traction BCs, body forces from ContinuumElementT */
	ContinuumElementT::RHSDriver();

	/* adaptive mass scaling: periodically recompute scale factors */
	if (fMassScalingType == kAdaptiveMassScaling && fMassScale) {
		static int step_count = 0;
		step_count++;
		if (step_count % fMassScaleInterval == 0)
			ApplyMassScaling();
	}

	/* fast path: batched internal force */
	if (fKernel && fBatchMaterial) {
		double constKd = 0.0;
		int formKd = fIntegrator->FormKd(constKd);
		if (formKd) {
			BatchedInternalForce(-constKd);
			return;
		}
	}

	/* fallback: use the generic per-element loop */
	SolidElementT::RHSDriver();
}

/*----------------------------------------------------------------------
 * BatchedInternalForce — the MVSIZ SoA loop
 *
 * Single-pass: computes both reference Jacobian (for F = dx/dX) and
 * current Jacobian (for dN/dx and detJ) inline, avoiding the double
 * kernel call overhead.
 *----------------------------------------------------------------------*/
void ExplicitElementT::BatchedInternalForce(double constKd)
{
	const int nsd = NumSD();
	const int nen = fKernel->NodesPerElement();
	const int nip = fKernel->NumIP();
	const int ndof = NumDOF();

	/* global coordinate arrays */
	const dArray2DT& init_coords = ElementSupport().InitialCoordinates();
	const dArray2DT& curr_coords = ElementSupport().CurrentCoordinates();
	int numnod = init_coords.MajorDim();

	/* timing instrumentation */
	static long long t_total = 0, t_assemble = 0;
	static int call_count = 0;
	call_count++;
	auto t0 = std::chrono::high_resolution_clock::now();

	/* global force accumulation — persistent array, zeroed each step */
	dArray2DT force(numnod, ndof);
	force = 0.0;

	/* raw coordinate pointers for direct indexing */
	const double* cur_ptr = curr_coords.Pointer();
	const double* ref_ptr = init_coords.Pointer();
	int coord_stride = ndof; /* nsd == ndof for displacement field */

	/* total batches across all blocks */
	int num_batches = (fTotalElements + MVSIZ - 1) / MVSIZ;

	/* OpenMP parallel loop over batches */
	#pragma omp parallel for schedule(dynamic, 1) if(num_batches > 1)
	for (int ibatch = 0; ibatch < num_batches; ibatch++) {
		int batch_start = ibatch * MVSIZ;
		int nel = batch_start + MVSIZ <= fTotalElements
		        ? MVSIZ : fTotalElements - batch_start;

		/* per-thread SoA workspace (stack allocated) */
		double xc[ExplicitKernelT::MAX_NEN][MVSIZ];
		double yc[ExplicitKernelT::MAX_NEN][MVSIZ];
		double zc[ExplicitKernelT::MAX_NEN][MVSIZ];
		double xr[ExplicitKernelT::MAX_NEN][MVSIZ];
		double yr[ExplicitKernelT::MAX_NEN][MVSIZ];
		double zr[ExplicitKernelT::MAX_NEN][MVSIZ];
		double dNdX[ExplicitKernelT::MAX_NEN][MVSIZ];
		double dNdY[ExplicitKernelT::MAX_NEN][MVSIZ];
		double dNdZ[ExplicitKernelT::MAX_NEN][MVSIZ];
		double detJ_cur[MVSIZ];
		double F2D[4][MVSIZ];
		double F3D[9][MVSIZ];
		double sig11[MVSIZ], sig22[MVSIZ], sig12[MVSIZ];
		double sig3D[6][MVSIZ];
		double fx[ExplicitKernelT::MAX_NEN][MVSIZ];
		double fy[ExplicitKernelT::MAX_NEN][MVSIZ];
		double fz[ExplicitKernelT::MAX_NEN][MVSIZ];
		int lconn[ExplicitKernelT::MAX_NEN][MVSIZ];

		/* 1. GATHER — flat connectivity, direct pointer indexing */
		const int* conn_base = fFlatConn + batch_start * nen;
		for (int i = 0; i < nel; i++) {
			const int* ec = conn_base + i * nen;
			for (int n = 0; n < nen; n++) {
				int node = ec[n];
				lconn[n][i] = node;
				xc[n][i] = cur_ptr[node * coord_stride + 0];
				yc[n][i] = cur_ptr[node * coord_stride + 1];
				xr[n][i] = ref_ptr[node * coord_stride + 0];
				yr[n][i] = ref_ptr[node * coord_stride + 1];
			}
			if (nsd == 3) {
				for (int n = 0; n < nen; n++) {
					int node = ec[n];
					zc[n][i] = cur_ptr[node * coord_stride + 2];
					zr[n][i] = ref_ptr[node * coord_stride + 2];
				}
			}
			for (int n = 0; n < nen; n++) {
				fx[n][i] = 0.0;
				fy[n][i] = 0.0;
			}
			if (nsd == 3)
				for (int n = 0; n < nen; n++)
					fz[n][i] = 0.0;
		}

			/* 2-5. SINGLE-PASS IP LOOP: one kernel call on reference
			 * coords gives dN/dX. Then compute F and detJ_cur inline. */
			for (int ip = 0; ip < nip; ip++) {
				double w_ref;

				/* reference config: dN/dX for F computation */
				fKernel->ComputeIPData(ip, nel,
					xr, yr, zr, dNdX, dNdY, dNdZ, detJ_cur, w_ref);
				/* Note: detJ_cur from ref coords = detJ_ref.
				 * We'll compute detJ_cur from F below. */

				if (nsd == 2) {
					/* Compute F = dx/dX and detJ_cur = det(F)*detJ_ref
					 * in a single vectorized loop */
					for (int i = 0; i < nel; i++) {
						double f11 = 0.0, f12 = 0.0, f21 = 0.0, f22 = 0.0;
						for (int n = 0; n < nen; n++) {
							f11 += xc[n][i] * dNdX[n][i];
							f12 += xc[n][i] * dNdY[n][i];
							f21 += yc[n][i] * dNdX[n][i];
							f22 += yc[n][i] * dNdY[n][i];
						}
						F2D[0][i] = f11; F2D[1][i] = f12;
						F2D[2][i] = f21; F2D[3][i] = f22;

						/* detJ_cur = det(F) * detJ_ref
						 * (since dx = F * dX, the current Jacobian = F * J_ref) */
						double detF = f11*f22 - f12*f21;
						detJ_cur[i] = detF * detJ_cur[i]; /* detJ_cur was detJ_ref */

						/* dN/dx = F^{-T} * dN/dX (push forward of shape derivs)
						 * F^{-1} = (1/detF) * [F22, -F12; -F21, F11]
						 * Reuse dNdX/dNdY arrays in-place for dN/dx */
					}

					/* Compute dN/dx from dN/dX via F^{-T} push-forward
					 * (separate loop for better vectorization) */
					for (int i = 0; i < nel; i++) {
						double f11 = F2D[0][i], f12 = F2D[1][i];
						double f21 = F2D[2][i], f22 = F2D[3][i];
						double invDetF = 1.0 / (f11*f22 - f12*f21);
						/* F^{-T} = (1/detF) * [F22, -F21; -F12, F11] */
						double FiT11 =  f22*invDetF, FiT12 = -f21*invDetF;
						double FiT21 = -f12*invDetF, FiT22 =  f11*invDetF;
						for (int n = 0; n < nen; n++) {
							double dX = dNdX[n][i], dY = dNdY[n][i];
							dNdX[n][i] = FiT11*dX + FiT12*dY; /* now dN/dx */
							dNdY[n][i] = FiT21*dX + FiT22*dY; /* now dN/dy */
						}
					}

					/* MATERIAL */
					fBatchMaterial->ComputeStress2D(nel, F2D,
						sig11, sig22, sig12, NULL);

					/* B^T * sigma (uses dNdX/dNdY which are now dN/dx) */
					for (int i = 0; i < nel; i++) {
						double scale = constKd * detJ_cur[i] * w_ref;
						for (int n = 0; n < nen; n++) {
							fx[n][i] += (dNdX[n][i]*sig11[i] + dNdY[n][i]*sig12[i]) * scale;
							fy[n][i] += (dNdX[n][i]*sig12[i] + dNdY[n][i]*sig22[i]) * scale;
						}
					}
				}
				else { /* nsd == 3 */
					/* F = dx/dX */
					for (int i = 0; i < nel; i++) {
						double f[9] = {0.0};
						for (int n = 0; n < nen; n++) {
							f[0] += xc[n][i]*dNdX[n][i]; f[1] += xc[n][i]*dNdY[n][i]; f[2] += xc[n][i]*dNdZ[n][i];
							f[3] += yc[n][i]*dNdX[n][i]; f[4] += yc[n][i]*dNdY[n][i]; f[5] += yc[n][i]*dNdZ[n][i];
							f[6] += zc[n][i]*dNdX[n][i]; f[7] += zc[n][i]*dNdY[n][i]; f[8] += zc[n][i]*dNdZ[n][i];
						}
						for (int k = 0; k < 9; k++) F3D[k][i] = f[k];

						/* detJ_cur = det(F) * detJ_ref */
						double detF = f[0]*(f[4]*f[8]-f[5]*f[7])
						            - f[1]*(f[3]*f[8]-f[5]*f[6])
						            + f[2]*(f[3]*f[7]-f[4]*f[6]);
						detJ_cur[i] = detF * detJ_cur[i];
					}

					/* dN/dx = F^{-T} * dN/dX (3D push-forward) */
					for (int i = 0; i < nel; i++) {
						double *f = &F3D[0][i]; /* stride = MVSIZ */
						double f0=F3D[0][i],f1=F3D[1][i],f2=F3D[2][i];
						double f3=F3D[3][i],f4=F3D[4][i],f5=F3D[5][i];
						double f6=F3D[6][i],f7=F3D[7][i],f8=F3D[8][i];
						double detF = f0*(f4*f8-f5*f7) - f1*(f3*f8-f5*f6) + f2*(f3*f7-f4*f6);
						double id = 1.0/detF;
						/* F^{-T} cofactors / detF */
						double c00=(f4*f8-f5*f7)*id, c01=(f5*f6-f3*f8)*id, c02=(f3*f7-f4*f6)*id;
						double c10=(f2*f7-f1*f8)*id, c11=(f0*f8-f2*f6)*id, c12=(f1*f6-f0*f7)*id;
						double c20=(f1*f5-f2*f4)*id, c21=(f2*f3-f0*f5)*id, c22=(f0*f4-f1*f3)*id;
						for (int n = 0; n < nen; n++) {
							double dX=dNdX[n][i], dY=dNdY[n][i], dZ=dNdZ[n][i];
							dNdX[n][i] = c00*dX + c01*dY + c02*dZ;
							dNdY[n][i] = c10*dX + c11*dY + c12*dZ;
							dNdZ[n][i] = c20*dX + c21*dY + c22*dZ;
						}
					}

					fBatchMaterial->ComputeStress3D(nel, F3D, sig3D, NULL);

					for (int i = 0; i < nel; i++) {
						double scale = constKd * detJ_cur[i] * w_ref;
						double ss11=sig3D[0][i], ss22=sig3D[1][i], ss33=sig3D[2][i];
						double ss23=sig3D[3][i], ss13=sig3D[4][i], ss12=sig3D[5][i];
						for (int n = 0; n < nen; n++) {
							fx[n][i] += (dNdX[n][i]*ss11 + dNdY[n][i]*ss12 + dNdZ[n][i]*ss13) * scale;
							fy[n][i] += (dNdX[n][i]*ss12 + dNdY[n][i]*ss22 + dNdZ[n][i]*ss23) * scale;
							fz[n][i] += (dNdX[n][i]*ss13 + dNdY[n][i]*ss23 + dNdZ[n][i]*ss33) * scale;
						}
					}
				}
			} /* end IP loop */

			/* HOURGLASS CONTROL (only for 1-IP reduced integration) */
			if (fHourglassType != kNoHourglass && nip == 1) {
				double hg_coeff = fHourglassCoeff;

				if (nsd == 2 && nen == 4) {
					/* Q4 hourglass: one mode, gamma = {1,-1,1,-1}
					 * Flanagan-Belytschko formulation:
					 *   q_hg_x = sum_n gamma_n * u_x_n (hourglass projection)
					 *   F_hg_x_n = coeff * kappa_hg * gamma_n * q_hg_x
					 * where kappa_hg = stiffness scale from element */
					static const double gamma[4] = {1.0, -1.0, 1.0, -1.0};
					const double w_1ip = 4.0; /* Q4 1-point weight */

					for (int i = 0; i < nel; i++) {
						/* element area (from detJ at 1-IP) */
						double area = fabs(detJ_cur[i]) * w_1ip;

						/* characteristic length */
						double h = sqrt(area);

						/* hourglass stiffness scale:
						 * For stiffness: kappa_hg = coeff * (kappa + 4mu/3) / area
						 * For viscous:   kappa_hg = coeff * rho * c * h / area
						 * We use material bulk+shear modulus from the batch material */
						double rho = fBatchMaterial->Density();
						double c = fBatchMaterial->WaveSpeed();
						double kappa_hg;
						if (fHourglassType == kStiffnessHG)
							kappa_hg = hg_coeff * rho * c * c / h;
						else /* kViscousHG */
							kappa_hg = hg_coeff * rho * c;

						/* compute hourglass mode projection */
						double qx = 0.0, qy = 0.0;
						for (int n = 0; n < 4; n++) {
							/* displacement = current - reference */
							double ux = xc[n][i] - xr[n][i];
							double uy = yc[n][i] - yr[n][i];
							qx += gamma[n] * ux;
							qy += gamma[n] * uy;
						}

						/* hourglass force: F_n = kappa_hg * gamma_n * q */
						for (int n = 0; n < 4; n++) {
							fx[n][i] += constKd * kappa_hg * gamma[n] * qx;
							fy[n][i] += constKd * kappa_hg * gamma[n] * qy;
						}
					}
				}
				else if (nsd == 3 && nen == 8) {
					/* Hex8 hourglass: 4 modes
					 * Flanagan-Belytschko base vectors for 8-node hex */
					static const double gamma[4][8] = {
						{ 1,-1, 1,-1,-1, 1,-1, 1},
						{ 1,-1,-1, 1, 1,-1,-1, 1},
						{ 1, 1,-1,-1,-1,-1, 1, 1},
						{-1, 1,-1, 1, 1,-1, 1,-1}
					};
					const double w_1ip = 8.0; /* Hex8 1-point weight */

					for (int i = 0; i < nel; i++) {
						double vol = fabs(detJ_cur[i]) * w_1ip;
						double h = cbrt(vol);
						double rho = fBatchMaterial->Density();
						double c = fBatchMaterial->WaveSpeed();

						double kappa_hg;
						if (fHourglassType == kStiffnessHG)
							kappa_hg = hg_coeff * rho * c * c / h;
						else
							kappa_hg = hg_coeff * rho * c;

						for (int m = 0; m < 4; m++) {
							double qx = 0.0, qy = 0.0, qz = 0.0;
							for (int n = 0; n < 8; n++) {
								qx += gamma[m][n] * (xc[n][i] - xr[n][i]);
								qy += gamma[m][n] * (yc[n][i] - yr[n][i]);
								qz += gamma[m][n] * (zc[n][i] - zr[n][i]);
							}
							for (int n = 0; n < 8; n++) {
								fx[n][i] += constKd * kappa_hg * gamma[m][n] * qx;
								fy[n][i] += constKd * kappa_hg * gamma[m][n] * qy;
								fz[n][i] += constKd * kappa_hg * gamma[m][n] * qz;
							}
						}
					}
				}
			} /* end hourglass */

			/* 6. SCATTER — atomic adds for thread safety */
			for (int i = 0; i < nel; i++) {
				for (int n = 0; n < nen; n++) {
					int node = lconn[n][i];
					#pragma omp atomic
					force(node, 0) += fx[n][i];
					#pragma omp atomic
					force(node, 1) += fy[n][i];
					if (nsd == 3) {
						#pragma omp atomic
						force(node, 2) += fz[n][i];
					}
				}
			}
	} /* end batch (parallel) loop */

	auto t1 = std::chrono::high_resolution_clock::now();
	t_total += std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count();

	/* assemble global forces */
	ElementSupport().AssembleRHS(Group(), force, Field().Equations());

	auto t2 = std::chrono::high_resolution_clock::now();
	t_assemble += std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();

	/* print timing summary */
	if (call_count == 10000) {
		std::cout << "ExplicitElementT timing (" << call_count << " calls):"
		          << " element_force=" << (t_total/1000.0) << "ms"
		          << " assembly=" << (t_assemble/1000.0) << "ms"
		          << " (" << (100.0*t_total/(t_total+t_assemble)) << "% / "
		          << (100.0*t_assemble/(t_total+t_assemble)) << "%)"
		          << std::endl;
		std::cout.flush();
	}
}
