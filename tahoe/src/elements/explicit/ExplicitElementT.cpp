/* $Id: ExplicitElementT.cpp,v 2.0 2026/05/08 samanseifi Exp $ */
/* ExplicitElementT.cpp — MVSIZ-batched explicit solid element. */
#include "ExplicitElementT.h"

#include "UpdatedLagrangianT.h"
#include "ExplicitKernelT.h"
#include "Q4KernelT.h"
#include "Hex8KernelT.h"
#include "Tet4KernelT.h"
#include "ExplicitMaterialT.h"
#include "ExplNeoHookeanT.h"
#include "ExplJ2PlasticityT.h"
#include "ANPHelperT.h"

#include "ElementSupportT.h"
#include "ModelManagerT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"
#include "eIntegratorT.h"
#include "MaterialListT.h"
#include "ContinuumMaterialT.h"
#include "SolidMaterialT.h"

#include <iostream>
#include <cstring>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif

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
	  fMassScale(NULL),
	  fHistory(NULL),
	  fNumHist(0),
	  fANPEnabled(false),
	  fANP(NULL),
	  fVrefE(NULL),
	  fJe(NULL),
	  fJbarE(NULL)
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
	delete[] fHistory;
	delete fANP;
	delete[] fVrefE;
	delete[] fJe;
	delete[] fJbarE;
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
	sub_list.AddSub("j2_plasticity", ParameterListT::ZeroOrOnce);
	sub_list.AddSub("anp_tet4", ParameterListT::ZeroOrOnce);
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
	if (name == "anp_tet4") {
		ParameterContainerT* anp = new ParameterContainerT(name);
		ParameterT enabled(ParameterT::Boolean, "enabled");
		enabled.SetDefault(true);
		anp->AddParameter(enabled);
		return anp;
	}
	if (name == "j2_plasticity") {
		ParameterContainerT* j2 = new ParameterContainerT(name);
		ParameterT sigY(ParameterT::Double, "sigma_Y");
		sigY.SetDefault(0.0);
		j2->AddParameter(sigY);
		ParameterT hard(ParameterT::Double, "hardening");
		hard.SetDefault(0.0);
		j2->AddParameter(hard);
		return j2;
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
	else if (nsd == 3 && nen == 4)
		fKernel = new Tet4KernelT;
	else {
		std::cout << "ExplicitElementT: unsupported topology nsd="
		          << nsd << " nen=" << nen
		          << " — using legacy element loop" << std::endl;
		fKernel = NULL;
	}

	/* create batch material from the first block's XML parameters.
	 * Scan recursively for mu/kappa/density in the material sub-tree.
	 * If <j2_plasticity> is present, build J2 plasticity; otherwise Neo-Hookean. */
	int num_blocks = list.NumLists("large_strain_element_block");
	if (num_blocks > 0) {
		const ParameterListT& block = list.GetList("large_strain_element_block");
		const ParameterListT& mat_choice = block.GetListChoice(*this, "large_strain_material_choice");

		double mu = 0.0, kappa = 0.0, density = 1.0;
		bool found_mu = false, found_kappa = false;

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
				const ArrayT<ParameterListT>& subs = p.Lists();
				for (int i = 0; i < subs.Length(); i++)
					Find(subs[i], mu, kappa, density, found_mu, found_kappa);
			}
		};

		ParamFinder::Find(mat_choice, mu, kappa, density, found_mu, found_kappa);

		if (found_mu && found_kappa) {
			if (list.NumLists("j2_plasticity") > 0) {
				const ParameterListT& j2 = list.GetList("j2_plasticity");
				double sigma_Y = j2.GetParameter("sigma_Y");
				double H = j2.GetParameter("hardening");
				fBatchMaterial = new ExplJ2PlasticityT(mu, kappa, sigma_Y, H, density);
				std::cout << "ExplicitElementT: J2 plasticity ("
				          << "mu=" << mu << " kappa=" << kappa
				          << " sigma_Y=" << sigma_Y << " H=" << H << ")" << std::endl;
			} else
				fBatchMaterial = new ExplNeoHookeanT(mu, kappa, density);
		}
	}

	/* hourglass control */
	if (list.NumLists("hourglass_control") > 0) {
		const ParameterListT& hg = list.GetList("hourglass_control");
		int hg_type = hg.GetParameter("type");
		fHourglassType = (HourglassTypeT)hg_type;
		fHourglassCoeff = hg.GetParameter("coefficient");
	}

	/* ANP / F-bar (Bonet-Burton): only meaningful for Tet4 (nen=4) */
	if (list.NumLists("anp_tet4") > 0) {
		const ParameterListT& anp = list.GetList("anp_tet4");
		fANPEnabled = anp.GetParameter("enabled");
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

		/* ANP/F-bar: only Tet4 (nen=4 in 3D) supported here.  Skip silently
		 * for hex/quad — not needed since those have full integration. */
		if (fANPEnabled && nsd == 3 && nen == 4) {
			BuildANPRefVolumes();
			fJe    = new double[fTotalElements]();
			fJbarE = new double[fTotalElements]();
			fANP = new ANPHelperT();
			int numnod = ElementSupport().InitialCoordinates().MajorDim();
			fANP->Init(fTotalElements, numnod, nen, fFlatConn, fVrefE);
			std::cout << "ExplicitElementT: ANP-Tet4 (LS-DYNA ELFORM=13) enabled "
			          << "for " << fTotalElements << " tets, "
			          << numnod << " nodes" << std::endl;
		} else if (fANPEnabled) {
			fANPEnabled = false;   /* only valid for Tet4 */
		}

		/* allocate and initialize history variables (plasticity, etc.) */
		fNumHist = fBatchMaterial->NumHistoryVars();
		if (fNumHist > 0) {
			int nip = fKernel->NumIP();
			long nwords = (long)fTotalElements * nip * fNumHist;
			fHistory = new double[nwords]();   /* zero-initialized */
			fBatchMaterial->InitializeHistory(nip, fTotalElements, fHistory);
			std::cout << "ExplicitElementT: history "
			          << fNumHist << " vars/IP x "
			          << nip << " IPs x "
			          << fTotalElements << " elems = "
			          << nwords << " doubles ("
			          << (nwords * 8) / (1024*1024) << " MB)" << std::endl;
		}

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
 * BuildANPRefVolumes — precompute reference tet volumes for F-bar averaging
 *----------------------------------------------------------------------*/
void ExplicitElementT::BuildANPRefVolumes(void)
{
	const int nen = fKernel->NodesPerElement();
	const dArray2DT& ref = ElementSupport().InitialCoordinates();

	delete[] fVrefE;
	fVrefE = new double[fTotalElements];
	for (int e = 0; e < fTotalElements; e++) {
		const int* ec = fFlatConn + e * nen;
		double x[4], y[4], z[4];
		for (int n = 0; n < 4; n++) {
			x[n] = ref(ec[n], 0); y[n] = ref(ec[n], 1); z[n] = ref(ec[n], 2);
		}
		double a1=x[0]-x[2], a2=y[0]-y[2], a3=z[0]-z[2];
		double b1=x[1]-x[2], b2=y[1]-y[2], b3=z[1]-z[2];
		double c1=x[3]-x[2], c2=y[3]-y[2], c3=z[3]-z[2];
		fVrefE[e] = std::fabs(a1*(b2*c3-b3*c2)
		                    - a2*(b1*c3-b3*c1)
		                    + a3*(b1*c2-b2*c1)) / 6.0;
	}
}

/*----------------------------------------------------------------------
 * ComputeAllJe — current-config tet volume / reference volume = J = det(F)
 *----------------------------------------------------------------------*/
void ExplicitElementT::ComputeAllJe(void)
{
	const int nen = fKernel->NodesPerElement();
	const dArray2DT& cur = ElementSupport().CurrentCoordinates();

	#pragma omp parallel for if(fTotalElements > 1024)
	for (int e = 0; e < fTotalElements; e++) {
		const int* ec = fFlatConn + e * nen;
		double x[4], y[4], z[4];
		for (int n = 0; n < 4; n++) {
			x[n] = cur(ec[n], 0); y[n] = cur(ec[n], 1); z[n] = cur(ec[n], 2);
		}
		double a1=x[0]-x[2], a2=y[0]-y[2], a3=z[0]-z[2];
		double b1=x[1]-x[2], b2=y[1]-y[2], b3=z[1]-z[2];
		double c1=x[3]-x[2], c2=y[3]-y[2], c3=z[3]-z[2];
		double V_cur = std::fabs(a1*(b2*c3-b3*c2)
		                       - a2*(b1*c3-b3*c1)
		                       + a3*(b1*c2-b2*c1)) / 6.0;
		fJe[e] = (fVrefE[e] > 0.0) ? V_cur / fVrefE[e] : 1.0;
	}
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
			/* 3D: characteristic length from element volume */
			double x[8], y[8], z[8];
			for (int n = 0; n < nen; n++) {
				x[n] = ref_coords(ec[n], 0);
				y[n] = ref_coords(ec[n], 1);
				z[n] = ref_coords(ec[n], 2);
			}
			double vol;
			if (nen == 8) {
				/* Hex8: 1/6 * |diag1 x diag2 . diag3| */
				double dx1 = x[6]-x[0], dy1 = y[6]-y[0], dz1 = z[6]-z[0];
				double dx2 = x[7]-x[1], dy2 = y[7]-y[1], dz2 = z[7]-z[1];
				double dx3 = x[5]-x[3], dy3 = y[5]-y[3], dz3 = z[5]-z[3];
				vol = fabs(dx1*(dy2*dz3-dz2*dy3)
				         - dy1*(dx2*dz3-dz2*dx3)
				         + dz1*(dx2*dy3-dy2*dx3)) / 6.0;
			} else {
				/* Tet4: 1/6 * |(x_0-x_2) x (x_1-x_2) . (x_3-x_2)|
				 * (uses Tahoe Tet4 node ordering — node 2 is the apex at
				 * parametric origin) */
				double a1 = x[0]-x[2], a2 = y[0]-y[2], a3 = z[0]-z[2];
				double b1 = x[1]-x[2], b2 = y[1]-y[2], b3 = z[1]-z[2];
				double c1 = x[3]-x[2], c2 = y[3]-y[2], c3 = z[3]-z[2];
				vol = fabs(a1*(b2*c3-b3*c2)
				         - a2*(b1*c3-b3*c1)
				         + a3*(b1*c2-b2*c1)) / 6.0;
			}
			double h = cbrt(vol);
			double c = fBatchMaterial->WaveSpeed();
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
			double vol;
			if (nen == 8) {
				double dx1=x[6]-x[0], dy1=y[6]-y[0], dz1=z[6]-z[0];
				double dx2=x[7]-x[1], dy2=y[7]-y[1], dz2=z[7]-z[1];
				double dx3=x[5]-x[3], dy3=y[5]-y[3], dz3=z[5]-z[3];
				vol = fabs(dx1*(dy2*dz3-dz2*dy3)
				         - dy1*(dx2*dz3-dz2*dx3)
				         + dz1*(dx2*dy3-dy2*dx3)) / 6.0;
			} else {
				/* Tet4 */
				double a1 = x[0]-x[2], a2 = y[0]-y[2], a3 = z[0]-z[2];
				double b1 = x[1]-x[2], b2 = y[1]-y[2], b3 = z[1]-z[2];
				double c1 = x[3]-x[2], c2 = y[3]-y[2], c3 = z[3]-z[2];
				vol = fabs(a1*(b2*c3-b3*c2)
				         - a2*(b1*c3-b3*c1)
				         + a3*(b1*c2-b2*c1)) / 6.0;
			}
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
 * LHSDriver — override to apply mass scaling to lumped mass
 *----------------------------------------------------------------------*/
void ExplicitElementT::LHSDriver(GlobalT::SystemTypeT sys_type)
{
	/* if no mass scaling, use the parent's assembly */
	if (!fMassScale || fMassScalingType == kNoMassScaling) {
		SolidElementT::LHSDriver(sys_type);
		return;
	}

	/* inherited: constraint equations from ContinuumElementT */
	ContinuumElementT::LHSDriver(sys_type);

	/* mass assembly with per-element scaling */
	double constM = 0.0;
	double constK = 0.0;
	int formM = fIntegrator->FormM(constM);
	int formK = fIntegrator->FormK(constK);
	if (fMassType == kNoMass) formM = 0;
	if (formM == 0 && formK == 0) return;

	bool axisymmetric = Axisymmetric();
	int elem_index = 0;
	Top();
	while (NextElement())
	{
		if (CurrentElement().Flag() != ElementCardT::kOFF)
		{
			fLHS = 0.0;
			SetGlobalShape();

			/* element mass with scale factor */
			double constMe = constM;
			if (fabs(constMe) > kSmall) {
				double density = fCurrMaterial->Density();
				double scaled_density = density * fMassScale[elem_index];
				FormMass(fMassType, constMe * scaled_density, axisymmetric, NULL);
			}

			/* stiffness (explicit: constK = 0, this block is skipped) */
			if (fabs(constK) > kSmall)
				FormStiffness(constK);

			AssembleLHS();
		}
		elem_index++;
	}
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

	/* ANP / F-bar pre-pass: compute J_e for all elements, average at nodes,
	 * scatter back to J_bar_e.  Done before the batched force loop so the
	 * J_bar value for each element is ready when its batch is processed. */
	if (fANPEnabled && fANP) {
		ComputeAllJe();
		fANP->ComputeJBar(fJe, fJbarE);
	}

	/* per-element RHS vector for assembly */
	dArrayT elem_rhs(nen * ndof);

	/* raw coordinate pointers for direct indexing */
	const double* cur_ptr = curr_coords.Pointer();
	const double* ref_ptr = init_coords.Pointer();
	int coord_stride = init_coords.MinorDim();


	/* total batches across all blocks */
	int num_batches = (fTotalElements + MVSIZ - 1) / MVSIZ;

	/* OpenMP parallel loop over batches.
	 *
	 * Threshold: only parallelise when there are enough batches to amortise
	 * per-task sync overhead.  Heuristic — at least 4 batches per available
	 * thread.  Below that, run serial (avoids the small-problem slowdown
	 * documented in #32, where 250-hex tests ran 2-3x SLOWER with 12 threads
	 * than with 1 thread).
	 *
	 * This is a static decision per call; the value of fThreadingActive is
	 * also stored for the timing-banner stats below. */
	int max_threads = 1;
#ifdef _OPENMP
	max_threads = omp_get_max_threads();
#endif
	bool use_omp = (num_batches >= 4 * max_threads);

	#pragma omp parallel for schedule(dynamic, 1) if(use_omp)
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

		/* per-batch history scratch buffer: [var * MVSIZ + i] layout */
		static const int MAX_HIST = 20;
		double hist_batch[MAX_HIST * MVSIZ];

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

			/* 2-5. IP LOOP: two kernel calls (reference + current) */
			for (int ip = 0; ip < nip; ip++) {
				double w_ref, w_cur;

				/* reference config: dN/dX for F computation */
				fKernel->ComputeIPData(ip, nel,
					xr, yr, zr, dNdX, dNdY, dNdZ, detJ_cur, w_ref);

				/* current config: dN/dx for B^T*sigma */
				double dNdx_cur[ExplicitKernelT::MAX_NEN][MVSIZ];
				double dNdy_cur[ExplicitKernelT::MAX_NEN][MVSIZ];
				double dNdz_cur[ExplicitKernelT::MAX_NEN][MVSIZ];
				double detJ_actual[MVSIZ];
				fKernel->ComputeIPData(ip, nel,
					xc, yc, zc, dNdx_cur, dNdy_cur, dNdz_cur, detJ_actual, w_cur);

				if (nsd == 2) {
					/* F = dx/dX = sum x_cur * dN/dX */
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
					}

					/* gather history slice for this batch at this IP */
					double* hist_ptr = NULL;
					if (fNumHist > 0 && fHistory) {
						hist_ptr = hist_batch;
						const double* hsrc = fHistory + ((long)ip * fNumHist * fTotalElements);
						for (int var = 0; var < fNumHist; var++) {
							const double* src = hsrc + (long)var * fTotalElements + batch_start;
							double* dst = hist_batch + var * MVSIZ;
							for (int i = 0; i < nel; i++) dst[i] = src[i];
						}
					}

					fBatchMaterial->ComputeStress2D(nel, F2D,
						sig11, sig22, sig12, hist_ptr);

					/* scatter history back */
					if (hist_ptr) {
						double* hdst = fHistory + ((long)ip * fNumHist * fTotalElements);
						for (int var = 0; var < fNumHist; var++) {
							double* dst = hdst + (long)var * fTotalElements + batch_start;
							const double* src = hist_batch + var * MVSIZ;
							for (int i = 0; i < nel; i++) dst[i] = src[i];
						}
					}

					/* B^T * sigma using current-config dN/dx */
					for (int i = 0; i < nel; i++) {
						double scale = constKd * detJ_actual[i] * w_cur;
						for (int n = 0; n < nen; n++) {
							fx[n][i] += (dNdx_cur[n][i]*sig11[i] + dNdy_cur[n][i]*sig12[i]) * scale;
							fy[n][i] += (dNdx_cur[n][i]*sig12[i] + dNdy_cur[n][i]*sig22[i]) * scale;
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
					}

					/* gather history slice for this batch at this IP */
					double* hist_ptr3 = NULL;
					if (fNumHist > 0 && fHistory) {
						hist_ptr3 = hist_batch;
						const double* hsrc = fHistory + ((long)ip * fNumHist * fTotalElements);
						for (int var = 0; var < fNumHist; var++) {
							const double* src = hsrc + (long)var * fTotalElements + batch_start;
							double* dst = hist_batch + var * MVSIZ;
							for (int i = 0; i < nel; i++) dst[i] = src[i];
						}
					}

					/* ANP / F-bar correction: scale F by (J_bar/J)^(1/3)
					 * before passing to material so the material sees the
					 * nodal-averaged dilatation.  Removes Tet4 volumetric
					 * locking under near-incompressibility. */
					if (fANPEnabled && fANP) {
						for (int i = 0; i < nel; i++) {
							double J = fJe[batch_start + i];
							double Jb = fJbarE[batch_start + i];
							double scale = (J > 0.0) ? std::cbrt(Jb / J) : 1.0;
							for (int k = 0; k < 9; k++)
								F3D[k][i] *= scale;
						}
					}

					fBatchMaterial->ComputeStress3D(nel, F3D, sig3D, hist_ptr3);

					/* scatter history back */
					if (hist_ptr3) {
						double* hdst = fHistory + ((long)ip * fNumHist * fTotalElements);
						for (int var = 0; var < fNumHist; var++) {
							double* dst = hdst + (long)var * fTotalElements + batch_start;
							const double* src = hist_batch + var * MVSIZ;
							for (int i = 0; i < nel; i++) dst[i] = src[i];
						}
					}

					/* B^T * sigma using current-config dN/dx */
					for (int i = 0; i < nel; i++) {
						double scale = constKd * detJ_actual[i] * w_cur;
						double ss11=sig3D[0][i], ss22=sig3D[1][i], ss33=sig3D[2][i];
						double ss23=sig3D[3][i], ss13=sig3D[4][i], ss12=sig3D[5][i];
						for (int n = 0; n < nen; n++) {
							fx[n][i] += (dNdx_cur[n][i]*ss11 + dNdy_cur[n][i]*ss12 + dNdz_cur[n][i]*ss13) * scale;
							fy[n][i] += (dNdx_cur[n][i]*ss12 + dNdy_cur[n][i]*ss22 + dNdz_cur[n][i]*ss23) * scale;
							fz[n][i] += (dNdx_cur[n][i]*ss13 + dNdy_cur[n][i]*ss23 + dNdz_cur[n][i]*ss33) * scale;
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

			/* 6. SCATTER — per-element assembly via equation numbers */
			for (int i = 0; i < nel; i++) {
				int global_elem = batch_start + i;
				const ElementCardT& card = ElementCard(global_elem);
				const iArrayT& eqnos = card.Equations();
				dArrayT erhs(nen * ndof);
				for (int n = 0; n < nen; n++) {
					erhs[n*ndof + 0] = fx[n][i];
					erhs[n*ndof + 1] = fy[n][i];
					if (nsd == 3)
						erhs[n*ndof + 2] = fz[n][i];
				}
				#pragma omp critical
				ElementSupport().AssembleRHS(Group(), erhs, eqnos);
			}
	} /* end batch (parallel) loop */

	auto t1 = std::chrono::high_resolution_clock::now();
	t_total += std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count();

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
