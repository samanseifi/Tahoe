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

using namespace Tahoe;

ExplicitElementT::ExplicitElementT(const ElementSupportT& support)
	: UpdatedLagrangianT(support),
	  fKernel(NULL),
	  fBatchMaterial(NULL)
{
	SetName("explicit_solid");
}

ExplicitElementT::~ExplicitElementT(void)
{
	delete fKernel;
	delete fBatchMaterial;
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
	/* use the same XML format as updated_lagrangian so the parent
	 * initialization (connectivity, shape functions, equation numbering)
	 * works without modification. The batch material is created from
	 * the block's material parameters in TakeParameterList. */
	UpdatedLagrangianT::DefineSubs(sub_list);
}

ParameterInterfaceT* ExplicitElementT::NewSub(const StringT& name) const
{
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

	if (fKernel && fBatchMaterial)
		std::cout << "ExplicitElementT: MVSIZ=" << MVSIZ
		          << " batched path active (nsd=" << nsd
		          << " nen=" << nen << ")" << std::endl;
	else
		std::cout << "ExplicitElementT: falling back to legacy element loop" << std::endl;
}

/*----------------------------------------------------------------------
 * ElementRHSDriver — override with batched loop
 *----------------------------------------------------------------------*/
void ExplicitElementT::ElementRHSDriver(void)
{
	/* fast path: if kernel and material are available */
	if (fKernel && fBatchMaterial) {
		/* get integrator coefficients */
		double constKd = 0.0;
		int formKd = fIntegrator->FormKd(constKd);
		if (formKd)
			BatchedInternalForce(-constKd);
		return;
	}

	/* fallback: use the generic per-element loop */
	UpdatedLagrangianT::ElementRHSDriver();
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

	/* global force accumulation */
	dArray2DT force(numnod, ndof);
	force = 0.0;

	/* SoA workspace */
	double xc[ExplicitKernelT::MAX_NEN][MVSIZ];
	double yc[ExplicitKernelT::MAX_NEN][MVSIZ];
	double zc[ExplicitKernelT::MAX_NEN][MVSIZ];
	double xr[ExplicitKernelT::MAX_NEN][MVSIZ];
	double yr[ExplicitKernelT::MAX_NEN][MVSIZ];
	double zr[ExplicitKernelT::MAX_NEN][MVSIZ];

	/* shape derivatives */
	double dNdX[ExplicitKernelT::MAX_NEN][MVSIZ]; /* reference */
	double dNdY[ExplicitKernelT::MAX_NEN][MVSIZ];
	double dNdZ[ExplicitKernelT::MAX_NEN][MVSIZ];
	double detJ_cur[MVSIZ];

	/* deformation gradient and stress */
	double F2D[4][MVSIZ];
	double F3D[9][MVSIZ];
	double sig11[MVSIZ], sig22[MVSIZ], sig12[MVSIZ];
	double sig3D[6][MVSIZ];

	/* force accumulators and connectivity */
	double fx[ExplicitKernelT::MAX_NEN][MVSIZ];
	double fy[ExplicitKernelT::MAX_NEN][MVSIZ];
	double fz[ExplicitKernelT::MAX_NEN][MVSIZ];
	int lconn[ExplicitKernelT::MAX_NEN][MVSIZ];

	int elem_offset = 0;
	for (int b = 0; b < fBlockData.Length(); b++) {
		int block_nel = fBlockData[b].Dimension();

		for (int batch_start = 0; batch_start < block_nel; batch_start += MVSIZ) {
			int nel = batch_start + MVSIZ <= block_nel
			        ? MVSIZ : block_nel - batch_start;

			/* 1. GATHER */
			for (int i = 0; i < nel; i++) {
				int global_elem = elem_offset + batch_start + i;
				const ElementCardT& card = ElementCard(global_elem);
				const iArrayT& nodes = card.NodesX();
				for (int n = 0; n < nen; n++) {
					int node = nodes[n];
					lconn[n][i] = node;
					xc[n][i] = curr_coords(node, 0);
					yc[n][i] = curr_coords(node, 1);
					xr[n][i] = init_coords(node, 0);
					yr[n][i] = init_coords(node, 1);
				}
				if (nsd == 3) {
					for (int n = 0; n < nen; n++) {
						int node = nodes[n];
						zc[n][i] = curr_coords(node, 2);
						zr[n][i] = init_coords(node, 2);
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

			/* 6. SCATTER */
			for (int i = 0; i < nel; i++) {
				for (int n = 0; n < nen; n++) {
					int node = lconn[n][i];
					force(node, 0) += fx[n][i];
					force(node, 1) += fy[n][i];
					if (nsd == 3)
						force(node, 2) += fz[n][i];
				}
			}
		} /* end batch loop */

		elem_offset += block_nel;
	} /* end block loop */

	/* assemble global forces */
	ElementSupport().AssembleRHS(Group(), force, Field().Equations());
}
