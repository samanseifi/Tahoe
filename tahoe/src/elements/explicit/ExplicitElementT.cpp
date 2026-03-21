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
 *----------------------------------------------------------------------*/
void ExplicitElementT::BatchedInternalForce(double constKd)
{
	int nsd = NumSD();
	int nen = fKernel->NodesPerElement();
	int nip = fKernel->NumIP();
	int ndof = NumDOF();

	/* global coordinate arrays */
	const dArray2DT& init_coords = ElementSupport().InitialCoordinates();
	const dArray2DT& curr_coords = ElementSupport().CurrentCoordinates();
	int numnod = init_coords.MajorDim();

	/* allocate global force accumulation array */
	dArray2DT force(numnod, ndof);
	force = 0.0;

	/* get total number of elements across all blocks */
	int total_elements = 0;
	for (int b = 0; b < fBlockData.Length(); b++)
		total_elements += fBlockData[b].Dimension();

	/* SoA workspace (stack allocated for cache locality) */
	double xc[ExplicitKernelT::MAX_NEN][MVSIZ];
	double yc[ExplicitKernelT::MAX_NEN][MVSIZ];
	double zc[ExplicitKernelT::MAX_NEN][MVSIZ];
	double xr[ExplicitKernelT::MAX_NEN][MVSIZ]; /* reference coords */
	double yr[ExplicitKernelT::MAX_NEN][MVSIZ];
	double zr[ExplicitKernelT::MAX_NEN][MVSIZ];

	/* shape derivatives and Jacobian */
	double dNdx[ExplicitKernelT::MAX_NEN][MVSIZ];
	double dNdy[ExplicitKernelT::MAX_NEN][MVSIZ];
	double dNdz[ExplicitKernelT::MAX_NEN][MVSIZ];
	double detJ[MVSIZ];

	/* reference config shape derivatives (for F computation) */
	double dNdX[ExplicitKernelT::MAX_NEN][MVSIZ];
	double dNdY[ExplicitKernelT::MAX_NEN][MVSIZ];
	double dNdZ[ExplicitKernelT::MAX_NEN][MVSIZ];
	double detJ_ref[MVSIZ];

	/* deformation gradient */
	double F2D[4][MVSIZ];  /* 2D: F11,F12,F21,F22 */
	double F3D[9][MVSIZ];  /* 3D: F11..F33 row-major */

	/* stress */
	double sig11[MVSIZ], sig22[MVSIZ], sig12[MVSIZ];
	double sig3D[6][MVSIZ];

	/* element force accumulators */
	double fx[ExplicitKernelT::MAX_NEN][MVSIZ];
	double fy[ExplicitKernelT::MAX_NEN][MVSIZ];
	double fz[ExplicitKernelT::MAX_NEN][MVSIZ];

	/* local connectivity */
	int lconn[ExplicitKernelT::MAX_NEN][MVSIZ];

	/* process elements in batches of MVSIZ */
	int elem_offset = 0;
	for (int b = 0; b < fBlockData.Length(); b++) {
		int block_nel = fBlockData[b].Dimension();

		for (int batch_start = 0; batch_start < block_nel; batch_start += MVSIZ) {
			int nel = batch_start + MVSIZ <= block_nel
			        ? MVSIZ : block_nel - batch_start;

			/* 1. GATHER: global AoS → local SoA */
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
					if (nsd == 3) {
						zc[n][i] = curr_coords(node, 2);
						zr[n][i] = init_coords(node, 2);
					}
				}
				for (int n = 0; n < nen; n++) {
					fx[n][i] = 0.0;
					fy[n][i] = 0.0;
					if (nsd == 3) fz[n][i] = 0.0;
				}
			}

			/* 2-5. LOOP OVER INTEGRATION POINTS */
			for (int ip = 0; ip < nip; ip++) {
				double w_cur, w_ref;

				/* current config: spatial derivatives (for B^T * sigma) */
				fKernel->ComputeIPData(ip, nel,
					xc, yc, zc, dNdx, dNdy, dNdz, detJ, w_cur);

				/* reference config: material derivatives (for F = dx/dX) */
				fKernel->ComputeIPData(ip, nel,
					xr, yr, zr, dNdX, dNdY, dNdZ, detJ_ref, w_ref);

				if (nsd == 2) {
					/* Compute F = dx/dX = sum_n x_n (dN_n/dX_j) */
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

					/* 4. MATERIAL: batch stress from F */
					fBatchMaterial->ComputeStress2D(nel, F2D,
						sig11, sig22, sig12, NULL);

					/* 5. B^T * sigma — force accumulation
					 * Uses current config spatial derivatives and detJ */
					for (int i = 0; i < nel; i++) {
						double scale = constKd * detJ[i] * w_cur;
						for (int n = 0; n < nen; n++) {
							fx[n][i] += (dNdx[n][i]*sig11[i] + dNdy[n][i]*sig12[i]) * scale;
							fy[n][i] += (dNdx[n][i]*sig12[i] + dNdy[n][i]*sig22[i]) * scale;
						}
					}
				}
				else { /* nsd == 3 */
					/* Compute F = dx/dX */
					for (int i = 0; i < nel; i++) {
						double f[9] = {0.0};
						for (int n = 0; n < nen; n++) {
							f[0] += xc[n][i]*dNdX[n][i]; f[1] += xc[n][i]*dNdY[n][i]; f[2] += xc[n][i]*dNdZ[n][i];
							f[3] += yc[n][i]*dNdX[n][i]; f[4] += yc[n][i]*dNdY[n][i]; f[5] += yc[n][i]*dNdZ[n][i];
							f[6] += zc[n][i]*dNdX[n][i]; f[7] += zc[n][i]*dNdY[n][i]; f[8] += zc[n][i]*dNdZ[n][i];
						}
						for (int k = 0; k < 9; k++) F3D[k][i] = f[k];
					}

					fBatchMaterial->ComputeStress3D(nel, F3D, sig3D, NULL);

					for (int i = 0; i < nel; i++) {
						double scale = constKd * detJ[i] * w_cur;
						double ss11=sig3D[0][i], ss22=sig3D[1][i], ss33=sig3D[2][i];
						double ss23=sig3D[3][i], ss13=sig3D[4][i], ss12=sig3D[5][i];
						for (int n = 0; n < nen; n++) {
							fx[n][i] += (dNdx[n][i]*ss11 + dNdy[n][i]*ss12 + dNdz[n][i]*ss13) * scale;
							fy[n][i] += (dNdx[n][i]*ss12 + dNdy[n][i]*ss22 + dNdz[n][i]*ss23) * scale;
							fz[n][i] += (dNdx[n][i]*ss13 + dNdy[n][i]*ss23 + dNdz[n][i]*ss33) * scale;
						}
					}
				}
			} /* end IP loop */

			/* 6. SCATTER: SoA forces → global force array */
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

	/* assemble global forces into the system RHS */
	ElementSupport().AssembleRHS(Group(), force, Field().Equations());
}
