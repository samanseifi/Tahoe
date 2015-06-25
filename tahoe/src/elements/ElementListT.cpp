/* $Id: ElementListT.cpp,v 1.158 2014/08/25 23:59:45 hspark Exp $ */

/* Revision 1.148  2011/11/28 15:26:08  hspark
/* correct 2D DE visco
/*
/* Revision 1.147  2011/11/28 14:29:06  hspark
/* Add 2D DE visco
/*
/* Revision 1.146  2011/11/02 01:05:43  hspark
/* add viscoelastic DE
/*
/* Revision 1.145.2.1  2011/10/29 06:09:07  bcyansfn
/* use c++ standard headers
/*
/* Revision 1.145  2010/11/08 15:33:56  hspark
/* Fixes for FSDielectricElastomerT
/*
/* Revision 1.144  2009/05/21 23:24:45  tdnguye
/* added optimization class
/*
/* Revision 1.143  2009/05/11 21:43:28  regueiro
/*
/* adding xfem development option
/*
/* Revision 1.142  2009/05/11 14:24:50  regueiro
/*
/* adding micromorphic element option
/*
/* Revision 1.141  2009/04/02 00:35:08  lxmota
/* Changed symbol name for mixed Hu-Washizu element.
/*
/* Revision 1.140  2008/12/12 01:00:38  lxmota
/* Add support for neohookean elasticity with simple damage model.
/*
/* Revision 1.139  2008/12/12 00:30:41  lxmota
/* Added conditional support for piezoelectric and HuWashizu elements.
/* */
/* created: paklein (04/20/1998) */
#include "ElementListT.h"
#include "ElementsConfig.h"

#ifdef __DEVELOPMENT__
#include "DevelopmentElementsConfig.h"
#endif

#include <iostream>
#include "ifstreamT.h"
#include "StringT.h"
#include "ElementT.h"
#include "ElementSupportT.h"
#include "GeometryT.h"

#include "ElementBaseT.h"

#ifdef ADHESION_ELEMENT
#include "AdhesionT.h"
#endif

#ifdef CONSTANT_VOLUME_ELEMENT
#include "ConstantVolumeT.h"
#endif

#ifdef COHESIVE_SURFACE_ELEMENT
#include "CSEIsoT.h"
#include "CSEAnisoT.h"
#include "CSESymAnisoT.h"
#include "MeshFreeCSEAnisoT.h"
#include "ThermalSurfaceT.h"

#ifdef COHESIVE_SURFACE_ELEMENT_DEV
#include "RigidCSEAnisoT.h"
#include "NodalRigidCSEAnisoT.h"
#include "CSEAnisoNodal.h"
#include "NodalRigidCSEAnisoMRT.h"
#endif
#endif

#ifdef CONTINUUM_ELEMENT
#include "ViscousDragT.h"
#include "SmallStrainT.h"
#include "SmallStrainAxiT.h"
#include "UpdatedLagrangianT.h"
#include "UpdatedLagrangianAxiT.h"
#include "TotalLagrangianT.h"
#include "TotalLagrangianAxiT.h"
#include "LocalizerT.h"
#include "SimoFiniteStrainT.h"
#include "SimoQ1P0.h"
#include "SimoQ1P0_inv.h"
#include "SimoQ1P0Axi.h"
#include "SimoQ1P0Axi_inv.h"
#include "DiffusionElementT.h"
#include "NLDiffusionElementT.h"
#include "HyperbolicDiffusionElementT.h"
#include "MeshFreeSSSolidT.h"
#include "MeshFreeFSSolidT.h"
#include "MeshFreeFSSolidAxiT.h"
#include "D2MeshFreeFSSolidT.h"
#include "SS_SCNIMFT.h"
#include "FS_SCNIMFT.h"
#include "SS_SCNIMF_AxiT.h"
#include "FS_SCNIMF_AxiT.h"
#include "UpLagr_ExternalFieldT.h"

#ifdef PIEZOELECTRIC
#include "FSPiezoElectricSolidT.h"
#endif

#ifdef DIELECTRIC_ELASTOMER
#include "FSDielectricElastomerT.h"
#endif

#ifdef DIELECTRIC_ELASTOMER_Q1P0
#include "FSDielectricElastomerQ1P0T.h"
#endif

#ifdef DIELECTRIC_ELASTOMER_Q1P0_ELASTOCAPILLARY
#include "FSDielectricElastomerQ1P0SurfaceT.h"
#endif

#ifdef DIELECTRIC_ELASTOMER_Q1P0_VISCO
#include "FSDielectricElastomerQ1P0ViscoT.h"
#endif

#ifdef DIELECTRIC_ELASTOMER_2D
#include "FSDielectricElastomer2DT.h"
#endif

#ifdef DIELECTRIC_ELASTOMER_Q1P02D
#include "FSDielectricElastomerQ1P02DT.h"
#endif

#ifdef DIELECTRIC_ELASTOMER_VISCO
#include "FSDielectricElastomerViscoT.h"
#endif

#ifdef DIELECTRIC_ELASTOMER_2D_VISCO
#include "FSDielectricElastomer2DViscoT.h"
#endif

#ifdef HUWASHIZU
#include "FSHuWashizuUSCT.h"
#endif

#ifdef SIMPLE_SOLID_DEV
#include "TotalLagrangianFlatT.h"
#endif
#ifdef COHESIVE_SURFACE_ELEMENT_DEV
#include "UpLagAdaptiveT.h"
#endif
#endif

#ifdef BRIDGING_ELEMENT
#include "BridgingScaleT.h"
#include "MeshfreeBridgingT.h"
#endif

#ifdef CONTACT_ELEMENT
#include "PenaltyContact2DT.h"
#include "PenaltyContact3DT.h"
#include "AugLagContact2DT.h"
#include "AugLagContact3DT.h"
#include "ACME_Contact3DT.h"
#include "PenaltyContactDrag2DT.h"
#include "PenaltyContactDrag3DT.h"
#ifdef CONTINUUM_ELEMENT /* need meshfree code */
#include "MFPenaltyContact2DT.h"
#endif
#endif

#ifdef PARTICLE_ELEMENT
#include "ParticlePairT.h"
#include "EAMT.h"
#include "ParticleThreeBodyT.h"
#include "ParticleTersoffT.h"
#endif

#ifdef SPRING_ELEMENT
#include "SWDiamondT.h"
#include "MixedSWDiamondT.h"
#include "VirtualSWDC.h"
#include "RodT.h"
#include "UnConnectedRodT.h"
#include "VirtualRodT.h"
#endif

#ifdef CONTACT_ELEMENT_DEV
#include "MultiplierContactElement3DT.h"
#include "MultiplierContactElement2DT.h"
#include "PenaltyContactElement2DT.h"
#include "PenaltyContactElement3DT.h"
#include "FrictionalContactElement2DT.h"
#endif

#ifdef BEM_ELEMENT_DEV
#include "BEMelement.h"
#endif

#ifdef MULTISCALE_ELEMENT_DEV
#include "StaggeredMultiScaleT.h"
#endif

#ifdef MULTISCALE_APS_DEV
#include "APS_AssemblyT.h"
#endif

#ifdef MULTISCALE_APS_V_DEV
#include "APS_V_AssemblyT.h"
#endif

#ifdef MESHFREE_GRAD_PLAST_DEV
#include "MFGPElementT.h"
#endif

#ifdef ENHANCED_STRAIN_LOC_DEV
#include "SmallStrainEnhLocT.h"
#endif

#ifdef ENHANCED_STRAIN_LOC_DEV_CRAIG
#include "SSEnhLocCraigT.h"
#include "SSEnhLocDieterichT.h"
#include "SSEnhLocOpenT.h"
#include "SSEnhLocLinearT.h"
#endif

#ifdef GRAD_SMALL_STRAIN_DEV
#include "GradSmallStrainT.h"
#endif

#ifdef SOLID_ELEMENT_DEV
#ifdef MATERIAL_FORCE_ELEMENT_DEV
#include "SSMF.h"
#endif /* MATERIAL_FORCE_ELEMENT_DEV */

#ifdef SPLIT_INTEGRATION_DEV
#include "SplitIntegrationT.h"
#endif
#endif

#ifdef MIXTURE_THEORY_DEV
#include "MixtureSpeciesT.h"
#include "CurrMixtureSpeciesT.h"
#include "UpdatedLagMixtureT.h"
#include "Q1P0MixtureT.h"
#endif

#ifdef SURFACE_CB_DEV
#include "TotalLagrangianCBSurfaceT.h"
#endif

#ifdef FLUID_ELEMENT_DEV
#include "FluidElementT.h"
#endif

#ifdef FIBER_COMP_DEV
#include "UpLagFiberCompT.h"
#endif

#ifdef SOLID_FLUID_MIX_DEV
#include "FSSolidFluidMixT.h"
#include "FSSolidFluidMixQ8P8T.h"
#endif

#ifdef SOLID_OPTIMIZATION_DEV
//#include "SS_Optimize_Primal.h";
#include "SS_Optimize_Dual.h";
#include "FSFiber_Optimize_Dual.h";
#include "FSFiber_OptSurf.h";
#include "FSFiber_OptNS.h";
#endif

#ifdef MICROMORPHIC_DEV
#include "FSMicromorphic2DT.h"
#include "FSMicromorphic3DT.h"
#endif

#ifdef MICROMORPHIC2_DEV
#include "FSMicromorphic2_3DT.h"
#endif

#ifdef MICROMORPHIC_CURR_CONFIG_DEV
#include "FSMicromorphic3DCurrConfigT.h"
#endif

#ifdef XFEM_DEV
#include "SSXfem2DT.h"
#endif

using namespace Tahoe;

/* constructors */
ElementListT::ElementListT(FEManagerT& fe):
	ParameterInterfaceT("element_list"),
	fHasContact(false)
{
	/* initialize element support */
	fSupport.SetFEManager(&fe);
}

/* destructor */
ElementListT::~ElementListT(void)
{
	/* activate all */
	if (fAllElementGroups.Length() > 0) {
		ArrayT<bool> mask(fAllElementGroups.Length());
		mask = true;
		SetActiveElementGroupMask(mask);
	}
}

/* returns true of ALL element groups have interpolant DOF's */
bool ElementListT::InterpolantDOFs(void) const
{
	bool are_interp = true;
	for (int i = 0; i < Length(); i++)
		if (!fArray[i]->InterpolantDOFs()) are_interp = false;

	return are_interp;
}

/* change the number of active element groups */
void ElementListT::SetActiveElementGroupMask(const ArrayT<bool>& mask)
{
	/* first time */
	if (fAllElementGroups.Length() == 0)
	{
		/* cache all pointers */
		fAllElementGroups.Dimension(Length());
		for (int i = 0; i < fAllElementGroups.Length(); i++)
		{
			ElementBaseT* element = (*this)[i];
			fAllElementGroups[i] = element;
		}
	}

	/* check */
	if (mask.Length() != fAllElementGroups.Length())
		ExceptionT::SizeMismatch("ElementListT::SetActiveElementGroupMask",
			"expecting mask length %d not %d", fAllElementGroups.Length(), mask.Length());

	/* reset active element groups */
	int num_active = 0;
	for (int i = 0; i < mask.Length(); i++)
		if (mask[i])
			num_active++;

	/* cast this to an ArrayT */
	ArrayT<ElementBaseT*>& element_list = *this;
	element_list.Dimension(num_active);
	num_active = 0;
	for (int i = 0; i < mask.Length(); i++)
		if (mask[i])
			element_list[num_active++] = fAllElementGroups[i];
}

/* information about subordinate parameter lists */
void ElementListT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ParameterInterfaceT::DefineSubs(sub_list);

	/* the element groups - an array of choices */
	sub_list.AddSub("element_groups", ParameterListT::OnePlus, true);
}

/* return the description of the given inline subordinate parameter list */
void ElementListT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order,
	SubListT& sub_lists) const
{
	if (name == "element_groups")
	{
		order = ParameterListT::Choice;

#ifdef COHESIVE_SURFACE_ELEMENT
		sub_lists.AddSub("isotropic_CSE");
		sub_lists.AddSub("anisotropic_CSE");
		sub_lists.AddSub("anisotropic_symmetry_CSE");
		sub_lists.AddSub("thermal_CSE");

#ifdef COHESIVE_SURFACE_ELEMENT_DEV
		sub_lists.AddSub("rigid_anisotropic_CSE");
		sub_lists.AddSub("nodal_rigid_anisotropic_CSE");
		sub_lists.AddSub("CSE_rigid_aniso_nodal_int");
		sub_lists.AddSub("nodal_rigid_anisotropic_CSE_MR");
#endif /* COHESIVE_SURFACE_ELEMENT_DEV */

#endif

#ifdef ADHESION_ELEMENT
		sub_lists.AddSub("adhesion");
#endif

#ifdef CONSTANT_VOLUME_ELEMENT
		sub_lists.AddSub("constant_volume");
#endif

#ifdef CONTACT_ELEMENT
		sub_lists.AddSub("contact_2D_penalty");
		sub_lists.AddSub("contact_3D_penalty");

		sub_lists.AddSub("contact_2D_multiplier");
		sub_lists.AddSub("contact_3D_multiplier");

		sub_lists.AddSub("contact_drag_2D_penalty");
		sub_lists.AddSub("contact_drag_3D_penalty");
#endif

#ifdef CONTACT_ELEMENT
#ifdef CONTINUUM_ELEMENT /* need meshfree code */
		sub_lists.AddSub("meshfree_contact_2D_penalty");
#endif
#endif

#ifdef PARTICLE_ELEMENT
		sub_lists.AddSub("particle_pair");
		sub_lists.AddSub("particle_EAM");
		sub_lists.AddSub("particle_three_body");
		sub_lists.AddSub("particle_tersoff");
#endif

#ifdef CONTINUUM_ELEMENT
		sub_lists.AddSub("diffusion");
		sub_lists.AddSub("viscous_drag");
		sub_lists.AddSub("nonlinear_diffusion");
		sub_lists.AddSub("hyperbolic_diffusion");
		sub_lists.AddSub("small_strain");
		sub_lists.AddSub("updated_lagrangian");
		sub_lists.AddSub("updated_lagrangian_Q1P0");
		sub_lists.AddSub("updated_lagrangian_Q1P0_inv");
		sub_lists.AddSub("total_lagrangian");
		sub_lists.AddSub("small_strain_meshfree");
		sub_lists.AddSub("large_strain_meshfree");
		sub_lists.AddSub("small_strain_axi");
		sub_lists.AddSub("updated_lagrangian_axi");
		sub_lists.AddSub("total_lagrangian_axi");
		sub_lists.AddSub("updated_lagrangian_Q1P0_axi");
		sub_lists.AddSub("updated_lagrangian_Q1P0_inv_axi");
		sub_lists.AddSub("large_strain_meshfree_axi");
		sub_lists.AddSub("ss_mfparticle");
		sub_lists.AddSub("fd_mfparticle");
		sub_lists.AddSub("ss_mfparticle_axi");
		sub_lists.AddSub("fd_mfparticle_axi");

#ifdef PIEZOELECTRIC
    sub_lists.AddSub("piezoelectric");
#endif

#ifdef DIELECTRIC_ELASTOMER
    sub_lists.AddSub("dielectric_elastomer");
#endif

#ifdef DIELECTRIC_ELASTOMER_Q1P0
    sub_lists.AddSub("dielectric_elastomer_Q1P0");
#endif

#ifdef DIELECTRIC_ELASTOMER_Q1P0_ELASTOCAPILLARY
    sub_lists.AddSub("dielectric_elastomer_Q1P0Elastocapillary");
#endif

#ifdef DIELECTRIC_ELASTOMER_Q1P0_VISCO
    sub_lists.AddSub("dielectric_elastomer_Q1P0_visco");
#endif

#ifdef DIELECTRIC_ELASTOMER_2D
    sub_lists.AddSub("dielectric_elastomer_2D");
#endif

#ifdef DIELECTRIC_ELASTOMER_Q1P02D
    sub_lists.AddSub("dielectric_elastomer_Q1P02D");
#endif

#ifdef DIELECTRIC_ELASTOMER_VISCO
    sub_lists.AddSub("dielectric_elastomer_visco");
#endif

#ifdef DIELECTRIC_ELASTOMER_2D_VISCO
    sub_lists.AddSub("dielectric_elastomer_2D_visco");
#endif

#ifdef HUWASHIZU
    sub_lists.AddSub("Hu_Washizu_USC");
#endif

#ifdef BRIDGING_ELEMENT
		sub_lists.AddSub("bridging");
		sub_lists.AddSub("meshfree_bridging");
#endif

#ifdef COHESIVE_SURFACE_ELEMENT_DEV
		sub_lists.AddSub("updated_lagrangian_adaptive_insertion");
#endif

#endif /* CONTINUUM_ELEMENT */

#ifdef SPRING_ELEMENT
		sub_lists.AddSub("spring_element");
#endif

#ifdef GRAD_SMALL_STRAIN_DEV
		sub_lists.AddSub("grad_small_strain");
#endif

/*
#ifdef MULTISCALE_ELEMENT_DEV
		sub_lists.AddSub("variational_multiscale");
#endif
*/

#ifdef MULTISCALE_APS_DEV
		sub_lists.AddSub("antiplane_shear_grad_plast");
#endif

/*
#ifdef MULTISCALE_APS_V_DEV
		sub_lists.AddSub("antiplane_shear_grad_plast_V");
#endif
*/

#ifdef MESHFREE_GRAD_PLAST_DEV
		sub_lists.AddSub("mfgp_element");
#endif

#ifdef ENHANCED_STRAIN_LOC_DEV
		sub_lists.AddSub("small_strain_enh_loc");
#endif

#ifdef ENHANCED_STRAIN_LOC_DEV_CRAIG
		sub_lists.AddSub("small_strain_enh_loc_craig");
		sub_lists.AddSub("small_strain_enh_loc_dieterich");
		sub_lists.AddSub("small_strain_enh_loc_open");
		sub_lists.AddSub("small_strain_enh_loc_linear");
#endif

#if defined(SOLID_ELEMENT_DEV) && defined(MATERIAL_FORCE_ELEMENT_DEV)
		sub_lists.AddSub("small_strain_material_force");
#endif

#ifdef MIXTURE_THEORY_DEV
		sub_lists.AddSub("updated_lagrangian_mixture");
		sub_lists.AddSub("Q1P0_mixture");
		sub_lists.AddSub("mixture_species");
		sub_lists.AddSub("current_mixture_species");
#endif

#ifdef CONTACT_ELEMENT_DEV
		sub_lists.AddSub("Jones_penalty_contact_2D");
		sub_lists.AddSub("Jones_penalty_contact_3D");
		sub_lists.AddSub("Jones_multiplier_contact_2D");
		sub_lists.AddSub("Jones_multiplier_contact_3D");
		sub_lists.AddSub("Jones_frictional_contact_2D");
#endif

#ifdef SURFACE_CB_DEV
		sub_lists.AddSub("total_lagrangian_CBsurface");
#endif

#ifdef FLUID_ELEMENT_DEV
		sub_lists.AddSub("incompressible_newtonian_fluid_element");
#endif

#ifdef FIBER_COMP_DEV
		sub_lists.AddSub("uplag_fiber_comp_planar");
#endif

#ifdef SOLID_FLUID_MIX_DEV
		sub_lists.AddSub("total_lagrangian_solid_fluid_mix");
		sub_lists.AddSub("total_lagrangian_solid_fluid_mix_q8p8");
#endif

#ifdef SOLID_OPTIMIZATION_DEV
	sub_lists.AddSub("small_strain_optimize_dual");
	sub_lists.AddSub("uplag_fiber_optimize_dual");
	sub_lists.AddSub("uplag_fiber_optsurf");
	sub_lists.AddSub("uplag_fiber_opt_ns");
#endif

#ifdef MICROMORPHIC_DEV
		sub_lists.AddSub("micromorphic_FS_2D");
		sub_lists.AddSub("micromorphic_FS_3D");
#endif

#ifdef MICROMORPHIC2_DEV
		sub_lists.AddSub("micromorphic2_FS_3D");
#endif

#ifdef MICROMORPHIC_CURR_CONFIG_DEV
		sub_lists.AddSub("micromorphic_FS_3D_inc_curr_config");
#endif

#ifdef XFEM_DEV
		sub_lists.AddSub("xfem_SS_2D");
#endif

	}
	else /* inherited */
		ParameterInterfaceT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* ElementListT::NewSub(const StringT& name) const
{
	/* try to construct element */
	ElementBaseT* element = NewElement(name);
	if (element)
		return element;
	else /* inherited */
		return ParameterInterfaceT::NewSub(name);
}

/* accept parameter list */
void ElementListT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	/* dimension */
	const ArrayT<ParameterListT>& subs = list.Lists();
	Dimension(subs.Length());
	for (int i = 0; i < Length(); i++) {

		/* construct element */
		ElementBaseT* element = NewElement(subs[i].Name());
		if (!element)
			ExceptionT::GeneralFail("ElementListT::TakeParameterList", "could not construct \"%s\"");

		/* store */
		fArray[i] = element;

		/* initialize */
		element->TakeParameterList(subs[i]);
	}
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* return a pointer to a new element group or NULL if the request cannot be completed */
ElementBaseT* ElementListT::NewElement(const StringT& name) const
{
	if (false) /* dummy */
		return NULL;

#ifdef COHESIVE_SURFACE_ELEMENT
	else if (name == "isotropic_CSE")
		return new CSEIsoT(fSupport);

	else if (name == "anisotropic_CSE")
		return new CSEAnisoT(fSupport);

	else if (name == "anisotropic_symmetry_CSE")
		return new CSESymAnisoT(fSupport);

	else if (name == "thermal_CSE")
		return new ThermalSurfaceT(fSupport);

#ifdef COHESIVE_SURFACE_ELEMENT_DEV
	else if (name == "rigid_anisotropic_CSE")
		return new RigidCSEAnisoT(fSupport);
	else if (name == "nodal_rigid_anisotropic_CSE")
		return new NodalRigidCSEAnisoT(fSupport);
	else if (name == "CSE_rigid_aniso_nodal_int")
		return new CSEAnisoNodal(fSupport);
	else if (name == "nodal_rigid_anisotropic_CSE_MR")
		return new NodalRigidCSEAnisoMRT(fSupport);
#endif /* COHESIVE_SURFACE_ELEMENT_DEV */

#endif

#ifdef ADHESION_ELEMENT
	else if (name == "adhesion")
		return new AdhesionT(fSupport);
#endif

#ifdef CONSTANT_VOLUME_ELEMENT
	else if (name == "constant_volume")
	        return new ConstantVolumeT(fSupport);
#endif

#ifdef CONTACT_ELEMENT
	else if (name == "contact_2D_penalty")
		return new PenaltyContact2DT(fSupport);
	else if (name == "contact_3D_penalty")
		return new PenaltyContact3DT(fSupport);

	else if (name == "contact_2D_multiplier")
		return new AugLagContact2DT(fSupport);
	else if (name == "contact_3D_multiplier")
		return new AugLagContact3DT(fSupport);

	else if (name == "contact_drag_2D_penalty")
		return new PenaltyContactDrag2DT(fSupport);
	else if (name == "contact_drag_3D_penalty")
		return new PenaltyContactDrag3DT(fSupport);
#endif

#ifdef CONTACT_ELEMENT
#ifdef CONTINUUM_ELEMENT /* need meshfree code */
	else if (name == "meshfree_contact_2D_penalty")
		return new MFPenaltyContact2DT(fSupport);
#endif
#endif

#ifdef PARTICLE_ELEMENT
	else if (name == "particle_pair")
		return new ParticlePairT(fSupport);
	else if (name == "particle_EAM")
		return new EAMT(fSupport);
	else if (name == "particle_three_body")
		return new ParticleThreeBodyT(fSupport);
	else if (name == "particle_tersoff")
		return new ParticleTersoffT(fSupport);
#endif

#ifdef CONTINUUM_ELEMENT
	else if (name == "diffusion")
		return new DiffusionElementT(fSupport);
	else if (name == "viscous_drag")
		return new ViscousDragT(fSupport);
	else if (name == "nonlinear_diffusion")
		return new NLDiffusionElementT(fSupport);
	else if (name == "hyperbolic_diffusion")
		return new HyperbolicDiffusionElementT(fSupport);
	else if (name == "small_strain")
		return new SmallStrainT(fSupport);
	else if (name == "updated_lagrangian")
		return new UpdatedLagrangianT(fSupport);
	else if (name == "updated_lagrangian_Q1P0")
		return new SimoQ1P0(fSupport);
	else if (name == "updated_lagrangian_Q1P0_inv")
		return new SimoQ1P0_inv(fSupport);
	else if (name == "total_lagrangian")
		return new TotalLagrangianT(fSupport);
	else if (name == "small_strain_meshfree")
		return new MeshFreeSSSolidT(fSupport);
	else if (name == "large_strain_meshfree")
		return new MeshFreeFSSolidT(fSupport);
	else if (name == "small_strain_axi")
		return new SmallStrainAxiT(fSupport);
	else if (name == "updated_lagrangian_axi")
		return new UpdatedLagrangianAxiT(fSupport);
	else if (name == "total_lagrangian_axi")
		return new TotalLagrangianAxiT(fSupport);
	else if (name == "updated_lagrangian_Q1P0_axi")
		return new SimoQ1P0Axi(fSupport);
	else if (name == "updated_lagrangian_Q1P0_inv_axi")
		return new SimoQ1P0Axi_inv(fSupport);
	else if (name == "large_strain_meshfree_axi")
		return new MeshFreeFSSolidAxiT(fSupport);
	else if (name == "ss_mfparticle")
		return new SS_SCNIMFT(fSupport);
	else if (name == "fd_mfparticle")
		return new FS_SCNIMFT(fSupport);
	else if (name == "ss_mfparticle_axi")
		return new SS_SCNIMF_AxiT(fSupport);
	else if (name == "fd_mfparticle_axi")
	  return new FS_SCNIMF_AxiT(fSupport);

#ifdef PIEZOELECTRIC
  else if (name == "piezoelectric")
    return new FSPiezoElectricSolidT(fSupport);
#endif

#ifdef DIELECTRIC_ELASTOMER
  else if (name == "dielectric_elastomer")
    return new FSDielectricElastomerT(fSupport);
#endif

#ifdef DIELECTRIC_ELASTOMER_Q1P0
  else if (name == "dielectric_elastomer_Q1P0")
    return new FSDielectricElastomerQ1P0T(fSupport);
#endif

#ifdef DIELECTRIC_ELASTOMER_Q1P0_ELASTOCAPILLARY
  else if (name == "dielectric_elastomer_Q1P0Elastocapillary")
    return new FSDielectricElastomerQ1P0SurfaceT(fSupport);
#endif

#ifdef DIELECTRIC_ELASTOMER_Q1P0_VISCO
  else if (name == "dielectric_elastomer_Q1P0_visco")
    return new FSDielectricElastomerQ1P0ViscoT(fSupport);
#endif

#ifdef DIELECTRIC_ELASTOMER_2D
  else if (name == "dielectric_elastomer_2D")
    return new FSDielectricElastomer2DT(fSupport);
#endif

#ifdef DIELECTRIC_ELASTOMER_Q1P02D
  else if (name == "dielectric_elastomer_Q1P02D")
    return new FSDielectricElastomerQ1P02DT(fSupport);
#endif

#ifdef DIELECTRIC_ELASTOMER_VISCO
  else if (name == "dielectric_elastomer_visco")
    return new FSDielectricElastomerViscoT(fSupport);
#endif

#ifdef DIELECTRIC_ELASTOMER_2D_VISCO
  else if (name == "dielectric_elastomer_2D_visco")
    return new FSDielectricElastomer2DViscoT(fSupport);
#endif

#ifdef HUWASHIZU
  else if (name == "Hu_Washizu_USC")
    return new FSHuWashizuUSCT(fSupport);
#endif

#ifdef BRIDGING_ELEMENT
	else if (name == "bridging")
		return new BridgingScaleT(fSupport);
	else if (name == "meshfree_bridging")
		return new MeshfreeBridgingT(fSupport);
#endif

#ifdef COHESIVE_SURFACE_ELEMENT_DEV
	else if (name == "updated_lagrangian_adaptive_insertion")
		return new UpLagAdaptiveT(fSupport);
#endif

#endif /* CONTINUUM_ELEMENT */

#ifdef SPRING_ELEMENT
	else if (name == "spring_element")
		return new RodT(fSupport);
#endif

#ifdef GRAD_SMALL_STRAIN_DEV
	else if (name == "grad_small_strain")
		return new GradSmallStrainT(fSupport);
#endif

/*
#ifdef MULTISCALE_ELEMENT_DEV
	else if (name == "variational_multiscale")
		return new StaggeredMultiScaleT(fSupport);
#endif
*/

#ifdef MULTISCALE_APS_DEV
	else if (name == "antiplane_shear_grad_plast")
		return new APS_AssemblyT(fSupport);
#endif

/*
#ifdef MULTISCALE_APS_V_DEV
	else if (name == "antiplane_shear_grad_plast_V")
		return new APS_V_AssemblyT(fSupport);
#endif
*/

#ifdef MESHFREE_GRAD_PLAST_DEV
	else if (name == "mfgp_element")
		return new MFGPElementT(fSupport);
#endif

#ifdef ENHANCED_STRAIN_LOC_DEV
	else if (name == "small_strain_enh_loc")
		return new SmallStrainEnhLocT(fSupport);
#endif

#ifdef ENHANCED_STRAIN_LOC_DEV_CRAIG
	else if (name == "small_strain_enh_loc_craig")
		return new SSEnhLocCraigT(fSupport);
	else if (name == "small_strain_enh_loc_dieterich")
		return new SSEnhLocDieterichT(fSupport);
	else if (name == "small_strain_enh_loc_open")
		return new SSEnhLocOpenT(fSupport);
	else if (name == "small_strain_enh_loc_linear")
		return new SSEnhLocLinearT(fSupport);
#endif

#ifdef MIXTURE_THEORY_DEV
	else if (name == "updated_lagrangian_mixture")
		return new UpdatedLagMixtureT(fSupport);
	else if (name == "Q1P0_mixture")
		return new Q1P0MixtureT(fSupport);
	else if (name == "mixture_species")
		return new MixtureSpeciesT(fSupport);
	else if (name == "current_mixture_species")
		return new CurrMixtureSpeciesT(fSupport);
#endif

#if defined(SOLID_ELEMENT_DEV) && defined(MATERIAL_FORCE_ELEMENT_DEV)
	else if (name == "small_strain_material_force")
		return new SSMF(fSupport);
#endif

#ifdef CONTACT_ELEMENT_DEV
	else if (name == "Jones_penalty_contact_2D")
		return new PenaltyContactElement2DT(fSupport);
	else if (name == "Jones_penalty_contact_3D")
		return new PenaltyContactElement3DT(fSupport);
	else if (name == "Jones_multiplier_contact_2D")
		return new MultiplierContactElement2DT(fSupport);
	else if (name == "Jones_multiplier_contact_3D")
		return new MultiplierContactElement3DT(fSupport);
	else if (name == "Jones_frictional_contact_2D")
		return new FrictionalContactElement2DT(fSupport);
#endif

#ifdef SURFACE_CB_DEV
	else if (name == "total_lagrangian_CBsurface")
		return new TotalLagrangianCBSurfaceT(fSupport);
#endif

#ifdef FLUID_ELEMENT_DEV
	else if (name == "incompressible_newtonian_fluid_element")
		return new FluidElementT(fSupport);
#endif

#ifdef FIBER_COMP_DEV
	else if (name == "uplag_fiber_comp_planar")
	  return new UpLagFiberCompT(fSupport);
#endif

#ifdef SOLID_FLUID_MIX_DEV
	else if (name == "total_lagrangian_solid_fluid_mix")
	  return new FSSolidFluidMixT(fSupport);
	else if (name == "total_lagrangian_solid_fluid_mix_q8p8")
	  return new FSSolidFluidMixQ8P8T(fSupport);
#endif

#ifdef SOLID_OPTIMIZATION_DEV
//	else if (name == "small_strain_optimize_primal")
//	  return new SS_Optimize_Primal(fSupport);
	else if (name == "small_strain_optimize_dual")
	  return new SS_Optimize_Dual(fSupport);
	  
	else if (name == "uplag_fiber_optimize_dual")
	  return new FSFiber_Optimize_Dual(fSupport);

	else if (name == "uplag_fiber_optsurf")
	  return new FSFiber_OptSurf(fSupport);

	else if (name == "uplag_fiber_opt_ns")
	  return new FSFiber_OptNS(fSupport);
#endif

#ifdef MICROMORPHIC_DEV
	else if (name == "micromorphic_FS_2D")
	  return new FSMicromorphic2DT(fSupport);
	else if (name == "micromorphic_FS_3D")
	  return new FSMicromorphic3DT(fSupport);
#endif

#ifdef MICROMORPHIC2_DEV
	else if (name == "micromorphic2_FS_3D")
	  return new FSMicromorphic2_3DT(fSupport);
#endif

#ifdef MICROMORPHIC_CURR_CONFIG_DEV
	else if (name == "micromorphic_FS_3D_inc_curr_config")
	  return new FSMicromorphic3DCurrConfigT(fSupport);
#endif

#ifdef XFEM_DEV
	else if (name == "xfem_SS_2D")
	  return new SSXfem2DT(fSupport);
#endif

	/* default */

	else
		return NULL;
}
