# $Id$
# This file was generated by MakeConfigMakeFile.pl from DevelopmentElementsConfig.conf
# created: Tue Jul 28 14:52:26 EDT 2015
#
# \file DevelopmentElementsConfig.make
# Configuration of optional components within Tahoe.
# Sections of the code are included or excluded in the build of Tahoe depending in 
# this flags in this file and in the file DevelopmentElementsConfig.h. Each option has
# a #define definition in DevelopmentElementsConfig.h and a corresponding directory definition
# in this file. The two items must be set consistently to enable or
# disable materials models. To enable an option:
# -# in this file, uncomment the macro.
# -# in DevelopmentElementsConfig.h, uncomment the #define statement
#
# The naming convention for the definitions in this file and the macros in
# DevelopmentElementsConfig.h are as follows. For the option [OPTION]:
# -# the macro below defining the corresponding source directory will be DIRECTORY_[OPTION]
# -# the symbol in DevelopmentElementsConfig.h file will be [OPTION]

# \def DIRECTORY_BEM_ELEMENT_DEV
# Boundary element method element.
# This option must be set in conjunction with #define BEM_ELEMENT_DEV
# in DevelopmentElementsConfig.h. */
#DIRECTORY_BEM_ELEMENT_DEV = BEM_element

# \def DIRECTORY_PML_ELEMENT_DEV
# Small strain, perfectly matched layer formulation.
# This option must be set in conjunction with #define PML_ELEMENT_DEV
# in DevelopmentElementsConfig.h. */
#DIRECTORY_PML_ELEMENT_DEV = PMLElement

# \def DIRECTORY_COHESIVE_SURFACE_ELEMENT_DEV
# Cohesive surface elements.
# This option must be set in conjunction with #define COHESIVE_SURFACE_ELEMENT_DEV
# in DevelopmentElementsConfig.h. */
DIRECTORY_COHESIVE_SURFACE_ELEMENT_DEV = cohesive_surface

# \def DIRECTORY_CONTACT_ELEMENT_DEV
# Methods for enforcing contact constraints.
# This option must be set in conjunction with #define CONTACT_ELEMENT_DEV
# in DevelopmentElementsConfig.h. */
#DIRECTORY_CONTACT_ELEMENT_DEV = contact

# \def DIRECTORY_MULTISCALE_ELEMENT_DEV
# Variational multi-scale method.
# This option must be set in conjunction with #define MULTISCALE_ELEMENT_DEV
# in DevelopmentElementsConfig.h. */
#DIRECTORY_MULTISCALE_ELEMENT_DEV = multiscale

# \def DIRECTORY_SOLID_ELEMENT_DEV
# Solid elements and constitutive models.
# This option must be set in conjunction with #define SOLID_ELEMENT_DEV
# in DevelopmentElementsConfig.h. */
DIRECTORY_SOLID_ELEMENT_DEV = solid

# \def DIRECTORY_MATERIAL_FORCE_ELEMENT_DEV
# Material force formulations.
# This option must be set in conjunction with #define MATERIAL_FORCE_ELEMENT_DEV
# in DevelopmentElementsConfig.h. */
#DIRECTORY_MATERIAL_FORCE_ELEMENT_DEV = solid/materialforce

# \def DIRECTORY_GRAD_SMALL_STRAIN_DEV
# Gradient enhanced C1 small strain element.
# This option must be set in conjunction with #define GRAD_SMALL_STRAIN_DEV
# in DevelopmentElementsConfig.h. */
#DIRECTORY_GRAD_SMALL_STRAIN_DEV = dorgan

# \def DIRECTORY_MULTISCALE_APS_DEV
# coupled FE implementation of anti-plane shear curl gradient plasticity.
# This option must be set in conjunction with #define MULTISCALE_APS_DEV
# in DevelopmentElementsConfig.h. */
#DIRECTORY_MULTISCALE_APS_DEV = APS_grad

# \def DIRECTORY_SPLIT_INTEGRATION_DEV
# split integration element.
# This option must be set in conjunction with #define SPLIT_INTEGRATION_DEV
# in DevelopmentElementsConfig.h. */
#DIRECTORY_SPLIT_INTEGRATION_DEV = solid/hspark

# \def DIRECTORY_SIMPLE_SOLID_DEV
# simplied implementation of finite strain solid.
# This option must be set in conjunction with #define SIMPLE_SOLID_DEV
# in DevelopmentElementsConfig.h. */
DIRECTORY_SIMPLE_SOLID_DEV = solid/simple_solid

# \def DIRECTORY_MULTISCALE_APS_V_DEV
# coupled vector FE implementation of anti-plane shear curl gradient plasticity.
# This option must be set in conjunction with #define MULTISCALE_APS_V_DEV
# in DevelopmentElementsConfig.h. */
#DIRECTORY_MULTISCALE_APS_V_DEV = APS_grad_vector

# \def DIRECTORY_MESHFREE_GRAD_PLAST_DEV
# meshfree implementation of gradient plasticity.
# This option must be set in conjunction with #define MESHFREE_GRAD_PLAST_DEV
# in DevelopmentElementsConfig.h. */
#DIRECTORY_MESHFREE_GRAD_PLAST_DEV = meshfree_grad_plast

# \def DIRECTORY_ENHANCED_STRAIN_LOC_DEV
# enhanced strain embedded discontinuity element.
# This option must be set in conjunction with #define ENHANCED_STRAIN_LOC_DEV
# in DevelopmentElementsConfig.h. */
#DIRECTORY_ENHANCED_STRAIN_LOC_DEV = enhanced_strain_loc

# \def DIRECTORY_ENHANCED_STRAIN_LOC_DEV_CRAIG
# enhanced strain embedded discontinuity element, 2nd version.
# This option must be set in conjunction with #define ENHANCED_STRAIN_LOC_DEV_CRAIG
# in DevelopmentElementsConfig.h. */
#DIRECTORY_ENHANCED_STRAIN_LOC_DEV_CRAIG = craig_enhanced_strain_loc

# \def DIRECTORY_MIXTURE_THEORY_DEV
# mixture theory.
# This option must be set in conjunction with #define MIXTURE_THEORY_DEV
# in DevelopmentElementsConfig.h. */
DIRECTORY_MIXTURE_THEORY_DEV = mixture

# \def DIRECTORY_BRIDGING_ELEMENT_DEV
# multiscale bridging methods.
# This option must be set in conjunction with #define BRIDGING_ELEMENT_DEV
# in DevelopmentElementsConfig.h. */
#DIRECTORY_BRIDGING_ELEMENT_DEV = NONE

# \def DIRECTORY_SURFACE_CB_DEV
# Cauchy-Born approach to capture surface effects using pair potentials.
# This option must be set in conjunction with #define SURFACE_CB_DEV
# in DevelopmentElementsConfig.h. */
#DIRECTORY_SURFACE_CB_DEV = surface_CB

# \def DIRECTORY_SURFACE_CB_EAM_DEV
# Cauchy-Born approach to capture surface effects using EAM potentials.
# This option must be set in conjunction with #define SURFACE_CB_EAM_DEV
# in DevelopmentElementsConfig.h. */
#DIRECTORY_SURFACE_CB_EAM_DEV = surface_CB_EAM

# \def DIRECTORY_SURFACE_CB_SI_DEV
# Cauchy-Born approach to capture surface effects using Tersoff potentials.
# This option must be set in conjunction with #define SURFACE_CB_SI_DEV
# in DevelopmentElementsConfig.h. */
#DIRECTORY_SURFACE_CB_SI_DEV = surface_CB_Si

# \def DIRECTORY_SURFACE_CB_SI_DIMER_DEV
# Cauchy-Born approach to capture {100} surface reconstructions in Silicon.
# This option must be set in conjunction with #define SURFACE_CB_SI_DIMER_DEV
# in DevelopmentElementsConfig.h. */
#DIRECTORY_SURFACE_CB_SI_DIMER_DEV = surface_CB_Si_Dimer

# \def DIRECTORY_SURFACE_CB_WURTZITE_DEV
# Cauchy-Born approach to capture surface effects using Tersoff-Abel potentials for Wurtzite lattices.
# This option must be set in conjunction with #define SURFACE_CB_WURTZITE_DEV
# in DevelopmentElementsConfig.h. */
#DIRECTORY_SURFACE_CB_WURTZITE_DEV = surface_CB_Wurtzite

# \def DIRECTORY_SURFACE_CB_ZB_DEV
# Cauchy-Born approach to capture surface effects using Tersoff-Abel potentials for Zinc Blende lattices.
# This option must be set in conjunction with #define SURFACE_CB_ZB_DEV
# in DevelopmentElementsConfig.h. */
#DIRECTORY_SURFACE_CB_ZB_DEV = surface_CB_ZB

# \def DIRECTORY_FLUID_ELEMENT_DEV
# Incompressible Newtonian fluid element (constant density/dynamic viscosity).
# This option must be set in conjunction with #define FLUID_ELEMENT_DEV
# in DevelopmentElementsConfig.h. */
#DIRECTORY_FLUID_ELEMENT_DEV = fluid_element

# \def DIRECTORY_FIBER_COMP_DEV
# Updated Langriagian element with element level decription of fiber orientations.
# This option must be set in conjunction with #define FIBER_COMP_DEV
# in DevelopmentElementsConfig.h. */
#DIRECTORY_FIBER_COMP_DEV = fiber_composite

# \def DIRECTORY_SOLID_FLUID_MIX_DEV
# 3D Total Lagrangian dynamic solid fluid mixture element.
# This option must be set in conjunction with #define SOLID_FLUID_MIX_DEV
# in DevelopmentElementsConfig.h. */
#DIRECTORY_SOLID_FLUID_MIX_DEV = solid_fluid_mix

# \def DIRECTORY_SOLID_OPTIMIZATION_DEV
# Implementation of the adjoint method for inverse elasticity problems.
# This option must be set in conjunction with #define SOLID_OPTIMIZATION_DEV
# in DevelopmentElementsConfig.h. */
#DIRECTORY_SOLID_OPTIMIZATION_DEV = optimization

# \def DIRECTORY_DEM_COUPLING_DEV
# 3D ellipsoidal particle discrete element for coupling.
# This option must be set in conjunction with #define DEM_COUPLING_DEV
# in DevelopmentElementsConfig.h. */
#DIRECTORY_DEM_COUPLING_DEV = DEM_coupling

# \def DIRECTORY_PIEZOELECTRIC
# piezoelectric element.
# This option must be set in conjunction with #define PIEZOELECTRIC
# in DevelopmentElementsConfig.h. */
#DIRECTORY_PIEZOELECTRIC = piezoelectric

# \def DIRECTORY_HUWASHIZU
# Hu-Washizu element.
# This option must be set in conjunction with #define HUWASHIZU
# in DevelopmentElementsConfig.h. */
#DIRECTORY_HUWASHIZU = HuWashizu

# \def DIRECTORY_MICROMORPHIC_DEV
# Micromorphic finite strain elastoplastic coupled element.
# This option must be set in conjunction with #define MICROMORPHIC_DEV
# in DevelopmentElementsConfig.h. */
#DIRECTORY_MICROMORPHIC_DEV = micromorphic

# \def DIRECTORY_MICROMORPHIC2_DEV
# Micromorphic finite strain elastoplastic coupled element, with changes.
# This option must be set in conjunction with #define MICROMORPHIC2_DEV
# in DevelopmentElementsConfig.h. */
#DIRECTORY_MICROMORPHIC2_DEV = micromorphic2

# \def DIRECTORY_MICROMORPHIC_CURR_CONFIG_DEV
# Micromorphic finite strain elastoplastic coupled element, current configuration.
# This option must be set in conjunction with #define MICROMORPHIC_CURR_CONFIG_DEV
# in DevelopmentElementsConfig.h. */
#DIRECTORY_MICROMORPHIC_CURR_CONFIG_DEV = micromorphic_curr_config

# \def DIRECTORY_XFEM_DEV
# extended finite element implementation for crack problems.
# This option must be set in conjunction with #define XFEM_DEV
# in DevelopmentElementsConfig.h. */
#DIRECTORY_XFEM_DEV = xfem

# \def DIRECTORY_DIELECTRIC_ELASTOMER
# Coupled electromechanical field equations for 3D dielectric elastomers.
# This option must be set in conjunction with #define DIELECTRIC_ELASTOMER
# in DevelopmentElementsConfig.h. */
DIRECTORY_DIELECTRIC_ELASTOMER = Dielectric_Elastomer

# \def DIRECTORY_DIELECTRIC_ELASTOMER_Q1P0
# 3D Dielectric Elastomers with (hopefully) locking-free Q1P0 formulation.
# This option must be set in conjunction with #define DIELECTRIC_ELASTOMER_Q1P0
# in DevelopmentElementsConfig.h. */
DIRECTORY_DIELECTRIC_ELASTOMER_Q1P0 = Dielectric_Elastomer_Q1P0

# \def DIRECTORY_DIELECTRIC_ELASTOMER_Q1P0_ELASTOCAPILLARY
# 3D Dielectric Elastomers with locking-free Q1P0 formulation and surface tension.
# This option must be set in conjunction with #define DIELECTRIC_ELASTOMER_Q1P0_ELASTOCAPILLARY
# in DevelopmentElementsConfig.h. */
DIRECTORY_DIELECTRIC_ELASTOMER_Q1P0_ELASTOCAPILLARY = Dielectric_Elastomer_Q1P0_Elastocapillary

# \def DIRECTORY_DIELECTRIC_ELASTOMER_Q1P0_VISCO
# 3D Dielectric Elastomers with (hopefully) locking-free Q1P0 formulation and viscoelasticity.
# This option must be set in conjunction with #define DIELECTRIC_ELASTOMER_Q1P0_VISCO
# in DevelopmentElementsConfig.h. */
DIRECTORY_DIELECTRIC_ELASTOMER_Q1P0_VISCO = Dielectric_Elastomer_Q1P0_Visco

# \def DIRECTORY_DIELECTRIC_ELASTOMER_2D
# Coupled electromechanical field equations for 2D dielectric elastomers.
# This option must be set in conjunction with #define DIELECTRIC_ELASTOMER_2D
# in DevelopmentElementsConfig.h. */
DIRECTORY_DIELECTRIC_ELASTOMER_2D = Dielectric_Elastomer_2D

# \def DIRECTORY_DIELECTRIC_ELASTOMER_Q1P02D
# 2D Dielectric Elastomers with (hopefully) locking-free Q1P0 formulation.
# This option must be set in conjunction with #define DIELECTRIC_ELASTOMER_Q1P02D
# in DevelopmentElementsConfig.h. */
DIRECTORY_DIELECTRIC_ELASTOMER_Q1P02D = Dielectric_Elastomer_Q1P02D

# \def DIRECTORY_DIELECTRIC_ELASTOMER_VISCO
# 3D dielectric elastomer with viscoelastic material behavior.
# This option must be set in conjunction with #define DIELECTRIC_ELASTOMER_VISCO
# in DevelopmentElementsConfig.h. */
#DIRECTORY_DIELECTRIC_ELASTOMER_VISCO = Dielectric_Elastomer_Visco

# \def DIRECTORY_DIELECTRIC_ELASTOMER_2D_VISCO
# 2D dielectric elastomer with viscoelastic material behavior.
# This option must be set in conjunction with #define DIELECTRIC_ELASTOMER_2D_VISCO
# in DevelopmentElementsConfig.h. */
#DIRECTORY_DIELECTRIC_ELASTOMER_2D_VISCO = Dielectric_Elastomer_2D_Visco

#############################################################
############## no configuration options below ###############
# Unless you are adding another option, you should not need #
# to edit the contents of below.                            #
#############################################################

DEVELOPMENT_ELEMENTS_CONFIG_DIRECTORIES = NONE \
	$(DIRECTORY_BEM_ELEMENT_DEV) \
	$(DIRECTORY_PML_ELEMENT_DEV) \
	$(DIRECTORY_COHESIVE_SURFACE_ELEMENT_DEV) \
	$(DIRECTORY_CONTACT_ELEMENT_DEV) \
	$(DIRECTORY_MULTISCALE_ELEMENT_DEV) \
	$(DIRECTORY_SOLID_ELEMENT_DEV) \
	$(DIRECTORY_MATERIAL_FORCE_ELEMENT_DEV) \
	$(DIRECTORY_GRAD_SMALL_STRAIN_DEV) \
	$(DIRECTORY_MULTISCALE_APS_DEV) \
	$(DIRECTORY_SPLIT_INTEGRATION_DEV) \
	$(DIRECTORY_SIMPLE_SOLID_DEV) \
	$(DIRECTORY_MULTISCALE_APS_V_DEV) \
	$(DIRECTORY_MESHFREE_GRAD_PLAST_DEV) \
	$(DIRECTORY_ENHANCED_STRAIN_LOC_DEV) \
	$(DIRECTORY_ENHANCED_STRAIN_LOC_DEV_CRAIG) \
	$(DIRECTORY_MIXTURE_THEORY_DEV) \
	$(DIRECTORY_BRIDGING_ELEMENT_DEV) \
	$(DIRECTORY_SURFACE_CB_DEV) \
	$(DIRECTORY_SURFACE_CB_EAM_DEV) \
	$(DIRECTORY_SURFACE_CB_SI_DEV) \
	$(DIRECTORY_SURFACE_CB_SI_DIMER_DEV) \
	$(DIRECTORY_SURFACE_CB_WURTZITE_DEV) \
	$(DIRECTORY_SURFACE_CB_ZB_DEV) \
	$(DIRECTORY_FLUID_ELEMENT_DEV) \
	$(DIRECTORY_FIBER_COMP_DEV) \
	$(DIRECTORY_SOLID_FLUID_MIX_DEV) \
	$(DIRECTORY_SOLID_OPTIMIZATION_DEV) \
	$(DIRECTORY_DEM_COUPLING_DEV) \
	$(DIRECTORY_PIEZOELECTRIC) \
	$(DIRECTORY_HUWASHIZU) \
	$(DIRECTORY_MICROMORPHIC_DEV) \
	$(DIRECTORY_MICROMORPHIC2_DEV) \
	$(DIRECTORY_MICROMORPHIC_CURR_CONFIG_DEV) \
	$(DIRECTORY_XFEM_DEV) \
	$(DIRECTORY_DIELECTRIC_ELASTOMER) \
	$(DIRECTORY_DIELECTRIC_ELASTOMER_Q1P0) \
	$(DIRECTORY_DIELECTRIC_ELASTOMER_Q1P0_ELASTOCAPILLARY) \
	$(DIRECTORY_DIELECTRIC_ELASTOMER_Q1P0_VISCO) \
	$(DIRECTORY_DIELECTRIC_ELASTOMER_2D) \
	$(DIRECTORY_DIELECTRIC_ELASTOMER_Q1P02D) \
	$(DIRECTORY_DIELECTRIC_ELASTOMER_VISCO) \
	$(DIRECTORY_DIELECTRIC_ELASTOMER_2D_VISCO) 
