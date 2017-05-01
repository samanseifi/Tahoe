# $Id$
# This file was generated by MakeConfigMakeFile.pl from DevelopmentMaterialsConfig.conf
# created: Mon May  1 17:36:17 EDT 2017
#
# \file DevelopmentMaterialsConfig.make
# Configuration of optional components within Tahoe.
# Sections of the code are included or excluded in the build of Tahoe depending in 
# this flags in this file and in the file DevelopmentMaterialsConfig.h. Each option has
# a #define definition in DevelopmentMaterialsConfig.h and a corresponding directory definition
# in this file. The two items must be set consistently to enable or
# disable materials models. To enable an option:
# -# in this file, uncomment the macro.
# -# in DevelopmentMaterialsConfig.h, uncomment the #define statement
#
# The naming convention for the definitions in this file and the macros in
# DevelopmentMaterialsConfig.h are as follows. For the option [OPTION]:
# -# the macro below defining the corresponding source directory will be DIRECTORY_[OPTION]
# -# the symbol in DevelopmentMaterialsConfig.h file will be [OPTION]

# \def DIRECTORY_FOSSUM_MATERIAL_DEV
# Cap plasticity model.
# This option must be set in conjunction with #define FOSSUM_MATERIAL_DEV
# in DevelopmentMaterialsConfig.h. */
#DIRECTORY_FOSSUM_MATERIAL_DEV = fossum

# \def DIRECTORY_PLASTICITY_MR_MATERIAL_DEV
# Geomaterial plasticity model.
# This option must be set in conjunction with #define PLASTICITY_MR_MATERIAL_DEV
# in DevelopmentMaterialsConfig.h. */
#DIRECTORY_PLASTICITY_MR_MATERIAL_DEV = plasticity_MR

# \def DIRECTORY_PLASTICITY_SMR_MATERIAL_DEV
# Simplified Geomaterial plasticity model.
# This option must be set in conjunction with #define PLASTICITY_SMR_MATERIAL_DEV
# in DevelopmentMaterialsConfig.h. */
#DIRECTORY_PLASTICITY_SMR_MATERIAL_DEV = plasticity_SMR

# \def DIRECTORY_GRAD_PLASTICITY_MR_MATERIAL_DEV
# Geomaterial gradient plasticity model.
# This option must be set in conjunction with #define GRAD_PLASTICITY_MR_MATERIAL_DEV
# in DevelopmentMaterialsConfig.h. */
#DIRECTORY_GRAD_PLASTICITY_MR_MATERIAL_DEV = grad_plasticity_MR

# \def DIRECTORY_PLASTICITY_DP_LOC_MATERIAL_DEV
# Drucker Prager plasticity model with localization.
# This option must be set in conjunction with #define PLASTICITY_DP_LOC_MATERIAL_DEV
# in DevelopmentMaterialsConfig.h. */
#DIRECTORY_PLASTICITY_DP_LOC_MATERIAL_DEV = plasticity_DP_loc

# \def DIRECTORY_VISCOELASTIC_MATERIALS_DEV
# viscoelasticity models.
# This option must be set in conjunction with #define VISCOELASTIC_MATERIALS_DEV
# in DevelopmentMaterialsConfig.h. */
DIRECTORY_VISCOELASTIC_MATERIALS_DEV = viscoelasticity

# \def DIRECTORY_ABAQUS_DEV
# Directory of ABAQUS UMAT's and VUMAT's.
# This option must be set in conjunction with #define ABAQUS_DEV
# in DevelopmentMaterialsConfig.h. */
#DIRECTORY_ABAQUS_DEV = ABAQUS

# \def DIRECTORY_ABAQUS_BCJ_MATERIAL_DEV
# implicit and explicit ABAQUS BCJ models.
# This option must be set in conjunction with #define ABAQUS_BCJ_MATERIAL_DEV
# in DevelopmentMaterialsConfig.h. */
#DIRECTORY_ABAQUS_BCJ_MATERIAL_DEV = ABAQUS_BCJ

# \def DIRECTORY_ABAQUS_XTAL_DEV
# ABAQUS UMAT for xtalplasticity models.
# This option must be set in conjunction with #define ABAQUS_XTAL_DEV
# in DevelopmentMaterialsConfig.h. */
#DIRECTORY_ABAQUS_XTAL_DEV = ABAQUS_Xtal

# \def DIRECTORY_ABAQUS_TI_DEV
# ABAQUS UMAT model for titanium.
# This option must be set in conjunction with #define ABAQUS_TI_DEV
# in DevelopmentMaterialsConfig.h. */
#DIRECTORY_ABAQUS_TI_DEV = jrmayeu

# \def DIRECTORY_BIO_MODELS
# models of biomaterials.
# This option must be set in conjunction with #define BIO_MODELS
# in DevelopmentMaterialsConfig.h. */
DIRECTORY_BIO_MODELS = bio_materials

# \def DIRECTORY_PIEZOELECTRIC
# simple piezoelectric material.
# This option must be set in conjunction with #define PIEZOELECTRIC
# in DevelopmentMaterialsConfig.h. */
#DIRECTORY_PIEZOELECTRIC = piezoelectricity

# \def DIRECTORY_DIELECTRIC_ELASTOMER
# dielectric elastomer and electromechanical coupling in 3D.
# This option must be set in conjunction with #define DIELECTRIC_ELASTOMER
# in DevelopmentMaterialsConfig.h. */
DIRECTORY_DIELECTRIC_ELASTOMER = dielectric_elastomer

# \def DIRECTORY_DIELECTRIC_ELASTOMER_Q1P0
# 3D dielectric elastomer with (hopefully) locking-free Q1P0 formulation.
# This option must be set in conjunction with #define DIELECTRIC_ELASTOMER_Q1P0
# in DevelopmentMaterialsConfig.h. */
DIRECTORY_DIELECTRIC_ELASTOMER_Q1P0 = dielectric_elastomer_Q1P0

# \def DIRECTORY_DIELECTRIC_ELASTOMER_Q1P0_ELASTOCAPILLARY
# 3D dielectric elastomer with (hopefully) locking-free Q1P0 formulation + elastocapillary effects.
# This option must be set in conjunction with #define DIELECTRIC_ELASTOMER_Q1P0_ELASTOCAPILLARY
# in DevelopmentMaterialsConfig.h. */
DIRECTORY_DIELECTRIC_ELASTOMER_Q1P0_ELASTOCAPILLARY = dielectric_elastomer_Q1P0_Elastocapillary

# \def DIRECTORY_DIELECTRIC_ELASTOMER_Q1P0_VISCO
# 3D dielectric elastomer with (hopefully) locking-free Q1P0 formulation and viscoelasticity.
# This option must be set in conjunction with #define DIELECTRIC_ELASTOMER_Q1P0_VISCO
# in DevelopmentMaterialsConfig.h. */
DIRECTORY_DIELECTRIC_ELASTOMER_Q1P0_VISCO = dielectric_elastomer_Q1P0_visco

# \def DIRECTORY_DIELECTRIC_ELASTOMER_2D
# dielectric elastomer and electromechanical coupling in 2D.
# This option must be set in conjunction with #define DIELECTRIC_ELASTOMER_2D
# in DevelopmentMaterialsConfig.h. */
DIRECTORY_DIELECTRIC_ELASTOMER_2D = dielectric_elastomer_2D

# \def DIRECTORY_DIELECTRIC_ELASTOMER_Q1P02D
# 2D dielectric elastomer with (hopefully) locking-free Q1P0 formulation.
# This option must be set in conjunction with #define DIELECTRIC_ELASTOMER_Q1P02D
# in DevelopmentMaterialsConfig.h. */
DIRECTORY_DIELECTRIC_ELASTOMER_Q1P02D = dielectric_elastomer_Q1P02D

# \def DIRECTORY_DIELECTRIC_ELASTOMER_VISCO
# 3D dielectric elastomer with viscoelastic material behavior.
# This option must be set in conjunction with #define DIELECTRIC_ELASTOMER_VISCO
# in DevelopmentMaterialsConfig.h. */
#DIRECTORY_DIELECTRIC_ELASTOMER_VISCO = dielectric_elastomer_visco

# \def DIRECTORY_DIELECTRIC_ELASTOMER_2D_VISCO
# 2D dielectric elastomer with viscoelastic material behavior.
# This option must be set in conjunction with #define DIELECTRIC_ELASTOMER_2D_VISCO
# in DevelopmentMaterialsConfig.h. */
#DIRECTORY_DIELECTRIC_ELASTOMER_2D_VISCO = dielectric_elastomer_2D_visco

# \def DIRECTORY_NEOHOOKEDAMAGE
# simple elastic+damage material.
# This option must be set in conjunction with #define NEOHOOKEDAMAGE
# in DevelopmentMaterialsConfig.h. */
#DIRECTORY_NEOHOOKEDAMAGE = damage

#############################################################
############## no configuration options below ###############
# Unless you are adding another option, you should not need #
# to edit the contents of below.                            #
#############################################################

DEVELOPMENT_MATERIALS_CONFIG_DIRECTORIES = NONE \
	$(DIRECTORY_FOSSUM_MATERIAL_DEV) \
	$(DIRECTORY_PLASTICITY_MR_MATERIAL_DEV) \
	$(DIRECTORY_PLASTICITY_SMR_MATERIAL_DEV) \
	$(DIRECTORY_GRAD_PLASTICITY_MR_MATERIAL_DEV) \
	$(DIRECTORY_PLASTICITY_DP_LOC_MATERIAL_DEV) \
	$(DIRECTORY_VISCOELASTIC_MATERIALS_DEV) \
	$(DIRECTORY_ABAQUS_DEV) \
	$(DIRECTORY_ABAQUS_BCJ_MATERIAL_DEV) \
	$(DIRECTORY_ABAQUS_XTAL_DEV) \
	$(DIRECTORY_ABAQUS_TI_DEV) \
	$(DIRECTORY_BIO_MODELS) \
	$(DIRECTORY_PIEZOELECTRIC) \
	$(DIRECTORY_DIELECTRIC_ELASTOMER) \
	$(DIRECTORY_DIELECTRIC_ELASTOMER_Q1P0) \
	$(DIRECTORY_DIELECTRIC_ELASTOMER_Q1P0_ELASTOCAPILLARY) \
	$(DIRECTORY_DIELECTRIC_ELASTOMER_Q1P0_VISCO) \
	$(DIRECTORY_DIELECTRIC_ELASTOMER_2D) \
	$(DIRECTORY_DIELECTRIC_ELASTOMER_Q1P02D) \
	$(DIRECTORY_DIELECTRIC_ELASTOMER_VISCO) \
	$(DIRECTORY_DIELECTRIC_ELASTOMER_2D_VISCO) \
	$(DIRECTORY_NEOHOOKEDAMAGE) 
