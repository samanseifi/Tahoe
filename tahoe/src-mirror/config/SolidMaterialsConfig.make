# $Id$
# This file was generated by MakeConfigMakeFile.pl from SolidMaterialsConfig.conf
# created: Mon Sep  4 14:52:04 EDT 2017
#
# \file SolidMaterialsConfig.make
# Configuration of optional components within Tahoe.
# Sections of the code are included or excluded in the build of Tahoe depending in 
# this flags in this file and in the file SolidMaterialsConfig.h. Each option has
# a #define definition in SolidMaterialsConfig.h and a corresponding directory definition
# in this file. The two items must be set consistently to enable or
# disable materials models. To enable an option:
# -# in this file, uncomment the macro.
# -# in SolidMaterialsConfig.h, uncomment the #define statement
#
# The naming convention for the definitions in this file and the macros in
# SolidMaterialsConfig.h are as follows. For the option [OPTION]:
# -# the macro below defining the corresponding source directory will be DIRECTORY_[OPTION]
# -# the symbol in SolidMaterialsConfig.h file will be [OPTION]

# \def DIRECTORY_CAUCHY_BORN_MATERIAL
# Cauchy-Born crystal elasticity models.
# This option must be set in conjunction with #define CAUCHY_BORN_MATERIAL
# in SolidMaterialsConfig.h. */
DIRECTORY_CAUCHY_BORN_MATERIAL = EAM CB

# \def DIRECTORY_MODCBSW_MATERIAL
# modified Cauchy-Born crystal elasticity models using the Stillinger-Weber potentials.
# This option must be set in conjunction with #define MODCBSW_MATERIAL
# in SolidMaterialsConfig.h. */
DIRECTORY_MODCBSW_MATERIAL = modCBSW

# \def DIRECTORY_VIB_MATERIAL
# Virtual Internal Bond models.
# This option must be set in conjunction with #define VIB_MATERIAL
# in SolidMaterialsConfig.h. */
DIRECTORY_VIB_MATERIAL = VIB

# \def DIRECTORY_PLASTICITY_CRYSTAL_MATERIAL
# Local and non-local crystal plasticity models.
# This option must be set in conjunction with #define PLASTICITY_CRYSTAL_MATERIAL
# in SolidMaterialsConfig.h. */
DIRECTORY_PLASTICITY_CRYSTAL_MATERIAL = plasticity_crystl

# \def DIRECTORY_PLASTICITY_MACRO_MATERIAL
# Local and non-local crystal plasticity models.
# This option must be set in conjunction with #define PLASTICITY_MACRO_MATERIAL
# in SolidMaterialsConfig.h. */
DIRECTORY_PLASTICITY_MACRO_MATERIAL = plasticity_macro

# \def DIRECTORY_ELASTICITY_CRYSTAL_MATERIAL
# crystal elasticity models.
# This option must be set in conjunction with #define ELASTICITY_CRYSTAL_MATERIAL
# in SolidMaterialsConfig.h. */
DIRECTORY_ELASTICITY_CRYSTAL_MATERIAL = elasticity_crystl

# \def DIRECTORY_PLASTICITY_J2_MATERIAL
# J2 plasticity models.
# This option must be set in conjunction with #define PLASTICITY_J2_MATERIAL
# in SolidMaterialsConfig.h. */
DIRECTORY_PLASTICITY_J2_MATERIAL = plasticity_J2 plasticity_J2NlHard

# \def DIRECTORY_PLASTICITY_DP_MATERIAL
# Drucker-Prager plasticity models.
# This option must be set in conjunction with #define PLASTICITY_DP_MATERIAL
# in SolidMaterialsConfig.h. */
#DIRECTORY_PLASTICITY_DP_MATERIAL = plasticity_DP

# \def DIRECTORY_GEOMODEL_MATERIAL
# Sandia Geomodel plasticity model.
# This option must be set in conjunction with #define GEOMODEL_MATERIAL
# in SolidMaterialsConfig.h. */
DIRECTORY_GEOMODEL_MATERIAL = geomodel

# \def DIRECTORY_ABAQUS_MATERIAL
# ABAQUS UMAT's and VUMAT's (requires the f2c module).
# This option must be set in conjunction with #define ABAQUS_MATERIAL
# in SolidMaterialsConfig.h. */
#DIRECTORY_ABAQUS_MATERIAL = ABAQUS_UMAT ABAQUS_VUMAT

# \def DIRECTORY_THERMO_VISCO_PLASTIC_MATERIAL
# Thermo-visco-plastic models.
# This option must be set in conjunction with #define THERMO_VISCO_PLASTIC_MATERIAL
# in SolidMaterialsConfig.h. */
#DIRECTORY_THERMO_VISCO_PLASTIC_MATERIAL = tevp

# \def DIRECTORY_VISCOELASTICITY
# Viscoelasticity models.
# This option must be set in conjunction with #define VISCOELASTICITY
# in SolidMaterialsConfig.h. */
DIRECTORY_VISCOELASTICITY = viscoelasticity

# \def DIRECTORY_SIERRA_MATERIAL
# Models using the Sierra materials interface.
# This option must be set in conjunction with #define SIERRA_MATERIAL
# in SolidMaterialsConfig.h. */
#DIRECTORY_SIERRA_MATERIAL = SIERRA

#############################################################
############## no configuration options below ###############
# Unless you are adding another option, you should not need #
# to edit the contents of below.                            #
#############################################################

SOLID_MATERIALS_CONFIG_DIRECTORIES = NONE \
	$(DIRECTORY_CAUCHY_BORN_MATERIAL) \
	$(DIRECTORY_MODCBSW_MATERIAL) \
	$(DIRECTORY_VIB_MATERIAL) \
	$(DIRECTORY_PLASTICITY_CRYSTAL_MATERIAL) \
	$(DIRECTORY_PLASTICITY_MACRO_MATERIAL) \
	$(DIRECTORY_ELASTICITY_CRYSTAL_MATERIAL) \
	$(DIRECTORY_PLASTICITY_J2_MATERIAL) \
	$(DIRECTORY_PLASTICITY_DP_MATERIAL) \
	$(DIRECTORY_GEOMODEL_MATERIAL) \
	$(DIRECTORY_ABAQUS_MATERIAL) \
	$(DIRECTORY_THERMO_VISCO_PLASTIC_MATERIAL) \
	$(DIRECTORY_VISCOELASTICITY) \
	$(DIRECTORY_SIERRA_MATERIAL) 
