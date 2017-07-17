# $Id$
# This file was generated by MakeConfigMakeFile.pl from ElementsConfig.conf
# created: Mon Jul 17 12:35:22 EDT 2017
#
# \file ElementsConfig.make
# Configuration of optional components within Tahoe.
# Sections of the code are included or excluded in the build of Tahoe depending in 
# this flags in this file and in the file ElementsConfig.h. Each option has
# a #define definition in ElementsConfig.h and a corresponding directory definition
# in this file. The two items must be set consistently to enable or
# disable materials models. To enable an option:
# -# in this file, uncomment the macro.
# -# in ElementsConfig.h, uncomment the #define statement
#
# The naming convention for the definitions in this file and the macros in
# ElementsConfig.h are as follows. For the option [OPTION]:
# -# the macro below defining the corresponding source directory will be DIRECTORY_[OPTION]
# -# the symbol in ElementsConfig.h file will be [OPTION]

# \def DIRECTORY_CONTACT_ELEMENT
# Methods for enforcing contact constraints.
# This option must be set in conjunction with #define CONTACT_ELEMENT
# in ElementsConfig.h. */
DIRECTORY_CONTACT_ELEMENT = contact

# \def DIRECTORY_COHESIVE_SURFACE_ELEMENT
# Cohesive surface elements.
# This option must be set in conjunction with #define COHESIVE_SURFACE_ELEMENT
# in ElementsConfig.h. */
DIRECTORY_COHESIVE_SURFACE_ELEMENT = cohesive_surface

# \def DIRECTORY_CONTINUUM_ELEMENT
# Continuum elements.
# This option must be set in conjunction with #define CONTINUUM_ELEMENT
# in ElementsConfig.h. */
DIRECTORY_CONTINUUM_ELEMENT = continuum

# \def DIRECTORY_SHAPE_FUNCTION_CLASSES
# Shape function classes.
# This option must be set in conjunction with #define SHAPE_FUNCTION_CLASSES
# in ElementsConfig.h. */
DIRECTORY_SHAPE_FUNCTION_CLASSES = shape_functions

# \def DIRECTORY_BRIDGING_ELEMENT
# Scale-bridging element, requires both CONTINUUM_ELEMENT and PARTICLE_ELEMENT.
# This option must be set in conjunction with #define BRIDGING_ELEMENT
# in ElementsConfig.h. */
#DIRECTORY_BRIDGING_ELEMENT = bridging_scale

# \def DIRECTORY_PARTICLE_ELEMENT
# Methods for particle simulations.
# This option must be set in conjunction with #define PARTICLE_ELEMENT
# in ElementsConfig.h. */
DIRECTORY_PARTICLE_ELEMENT = particle

# \def DIRECTORY_ADHESION_ELEMENT
# Surface-to-surface interaction.
# This option must be set in conjunction with #define ADHESION_ELEMENT
# in ElementsConfig.h. */
#DIRECTORY_ADHESION_ELEMENT = adhesion

# \def DIRECTORY_CONSTANT_VOLUME_ELEMENT
# Constant volume constraint elements.
# This option must be set in conjunction with #define CONSTANT_VOLUME_ELEMENT
# in ElementsConfig.h. */
#DIRECTORY_CONSTANT_VOLUME_ELEMENT = constant_volume

# \def DIRECTORY_SPRING_ELEMENT
# Linear and diamond cubic lattice elements.
# This option must be set in conjunction with #define SPRING_ELEMENT
# in ElementsConfig.h. */
DIRECTORY_SPRING_ELEMENT = spring

#############################################################
############## no configuration options below ###############
# Unless you are adding another option, you should not need #
# to edit the contents of below.                            #
#############################################################

ELEMENTS_CONFIG_DIRECTORIES = NONE \
	$(DIRECTORY_CONTACT_ELEMENT) \
	$(DIRECTORY_COHESIVE_SURFACE_ELEMENT) \
	$(DIRECTORY_CONTINUUM_ELEMENT) \
	$(DIRECTORY_SHAPE_FUNCTION_CLASSES) \
	$(DIRECTORY_BRIDGING_ELEMENT) \
	$(DIRECTORY_PARTICLE_ELEMENT) \
	$(DIRECTORY_ADHESION_ELEMENT) \
	$(DIRECTORY_CONSTANT_VOLUME_ELEMENT) \
	$(DIRECTORY_SPRING_ELEMENT) 
