/* $Id$ */
#ifndef _ELEMENTS_CONFIG_H_
#define _ELEMENTS_CONFIG_H_
/* This file was generated by MakeConfigHeaderFile.pl from ElementsConfig.conf */
/* created: Mon Aug 28 18:23:33 EDT 2017 */

/** \file ElementsConfig.h
 * Configuration of optional components within Tahoe.
 * Sections of the code are included or excluded in the build of Tahoe depending in 
 * this flags in this file and in the file ElementsConfig.make. Each option has
 * a #define definition in this file and a corresponding directory definition
 * in ElementsConfig.make. The two items must be set consistently to enable or
 * disable materials models. To enable an option:
 * -# in this file, uncomment the #define statement
 * -# in ElementsConfig.make, uncomment the macro.
 *
 * The naming convention for the definitions in this file and the macros in
 * ElementsConfig.make are as follows. For the option [OPTION]:
 * -# the symbol in this file will be [OPTION]
 * -# the macro defining the corresponding source directory will be DIRECTORY_[OPTION]
 */

/** \def CONTACT_ELEMENT
 * Methods for enforcing contact constraints.
 * This option must be set in conjunction with the DIRECTORY_CONTACT_ELEMENT macro
 * in ElementsConfig.make. */
#define CONTACT_ELEMENT 1

/** \def COHESIVE_SURFACE_ELEMENT
 * Cohesive surface elements.
 * This option must be set in conjunction with the DIRECTORY_COHESIVE_SURFACE_ELEMENT macro
 * in ElementsConfig.make. */
#define COHESIVE_SURFACE_ELEMENT 1

/** \def CONTINUUM_ELEMENT
 * Continuum elements.
 * This option must be set in conjunction with the DIRECTORY_CONTINUUM_ELEMENT macro
 * in ElementsConfig.make. */
#define CONTINUUM_ELEMENT 1

/** \def SHAPE_FUNCTION_CLASSES
 * Shape function classes.
 * This option must be set in conjunction with the DIRECTORY_SHAPE_FUNCTION_CLASSES macro
 * in ElementsConfig.make. */
#define SHAPE_FUNCTION_CLASSES 1

/** \def BRIDGING_ELEMENT
 * Scale-bridging element, requires both CONTINUUM_ELEMENT and PARTICLE_ELEMENT.
 * This option must be set in conjunction with the DIRECTORY_BRIDGING_ELEMENT macro
 * in ElementsConfig.make. */
/* #define BRIDGING_ELEMENT 1 */

/** \def PARTICLE_ELEMENT
 * Methods for particle simulations.
 * This option must be set in conjunction with the DIRECTORY_PARTICLE_ELEMENT macro
 * in ElementsConfig.make. */
#define PARTICLE_ELEMENT 1

/** \def ADHESION_ELEMENT
 * Surface-to-surface interaction.
 * This option must be set in conjunction with the DIRECTORY_ADHESION_ELEMENT macro
 * in ElementsConfig.make. */
/* #define ADHESION_ELEMENT 1 */

/** \def CONSTANT_VOLUME_ELEMENT
 * Constant volume constraint elements.
 * This option must be set in conjunction with the DIRECTORY_CONSTANT_VOLUME_ELEMENT macro
 * in ElementsConfig.make. */
/* #define CONSTANT_VOLUME_ELEMENT 1 */

/** \def SPRING_ELEMENT
 * Linear and diamond cubic lattice elements.
 * This option must be set in conjunction with the DIRECTORY_SPRING_ELEMENT macro
 * in ElementsConfig.make. */
#define SPRING_ELEMENT 1

#endif
