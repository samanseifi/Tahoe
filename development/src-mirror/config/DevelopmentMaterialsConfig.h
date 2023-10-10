/* $Id$ */
#ifndef _DEVELOPMENT_MATERIALS_CONFIG_H_
#define _DEVELOPMENT_MATERIALS_CONFIG_H_
/* This file was generated by MakeConfigHeaderFile.pl from DevelopmentMaterialsConfig.conf */
/* created: Sun Jan 12 19:24:10 EST 2020 */

/** \file DevelopmentMaterialsConfig.h
 * Configuration of optional components within Tahoe.
 * Sections of the code are included or excluded in the build of Tahoe depending in 
 * this flags in this file and in the file DevelopmentMaterialsConfig.make. Each option has
 * a #define definition in this file and a corresponding directory definition
 * in DevelopmentMaterialsConfig.make. The two items must be set consistently to enable or
 * disable materials models. To enable an option:
 * -# in this file, uncomment the #define statement
 * -# in DevelopmentMaterialsConfig.make, uncomment the macro.
 *
 * The naming convention for the definitions in this file and the macros in
 * DevelopmentMaterialsConfig.make are as follows. For the option [OPTION]:
 * -# the symbol in this file will be [OPTION]
 * -# the macro defining the corresponding source directory will be DIRECTORY_[OPTION]
 */

/** \def FOSSUM_MATERIAL_DEV
 * Cap plasticity model.
 * This option must be set in conjunction with the DIRECTORY_FOSSUM_MATERIAL_DEV macro
 * in DevelopmentMaterialsConfig.make. */
/* #define FOSSUM_MATERIAL_DEV 1 */

/** \def PLASTICITY_MR_MATERIAL_DEV
 * Geomaterial plasticity model.
 * This option must be set in conjunction with the DIRECTORY_PLASTICITY_MR_MATERIAL_DEV macro
 * in DevelopmentMaterialsConfig.make. */
/* #define PLASTICITY_MR_MATERIAL_DEV 1 */

/** \def PLASTICITY_SMR_MATERIAL_DEV
 * Simplified Geomaterial plasticity model.
 * This option must be set in conjunction with the DIRECTORY_PLASTICITY_SMR_MATERIAL_DEV macro
 * in DevelopmentMaterialsConfig.make. */
/* #define PLASTICITY_SMR_MATERIAL_DEV 1 */

/** \def GRAD_PLASTICITY_MR_MATERIAL_DEV
 * Geomaterial gradient plasticity model.
 * This option must be set in conjunction with the DIRECTORY_GRAD_PLASTICITY_MR_MATERIAL_DEV macro
 * in DevelopmentMaterialsConfig.make. */
/* #define GRAD_PLASTICITY_MR_MATERIAL_DEV 1 */

/** \def PLASTICITY_DP_LOC_MATERIAL_DEV
 * Drucker Prager plasticity model with localization.
 * This option must be set in conjunction with the DIRECTORY_PLASTICITY_DP_LOC_MATERIAL_DEV macro
 * in DevelopmentMaterialsConfig.make. */
/* #define PLASTICITY_DP_LOC_MATERIAL_DEV 1 */

/** \def VISCOELASTIC_MATERIALS_DEV
 * viscoelasticity models.
 * This option must be set in conjunction with the DIRECTORY_VISCOELASTIC_MATERIALS_DEV macro
 * in DevelopmentMaterialsConfig.make. */
#define VISCOELASTIC_MATERIALS_DEV 1

/** \def ABAQUS_DEV
 * Directory of ABAQUS UMAT's and VUMAT's.
 * This option must be set in conjunction with the DIRECTORY_ABAQUS_DEV macro
 * in DevelopmentMaterialsConfig.make. */
/* #define ABAQUS_DEV 1 */

/** \def ABAQUS_BCJ_MATERIAL_DEV
 * implicit and explicit ABAQUS BCJ models.
 * This option must be set in conjunction with the DIRECTORY_ABAQUS_BCJ_MATERIAL_DEV macro
 * in DevelopmentMaterialsConfig.make. */
/* #define ABAQUS_BCJ_MATERIAL_DEV 1 */

/** \def ABAQUS_XTAL_DEV
 * ABAQUS UMAT for xtalplasticity models.
 * This option must be set in conjunction with the DIRECTORY_ABAQUS_XTAL_DEV macro
 * in DevelopmentMaterialsConfig.make. */
/* #define ABAQUS_XTAL_DEV 1 */

/** \def ABAQUS_TI_DEV
 * ABAQUS UMAT model for titanium.
 * This option must be set in conjunction with the DIRECTORY_ABAQUS_TI_DEV macro
 * in DevelopmentMaterialsConfig.make. */
/* #define ABAQUS_TI_DEV 1 */

/** \def BIO_MODELS
 * models of biomaterials.
 * This option must be set in conjunction with the DIRECTORY_BIO_MODELS macro
 * in DevelopmentMaterialsConfig.make. */
#define BIO_MODELS 1

/** \def PIEZOELECTRIC
 * simple piezoelectric material.
 * This option must be set in conjunction with the DIRECTORY_PIEZOELECTRIC macro
 * in DevelopmentMaterialsConfig.make. */
/* #define PIEZOELECTRIC 1 */

/** \def DIELECTRIC_ELASTOMER
 * dielectric elastomer and electromechanical coupling in 3D.
 * This option must be set in conjunction with the DIRECTORY_DIELECTRIC_ELASTOMER macro
 * in DevelopmentMaterialsConfig.make. */
#define DIELECTRIC_ELASTOMER 1

/** \def DIELECTRIC_ELASTOMER_Q1P0
 * 3D dielectric elastomer with (hopefully) locking-free Q1P0 formulation.
 * This option must be set in conjunction with the DIRECTORY_DIELECTRIC_ELASTOMER_Q1P0 macro
 * in DevelopmentMaterialsConfig.make. */
#define DIELECTRIC_ELASTOMER_Q1P0 1

/** \def DIELECTRIC_ELASTOMER_Q1P0_ELASTOCAPILLARY
 * 3D dielectric elastomer with (hopefully) locking-free Q1P0 formulation + elastocapillary effects.
 * This option must be set in conjunction with the DIRECTORY_DIELECTRIC_ELASTOMER_Q1P0_ELASTOCAPILLARY macro
 * in DevelopmentMaterialsConfig.make. */
#define DIELECTRIC_ELASTOMER_Q1P0_ELASTOCAPILLARY 1

/** \def DIELECTRIC_ELASTOMER_Q1P0_VISCO
 * 3D dielectric elastomer with (hopefully) locking-free Q1P0 formulation and viscoelasticity.
 * This option must be set in conjunction with the DIRECTORY_DIELECTRIC_ELASTOMER_Q1P0_VISCO macro
 * in DevelopmentMaterialsConfig.make. */
#define DIELECTRIC_ELASTOMER_Q1P0_VISCO 1

/** \def DIELECTRIC_ELASTOMER_2D
 * dielectric elastomer and electromechanical coupling in 2D.
 * This option must be set in conjunction with the DIRECTORY_DIELECTRIC_ELASTOMER_2D macro
 * in DevelopmentMaterialsConfig.make. */
#define DIELECTRIC_ELASTOMER_2D 1

/** \def DIELECTRIC_ELASTOMER_Q1P02D
 * 2D dielectric elastomer with (hopefully) locking-free Q1P0 formulation.
 * This option must be set in conjunction with the DIRECTORY_DIELECTRIC_ELASTOMER_Q1P02D macro
 * in DevelopmentMaterialsConfig.make. */
#define DIELECTRIC_ELASTOMER_Q1P02D 1

/** \def DIELECTRIC_ELASTOMER_VISCO
 * 3D dielectric elastomer with viscoelastic material behavior.
 * This option must be set in conjunction with the DIRECTORY_DIELECTRIC_ELASTOMER_VISCO macro
 * in DevelopmentMaterialsConfig.make. */
/* #define DIELECTRIC_ELASTOMER_VISCO 1 */

/** \def DIELECTRIC_ELASTOMER_2D_VISCO
 * 2D dielectric elastomer with viscoelastic material behavior.
 * This option must be set in conjunction with the DIRECTORY_DIELECTRIC_ELASTOMER_2D_VISCO macro
 * in DevelopmentMaterialsConfig.make. */
/* #define DIELECTRIC_ELASTOMER_2D_VISCO 1 */

/** \def NEOHOOKEDAMAGE
 * simple elastic+damage material.
 * This option must be set in conjunction with the DIRECTORY_NEOHOOKEDAMAGE macro
 * in DevelopmentMaterialsConfig.make. */
/* #define NEOHOOKEDAMAGE 1 */

#endif
