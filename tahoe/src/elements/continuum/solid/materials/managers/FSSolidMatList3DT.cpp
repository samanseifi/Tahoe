/* $Id: FSSolidMatList3DT.cpp,v 1.48 2013/11/27 18:52:48 xiaorui Exp $ */
/*Modified by RXiao 2013/11/26 by adding hydrogel and SMP-solvent  model*/
#include "FSSolidMatList3DT.h"

#include "SolidMaterialsConfig.h"


#ifdef __DEVELOPMENT__
#include "DevelopmentMaterialsConfig.h"
#include "DevelopmentElementsConfig.h"
#endif

#include "FDKStV.h"
#include "FDCubicT.h"
#include "QuadLog3D.h"
#include "SimoIso3D.h"
#include "QuadLogOgden3DT.h"

#ifdef CAUCHY_BORN_MATERIAL
#include "EAMFCC3DMatT.h"
#include "FCC3D.h"
#endif

#ifdef MODCBSW_MATERIAL
#include "ModCB3DT.h"
#endif

#ifdef VIB_MATERIAL
#include "VIB3D.h"
#include "IsoVIB3D.h"
#include "J2IsoVIB3DLinHardT.h"
#include "OgdenIsoVIB3D.h"
#endif

#ifdef VISCOELASTICITY
#include "RGSplitT.h"
#include "RGSplitT2.h"
#include "BoyceViscoPlasticity.h"
#endif

#ifdef FINITE_ANISOTROPY
#include "WLC.h"
#endif

#ifdef BIO_MODELS
#include "VerondaWestmannT.h"
#endif

#ifdef PLASTICITY_CRYSTAL_MATERIAL
#include "LocalCrystalPlast.h"
#include "LocalCrystalPlast_C.h"
#include "GradCrystalPlast.h"
#include "LocalCrystalPlastFp.h"
#include "LocalCrystalPlastFp_C.h"
#include "GradCrystalPlastFp.h"
#endif

#ifdef PLASTICITY_MACRO_MATERIAL
#include "HyperEVP3D.h"
#include "BCJHypo3D.h"
#include "BCJHypoIsoDamageKE3D.h"
#include "BCJHypoIsoDamageYC3D.h"
#endif

#ifdef PLASTICITY_J2_MATERIAL
#include "J2Simo3D.h"
#include "J2QLLinHardT.h"
#endif

#ifdef ELASTICITY_CRYSTAL_MATERIAL
#include "FDCrystalElast.h"
#endif

#ifdef THERMO_VISCO_PLASTIC_MATERIAL
#include "tevp3D.h"
#endif

#ifdef SIERRA_MATERIAL
#include "SIERRA_HypoElasticT.h"

#ifdef __FOSSUM__
#include "SIERRA_Isotropic_Geomaterial.h"
#endif /* __FOSSUM__ */

#ifdef __SIERRA__
#include "SIERRA_BCJ.h"
#include "SIERRA_EMMI.h"
#include "SIERRA_EMMI_iso.h"
#endif /* __SIERRA__ */

#endif /* SIERRA_MATERIAL */

#ifdef MIXTURE_THEORY_DEV
#include "FSSolidMixtureT.h"
#endif

#ifdef SURFACE_CB_SI_DEV
#include "CB_TersoffT.h"
#endif

#ifdef SURFACE_CB_SI_DIMER_DEV
#include "CB_TersoffDimerT.h"
#endif

#ifdef SURFACE_CB_WURTZITE_DEV
#include "CB_WurtziteT.h"
#endif

#ifdef SURFACE_CB_ZB_DEV
#include "CB_ZBT.h"
#endif

#ifdef ARRUDA_BOYCE
#include "Arruda_Boyce3D.h"
#endif

#ifdef PIEZOELECTRIC
#include "FSNeoHookePZLinT.h"
#include "FSPZMatSupportT.h"
#endif

#ifdef NEOHOOKEDAMAGE
#include "FSNeoHookeDamageT.h"
#endif

#ifdef DIELECTRIC_ELASTOMER
#include "FSDielectricElastomerT.h"
#include "FSDEMatSupportT.h"
#endif

#ifdef DIELECTRIC_ELASTOMER_Q1P0
#include "FSDielectricElastomerQ1P0T.h"
#include "FSDEMatSupportQ1P0T.h"
#endif

// #ifdef DIELECTRIC_ELASTOMER_Q1P0SURFACE
// #include "FSDielectricElastomerQ1P0SurfaceT.h"
// #include "FSDEMatSupportQ1P0SurfaceT.h"
// #endif

#ifdef DIELECTRIC_ELASTOMER_Q1P0_VISCO
#include "FSDielectricElastomerQ1P0ViscoT.h"
#include "FSDEMatSupportQ1P0ViscoT.h"
#endif

#ifdef DIELECTRIC_ELASTOMER_VISCO
#include "FSDielectricElastomerViscoT.h"
#include "FSDEMatSupportViscoT.h"
#endif

/* development module materials require solid element development to be enabled */
#ifdef SOLID_ELEMENT_DEV

#ifdef VISCOELASTIC_MATERIALS_DEV
#include "FDSV_KStV3D.h"
#include "SMP_simple.h"
#include "SMP_multi.h"
#include "SMP_solvent.h"
#include "SMP_multisolvent.h"
#include "ModBoyceVisco.h"
#include "BergstromBoyce.h"
#endif

#ifdef HYDROGEL
#include "ElasticHydrogelT.h"
#include "ElasticHydrogelSuo.h"
#include "ElasticHydrogelSuoT.h"
#endif

#ifdef ABAQUS_MATERIAL
#ifdef ABAQUS_BCJ_MATERIAL_DEV
#include "ABAQUS_BCJ.h"
#include "ABAQUS_BCJ_ISO.h"
#include "ABAQUS_VUMAT_BCJ.h"
#endif
#ifdef ABAQUS_TI_DEV
#include "ABAQUS_Ti.h"
#endif
#ifdef ABAQUS_XTAL_DEV
#include "ABAQUS_UMAT_Xtal.h"
#endif
#endif

#endif /* SOLID_ELEMENT_DEV */

using namespace Tahoe;

/* constructors */
FSSolidMatList3DT::FSSolidMatList3DT(int length, const FSMatSupportT& support):
	SolidMatListT(length, support),
	fFSMatSupport(&support)
{
	SetName("large_strain_material_3D");
}

FSSolidMatList3DT::FSSolidMatList3DT(void):
	fFSMatSupport(NULL)
{
	SetName("large_strain_material_3D");
}

/* information about subordinate parameter lists */
void FSSolidMatList3DT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SolidMatListT::DefineSubs(sub_list);

	/* choice of 2D materials */
	sub_list.AddSub("fs_material_list_3D", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void FSSolidMatList3DT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order,
		SubListT& sub_lists) const
{
	if (name == "fs_material_list_3D")
	{
		order = ParameterListT::Choice;

		sub_lists.AddSub("large_strain_Hookean");
		sub_lists.AddSub("large_strain_cubic");
		sub_lists.AddSub("large_strain_StVenant");
		sub_lists.AddSub("Simo_isotropic");
		sub_lists.AddSub("quad_log");
		sub_lists.AddSub("quad_log_Ogden");

#ifdef PLASTICITY_J2_MATERIAL
		sub_lists.AddSub("Simo_J2");
		sub_lists.AddSub("quad_log_J2");
#endif

#ifdef PLASTICITY_CRYSTAL_MATERIAL
		sub_lists.AddSub("local_crystal_plasticity");
		sub_lists.AddSub("local_crystal_plasticity_Fp");
		sub_lists.AddSub("local_crystal_plasticity_C");
		sub_lists.AddSub("local_crystal_plasticity_Fp_C");
		sub_lists.AddSub("gradient_crystal_plasticity_Fp");
#endif

#ifdef CAUCHY_BORN_MATERIAL
		sub_lists.AddSub("FCC_3D");
		sub_lists.AddSub("FCC_EAM");
#endif

#ifdef MODCBSW_MATERIAL
		sub_lists.AddSub("Cauchy-Born_diamond");
#endif

#ifdef VIB_MATERIAL
		sub_lists.AddSub("VIB");
		sub_lists.AddSub("isotropic_VIB");
		sub_lists.AddSub("Ogden_isotropic_VIB");
#endif

#ifdef VISCOELASTICITY
		sub_lists.AddSub("Reese-Govindjee_split");
		sub_lists.AddSub("RG_split_general");
                sub_lists.AddSub("boyce_viscoplasticity");
#endif

#ifdef VISCOELASTIC_MATERIALS_DEV
		sub_lists.AddSub("SMP_simple");
	    sub_lists.AddSub("SMP_multi");
		sub_lists.AddSub("SMP_solvent");
		sub_lists.AddSub("SMP_multisolvent");
		sub_lists.AddSub("ModBoyceVisco");
		sub_lists.AddSub("BergstromBoyce");	
#endif

#ifdef HYDROGEL
        sub_lists.AddSub("ElasticHydrogelT");
        sub_lists.AddSub("ElasticHydrogelSuo");
        sub_lists.AddSub("ElasticHydrogelSuoT");
#endif

#ifdef BIO_MODELS
		sub_lists.AddSub("veronda_westmann_potential");
#endif

#ifdef FINITE_ANISOTROPY
                sub_lists.AddSub("Bischoff-Arruda_WLC");
#endif

#ifdef ABAQUS_MATERIAL
#ifdef ABAQUS_BCJ_MATERIAL_DEV
		sub_lists.AddSub("ABAQUS_UMAT_BCJ");
		sub_lists.AddSub("ABAQUS_VUMAT_BCJ");
		sub_lists.AddSub("ABAQUS_UMAT_BCJ_iso-damage");
#endif
#ifdef ABAQUS_TI_DEV
		sub_lists.AddSub("ABAQUS_UMAT_Ti");
#endif
#ifdef ABAQUS_XTAL_DEV
		sub_lists.AddSub("ABAQUS_UMAT_Xtal");
#endif
#endif

#ifdef SIERRA_MATERIAL
		sub_lists.AddSub("SIERRA_hypoelastic");

#ifdef __SIERRA__
		sub_lists.AddSub("SIERRA_BCJ");
		sub_lists.AddSub("SIERRA_EMMI");
		sub_lists.AddSub("SIERRA_EMMI_iso");
#endif /* __SIERRA__ */

#endif

#ifdef MIXTURE_THEORY_DEV
		sub_lists.AddSub("large_strain_solid_mixture");
#endif

#ifdef SURFACE_CB_SI_DEV
		sub_lists.AddSub("Tersoff_CB");
#endif
#ifdef SURFACE_CB_SI_DIMER_DEV
		sub_lists.AddSub("TersoffDimer_CB");
#endif
#ifdef SURFACE_CB_WURTZITE_DEV
		sub_lists.AddSub("Wurtzite_CB");
#endif
#ifdef SURFACE_CB_ZB_DEV
		sub_lists.AddSub("ZB_CB");
#endif

#ifdef ARRUDA_BOYCE
		sub_lists.AddSub("Arruda_Boyce_3D");
#endif

#ifdef DIELECTRIC_ELASTOMER
		sub_lists.AddSub(FSDEMatT::Name);
#endif

#ifdef DIELECTRIC_ELASTOMER_Q1P0
		sub_lists.AddSub(FSDEMatQ1P0T::Name);
#endif

// #ifdef DIELECTRIC_ELASTOMER_Q1P0SURFACE
// 		sub_lists.AddSub(FSDEMatQ1P0SurfaceT::Name);
// #endif

#ifdef DIELECTRIC_ELASTOMER_Q1P0_VISCO
		sub_lists.AddSub(FSDEMatQ1P0ViscoT::Name);
#endif

#ifdef DIELECTRIC_ELASTOMER_VISCO
		sub_lists.AddSub(FSDEMatViscoT::Name);
#endif

#ifdef PIEZOELECTRIC
		sub_lists.AddSub(FSNeoHookePZLinT::Name);
#endif

#ifdef NEOHOOKEDAMAGE
    sub_lists.AddSub(FSNeoHookeDamageT::Name);
#endif

	}
	else /* inherited */
		SolidMatListT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FSSolidMatList3DT::NewSub(const StringT& name) const
{
	/* try to construct material */
	FSSolidMatT* fs_solid_mat = NewFSSolidMat(name);
	if (fs_solid_mat)
		return fs_solid_mat;
	else /* inherited */
		return SolidMatListT::NewSub(name);
}

/* accept parameter list */
void FSSolidMatList3DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SolidMatListT::TakeParameterList(list);

	/* construct materials - NOTE: subs have been defined as a choice, but
	 * here we construct as many materials as are passed in */
	const ArrayT<ParameterListT>& subs = list.Lists();
	int count = 0;
	for (int i = 0; i < subs.Length(); i++) {
		const ParameterListT& sub = subs[i];
		FSSolidMatT* mat = NewFSSolidMat(sub.Name());
		if (mat) {

			/* store pointer */
			(*this)[count++] = mat;

			/* initialize material */
			mat->TakeParameterList(sub);

			/* set flags */
			if (mat->HasHistory()) fHasHistory = true;
			if (mat->HasThermalStrain()) fHasThermal = true;
			if (mat->HasLocalization()) fHasLocalizers = true;
		}
	}
}

/* construct the specified material or NULL if the request cannot be completed */
FSSolidMatT* FSSolidMatList3DT::NewFSSolidMat(const StringT& name) const
{
	FSSolidMatT* mat = NULL;

	if (name == "large_strain_Hookean")
		mat = new FDHookeanMatT;
	else if (name == "large_strain_cubic")
		mat = new FDCubicT;
	else if (name == "large_strain_StVenant")
		mat = new FDKStV;
	else if (name == "Simo_isotropic")
		mat = new SimoIso3D;
	else if (name == "quad_log")
		mat = new QuadLog3D;
	else if (name == "quad_log_Ogden")
		mat = new QuadLogOgden3DT;

#ifdef PLASTICITY_J2_MATERIAL
	else if (name == "Simo_J2")
		mat = new J2Simo3D;
	else if (name == "quad_log_J2")
		mat = new J2QLLinHardT;
#endif

#ifdef PLASTICITY_CRYSTAL_MATERIAL
	else if (name == "local_crystal_plasticity")
		mat = new LocalCrystalPlast;
	else if (name == "local_crystal_plasticity_Fp")
		mat = new LocalCrystalPlastFp;
	else if (name == "gradient_crystal_plasticity_Fp")
		mat = new GradCrystalPlastFp;
	else if (name == "local_crystal_plasticity_C")
		mat = new LocalCrystalPlast_C;
	else if (name == "local_crystal_plasticity_Fp_C")
		mat = new LocalCrystalPlastFp_C;
#endif

#ifdef CAUCHY_BORN_MATERIAL
	else if (name == "FCC_3D")
		mat = new FCC3D;
	else if (name == "FCC_EAM")
		mat = new EAMFCC3DMatT;
#endif

#ifdef MODCBSW_MATERIAL
	else if (name == "Cauchy-Born_diamond")
		mat= new ModCB3DT;
#endif

#ifdef VIB_MATERIAL
	else if (name == "VIB")
		mat= new VIB3D;
	else if (name == "isotropic_VIB")
		mat= new IsoVIB3D;
	else if (name == "Ogden_isotropic_VIB")
		mat= new OgdenIsoVIB3D;
#endif

#ifdef VISCOELASTICITY
	else if (name == "Reese-Govindjee_split")
		mat= new RGSplitT;
	else if (name == "RG_split_general")
		mat= new RGSplitT2;
        else if (name == "boyce_viscoplasticity")
                mat = new BoyceViscoPlasticity;
#endif

#ifdef VISCOELASTIC_MATERIALS_DEV
	else if (name == "SMP_simple")
		mat= new SMP_simple;
	else if (name == "SMP_multi")
	   mat= new SMP_multi;
	else if (name == "SMP_solvent")
		mat= new SMP_solvent;
	else if (name == "SMP_multisolvent")
		mat= new SMP_multisolvent;
		/*add the following sentence*/
	else if (name == "ModBoyceVisco")
		mat= new ModBoyceVisco;
	else if (name == "BergstromBoyce")
		mat= new BergstromBoyce;
#endif

#ifdef HYDROGEL
       else if (name =="ElasticHydrogelT")
		mat= new ElasticHydrogelT;
        else if (name =="ElasticHydrogelSuo")
        mat= new ElasticHydrogelSuo;
        else if (name =="ElasticHydrogelSuoT")
            mat= new ElasticHydrogelSuoT;
#endif

#ifdef FINITE_ANISOTROPY
	else if (name == "Bischoff-Arruda_WLC")
	  mat= new WLC;
#endif

#ifdef ABAQUS_MATERIAL
#ifdef ABAQUS_BCJ_MATERIAL_DEV
	else if (name == "ABAQUS_UMAT_BCJ")
		mat= new ABAQUS_BCJ;
	else if (name == "ABAQUS_VUMAT_BCJ")
		mat= new ABAQUS_VUMAT_BCJ;
	else if (name == "ABAQUS_UMAT_BCJ_iso-damage")
		mat= new ABAQUS_BCJ_ISO;
#endif
#ifdef ABAQUS_TI_DEV
	else if (name == "ABAQUS_UMAT_Ti")
		mat= new ABAQUS_Ti;
#endif
#ifdef ABAQUS_XTAL_DEV
	else if (name == "ABAQUS_UMAT_Xtal")
		mat= new ABAQUS_UMAT_Xtal;
#endif
#endif

#ifdef SIERRA_MATERIAL
	else if (name == "SIERRA_hypoelastic")
		mat= new SIERRA_HypoElasticT;

#ifdef __SIERRA__
	else if (name == "SIERRA_BCJ")
		mat= new SIERRA_BCJ;
	else if (name == "SIERRA_EMMI")
		mat= new SIERRA_EMMI;
	else if (name == "SIERRA_EMMI_iso")
		mat= new SIERRA_EMMI_iso;
#endif /* __SIERRA__ */

#endif

#ifdef MIXTURE_THEORY_DEV
	else if (name == "large_strain_solid_mixture")
		mat= new FSSolidMixtureT;
#endif

#ifdef BIO_MODELS
	else if (name == "veronda_westmann_potential")
		mat = new VerondaWestmannT;
#endif

#ifdef SURFACE_CB_SI_DEV
	else if (name == "Tersoff_CB")
	  mat= new CB_TersoffT;
#endif
#ifdef SURFACE_CB_SI_DIMER_DEV
	else if (name == "TersoffDimer_CB")
	  mat= new CB_TersoffDimerT;
#endif
#ifdef SURFACE_CB_WURTZITE_DEV
	else if (name == "Wurtzite_CB")
	  mat= new CB_WurtziteT;
#endif
#ifdef SURFACE_CB_ZB_DEV
	else if (name == "ZB_CB")
	  mat= new CB_ZBT;
#endif

#ifdef ARRUDA_BOYCE
	else if (name == "Arruda_Boyce_3D")
	  mat= new Arruda_Boyce3D;
#endif

#ifdef DIELECTRIC_ELASTOMER
	else if (name == FSDEMatT::Name) {
	  FSDEMatT* demat = new FSDEMatT;
	  if (demat != 0) {
	    FSMatSupportT* dems = const_cast<FSMatSupportT*>(fFSMatSupport);
	    const FSDEMatSupportT* ddems = dynamic_cast<FSDEMatSupportT*>(dems);
	    demat->SetFSDEMatSupport(ddems);
	    mat = demat;
	  }
	}
#endif

#ifdef DIELECTRIC_ELASTOMER_Q1P0
	else if (name == FSDEMatQ1P0T::Name) {
	  FSDEMatQ1P0T* dematqp = new FSDEMatQ1P0T;
	  if (dematqp != 0) {
	    FSMatSupportT* demsqp = const_cast<FSMatSupportT*>(fFSMatSupport);
	    const FSDEMatSupportQ1P0T* ddemsqp = dynamic_cast<FSDEMatSupportQ1P0T*>(demsqp);
	    dematqp->SetFSDEMatSupportQ1P0(ddemsqp);
	    mat = dematqp;
	  }
	}
#endif

// #ifdef DIELECTRIC_ELASTOMER_Q1P0SURFACE
// 	else if (name == FSDEMatQ1P0SurfaceT::Name) {
// 	  FSDEMatQ1P0SurfaceT* dematqps = new FSDEMatQ1P0SurfaceT;
// 	  if (dematqps != 0) {
// 	    FSMatSupportT* demsqps = const_cast<FSMatSupportT*>(fFSMatSupport);
// 	    const FSDEMatSupportQ1P0SurfaceT* ddemsqps = dynamic_cast<FSDEMatSupportQ1P0SurfaceT*>(demsqps);
// 	    dematqps->SetFSDEMatSupportQ1P0Surface(ddemsqps);
// 	    mat = dematqps;
// 	  }
// 	}
// #endif

#ifdef DIELECTRIC_ELASTOMER_Q1P0_VISCO
	else if (name == FSDEMatQ1P0ViscoT::Name) {
	  FSDEMatQ1P0ViscoT* dematqpv = new FSDEMatQ1P0ViscoT;
	  if (dematqpv != 0) {
	    FSMatSupportT* demsqpv = const_cast<FSMatSupportT*>(fFSMatSupport);
	    const FSDEMatSupportQ1P0ViscoT* ddemsqpv = dynamic_cast<FSDEMatSupportQ1P0ViscoT*>(demsqpv);
	    dematqpv->SetFSDEMatSupportQ1P0Visco(ddemsqpv);
	    mat = dematqpv;
	  }
	}
#endif

#ifdef DIELECTRIC_ELASTOMER_VISCO
	else if (name == FSDEMatViscoT::Name) {
	  FSDEMatViscoT* dematvisco = new FSDEMatViscoT;
	  if (dematvisco != 0) {
	    FSMatSupportT* demsvisco = const_cast<FSMatSupportT*>(fFSMatSupport);
	    const FSDEMatSupportViscoT* ddemsvisco = dynamic_cast<FSDEMatSupportViscoT*>(demsvisco);
	    dematvisco->SetFSDEMatSupportVisco(ddemsvisco);
	    mat = dematvisco;
	  }
	}
#endif

#ifdef PIEZOELECTRIC
	else if (name == FSNeoHookePZLinT::Name) {
	  FSNeoHookePZLinT* pzmat = new FSNeoHookePZLinT;
	  if (pzmat != 0) {
	    FSMatSupportT* pms = const_cast<FSMatSupportT*>(fFSMatSupport);
	    const FSPZMatSupportT* ppzms = dynamic_cast<FSPZMatSupportT*>(pms);
	    pzmat->SetFSPZMatSupport(ppzms);
	    mat = pzmat;
	  }
	}
#endif

#ifdef NEOHOOKEDAMAGE
  else if (name == FSNeoHookeDamageT::Name) {
    mat = new FSNeoHookeDamageT;
  }
#endif

  /* set support */
	if (mat) mat->SetFSMatSupport(fFSMatSupport);

	return mat;
}
