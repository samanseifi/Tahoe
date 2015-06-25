/* $Id: FSSolidMatList2DT.cpp,v 1.23 2014/08/26 00:00:53 hspark Exp $ */
#include "FSSolidMatList2DT.h"
#include "FSMatSupportT.h"

#include "SolidMaterialsConfig.h"

#ifdef __DEVELOPMENT__
#include "DevelopmentMaterialsConfig.h"
#include "DevelopmentElementsConfig.h"
#endif

#include "FDHookeanMat2DT.h"
#include "FDKStV2D.h"
#include "FDCubic2DT.h"
#include "SimoIso2D.h"
#include "QuadLog2D.h"
#include "QuadLogOgden2DT.h"

#ifdef CAUCHY_BORN_MATERIAL
#include "EAMFCC2D.h"
#include "LJTr2D.h"
#include "Hex2D.h"
#endif

#ifdef MODCBSW_MATERIAL
#include "ModCB2DT.h"
#endif

#ifdef VIB_MATERIAL
#include "IsoVIB2D.h"
#include "J2IsoVIB2DLinHardT.h"
#include "VIB2D.h"
#include "D2VIB2D_a.h"
#include "OgdenIsoVIB2D.h"
#endif

#ifdef PLASTICITY_MACRO_MATERIAL
#include "HyperEVP2D.h"
#include "BCJHypo2D.h"
#include "BCJHypoIsoDamageKE2D.h"
#include "BCJHypoIsoDamageYC2D.h"
#endif

#ifdef PLASTICITY_CRYSTAL_MATERIAL
#include "LocalCrystalPlast2D.h"
#include "GradCrystalPlast2D.h"
#include "LocalCrystalPlastFp2D.h"
#include "GradCrystalPlastFp2D.h"
#endif

#ifdef PLASTICITY_J2_MATERIAL
#include "J2Simo2D.h"
#include "J2QL2DLinHardT.h"
#endif

#ifdef VISCOELASTICITY
#include "RGSplitT.h"
#include "RGSplitT2.h"
#include "BoyceViscoPlasticity.h"
#endif

#ifdef THERMO_VISCO_PLASTIC_MATERIAL
#include "tevp2D.h"
#include "povirk2D.h"
#endif

#ifdef FINITE_ANISOTROPY
#include "WLC.h"
#endif

/* development module materials require solid element development to be enabled */
#ifdef SOLID_ELEMENT_DEV

#ifdef VISCOELASTIC_MATERIALS_DEV
#include "FDSV_KStV2D.h"
#include "SMP_simple.h"
#include "SMP_multi.h"
//#include "SMP_solvent.h"
//#include "SMP_multisolvent.h"
#include "ModBoyceVisco.h"
#include "BergstromBoyce.h"
// #include "ElasticHydrogelT.h"
// #include "ElasticHydrogelSuo.h"
#endif

#ifdef DIELECTRIC_ELASTOMER_2D
#include "FSDielectricElastomer2DT.h"
#include "FSDEMatSupport2DT.h"
#endif

#ifdef DIELECTRIC_ELASTOMER_Q1P02D
#include "FSDielectricElastomerQ1P02DT.h"
#include "FSDEMatSupportQ1P02DT.h"
#endif

#ifdef DIELECTRIC_ELASTOMER_2D_VISCO
#include "FSDielectricElastomer2DViscoT.h"
#include "FSDEMatSupport2DViscoT.h"
#endif

#ifdef ABAQUS_MATERIAL
#ifdef ABAQUS_BCJ_MATERIAL_DEV
#include "ABAQUS_BCJ.h"
#include "ABAQUS_BCJ_ISO.h"
#include "ABAQUS_VUMAT_BCJ.h"
#endif /* ABAQUS_BCJ_MATERIAL_DEV */
#ifdef ABAQUS_TI_DEV
#include "ABAQUS_Ti.h"
#endif
#endif

#endif /* SOLID_ELEMENT_DEV */

using namespace Tahoe;

/* constructor */
FSSolidMatList2DT::FSSolidMatList2DT(int length, const FSMatSupportT& support):
	SolidMatListT(length, support),
	fFSMatSupport(&support)
{
	SetName("large_strain_material_2D");
}

FSSolidMatList2DT::FSSolidMatList2DT(void):
	fFSMatSupport(NULL)
{
	SetName("large_strain_material_2D");
}

/* return true if the list contains plane stress models */
bool FSSolidMatList2DT::HasPlaneStress(void) const
{
	/* check materials */
	for (int i = 0; i < Length(); i++)
	{
		/* get pointer to Material2DT */
		const ContinuumMaterialT* cont_mat = fArray[i];
		const SolidMaterialT* sol_mat = TB_DYNAMIC_CAST(const SolidMaterialT*, cont_mat);
		
		/* assume materials that don't have Material2DT are plane strain */
		if (sol_mat && sol_mat->Constraint() == SolidMaterialT::kPlaneStress) 
			return true;
	}
	return false;
}

/* information about subordinate parameter lists */
void FSSolidMatList2DT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SolidMatListT::DefineSubs(sub_list);
	
	/* choice of 2D materials */
	sub_list.AddSub("fs_material_list_2D", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void FSSolidMatList2DT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const
{
	if (name == "fs_material_list_2D")
	{
		order = ParameterListT::Choice;
	
		sub_lists.AddSub("large_strain_Hookean_2D");
		sub_lists.AddSub("large_strain_cubic_2D");
		sub_lists.AddSub("large_strain_StVenant_2D");
		sub_lists.AddSub("Simo_isotropic_2D");
		sub_lists.AddSub("quad_log_2D");
		sub_lists.AddSub("quad_log_Ogden_2D");

#ifdef PLASTICITY_J2_MATERIAL
		sub_lists.AddSub("Simo_J2_2D");
		sub_lists.AddSub("quad_log_J2_2D");
#endif

#ifdef PLASTICITY_CRYSTAL_MATERIAL
		sub_lists.AddSub("local_crystal_plasticity_2D");
		sub_lists.AddSub("local_crystal_plasticity_Fp_2D");
		sub_lists.AddSub("gradient_crystal_plasticity_Fp_2D");
#endif

#ifdef CAUCHY_BORN_MATERIAL
		sub_lists.AddSub("LJ_triangular_2D");
		sub_lists.AddSub("hex_2D");
		sub_lists.AddSub("FCC_EAM_2D");
#endif

#ifdef MODCBSW_MATERIAL
		sub_lists.AddSub("Cauchy-Born_diamond_2D");
#endif

#ifdef VIB_MATERIAL
		sub_lists.AddSub("VIB_2D");
		sub_lists.AddSub("isotropic_VIB_2D");
		sub_lists.AddSub("Ogden_isotropic_VIB_2D");
#endif

#ifdef VISCOELASTICITY
		sub_lists.AddSub("Reese-Govindjee_split");
		sub_lists.AddSub("RG_split_general");
		sub_lists.AddSub("boyce_viscoplasticity");
#endif
#ifdef VISCOELASTIC_MATERIALS_DEV
		sub_lists.AddSub("SMP_simple");
		sub_lists.AddSub("SMP_multi");
		//        sub_lists.AddSub("SMP_solvent");
		//		sub_lists.AddSub("SMP_multisolvent");
		sub_lists.AddSub("ModBoyceVisco");
		sub_lists.AddSub("BergstromBoyce");
//        sub_lists.AddSub("ElasticHydrogelT");
//		sub_lists.AddSub("ElasticHydrogelSuo");
#endif

#ifdef DIELECTRIC_ELASTOMER_2D
		sub_lists.AddSub(FSDEMat2DT::Name);
#endif

#ifdef DIELECTRIC_ELASTOMER_Q1P02D
		sub_lists.AddSub(FSDEMatQ1P02DT::Name);
#endif

#ifdef DIELECTRIC_ELASTOMER_2D_VISCO
		sub_lists.AddSub(FSDEMat2DViscoT::Name);
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
#endif
	}
	else /* inherited */
		SolidMatListT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FSSolidMatList2DT::NewSub(const StringT& name) const
{
	/* try to construct material */
	FSSolidMatT* fs_solid_mat = NewFSSolidMat(name);
	if (fs_solid_mat)
		return fs_solid_mat;
	else /* inherited */
		return SolidMatListT::NewSub(name);
}

/* accept parameter list */
void FSSolidMatList2DT::TakeParameterList(const ParameterListT& list)
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
FSSolidMatT* FSSolidMatList2DT::NewFSSolidMat(const StringT& name) const
{
	FSSolidMatT* mat = NULL;

	if (name == "large_strain_Hookean_2D")
		mat = new FDHookeanMat2DT;
	else if (name == "large_strain_cubic_2D")
		mat = new FDCubic2DT;
	else if (name == "large_strain_StVenant_2D")
		mat = new FDKStV2D;
	else if (name == "Simo_isotropic_2D")
		mat = new SimoIso2D;
	else if (name == "quad_log_2D")
		mat = new QuadLog2D;
	else if (name == "quad_log_Ogden_2D")
		mat = new QuadLogOgden2DT;

#ifdef PLASTICITY_J2_MATERIAL
	else if (name == "Simo_J2_2D")
		mat = new J2Simo2D;
	else if (name == "quad_log_J2_2D")
		mat = new J2QL2DLinHardT;
#endif

#ifdef PLASTICITY_CRYSTAL_MATERIAL
	else if (name == "local_crystal_plasticity_2D")
		mat = new LocalCrystalPlast2D;
	else if (name == "gradient_crystal_plasticity_Fp_2D")
		mat = new GradCrystalPlastFp2D;
	else if (name == "local_crystal_plasticity_Fp_2D")
		mat = new LocalCrystalPlastFp2D;
#endif

#ifdef CAUCHY_BORN_MATERIAL
	else if (name == "LJ_triangular_2D")
		mat = new LJTr2D;
	else if (name == "hex_2D")
		mat = new Hex2D;
	else if (name == "FCC_EAM_2D")
		mat = new EAMFCC2D;
#endif

#ifdef MODCBSW_MATERIAL
	else if (name == "Cauchy-Born_diamond_2D")
		mat = new ModCB2DT;
#endif

#ifdef VIB_MATERIAL
	else if (name == "VIB_2D")
		mat = new VIB2D;
	else if (name == "isotropic_VIB_2D")
		mat = new IsoVIB2D;
	else if (name == "Ogden_isotropic_VIB_2D")
		mat = new OgdenIsoVIB2D;
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
	//    else if (name == "SMP_solvent")
	//		mat= new SMP_solvent;
	//	else if (name == "SMP_multisolvent")
	//		mat= new SMP_multisolvent;
	else if (name == "ModBoyceVisco")
		mat= new ModBoyceVisco;
	else if (name == "BergstromBoyce")
		mat= new BergstromBoyce;
//   else if (name =="ElasticHydrogelT")
//		mat= new ElasticHydrogelT;
//    else if (name =="ElasticHydrogelSuo")
//		mat= new ElasticHydrogelSuo;
#endif

#ifdef FINITE_ANISOTROPY
	else if (name == "Bischoff-Arruda_WLC")
		mat= new WLC;
#endif

#ifdef DIELECTRIC_ELASTOMER_2D
	else if (name == FSDEMat2DT::Name) {
	  FSDEMat2DT* demat = new FSDEMat2DT;
	  if (demat != 0) {
	    FSMatSupportT* dems = const_cast<FSMatSupportT*>(fFSMatSupport);
	    const FSDEMatSupport2DT* ddems = dynamic_cast<FSDEMatSupport2DT*>(dems);
	    demat->SetFSDEMatSupport2D(ddems);
	    mat = demat;
	  }
	}
#endif

#ifdef DIELECTRIC_ELASTOMER_Q1P02D
	else if (name == FSDEMatQ1P02DT::Name) {
	  FSDEMatQ1P02DT* dematqp2d = new FSDEMatQ1P02DT;
	  if (dematqp2d != 0) {
	    FSMatSupportT* demqp2d = const_cast<FSMatSupportT*>(fFSMatSupport);
	    const FSDEMatSupportQ1P02DT* ddemqp2d = dynamic_cast<FSDEMatSupportQ1P02DT*>(demqp2d);
	    dematqp2d->SetFSDEMatSupportQ1P02D(ddemqp2d);
	    mat = dematqp2d;
	  }
	}
#endif

#ifdef DIELECTRIC_ELASTOMER_2D_VISCO
	else if (name == FSDEMat2DViscoT::Name) {
	  FSDEMat2DViscoT* dematv = new FSDEMat2DViscoT;
	  if (dematv != 0) {
	    FSMatSupportT* demsv = const_cast<FSMatSupportT*>(fFSMatSupport);
	    const FSDEMatSupport2DViscoT* ddemsv = dynamic_cast<FSDEMatSupport2DViscoT*>(demsv);
	    dematv->SetFSDEMatSupport2DVisco(ddemsv);
	    mat = dematv;
	  }
	}
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
#endif

	/* set support */
	if (mat) mat->SetFSMatSupport(fFSMatSupport);

	return mat;

}
