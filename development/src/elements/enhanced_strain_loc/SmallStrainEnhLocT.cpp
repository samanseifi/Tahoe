/* $Id: SmallStrainEnhLocT.cpp,v 1.51 2010/07/22 15:17:41 regueiro Exp $ */
#include "SmallStrainEnhLocT.h"
#include "ShapeFunctionT.h"
#include "SSSolidMatT.h"

/* materials lists */
#include "SSSolidMatList1DT.h"
#include "SSSolidMatList2DT.h"
#include "SSSolidMatList3DT.h"
#include "SSEnhLocMatSupportT.h"
#include "ParameterContainerT.h"
#include "ModelManagerT.h"

#include "ParameterUtils.h"

using namespace Tahoe;

/* parameters */
const double Pi = acos(-1.0);
const double smallnum = 1.0e-4;
const double verysmallnum = 1.0e-40;

/* initialize static variables */
bool SmallStrainEnhLocT::fFirstPass = true;

#ifndef __MWERKS__ // for compilation outside CodeWarrior
#ifdef NDEBUG
bool SmallStrainEnhLocT::fDeBug = false; // no output
#else
bool SmallStrainEnhLocT::fDeBug = true; // output info for debugging
#endif
#else
bool SmallStrainEnhLocT::fDeBug = true; // output info for debugging
#endif // __MWERKS__

//bool SmallStrainEnhLocT::fFirstTrace = false;

/* constructor */
SmallStrainEnhLocT::SmallStrainEnhLocT(const ElementSupportT& support):
	SolidElementT(support),
	fNeedsOffset(-1),
	fSSEnhLocMatSupport(NULL)
{
	SetName("small_strain_enh_loc");
}

/* destructor */
SmallStrainEnhLocT::~SmallStrainEnhLocT(void)
{
	delete fSSEnhLocMatSupport;
}


/* initialize current step */
void SmallStrainEnhLocT::InitStep(void)
{
	/* inherited */
	SolidElementT::InitStep();
	
	/* initialize to converged previous solution */
	fElementLocScalars = fElementLocScalars_last;
	fElementLocSlipDir = fElementLocSlipDir_last;
	fElementLocMuDir = fElementLocMuDir_last;
	fElementLocInternalVars = fElementLocInternalVars_last;
	fElementVolume = fElementVolume_last;
	
	fElementStress = fElementStress_last;
	
	if (fDeBug)
	{
		ss_enh_out	<< endl << "**********************************************************************************************";
		ss_enh_out	<< endl << "**********************************************************************************************" << endl;
		ss_enh_out	<< endl 
					<< setw(outputFileWidth) << "time_step"
					<< endl;
		int step = fSSEnhLocMatSupport->StepNumber();
		ss_enh_out	<< setw(outputFileWidth) << step
					<< endl;	
	}
	
	// initialize coupling stiffness matrices
	fElementLocKzetad = 0.0;
	fElementLocKdzeta = 0.0;
	int elem_num;
	Top();
	while (NextElement())
	{
		elem_num = CurrElementNumber();
		fElementLocScalars[kNUM_SCALAR_TERMS*elem_num + kKzetazeta] = 0.0;
		//fElementLocScalars[kNUM_SCALAR_TERMS*elem_num + kr_S] = 0.0;
	}
		
}


/* finalize current step - step is solved */
void SmallStrainEnhLocT::CloseStep(void)
{
	/* inherited */
	SolidElementT::CloseStep();
	
	int nen, elem_num;
	
	Top();
	while (NextElement())
	{
		/* shape function derivatives, jacobians, local coords */
		SetGlobalShape();
		
		/* get displacements */
		SetLocalU(fLocDisp);
        
		/* get reference geometry */
		SetLocalX(fLocInitCoords);
		/* get current geometry */
		//SetLocalX(fLocCurrCoords); // same for small strain
		
		elem_num = CurrElementNumber();
		nen = NumElementNodes();
		loc_flag = fElementLocFlag[elem_num];

		if ( loc_flag == 1 )
		{
			/* determine active nodes and band trace */
			DetermineActiveNodesTrace(fLocInitCoords, elem_num, nen);
		}
		
		/*
		if ( fDeBug ) ss_enh_out << elem_num 
								<< setw(outputFileWidth) << fElementLocFlag[elem_num]
								<< endl;
								*/
		
		if ( fDeBug &&  loc_flag == 2 )
		{
			ss_enh_isv	<< elem_num << setw(outputFileWidth) << fElementLocFlag[elem_num] 
						<< setw(outputFileWidth) << fElementLocScalars[kNUM_SCALAR_TERMS*elem_num + kzeta] 
						<< setw(outputFileWidth) << fElementLocScalars[kNUM_SCALAR_TERMS*elem_num + kgamma_delta] 
						<< setw(outputFileWidth) << fElementLocScalars[kNUM_SCALAR_TERMS*elem_num + kQ_S]
						<< setw(outputFileWidth) << fElementLocScalars[kNUM_SCALAR_TERMS*elem_num + kP_S] 
						<< setw(outputFileWidth) << fElementLocScalars[kNUM_SCALAR_TERMS*elem_num + kq_St]
			 			<< setw(outputFileWidth) << fElementLocInternalVars[kNUM_ISV_TERMS*elem_num + kCohesion] 
			 			<< setw(outputFileWidth) << fElementLocInternalVars[kNUM_ISV_TERMS*elem_num + kFriction] 
						<< setw(outputFileWidth) << fElementLocInternalVars[kNUM_ISV_TERMS*elem_num + kDilation]
						<< endl;
		}		
		
	} // while next element
	
	/* store converged solution */
	fElementLocScalars_last = fElementLocScalars;
	fElementLocSlipDir_last = fElementLocSlipDir;
	fElementLocMuDir_last = fElementLocMuDir;
	fElementLocInternalVars_last = fElementLocInternalVars;
	fElementVolume_last = fElementVolume;
	
	fElementStress_last = fElementStress;
}
	
/* restore last converged state */
GlobalT::RelaxCodeT SmallStrainEnhLocT::ResetStep(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = SolidElementT::ResetStep();
	
	/* restore converged solution */
	fElementLocScalars = fElementLocScalars_last;
	fElementLocSlipDir = fElementLocSlipDir_last;
	fElementLocMuDir = fElementLocMuDir_last;
	fElementLocInternalVars = fElementLocInternalVars_last;
	fElementVolume = fElementVolume_last;
	
	fElementStress = fElementStress_last;
	
	return relax;
}

/* read restart information from stream */
void SmallStrainEnhLocT::ReadRestart(istream& in)
{
	/* inherited */
	SolidElementT::ReadRestart(in);
	
	/* read restart data */
	in >> fElementLocScalars;
	in >> fElementLocSlipDir;
	in >> fElementLocMuDir;
	in >> fElementLocInternalVars;
	in >> fElementVolume;
	
	in >> fElementStress;
	
	/* reset last state */
	fElementLocScalars_last = fElementLocScalars;
	fElementLocSlipDir_last = fElementLocSlipDir;
	fElementLocMuDir_last = fElementLocMuDir;
	fElementLocInternalVars_last = fElementLocInternalVars;
	fElementVolume_last = fElementVolume;
	
	fElementStress_last = fElementStress;
}

/* write restart information from stream */
void SmallStrainEnhLocT::WriteRestart(ostream& out) const
{
	/* inherited */
	SolidElementT::WriteRestart(out);
	
	/* write restart data */
	out << fElementLocScalars << '\n';
	out << fElementLocSlipDir << '\n';
	out << fElementLocMuDir << '\n';
	out << fElementLocInternalVars << '\n';
	out << fElementVolume << '\n';
	out << fElementStress << '\n';
}

/* check for localization at end of step */
GlobalT::RelaxCodeT SmallStrainEnhLocT::RelaxSystem(void)
{
	/* inherited */
	GlobalT::RelaxCodeT code = SolidElementT::RelaxSystem();
	
	int nen, elem_num;
	
	Top();
	while (NextElement())
	{
		/* shape function derivatives, jacobians, local coords */
		SetGlobalShape();
		
		/* get displacements */
        SetLocalU(fLocDisp);
		
		elem_num = CurrElementNumber();
		nen = NumElementNodes();
		loc_flag = fElementLocFlag[elem_num];
		const ElementCardT& element = CurrentElement();

		/* until element is traced, i.e. loc_flag = 2, check for localization */
		if ( loc_flag < 2 && element.IsAllocated() )
		{
			/* initialize element localization data */
			fElementLocFlag[elem_num] = 0;
			
			fElementLocScalars[kNUM_SCALAR_TERMS*elem_num + kdetAmin] = 1.0e99;
			
			fElementLocScalars[kNUM_SCALAR_TERMS*elem_num + kdissip_max] = 0.0;
			
			normal_chosen = 0.0;
			fElementLocNormal.SetRow(elem_num, normal_chosen);
			slipdir_chosen = 0.0;
			fElementLocSlipDir.SetRow(elem_num, slipdir_chosen);
			tangent_chosen = 0.0;
			fElementLocTangent.SetRow(elem_num, tangent_chosen);

			/* check for localization */
			CheckLocalization(elem_num, fLocDisp);
			
			/* choose normal and slipdir */
			loc_flag = fElementLocFlag[elem_num];
			if ( loc_flag == 1 )
			{
				/* choose normal and slip direction based on current element deformation */
				ChooseNormalAndSlipDir(fLocDisp, elem_num, nen);
			}
	
		}
	
	} // while next element

	return code;
}



/* implementation of the ParameterInterfaceT interface */
void SmallStrainEnhLocT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SolidElementT::DefineParameters(list);

	/* strain-displacement relation */
	ParameterT strain_displacement(ParameterT::Enumeration, "strain_displacement");
	strain_displacement.AddEnumeration("standard", kStandardB);
    strain_displacement.AddEnumeration("B-bar", kMeanDilBbar);
    strain_displacement.SetDefault(kStandardB);
	list.AddParameter(strain_displacement);
	
	// select model type
	int modeltype;
	list.AddParameter(modeltype, "model_type");
	
	/* cohesive surface model constants */
	// eventually put this in a separate class to try different cohesive surface models
	double cohesion_r, cohesion_p, alpha_c, phi_r, phi_p, alpha_phi, psi_p, alpha_psi;
	list.AddParameter(cohesion_r, "residual_cohesion");
	list.AddParameter(cohesion_p, "peak_cohesion");
	list.AddParameter(alpha_c, "cohesion_softening_coefficient");
	list.AddParameter(phi_r, "residual_friction_angle_rad");
	list.AddParameter(phi_p, "peak_friction_angle_rad");
	list.AddParameter(alpha_phi, "friction_softening_coefficient");
	list.AddParameter(psi_p, "peak_dilation_angle_rad");
	list.AddParameter(alpha_psi, "dilation_softening_coefficient");
	
	int choosenormal;
	list.AddParameter(choosenormal, "choose_normal");
	
	int chooseelement;
	list.AddParameter(chooseelement, "choose_element");	
}

/* information about subordinate parameter lists */
void SmallStrainEnhLocT::DefineSubs(SubListT& sub_list) const
{	
	/* inherited */
	SolidElementT::DefineSubs(sub_list);

	/* element block/material specification */
	sub_list.AddSub("small_strain_enh_loc_element_block", ParameterListT::OnePlus);
	
	/* coordinates of point from which to start discontinuity surface */
	sub_list.AddSub("start_surface_coord", ParameterListT::ZeroOrOnce);
}

/* return the description of the given inline subordinate parameter list */
ParameterInterfaceT* SmallStrainEnhLocT::NewSub(const StringT& name) const
{
	if (name == "small_strain_enh_loc_element_block")
	{
		ParameterContainerT* block = new ParameterContainerT(name);
		
		/* list of element block ID's (defined by ElementBaseT) */
		block->AddSub("block_ID_list", ParameterListT::Once);
	
		/* choice of materials lists (inline) */
		block->AddSub("small_strain_enh_loc_material_choice", ParameterListT::Once, true);
	
		/* set this as source of subs */
		block->SetSubSource(this);
		
		return block;
	}
	else if (name == "start_surface_coord") {

		ParameterContainerT* s_choice = new ParameterContainerT(name);
		
		/* by dimension */
		s_choice->SetListOrder(ParameterListT::Choice);
		s_choice->AddSub("Vector_2");
		s_choice->AddSub("Vector_3");

		return s_choice;	
	}
	else /* inherited */
		return SolidElementT::NewSub(name);
}

/* return the description of the given inline subordinate parameter list. */
void SmallStrainEnhLocT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	if (name == "small_strain_enh_loc_material_choice")
	{
		order = ParameterListT::Choice;
		
		/* list of choices */
		sub_lists.AddSub("small_strain_material_1D");
		sub_lists.AddSub("small_strain_material_2D");
		sub_lists.AddSub("small_strain_material_3D");
	}
	else /* inherited */
		SolidElementT::DefineInlineSub(name, order, sub_lists);
}

void SmallStrainEnhLocT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "SmallStrainEnhLocT::TakeParameterList";

	/* strain displacement option before calling SolidElementT::TakeParameterList */
	int b = list.GetParameter("strain_displacement");
	fStrainDispOpt = (b == kStandardB) ? kStandardB : kMeanDilBbar;

	/* inherited */
	SolidElementT::TakeParameterList(list);
	
	/* get coordinates of start surface */
	const ParameterListT* start_check = list.List("start_surface_coord");
	if (start_check)
	{
		int nsd = NumSD();
		const ParameterListT& start_surf = list.GetListChoice(*this, "start_surface_coord");
		VectorParameterT::Extract(start_surf, start_surface_coord_read);
		if (start_surface_coord_read.Length() != nsd) 
			ExceptionT::GeneralFail(caller, "\"start_surface_coord\" should be length %d not %d", nsd, start_surface_coord_read.Length());
	}
	
	/* dimension workspace */
	fGradU.Dimension(NumSD());	
	if (fStrainDispOpt == kMeanDilBbar) {
		fLocDispTranspose.Dimension(fLocDisp.Length());
		//fMeanGradient.Dimension(NumSD(), NumElementNodes());
	}	
	fMeanGradient.Dimension(NumSD(), NumElementNodes());

	/* offset to class needs flags */
	fNeedsOffset = fMaterialNeeds[0].Length();
	
	/* set material needs */
	for (int i = 0; i < fMaterialNeeds.Length(); i++)
	{
		/* needs array */
		ArrayT<bool>& needs = fMaterialNeeds[i];

		/* resize array */
		needs.Resize(needs.Length() + 2, true);

		/* casts are safe since class contructs materials list */
		ContinuumMaterialT* pcont_mat = (*fMaterialList)[i];
		SSSolidMatT* mat = (SSSolidMatT*) pcont_mat;

		/* collect needs */
		needs[fNeedsOffset + kstrain     ] = mat->Need_Strain();
		needs[fNeedsOffset + kstrain_last] = mat->Need_Strain_last();
		
		/* consistency */
		needs[kNeedDisp] = needs[kNeedDisp] || needs[fNeedsOffset + kstrain];
		needs[KNeedLastDisp] = needs[KNeedLastDisp] || needs[fNeedsOffset + kstrain_last];
	}

	/* what's needed */
	bool need_strain = false;
	bool need_strain_last = false;
	for (int i = 0; i < fMaterialNeeds.Length(); i++) {
		const ArrayT<bool>& needs = fMaterialNeeds[i];
		need_strain = need_strain || needs[fNeedsOffset + kstrain];
		need_strain_last = need_strain_last || needs[fNeedsOffset + kstrain_last];
	}

	/* allocate strain list */
	if (need_strain) {
		fStrain_List.Dimension(NumIP());
		for (int i = 0; i < NumIP(); i++)
			fStrain_List[i].Dimension(NumSD());
	}
	
	/* allocate "last" strain list */
	if (need_strain_last) {
		fStrain_last_List.Dimension(NumIP());
		for (int i = 0; i < NumIP(); i++)
			fStrain_last_List[i].Dimension(NumSD());
	}
	
	
	/* allocate stress list */
	fStress_List.Dimension(NumIP());
	for (int j = 0; j < NumIP(); j++)
		fStress_List[j].Dimension(NumSD());
		
	fElementStress_List.Dimension(NumElements(),NumIP());
	for (int i = 0; i < NumElements(); i++)
		for (int j = 0; j < NumIP(); j++)
			fElementStress_List[i,j].Dimension(NumSD());

	/* allocate "last" stress list */
	/*
	fStress_last_List.Dimension(NumIP());
	for (int j = 0; j < NumIP(); j++)
		fStress_last_List[j].Dimension(NumSD());
		*/
	
	/** \name localization element info storage */
	/*@{*/
	/** fixed values */
	fElementLocNormal.Dimension(NumElements(),NumSD());
	fElementLocNormal = 0.0;
	fElementLocTangent.Dimension(NumElements(),NumSD());
	fElementLocTangent = 0.0;
	fElementLocPsi.Dimension(NumElements());
	fElementLocPsi = 0.0;
	fElementLocGradEnh.Dimension(NumElements(),NumSD()*NumIP());
	fElementLocGradEnh = 0.0;
	fElementLocGradEnhIP.Dimension(NumIP(),NumSD());
	fElementLocGradEnhIP = 0.0;
	/*
	hardcode element edge number for hexahedral element for now;
	should generalize to handle quadrilaterals as well, in the future;
	band tracing NOT currently implemented
	*/
	int numedges = 12;
	fElementLocEdgeIntersect.Dimension(NumElements(),numedges);
	fElementLocEdgeIntersect = 0.0;
	
	fElementLocStartSurface.Dimension(NumElements(),NumSD());
	fElementLocStartSurface = 0.0;
	
	fElementLocISVSoften.Dimension(NumElements(),kNUM_ISV_TERMS);
	fElementLocISVSoften = 0;
	
	fElementLocNodesActive.Dimension(NumElements(),NumElementNodes());
	fElementLocNodesActive = 0;
	
	/** variable from time step to time step */
	fElementLocScalars.Dimension(NumElements(),kNUM_SCALAR_TERMS);	
	fElementLocScalars = 0.0;
	fElementLocScalars_last.Dimension(NumElements(),kNUM_SCALAR_TERMS);
	fElementLocScalars_last = 0.0;
	fElementLocSlipDir.Dimension(NumElements(),NumSD());
	fElementLocSlipDir = 0.0;
	fElementLocSlipDir_last.Dimension(NumElements(),NumSD());
	fElementLocSlipDir_last = 0.0;
	fElementLocMuDir.Dimension(NumElements(),NumSD());
	fElementLocMuDir = 0.0;
	fElementLocMuDir_last.Dimension(NumElements(),NumSD());
	fElementLocMuDir_last = 0.0;
	fElementLocInternalVars.Dimension(NumElements(),kNUM_ISV_TERMS);	
	fElementLocInternalVars = 0.0;
	fElementLocInternalVars_last.Dimension(NumElements(),kNUM_ISV_TERMS);
	fElementLocInternalVars_last = 0.0;
	fElementVolume.Dimension(NumElements());
	fElementVolume = 0.0;
	fElementVolume_last.Dimension(NumElements());
	fElementVolume_last = 0.0;
	fElementLocFlag.Dimension(NumElements());
	fElementLocFlag = 0;
	fElementYieldTrial.Dimension(NumElements());
	fElementYieldTrial = 0.0;
	fElementStress.Dimension(NumElements(),NumIP()*dSymMatrixT::NumValues(NumSD()));
	fElementStress = 0.0;
	fElementStress_last.Dimension(NumElements(),NumIP()*dSymMatrixT::NumValues(NumSD()));
	fElementStress_last = 0.0;
	fElementLocKzetad.Dimension(NumElements(),NumElementNodes()*NumDOF());
	fElementLocKzetad = 0.0;
	fElementLocKdzeta.Dimension(NumElements(),NumElementNodes()*NumDOF());
	fElementLocKdzeta = 0.0;
	/*@}*/
	
	// increment and last iterated displacement
	fLocdeltaDisp.Dimension(fLocDisp.Length());
	fLocdeltaDisp = 0.0;
	fLocLastIterateDisp.Dimension(fLocDisp.Length());
	fLocLastIterateDisp = 0.0;
	fElementLastIterateDisp.Dimension(NumElements(), fLocDisp.Length());
	fElementLastIterateDisp = 0.0;
	
	/* allocate element centroid space */
	fElementCentroid.Dimension(NumElements(),NumSD());
	fElementCentroid = 0.0;
	
	/* localization element arrays */
	normals.Dimension(NumSD());
	slipdirs.Dimension(NumSD());
	
	normal_tmp.Dimension(NumSD());
	slipdir_tmp.Dimension(NumSD());
	tangent_tmp.Dimension(NumSD());
	
	normal_chosen.Dimension(NumSD());
	slipdir_chosen.Dimension(NumSD());
	tangent_chosen.Dimension(NumSD());
	
	mu_dir.Dimension(NumSD());
	mu_dir_last.Dimension(NumSD());
		
	start_surface_coord_read.Dimension(NumSD());
	
	fDe.Dimension(dSymMatrixT::NumValues(NumSD()));
	
	fK_dd.Dimension(NumElementNodes()*NumDOF());
	fK_dd = 0.0;
	fK_dzeta_x_Kzetad.Dimension(NumElementNodes()*NumDOF());
	fK_dzeta_x_Kzetad = 0.0;
	
	fK_dzeta.Dimension(NumElementNodes()*NumDOF());
	fK_dzeta = 0.0;
	fK_zetad.Dimension(NumElementNodes()*NumDOF());
	fK_zetad = 0.0;
	
	/* to NOT let localized element be traced, set fFirstTrace = true */
	//fFirstTrace = false;
	fFirstTrace = true;
	
	// set model type
	model_type = list.GetParameter("model_type");
	
	/* cohesive surface material parameters */
	fCohesiveSurface_Params.Dimension ( kNUM_CS_TERMS );
	fCohesiveSurface_Params[kc_r] = list.GetParameter("residual_cohesion");
	// c_p is determined at localization
	//fCohesiveSurface_Params[kc_p] = list.GetParameter("peak_cohesion");
	fCohesiveSurface_Params[kc_p] = 0.0;
	// store and check element by element, since can vary between multiple elements
	fCohesionParam.Dimension(NumElements());
	fCohesionParam = 0.0;
	fCohesiveSurface_Params[kalpha_c] = list.GetParameter("cohesion_softening_coefficient");
	fCohesiveSurface_Params[kphi_r] = list.GetParameter("residual_friction_angle_rad");
	fCohesiveSurface_Params[kphi_p] = list.GetParameter("peak_friction_angle_rad");
	fCohesiveSurface_Params[kalpha_phi] = list.GetParameter("friction_softening_coefficient");
	// psi_p is determined at localization
	//fCohesiveSurface_Params[kpsi_p] = list.GetParameter("peak_dilation_angle_rad");
	fCohesiveSurface_Params[kpsi_p] = 0.0;
	fCohesiveSurface_Params[kalpha_psi] = list.GetParameter("dilation_softening_coefficient");
	
	choose_normal = list.GetParameter("choose_normal");
	
	choose_element = list.GetParameter("choose_element");
	
	/* need to initialize previous volume */
	Top();
	while (NextElement())
	{
		/* inherited - computes gradients and standard 
		 * deformation gradients */
		SolidElementT::SetGlobalShape();

		/* compute mean of shape function gradients */
		double& vol = fElementVolume_last[CurrElementNumber()];
		SetMeanGradient(fMeanGradient, vol);
		
		// initialize phi as phi_p and psi as psi_p
		int elem = CurrElementNumber();
		fElementLocInternalVars[kNUM_ISV_TERMS*elem + kFriction] = fCohesiveSurface_Params[kphi_p];
		fElementLocInternalVars_last[kNUM_ISV_TERMS*elem + kFriction] = fCohesiveSurface_Params[kphi_p];
		// psi will be determined by LocPsi
		fElementLocInternalVars[kNUM_ISV_TERMS*elem + kDilation] = fCohesiveSurface_Params[kpsi_p];
		fElementLocInternalVars_last[kNUM_ISV_TERMS*elem + kDilation] = fCohesiveSurface_Params[kpsi_p];
	}
	
	
	/* create output file for debugging */
	if (fDeBug)
	{
		outputPrecision = 10;
		outputFileWidth = outputPrecision + 8;
		
		/*
		if (fFirstPass) 
		{
			ss_enh_out.open("ss_enh.info");
			fFirstPass = false;
		}
		else ss_enh_out.open_append("ss_enh.info");
		*/
		
		ss_enh_out.open("ss_enh.info");
		ss_enh_isv.open("ss_enh_isv.txt");
	}
	
}

/* extract the list of material parameters */
void SmallStrainEnhLocT::CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const
{
	const char caller[] = "SmallStrainEnhLocT::CollectMaterialInfo";
	
	/* initialize */
	mat_params.Clear();
	
	/* collected material parameters */
	int num_blocks = all_params.NumLists("small_strain_enh_loc_element_block");
	for (int i = 0; i < num_blocks; i++) {

		/* block information */	
		const ParameterListT& block = all_params.GetList("small_strain_enh_loc_element_block", i);
		
		/* resolve material list name */
		if (i == 0) {
			const ParameterListT& mat_list_params = block.GetListChoice(*this, "small_strain_enh_loc_material_choice");
			mat_params.SetName(mat_list_params.Name());
		}
		
		/* collect material parameters */
		const ParameterListT& mat_list = block.GetList(mat_params.Name());
		const ArrayT<ParameterListT>& mat = mat_list.Lists();
		mat_params.AddList(mat[0]);
	}
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* construct a new material support and return a pointer */
MaterialSupportT* SmallStrainEnhLocT::NewMaterialSupport(MaterialSupportT* p) const
{
	
	// allocate
	if (!p) p = new SSEnhLocMatSupportT(NumDOF(), NumIP());

	// inherited initializations
	SolidElementT::NewMaterialSupport(p);
	
	// set SSEnhLocMatSupportT fields
	SSEnhLocMatSupportT* ps = TB_DYNAMIC_CAST(SSEnhLocMatSupportT*, p);
	if (ps) {
		ps->SetLinearStrain(&fStrain_List);
		ps->SetLinearStrain_last(&fStrain_last_List);
		/*ps->SetElementStress(&fStress_List);
		ps->SetElementLocFlag(&fLocFlag);*/
		ps->SetElementStress(&fElementStress_List);
		ps->SetElementLocFlag(&fElementLocFlag);
	}

	return p;
}


/* return a pointer to a new material list */
MaterialListT* SmallStrainEnhLocT::NewMaterialList(const StringT& name, int size)
{
	/* resolve dimension */
	int nsd = -1;
	if (name == "small_strain_material_1D") nsd = 1;
	else if (name == "small_strain_material_2D") nsd = 2;
	else if (name == "small_strain_material_3D") nsd = 3;
	
	/* no match */
	if (nsd == -1) return NULL;

	/* full list */
	if (size > 0)
	{
		/* material support */
		if (!fSSEnhLocMatSupport) {
			fSSEnhLocMatSupport = TB_DYNAMIC_CAST(SSEnhLocMatSupportT*, NewMaterialSupport());
			if (!fSSEnhLocMatSupport)
				ExceptionT::GeneralFail("SmallStrainEnhLocT::NewMaterialList");
		}
		
		if (nsd == 1)
			return new SSSolidMatList1DT(size, *fSSEnhLocMatSupport);
		else if (nsd == 2)
			return new SSSolidMatList2DT(size, *fSSEnhLocMatSupport);
		else if (nsd == 3)
			return new SSSolidMatList3DT(size, *fSSEnhLocMatSupport);
	}
	else
	{
		if (nsd == 1)
			return new SSSolidMatList1DT;
		else if (nsd == 2)
			return new SSSolidMatList2DT;
		else if (nsd == 3)
			return new SSSolidMatList3DT;
	}
	
	/* no match */
	return NULL;
}


/* check for localization */
void SmallStrainEnhLocT::CheckLocalization(int& elem, LocalArrayT& displ_elem)
{
	dMatrixT grad_displ;
	grad_displ.Dimension(NumSD(),NumSD());
	grad_displ = 0.0;
	double grad_displ_mn_scalar = 0.0;
	
	fShapes->TopIP();
	while ( fShapes->NextIP() )
	{
		/* determine localization */
		//bool checkloc = fCurrMaterial->IsLocalized(normals,slipdirs);
		//bool checkloc = fCurrMaterial->IsLocalized(normals,slipdirs,detAmin);
		bool checkloc = fCurrMaterial->IsLocalized(normals,slipdirs,detAs,dissipations_fact);
		//bool checkloc = fCurrMaterial->IsLocalized_DB(normals,slipdirs,detAs,dissipations_fact);
		
		/* check if localized and if has minimum determinant for IPs */
		//if (checkloc && detAmin < fElementLocScalars[kNUM_SCALAR_TERMS*elem + kdetAmin])
		/* check if localized */
		if (checkloc)
		{
			/* determine minimum determinant among normals at different IPs */
			double detAmin = 1.0e99;
			double detA_tmp;
			detAs.Top();
			while (detAs.Next())
			{
				detA_tmp = detAs.Current();
				if (detA_tmp < detAmin) detAmin = detA_tmp;
			}
			
			/* determine maximum dissipation among normals at different IPs */
			double dissip_max = 0.0;
			double dissip_tmp;
			dissipations_fact.Top();
			while (dissipations_fact.Next())
			{
				dissip_tmp = dissipations_fact.Current();
				if (dissip_tmp > dissip_max) dissip_max = dissip_tmp;
			}
			
			/* if maximum dissipation, do further localization calcs */
			/*
			if ( dissip_max > fElementLocScalars[kNUM_SCALAR_TERMS*elem + kdissip_max] || 
				(dissip_max/fElementLocScalars[kNUM_SCALAR_TERMS*elem + kdissip_max] < 1.1 &&
				dissip_max/fElementLocScalars[kNUM_SCALAR_TERMS*elem + kdissip_max] > 0.0) )
			*/
			/* if minimum determinant, do further localization calcs */
			if ( detAmin < fElementLocScalars[kNUM_SCALAR_TERMS*elem + kdetAmin] || 
				(detAmin/fElementLocScalars[kNUM_SCALAR_TERMS*elem + kdetAmin] < 1.1 &&
				detAmin/fElementLocScalars[kNUM_SCALAR_TERMS*elem + kdetAmin] > 0.0) )
			{
				if (fDeBug)
				{
					const int ip = fShapes->CurrIP()+1;
					ss_enh_out	<< endl 
								<< setw(outputFileWidth) << "element " << elem 
								<< setw(outputFileWidth) <<  "IP" << ip;
				}
				
				loc_flag = 1; // localized, not traced
				fElementLocFlag[elem] = loc_flag;
				
				fElementLocScalars[kNUM_SCALAR_TERMS*elem + kdetAmin] = detAmin;
				
				fElementLocScalars[kNUM_SCALAR_TERMS*elem + kdissip_max] = dissip_max;
				
				/* calculate tangent unit vectors and write output for debugging */
				normals.Top();
				slipdirs.Top();
				tangents.Free();
				psis.Free();
				grad_displ_mns.Free();
				detAs.Top();
				dissipations_fact.Top();
				int num_normals = normals.Length();
				dArrayT dummyt(3);
				while (normals.Next())
				{
					normal_tmp = normals.Current();
					slipdirs.Next();
					slipdir_tmp = slipdirs.Current();
					double product = dArrayT::Dot(normal_tmp,slipdir_tmp);
					double psi_tmp = asin(product);
					psis.Append(psi_tmp);
					if (fabs(psi_tmp-Pi/2.0) < smallnum)
					{
						tangent_tmp = 0.0;
					}
					else 
					{
						double sec = 1.0/cos(psi_tmp);
						tangent_tmp = slipdir_tmp;
						tangent_tmp *= sec;
						dummyt = normal_tmp;
						dummyt *= tan(psi_tmp);
						tangent_tmp -= dummyt;
					}
					tangents.Append(tangent_tmp);
					
					// calculate inner product \bnablabu : \bm \otimes \bn at current IP
					IP_ComputeGradient(displ_elem, grad_displ);
					grad_displ_mn_scalar = grad_displ.MultmBn(slipdir_tmp, normal_tmp);
					grad_displ_mns.Append(grad_displ_mn_scalar);
	
					if (fDeBug)
					{
						detAs.Next();
						detA_tmp = detAs.Current();
						dissipations_fact.Next();
						dissip_tmp = dissipations_fact.Current();
						ss_enh_out	<< endl << "detA_min" << setw(outputFileWidth) << "dissip_max"  << setw(outputFileWidth) << "psi (rad)"
									<< setw(outputFileWidth) << "grad_displ_mn";
						ss_enh_out	<< endl << detA_tmp << setw(outputFileWidth) << dissip_tmp << setw(outputFileWidth) << psi_tmp
									<< setw(outputFileWidth) << grad_displ_mn_scalar; 			
						if (NumSD() == 2)
						{
							ss_enh_out	<< endl << "normal: " << setw(outputFileWidth) << normal_tmp[0] 
									<< setw(outputFileWidth) << normal_tmp[1]
									<< setw(outputFileWidth) << "slipdir: " << setw(outputFileWidth) << slipdir_tmp[0] 
									<< setw(outputFileWidth) << slipdir_tmp[1]; 
							ss_enh_out	<< endl << "tangent:" << setw(outputFileWidth) << tangent_tmp[0] 
									<< setw(outputFileWidth) << tangent_tmp[1];
						}
						else if (NumSD() == 3)
						{
							ss_enh_out	<< endl << "normal: " << setw(outputFileWidth) << normal_tmp[0] 
									<< setw(outputFileWidth) << normal_tmp[1] <<  setw(outputFileWidth) << normal_tmp[2]
									<< setw(outputFileWidth) << "slipdir: " << setw(outputFileWidth) << slipdir_tmp[0] 
									<< setw(outputFileWidth) << slipdir_tmp[1] <<  setw(outputFileWidth) << slipdir_tmp[2]; 
							ss_enh_out	<< endl << "tangent:" << setw(outputFileWidth) << tangent_tmp[0] 
									<< setw(outputFileWidth) << tangent_tmp[1] <<  setw(outputFileWidth) << tangent_tmp[2];
						}
					}
				} // while normals.Next
				
				/* assign current normals, slipdirs, tangents, and detAs 
				 * as values associated with minimum determinant for element */
				normals_min.Free();
				slipdirs_min.Free();
				tangents_min.Free();
				psis_min.Free();
				detAs_min.Free();
				dissipations_fact_min.Free();
				normals_min = normals;
				slipdirs_min = slipdirs;
				tangents_min = tangents;
				psis_min = psis;
				detAs_min = detAs;
				dissipations_fact_min = dissipations_fact;
				grad_displ_mns_min = grad_displ_mns;
			
			} // if (detAmin < fElementLocScalars[kNUM_SCALAR_TERMS*elem + kdetAmin])

		} // if (checkloc)
		
	} // while next IP
	
}


/* choose the normal and slipdir given normals and slipdirs from bifurcation condition */
void SmallStrainEnhLocT::ChooseNormalAndSlipDir(LocalArrayT& displ_elem, int& elem, int& nen)
{
	/* determine normal by using minimum determinant associated with normal */
	/*
	normals_min.Top();
	slipdirs_min.Top();
	tangents_min.Top();
	psis_min.Top();
	detAs_min.Top();
	detAmin = 1.0e99;
	while (normals_min.Next())
	{
		detAs_min.Next();
		detA_tmp = detAs_min.Current();
		slipdirs_min.Next();
		tangents_min.Next();
		psis_min.Next();
		if (detA_tmp < detAmin)
		{
			detAmin = detA_tmp;
			normal_chosen = normals_min.Current();
			slipdir_chosen = slipdirs_min.Current();
			tangent_chosen = tangents_min.Current();
			psi_chosen = psis_min.Current();
		}
	}
	*/

	/* determine normal based on maximum dissipation calculation */
	/*
	normals_min.Top();
	slipdirs_min.Top();
	tangents_min.Top();
	psis_min.Top();
	dissipations_fact_min.Top();
	dissip_max = 0.0;
	while (normals_min.Next())
	{
		dissipations_fact_min.Next();
		dissip_tmp = dissipations_fact_min.Current();
		slipdirs_min.Next();
		tangents_min.Next();
		psis_min.Next();
		if (dissip_tmp > dissip_max)
		{
			dissip_max = dissip_tmp;
			normal_chosen = normals_min.Current();
			slipdir_chosen = slipdirs_min.Current();
			tangent_chosen = tangents_min.Current();
			psi_chosen = psis_min.Current();
		}
	}
	*/

	/* determines normal for shear loading,
	 * which is not applicable for cracking under pure tension, 
	 * or compaction banding under compression */
	/*
	
	double product, sum, sum_min, sum_tmp;
	AutoArrayT <double> sums, sum_products;
	
	dArrayT node_displ(3);
	
	//initialize sums
	sums.Free();
	normals_min.Top();
	while (normals_min.Next())
	{
		sum = 0.0;
		sums.Append(sum);
	}
		
	for (int i=0; i < nen; i++)
	{
		// grab nodal displacement vector
		if (NumSD() == 2)
		{
			node_displ[0] = displ_elem[i];
			node_displ[1] = displ_elem[i+nen];
		}
		else if (NumSD() == 3)
		{
			node_displ[0] = displ_elem[i];
			node_displ[1] = displ_elem[i+nen];
			node_displ[2] = displ_elem[i+2*nen];
		}
		
		// take inner product of nodal displacement vector and normal;
		// then sum up inner products over nodes
		normals_min.Top();
		sum_products.Free();
		sums.Top();
		while (normals_min.Next())
		{
			normal_tmp = normals_min.Current();
			product = dArrayT::Dot(normal_tmp,node_displ);
			sums.Next();
			sum = sums.Current();
			sum += product;
			sum_products.Append(sum);
		}
		// put back sum on sums list
		sums.Free();
		sum_products.Top();
		while (sum_products.Next())
		{
			sum = sum_products.Current();
			sums.Append(sum);
		}
	}
	
	normals_min.Top();
	slipdirs_min.Top();
	tangents_min.Top();
	psis_min.Top();
	sums.Top();
	sum_min = 1.0e99;
	while (normals_min.Next())
	{
		sums.Next();
		sum_tmp = sums.Current();
		// this check works for shear localization, not for pure tension or compression
		slipdirs_min.Next();
		tangents_min.Next();
		psis_min.Next();
		if (sum_tmp < sum_min)
		{
			sum_min = sum_tmp;
			normal_chosen = normals_min.Current();
			slipdir_chosen = slipdirs_min.Current();
			tangent_chosen = tangents_min.Current();
			psi_chosen = psis_min.Current();
		}
	}
	*/
	
	
	/* determine normal based on maximum \bnablabu : \bm \otimes \bn */
	normals_min.Top();
	slipdirs_min.Top();
	tangents_min.Top();
	psis_min.Top();
	grad_displ_mns_min.Top();
	double grad_displ_mn_max = 0.0, grad_displ_mn_tmp;
	while (normals_min.Next())
	{
		grad_displ_mns_min.Next();
		grad_displ_mn_tmp = grad_displ_mns_min.Current();
		slipdirs_min.Next();
		tangents_min.Next();
		psis_min.Next();
		if (grad_displ_mn_tmp > grad_displ_mn_max)
		{
			grad_displ_mn_max = grad_displ_mn_tmp;
			normal_chosen = normals_min.Current();
			slipdir_chosen = slipdirs_min.Current();
			tangent_chosen = tangents_min.Current();
			psi_chosen = psis_min.Current();
		}
	}
	
	
	/* determine normal based on arbitrary choice */
	if ( choose_normal > 0 )
	{
		normals_min.Top();
		slipdirs_min.Top();
		tangents_min.Top();
		psis_min.Top();
		while (normals_min.Next())
		{
			int position = normals_min.Position() + 1;
			slipdirs_min.Next();
			tangents_min.Next();
			psis_min.Next();
			if (position == choose_normal)
			{
				normal_chosen = normals_min.Current();
				slipdir_chosen = slipdirs_min.Current();
				tangent_chosen = tangents_min.Current();
				psi_chosen = psis_min.Current();
			}
		}
	}
	
	/* store chosen normal and slip direction vectors */
	fElementLocNormal.SetRow(elem, normal_chosen);
	fElementLocSlipDir.SetRow(elem, slipdir_chosen);
	fElementLocTangent.SetRow(elem, tangent_chosen);
	fElementLocPsi[elem] = psi_chosen;
	//fElementLocInternalVars[kNUM_ISV_TERMS*elem + kDilation] = psi_chosen;
	//fElementLocInternalVars_last[kNUM_ISV_TERMS*elem + kDilation] = psi_chosen;
	//fCohesiveSurface_Params[kpsi_p] = psi_chosen;
}

/* given the normal and one point, determine active nodes */
void SmallStrainEnhLocT::DetermineActiveNodesTrace(LocalArrayT& coords_elem, int& elem, int& nen)
{
	/* fetch chosen normal */
	fElementLocNormal.RowCopy(elem, normal_chosen);
	
	dArrayT start_surface_coords(NumSD()), elem_centroid(NumSD());
	
	/* calculate element centroid */
	coords_elem.Average(elem_centroid);
	fElementCentroid.SetRow(elem, elem_centroid);
	
	/* fetch slip surface starting intersection point */
	// for one chosen element
	if (!fFirstTrace && elem == choose_element-1)
	{
		fElementLocStartSurface.SetRow(elem, start_surface_coord_read);
		if (start_surface_coord_read.Magnitude() > verysmallnum) 
		{
			start_surface_coords = start_surface_coord_read;
		}
		else 
		{
			start_surface_coords = elem_centroid;
		}
		fFirstTrace = true;
		fElementLocFlag[elem] = 2;
	}
	// for all localized elements
	else if (!fFirstTrace)
	{
		start_surface_coords = elem_centroid;
		fElementLocFlag[elem] = 2;
	}
	else
	{
		//determine whether element is adjacent to traced element to see whether it can be traced
		
		//if found adjacent, set start_surface_coords
		//fElementLocStartSurface.SetRow(elem, start_surface_coords);
		
		//fElementLocFlag[elem] = 2;
		fElementLocFlag[elem] = 0; // temporary
	}
	
	/*
	read first start point from input file;
	calculate others based on slip surface intersection with adjacent elements
	*/
	
	dArrayT diff_vector(NumSD()), node_coords(NumSD());
	double product;
	
	/* loop through nodes and determine active nodes */
	for (int i=0; i < nen; i++)
	{
		if (NumSD() == 2)
		{
			node_coords[0] = coords_elem[i];
			node_coords[1] = coords_elem[i+nen];
		}
		else if (NumSD() == 3)
		{
			node_coords[0] = coords_elem[i];
			node_coords[1] = coords_elem[i+nen];
			node_coords[2] = coords_elem[i+2*nen];
		}
		diff_vector.DiffOf(node_coords,start_surface_coords);
		product = dArrayT::Dot(normal_chosen,diff_vector);
		if (product > 0.0)
		{
			fElementLocNodesActive[NumElementNodes()*elem + i] = 1;
			//calculate enhancement in FormKd
		}
	}
	
	/* loop through nodes and determine edge intersectons */
	/*
	for (int i = 0; i < numedges; i++)
	{
		//grab associated nodes with edges??
		fElementLocEdgeIntersect[elem,i] = ?;
	}
	*/

}

/* calculate the internal force contribution ("-k*d") 
   plus other terms due to jump displacement */
void SmallStrainEnhLocT::FormKd(double constK)
{
	/* collect incremental heat */
	bool need_heat = fElementHeat.Length() == fShapes->NumIP();

	/* current element number */
	int elem = CurrElementNumber();
	loc_flag = fElementLocFlag[elem];
	fLocFlag = fElementLocFlag[elem];
	int nen = NumElementNodes();
	double vol = fElementVolume[elem];
	double divide_vol = 0.0;
	if ( fabs(vol) > verysmallnum ) divide_vol = 1.0/vol;
	//int iter = fSSEnhLocMatSupport->IterationNumber();
	//int iter = fSSMatSupport->IterationNumber();
	int iter = IterationNumber();
	
	dArray2DT stress_IPs(NumIP(),dSymMatrixT::NumValues(NumSD()));
	dArray2DT stress_last_IPs(NumIP(),dSymMatrixT::NumValues(NumSD()));
	
	dArrayT tmp_array(NumSD()), tmp_Dslipdir_array(NumSD());
	
	dArrayT DqDzeta(kNUM_ISV_TERMS), Dh_qDzeta(kNUM_ISV_TERMS), Dh_qDpsi(kNUM_ISV_TERMS), DfDq(kNUM_ISV_TERMS);
	dMatrixT Dh_qDq(kNUM_ISV_TERMS), A_matrix(kNUM_ISV_TERMS), tmp_A_matrix(kNUM_ISV_TERMS);
	
	dArrayT tmp_h_q(kNUM_ISV_TERMS), tmp_q_zeta(kNUM_ISV_TERMS);
	dArrayT q_isv(kNUM_ISV_TERMS), q_isv_last(kNUM_ISV_TERMS), h_q(kNUM_ISV_TERMS);
	
	dArrayT DslipdirDpsi(NumSD()), DslipdirDzeta(NumSD());
	
	dMatrixT DslipdirDq(NumSD(), kNUM_ISV_TERMS);
	
	dArrayT grad_enh_IP(NumSD());
		
	dSymMatrixT tmp_sym_matrix(NumSD()), inner_matrix(NumSD());
	dSymMatrixT DF_munDzeta(NumSD()), DGenhDzeta(NumSD());
	
	dSymMatrixT strain_inc(NumSD()), stress_inc(NumSD()), stress_last(NumSD()), stress_return(NumSD());

	double inner, inner1, inner2, inner_abs;

	dSymMatrixT fStressTrial(NumSD()), fStressCurr(NumSD()), DsigDzeta(NumSD());
	
	dArrayT tmp_vec1(NumElementNodes()*NumDOF());
	
	dSymMatrixT F_nn(NumSD()), F_tn(NumSD()), F_mun(NumSD());
	dSymMatrixT F_mun_last(NumSD()), G_enh(NumSD());
	dMatrixT F_mun_nonsym(NumSD()), F_nn_nonsym(NumSD());
	dMatrixT tmp_matrix(NumSD());
	
	double Q_S_last, Q_Sn_trial, Q_S;
	double P_S_last, P_S_trial, P_S;
	double q_St_last, q_St_trial, q_St, sign_q_St;
	double q_St_abs_last, q_St_abs_trial, q_St_abs;
	double DQ_SDzeta, DQ_SnDzeta, DP_SDzeta, Dr_SDzeta;
	
	double zeta, zeta_last, delta_zeta, Delta_zeta, gamma_delta, sign_Delta_zeta;
	double K_zetazeta;
	
	double tanphi_n;
	
	double cospsi, sinpsi, tanphi, secphi2, cospsi2;
	
	double h_c, h_phi, h_psi;
	double Dh_cDzeta, Dh_phiDzeta, Dh_psiDzeta;
	double Dh_cDpsi, Dh_phiDpsi, Dh_psiDpsi;
	
	double DpsiDzeta, DphiDzeta;
	double DcospsiDzeta, DgammadeltaDzeta;
	
	double c_r, c_p, alpha_c, phi_r, phi_p, alpha_phi, psi_p, alpha_psi, psi;
	
	if ( loc_flag == 2 )
	{
		/* fetch material constants */
		c_r = fCohesiveSurface_Params[kc_r];
		//c_p = fCohesiveSurface_Params[kc_p];
		//check per element rather than for each element, because c_p will be calculated based on element average shear stress
		c_p = fCohesionParam[elem];
		alpha_c = fCohesiveSurface_Params[kalpha_c];
		phi_r = fCohesiveSurface_Params[kphi_r];
		phi_p = fCohesiveSurface_Params[kphi_p];
		alpha_phi = fCohesiveSurface_Params[kalpha_phi];
		//psi_p = fCohesiveSurface_Params[kpsi_p];
		psi_p = fElementLocPsi[elem];
		alpha_psi = fCohesiveSurface_Params[kalpha_psi];
	
		/* fetch normal and MuDir_last for element */
		fElementLocNormal.RowCopy(elem, normal_chosen);
		fElementLocTangent.RowCopy(elem, tangent_chosen);
		fElementLocMuDir_last.RowCopy(elem, mu_dir_last);
		
		tmp_matrix.Outer(mu_dir_last, normal_chosen);
		F_mun_last.Symmetrize(tmp_matrix);
		
		F_nn_nonsym.Outer(normal_chosen, normal_chosen);
		F_nn.Symmetrize(F_nn_nonsym);
		
		tmp_matrix.Outer(tangent_chosen, normal_chosen);
		F_tn.Symmetrize(tmp_matrix);
		
		// calc tanphi from phi_last
		double phi_tmp = fElementLocInternalVars_last[kNUM_ISV_TERMS*elem + kFriction];		
		tanphi_n = tan(phi_tmp);
		
		// initialize resolved stresses and derivatives
		Q_S_last = 0.0;
		Q_Sn_trial = 0.0;
		Q_S = 0.0;
		P_S_last = 0.0;
		P_S_trial = 0.0;
		P_S = 0.0;
		q_St_last = 0.0;
		q_St_trial = 0.0;
		q_St = 0.0;
		q_St_abs_last = 0.0;
		q_St_abs_trial = 0.0;
		q_St_abs = 0.0;
		DQ_SDzeta = 0.0;
		DQ_SnDzeta = 0.0;
		DP_SDzeta = 0.0;
		
		fK_dzeta = 0.0;
		
		DqDzeta = 0.0;
		Dh_qDzeta = 0.0;
		Dh_qDq = 0.0;
		Dh_qDpsi = 0.0;
		DfDq = 0.0;
		A_matrix = 0.0;
		
		DslipdirDpsi = 0.0;
		DslipdirDzeta = 0.0;
		DslipdirDq = 0.0;
		
		/* get displacements and displacement increments, 
		   and update delta_zeta;
		   delta d = d^{k+1} - d^k
		   */
	    SetLocalU(fLocDisp);
    	// remap nodal displacements
	    dArrayT nodal_displ;
	    nodal_displ.Dimension(fLocDisp.Length());
	    nodal_displ = 0.0;
	    for (int i=0; i < nen; i++)
		{
			int ii = NumSD()*i;
			if (NumSD() == 2)
			{
				nodal_displ[ii] = fLocDisp[i];
				nodal_displ[ii+1] = fLocDisp[i+nen];
			}
			else if (NumSD() == 3)
			{
				nodal_displ[ii] = fLocDisp[i];
				nodal_displ[ii+1] = fLocDisp[i+nen];
				nodal_displ[ii+2] = fLocDisp[i+2*nen];
			}
		}
		
		// calculate displacement increment
    	fElementLastIterateDisp.RowCopy(elem, fLocLastIterateDisp);
    	fLocdeltaDisp = nodal_displ;
    	fLocdeltaDisp -= fLocLastIterateDisp;
    	fElementLastIterateDisp.SetRow(elem, nodal_displ);
    	
    	// calculate delta_zeta based on previous residual and matrices
    	fElementLocKzetad.RowCopy(elem, fK_zetad);
		delta_zeta = dArrayT::Dot(fK_zetad,fLocdeltaDisp);
		delta_zeta += fElementLocScalars[kNUM_SCALAR_TERMS*elem + kr_S];
		delta_zeta *= -1.0;
		K_zetazeta = fElementLocScalars[kNUM_SCALAR_TERMS*elem + kKzetazeta];
		if ( fabs(K_zetazeta) > verysmallnum ) 
		{
			delta_zeta /= K_zetazeta;
		}
		else
		{
			delta_zeta = 0.0;
		}
		//delta_zeta = 0.5*(delta_zeta+fabs(delta_zeta));
	    
		// update zeta and calculate Delta_zeta
		// zeta^{k+1} = zeta^k + delta zeta
		double zeta_tmp = fElementLocScalars[kNUM_SCALAR_TERMS*elem + kzeta];
		zeta_tmp += delta_zeta;
		// Delta zeta = zeta^{k+1} - zeta_n
		Delta_zeta = zeta_tmp;
		zeta_last = fElementLocScalars_last[kNUM_SCALAR_TERMS*elem + kzeta];
		Delta_zeta -= zeta_last;
		sign_Delta_zeta = 1.0;
		
		if (fabs(Delta_zeta) > verysmallnum) sign_Delta_zeta = Delta_zeta/fabs(Delta_zeta);
		sign_Delta_zeta = 0.5*(sign_Delta_zeta+fabs(sign_Delta_zeta));
		Delta_zeta = 0.5*(Delta_zeta+fabs(Delta_zeta));
		
		//if (zeta_last > verysmallnum) Delta_zeta = 0.5*(Delta_zeta+fabs(Delta_zeta));
		/*
		if (zeta_last > verysmallnum)
		{
			if (fabs(Delta_zeta) > verysmallnum) sign_Delta_zeta = Delta_zeta/fabs(Delta_zeta);
			sign_Delta_zeta = 0.5*(sign_Delta_zeta+fabs(sign_Delta_zeta));
			Delta_zeta = 0.5*(Delta_zeta+fabs(Delta_zeta));
		}
		*/
		zeta = zeta_last + Delta_zeta;
		
		// store zeta and Delta_zeta
		fElementLocScalars[kNUM_SCALAR_TERMS*elem + kzeta] = zeta;
		fElementLocScalars[kNUM_SCALAR_TERMS*elem + kDelta_zeta] = Delta_zeta;
		
		// update internal variables using Delta_zeta
		if ( model_type == 0 || model_type == 1 ) // M-C model
		{
			psi = fElementLocInternalVars[kNUM_ISV_TERMS*elem + kDilation];
			cospsi = cos(psi);
			sinpsi = sin(psi);
			gamma_delta = fElementLocScalars_last[kNUM_SCALAR_TERMS*elem + kgamma_delta] + cospsi*Delta_zeta;
			fElementLocScalars[kNUM_SCALAR_TERMS*elem + kgamma_delta] = gamma_delta;
		}
		
		// calculate ISV hardening functions
		if ( model_type == 0 ) // linear softening
		{
			h_c = -alpha_c*cospsi;
			h_phi = -alpha_phi*cospsi;
			h_psi = -alpha_psi*cospsi;
		}
		else if  ( model_type == 1 ) // exponential softening
		{
			h_c = -alpha_c*(c_p-c_r)*exp(-alpha_c*gamma_delta)*cospsi;
			h_phi = -alpha_phi*(phi_p-phi_r)*exp(-alpha_phi*gamma_delta)*cospsi;
			h_psi = -alpha_psi*psi_p*exp(-alpha_psi*gamma_delta)*cospsi;
		}
		else // default to linear
		{
			h_c = -alpha_c*cospsi;
			h_phi = -alpha_phi*cospsi;
			h_psi = -alpha_psi*cospsi;
		}
		fElementLocScalars[kNUM_SCALAR_TERMS*elem + kh_phi] = h_phi;
				
		h_q[0] = h_c;
		h_q[1] = h_phi;
		h_q[2] = h_psi;
		tmp_h_q = h_q;
		tmp_h_q *= Delta_zeta;
		
		// update ISVs
		q_isv_last = 0.0;
		if (fElementLocISVSoften[kNUM_ISV_TERMS*elem + kCohesion] < 1)
			q_isv_last[0] = fElementLocInternalVars_last[kNUM_ISV_TERMS*elem + kCohesion];
		if (fElementLocISVSoften[kNUM_ISV_TERMS*elem + kFriction] < 1)
			q_isv_last[1] = fElementLocInternalVars_last[kNUM_ISV_TERMS*elem + kFriction];
		if (fElementLocISVSoften[kNUM_ISV_TERMS*elem + kDilation] < 1)
			q_isv_last[2] = fElementLocInternalVars_last[kNUM_ISV_TERMS*elem + kDilation];
		q_isv = q_isv_last;
		q_isv += tmp_h_q;
		
		// check that do not soften less than residual values, for linear softening
		if ( model_type == 0 )
		{
			if ( q_isv[0] < c_r && c_p > verysmallnum )
			{
				q_isv[0] = c_r;
				alpha_c = 0.0;
				fCohesiveSurface_Params[kalpha_c] = alpha_c;
				fElementLocISVSoften[kNUM_ISV_TERMS*elem + kCohesion] = 1;
			}
			if ( q_isv[1] < phi_r )
			{
				q_isv[1] = phi_r;
				alpha_phi = 0.0;
				fCohesiveSurface_Params[kalpha_phi] = alpha_phi;
				fElementLocISVSoften[kNUM_ISV_TERMS*elem + kFriction] = 1;
			}
			if ( q_isv[2] < verysmallnum )
			{
				q_isv[2] = verysmallnum;
				alpha_psi = 0.0;
				fCohesiveSurface_Params[kalpha_psi] = alpha_psi;
				fElementLocISVSoften[kNUM_ISV_TERMS*elem + kDilation] = 1;
			}
			
			// update cospsi
			psi = q_isv[2];
			cospsi = cos(psi);
			sinpsi = sin(psi);
			
			// recalculate hardening functions
			h_c = -alpha_c*cospsi;
			h_phi = -alpha_phi*cospsi;
			fElementLocScalars[kNUM_SCALAR_TERMS*elem + kh_phi] = h_phi;
			h_psi = -alpha_psi*cospsi;
			h_q[0] = h_c;
			h_q[1] = h_phi;
			h_q[2] = h_psi;
		}
		
		// update cospsi
		psi = q_isv[2];
		cospsi = cos(psi);
		sinpsi = sin(psi);
		
		// store updated internal variables
		fElementLocInternalVars[kNUM_ISV_TERMS*elem + kCohesion] = q_isv[0];
		fElementLocInternalVars[kNUM_ISV_TERMS*elem + kFriction] = q_isv[1];
		fElementLocInternalVars[kNUM_ISV_TERMS*elem + kDilation] = q_isv[2];
		
		// update tanphi
		tanphi = tan(q_isv[1]);
		double secphi = 1.0;
		if (fabs(cos(q_isv[1])) > verysmallnum) secphi /= cos(q_isv[1]);
		secphi2 = secphi*secphi;
		fElementLocScalars[kNUM_SCALAR_TERMS*elem + ksecphi2] = secphi2;
		
		// calculate the trial direction sign_q_St
		const double* Det    = fShapes->IPDets();
		const double* Weight = fShapes->IPWeights();
		fShapes->TopIP();
		while (fShapes->NextIP())
		{
			double scale_no_constK = (*Det++)*(*Weight++);
			double scale_vol = scale_no_constK*divide_vol;
			
			const int ip = fShapes->CurrIP();
		
			// fetch last stress at current IP
			stress_last_IPs.RowCopy(ip, stress_last);
			
			// calculate trial stress
			strain_inc=0.0;
			strain_inc.DiffOf(LinearStrain(), LinearStrain_last());
			strain_inc.ScaleOffDiagonal(2.0);
			fStressTrial = stress_last;
			fDe = fCurrMaterial->ce_ijkl();
			fDe.Multx(strain_inc, stress_inc);
			fStressTrial += stress_inc;
			
			// calc q_St_trial
			inner_matrix = 0.0;
			inner_matrix.MultAB(fStressTrial, F_tn);
			inner = inner_matrix.Trace();
			inner *= scale_vol;
			q_St_trial += inner;			
		}
		sign_q_St = 1.0;
		//if (fabs(q_St_trial) > verysmallnum) sign_q_St = q_St_trial/fabs(q_St_trial);
		// change dir if large relative to fraction of elastic modulus??
		if (fabs(q_St_trial) > smallnum) sign_q_St = q_St_trial/fabs(q_St_trial);
		fElementLocScalars[kNUM_SCALAR_TERMS*elem + ksign_q_St] = sign_q_St;
		//reinitialize
		q_St_trial = 0.0;
		
		// update mu_dir using last direction
		if ( model_type == 0 || model_type == 1 ) // M-C model
		{
			//sign_q_St = fElementLocScalars[kNUM_SCALAR_TERMS*elem + ksign_q_St];
			//if ( fabs(sign_q_St) < verysmallnum ) sign_q_St = 1.0;
			tmp_array = normal_chosen;
			tmp_array *= tanphi;
			mu_dir = tangent_chosen;
			mu_dir *= sign_q_St;
			mu_dir += tmp_array;
			// update F_mun
			F_mun_nonsym.Outer(mu_dir, normal_chosen);
			F_mun.Symmetrize(F_mun_nonsym);
			
			// update slipdir \bm
			tmp_array = normal_chosen;
			tmp_array *= sinpsi;
			slipdir_tmp = tangent_chosen;
			slipdir_tmp *= cospsi;
			// be careful that stays in direction of loading
			//slipdir_tmp *= sign_q_St;
			slipdir_tmp += tmp_array;
			fElementLocSlipDir.SetRow(elem, slipdir_tmp);
			slipdir_chosen = slipdir_tmp;
		}
		
		// calculate derivatives of h_q w.r.t. ISVs
		if ( model_type == 0 ) // linear softening
		{
			Dh_cDpsi = alpha_c*sinpsi;
			Dh_phiDpsi = alpha_phi*sinpsi;
			Dh_psiDpsi = alpha_psi*sinpsi;
			Dh_cDzeta = 0.0;
			Dh_phiDzeta = 0.0;
			Dh_psiDzeta = 0.0;
		}
		else if  ( model_type == 1 ) // exponential softening
		{
			Dh_cDpsi = alpha_c*(c_p-c_r)*exp(-alpha_c*gamma_delta)*sinpsi*(1.0-alpha_c*cospsi*Delta_zeta);
			Dh_phiDpsi = alpha_phi*(phi_p-phi_r)*exp(-alpha_phi*gamma_delta)*sinpsi*(1.0-alpha_phi*cospsi*Delta_zeta);
			Dh_psiDpsi = alpha_psi*psi_p*exp(-alpha_psi*gamma_delta)*sinpsi*(1.0-alpha_psi*cospsi*Delta_zeta);
			Dh_cDzeta = alpha_c*alpha_c*cospsi2*(c_p-c_r)*exp(-alpha_c*gamma_delta)*sign_Delta_zeta;
			Dh_phiDzeta = alpha_phi*alpha_phi*cospsi2*(phi_p-phi_r)*exp(-alpha_phi*gamma_delta)*sign_Delta_zeta;
			Dh_psiDzeta = alpha_psi*alpha_psi*cospsi2*psi_p*exp(-alpha_psi*gamma_delta)*sign_Delta_zeta;
		}
		else // default to linear
		{
			Dh_cDpsi = alpha_c*sinpsi;
			Dh_phiDpsi = alpha_phi*sinpsi;
			Dh_psiDpsi = alpha_psi*sinpsi;
			Dh_cDzeta = 0.0;
			Dh_phiDzeta = 0.0;
			Dh_psiDzeta = 0.0;
		}
		
		if ( model_type == 0 || model_type == 1 ) // M-C model
		{
			Dh_qDpsi[0] = Dh_cDpsi;
			Dh_qDpsi[1] = Dh_phiDpsi;
			Dh_qDpsi[2] = Dh_psiDpsi;
			Dh_qDq.SetCol(kDilation, Dh_qDpsi);
			Dh_qDq *= Delta_zeta;
			
			// construct A_matrix
			A_matrix = 0.0;
			A_matrix.PlusIdentity(1.0);
			A_matrix -= Dh_qDq;
			A_matrix.Inverse();
			
			// store derivatives of ISVs w.r.t zeta (jump displ.)
			Dh_qDzeta[0] = Dh_cDzeta;
			Dh_qDzeta[1] = Dh_phiDzeta;
			Dh_qDzeta[2] = Dh_psiDzeta;
			
			// form DqDzeta
			tmp_h_q = h_q;
			tmp_h_q *= sign_Delta_zeta;
			tmp_q_zeta = Dh_qDzeta;
			tmp_q_zeta *= Delta_zeta;
			tmp_q_zeta += tmp_h_q;
			A_matrix.Multx(tmp_q_zeta, DqDzeta);
			
			// form DslipdirDq
			tmp_array = normal_chosen;
			tmp_array *= cospsi;
			DslipdirDpsi = tangent_chosen;
			DslipdirDpsi *= -sinpsi;
			DslipdirDpsi *= sign_q_St;
			DslipdirDpsi += tmp_array;
			DslipdirDq.SetCol(kDilation, DslipdirDpsi);
			
			// calculate DslipdirDzeta
			DslipdirDq.Multx(DqDzeta, DslipdirDzeta);
			
			// calc for use in forming DsigDzeta
			tmp_Dslipdir_array = DslipdirDzeta;
			tmp_Dslipdir_array *= Delta_zeta;
			slipdir_tmp *= sign_Delta_zeta;
			tmp_Dslipdir_array += slipdir_tmp;
		}
		
		// fetch last stress
		fElementStress_last.RowCopy(elem, stress_last_IPs);
	}
	
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	fShapes->TopIP();
	while (fShapes->NextIP())
	{
		double scale_no_constK = (*Det++)*(*Weight++);
		double scale_vol = scale_no_constK*divide_vol;
		double scale = constK*scale_no_constK;
		
		const int ip = fShapes->CurrIP();
		
		/* strain displacement matrix */
		if (fStrainDispOpt == kMeanDilBbar)
			Set_B_bar(fShapes->Derivatives_U(), fMeanGradient, fB);
		else
			Set_B(fShapes->Derivatives_U(), fB);
		
		if ( loc_flag == 2 )
		{
			// fetch last stress at current IP
			stress_last_IPs.RowCopy(ip, stress_last);
			
			// calc q_St_last
			inner_matrix = 0.0;
			inner_matrix.MultAB(stress_last, F_tn);
			inner = inner_matrix.Trace();
			inner_abs = fabs(inner);
			inner *= scale_vol;
			inner_abs *= scale_vol;
			q_St_last += inner;
			q_St_abs_last += inner_abs;
			
			// calc P_S_last
			inner_matrix = 0.0;
			inner_matrix.MultAB(stress_last, F_nn);
			inner = inner_matrix.Trace();
			inner *= scale_vol;
			P_S_last += inner;
			
			// calculate trial stress
			strain_inc=0.0;
			strain_inc.DiffOf(LinearStrain(), LinearStrain_last());
			strain_inc.ScaleOffDiagonal(2.0);
			fStressTrial = stress_last;
			fDe = fCurrMaterial->ce_ijkl();
			fDe.Multx(strain_inc, stress_inc);
			fStressTrial += stress_inc;
			
			// calc q_St_trial
			inner_matrix = 0.0;
			inner_matrix.MultAB(fStressTrial, F_tn);
			inner = inner_matrix.Trace();
			inner_abs = fabs(inner);
			inner *= scale_vol;
			inner_abs *= scale_vol;
			q_St_trial += inner;
			q_St_abs_trial += inner_abs;
			
			// calc P_S_trial
			inner_matrix = 0.0;
			inner_matrix.MultAB(fStressTrial, F_nn);
			inner = inner_matrix.Trace();
			inner *= scale_vol;
			P_S_trial += inner;
			
			/* loop through nodes and calculate enhancement function */
			const dArray2DT& DNa = fShapes->Derivatives_U();			
			dArrayT node_shape_deriv(NumSD());
			grad_enh_IP = 0.0;
			for (int i=0; i < nen; i++)
			{
				DNa.ColumnCopy(i,node_shape_deriv);
				if ( fElementLocNodesActive[NumElementNodes()*elem + i] == 1 )
				{
					grad_enh_IP += node_shape_deriv;
				}
			}
			fElementLocGradEnhIP.SetRow(ip, grad_enh_IP);
			
			// update stress
			stress_return = 0.0;
			tmp_matrix.Outer(slipdir_chosen, grad_enh_IP);
			G_enh.Symmetrize(tmp_matrix);
			G_enh.ScaleOffDiagonal(2.0);
			fDe.Multx(G_enh, stress_return);
			stress_return *= Delta_zeta;
			fStressCurr = fStressTrial;
			fStressCurr -= stress_return;
			stress_IPs.SetRow(ip, fStressCurr);
			fStress_List[ip] = fStressCurr;
			fElementStress_List[elem,ip] = fStressCurr;
	
			// calc q_St
			inner_matrix = 0.0;
			inner_matrix.MultAB(fStressCurr, F_tn);
			inner = inner_matrix.Trace();
			inner_abs = fabs(inner);
			inner *= scale_vol;
			inner_abs *= scale_vol;
			q_St += inner;
			q_St_abs += inner_abs;
			
			// calc P_S
			inner_matrix = 0.0;
			inner_matrix.MultAB(fStressCurr, F_nn);
			inner = inner_matrix.Trace();
			inner *= scale_vol;
			P_S += inner;
			
			// calc DsigDzeta
			tmp_matrix.Outer(tmp_Dslipdir_array, grad_enh_IP);
			tmp_sym_matrix.Symmetrize(tmp_matrix);
			tmp_sym_matrix.ScaleOffDiagonal(2.0);
			tmp_sym_matrix *= -1.0;
			fDe.Multx(tmp_sym_matrix, DsigDzeta);
			
			// calc DQ_SDzeta
			inner_matrix = 0.0;
			inner_matrix.MultAB(DsigDzeta, F_mun);
			inner = inner_matrix.Trace();
			inner *= scale_vol;
			DQ_SDzeta += inner;
			
			/* B^T * Cauchy stress */
			fB.MultTx(fStressCurr, fNEEvec);
			
			/* B^T * derivative Cauchy stress w.r.t zeta */
			fB.MultTx(DsigDzeta, tmp_vec1);

			/* accumulate */
			fRHS.AddScaled(scale, fNEEvec);
			fK_dzeta.AddScaled(scale_no_constK, tmp_vec1);
		}
		else 
		{
			stress_IPs.SetRow(ip, fCurrMaterial->s_ij());
			
			/* B^T * Cauchy stress */
			fB.MultTx(fCurrMaterial->s_ij(), fNEEvec);
			
			/* accumulate */
			fRHS.AddScaled(scale, fNEEvec);
		}
		
		/* incremental heat generation */
		if (need_heat) 
			fElementHeat[fShapes->CurrIP()] += fCurrMaterial->IncrementalHeat();
	}
	
	// store current stress
	fElementStress.SetRow(elem, stress_IPs);
	
	if ( loc_flag == 2 )
	{	
		if ( model_type == 0 || model_type == 1 ) // M-C model
		{
			// ensure that P_S is negative (in compression) or zero along \cS
			double P_S_abs;
			P_S_abs = fabs(P_S_last);
			P_S_last -= P_S_abs;
			P_S_last *= 0.5;
			P_S_abs = fabs(P_S_trial);
			P_S_trial -= P_S_abs;
			P_S_trial *= 0.5;
			P_S_abs = fabs(P_S);
			P_S -= P_S_abs;
			P_S *= 0.5;
		}
			
		// store volume averaged resolved stresses
		fElementLocScalars[kNUM_SCALAR_TERMS*elem + kP_S] = P_S;
		fElementLocScalars[kNUM_SCALAR_TERMS*elem + kq_St] = q_St;
		
		// store enhancement function
		fElementLocGradEnh.SetRow(elem, fElementLocGradEnhIP);
		
		// initialize fYieldTrial and r_S
		double fYieldTrial = 0.0;
		double r_S = 0.0;
		
		if ( model_type == 0 || model_type == 1 ) // M-C model
		{
			// calc sign_q_St and update MuDir
			//sign_q_St = 1.0;
			//if (fabs(q_St) > verysmallnum) sign_q_St = q_St/fabs(q_St);
			//fElementLocScalars[kNUM_SCALAR_TERMS*elem + ksign_q_St] = sign_q_St;
			sign_q_St = fElementLocScalars[kNUM_SCALAR_TERMS*elem + ksign_q_St];
			tmp_array = normal_chosen;
			tmp_array *= tanphi;
			mu_dir = tangent_chosen;
			mu_dir *= sign_q_St;
			mu_dir += tmp_array;
			fElementLocMuDir.SetRow(elem, mu_dir);
		
			// use abs of integral?
			q_St = fabs(q_St);
			q_St_trial = fabs(q_St_trial);
			q_St_last = fabs(q_St_last);
			// or
			// q_St = q_St_abs;
			// q_St_trial = q_St_abs_trial;
			// q_St_last = q_St_abs_last;
			
			// calc Q_S
			Q_S = q_St + tanphi*P_S;
			// store volume averaged resolved stress
			fElementLocScalars[kNUM_SCALAR_TERMS*elem + kQ_S] = Q_S;
			
			// calc Q_Sn_trial
			Q_Sn_trial = q_St_trial + tanphi_n*P_S_trial;
			
			// calc Q_S_last
			Q_S_last = q_St_last + tanphi*P_S_last;
			
			// initialize cohesive strength parameter just after localization
			if ( fabs(zeta) < verysmallnum && c_p < verysmallnum ) 
			{
				c_p = Q_S_last;
				fElementLocInternalVars[kNUM_ISV_TERMS*elem + kCohesion] = c_p;
				fElementLocInternalVars_last[kNUM_ISV_TERMS*elem + kCohesion] = c_p;
				if (c_r > c_p) 
				{
					c_r = 0.1*c_p;
					fCohesiveSurface_Params[kc_r] = c_r;
				}
				//fCohesiveSurface_Params[kc_p] = c_p;
				fCohesionParam[elem] = c_p;
				
				fElementLocInternalVars[kNUM_ISV_TERMS*elem + kDilation] = fElementLocPsi[elem];
				fElementLocInternalVars_last[kNUM_ISV_TERMS*elem + kDilation] = fElementLocPsi[elem];
				
				psi = fElementLocInternalVars[kNUM_ISV_TERMS*elem + kDilation];
				cospsi = cos(psi);
				sinpsi = sin(psi);
				
				if (model_type == 1)
				{
					// recalc h_c
					h_c = -alpha_c*(c_p-c_r)*exp(-alpha_c*gamma_delta)*cospsi;
					h_q[0] = h_c;
					
					// recalculate derivative of h_c w.r.t. ISVs
					Dh_cDpsi = alpha_c*(c_p-c_r)*exp(-alpha_c*gamma_delta)*sinpsi*(1.0-alpha_c*cospsi*Delta_zeta);
					Dh_qDpsi[0] = Dh_cDpsi;
					Dh_qDq = 0.0;
					Dh_qDq.SetCol(kDilation, Dh_qDpsi);
					Dh_qDq *= Delta_zeta;
					// construct A_matrix
					A_matrix = 0.0;
					A_matrix.PlusIdentity(1.0);
					A_matrix -= Dh_qDq;
					A_matrix.Inverse();
					
					// recalculate derivative of ISV w.r.t zeta (jump displ.)
					Dh_cDzeta = alpha_c*alpha_c*cospsi2*(c_p-c_r)*exp(-alpha_c*gamma_delta)*sign_Delta_zeta;
					Dh_qDzeta[0] = Dh_cDzeta;
					
					// reform DqDzeta
					tmp_h_q = h_q;
					tmp_h_q *= sign_Delta_zeta;
					tmp_q_zeta = Dh_qDzeta;
					tmp_q_zeta *= Delta_zeta;
					tmp_q_zeta += tmp_h_q;
					A_matrix.Multx(tmp_q_zeta, DqDzeta);
				}
			}
			
			// calc residual
			r_S = Q_S - fElementLocInternalVars[kNUM_ISV_TERMS*elem + kCohesion];
			fElementLocScalars[kNUM_SCALAR_TERMS*elem + kr_S] = r_S;
			
			// calc DfDq
			DfDq[0] = -1.0;
			DfDq[1] = P_S*secphi2;
			DfDq[2] = 0.0;
			
			// calc Dr_SDzeta
			double dotprod = dArrayT::Dot(DfDq,DqDzeta);
			double Dr_SDzeta = DQ_SDzeta + dotprod;
			K_zetazeta = Dr_SDzeta;
			
			// calculate yield on discontinuity surface
			fYieldTrial = Q_Sn_trial - fElementLocInternalVars_last[kNUM_ISV_TERMS*elem + kCohesion];
			//fYieldTrial = 1.0;
			fElementYieldTrial[elem] = fYieldTrial;
		}
		
		// modify fRHS if yielding
		if (fYieldTrial > verysmallnum)
		{
			// store fK_dzeta
			fElementLocKdzeta.SetRow(elem, fK_dzeta);
			
			// store K_zetazeta
			fElementLocScalars[kNUM_SCALAR_TERMS*elem + kKzetazeta] = K_zetazeta;

			if ( fabs(K_zetazeta) > verysmallnum )
			{
				fK_dzeta /= K_zetazeta;
			}
			else
			{
				fK_dzeta = 0.0;
			}
			fK_dzeta *= r_S;
			fRHS += fK_dzeta;
		}
		else
		{
			// store fK_dzeta
			fK_dzeta = 0.0;
			fElementLocKdzeta.SetRow(elem, fK_dzeta);
			
			// store K_zetazeta
			fElementLocScalars[kNUM_SCALAR_TERMS*elem + kKzetazeta] = 0.0;
		}
		
		if (fDeBug)
		{			
			ss_enh_out	<< endl << endl << "----------------------------------------------------------------------------------------------" << endl;
			ss_enh_out	<< endl 
						<< setw(outputFileWidth) << "element " << elem 
						<< setw(outputFileWidth) << "iteration " << iter;
						
			ss_enh_out	<< endl << endl << "element volume:" << setw(outputFileWidth) << vol; 
						
			ss_enh_out	<< endl << endl << "detA_min" << setw(outputFileWidth) << "dissip_max"  << setw(outputFileWidth) << "psi (rad)";
			ss_enh_out	<< endl << fElementLocScalars[kNUM_SCALAR_TERMS*elem + kdetAmin] << setw(outputFileWidth) << fElementLocScalars[kNUM_SCALAR_TERMS*elem + kdissip_max]
						<< setw(outputFileWidth) << fElementLocPsi[elem]; 
			
			if (NumSD() == 2)	
			{
				ss_enh_out	<< endl << " normal_chosen: " << setw(outputFileWidth) << normal_chosen[0] 
						<< setw(outputFileWidth) << normal_chosen[1]; 
				ss_enh_out	<< endl << "slipdir_chosen: " << setw(outputFileWidth) << slipdir_chosen[0] 
							<< setw(outputFileWidth) << slipdir_chosen[1]; 
				ss_enh_out	<< endl << "tangent_chosen: " << setw(outputFileWidth) << tangent_chosen[0] 
							<< setw(outputFileWidth) << tangent_chosen[1]; 
				if ( model_type == 0 || model_type == 1 ) // M-C model
				{
					ss_enh_out	<< endl << "mu_dir: " << setw(outputFileWidth) << mu_dir[0] 
							<< setw(outputFileWidth) << mu_dir[1];
				}
			}
			else if (NumSD() == 3)	
			{
				ss_enh_out	<< endl << " normal_chosen: " << setw(outputFileWidth) << normal_chosen[0] 
						<< setw(outputFileWidth) << normal_chosen[1] <<  setw(outputFileWidth) << normal_chosen[2]; 
				ss_enh_out	<< endl << "slipdir_chosen: " << setw(outputFileWidth) << slipdir_chosen[0] 
							<< setw(outputFileWidth) << slipdir_chosen[1] <<  setw(outputFileWidth) << slipdir_chosen[2]; 
				ss_enh_out	<< endl << "tangent_chosen: " << setw(outputFileWidth) << tangent_chosen[0] 
							<< setw(outputFileWidth) << tangent_chosen[1] <<  setw(outputFileWidth) << tangent_chosen[2]; 
				if ( model_type == 0 || model_type == 1 ) // M-C model
				{
					ss_enh_out	<< endl << "mu_dir: " << setw(outputFileWidth) << mu_dir[0] 
							<< setw(outputFileWidth) << mu_dir[1] <<  setw(outputFileWidth) << mu_dir[2];
				}
			}
			
			if ( model_type == 0 || model_type == 1 ) // M-C model
			{
				ss_enh_out	<< endl << endl << "loc_flag" 
							<< setw(outputFileWidth) << "jump_displ" 
							<< setw(outputFileWidth) << "gamma_delta" 
							<< setw(outputFileWidth) << "Q_S"
							<< setw(outputFileWidth) << "P_S" 
							<< setw(outputFileWidth) << "q_St"
							<< setw(outputFileWidth) << "sign_q_St"
							<< setw(outputFileWidth) << "YieldTrial"; 
				ss_enh_out	<< endl << fElementLocFlag[elem] 
							<< setw(outputFileWidth) << fElementLocScalars[kNUM_SCALAR_TERMS*elem + kzeta] 
							<< setw(outputFileWidth) << fElementLocScalars[kNUM_SCALAR_TERMS*elem + kgamma_delta] 
							<< setw(outputFileWidth) << fElementLocScalars[kNUM_SCALAR_TERMS*elem + kQ_S]
							<< setw(outputFileWidth) << fElementLocScalars[kNUM_SCALAR_TERMS*elem + kP_S] 
							<< setw(outputFileWidth) << fElementLocScalars[kNUM_SCALAR_TERMS*elem + kq_St]
							<< setw(outputFileWidth) << fElementLocScalars[kNUM_SCALAR_TERMS*elem + ksign_q_St] 
							<< setw(outputFileWidth) << fElementYieldTrial[elem]; 
				
				ss_enh_out	<< endl << endl << "cohesion" 
							<< setw(outputFileWidth) << "friction (rad)" 
							<< setw(outputFileWidth) << "dilation (rad)"; 
				ss_enh_out	<< endl << fElementLocInternalVars[kNUM_ISV_TERMS*elem + kCohesion] 
							<< setw(outputFileWidth) << fElementLocInternalVars[kNUM_ISV_TERMS*elem + kFriction] 
							<< setw(outputFileWidth) << fElementLocInternalVars[kNUM_ISV_TERMS*elem + kDilation] << endl;	
			}
			
			for (int i=0; i < NumElementNodes(); i++)
			{		
				ss_enh_out	<< endl << "activity of node " << i+1
							<< setw(outputFileWidth)<< fElementLocNodesActive[NumElementNodes()*elem + i];
			}			
			dArrayT elem_centroid(NumSD());
			fElementCentroid.RowCopy(elem, elem_centroid);
			if (NumSD() == 2)	
			{
				ss_enh_out	<< endl << "elem_centroid: " << setw(outputFileWidth) << elem_centroid[0] 
							<< setw(outputFileWidth) << elem_centroid[1] << endl;
			}
			else if (NumSD() == 3)	
			{
				ss_enh_out	<< endl << "elem_centroid: " << setw(outputFileWidth) << elem_centroid[0] 
							<< setw(outputFileWidth) << elem_centroid[1] 
							<< setw(outputFileWidth) << elem_centroid[2] << endl;
			}
		}
				
	} // if loc_flag == 2

}

/* form the element stiffness matrix */
void SmallStrainEnhLocT::FormStiffness(double constK)
{
	/* matrix format */
	dMatrixT::SymmetryFlagT format =
		(fLHS.Format() == ElementMatrixT::kNonSymmetric) ?
		dMatrixT::kWhole :
		dMatrixT::kUpperOnly;

	/* integrate element stiffness */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	/********DEBUG*******/
	bool print = false; 
	int pos = fElementCards.Position(); 
	if (pos == 1&&0) print = true; 
	/*******************/
	
	/* current element info */
	int elem = CurrElementNumber();
	loc_flag = fElementLocFlag[elem];
	double vol = fElementVolume[elem];
	double divide_vol = 0.0;
	if ( vol > verysmallnum ) divide_vol = 1.0/vol;
	//int iter = fSSEnhLocMatSupport->IterationNumber();
	//int iter = fSSMatSupport->IterationNumber();
	int iter = IterationNumber();
	
	dArrayT tmp_vec(NumElementNodes()*NumDOF());

	dSymMatrixT F_mun(NumSD());
	dMatrixT tmp_matrix(NumSD());
	
	double K_zetazeta, fYieldTrial;
	
	dSymMatrixT tmp_stress(NumSD());
	
	dArray2DT grad_enh_IPs(NumIP(),NumSD());

	/* element has localized and has been traced, thus fetch data to modify the stiffness matrix */
	if ( loc_flag == 2 ) 
	{
		K_zetazeta = fElementLocScalars[kNUM_SCALAR_TERMS*elem + kKzetazeta];
		
		fElementLocKdzeta.RowCopy(elem, fK_dzeta);
		
		fK_zetad = 0.0;
		
		fYieldTrial = fElementYieldTrial[elem];
		
		/* fetch normal for element */
		fElementLocNormal.RowCopy(elem, normal_chosen);
		
		/* fetch MuDir for element */
		if ( model_type == 0 || model_type == 1 ) // M-C model
		{
			fElementLocMuDir.RowCopy(elem, mu_dir);
			tmp_matrix.Outer(mu_dir, normal_chosen);
			F_mun.Symmetrize(tmp_matrix);
			F_mun.ScaleOffDiagonal(2.0);
		}
		
		if (fDeBug) fElementLocGradEnh.RowCopy(elem, grad_enh_IPs);
	}
	
	fShapes->TopIP();
	while ( fShapes->NextIP() )
	{
		double scale_no_constK = (*Det++)*(*Weight++);
		double scale = constK*scale_no_constK;
	
		/* strain displacement matrix */
		if (fStrainDispOpt == kMeanDilBbar)
			Set_B_bar(fShapes->Derivatives_U(), fMeanGradient, fB);
		else
			Set_B(fShapes->Derivatives_U(), fB);
			
		/*
		check if element has traced and localized
		*/
		if ( loc_flag == 2 ) 
		{			
			/* get De matrix */
			fDe.SetToScaled(scale_no_constK, fCurrMaterial->ce_ijkl());
			if (print) cout << "\nelas_modulus: "<<fCurrMaterial->ce_ijkl();
			
			/* multiply b(transpose) * db, taking account of symmetry, */
			/* and accumulate in elstif */
			fLHS.MultQTBQ(fB, fDe, format, dMatrixT::kAccumulate);
						
			// calculate separate coupling stiffness matrices
			if ( fYieldTrial > verysmallnum)
			{
				// calc fK_zetad
				//fK_zetad.MultATBC(F_mun_nonsym, fDe, fB, format, dMatrixT::kAccumulate);
				tmp_stress = 0.0;
				fDe.Multx(F_mun, tmp_stress);
				tmp_vec = 0.0;
				fB.MultTx(tmp_stress, tmp_vec);
				tmp_vec *= divide_vol;
				fK_zetad += tmp_vec;

				if (fDeBug)
				{
					const int ip = fShapes->CurrIP()+1;
					int array_location = fShapes->CurrIP()*NumSD();
					if ( NumSD() == 2 )
					{
						ss_enh_out	<< endl << "GradEnh at IP" << ip
								<< setw(outputFileWidth) << grad_enh_IPs[array_location] 
								<< setw(outputFileWidth) << grad_enh_IPs[array_location+1];
					}
					else if ( NumSD() == 3 )
					{
						ss_enh_out	<< endl << "GradEnh at IP" << ip
								<< setw(outputFileWidth) << grad_enh_IPs[array_location] 
								<< setw(outputFileWidth) << grad_enh_IPs[array_location+1] 
								<< setw(outputFileWidth) << grad_enh_IPs[array_location+2];
					}
				}
			}
			
		}
		else // update for standard plasticity
		{
			/* get D matrix */
			fD.SetToScaled(scale, fCurrMaterial->c_ijkl());
			if (print) cout << "\nmodulus: "<<fCurrMaterial->c_ijkl();
			
			/* multiply b(transpose) * db, taking account of symmetry, */
			/* and accumulate in elstif */
			fLHS.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);	
		}
		
	} // while next IP
	
	if ( loc_flag == 2 )
	{
		// store fK_zetad
		fElementLocKzetad.SetRow(elem, fK_zetad);
		
		if ( fYieldTrial > verysmallnum )
		{
			// modify stiffness matrix
			fK_dzeta_x_Kzetad.Outer(fK_dzeta, fK_zetad);
			
			if ( fabs(K_zetazeta) > verysmallnum )
			{
				fK_dzeta_x_Kzetad /= K_zetazeta;
			}
			else
			{
				fK_dzeta_x_Kzetad = 0.0;
			}
			
			fLHS -= fK_dzeta_x_Kzetad;
		}
		
	}
	
}

/* compute the measures of strain/deformation over the element */
void SmallStrainEnhLocT::SetGlobalShape(void)
{
	/* current element number */
	int elem = CurrElementNumber();
	
	/* inherited */
	SolidElementT::SetGlobalShape();

	/* compute mean of shape function gradients */
	double& vol = fElementVolume[elem];
	SetMeanGradient(fMeanGradient, vol);
	
	/* last deformed volume */
	double& vol_last = fElementVolume_last[elem];
	
	/* material information */
	int material_number = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[material_number];
	
	/* using B-bar */
	if (fStrainDispOpt == kMeanDilBbar)
	{
		/* compute mean of shape function gradients */
		//SetMeanGradient(fMeanGradient, vol);

		/* loop over integration points */
		fShapes->TopIP();
		while (fShapes->NextIP())
		{
			/* set B-bar */
			int ip = fShapes->CurrIP();
			Set_B_bar(fShapes->Derivatives_U(ip), fMeanGradient, fB);
	
			/* deformation gradient */
			if (needs[fNeedsOffset + kstrain])
			{
				/* transpose displacement array */
				fLocDisp.ReturnTranspose(fLocDispTranspose);

				/* compute strain using B-bar */
				dSymMatrixT& strain = fStrain_List[ip];
				fB.Multx(fLocDispTranspose, strain);
				strain.ScaleOffDiagonal(0.5);
			}

			/* "last" deformation gradient */
			if (needs[fNeedsOffset + kstrain_last])
			{
				/* transpose displacement array */
				fLocLastDisp.ReturnTranspose(fLocDispTranspose);

				/* compute strain using B-bar */
				dSymMatrixT& strain = fStrain_last_List[ip];
				fB.Multx(fLocDispTranspose, strain);
				strain.ScaleOffDiagonal(0.5);
			}
		}		
	}
	else
	{
		/* loop over integration points */
		for (int i = 0; i < NumIP(); i++)
		{
			/* deformation gradient */
			if (needs[fNeedsOffset + kstrain])
			{
				/* displacement gradient */
				fShapes->GradU(fLocDisp, fGradU, i);

				/* symmetric part */
				fStrain_List[i].Symmetrize(fGradU);
			}

			/* "last" deformation gradient */
			if (needs[fNeedsOffset + kstrain_last])
			{
				/* displacement gradient */
				fShapes->GradU(fLocLastDisp, fGradU, i);

				/* symmetric part */
				fStrain_last_List[i].Symmetrize(fGradU);
			}
		}
	}
}

/***********************************************************************
 * Private
 ***********************************************************************/

/* compute mean shape function gradient, Hughes (4.5.23) */
void SmallStrainEnhLocT::SetMeanGradient(dArray2DT& mean_gradient, double& vol) const
{
	int nip = NumIP();
	const double* det = fShapes->IPDets();
	const double*   w = fShapes->IPWeights();

	/* volume */
	vol = 0.0;
	for (int i = 0; i < nip; i++)
		vol += w[i]*det[i];

	/* initialize */
	mean_gradient = 0.0;			

	/* integrate */
	for (int i = 0; i < nip; i++)
		mean_gradient.AddScaled(w[i]*det[i]/vol, fShapes->Derivatives_U(i));
}

