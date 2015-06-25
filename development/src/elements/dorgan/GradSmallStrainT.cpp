/* $id: GradSmallStrainT.cpp,v 1.10 2004/06/17 20:45:13 rdorgan Exp $ */ 
#include "GradSmallStrainT.h"
#include "ShapeFunctionT.h"
#include "GradSSSolidMatT.h"

/* shape functions */
#include "ShapeTools.h"
#include "C0_LineT.h"
#include "C1_LineT.h"
#include "C0_QuadT.h"
#include "C1_QuadT.h"
#include "C0_HexahedronT.h"

#include <cmath>
#include "ifstreamT.h"
#include "ofstreamT.h"

/* materials lists */
#include "GradSSSolidMatList1DT.h"
#include "GradSSSolidMatList2DT.h"
#include "GradSSSolidMatList3DT.h"
#include "GradSSMatSupportT.h"
#include "ParameterContainerT.h"
#include "ModelManagerT.h"
#include "NodeManagerT.h"
#include "KBC_ControllerT.h"

#include "FieldT.h"
#include "FieldSupportT.h"

using namespace Tahoe;

/* parameters */
const double kYieldTol    = 1.0e-10;

/* constructor */
GradSmallStrainT::GradSmallStrainT(const ElementSupportT& support):
	SmallStrainT(support), //pass the displacement field to the base class
	fGradSSMatSupport(NULL),

	fMaxWeakened(-1),
	fWeakened(0),
	fHoldTime(false),

	fK_bb(ElementMatrixT::kNonSymmetric),
	fK_bh(ElementMatrixT::kNonSymmetric),
	fK_hb(ElementMatrixT::kNonSymmetric),
	fK_hh(ElementMatrixT::kNonSymmetric),
	fK_hp(ElementMatrixT::kNonSymmetric),
	fK_hq(ElementMatrixT::kNonSymmetric),
	fK_ct(ElementMatrixT::kNonSymmetric),

	fLocPMultiplier(LocalArrayT::kDisp),
	fLocLastPMultiplier(LocalArrayT::kLastDisp),

	fDisplacement(NULL),
	fPMultiplier(NULL),
	
	fShapes_PMultiplier(NULL),
	
	fNumSD(-1),
	fNumIP_Disp(-1),
	fNumElementNodes_Disp(-1),
	fNumDOF_Disp(-1),
	fNumDOF_PMultiplier(-1),
	fNumEQ_Total(-1),

	fNumIP_PMultiplier(-1),
	fNumElementNodes_PMultiplier(-1),
	fDegreeOfContinuity_PMultiplier(-1),
	fNodalConstraint(-1),

	fprint_GlobalShape(false),
	fprint_Kd(false),
	fprint_KdMatrix(false),
	fprint_Stiffness(false),
	fprint_StiffnessMatrix(false),
	fprint_All(false),

	fI(1)
{
	SetName("grad_small_strain");
}

/* destructor */
GradSmallStrainT::~GradSmallStrainT(void) 
{  
	delete fGradSSMatSupport;
	delete fShapes_PMultiplier;
}

void GradSmallStrainT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
								 AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
#pragma unused(eq_2)

	const char caller[] = "GradSmallStrainT::Equations";
	
#if __option(extended_errorcheck)
	if (fConnectivities.Length() != fEqnos.Length()) throw ExceptionT::kSizeMismatch;
	if (fConnectivities_PMultiplier.Length() != fEqnos.Length()) throw ExceptionT::kSizeMismatch;
#endif

	/* loop over connectivity blocks */
	for (int i = 0; i < fEqnos.Length(); i++)
	{
		/* connectivities */
		const iArray2DT& connects = *(fConnectivities[i]);
		const iArray2DT& connects_pmultiplier = fConnectivities_PMultiplier[i];
		int fNumElements = connects.MajorDim();

		/* dimension */
		fEqnos[i].Dimension   (fNumElements, fNumEQ_Total);
		iArray2DT fEqnos_Disp (fNumElements, fNumElementNodes_Disp * fNumDOF_Disp );
		iArray2DT fEqnos_PMultiplier(fNumElements, fNumElementNodes_PMultiplier*fNumDOF_PMultiplier);

		/* get equation numbers */
		fDisplacement->SetLocalEqnos(connects, fEqnos_Disp);
		fPMultiplier->SetLocalEqnos(connects_pmultiplier, fEqnos_PMultiplier);

		/* write into one array */
		fEqnos[i].BlockColumnCopyAt(fEqnos_Disp, 0);
		fEqnos[i].BlockColumnCopyAt(fEqnos_PMultiplier, fEqnos_Disp.MinorDim());

		/* add to list of equation numbers */
		eq_1.Append(&fEqnos[i]);
	}

	/* reset pointers to element cards */
	SetElementCards(fBlockData, fConnectivities, fEqnos, fElementCards);
}

/* implementation of the ParameterInterfaceT interface */
void GradSmallStrainT::DefineParameters(ParameterListT& list) const
{
	const char caller[] = "GradSmallStrainT::DefineParameters";

	/* inherited */
	SmallStrainT::DefineParameters(list);

	/* associate field */
	list.AddParameter(ParameterT::Word, "grad_plasticity_field_name");

	ParameterT degree_of_continuity_pmultiplier(ParameterT::Enumeration, "degree_of_continuity_pmultiplier");
	degree_of_continuity_pmultiplier.AddEnumeration("C0", kC0);
	degree_of_continuity_pmultiplier.AddEnumeration("C1", kC1);
    degree_of_continuity_pmultiplier.SetDefault(kC1);
	list.AddParameter(degree_of_continuity_pmultiplier);

	ParameterT max_weakened_ip(fMaxWeakened, "max_weakened_ip");
	max_weakened_ip.SetDefault(1000);
	list.AddParameter(max_weakened_ip);

	ParameterT nodal_constraint(fNodalConstraint, "nodal_constraint");
	nodal_constraint.SetDefault(0.0);
	list.AddParameter(nodal_constraint);

	ParameterT print_globalShape(fprint_GlobalShape, "print_globalShape");
	print_globalShape.SetDefault(fprint_GlobalShape);
	list.AddParameter(print_globalShape);

	ParameterT print_kd(fprint_Kd, "print_kd");
	print_kd.SetDefault(fprint_Kd);
	list.AddParameter(print_kd);

	ParameterT print_kd_matrix(fprint_KdMatrix, "print_kd_matrix");
	print_kd_matrix.SetDefault(fprint_KdMatrix);
	list.AddParameter(print_kd_matrix);

	ParameterT print_stiffness(fprint_Stiffness, "print_stiffness");
	print_stiffness.SetDefault(fprint_Stiffness);
	list.AddParameter(print_stiffness);

	ParameterT print_stiffness_matrix(fprint_StiffnessMatrix, "print_stiffness_matrix");
	print_stiffness_matrix.SetDefault(fprint_StiffnessMatrix);
	list.AddParameter(print_stiffness_matrix);

	ParameterT print_all(fprint_All, "print_all");
	print_all.SetDefault(fprint_All);
	list.AddParameter(print_all);
}

/* information about subordinate parameter lists */
void GradSmallStrainT::DefineSubs(SubListT& sub_list) const
{
	const char caller[] = "GradSmallStrainT::DefineSubs";

	/* inherited */
	SolidElementT::DefineSubs(sub_list);	

	/* element block/material specification */
	sub_list.AddSub("grad_small_strain_element_block", ParameterListT::OnePlus);
}

/* return the description of the given inline subordinate parameter list */
ParameterInterfaceT* GradSmallStrainT::NewSub(const StringT& name) const
{
	const char caller[] = "GradSmallStrainT::NewSub";

	if (name == "grad_small_strain_element_block")
	{
		ParameterContainerT* block = new ParameterContainerT(name);
		
		/* list of element block ID's (defined by ElementBaseT) */
		block->AddSub("block_ID_list", ParameterListT::Once);
	
		/* choice of materials lists (inline) */
		block->AddSub("grad_small_strain_material_choice", ParameterListT::Once, true);
	
		/* set this as source of subs */
		block->SetSubSource(this);
		
		return block;
	}
	else /* inherited */
		return SolidElementT::NewSub(name);
}

/* return the description of the given inline subordinate parameter list. */
void GradSmallStrainT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	const char caller[] = "GradSmallStrainT::DefineInlineSub";

	if (name == "grad_small_strain_material_choice")
	{
		order = ParameterListT::Choice;
		
		/* list of choices */
		sub_lists.AddSub("grad_small_strain_material_1D");
		sub_lists.AddSub("grad_small_strain_material_2D");
		sub_lists.AddSub("grad_small_strain_material_3D");
	}
	else /* inherited */
		SolidElementT::DefineInlineSub(name, order, sub_lists);
}

void GradSmallStrainT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "GradSmallStrainT::TakeParameterList";

	/* get the displacement field */
	const StringT& field_name = list.GetParameter("field_name");
	fDisplacement = ElementSupport().Field(field_name);
	if (!fDisplacement)
		ExceptionT::GeneralFail(caller, "could not resolve \"%s\" field", field_name.Pointer());

	/* get the plastic multiplier field */
	const StringT& grad_plasticity_field_name = list.GetParameter("grad_plasticity_field_name");
	fPMultiplier = ElementSupport().Field(grad_plasticity_field_name);
	if (!fPMultiplier)
		ExceptionT::GeneralFail(caller, "could not resolve \"%s\" grad_plasticity_field", grad_plasticity_field_name.Pointer());

	/* get the element parameters */
	int b = list.GetParameter("degree_of_continuity_pmultiplier");
	fDegreeOfContinuity_PMultiplier = (b == kC1) ? kC1 : kC0;
	fMaxWeakened               = list.GetParameter("max_weakened_ip");
	fNodalConstraint           = list.GetParameter("nodal_constraint");
	fprint_GlobalShape         = list.GetParameter("print_globalShape");
	fprint_Kd                  = list.GetParameter("print_kd");
	fprint_KdMatrix            = list.GetParameter("print_kd_matrix");
	fprint_Stiffness           = list.GetParameter("print_stiffness");
	fprint_StiffnessMatrix     = list.GetParameter("print_stiffness_matrix");
	fprint_All                 = list.GetParameter("print_all");

	if (fprint_All == true)
		fprint_GlobalShape = fprint_Kd = fprint_KdMatrix = fprint_Stiffness = fprint_StiffnessMatrix = true;

	/* dimensions */
	fNumSD = NumSD();
	fNumDOF_Disp = fDisplacement->NumDOF();
	fNumDOF_PMultiplier = fPMultiplier->NumDOF();

	/* vertex nodes only */
	if (fNumSD == 1)
		fNumElementNodes_PMultiplier = 2;  // line
	else if (fNumSD == 2)
		fNumElementNodes_PMultiplier = 4;  // quad
	else 
		fNumElementNodes_PMultiplier = 8;  // hex
	
	/* inherited */
	SmallStrainT::TakeParameterList(list);
	
	/* dimensions */
	fNumElementNodes_Disp = NumElementNodes();
	fNumIP_Disp = NumIP();
	fNumIP_PMultiplier = fNumIP_Disp;
	
	/* allocate lists */
	fPMultiplier_List.Dimension(fNumIP_PMultiplier);           // pmultiplier
	fPMultiplier_last_List.Dimension(fNumIP_PMultiplier);      // "last" pmultiplier
	fGradPMultiplier_List.Dimension(fNumIP_PMultiplier);       // gradient pmultiplier
	fGradPMultiplier_last_List.Dimension(fNumIP_PMultiplier);  // "last" gradient pmultiplier
	fLapPMultiplier_List.Dimension(fNumIP_PMultiplier);        // Laplacian pmultiplier
	fLapPMultiplier_last_List.Dimension(fNumIP_PMultiplier);   // "last" Laplacian pmultiplier
	fYield_List.Dimension(fNumIP_Disp);                        // yield condition pmultiplier

	/* additional dimensions for PMultiplier */
	for (int i = 0; i < NumIP(); i++)
	{
		fPMultiplier_List[i].Dimension(1);
		fPMultiplier_last_List[i].Dimension(1);
		fGradPMultiplier_List[i].Dimension(fNumSD,1);
		fGradPMultiplier_last_List[i].Dimension(fNumSD,1);
		fLapPMultiplier_List[i].Dimension(1);
		fLapPMultiplier_last_List[i].Dimension(1);
	}

	/* dimension moduli in workspace */
	fDM_bb.Dimension(dSymMatrixT::NumValues(fNumSD));
	fODM_bh.Dimension(dSymMatrixT::NumValues(fNumSD), 1);
	fODM_hb.Dimension(1, dSymMatrixT::NumValues(fNumSD));
	fGM_hh.Dimension(1);
	fGM_hp.Dimension(1, fNumSD);
	fGM_hq.Dimension(1);

	/* dimension stiffness matrices in workspace */
	fK_bb.Dimension(fNumElementNodes_Disp * fNumDOF_Disp, fNumElementNodes_Disp * fNumDOF_Disp );
	fK_bh.Dimension(fNumElementNodes_Disp * fNumDOF_Disp, fNumElementNodes_PMultiplier * fNumDOF_PMultiplier);
	fK_hb.Dimension(fNumElementNodes_PMultiplier * fNumDOF_PMultiplier, fNumElementNodes_Disp * fNumDOF_Disp );
	fK_hh.Dimension(fNumElementNodes_PMultiplier * fNumDOF_PMultiplier, fNumElementNodes_PMultiplier * fNumDOF_PMultiplier);
	fK_hp.Dimension(fNumElementNodes_PMultiplier * fNumDOF_PMultiplier, fNumElementNodes_PMultiplier * fNumDOF_PMultiplier);
	fK_hq.Dimension(fNumElementNodes_PMultiplier * fNumDOF_PMultiplier, fNumElementNodes_PMultiplier * fNumDOF_PMultiplier);
	fK_ct.Dimension(fNumElementNodes_PMultiplier * fNumDOF_PMultiplier, fNumElementNodes_PMultiplier * fNumDOF_PMultiplier);

	/* resize work arrays */
	fNumEQ_Total = fNumElementNodes_Disp*fNumDOF_Disp+fNumElementNodes_PMultiplier*fNumDOF_PMultiplier;
	fRHS.Dimension(fNumEQ_Total);
	fLHS.Dimension(fRHS.Length());

	/* allocate shape functions */
	fh.Dimension (1, fNumElementNodes_PMultiplier*fNumDOF_PMultiplier);
	fhT.Dimension(fNumElementNodes_PMultiplier*fNumDOF_PMultiplier, 1);
	fp.Dimension (fNumSD, fNumElementNodes_PMultiplier*fNumDOF_PMultiplier);
	fq.Dimension (1, fNumElementNodes_PMultiplier*fNumDOF_PMultiplier);
}

/* extract the list of material parameters */
void GradSmallStrainT::CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const
{
	const char caller[] = "GradSmallStrainT::CollectMaterialInfo";

	/* initialize */
	mat_params.Clear();
	
	/* collected material parameters */
	int num_blocks = all_params.NumLists("grad_small_strain_element_block");
	for (int i = 0; i < num_blocks; i++) {

		/* block information */	
		const ParameterListT& block = all_params.GetList("grad_small_strain_element_block", i);
		
		/* resolve material list name */
		if (i == 0) {
			const ParameterListT& mat_list_params = block.GetListChoice(*this, "grad_small_strain_material_choice");
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
MaterialSupportT* GradSmallStrainT::NewMaterialSupport(MaterialSupportT* p) const
{
	const char caller[] = "GradSmallStrainT::NewMaterialSupport";

	/* allocate */
	if (!p) p = new GradSSMatSupportT(fDisplacement->NumDOF(), fPMultiplier->NumDOF(), NumIP(), NumIP());

	/* inherited initializations */
	SmallStrainT::NewMaterialSupport(p);
	
	/* set SolidMatSupportT fields */
	GradSSMatSupportT* ps = TB_DYNAMIC_CAST(GradSSMatSupportT*, p);
	if (ps) {
		ps->SetLinearPMultiplier(&fPMultiplier_List);
		ps->SetLinearPMultiplier_last(&fPMultiplier_last_List);

		ps->SetLinearGradPMultiplier(&fGradPMultiplier_List);
		ps->SetLinearGradPMultiplier_last(&fGradPMultiplier_last_List);

		ps->SetLinearLapPMultiplier(&fLapPMultiplier_List);
		ps->SetLinearLapPMultiplier_last(&fLapPMultiplier_last_List);
	}

	return p;
}

/* return a pointer to a new material list */
MaterialListT* GradSmallStrainT::NewMaterialList(const StringT& name, int size)
{
	const char caller[] = "GradSmallStrainT::NewMaterialList";

	/* resolve dimension */
	int nsd = -1;
	if (name == "grad_small_strain_material_1D") nsd = 1;
	else if (name == "grad_small_strain_material_2D") nsd = 2;
	else if (name == "grad_small_strain_material_3D") nsd = 3;
	
	/* no match */
	if (nsd == -1) return NULL;

	/* full list */
	if (size > 0)
	{
		/* material support */
		if (!fGradSSMatSupport)
		{
			fGradSSMatSupport = TB_DYNAMIC_CAST(GradSSMatSupportT*, NewMaterialSupport());
			if (!fGradSSMatSupport)
				ExceptionT::GeneralFail("GradSmallStrainT::NewMaterialList");
		}

		if (nsd == 1)
			return new GradSSSolidMatList1DT(size, *fGradSSMatSupport);
		else if (nsd == 2)
			return new GradSSSolidMatList2DT(size, *fGradSSMatSupport);
 		else if (nsd == 3)
 			return new GradSSSolidMatList3DT(size, *fGradSSMatSupport);
	}
else
	{
		if (nsd == 1)
			return new GradSSSolidMatList1DT;
		else if (nsd == 2)
			return new GradSSSolidMatList2DT;
		else if (nsd == 3)
 			return new GradSSSolidMatList3DT;
	}
	
	/* no match */
	return NULL;
}

/* define the elements blocks for the element group */
void GradSmallStrainT::DefineElements(const ArrayT<StringT>& block_ID, const ArrayT<int>& mat_index)
{
	const char caller[] = "GradSmallStrainT::DefineElements";

	/* initialize displacement connectivities */
	SmallStrainT::DefineElements(block_ID, mat_index);

	fConnectivities_PMultiplier.Dimension(fConnectivities.Length());

	int num_vertex_nodes = -1;

	/* number of vertex nodes */
	if (fNumSD == 1)
		num_vertex_nodes = 2;  // line
	else if (fNumSD == 2)
		num_vertex_nodes = 4;  // quad
	else
		num_vertex_nodes = 8;  // quad
		
	/* list of the vertex nodes */
	iArrayT vertex_nodes(num_vertex_nodes);

	for (int nd = 0; nd < num_vertex_nodes; nd ++)
		vertex_nodes[nd] = nd;

	fConnectivities_All.Dimension(NumElements(), num_vertex_nodes);
	fFixedPMultiplier.Dimension(fConnectivities.Length());
	fFixedPMultiplier = NULL;

	/* translate blocks */
	int count = 0;
	for (int i = 0; i < fConnectivities.Length(); i++)
	{
		const iArray2DT& connects = *(fConnectivities[i]);
		iArray2DT& connects_pmultiplier = fConnectivities_PMultiplier[i];
		connects_pmultiplier.Alias(connects.MajorDim(), num_vertex_nodes, fConnectivities_All(count));
		
		/* extract */
		for (int j = 0; j < vertex_nodes.Length(); j++)
			connects_pmultiplier.ColumnCopy(j, connects, vertex_nodes[j]);
				
		/* next block */
		count += connects.MajorDim();

		if (NumElementNodes() > num_vertex_nodes)
		{
			/* prescribe fixed multiplier field at center node */
			FieldT* non_constPMultiplier = const_cast<FieldT*>(fPMultiplier);

			/* construct new contoller */
			fFixedPMultiplier[i] = ElementSupport().NodeManager().NewKBC_Controller(*non_constPMultiplier, KBC_ControllerT::kPrescribed);

			/* add to field */
			non_constPMultiplier->AddKBCController(fFixedPMultiplier[i]);

			/* define fixed conditions */
			ArrayT<KBC_CardT>& KBC_cards = fFixedPMultiplier[i]->KBC_Cards();
			KBC_cards.Dimension(connects.MajorDim()*fNumDOF_PMultiplier);
			for (int fixed_node = num_vertex_nodes; fixed_node < NumElementNodes(); fixed_node ++)
			{
				int dex = 0;
				for (int j = 0; j < fNumDOF_PMultiplier; j++)
					for (int i = 0; i < connects.MajorDim(); i++)
						KBC_cards[dex++].SetValues(connects(i,fixed_node), j, KBC_CardT::kFix, NULL, 0.0);
			}
		}
		else
			/* same connectivities - make aliases */
			fConnectivities_PMultiplier[i].Alias(*(fConnectivities[i]));
	}
}

/* initialize local arrays */
void GradSmallStrainT::SetLocalArrays(void)
{
	const char caller[] = "GradSmallStrainT::SetLocalArrays";

	/* inherited */
	SmallStrainT::SetLocalArrays();

	/* dimension local arrays */
	fLocPMultiplier.Dimension(fNumElementNodes_PMultiplier, fPMultiplier->NumDOF());
	fLocLastPMultiplier.Dimension(fNumElementNodes_PMultiplier, fPMultiplier->NumDOF());

	fLocPMultiplierTranspose.Dimension(fLocPMultiplier.Length());
	
	/* register local arrays */
	fPMultiplier->RegisterLocal(fLocPMultiplier);
	fPMultiplier->RegisterLocal(fLocLastPMultiplier);
}

/* initialization functions */
void GradSmallStrainT::SetShape(void)
{
	const char caller[] = "GradSmallStrainT::SetShape";

	/* inherited */
	SmallStrainT::SetShape();

	/* select shape function tools to use */
	if (fNumSD == 1)
		if (fDegreeOfContinuity_PMultiplier == kC0)
			/* 1D, C0 */
			fShapes_PMultiplier = new C0_LineT(GeometryCode(), NumIP(), fNumElementNodes_PMultiplier, fPMultiplier->NumDOF(), fLocInitCoords);
		else
			/* 1D, C1 */
			fShapes_PMultiplier = new C1_LineT(GeometryCode(), NumIP(), fNumElementNodes_PMultiplier, fPMultiplier->NumDOF(), fLocInitCoords);
	else if (fNumSD == 2)
		if (fDegreeOfContinuity_PMultiplier == kC0)
			/* 2D, C0 */
			fShapes_PMultiplier = new C0_QuadT(GeometryCode(), NumIP(), fNumElementNodes_PMultiplier, fPMultiplier->NumDOF(), fLocInitCoords);
		else
			/* 2D, C1 */
			fShapes_PMultiplier = new C1_QuadT(GeometryCode(), NumIP(), fNumElementNodes_PMultiplier, fPMultiplier->NumDOF(), fLocInitCoords);
 	else
 		if (fDegreeOfContinuity_PMultiplier == kC0)
			/* 3D, C0 */
  			fShapes_PMultiplier = new C0_HexahedronT(GeometryCode(), NumIP(), fNumElementNodes_PMultiplier, fPMultiplier->NumDOF(), fLocInitCoords);
 		else
			/* 3D, C1 */
 			ExceptionT::GeneralFail(caller, "C1_HexahedronT not implemented");		

	if (!fShapes_PMultiplier) throw ExceptionT::kOutOfMemory;

	/* initialize */
	fShapes_PMultiplier->Initialize();
}

/* form shape functions and derivatives */
void GradSmallStrainT::SetGlobalShape(void)
{
	const char caller[] = "GradSmallStrainT::SetGlobalShape";

	/* inherited */
	SmallStrainT::SetGlobalShape();

	/* collect element values of PMultiplier */
	fConnectivities_All.RowAlias(CurrElementNumber(), fNodesPMultiplier);
	fLocPMultiplier.SetLocal(fNodesPMultiplier);
	fLocLastPMultiplier.SetLocal(fNodesPMultiplier);

	/* compute shape function derivatives */
	fShapes_PMultiplier->SetDerivatives();

	/********DEBUG*******/
	if (fprint_GlobalShape)
	{
		cout << caller << endl;
		cout << "  element: " << CurrElementNumber() << endl;
	}
	/*******************/

	fShapes->TopIP();
	/* loop over integration points */
	while(fShapes->NextIP())
	{
		int ip = fShapes->CurrIP();
			
		/* accumulate shape functions */
		Set_h(fh);
		Set_p(fp);
		Set_q(fq);

		/* compute current fields using pmultiplier shape functions */
		fLocPMultiplier.ReturnTranspose(fLocPMultiplierTranspose);

		dMatrixT& pmultiplier     = fPMultiplier_List[ip];
		dMatrixT& gradpmultiplier = fGradPMultiplier_List[ip];
		dMatrixT& plapmultiplier  = fLapPMultiplier_List[ip];

		fh.Multx(fLocPMultiplierTranspose, pmultiplier);
		fp.Multx(fLocPMultiplierTranspose, gradpmultiplier);
		fq.Multx(fLocPMultiplierTranspose, plapmultiplier);

		/* compute "last" fields using pmultiplier shape functions */
		fLocLastPMultiplier.ReturnTranspose(fLocPMultiplierTranspose);

		dMatrixT& pmultiplier_last     = fPMultiplier_last_List[ip];
		dMatrixT& gradpmultiplier_last = fGradPMultiplier_last_List[ip];
		dMatrixT& plapmultiplier_last  = fLapPMultiplier_last_List[ip];

		fh.Multx(fLocPMultiplierTranspose, pmultiplier_last);
		fp.Multx(fLocPMultiplierTranspose, gradpmultiplier_last);
		fq.Multx(fLocPMultiplierTranspose, plapmultiplier_last);
	}

	/********DEBUG*******/
	if (fprint_GlobalShape)
	{
		for (int nd_dof = 0; nd_dof < fNumElementNodes_PMultiplier*fNumDOF_PMultiplier; nd_dof++)
			cout << "                 fLocPMultiplier[" << nd_dof << "]  : " << fLocPMultiplier[nd_dof] << " * " << fh[nd_dof] << endl;
		cout << endl;
		for (int ip = 0; ip < fNumIP_Disp; ip++)
			cout << "                 fPMultiplier_List[" << ip << "]: " << fPMultiplier_List[ip] << endl;
		for (int ip = 0; ip < fNumIP_Disp; ip++)
			cout << "                 fLapPMultiplier_List[" << ip << "]: " << fLapPMultiplier_List[ip] << endl;
	}
	/*******************/
}

/* form the element stiffness matrix */
void GradSmallStrainT::FormStiffness(double constK)
{
	const char caller[] = "GradSmallStrainT::FormStiffness";
	int curr_group = ElementSupport().CurrentGroup();

	/* dmatrix format */
	dMatrixT::SymmetryFlagT format =
		(fLHS.Format() == ElementMatrixT::kNonSymmetric) ?
		dMatrixT::kWhole:
		dMatrixT::kUpperOnly;

	if(format != dMatrixT::kWhole)
		ExceptionT::BadInputValue(caller, "LHS must be NonSymmetric");

	/* integration rules (same for both fields)*/
	const double* Det = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	/********DEBUG*******/
	if (fprint_Stiffness)
	{
		cout << caller << endl;
		cout << "  element: " << CurrElementNumber() << endl;
	}
	/*******************/

	/* initialize */
	fK_bb = 0.0;
	fK_bh = 0.0;
	fK_hb = 0.0;
	fK_hh = 0.0;
	fK_hp = 0.0;
	fK_hq = 0.0;
	fK_ct = 0.0;

	fShapes->TopIP();
	/* integrate stiffness over the elements */
	while(fShapes->NextIP())
	{
		/* integration weight */
		double scale = constK*(*Det++)*(*Weight++);

		/* accumulate strain displacement matrix */
		if (fStrainDispOpt == kMeanDilBbar)
			Set_B_bar(fShapes->Derivatives_U(), fMeanGradient, fB);
		else
			Set_B(fShapes->Derivatives_U(), fB);

		/* accumulate shape functions */
		Set_h(fh);
		Set_p(fp);
		Set_q(fq);

		/* compute elastic stiffness matrix */
		fDM_bb.SetToScaled(scale, fCurrMaterial_Grad->c_ijkl());
		fK_bb.MultQTBQ(fB, fDM_bb, format, dMatrixT::kAccumulate);

		/* compute off-diagonal matrics */
		fODM_bh.SetToScaled(scale, fCurrMaterial_Grad->odm_bh_ij());
		fK_bh.MultATBC(fB, fODM_bh, fh, format, dMatrixT::kAccumulate);

		fODM_hb.SetToScaled(scale, fCurrMaterial_Grad->odm_hb_ij());
		fK_hb.MultATBC(fh, fODM_hb, fB, format, dMatrixT::kAccumulate);

		/* compute non-symmetric, gradient dependent matrix */
		fGM_hh.SetToScaled(scale, fCurrMaterial_Grad->gm_hh());
		fK_hh.MultQTBQ(fh, fGM_hh, format, dMatrixT::kAccumulate);

		fGM_hp.SetToScaled(scale, fCurrMaterial_Grad->gm_hp());
		fK_hp.MultATBC(fh, fGM_hp, fp, format, dMatrixT::kAccumulate);

		fGM_hq.SetToScaled(scale, fCurrMaterial_Grad->gm_hq());
		fK_hq.MultATBC(fh, fGM_hq, fq, format, dMatrixT::kAccumulate);

		/* obtain yield condition residual */
		fYield_List[CurrIP()] = fCurrMaterial_Grad->yc();

		/********DEBUG*******/
		if (fprint_Stiffness)
		{
			cout<<"            ip: "             << CurrIP()                                                 <<endl;
			cout<<"                 strain         : " << sqrt(fStrain_List[CurrIP()].ScalarProduct())       <<endl;
			cout<<"                 stress         : " << sqrt(fCurrMaterial_Grad->s_ij().ScalarProduct())   <<endl;
			cout<<"                 yc             : " << fCurrMaterial_Grad->yc()                           <<endl;
			cout<<"                 PMultiplier    : " << fPMultiplier_List[CurrIP()]                        <<endl;
			cout<<"                 del_PMultiplier: " << fPMultiplier_List[CurrIP()][0] -  fPMultiplier_last_List[CurrIP()][0]                      <<endl;
			cout<<"                 LapPMultiplier : " << fLapPMultiplier_List[CurrIP()]                     <<endl;
			cout<<"                 odm_bh_ij[0]   : " << (fCurrMaterial_Grad->odm_hb_ij())[0]               <<endl;
			cout<<"                 odm_hb_ij[0]   : " << (fCurrMaterial_Grad->odm_bh_ij())[0]               <<endl;
			cout<<"                 gm_hh          : " << fCurrMaterial_Grad->gm_hh()                        <<endl;
			cout<<"                 gm_hp[0]       : " << fCurrMaterial_Grad->gm_hp()[0]                     <<endl;
			cout<<"                 gm_hq          : " << fCurrMaterial_Grad->gm_hq()                        <<endl;
		}
		/*******************/
	}

	/* assemble into element stiffness matrix */
	fLHS.AddBlock(0,            0,            fK_bb);
	fLHS.AddBlock(0,            fK_bb.Cols(), fK_bh);
	fLHS.AddBlock(fK_bb.Rows(), 0,            fK_hb);
	fLHS.AddBlock(fK_bb.Rows(), fK_bb.Cols(), fK_hh);
	fLHS.AddBlock(fK_bb.Rows(), fK_bb.Cols(), fK_hp);
	fLHS.AddBlock(fK_bb.Rows(), fK_bb.Cols(), fK_hq);

	dArrayT fLocYield(fNumElementNodes_Disp);
	IP_ExtrapolateAll(fYield_List, fLocYield);

	/* add constraint to Krr if elastic */
	fK_ct = 0.;
	for (int nd = 0; nd < fNumElementNodes_PMultiplier ; nd++)
		if (fLocPMultiplier[nd] <= fLocLastPMultiplier[nd] && false)
			fK_ct(nd)[nd] = fNodalConstraint*fNodalConstraint;
	fLHS.AddBlock(fK_bb.Rows(), fK_bb.Cols(), fK_ct);

	/********DEBUG*******/
	if (fprint_All)
	{
		cout << caller << " Stiffness Matrix" << endl;
		cout << "  element: " << CurrElementNumber() << endl;
		for (int i=0; i < fNumElementNodes_Disp * fNumDOF_Disp; i++)
		{
			cout<<"                 K("<< i << ",j): ";
			for (int j=0; j < fNumElementNodes_Disp * fNumDOF_Disp; j++)
				cout << fK_bb(j)[i] << "  ";
			for (int j=0; j < fNumElementNodes_PMultiplier * fNumDOF_PMultiplier; j++)
				cout << fK_bh(j)[i] << "      ";
			cout << endl;
		}
		for (int i=0; i < fNumElementNodes_PMultiplier * fNumDOF_PMultiplier; i++)
		{
			cout<<"                 K("<< i << ",j): ";
			for (int j=0; j < fNumElementNodes_Disp * fNumDOF_Disp; j++)
				cout << fK_hb(j)[i] << "  ";
			for (int j=0; j < fNumElementNodes_PMultiplier * fNumDOF_PMultiplier; j++)
			{
				cout << fK_hh(j)[i] << "+";
				if (fK_ct(j)[i] != 0)
					cout << "cst" << "  ";
				else
					cout << "0.0" << "  ";
			}
			cout << endl;
		}
	}
	/*******************/
}

/* calculate the internal force contribution ("-k*d") */
void GradSmallStrainT::FormKd(double constK)
{
	const char caller[] = "GradSmallStrainT::FormKd";
	int curr_group = ElementSupport().CurrentGroup();

	/* partition residual force vector */
	int neq_Disp = fNumElementNodes_Disp*fNumDOF_Disp;
	int neq_PMultiplier = fNumElementNodes_PMultiplier*fNumDOF_PMultiplier;

	dArrayT RHS_Disp(neq_Disp, fRHS.Pointer());
	dArrayT RHS_PMultiplier(neq_PMultiplier, fRHS.Pointer(neq_Disp));

	/* integration rules (same for both fields)*/
	const double* Det = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	/********DEBUG*******/
	if (fprint_Kd)
	{
		cout << caller << endl;
		cout<<"	 element: "<<CurrElementNumber()<<endl;
	}
	/*******************/

	/* initialize */
	fShapes->TopIP();

	/* integrate residuals over the elements */
	while(fShapes->NextIP())
	{
		/* integration weight */
		double scale = constK*(*Det++)*(*Weight++);

		/* accumulate strain displacement matrix */
		if (fStrainDispOpt == kMeanDilBbar)
			Set_B_bar(fShapes->Derivatives_U(), fMeanGradient, fB);
		else
			Set_B(fShapes->Derivatives_U(), fB);

		/* compute internal stresses residual */
		fB.MultTx(fCurrMaterial_Grad->s_ij(), fNEEvec);

		/* accumulate in rhs dislacement equations */
		RHS_Disp.AddScaled(scale, fNEEvec);
		
		/* accumulate shape functions */
		Set_h(fh);
		Set_p(fp);
		Set_q(fq);

		/* obtain yield condition residual */
		double yield = fCurrMaterial_Grad->yc();

		/* accumulate in rhs pmultiplier equations */
		RHS_PMultiplier.AddScaled(-scale*yield, fhT.Transpose(fh));

		/* increment if the ip leaves the softening branch */
		fWeakened += fCurrMaterial_Grad->weakened();

		/********DEBUG*******/
		if (fprint_Kd)
		{
			cout<<"            ip: "             << CurrIP()                                                 <<endl;
			cout<<"                 strain         : " << sqrt(fStrain_List[CurrIP()].ScalarProduct())       <<endl;
			cout<<"                 stress         : " << sqrt(fCurrMaterial_Grad->s_ij().ScalarProduct())   <<endl;
			cout<<"                 yc             : " << fCurrMaterial_Grad->yc()                           <<endl;
			cout<<"                 PMultiplier    : " << fPMultiplier_List[CurrIP()]                        <<endl;
			cout<<"                 del_PMultiplier: " << fPMultiplier_List[CurrIP()][0] -  fPMultiplier_last_List[CurrIP()][0]                      <<endl;
			cout<<"                 LapPMultiplier : " << fLapPMultiplier_List[CurrIP()]                     <<endl;
			cout<<"                 odm_bh_ij[0]   : " << (fCurrMaterial_Grad->odm_hb_ij())[0]               <<endl;
			cout<<"                 odm_hb_ij[0]   : " << (fCurrMaterial_Grad->odm_bh_ij())[0]               <<endl;
			cout<<"                 gm_hh          : " << fCurrMaterial_Grad->gm_hh()                        <<endl;
			cout<<"                 gm_hp[0]       : " << fCurrMaterial_Grad->gm_hp()[0]                     <<endl;
			cout<<"                 gm_hq          : " << fCurrMaterial_Grad->gm_hq()                        <<endl;
		}
		/*******************/
		/********DEBUG*******/
		if (false)
		{
			cout<<"            ip: "             << CurrIP()                                                 <<endl;
			if (NumSD() == 2)
			{
			for (int sd = 0; sd < 3; sd ++)
				cout<<"                 strain["<<sd<<"]      : " << fStrain_List[CurrIP()][sd]                         <<endl;
			for (int sd = 0; sd < 3; sd ++)
			cout<<"                 stress["<<sd<<"]      : " << fCurrMaterial_Grad->s_ij()[sd]                     <<endl;
			}
			else if (NumSD() == 3)
			{
			for (int sd = 0; sd < 6; sd ++)
				cout<<"                 strain["<<sd<<"]      : " << fStrain_List[CurrIP()][sd]                         <<endl;
			for (int sd = 0; sd < 6; sd ++)
			cout<<"                 stress["<<sd<<"]      : " << fCurrMaterial_Grad->s_ij()[sd]                     <<endl;
			}
			
			cout<<"                 yc             : " << fCurrMaterial_Grad->yc()                           <<endl;
			cout<<"                 PMultiplier    : " << fPMultiplier_List[CurrIP()]                        <<endl;
			cout<<"                 LapPMultiplier : " << fLapPMultiplier_List[CurrIP()]                     <<endl;
			cout<<"                 odm_bh_ij[0]   : " << (fCurrMaterial_Grad->odm_hb_ij())[0]               <<endl;
			cout<<"                 del_PMultiplier: " << fPMultiplier_List[CurrIP()][0] -  fPMultiplier_last_List[CurrIP()][0]                      <<endl;
			cout<<"                 odm_hb_ij[0]   : " << (fCurrMaterial_Grad->odm_bh_ij())[0]               <<endl;
			cout<<"                 gm_hh          : " << fCurrMaterial_Grad->gm_hh()                        <<endl;
			cout<<"                 gm_hp[0]       : " << fCurrMaterial_Grad->gm_hp()[0]                     <<endl;
			cout<<"                 gm_hq          : " << fCurrMaterial_Grad->gm_hq()                        <<endl;
		}
		/*******************/
	}
	/********DEBUG*******/
	if (fprint_KdMatrix)
	{
		cout << caller << " Kd Matrix" << endl;
		cout << "  element: " << CurrElementNumber() << endl;
		for (int i=0; i < neq_Disp; i++)
			cout<<"                 RHS_Disp["<< i << "]: " << RHS_Disp[i] << endl;
		for (int i=0; i < neq_PMultiplier; i++)
			cout<<"                 RHS_PMultiplier["<< i << "]: " << RHS_PMultiplier[i] << endl;
		cout << endl;
	}
	/*******************/
}

/* current element operations */
bool GradSmallStrainT::NextElement(void)
{
	const char caller[] = "GradSmallStrainT::NextElement";

	/* inherited */
	bool result = SmallStrainT::NextElement();
	
	/* get material pointer */
	if (result)
	{
		ContinuumMaterialT* pcont_mat = (*fMaterialList)[CurrentElement().MaterialNumber()];
	
		/* cast is safe since class contructs materials list */
		fCurrMaterial_Grad = (GradSSSolidMatT*) pcont_mat;
	}
	
	return result;
}

/* form of tangent matrix */
GlobalT::SystemTypeT GradSmallStrainT::TangentType(void) const
{
	const char caller[] = "GradSmallStrainT::TangentType";

	return GlobalT::kNonSymmetric;
}

/***********************************************************************
* Private
***********************************************************************/

/* accumulate shape function matrix */
void GradSmallStrainT::Set_h(dMatrixT& h) const
{
	const char caller[] = "GradSmallStrainT::Set_h";

	const double* pShapes_PMultiplier = fShapes_PMultiplier->IPShapeU(fShapes->CurrIP());
	double* ph = h.Pointer();

	for (int i = 0; i < fNumElementNodes_PMultiplier*fNumDOF_PMultiplier; i++)
		*ph++ = *pShapes_PMultiplier++;
}

/* accumulate shape function matrix */
void GradSmallStrainT::Set_p(dMatrixT& p) const
{
	const char caller[] = "GradSmallStrainT::Set_p";

	double* pp = p.Pointer();

	/* 1D */
	if (fNumSD == 1)
	{
		const double* pNax = fShapes_PMultiplier->IPDShapeU(fShapes->CurrIP())(0);
		for (int i = 0; i < fNumElementNodes_PMultiplier*fNumDOF_PMultiplier; i++)
			*pp++ = *pNax++;
	}
	/* 2D */
	else if (fNumSD == 2)
	{
		const double* pNax = fShapes_PMultiplier->IPDShapeU(fShapes->CurrIP())(0);
		const double* pNay = fShapes_PMultiplier->IPDShapeU(fShapes->CurrIP())(1);
		for (int i = 0; i < fNumElementNodes_PMultiplier*fNumDOF_PMultiplier; i++)
		{
			*pp++ = *pNax++;
			*pp++ = *pNay++;
		}
	}
	/* 3D */
	else
	{
		const double* pNax = fShapes_PMultiplier->IPDShapeU(fShapes->CurrIP())(0);
		const double* pNay = fShapes_PMultiplier->IPDShapeU(fShapes->CurrIP())(1);
		const double* pNaz = fShapes_PMultiplier->IPDShapeU(fShapes->CurrIP())(2);
		for (int i = 0; i < fNumElementNodes_PMultiplier*fNumDOF_PMultiplier; i++)
		{
			*pp++ = *pNax++;
			*pp++ = *pNay++;
			*pp++ = *pNaz++;
		}
	}
}

/* accumulate Laplacian C1 shape function matrix */
void GradSmallStrainT::Set_q(dMatrixT& q) const
{
	const char caller[] = "GradSmallStrainT::Set_q";

	double* pq = q.Pointer();

	/* 1D */
	if (fNumSD == 1)
	{
		const double* pNaxx = fShapes_PMultiplier->IPDDShapeU(fShapes->CurrIP())(0);
		for (int i = 0; i < fNumElementNodes_PMultiplier*fNumDOF_PMultiplier; i++)
			*pq++ = *pNaxx++;
	}
	/* 2D */
	else if (fNumSD == 2)
	{
		const double* pNaxx = fShapes_PMultiplier->IPDDShapeU(fShapes->CurrIP())(0);
		const double* pNayy = fShapes_PMultiplier->IPDDShapeU(fShapes->CurrIP())(1);
		for (int i = 0; i < fNumElementNodes_PMultiplier*fNumDOF_PMultiplier; i++)
			*pq++ = *pNaxx++ + *pNayy++;
	}
	/* 3D */
	else
	{
		const double* pNaxx = fShapes_PMultiplier->IPDDShapeU(fShapes->CurrIP())(0);
		const double* pNayy = fShapes_PMultiplier->IPDDShapeU(fShapes->CurrIP())(1);
		const double* pNazz = fShapes_PMultiplier->IPDDShapeU(fShapes->CurrIP())(2);
		for (int i = 0; i < fNumElementNodes_PMultiplier*fNumDOF_PMultiplier; i++)
			*pq++ = *pNaxx++ + *pNayy++ + *pNazz++;
	}
}
