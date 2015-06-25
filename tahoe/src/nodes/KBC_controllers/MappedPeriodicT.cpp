/* $Id: MappedPeriodicT.cpp,v 1.10 2005/06/07 07:32:07 paklein Exp $ */
/* created: paklein (04/07/1997) */
#include "MappedPeriodicT.h"

#include "FEManagerT.h"
#include "BasicFieldT.h"
#include "ParameterUtils.h"
#include "ParameterContainerT.h"

using namespace Tahoe;

/* column indeces */
const int kMaster = 0;
const int kSlave  = 1;

/* constructor */
MappedPeriodicT::MappedPeriodicT(const BasicSupportT& support, BasicFieldT& field):
	KBC_ControllerT(support),
	fField(field),
	fSchedule(NULL),
	fDummySchedule(1.0) //still need this?
{
	SetName("mapped_nodes");
}

/* initial condition */
void MappedPeriodicT::InitialCondition(void)
{
	/* reference coordinates */
	const dArray2DT& init_coords = fSupport.InitialCoordinates();

	/* compute mapping: F - 1 */
	fF = fFperturb;

	int ndof = fF.Rows();
	dArrayT	X;
	dArrayT d(ndof);
	int mappedcount = fMappedNodeList.Length();
	int dex = 0;
	for (int i = 0; i < mappedcount; i++)
	{
		int node = fMappedNodeList[i];
	
		/* fetch pointers */
		init_coords.RowAlias(node, X);
	
		/* map current coord */
		fF.Multx(X, d);
		
		/* set prescribed displacements */
		for (int j = 0; j < ndof; j++)
		{
			fMappedCards[dex].SetValues(node, j, KBC_CardT::kDsp, fSchedule, d[j]);
			dex++;
		}
	}
}

/* set BC cards for current step */
void MappedPeriodicT::InitStep(void)
{
	/* inherited */
	KBC_ControllerT::InitStep();
	
	/* just compute mapping */
	if (fSlaveMasterPairs.MajorDim() == 0)
	{
		/* compute F - 1  = fFperturb */
		fF.SetToScaled(fSchedule->Value(), fFperturb);
	}
	else /* apply mapping */
	{
		/* nodal information */
		const dArray2DT& init_coords = fSupport.InitialCoordinates();
		const dArray2DT& disp = fField[0];	

		/* compute F - 1  = fFperturb */
		fF.SetToScaled(fSchedule->Value(), fFperturb);

		int dex = 0;
		int ndof = fF.Rows();
		int num_pairs = fSlaveMasterPairs.MajorDim();	
		dArrayT	X_m, d_m;      //master data
		dArrayT X_s, d_s(ndof); //slave data
		for (int i = 0; i < num_pairs; i++)
		{
			int node_m = fSlaveMasterPairs(i, kMaster);
			int node_s = fSlaveMasterPairs(i, kSlave);

			/* fetch pointers */
			init_coords.RowAlias(node_m, X_m);
			disp.RowAlias(node_m, d_m);		
			init_coords.RowAlias(node_s, X_s);

			/* periodic displacements */		
			fD_sm.DiffOf(X_s, X_m);
			fF.Multx(fD_sm, d_s);
		
			/* set cards */
			for (int j = 0; j < ndof; j++)
				fSlaveCards[dex++].SetValues(node_s, j, KBC_CardT::kDsp, NULL, d_s[j] + d_m[j]);
		}
	}
}

/* output current configuration */
void MappedPeriodicT::WriteOutput(ostream& out) const
{
	/* inherited */
	KBC_ControllerT::WriteOutput(out);

	/* write output */
	out << "\n Field mapping:\n";
	out << fF << endl;
}

/* describe the parameters needed by the interface */
void MappedPeriodicT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	KBC_ControllerT::DefineParameters(list);

	/* schedule */
	list.AddParameter(ParameterT::Integer, "schedule");
}

/* information about subordinate parameter lists */
void MappedPeriodicT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	KBC_ControllerT::DefineSubs(sub_list);

	/* perturbation */
	sub_list.AddSub("F_perturb_choice", ParameterListT::Once, true);
	
	/* list of mapped nodes */
	sub_list.AddSub("mapped_node_ID_list");

	/* leader-follower node pairs */
	sub_list.AddSub("leader_follower_node_ID_list", ParameterListT::ZeroOrOnce, true);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* MappedPeriodicT::NewSub(const StringT& name) const
{
	if (name == "F_perturb_choice") {

		ParameterContainerT* F_choice = new ParameterContainerT(name);
		F_choice->SetSubSource(this);
		
		/* by dimension */
		F_choice->SetListOrder(ParameterListT::Choice);
		F_choice->AddSub("Matrix_1x1");
		F_choice->AddSub("Matrix_2x2");
		F_choice->AddSub("Matrix_3x3");
		F_choice->AddSub("Matrix_Nx1");
		F_choice->AddSub("Matrix_Nx2");
		F_choice->AddSub("Matrix_Nx3");
	
		return F_choice;
	}
	else if (name == "leader_follower_node_ID_list") {

		ParameterContainerT* node_pairs = new ParameterContainerT(name);
		
		/* by dimension */
		node_pairs->SetListOrder(ParameterListT::Sequence);
		node_pairs->AddSub("leader_node_ID_list");
		node_pairs->AddSub("follower_node_ID_list");
	
		return node_pairs;	
	}
	else if (name.StringMatch("Matrix_Nx"))
	{
		ParameterContainerT* matrix = new ParameterContainerT(name);
	
		/* row vectors */
		StringT vec_name = "Vector_";
		vec_name.Append(name[name.StringLength() - 1]);
		matrix->AddSub(vec_name, ParameterListT::OnePlus);
	
		return matrix;
	}
	else /* inherited */
		return KBC_ControllerT::NewSub(name);
}

/* accept parameter list */
void MappedPeriodicT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "MappedPeriodicT::TakeParameterList"; 

	/* inherited */
	KBC_ControllerT::TakeParameterList(list);

	/* dimension workspace */
	int nsd = fSupport.NumSD();
	int ndof = fField.NumDOF();
	fFperturb.Dimension(ndof,nsd);
	fF.Dimension(ndof,nsd);
	fD_sm.Dimension(nsd);

	/* schedule */
	int schedule = list.GetParameter("schedule");
	fSchedule = fSupport.Schedule(--schedule);
	if (!fSchedule) ExceptionT::GeneralFail(caller, "could not resolve schedule %d", schedule+1);
	
	/* perturbation matrix */
	const char *mat_names[] = {"", "Matrix_1x1", "Matrix_2x2", "Matrix_3x3"};
	const ParameterListT* matrix_parameters = list.List(mat_names[nsd]);
	if (matrix_parameters) /* square matrix */
		MatrixParameterT::Extract(*matrix_parameters, fFperturb);
	else /* general matrix */
	{
		const char *mat_names[] = {"", "Matrix_Nx1", "Matrix_Nx2", "Matrix_Nx3"};
		const ParameterListT& matrix_parameters = list.GetList(mat_names[nsd]);
		StringT vec_name = "Vector_";
		vec_name.Append(nsd);
		int num_rows = matrix_parameters.NumLists(vec_name);
		if (num_rows != ndof)
			ExceptionT::GeneralFail(caller, "expecting %d not %d occurrences of \"%s\" in \"%s\"",
				ndof, num_rows, vec_name.Pointer(), mat_names[nsd]);

		/* extract rows */
		dArrayT row(nsd);
		for (int i = 0; i < num_rows; i++) {
			VectorParameterT::Extract(matrix_parameters.GetList(vec_name, i), row);
			fFperturb.SetRow(i, row);
		}
	}

	/* list of mapped nodes */
	ArrayT<StringT> id_list;
	StringListT::Extract(list.GetList("mapped_node_ID_list"), id_list);
	GetNodes(id_list, fMappedNodeList);
	
	/* read leader-follower node pairs */
	const ParameterListT* leader_nodes = list.List("leader_node_ID_list");
	const ParameterListT* follower_nodes = list.List("follower_node_ID_list");
	if (leader_nodes && follower_nodes) {
	
		/* leader nodes */
		StringListT::Extract(*leader_nodes, id_list);
		iArrayT tmp;
		GetNodes(id_list, tmp);
		fSlaveMasterPairs.Dimension(tmp.Length(), 2);
		fSlaveMasterPairs.SetColumn(kMaster, tmp);		

		/* follower nodes */
		StringListT::Extract(*follower_nodes, id_list);
		if (tmp.Length() != fSlaveMasterPairs.MajorDim())
			ExceptionT::GeneralFail(caller, "%d follower nodes does not match %d leader nodes",
				tmp.Length(), fSlaveMasterPairs.MajorDim());
		fSlaveMasterPairs.SetColumn(kSlave, tmp);
	}
	else if (leader_nodes || follower_nodes)
		ExceptionT::BadInputValue(caller, "\"leader_node_ID_list\" and \"follower_node_ID_list\" must be defined together");
	
	/* generate BC cards */
	int num_BC = fMappedNodeList.Length() + fSlaveMasterPairs.MajorDim();
	fKBC_Cards.Dimension(num_BC*ndof);
	fMappedCards.Set(fMappedNodeList.Length()*ndof, fKBC_Cards.Pointer());
	fSlaveCards.Set(fSlaveMasterPairs.MajorDim()*ndof,
		fKBC_Cards.Pointer(fMappedCards.Length()));

	/* mapped nodes */
	int dex = 0;
	for (int i = 0; i < fMappedNodeList.Length(); i++)
		for (int j = 0; j < ndof; j++)
			fMappedCards[dex++].SetValues(fMappedNodeList[i], j, KBC_CardT::kDsp, NULL, 0.0);	

	/* slave nodes */
	dex = 0;
	for (int ii = 0; ii < fSlaveMasterPairs.MajorDim(); ii++)
		for (int jj = 0; jj < ndof; jj++)
		{
			/* set values */
			fSlaveCards[dex].SetValues(fSlaveMasterPairs(ii, kSlave), jj, KBC_CardT::kDsp, &fDummySchedule, 0.0);
			dex++;
		}	
}
