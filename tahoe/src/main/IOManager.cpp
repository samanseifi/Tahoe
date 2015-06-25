/* $Id: IOManager.cpp,v 1.25 2007/01/15 05:55:33 paklein Exp $ */
/* created: sawimme (10/12/1999) */
#include "IOManager.h"

#include "ifstreamT.h"
#include "OutputSetT.h"
#include "dArrayT.h"
#include "OutputBaseT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"

using namespace Tahoe;

IOManager::IOManager(ostream& log, const StringT& program_name,
	const StringT& version, const StringT& title, const StringT& input_file,
	IOBaseT::FileTypeT output_format):
	fLog(log),
	fOutputFormat(output_format),
	fOutput(NULL),
	fEcho (false),
	fOutputTime(0.0),
	fOutput_tmp(NULL),
	fChangingFlag(kNoChangingFlag)
{
	/* construct output formatter */
	fOutput = IOBaseT::NewOutput(program_name, version, title, input_file, fOutputFormat, fLog);
}

IOManager::~IOManager(void)
{
	/* in case output is diverted */
	RestoreOutput();
	
	delete fOutput;
	fOutput = NULL;
}

/* how to override changing geometry flags */
void IOManager::SetChangingFlag(ChangingFlagT changing_flag)
{
	const char caller[] = "ModelManagerT::SetChangingFlag";

	/* can only override once */
	if (fChangingFlag != kNoChangingFlag) 
		ExceptionT::GeneralFail(caller, "changing flag has already been overridden");
	
	/* can only override before sets register */
	if (fOutput->ElementSets().Length() != 0)
		ExceptionT::GeneralFail(caller, "can only override flag before sets register");
	
	fChangingFlag = changing_flag;
}

void IOManager::EchoData (ostream& o) const
{
  IOBaseT temp (o);
  o << " Output format . . . . . . . . . . . . . . . . . = " << fOutputFormat  << '\n';
  temp.OutputFormats (o);
}

/*********** OUTPUT **************/


void IOManager::NextTimeSequence(int sequence_number)
{
	fOutput->NextTimeSequence(sequence_number);
}

/* set model coordinates */
void IOManager::SetCoordinates(const dArray2DT& coordinates, const iArrayT* node_id)
{
	fOutput->SetCoordinates(coordinates, node_id); 
}
	
/* register the output for an element set. returns the output ID */
int IOManager::AddElementSet(const OutputSetT& output_set)
{
	/* attempt to override changing geometry flag */
	bool changing = output_set.Changing();
	if (changing && fChangingFlag == kForceNotChanging)
		ExceptionT::GeneralFail("IOManager::AddElementSet", "cannot override changing flag");
	else if (fChangingFlag == kForceChanging)
	{
		/* not so const */
		OutputSetT* os = (OutputSetT*) &output_set;
		os->SetChanging(true);
	}

	return fOutput->AddElementSet(output_set);
}

const ArrayT<OutputSetT*>& IOManager::ElementSets(void) const
{
	return fOutput->ElementSets();
}

void IOManager::AddNodeSet(const iArrayT& nodeset, const StringT& setID)
{
	fOutput->AddNodeSet(nodeset, setID);
}

void IOManager::AddSideSet(const iArray2DT& sideset, const StringT& setID, const StringT& group_ID)
{
	fOutput->AddSideSet(sideset, setID, group_ID);
}

/* output functions */
void IOManager::WriteGeometry(void)
{
	fOutput->WriteGeometry();
}

void IOManager::WriteGeometryFile(const StringT& file_name,
	IOBaseT::FileTypeT format) const
{
	if (!fOutput)
		ExceptionT::GeneralFail("IOManager::WriteGeometryFile", 
			"output must be configured");

	fOutput->WriteGeometryFile(file_name, format);
}

void IOManager::WriteOutput(int ID, const dArray2DT& n_values, const dArray2DT& e_values) {
	if (fInsertNodalLabels.Length() == 0) {
		fOutput->WriteOutput(fOutputTime, ID, n_values, e_values);
	} else {
		/* initialize by copying in */
		dArray2DT node_values = n_values;

		/* get output set data */
		const iArrayT& nodes_used = fOutput->NodesUsed(ID);
		const OutputSetT& output_set = fOutput->OutputSet(ID);
		const ArrayT<StringT>& node_output_labels = output_set.NodeOutputLabels();

		/* fill output data */
		dArrayT vals(nodes_used.Length());
		dArrayT column(fInsertNodalData.MajorDim());
		for (int i = 0; i < fInsertNodalLabels.Length(); i++) {
			int col = -1;
			for (int j = 0; j < node_output_labels.Length() && col == -1; j++)
				if (fInsertNodalLabels[i] == node_output_labels[j])
					col = j;

			if (col != -1) {
				/* collect values */
				fInsertNodalData.ColumnCopy(col, column);
				vals.Collect(nodes_used, column);
			
				/* write into data */
				node_values.SetColumn(col, vals);
			}
		}

		/* send to output */
		fOutput->WriteOutput(fOutputTime, ID, node_values, e_values);	
	}
}

void IOManager::WriteOutput(int ID, const dArray2DT& n_values) {
	dArray2DT e_values;
	WriteOutput(ID, n_values, e_values);	
}

void IOManager::WriteOutput(const StringT& file, const dArray2DT& coords, const iArrayT& node_map,
	const dArray2DT& values, const ArrayT<StringT>& labels) const
{
	StringT junk = "N/A";
 	OutputBaseT* output = IOBaseT::NewOutput(junk, junk, junk, file, fOutputFormat, fLog);
	output->SetCoordinates(coords, &node_map);
	iArray2DT connectivities(node_map.Length(), 1);
	connectivities.SetValueToPosition();
	OutputSetT output_set(GeometryT::kPoint, connectivities, labels);
	int id = output->AddElementSet(output_set);
	dArray2DT e_values;
	output->WriteOutput(0.0, id, values, e_values);
	delete output;
}

/* return the list of nodes used by the output set */
const iArrayT& IOManager::NodesUsed(int ID) const { return fOutput->NodesUsed(ID); }

/* (temporarily) re-route output */
void IOManager::DivertOutput(const StringT& outfile)
{
	/* can only divert once */
	if (fOutput_tmp != NULL)
		cout << "\n IOManager::DivertOutput: cannot divert output to \""
		     << outfile << "\".\n"
		     <<   "     Output is already diverted to \"" <<
		     fOutput->OutputRoot() << "\"" << endl;
	else
	{
		/* store main out */
		fOutput_tmp = fOutput;
	
		/* construct temporary output formatter */
		StringT tmp(outfile);
		tmp.Append(".ext"); //OutputBaseT takes root of name passed in
		fOutput = IOBaseT::NewOutput(fOutput_tmp->CodeName(), fOutput_tmp->Version(),
			fOutput_tmp->Title(), tmp, fOutputFormat, fLog);
		
		/* add all output sets */
		const ArrayT<OutputSetT*>& element_sets = fOutput_tmp->ElementSets();
		for (int i = 0; i < element_sets.Length(); i++)
			AddElementSet(*(element_sets[i]));
			
		/* set coordinate data */
		fOutput->SetCoordinates(fOutput_tmp->Coordinates(), fOutput_tmp->NodeID());
	}
}

void IOManager::RestoreOutput(void)
{
	if (fOutput_tmp != NULL)
	{
		/* delete temp output formatter */
		delete fOutput;
	
		/* restore main out */
		fOutput = fOutput_tmp;
		fOutput_tmp = NULL;
	}
}

const OutputSetT& IOManager::OutputSet(int ID) const
{
	return fOutput->OutputSet(ID);
}

void IOManager::InsertNodalData(const ArrayT<StringT>& labels, const dArray2DT& data) {
	const char caller[] = "IOManager::InsertNodalData";
	if (!fOutput) ExceptionT::GeneralFail(caller, "output not set");
	if (labels.Length() != data.MinorDim())
		ExceptionT::SizeMismatch(caller, "number of labels does not match data dimension");
	if (data.MajorDim() != fOutput->Coordinates().MajorDim())
		ExceptionT::SizeMismatch(caller, "major dimension");
	
	 /* copy */
	fInsertNodalLabels = labels;
	fInsertNodalData = data;
}

void IOManager::ClearInsertNodalData(void) {
	fInsertNodalLabels.Dimension(0);
	fInsertNodalData.Dimension(0,0);
}

