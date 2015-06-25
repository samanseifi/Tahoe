/* $Id: SCNIMFT.cpp,v 1.65 2005/11/16 22:23:37 jzimmer Exp $ */
#include "SCNIMFT.h"

#include "ArrayT.h"
#include "ofstreamT.h"
#include "ifstreamT.h"
#include "fstreamT.h"
#include "eIntegratorT.h"
#include "OutputSetT.h"
#include "ElementSupportT.h"
#include "ModelManagerT.h"
#include "iGridManagerT.h"
#include "iNodeT.h"
#include "CommManagerT.h"
#include "CommunicatorT.h"
#include "BasicFieldT.h"
#include "ParameterUtils.h"
#include "ParameterContainerT.h"
#include "InverseMapT.h"
#include "dSymMatrixT.h"

#include "MeshFreeSupport2DT.h"
#include "MeshFreeSupport3DT.h"
#include "MeshFreeNodalShapeFunctionT.h"
#include "ContinuumMaterialT.h"
#include "SolidMaterialT.h"
#include "SolidMatSupportT.h"

/* cell geometries */
#include "CellGeometryT.h"
#include "VoronoiDiagramT.h"
#include "CellFromMeshT.h"

#include "Traction_CardT.h"

/* materials lists */
#include "SSSolidMatList1DT.h"
#include "SSSolidMatList2DT.h"
#include "SSSolidMatList3DT.h"

/* for manufactured solutions */
#include "IsotropicT.h"

using namespace Tahoe;

const int kNoTractionVector = -1;

const int kNumOutput = 6;
static const char* OutputNames[kNumOutput] = {
  "coordinates",
  "displacement",
  "mass",
  "strain",
  "stress",
  "material_output",	
};

/* constructors */
SCNIMFT::SCNIMFT(const ElementSupportT& support, const FieldT& field):
  ElementBaseT(support),
  fMaterialList(NULL),
  fNodalShapes(NULL),
  fCellGeometry(NULL),
  qIsAxisymmetric(false),
  fBodySchedule(NULL)
{
#pragma unused(field)

  SetName("mfparticle");

  /* set matrix format */
  fLHS.SetFormat(ElementMatrixT::kNonSymmetric);
}

SCNIMFT::SCNIMFT(const ElementSupportT& support):
  ElementBaseT(support),
  fMaterialList(NULL),
  fNodalShapes(NULL),
  fCellGeometry(NULL),
  qIsAxisymmetric(false),
  fBodySchedule(NULL)
{
  SetName("mfparticle");

  /* set matrix format */
  fLHS.SetFormat(ElementMatrixT::kNonSymmetric);
}

/* destructor */
SCNIMFT::~SCNIMFT(void)
{
  delete fCellGeometry;
  delete fMaterialList;
  delete fNodalShapes;
}

/* initialization */
void SCNIMFT::TakeParameterList(const ParameterListT& list)
{
  const char caller[] = "SCNIMFT::TakeParameterList";
	
  /* dimension */
  int nsd = NumSD();
	
  /* get parameters needed to construct shape functions */
  fMeshfreeParameters = list.ListChoice(*this, "meshfree_support_choice");
	
  /* access to the model database */
  ModelManagerT& model = ElementSupport().ModelManager();

  /* extract particle ID's */
  const ParameterListT& particle_ID_params = list.GetList("mf_particle_ID_list");
  ArrayT<StringT> particle_ID_list;
  StringListT::Extract(particle_ID_params, particle_ID_list);

  //get nodes from ModelManagerT
  model.ManyNodeSets(particle_ID_list, fNodes);
	
  /** This class and its derived classes assume the list of nodes is sorted */
  fNodes.SortAscending();
		
  /* set inverse map */
  fNodes_inv.SetOutOfRange(InverseMapT::MinusOne);
  fNodes_inv.SetMap(fNodes);
	
  // get the cell geometry parameters
  const ParameterListT* geometry_params = list.ListChoice(*this, "cell_geometry_choice");
  if (geometry_params->Name() == "voronoi_diagram")
    fCellGeometry = new VoronoiDiagramT(ElementSupport(), qIsAxisymmetric);
  else if (geometry_params->Name() == "cell_from_mesh")
    fCellGeometry = new CellFromMeshT(ElementSupport(), qIsAxisymmetric);
  else
    ExceptionT::GeneralFail(caller,"Cannot get valid cell geometry from input file\n");
  fCellGeometry->TakeParameterList(*geometry_params);
  fCellGeometry->SetNodalElements(this);
	
  /* inherited */
  ElementBaseT::TakeParameterList(list);

  /* re-dimension "element" force and stiffness contributions */
  fLHS.Dimension(nsd);
	
  /* allocate work space */
  fForce_man.SetWard(0, fForce, nsd);
  fForce_man.SetMajorDimension(ElementSupport().NumNodes(), false);

  /* write parameters */
  ostream& out = ElementSupport().Output();
	
  /* shape functions */
  /* only support single list of integration cells for now */
  if (fElementConnectivities.Length() > 1) {
    ExceptionT::GeneralFail(caller,"Multiple ElementConnectivities not yet supported\n");
  }

  /* construct shape functions */
  fNodalShapes = new MeshFreeNodalShapeFunctionT(nsd,
						 ElementSupport().InitialCoordinates(), *fElementConnectivities[0], 
						 fNodalCoordinates, *fMeshfreeParameters);
  if (!fNodalShapes) throw ExceptionT::kOutOfMemory;
	
  /* echo parameters */
  fNodalShapes->WriteParameters(ElementSupport().Output());

  /* MLS stuff */
  fNodalShapes->SetSupportSize();

  /* exchange nodal parameters (only Dmax for now) */
  const ArrayT<int>* p_nodes_in = ElementSupport().ExternalNodes();
  if (p_nodes_in) {

    /* skip MLS fit at external nodes */
    iArrayT nodes_in;
    nodes_in.Alias(*p_nodes_in);
    fNodalShapes->SetSkipNodes(nodes_in);
		
    /* exchange */
    CommManagerT& comm = ElementSupport().CommManager();

    /* send all */
    dArray2DT& nodal_params = fNodalShapes->NodalParameters();

    /* initialize the exchange */
    int id = comm.Init_AllGather(nodal_params);
		
    /* do the exchange */
    comm.AllGather(id, nodal_params);
		
    /* clear the communication */
    comm.Clear_AllGather(id);
  }
	
  /* set nodal neighborhoods */
  fNodalShapes->SetNeighborData();
	
  /* final MLS initializations */
  fNodalShapes->WriteStatistics(ElementSupport().Output());
	
  /* initialize workspace for strain smoothing */
  fCellGeometry->SetNodesAndShapes(&fNodes, &fNodalCoordinates, fNodalShapes);
  fCellGeometry->ComputeBMatrices(nodalCellSupports, bVectorArray, fCellVolumes, 
				  fCellCentroids, circumferential_B);

  /* store shape functions at nodes */
  int nNodes = fNodes.Length();
  dArrayT nodalCoords;
  const RaggedArray2DT<int>& nodeSupport = fNodalShapes->NodeNeighbors();
  ArrayT<int> neighbor_list;
  ArrayT< LinkedListT<double> > nodal_phi;
  ArrayT< LinkedListT<int> > nodal_supports;
  nodal_phi.Dimension(nNodes);
  nodal_supports.Dimension(nNodes);
  for (int i = 0; i < nNodes; i++) {
    int node_i = fNodes[i];
    nodalCoords.Alias(nsd, fNodalCoordinates(i));
	
    neighbor_list.Dimension(nodeSupport.MinorDim(node_i));
    neighbor_list.Copy(nodeSupport(node_i));
    if (!fNodalShapes->SetFieldUsing(nodalCoords, neighbor_list)) // shift = 0 or not ?
      ExceptionT::GeneralFail("SCNIMFT::TakeParameterList","Shape Function evaluation"
			      "failed at node %d\n",fNodes[i]);		
		
    nodal_phi[i].AppendArray(fNodalShapes->FieldAt().Length(),
			     const_cast <double *> (fNodalShapes->FieldAt().Pointer()));
    nodal_supports[i].AppendArray(fNodalShapes->Neighbors().Length(),
				  const_cast <int *> (fNodalShapes->Neighbors().Pointer()));
  }
	
  // move into RaggedArray2DT's
  // move into more efficient storage for computation
  fNodalPhi.Configure(nodal_phi);
  fNodalSupports.Configure(nodal_supports);
	
  if (nodal_supports.Length() != nodal_phi.Length())
    ExceptionT::GeneralFail(caller,"nodal support indices and shape function values do not match\n");
		
  for (int i = 0; i < nodal_supports.Length(); i++) {
    int* irow_i = fNodalSupports(i);
    double* drow_i = fNodalPhi(i);
    LinkedListT<int>& ilist = nodal_supports[i];
    LinkedListT<double>& dlist = nodal_phi[i];
    ilist.Top(); dlist.Top();
    while (ilist.Next() && dlist.Next()) {
      *irow_i++ = *(ilist.CurrentValue());
      *drow_i++ = *(dlist.CurrentValue());
    }
  }

  /* output nodal shape function information */
  if (ElementSupport().Logging() == GlobalT::kVerbose)
    {
      /* output file root */
      StringT root;
      root.Root(ElementSupport().InputFile());
      ofstreamT out;

      /* nodal neighbors */
      StringT neighbor_file = root;
      neighbor_file.Append(".", Name(), ".nodal_neighbors");
      out.open(neighbor_file);
      fNodalShapes->MeshFreeSupport().WriteNodalNeighbors(out);
      out.close();

      /* nodal shape functions */
      StringT shape_file = root;
      shape_file.Append(".", Name(), ".nodal_phi");
      out.open(shape_file);
      fNodalShapes->MeshFreeSupport().WriteNodalShapes(out);
      out.close();
    }

  // store shape function information for boundary integration
  fCellGeometry->BoundaryShapeFunctions(fBoundaryPhi, fBoundarySupports, fBoundaryFacetNormals);
	
  /* material Data */
  ParameterListT mat_params;
  CollectMaterialInfo(list, mat_params);
  fMaterialList = NewMaterialList(mat_params.Name(), mat_params.NumLists());
 
  if (!fMaterialList)
    ExceptionT::GeneralFail(caller,"could not construct material list \"%s\"", mat_params.Name().Pointer());
  fMaterialList->TakeParameterList(mat_params);
	
  /* body force */
  const ParameterListT* body_force = list.List("body_force");
  if (body_force) {
    int schedule = body_force->GetParameter("schedule");
    fBodySchedule = ElementSupport().Schedule(--schedule);
    /* body force vector */
    const ArrayT<ParameterListT>& body_force_vector = body_force->Lists();
    if (body_force_vector.Length() != NumDOF())
      ExceptionT::BadInputValue(caller, "body force is length %d not %d",
				body_force_vector.Length(), NumDOF());
    fBody.Dimension(NumDOF());
    for (int i = 0; i < fBody.Length(); i++)
      fBody[i] = body_force_vector[i].GetParameter("value");
  }
	
  /* extract natural boundary conditions */
  TakeNaturalBC(list);

  /* output variables */
  fOutputFlags.Dimension(kNumOutput);
  fOutputFlags = 0;
  const ParameterListT* output = list.List("scni_output");
  if (output) 
    {
      /* set flags */
      for (int i = 0; i < kNumOutput; i++)
	{
	  /* look for entry */
	  const ParameterT* value = output->Parameter(OutputNames[i]);
	  if (value) {
	    int do_write = *value;
	    if (do_write)
	      fOutputFlags[i] = 1;
	  }
	}
    }
}

/* extract natural boundary condition information */
void SCNIMFT::TakeNaturalBC(const ParameterListT& list)
{
  const char caller[] = "SCNIMFT::TakeNaturalBC";

  int num_natural_bc = list.NumLists("natural_bc");
  int nsd = NumSD();
	
  // allocate data structures
  fTractionVectors.Dimension(num_natural_bc, nsd == 2 ? 3 : 6);
  fTractionVectors = 0.;
	
  // 
  fTractionBoundaryCondition.Dimension(fBoundaryFacetNormals.MajorDim());
  fTractionBoundaryCondition = kNoTractionVector;
	
  if (num_natural_bc > 0)
    {
      /* TEMP - turn on traction boundary condition for all boundary nodes */
      fTractionBoundaryCondition = 0;
		
      /* model manager */
      ModelManagerT& model = ElementSupport().ModelManager();
	
      /* temp space */
      ArrayT<StringT> block_ID(num_natural_bc);
      ArrayT<iArray2DT> localsides(num_natural_bc);
      iArrayT LTf(num_natural_bc);
      ArrayT<Traction_CardT::CoordSystemT> coord_sys(num_natural_bc);
      ArrayT<dArray2DT> values(num_natural_bc);
	    
      /* loop over natural BC's */
      int tot_num_sides = 0;
      for (int i = 0; i < num_natural_bc; i++) 
	{
	  const ParameterListT& natural_bc = list.GetList("natural_bc", i);
	    
	  /* side set */
	  const StringT& ss_ID = natural_bc.GetParameter("side_set_ID");
	  localsides[i] = model.SideSet(ss_ID);
	  int num_sides = localsides[i].MajorDim();
	  tot_num_sides += num_sides;
	  if (num_sides > 0)
	    {
	      block_ID[i] = model.SideSetGroupID(ss_ID);
	      LTf[i] = natural_bc.GetParameter("schedule");
	      coord_sys[i] = Traction_CardT::int2CoordSystemT(natural_bc.GetParameter("coordinate_system"));
	      int num_traction_vectors = natural_bc.NumLists("DoubleList");

	      /* switch to elements numbering within the group */
	      iArray2DT& side_set = localsides[i];
	      iArrayT elems(num_sides);
				
	      /* constant over the face */
	      if (num_traction_vectors == 1) {
		const ParameterListT& traction_vector = natural_bc.GetList("DoubleList");
		int dim = traction_vector.NumLists("Double");
						
		int minor_dim = nsd == 2 ? 3 : 6;
		if (dim != minor_dim)
		  ExceptionT::GeneralFail(caller, "expecting traction vector length %d not %d",
					  NumDOF(), dim);
		dArrayT t(minor_dim);
		for (int j = 0; j < minor_dim; j++)
		  t[j] = traction_vector.GetList("Double", j).GetParameter("value");	
							
		fTractionVectors.SetRow(0, t);
	      }
	      else
		{
		  ExceptionT::GeneralFail(caller,"Only constant traction over sideset implemented\n");
					
		}
	    }
	}
    }
}


/* form of tangent matrix */
GlobalT::SystemTypeT SCNIMFT::TangentType(void) const
{
  return GlobalT::kSymmetric;
}

void SCNIMFT::InitStep(void)
{
  /* inherited */
  ElementBaseT::InitStep();

  /* set material variables */
  if (fMaterialList)  fMaterialList->InitStep();
}

void SCNIMFT::CloseStep(void)
{
  /* inherited */
  ElementBaseT::CloseStep();

  if (fMaterialList) 
    {
      /* set material variables */
      fMaterialList->CloseStep();

      /* update element level internal variables */
      if (fMaterialList->HasHistoryMaterials())
	{
	  Top();
	  while (NextElement())
	    {
	      const ElementCardT& element = CurrentElement();
	      if (element.IsAllocated())
		{
		  ContinuumMaterialT* pmat = (*fMaterialList)[element.MaterialNumber()];

		  /* material update function */
		  pmat->UpdateHistory();
		}
	    }
	}
    }

  //---------------------
  // Compute error norms
  //---------------------

   // Create list of nodes
   const dArray2DT& initial_coordinates = fNodalCoordinates;
   int nNodes = fNodalCoordinates.MajorDim();
   ArrayT<int> nodeList(nNodes); 
   for (int i = 0; i < nNodes; i++) nodeList[i] = i;

   // Interpolate field
   dArray2DT fieldAtNodes(nNodes,NumSD());
   InterpolatedFieldAtNodes(nodeList, fieldAtNodes);

   // Loop over nodes
   double L2u = 0.0;
   double L2v = 0.0;
   double L1u = 0.0;
   double L1v = 0.0;
   double Linfu = 0.0;
   double Linfv = 0.0;
   double unorm = 0.0;
   double vnorm = 0.0;
   for (int i = 0; i < nNodes; i++) {

     double x = initial_coordinates(i,0);
     double y = initial_coordinates(i,1);
     double u = fieldAtNodes(i,0);
     double v = fieldAtNodes(i,1);

     // Linear field
//      double A1 =  0.010;
//      double B1 =  0.002;
//      double C1 =  0.000;
//      double A2 =  0.015;
//      double B2 = -0.008;
//      double C2 =  0.000;
//      double u_ex =  A1*x + B1*y + C1*x*y;
//      double v_ex =  A2*x + B2*y + C2*x*y;

     // Manufactured solution 1
     double A = 1.0;
     double B = 1.0;
     double C = 1.0;
     double D = 1.0;
     double Pi = acos(-1.0);
     double u_ex = sin(A*Pi*x)*cos(B*Pi*y);
     double v_ex = sin(C*Pi*x)*cos(D*Pi*y);


     double eu = fabs(u_ex - u);
     double ev = fabs(v_ex - v);
     unorm += u*u;
     vnorm += v*v;
     L2u += eu*eu;
     L2v += ev*ev;
     L1u += eu;
     L1v += ev;
     Linfu = (eu > Linfu) ? eu : Linfu; // max(eu,Linfu);
     Linfv = (ev > Linfv) ? ev : Linfv; // max(ev,Linfv);
   }
   unorm = sqrt(unorm/nNodes);
   vnorm = sqrt(vnorm/nNodes);
   L2u = sqrt(L2u/nNodes)/(unorm);
   L2v = sqrt(L2v/nNodes)/(vnorm);
   L1u = L1u/(nNodes*unorm);
   L1v = L1v/(nNodes*vnorm);
   Linfu = Linfu/unorm;
   Linfv = Linfv/vnorm;
/*   cout << endl;
   cout << "L2u   = " << L2u << endl;
   cout << "L2v   = " << L2v << endl;
   cout << "L1u   = " << L1u << endl;
   cout << "L1v   = " << L1v << endl;
   cout << "Linfu = " << Linfu << endl;
   cout << "Linfv = " << Linfv << endl;
   cout << "unorm = " << unorm << endl;
   cout << "vnorm = " << vnorm << endl;
   cout << endl; */

}

GlobalT::RelaxCodeT SCNIMFT::ResetStep(void)
{
  /* inherited */
  GlobalT::RelaxCodeT relax = ElementBaseT::ResetStep();

  /* update material internal variables */
  if (fMaterialList && fMaterialList->HasHistoryMaterials())
    {
      Top();
      while (NextElement())
	{
	  const ElementCardT& element = CurrentElement();		
	  if (element.IsAllocated())
	    {
	      ContinuumMaterialT* pmat = (*fMaterialList)[element.MaterialNumber()];

	      /* material reset function */
	      pmat->ResetHistory();
	    }
	}
    }

  return relax;
}

/* NOT implemented. Returns an zero force vector */
void SCNIMFT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
#pragma unused(field)
#pragma unused(node)
#pragma unused(force)
}

/* writing output */
void SCNIMFT::RegisterOutput(void)
{
  /* "point connectivities" needed for output */
  fPointConnectivities.Alias(fNodes.Length(), 1, fNodes.Pointer());

  /* block ID's */
  ArrayT<StringT> block_ID(fBlockData.Length());
  for (int i = 0; i < block_ID.Length(); i++)
    block_ID[i] = fBlockData[i].ID();

  /* get output labels (per node) */
  ArrayT<StringT> n_labels, e_labels;
  GenerateOutputLabels(n_labels);

  /* set output specifier */
  StringT set_ID;
  set_ID.Append(ElementSupport().ElementGroupNumber(this) + 1);
  OutputSetT output_set(GeometryT::kPoint, fPointConnectivities, n_labels, ChangingGeometry());
		
  /* register and get output ID */
  fOutputID = ElementSupport().RegisterOutput(output_set);
}

/* generate labels for output data */
void SCNIMFT::GenerateOutputLabels(ArrayT<StringT>& labels)
{
  const char caller[] = "SCNIMFT::GenerateOutputLabels";
  int ndof = NumDOF();
  if (ndof > 3) ExceptionT::GeneralFail(caller);

  /* number of output variables */
  iArrayT counts;
  SetOutputCount(fOutputFlags, counts);
  int num_output = counts.Sum();

  /* offsets to the different output values */
  iArrayT offsets(fOutputFlags.Length());
  offsets = 0;
  for (int i = 1; i < offsets.Length(); i++)
    offsets[i] = offsets[i-1] + counts[i-1];

  /* initialize */
  labels.Dimension(num_output);

  /* coordinates */
  if (fOutputFlags[kCoordinates]) {
    const char* ref[3] = {"X", "Y", "Z"};
    int index = offsets[kCoordinates];
    for (int i = 0; i < ndof; i++)
      labels[index++] = ref[i];
  }

  /* displacements */
  if (fOutputFlags[kDisplacement]) {

    /* labels from the field */
    const ArrayT<StringT>& field_labels = Field().Labels();

    int index = offsets[kDisplacement];
    for (int i = 0; i < ndof; i++)
      labels[index++] = field_labels[i];
  }

  /* mass */
  if (fOutputFlags[kMass])
    labels[offsets[kMass]] = "mass";

  /* strain */
  if (fOutputFlags[kStrain]) {
    const char* e1D[1] = {"e11"};
    const char* e2D[3] = {"e11", "e22", "e12"};
    const char* e3D[6] = {"e11", "e22", "e33", "e23", "e13", "e12"};

    const char** elabels = NULL;
    if (ndof == 1) elabels = e1D;
    else if (ndof == 2) elabels = e2D;
    else if (ndof == 3) elabels = e3D;
    else ExceptionT::GeneralFail(caller);	
		
    int nstrs = dSymMatrixT::NumValues(ndof);
    int index = offsets[kStrain];
    for (int i = 0; i < nstrs; i++)
      labels[index++] = elabels[i];
  }

  /* stress */
  if (fOutputFlags[kStress]) {
    const char* s1D[1] = {"s11"};
    const char* s2D[3] = {"s11", "s22", "s12"};
    const char* s3D[6] = {"s11", "s22", "s33", "s23", "s13", "s12"};

    const char** slabels = NULL;
    if (ndof == 1) slabels = s1D;
    else if (ndof == 2) slabels = s2D;
    else if (ndof == 3) slabels = s3D;
    else ExceptionT::GeneralFail(caller);	
		
    int nstrs = dSymMatrixT::NumValues(ndof);
    int index = offsets[kStress];
    for (int i = 0; i < nstrs; i++)
      labels[index++] = slabels[i];
  }

  /* material output labels */
  if (fOutputFlags[kMaterialOutput])
    {
      ArrayT<StringT> mat_labels;
      (*fMaterialList)[0]->OutputLabels(mat_labels);	
		
      int index = offsets[kMaterialOutput];
      for (int i = 0; i < mat_labels.Length(); i++)
	labels[index++] = mat_labels[i];
    }
}

/* compute specified output parameter(s) */
void SCNIMFT::SendOutput(int kincode)
{
#pragma unused(kincode)
  //TEMP: for now, do nothing
}

/* trigger reconfiguration */
GlobalT::RelaxCodeT SCNIMFT::RelaxSystem(void)
{
  /* inherited */
  GlobalT::RelaxCodeT relax = ElementBaseT::RelaxSystem();

  return relax;
}

/* construct field */
void SCNIMFT::NodalDOFs(const iArrayT& nodes, dArray2DT& DOFs) const
{
  /* get local numbers */
  iArrayT nodes_local = nodes;
  if (!GlobalToLocalNumbering(nodes_local))
    ExceptionT::GeneralFail("SCNIMFT::NodalDOFs", "map to local numbering failed");

  /* compute field at nodes */
  InterpolatedFieldAtNodes(nodes_local, DOFs);
}

/* write restart data to the output stream */
void SCNIMFT::WriteRestart(ostream& out) const
{
  /* inherited */
  ElementBaseT::WriteRestart(out);

  /* update element level internal variables */
  if (fMaterialList && fMaterialList->HasHistoryMaterials())
    {
      for (int i = 0; i < fElementCards.Length(); i++)
	{
	  const ElementCardT& element = fElementCards[i];
	  out << element.IsAllocated() << '\n';
	  if (element.IsAllocated()) element.WriteRestart(out);
	}
    }
}

/* read restart data to the output stream */
void SCNIMFT::ReadRestart(istream& in)
{
  /* inherited */
  ElementBaseT::ReadRestart(in);

  /* update element level internal variables */
  if (fMaterialList && fMaterialList->HasHistoryMaterials())
    {
      for (int i = 0; i < fElementCards.Length(); i++)
	{
	  int isallocated;
	  in >> isallocated;
	  if (isallocated) fElementCards[i].ReadRestart(in);
	}
    }
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* return number of values for each output variable */
void SCNIMFT::SetOutputCount(const iArrayT& flags, iArrayT& counts) const
{
  /* dimension check */
  if (flags.Length() != kNumOutput)
    ExceptionT::SizeMismatch("SCNIMFT::SetOutputCount");
	
  /* initialize */
  counts.Dimension(flags.Length());
  counts = 0;

  /* set output flags */
  if (flags[kCoordinates]) counts[kCoordinates] = NumSD();
  if (flags[kDisplacement]) counts[kDisplacement] = NumDOF();
  if (flags[kMass]) counts[kMass] = 1;
  if (flags[kStrain]) counts[kStrain] = dSymMatrixT::NumValues(NumSD());
  if (flags[kStress]) counts[kStress] = dSymMatrixT::NumValues(NumSD());

  /* material output variables */
  if (flags[kMaterialOutput])
    counts[kMaterialOutput] = (*fMaterialList)[0]->NumOutputVariables();
}

/* return true if connectivities are changing */
bool SCNIMFT::ChangingGeometry(void) const
{
  return ElementSupport().CommManager().PartitionNodesChanging();
}

/* echo element connectivity data */
void SCNIMFT::DefineElements(const ArrayT<StringT>& block_ID, const ArrayT<int>& mat_index) 
{
  const char caller[] = "SCNIMFT::DefineElements";

  //TEMP
  if (block_ID.Length() > 1)
    ExceptionT::GeneralFail(caller, "mutliple block ID's not supported %d",
			    block_ID.Length());
	
  /* access to the model database */
  ModelManagerT& model = ElementSupport().ModelManager();

  fElementConnectivities.Dimension(1);

  // NB THIS IS SPECIALIZED TO ONLY ONE ELEMENT BLOCK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  model.ReadConnectivity(block_ID[0]);

  /* set pointer to connectivity list */
  fElementConnectivities[0] = model.ElementGroupPointer(block_ID[0]);

  // Get nodal coordinates 
  int nsd = NumSD();
  fNodalCoordinates.Dimension(fNodes.Length(), nsd);
  fNodalCoordinates.RowCollect(fNodes, model.Coordinates());

  /* set up element cards for state variable storage */
  fElementCards.Dimension(fNodes.Length()); /* one card per node */
  for (int i = 0; i < fElementCards.Length(); i++)
    fElementCards[i].SetMaterialNumber(mat_index[0]);
	
  fCellGeometry->DefineElements(block_ID, mat_index);
}

/* collecting element group equation numbers */
void SCNIMFT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
			AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
#pragma unused(eq_1)

  /* dimension equations array */
  fEqnos.Configure(nodalCellSupports, NumDOF());

  /* get local equations numbers */
  Field().SetLocalEqnos(nodalCellSupports, fEqnos);

  /* add to list of equation numbers */
  eq_2.Append(&fEqnos);
}

/* assemble particle mass matrix into LHS of global equation system */
void SCNIMFT::AssembleParticleMass(const double rho)
{
  int nsd = NumSD();
  fForce = 0.0;
  int* nodes = fNodes.Pointer();
  double* volume = fCellVolumes.Pointer();
  for (int i = 0; i < fNodes.Length(); i++) {
    double* m = fForce(*nodes++);
    for (int j = 0; j < nsd; j++)
      *m++ = *volume*rho;
    volume++;
  }
	
  /* assemble all */
  ElementSupport().AssembleLHS(Group(), fForce, Field().Equations());
}

/* contribution from natural BCs */
void SCNIMFT::RHSDriver(void) 
{
  const char caller[] = "SCNIMFT::RHSDriver";

  int nsd = NumSD();
  fForce = 0.0;

  if (fTractionVectors.MajorDim()) {
    int nsd = NumSD();
    int numBoundaryFacets = fBoundaryFacetNormals.MajorDim();
    dArrayT traction_vector(nsd), normal_vector;
    dSymMatrixT workspace(nsd); 
    for (int i = 0; i < numBoundaryFacets; i++) 
      if (fTractionBoundaryCondition[i] != kNoTractionVector) {
	normal_vector.Alias(nsd, fBoundaryFacetNormals(i));
	fTractionVectors.RowCopy(fTractionBoundaryCondition[i], workspace.Pointer());
	//compute sigma dot n
	workspace.Multx(normal_vector, traction_vector);
				
	int* supp_i = fBoundarySupports(i) ;
	double* phi_i = fBoundaryPhi(i);
	int n_supp = fBoundaryPhi.MinorDim(i);
	for (int j = 0; j < n_supp; j++) {
	  double* fint = fForce(*supp_i++);
	  for (int k = 0; k < nsd; k++) 
	    *fint++ += traction_vector[k]*(*phi_i); 
	  phi_i++;
	}
      }

    // fForce gets multiplied by constKd?
  }
	
  double constMa = 0.0;
  int formMa = fIntegrator->FormMa(constMa); 
  ContinuumMaterialT *mat;
  SolidMaterialT* solid_material = NULL;
  IsotropicT* iso_material = NULL;
  int nnd;
  if (fBodySchedule || formMa)
    {
      /* just one material for now */
      mat = (*fMaterialList)[0];
      solid_material = TB_DYNAMIC_CAST(SolidMaterialT*, mat);
      if (!solid_material) ExceptionT::GeneralFail(caller, "cannot get material");
      nnd = fNodes.Length();

      /* cast to isotropic class */
      iso_material = TB_DYNAMIC_CAST(IsotropicT*, mat);
    }
	
  /* contribution from body force source */
  bool body_force_source = false;
  if (fBodySchedule || body_force_source) {

    /* nodal coordinates */
    const dArray2DT& initial_coordinates = ElementSupport().InitialCoordinates();

    /* work space */	
    double load_factor = fBodySchedule->Value();
    dArrayT bf_source(nsd);
    bf_source = 0.0;

    /* compute body force */
    double *f, *phi_i;
    int* supp_i, n_supp;
    int* nodes = fNodes.Pointer();
    double* volume = fCellVolumes.Pointer();
    double density = solid_material->Density();
    double twoPi = 2.0*acos(-1.0);
    for (int i = 0; i < nnd; i++) {
		  
      bf_source = 0.0;
		
      /* compute body force at node i */
      if (fBodySchedule) bf_source = fBody;
					
      //add additional contribution to the body force
      // Manufactured solution 1
      double x = initial_coordinates(i,0);
      double y = initial_coordinates(i,1);
      double A = 1.0;
      double B = 1.0;
      double C = 1.0;
      double D = 1.0;
      double Pi = acos(-1.0);
      double E = 1.0;
      double nu = 0.3;
      double lambda = nu*E/((1+nu)*(1-2*nu));
      double mu = E/(2*(1+nu));
      //bf_source[0] = Pi*Pi * ( (A*A*(lambda+2*mu)+B*B*mu)*sin(A*Pi*x)*cos(B*Pi*y)
 		//	       + C*D*(lambda+mu)*cos(C*Pi*x)*sin(D*Pi*y));
      //bf_source[1] = Pi*Pi * ( (C*C*mu + D*D*(lambda+2*mu))*sin(C*Pi*x)*cos(D*Pi*y)
		//	       + A*B*(lambda+mu)*cos(A*Pi*x)*sin(B*Pi*y));
			
      bf_source *= load_factor*density*(*volume++);
      if (qIsAxisymmetric) 
	bf_source *= twoPi*fCellCentroids(i,0);
		
      //			supp_i = nodalCellSupports(i);
      //			n_supp = nodalCellSupports.MinorDim(i);
      supp_i = fNodalSupports(i);
      n_supp = fNodalSupports.MinorDim(i);
      phi_i = fNodalPhi(i);
      for (int j = 0; j < n_supp; j++) {
	f = fForce(*supp_i++);
	for (int k = 0; k < nsd; k++)
	  // 			    *f++ -= bf_source[k]*(*phi_i); // This created - disp. for + body force
	  *f++ += bf_source[k]*(*phi_i); // This creates + disp. for + body force
	phi_i++;
      }
    }	
  }
	
  if (fIntegrator->FormMa(constMa)) {
    if (Field().Order() < 2)
      ExceptionT::GeneralFail(caller,"Field's Order does not have accelerations\n");
	
    const dArray2DT& a = Field()(0,2); 
    double* ma = fForce.Pointer();
    const double* acc;
    int* nodes = fNodes.Pointer();
    double* volume = fCellVolumes.Pointer();
    double density = solid_material->Density();
    if (!qIsAxisymmetric) {
      for (int i = 0; i < nnd; i++)
	{
	  acc = a(*nodes++);
	  for (int j = 0; j < nsd; j++)
	    *ma++ += density*(*volume)*(*acc++);
	  volume++;
	}
    } else {
      double twoPi = 2.0*acos(-1.0);
      double w_i;
      for (int i = 0; i < nnd; i++)
	{
	  acc = a(*nodes++);
	  w_i = density*(*volume)*twoPi*fCellCentroids(i,0);
	  for (int j = 0; j < nsd; j++)
	    *ma++ += w_i*(*acc++);
	  volume++;
	}
    }
				
  }
	
  if (fTractionVectors.MajorDim() || fBodySchedule || formMa) {
    /* Assemble */
    ElementSupport().AssembleRHS(Group(),fForce,Field().Equations());	
  }
	
}

int SCNIMFT::GlobalToLocalNumbering(ArrayT<int>& nodes) const
{
  int OK = 1;
  for (int i = 0; i < nodes.Length(); i++) {
    int local = fNodes_inv.Map(nodes[i]);
    if (local == -1) OK = 0;
    nodes[i] = local;
  }
  return OK;
}

int SCNIMFT::GlobalToLocalNumbering(RaggedArray2DT<int>& nodes)
{
  iArrayT row_i;
  for (int i = 0; i < nodes.MajorDim(); i++) {
    row_i.Set(nodes.MinorDim(i), nodes(i)); 
    if (!GlobalToLocalNumbering(row_i))
      return 0;
  }
	
  return 1;
}

void SCNIMFT::InterpolatedFieldAtNodes(const ArrayT<int>& nodes, dArray2DT& fieldAtNodes) const
{
  /* displacements */
  int nsd = NumSD();
  const dArray2DT& u = Field()(0,0);
  dArrayT vec, values_i;
  for (int i = 0; i < nodes.Length(); i++) {
    /* copy in */
    vec.Alias(nsd, fieldAtNodes.Pointer(i*nsd));
    vec = 0.;	
		
    int node_i = nodes[i];
    const int* nodal_supp = fNodalSupports(node_i);
    const double* phi_i = fNodalPhi(node_i);
    for (int j = 0; j < fNodalPhi.MinorDim(node_i); j++)
      vec.AddScaled(*phi_i++, u(*nodal_supp++));
  }
}

/** localNodes are local Numbers, so GlobalToLocalNumbering needs to have been called in whatever class 
 * calls this function. The node numbers returned in support are global. 
 */
void SCNIMFT::NodalSupportAndPhi(const iArrayT& localNodes, RaggedArray2DT<int>& support, 
				 RaggedArray2DT<double>& phi) const
{
  int nlnd = localNodes.Length();
  iArrayT minorDims(localNodes.Length());
  for (int i = 0; i < nlnd; i++) 
    minorDims[i] = fNodalPhi.MinorDim(localNodes[i]);
		
  support.Configure(minorDims);
  phi.Configure(minorDims);
	
  const int *lndi = localNodes.Pointer();
  for (int i = 0; i < nlnd; i++) {
    support.SetRow(i, fNodalSupports(*lndi));
    phi.SetRow(i, fNodalPhi(*lndi++));
  }
}

int SCNIMFT::SupportSize(int localNode) const {
  return fNodalPhi.MinorDim(localNode);
}

// XML stuff below

/* information about subordinate parameter lists */
void SCNIMFT::DefineSubs(SubListT& sub_list) const
{
  /* inherited */
  ElementBaseT::DefineSubs(sub_list);
	
  /* parameters for the meshfree support */
  sub_list.AddSub("cell_geometry_choice", ParameterListT::Once, true);
	
  /* parameters for the meshfree support */
  sub_list.AddSub("meshfree_support_choice", ParameterListT::Once, true);

  /* list of node set ID's defining which nodes get integrated */
  sub_list.AddSub("mf_particle_ID_list", ParameterListT::OnePlus);
	
  /* optional body force */
  sub_list.AddSub("body_force", ParameterListT::ZeroOrOnce);
	
  /* tractions */
  sub_list.AddSub("natural_bc", ParameterListT::Any);

  /* output variables */
  sub_list.AddSub("scni_output", ParameterListT::ZeroOrOnce);
}

/* return the description of the given inline subordinate parameter list */
void SCNIMFT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
			      SubListT& sub_lists) const
{
  /* inherited */
  ElementBaseT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SCNIMFT::NewSub(const StringT& name) const
{
  SCNIMFT* non_const_this = const_cast<SCNIMFT*>(this);
	
  /* try material list */
  MaterialListT* material_list = non_const_this->NewMaterialList(name, 0);

  if (material_list)
    return material_list;
  else if (name == "meshfree_support_choice") {
    ParameterContainerT* mf_choice = new ParameterContainerT(name);
    mf_choice->SetSubSource(this);
    mf_choice->SetListOrder(ParameterListT::Choice);
	  
    mf_choice->AddSub("meshfree_support_2D");
    mf_choice->AddSub("meshfree_support_3D");
	  
    return mf_choice;
  } else if (name == "meshfree_support_2D")
    return new MeshFreeSupport2DT;	
  else if (name == "meshfree_support_3D")
    return new MeshFreeSupport3DT;
  else if (name == "cell_geometry_choice") {
    ParameterContainerT* cg_choice = new ParameterContainerT(name);
    cg_choice->SetSubSource(this);
    cg_choice->SetListOrder(ParameterListT::Choice);
	  
    cg_choice->AddSub("voronoi_diagram");
    cg_choice->AddSub("cell_from_mesh");
	  
    return cg_choice;
  } else if (name == "voronoi_diagram")
    return new VoronoiDiagramT;
  else if (name == "cell_from_mesh")
    return new CellFromMeshT;
  else if (name == "body_force") { // body force
    ParameterContainerT* body_force = new ParameterContainerT(name);
	
    /* schedule number */
    body_force->AddParameter(ParameterT::Integer, "schedule");
	
    /* body force vector */
    body_force->AddSub("Double", ParameterListT::OnePlus); 		
		
    return body_force;
  } else if (name == "natural_bc") { /* traction bc */
    ParameterContainerT* natural_bc = new ParameterContainerT(name);

    natural_bc->AddParameter(ParameterT::Word, "side_set_ID");
    natural_bc->AddParameter(ParameterT::Integer, "schedule");

    ParameterT coord_sys(ParameterT::Enumeration, "coordinate_system");
    coord_sys.AddEnumeration("global", Traction_CardT::kCartesian);
    coord_sys.AddEnumeration( "local", Traction_CardT::kLocal);
    coord_sys.SetDefault(Traction_CardT::kCartesian);
    natural_bc->AddParameter(coord_sys);

    natural_bc->AddSub("DoubleList", ParameterListT::OnePlus); 		
		
    return natural_bc;
  } 
  else if (name == "scni_output") 
    {
      ParameterContainerT* output = new ParameterContainerT(name);
		
      /* all true by default */
      for (int i = 0; i < kNumOutput; i++) {
	ParameterT var(ParameterT::Integer, OutputNames[i]);
	var.SetDefault(1);
	output->AddParameter(var, ParameterListT::ZeroOrOnce);
      }

      return output;	
    }
  else /* inherited */
    return ElementBaseT::NewSub(name);
}

