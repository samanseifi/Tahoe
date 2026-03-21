/* $Id: SimoQ1P0.cpp,v 1.15 2026/03/01 23:09:12 samanseifi Exp $ */
#include "SimoQ1P0_Surface.h"

#include "FSDEMatQ1P02DT.h"
#include "FSDEMatSupportQ1P02DT.h"
#include "ShapeFunctionT.h"

/* From TLCBSurfaceT */
#include "ModelManagerT.h"
#include "MaterialListT.h"
#include "eIntegratorT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"
#include "InverseMapT.h"
#include "RowAutoFill2DT.h"
#include "OutputSetT.h"

#include "FEManagerT.h"
#include "TimeManagerT.h"

#include <iostream>
#include <algorithm>

using namespace std;

namespace Tahoe {

SimoQ1P0_Surface::SimoQ1P0_Surface(const ElementSupportT& support) :
    SimoQ1P0(support),
    fLocCurrCoords(LocalArrayT::kCurrCoords), fSurfaceCBSupport(NULL)
{
      SetName("updated_lagrangian_Q1P0_surface");
}

/* Destructor */
SimoQ1P0_Surface::~SimoQ1P0_Surface()
{
      /* TLCBSurface Stuff */
	/* free surface models */
	delete fSurfaceCBSupport;

}

// ------------------------------------------------specify parameters needed by the interface--------------------------------------------------------- //
void SimoQ1P0_Surface::DefineParameters(ParameterListT& list) const
{
      /* Inherited from the most similar element */
      SimoQ1P0::DefineParameters(list);

      /* Taking the output_surface */
      ParameterT output_surface(ParameterT::Boolean, "output_surface");
      output_surface.SetDefault(false);
      list.AddParameter(output_surface);
}

  /* ----------------------------- information about subordinate parameter lists ------------------------------------------------------------------------*/
void SimoQ1P0_Surface::DefineSubs(SubListT& sub_list) const
{

      /* inherited from the most similar element */
      SimoQ1P0::DefineSubs(sub_list);

      /* list of passivated surfaces - side set list */
      sub_list.AddSub("passivated_surface_ID_list", ParameterListT::ZeroOrOnce);

      /* one surface_tension block per sideset/value pair; zero or more allowed */
      sub_list.AddSub("surface_tension", ParameterListT::Any);

}

/* ----------------------------- construct sub-list parameter interface -----------------------------------------------------------------------*/
ParameterInterfaceT* SimoQ1P0_Surface::NewSub(const StringT& name) const
{
      if (name == "surface_tension") {
            ParameterContainerT* st = new ParameterContainerT(name);
            /* side set ID this surface tension is applied to */
            st->AddParameter(ParameterT::Word, "side_set_ID");
            /* surface tension coefficient gamma */
            st->AddParameter(ParameterT::Double, "gamma");
            /* ramp-up time; set to 0.0 for instant application */
            st->AddParameter(0.0, "t_0");
            return st;
      }
      return SimoQ1P0::NewSub(name);
}


  // ---------------------------------------------------- accept parameter list ------------------------------------------------------------------------ //
void SimoQ1P0_Surface::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "SimoQ1P0_Surface::TakeParameterList";

      /* inherited from the most similar element for bulk element */
      SimoQ1P0::TakeParameterList(list);

      /* Shape Functions and related values */
      const ShapeFunctionT& shape = ShapeFunction();
      int nsd = shape.NumSD();                          // # of spatial dimensions in problem
      int nfs = shape.NumFacets();                      // # of total possible surface facets?
      int nsi = shape.FacetShapeFunction(0).NumIP();    // # IPs per surface face (2 x 1 = 2 for 1D surface)
      int nfn = shape.FacetShapeFunction(0).NumNodes(); // # nodes on each surface face?
      int nen = NumElementNodes();                      // # nodes on each element
      int nme = nen * nsd;                              // # of mechanical degree of freedoms per element

      /* mechanical due to surface */
      fAmm_mat2.Dimension(nme, nme);

      /* support for the surface model */
      fF_Surf_List.Dimension(nsi);

      /* Need to actually place values into fF_Surf_List when testing (identity) */
      for (int i = 0; i < fF_Surf_List.Length(); i++)
            fF_Surf_List[i].Dimension(nsd);

      /* DUMMY INITIALIZE fF_Surf_List - SPECIFY DEFORMATION GRADIENT (copied from TotalLagrangianCBSurfaceT) */
      fF_Surf_List[0].Identity();

      /* Back to normal flow */
      /* Support for surface models */
      fSurfaceCBSupport = new FSMatSupportT(nsd, nsi);
      fSurfaceCBSupport->SetContinuumElement(this);
      fSurfaceCBSupport->SetDeformationGradient(&fF_Surf_List);



      /* Do we need to redefine this in "canonical" normal order? */
      /* 2D */
      double normals_dat[4*2] = {
            1.0, 0.0,
            -1.0, 0.0,
            0.0, 1.0,
            0.0,-1.0,
      };
      dArray2DT normals(4, 2, normals_dat);

      /* Allocating the normal vectors */
      fNormal.Dimension(nfs);
      for (int i = 0; i < nfs; i++)
      {
            fNormal[i].Dimension(nsd);
            fNormal[i] = normals(i);
      }

      /* surface tension is now specified per side set — see surface_tension sub-lists below */

      bool output_surface = list.GetParameter("output_surface");

      /* collect surface element information */
      ArrayT<StringT> block_ID;
      ElementBlockIDs(block_ID);
      ModelManagerT& model_manager = ElementSupport().ModelManager();
      model_manager.BoundingElements(block_ID, fSurfaceElements, fSurfaceElementNeighbors);

      ArrayT<const iArray2DT*> connects;
      model_manager.ElementGroupPointers(block_ID, connects);

      /* Pre-scan surface_tension sidesets: elements on internal interfaces
       * (multi-block) are not found by BoundingElements, so merge them in. */
      int n_st = list.NumLists("surface_tension");
      if (n_st > 0) {
            /* collect all sideset element/face pairs (group-numbered) */
            AutoArrayT<int> extra_elems;
            for (int k = 0; k < n_st; k++) {
                  const ParameterListT& st = list.GetList("surface_tension", k);
                  StringT ss_id;
                  st.GetParameter("side_set_ID", ss_id);
                  iArray2DT sides = model_manager.SideSet(ss_id);
                  const StringT& blk_id = model_manager.SideSetGroupID(ss_id);
                  if (sides.MajorDim() > 0) {
                        iArrayT elems(sides.MajorDim());
                        sides.ColumnCopy(0, elems);
                        BlockToGroupElementNumbers(elems, blk_id);
                        for (int j = 0; j < elems.Length(); j++)
                              extra_elems.AppendUnique(elems[j]);
                  }
            }

            /* check if any sideset elements are missing from fSurfaceElements */
            int nfs = fSurfaceElementNeighbors.MinorDim();
            iArrayT old_surf;
            old_surf.Dimension(fSurfaceElements.Length());
            old_surf.CopyPart(0, fSurfaceElements, 0, fSurfaceElements.Length());

            AutoArrayT<int> merged;
            for (int i = 0; i < old_surf.Length(); i++)
                  merged.AppendUnique(old_surf[i]);
            int old_len = merged.Length();
            for (int i = 0; i < extra_elems.Length(); i++)
                  merged.AppendUnique(extra_elems[i]);

            if (merged.Length() != old_len) {
                  /* rebuild: expand fSurfaceElements and fSurfaceElementNeighbors */
                  int new_len = merged.Length();
                  iArray2DT old_neighbors(fSurfaceElementNeighbors);

                  fSurfaceElements.Dimension(new_len);
                  fSurfaceElementNeighbors.Dimension(new_len, nfs);
                  fSurfaceElementNeighbors = -1; /* default: no neighbor (free surface) */

                  /* copy old data */
                  for (int i = 0; i < old_len; i++) {
                        fSurfaceElements[i] = merged[i];
                        for (int j = 0; j < nfs; j++)
                              fSurfaceElementNeighbors(i, j) = old_neighbors(i, j);
                  }
                  /* new entries: interior interface elements — mark all faces as free
                   * so surface tension can act on them */
                  for (int i = old_len; i < new_len; i++)
                        fSurfaceElements[i] = merged[i];
            }
      }

      /* initialize per-face surface tension arrays (zero = no contribution) */
      fSurfaceGamma.Dimension(fSurfaceElements.Length(), fSurfaceElementNeighbors.MinorDim());
      fSurfaceGamma = 0.0;
      fSurfaceT0.Dimension(fSurfaceElements.Length(), fSurfaceElementNeighbors.MinorDim());
      fSurfaceT0 = 0.0;

      /* parse surface_tension sub-lists and assign gamma/t_0 per face */
      if (n_st > 0) {
            InverseMapT st_elem_map;
            st_elem_map.SetOutOfRange(InverseMapT::Throw);
            st_elem_map.SetMap(fSurfaceElements);

            for (int k = 0; k < n_st; k++) {
                  const ParameterListT& st = list.GetList("surface_tension", k);
                  StringT ss_id;
                  double gamma, t_0;
                  st.GetParameter("side_set_ID", ss_id);
                  st.GetParameter("gamma", gamma);
                  st.GetParameter("t_0", t_0);

                  iArray2DT sides = model_manager.SideSet(ss_id);
                  const StringT& blk_id = model_manager.SideSetGroupID(ss_id);

                  if (sides.MajorDim() > 0) {
                        iArrayT elems(sides.MajorDim());
                        sides.ColumnCopy(0, elems);
                        BlockToGroupElementNumbers(elems, blk_id);
                        sides.SetColumn(0, elems);

                        for (int j = 0; j < sides.MajorDim(); j++) {
                              int elem = sides(j, 0);
                              int s_elem = st_elem_map.Map(elem);
                              int side = sides(j, 1);
                              fSurfaceGamma(s_elem, side) = gamma;
                              fSurfaceT0(s_elem, side) = t_0;
                        }
                  }
            }
      }

      /* determine normal type of each face */
      dMatrixT Q(nsd);
      dMatrixT jacobian(nsd, nsd-1);
      LocalArrayT face_coords(LocalArrayT::kCurrCoords, nfn, nsd); // kCurrCoords: current, kInitCoords: initial
      iArrayT face_nodes(nfn), face_nodes_index(nfn);
      ElementSupport().RegisterCoordinates(face_coords);
      fSurfaceElementFacesType = fSurfaceElementNeighbors;
      fSurfaceElementFacesType = -1;

      /* nodes per surface type */
      RowAutoFill2DT<int> nodes_on_surfaces(nfs, 25);

      for (int i = 0; i < fSurfaceElements.Length(); i++) /* loop over the surface elements for one element fSurfaceElements.Length() = 1 then i = 0*/
      {
            for (int j = 0; j < fSurfaceElementNeighbors.MinorDim(); j++) /* loop over faces */
            {
                  if (fSurfaceElementNeighbors(i,j) == -1) /* no neighbor => surface for one element (i,j) = (0,0), (0,1), (0,2) and (0,3)*/
                  {
                        /* face parent domain */
                        const ParentDomainT& surf_shape = shape.FacetShapeFunction(j);

                        /* collect coordinates of face nodes */
                        ElementCardT& element_card = ElementCard(fSurfaceElements[i]);
                        shape.NodesOnFacet(j, face_nodes_index);  // fni = 4 nodes of surface face
                        face_nodes.Collect(face_nodes_index, element_card.NodesX());
                        face_coords.SetLocal(face_nodes);

                        /* face normal (using 1st integration point) */
                        surf_shape.DomainJacobian(face_coords, 0, jacobian);
                        surf_shape.SurfaceJacobian(jacobian, Q);  // Last column of Q is normal vector to surface face

                        /* match to face normal types, i.e. normals_dat */
                        int normal_type = -1;
                        for (int k = 0; normal_type == -1 && k < fNormal.Length(); k++)
                        {
                              if ((Q.DotCol(nsd-1, fNormal[k]) - 1.0) < -kSmall)
                                    normal_type = -1;
                              else
                                    normal_type = k;
                        }
                        /* no match */
                        if (normal_type == -1)
                              cout << "could not classify normal on face " << j+1 << " of element " <<  fSurfaceElements[i]+1 << ". But that's ok" << endl;

                        /* store */
                        fSurfaceElementFacesType(i,j) = normal_type;

                        /* collect face nodes */
                        if (output_surface) nodes_on_surfaces.AppendUnique(normal_type, face_nodes);
                  }
            }
      }

      /* copy nodes on face types data */
      if (output_surface) {
            fSurfaceNodes.Dimension(nodes_on_surfaces.MajorDim());
            for (int i = 0; i < fSurfaceNodes.Length(); i++) {
                  fSurfaceNodes[i].Dimension(nodes_on_surfaces.MinorDim(i));
                  fSurfaceNodes[i].Copy(nodes_on_surfaces(i));
            }
      }

      /* process passivated surfaces */
      const ParameterListT* passivated_surfaces = list.List("passivated_surface_ID_list");
      if (passivated_surfaces) {
            ArrayT<StringT> ss_ID;
            StringListT::Extract(*passivated_surfaces, ss_ID);
            if (ss_ID.Length() > 0) {

                  /* need map into the surface element list */
                  InverseMapT surf_elem_map;
                  surf_elem_map.SetOutOfRange(InverseMapT::Throw);
                  surf_elem_map.SetMap(fSurfaceElements);

                  /* model manager */
                  ModelManagerT& model = ElementSupport().ModelManager();

                  /* loop over side set ID's */
                  for (int i = 0; i < ss_ID.Length(); i++) {
                        const StringT& id = ss_ID[i];

                          /* side set parameters */
                        iArray2DT sides = model.SideSet(id);
                        const StringT& block_id = model.SideSetGroupID(id);

                          if (sides.MajorDim() > 0) {

                              /* convert to element numbering within the group */
                              iArrayT elems(sides.MajorDim());
                              sides.ColumnCopy(0, elems);
                              BlockToGroupElementNumbers(elems, block_id);
                              sides.SetColumn(0, elems);

                              /* mark passivated faces */
                              for (int j = 0; j < sides.MajorDim(); j++) {
                                    int elem = sides(j,0);
                                    int s_elem = surf_elem_map.Map(elem);
                                    int side = sides(j,1);

                                    /* mark as non-surface (by giving negative neighbor id) */
                                    fSurfaceElementNeighbors(s_elem,side) = -2;
                              }
                        }
                  }
            }
      }

}

/***********************************************************************
 * Private
 ***********************************************************************/

/* ------------------------------------------------ calculate the LHS of residual, or element stiffness matrix ------------------------------------------------*/
void SimoQ1P0_Surface::FormStiffness(double constK)
{
	const char caller[] = "SimoQ1P0_Surface::FormStiffness";

	/* Inherited */
	SimoQ1P0::FormStiffness(constK);
	dMatrixT::SymmetryFlagT format = (fLHS.Format()
			== ElementMatrixT::kNonSymmetric)
            ? dMatrixT::kWhole
            : dMatrixT::kUpperOnly;

      const ShapeFunctionT& shape = ShapeFunction();
      int nsd = shape.NumSD();                          // # of spatial dimensions in problem
      int nfs = shape.NumFacets();                      // # of total possible element faces
      int nsi = shape.FacetShapeFunction(0).NumIP();    // # IPs per surface face
      int nfn = shape.FacetShapeFunction(0).NumNodes(); // # nodes on each surface face
      int nen = NumElementNodes();                      // # nodes in bulk element
      int nme = nen * nsd;

      iArrayT dof_indicies(nen);

      /* loop over surface elements */
      dMatrixT jacobian(nsd, nsd-1);
      iArrayT face_nodes(nfn), face_nodes_index(nfn);
      LocalArrayT face_coords(LocalArrayT::kCurrCoords, nfn, nsd);
      ElementSupport().RegisterCoordinates(face_coords);

      dMatrixT K1;
      dMatrixT K2;
      dMatrixT K_Total;
      dArrayT fB;
      dMatrixT K3(nen);

      K1.Dimension(nen, nen);
      K2.Dimension(nen, nen);
      fB.Dimension(nen);
      K_Total.Dimension(nen, nen);

      double CurrTime = ElementSupport().Time();

      fGrad_U.Dimension(2, NumSD());
      fGrad_U = 0.0;

 	 for (int i = 0; i < fSurfaceElements.Length(); i++)
 	 {
 		 int element = fSurfaceElements[i];
 		 if (element != CurrElementNumber()) continue; /* Skip if not current element */

 		 const ElementCardT& element_card = ElementCard(element);
 		 fLocInitCoords.SetLocal(element_card.NodesX());
 		 fLocDisp.SetLocal(element_card.NodesU());
 		 fAmm_mat2 = 0.0;

 		 for (int j = 0; j < fSurfaceElementNeighbors.MinorDim(); j++) /* loop over faces */
 		 {
 			 /* per-face surface tension with optional ramp-up */
 			 double gamma = fSurfaceGamma(i, j);
 			 if (gamma == 0.0) continue; /* Skip faces without surface tension */
 			 
 			 double t_0 = fSurfaceT0(i, j);
 			 double scaled_gamma = (t_0 > 0.0) ? gamma * min(1.0, CurrTime / t_0) : gamma;

 			 /* collect coordinates of face nodes */
 			 shape.NodesOnFacet(j, face_nodes_index);
 			 face_nodes.Collect(face_nodes_index, element_card.NodesX());
 			 face_coords.SetLocal(face_nodes);

 			 double x_1 = face_coords[0];
 			 double x_2 = face_coords[1];
 			 double y_1 = face_coords[2];
 			 double y_2 = face_coords[3];

 			 double dx = x_1 - x_2;
 			 double dy = y_1 - y_2;
 			 double edge_length_sq = dx*dx + dy*dy;
 			 double edge_length = sqrt(edge_length_sq);

 			 /* Stiffness matrix: K = coeff1 * I + coeff2 * (b ⊗ b) */
 			 double coeff1 = scaled_gamma / edge_length;
 			 double coeff2 = -scaled_gamma / (edge_length * edge_length_sq);

 			 /* Build K_Total without intermediate matrices */
 			 K1 = 0.0;
 			 K1(0,0) = 1.0;   K1(1,1) = 1.0;
 			 K1(2,2) = 1.0;   K1(3,3) = 1.0;
 			 K1(0,2) = -1.0;  K1(2,0) = -1.0;
 			 K1(1,3) = -1.0;  K1(3,1) = -1.0;
 			 K1 *= coeff1;

 			 K2 = 0.0;
 			 K2(0,0) = dx*dx;    K2(0,1) = dx*dy;    K2(0,2) = -dx*dx;   K2(0,3) = -dx*dy;
 			 K2(1,0) = dy*dx;    K2(1,1) = dy*dy;    K2(1,2) = -dy*dx;   K2(1,3) = -dy*dy;
 			 K2(2,0) = -dx*dx;   K2(2,1) = -dx*dy;   K2(2,2) = dx*dx;    K2(2,3) = dx*dy;
 			 K2(3,0) = -dy*dx;   K2(3,1) = -dy*dy;   K2(3,2) = dy*dx;    K2(3,3) = dy*dy;
 			 K2 *= coeff2;

 			 K_Total = K1;
 			 K_Total += K2;

 			 dof_indicies = CanonicalNodes(face_nodes_index[0], face_nodes_index[1]);

 			 for (int n = 0; n < 4; n++) {
 			 	 for (int m = 0; m < 4; m++) {
 			 	 	 fAmm_mat2(dof_indicies[n], dof_indicies[m]) += K_Total(n, m);
 			 	 }
 			 }
 		 } /* End of surface edge loop */

 		 /* scale by integrator coefficient and add to LHS */
 		 fAmm_mat2 *= constK;
 		 fLHS.AddBlock(0, 0, fAmm_mat2);
 	 } /* End of element loop */

} /* End of Function FormStiffness */

/* ------------------------------------------------ Compute RHS, or residual of element equations ------------------------------------------------*/
void SimoQ1P0_Surface::FormKd(double constK)
{
      const char caller[] = "SimoQ1P0_Surface::FormKd";

      /* Inherited */
      SimoQ1P0::FormKd(constK);

      /* dimensions */
      const ShapeFunctionT& shape = ShapeFunction();
      int nsd = shape.NumSD();                          // # of spatial dimensions in problem
      int nfs = shape.NumFacets();                      // # of total possible element faces
      int nsi = shape.FacetShapeFunction(0).NumIP();    // # IPs per surface face
      int nfn = shape.FacetShapeFunction(0).NumNodes(); // # nodes on each surface face
      int nen = NumElementNodes();                      // # nodes in bulk element
      int nme = nen * nsd;

      iArrayT dof_indicies(nen);

      /* loop over surface elements */
      dMatrixT jacobian(nsd, nsd-1);
      LocalArrayT face_coords(LocalArrayT::kCurrCoords, nfn, nsd); // kCurrCoords = current coordinates
      iArrayT face_nodes(nfn), face_nodes_index(nfn);
      ElementSupport().RegisterCoordinates(face_coords);

      dArrayT fD;
      dArrayT R_Total;
      dArrayT R;

      fD.Dimension(nen);
      R_Total.Dimension(nme);
      R.Dimension(fRHS.Length());

      double CurrTime = ElementSupport().Time();

      for (int i = 0; i < fSurfaceElements.Length(); i++)
      {
            int element = fSurfaceElements[i];
            if (element != CurrElementNumber()) continue; /* Skip if not current element */

            const ElementCardT& element_card = ElementCard(element);
            fLocInitCoords.SetLocal(element_card.NodesX());
            fLocDisp.SetLocal(element_card.NodesU());

            fD = 0.0;
            R_Total = 0.0;

            for (int j = 0; j < fSurfaceElementNeighbors.MinorDim(); j++) /* loop over faces */
            {
                  /* per-face surface tension with optional ramp-up */
                  double gamma = fSurfaceGamma(i, j);
                  if (gamma == 0.0) continue; /* Skip faces without surface tension */
                  
                  double t_0 = fSurfaceT0(i, j);
                  double scaled_gamma = (t_0 > 0.0) ? gamma * min(1.0, CurrTime / t_0) : gamma;

                  /* collect coordinates of face nodes */
                  shape.NodesOnFacet(j, face_nodes_index);
                  face_nodes.Collect(face_nodes_index, element_card.NodesX());
                  face_coords.SetLocal(face_nodes);

                  double x_1 = face_coords[0];
                  double x_2 = face_coords[1];
                  double y_1 = face_coords[2];
                  double y_2 = face_coords[3];

                  double dx = x_1 - x_2;
                  double dy = y_1 - y_2;
                  double edge_length = sqrt(dx*dx + dy*dy);
                  double scale_factor = -scaled_gamma / edge_length;

                  fD[0] = scale_factor * dx;
                  fD[1] = scale_factor * dy;
                  fD[2] = -fD[0];
                  fD[3] = -fD[1];

                  dof_indicies = CanonicalNodes(face_nodes_index[0], face_nodes_index[1]);

                  for (int k = 0; k < 4; k++)
                        R_Total[dof_indicies[k]] += fD[k];
            }

            R.CopyIn(0, R_Total);
            fRHS += R;
      }
}
/***********************************************************************
 * Protected
 ***********************************************************************/
iArrayT SimoQ1P0_Surface::CanonicalNodes(const int node_index0, const int node_index1)
{
	const char caller[] = "SimoQ1P0_Surface::CanonicalNodes";

	/* Tahoe's fRHS/fLHS use INTERLEAVED DOF ordering (confirmed via FieldT::SetLocalEqnos):
	 * [u_x0, u_y0, u_x1, u_y1, u_x2, u_y2, u_x3, u_y3]
	 * Face (n0,n1) DOFs: [x_n0=2*n0, y_n0=2*n0+1, x_n1=2*n1, y_n1=2*n1+1] */
	iArrayT dof_indicies(4);   /* 4 entries: [x_n0, y_n0, x_n1, y_n1] */

	if (node_index0 == 0 && node_index1 == 1) {
		dof_indicies[0] = 0; dof_indicies[1] = 1; dof_indicies[2] = 2; dof_indicies[3] = 3;
	} else if (node_index0 == 1 && node_index1 == 2) {
		dof_indicies[0] = 2; dof_indicies[1] = 3; dof_indicies[2] = 4; dof_indicies[3] = 5;
	} else if (node_index0 == 2 && node_index1 == 3) {
		dof_indicies[0] = 4; dof_indicies[1] = 5; dof_indicies[2] = 6; dof_indicies[3] = 7;
	} else if (node_index0 == 3 && node_index1 == 0) {
		dof_indicies[0] = 6; dof_indicies[1] = 7; dof_indicies[2] = 0; dof_indicies[3] = 1;
	} else {
		ExceptionT::GeneralFail(caller, "could not classify face with node index %d and %d", node_index0, node_index1);
	}

	return dof_indicies;
}


/***********************************************************************
 * Private
 ***********************************************************************/

  // ------------------------------------------------ extrapolate from integration points and compute output nodal/element values ------------------------ //
void SimoQ1P0_Surface::ComputeOutput(const iArrayT& n_codes,
      dArray2DT& n_values, const iArrayT& e_codes, dArray2DT& e_values)
{
      // cout << "\033[1;31mCan you see this?\033[0m" << endl;
      /* inherited - bulk contribution */
      SimoQ1P0::ComputeOutput(n_codes, n_values, e_codes, e_values);

}

} // namespace Tahoe
