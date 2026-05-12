#include "FSDielectricElastomerQ1P0ElastocapillaryT.h"
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

FSDielectricElastomerQ1P0ElastocapillaryT::FSDielectricElastomerQ1P0ElastocapillaryT(const ElementSupportT& support) :
    FSDielectricElastomerQ1P02DT(support),
    fLocCurrCoords(LocalArrayT::kCurrCoords), fSurfTension(0), fSurfaceCBSupport(NULL)
{
      SetName("dielectric_elastomer_Q1P0Elastocapillary");
}


/* Destructor */
FSDielectricElastomerQ1P0ElastocapillaryT::~FSDielectricElastomerQ1P0ElastocapillaryT()
{
      /* TLCBSurface Stuff */
	/* free surface models */
	delete fSurfaceCBSupport;

}


// ------------------------------------------------specify parameters needed by the interface--------------------------------------------------------- //
void FSDielectricElastomerQ1P0ElastocapillaryT::DefineParameters(ParameterListT& list) const
{
      /* Inherited from the most similar element */
      FSDielectricElastomerQ1P02DT::DefineParameters(list);

      /* Taking the output_surface */
      ParameterT output_surface(ParameterT::Boolean, "output_surface");
      output_surface.SetDefault(false);
      list.AddParameter(output_surface);
}

  /* ----------------------------- information about subordinate parameter lists ------------------------------------------------------------------------*/
void FSDielectricElastomerQ1P0ElastocapillaryT::DefineSubs(SubListT& sub_list) const
{

      /* inherited from the most similar element */
      FSDielectricElastomerQ1P02DT::DefineSubs(sub_list);

      /* list of passivated surfaces - side set list */
      sub_list.AddSub("passivated_surface_ID_list", ParameterListT::ZeroOrOnce);

}


  // ---------------------------------------------------- accept parameter list ------------------------------------------------------------------------ //
void FSDielectricElastomerQ1P0ElastocapillaryT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "FSDielectricElastomerQ1P0ElastocapillaryT::TakeParameterList";

      /* inherited from the most similar element for bulk element */
      FSDielectricElastomerQ1P02DT::TakeParameterList(list);

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

      /* find parameter list for the bulk material */
      int num_blocks = list.NumLists("large_strain_element_block");
      if (num_blocks > 1)
      ExceptionT::GeneralFail(caller, "expecting only 1 not %d element blocks", num_blocks);
      const ParameterListT& block = list.GetList("large_strain_element_block");
      const ParameterListT& mat_list = block.GetListChoice(*this, "large_strain_material_choice");
      const ArrayT<ParameterListT>& mat = mat_list.Lists();
      const ParameterListT& bulk_params = mat[0]; /* getting bulk parameters */

      fSurfTension = bulk_params.GetParameter("gamma"); /* Reading the value of gamma (surface tension) */


      fT_0 = bulk_params.GetParameter("t_0");

      bool output_surface = list.GetParameter("output_surface");

      /* collect surface element information */
      ArrayT<StringT> block_ID;
      ElementBlockIDs(block_ID);
      ModelManagerT& model_manager = ElementSupport().ModelManager();
      model_manager.BoundingElements(block_ID, fSurfaceElements, fSurfaceElementNeighbors);



      ArrayT<const iArray2DT*> connects;
      model_manager.ElementGroupPointers(block_ID, connects);
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

/* Surface-tension contribution to the element stiffness.
 *
 * For each free surface face (a line in 2D between face nodes A and B):
 *     K_face = (gamma/L) * P  -  (gamma/L^3) * (b ⊗ b)
 * where L = |x_A - x_B|, b = (dx, dy, -dx, -dy), and P is the 4x4 line
 * stiffness pattern [[I, -I], [-I, I]].  The 4x4 block is scattered into
 * the bulk 8x8 mechanical block using the flat DOF indices of the 2 face
 * nodes inside the 4-node bulk element.
 */
void FSDielectricElastomerQ1P0ElastocapillaryT::FormStiffness(double constK)
{
    /* Inherited bulk contribution */
    FSDielectricElastomerQ1P02DT::FormStiffness(constK);

    /* skip elements that aren't on the free surface */
    const int curr_elem = CurrElementNumber();
    int s_idx = -1;
    for (int i = 0; i < fSurfaceElements.Length(); i++)
        if (fSurfaceElements[i] == curr_elem) { s_idx = i; break; }
    if (s_idx == -1) return;

    const ShapeFunctionT& shape = ShapeFunction();
    const int nsd = shape.NumSD();
    const int nfn = shape.FacetShapeFunction(0).NumNodes();

    /* time-ramped surface tension */
    const double t = ElementSupport().Time();
    const double gamma = min(fSurfTension, t * fSurfTension / fT_0);

    /* P = [[I, -I], [-I, I]]: constant 4x4 line-stiffness pattern */
    dMatrixT P(4, 4);
    P = 0.0;
    P(0,0) = P(1,1) = P(2,2) = P(3,3) =  1.0;
    P(0,2) = P(2,0) = P(1,3) = P(3,1) = -1.0;

    iArrayT face_nodes(nfn), face_nodes_index(nfn);
    LocalArrayT face_coords(LocalArrayT::kCurrCoords, nfn, nsd);
    ElementSupport().RegisterCoordinates(face_coords);

    dArrayT b(4);
    dMatrixT bb(4, 4), K_face(4, 4);
    fAmm_mat2 = 0.0;

    for (int j = 0; j < fSurfaceElementNeighbors.MinorDim(); j++) {
        if (fSurfaceElementNeighbors(s_idx, j) != -1) continue;  /* not a free face */

        ElementCardT& element_card = ElementCard(curr_elem);
        shape.NodesOnFacet(j, face_nodes_index);
        face_nodes.Collect(face_nodes_index, element_card.NodesX());
        face_coords.SetLocal(face_nodes);

        const double dx = face_coords[0] - face_coords[1];
        const double dy = face_coords[2] - face_coords[3];
        const double L_e = sqrt(dx*dx + dy*dy);

        b[0] =  dx;  b[1] =  dy;
        b[2] = -dx;  b[3] = -dy;
        bb.Outer(b, b);

        K_face = P;
        K_face *= gamma / L_e;
        bb     *= -gamma / (L_e * L_e * L_e);
        K_face += bb;

        const iArrayT face_dofs = FaceDOFIndices(face_nodes_index[0], face_nodes_index[1]);
        for (int n = 0; n < 4; n++)
            for (int m = 0; m < 4; m++)
                fAmm_mat2(face_dofs[n], face_dofs[m]) += K_face(n, m);
    }

    /* dynamic formulation: scale once after all faces are accumulated */
    if (fIntegrator->Order() == 2)
        fAmm_mat2 *= constK;

    fLHS.AddBlock(0, 0, fAmm_mat2);
}

/* Surface-tension contribution to the element residual.
 *
 * For each free surface face (a line in 2D from node A to node B):
 *     R_face = -(gamma/L) * (x_A - x_B, y_A - y_B, x_B - x_A, y_B - y_A)
 * which is the gradient of the line-energy (gamma*L) w.r.t. the 4 face DOFs.
 */
void FSDielectricElastomerQ1P0ElastocapillaryT::FormKd(double constK)
{
    /* Inherited bulk contribution */
    FSDielectricElastomerQ1P02DT::FormKd(constK);

    /* skip elements that aren't on the free surface */
    const int curr_elem = CurrElementNumber();
    int s_idx = -1;
    for (int i = 0; i < fSurfaceElements.Length(); i++)
        if (fSurfaceElements[i] == curr_elem) { s_idx = i; break; }
    if (s_idx == -1) return;

    const ShapeFunctionT& shape = ShapeFunction();
    const int nsd = shape.NumSD();
    const int nfn = shape.FacetShapeFunction(0).NumNodes();
    const int nen = NumElementNodes();
    const int nme = nen * nsd;

    const double t = ElementSupport().Time();
    const double gamma = min(fSurfTension, t * fSurfTension / fT_0);

    iArrayT face_nodes(nfn), face_nodes_index(nfn);
    LocalArrayT face_coords(LocalArrayT::kCurrCoords, nfn, nsd);
    ElementSupport().RegisterCoordinates(face_coords);

    dArrayT R_mech(nme), R((nsd + 1) * nen);
    R_mech = 0.0;
    R      = 0.0;

    for (int j = 0; j < fSurfaceElementNeighbors.MinorDim(); j++) {
        if (fSurfaceElementNeighbors(s_idx, j) != -1) continue;  /* not a free face */

        ElementCardT& element_card = ElementCard(curr_elem);
        shape.NodesOnFacet(j, face_nodes_index);
        face_nodes.Collect(face_nodes_index, element_card.NodesX());
        face_coords.SetLocal(face_nodes);

        const double dx = face_coords[0] - face_coords[1];
        const double dy = face_coords[2] - face_coords[3];
        const double L_e = sqrt(dx*dx + dy*dy);
        const double scale = -gamma / L_e;

        const iArrayT face_dofs = FaceDOFIndices(face_nodes_index[0], face_nodes_index[1]);
        R_mech[face_dofs[0]] += scale *  dx;
        R_mech[face_dofs[1]] += scale *  dy;
        R_mech[face_dofs[2]] += scale * -dx;
        R_mech[face_dofs[3]] += scale * -dy;
    }

    R.CopyIn(0, R_mech);
    fRHS += R;
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* Flat DOF indices for the 2 nodes of a quad4 face in the 8-DOF mechanical
 * block (layout: [u_x_0, u_y_0, u_x_1, u_y_1, ..., u_y_3]).
 *
 *   node A -> DOFs (2A, 2A+1)
 *   node B -> DOFs (2B, 2B+1)
 */
iArrayT FSDielectricElastomerQ1P0ElastocapillaryT::FaceDOFIndices(int node_a, int node_b)
{
    iArrayT face_dofs(4);
    face_dofs[0] = 2 * node_a;
    face_dofs[1] = 2 * node_a + 1;
    face_dofs[2] = 2 * node_b;
    face_dofs[3] = 2 * node_b + 1;
    return face_dofs;
}


/***********************************************************************
 * Private
 ***********************************************************************/

  // ------------------------------------------------ extrapolate from integration points and compute output nodal/element values ------------------------ //
void FSDielectricElastomerQ1P0ElastocapillaryT::ComputeOutput(const iArrayT& n_codes,
      dArray2DT& n_values, const iArrayT& e_codes, dArray2DT& e_values)
{
      /* inherited - bulk contribution */
      FSDielectricElastomerQ1P02DT::ComputeOutput(n_codes, n_values, e_codes, e_values);

}

} // namespace Tahoe
