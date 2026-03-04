/* $Id: SimoQ1P0_3D_Surface.cpp,v 1.1 2026/03/02 00:00:00 samanseifi Exp $ */
/*
 * Updated Lagrangian Q1P0 3D hex element with surface tension.
 *
 * Surface tension formulation follows the covariant metric tensor approach:
 *   W_s = gamma * int H dxi1 dxi2,  H = sqrt(EG - F^2)
 * using 2x2 Gauss quadrature on each free hex face.
 *
 * Reference: Henann & Bertoldi (2014) uel_surften_3D_hex.for
 */
#include "SimoQ1P0_3D_Surface.h"

#include "ShapeFunctionT.h"
#include "ParentDomainT.h"

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

#include <cmath>
#include <iostream>
#include <algorithm>

using namespace std;

namespace Tahoe {

SimoQ1P0_3D_Surface::SimoQ1P0_3D_Surface(const ElementSupportT& support) :
    SimoQ1P0(support),
    fLocCurrCoords(LocalArrayT::kCurrCoords), fSurfaceCBSupport(NULL)
{
    SetName("updated_lagrangian_Q1P0_3D_surface");
}

SimoQ1P0_3D_Surface::~SimoQ1P0_3D_Surface()
{
    delete fSurfaceCBSupport;
}

void SimoQ1P0_3D_Surface::DefineParameters(ParameterListT& list) const
{
    SimoQ1P0::DefineParameters(list);

    ParameterT output_surface(ParameterT::Boolean, "output_surface");
    output_surface.SetDefault(false);
    list.AddParameter(output_surface);
}

void SimoQ1P0_3D_Surface::DefineSubs(SubListT& sub_list) const
{
    SimoQ1P0::DefineSubs(sub_list);
    sub_list.AddSub("passivated_surface_ID_list", ParameterListT::ZeroOrOnce);
    sub_list.AddSub("surface_tension", ParameterListT::Any);
}

ParameterInterfaceT* SimoQ1P0_3D_Surface::NewSub(const StringT& name) const
{
    if (name == "surface_tension") {
        ParameterContainerT* st = new ParameterContainerT(name);
        st->AddParameter(ParameterT::Word, "side_set_ID");
        st->AddParameter(ParameterT::Double, "gamma");
        st->AddParameter(0.0, "t_0");
        return st;
    }
    return SimoQ1P0::NewSub(name);
}

void SimoQ1P0_3D_Surface::TakeParameterList(const ParameterListT& list)
{
    const char caller[] = "SimoQ1P0_3D_Surface::TakeParameterList";

    SimoQ1P0::TakeParameterList(list);

    const ShapeFunctionT& shape = ShapeFunction();
    int nsd = shape.NumSD();                           // should be 3
    int nfs = shape.NumFacets();                       // should be 6 for hex
    int nsi = shape.FacetShapeFunction(0).NumIP();     // should be 4 (2x2 Gauss)
    int nfn = shape.FacetShapeFunction(0).NumNodes();  // should be 4 (quad face)
    int nen = NumElementNodes();                       // should be 8 for hex
    int nme = nen * nsd;                               // 24

    if (nsd != 3)
        ExceptionT::GeneralFail(caller, "expected 3D element, got nsd=%d", nsd);

    fAmm_mat2.Dimension(nme, nme);

    /* surface model support (needed to satisfy parent interface) */
    fF_Surf_List.Dimension(nsi);
    for (int i = 0; i < fF_Surf_List.Length(); i++)
        fF_Surf_List[i].Dimension(nsd);
    fF_Surf_List[0].Identity();

    fSurfaceCBSupport = new FSMatSupportT(nsd, nsi);
    fSurfaceCBSupport->SetContinuumElement(this);
    fSurfaceCBSupport->SetDeformationGradient(&fF_Surf_List);

    /* canonical normals for axis-aligned hex faces (+/-z, +/-y, +/-x) */
    double normals_dat[6*3] = {
         0.0,  0.0, -1.0,  /* Face 0: bottom (z-) */
         0.0,  0.0,  1.0,  /* Face 1: top    (z+) */
         0.0, -1.0,  0.0,  /* Face 2: front  (y-) */
         1.0,  0.0,  0.0,  /* Face 3: right  (x+) */
         0.0,  1.0,  0.0,  /* Face 4: back   (y+) */
        -1.0,  0.0,  0.0,  /* Face 5: left   (x-) */
    };
    dArray2DT normals(6, 3, normals_dat);
    fNormal.Dimension(nfs);
    for (int i = 0; i < nfs; i++) {
        fNormal[i].Dimension(nsd);
        fNormal[i] = normals(i);
    }

    /* collect surface elements */
    ArrayT<StringT> block_ID;
    ElementBlockIDs(block_ID);
    ModelManagerT& model_manager = ElementSupport().ModelManager();
    model_manager.BoundingElements(block_ID, fSurfaceElements, fSurfaceElementNeighbors);

    /* initialize per-face surface tension arrays */
    fSurfaceGamma.Dimension(fSurfaceElements.Length(), fSurfaceElementNeighbors.MinorDim());
    fSurfaceGamma = 0.0;
    fSurfaceT0.Dimension(fSurfaceElements.Length(), fSurfaceElementNeighbors.MinorDim());
    fSurfaceT0 = 0.0;

    /* parse surface_tension sub-lists */
    int n_st = list.NumLists("surface_tension");
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

    /* classify faces by normal direction (for output_surface feature) */
    bool output_surface = list.GetParameter("output_surface");
    dMatrixT Q(nsd);
    dMatrixT jacobian(nsd, nsd-1);
    LocalArrayT face_coords(LocalArrayT::kCurrCoords, nfn, nsd);
    iArrayT face_nodes(nfn), face_nodes_index(nfn);
    ElementSupport().RegisterCoordinates(face_coords);
    fSurfaceElementFacesType = fSurfaceElementNeighbors;
    fSurfaceElementFacesType = -1;

    RowAutoFill2DT<int> nodes_on_surfaces(nfs, 25);

    for (int i = 0; i < fSurfaceElements.Length(); i++) {
        for (int j = 0; j < fSurfaceElementNeighbors.MinorDim(); j++) {
            if (fSurfaceElementNeighbors(i,j) == -1) {
                const ParentDomainT& surf_shape = shape.FacetShapeFunction(j);

                ElementCardT& element_card = ElementCard(fSurfaceElements[i]);
                shape.NodesOnFacet(j, face_nodes_index);
                face_nodes.Collect(face_nodes_index, element_card.NodesX());
                face_coords.SetLocal(face_nodes);

                surf_shape.DomainJacobian(face_coords, 0, jacobian);
                surf_shape.SurfaceJacobian(jacobian, Q);

                int normal_type = -1;
                for (int k = 0; normal_type == -1 && k < fNormal.Length(); k++) {
                    if ((Q.DotCol(nsd-1, fNormal[k]) - 1.0) >= -kSmall)
                        normal_type = k;
                }
                if (normal_type == -1)
                    cout << "could not classify normal on face " << j+1
                         << " of element " << fSurfaceElements[i]+1
                         << ". But that's ok" << endl;

                fSurfaceElementFacesType(i,j) = normal_type;
                if (output_surface) nodes_on_surfaces.AppendUnique(normal_type, face_nodes);
            }
        }
    }

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
            InverseMapT surf_elem_map;
            surf_elem_map.SetOutOfRange(InverseMapT::Throw);
            surf_elem_map.SetMap(fSurfaceElements);

            ModelManagerT& model = ElementSupport().ModelManager();
            for (int i = 0; i < ss_ID.Length(); i++) {
                const StringT& id = ss_ID[i];
                iArray2DT sides = model.SideSet(id);
                const StringT& block_id = model.SideSetGroupID(id);

                if (sides.MajorDim() > 0) {
                    iArrayT elems(sides.MajorDim());
                    sides.ColumnCopy(0, elems);
                    BlockToGroupElementNumbers(elems, block_id);
                    sides.SetColumn(0, elems);

                    for (int j = 0; j < sides.MajorDim(); j++) {
                        int elem = sides(j,0);
                        int s_elem = surf_elem_map.Map(elem);
                        int side = sides(j,1);
                        fSurfaceElementNeighbors(s_elem,side) = -2;
                    }
                }
            }
        }
    }
}

/* ---- FormKd: surface tension residual contribution ---- */
void SimoQ1P0_3D_Surface::FormKd(double constK)
{
    /* inherited bulk contribution */
    SimoQ1P0::FormKd(constK);

    const ShapeFunctionT& shape = ShapeFunction();
    int nsd = 3;
    int nfn = 4;          /* nodes per quad face */
    int nen = NumElementNodes();   /* 8 for hex */
    int nme = nen * nsd;           /* 24 */

    dMatrixT jacobian(nsd, nsd-1);  /* 3x2: columns = g1, g2 */
    iArrayT face_nodes(nfn), face_nodes_index(nfn);
    LocalArrayT face_coords(LocalArrayT::kCurrCoords, nfn, nsd);
    ElementSupport().RegisterCoordinates(face_coords);

    iArrayT dof_indices(nfn * nsd);  /* 12 element-level DOF indices */
    dArrayT R_Total(nme);
    dArrayT R(fRHS.Length());

    double g1[3], g2[3];

    double CurrTime = ElementSupport().Time();

    for (int i = 0; i < fSurfaceElements.Length(); i++) {
        int element = fSurfaceElements[i];
        if (element != CurrElementNumber()) continue;

        const ElementCardT& element_card = ElementCard(element);
        fLocInitCoords.SetLocal(element_card.NodesX());
        fLocDisp.SetLocal(element_card.NodesU());

        R_Total = 0.0;

        for (int j = 0; j < fSurfaceElementNeighbors.MinorDim(); j++) {
            if (fSurfaceElementNeighbors(i,j) != -1) continue;

            double gamma = fSurfaceGamma(i, j);
            if (gamma == 0.0) continue;

            double t_0 = fSurfaceT0(i, j);
            double scaled_gamma = (t_0 > 0.0) ? gamma * min(1.0, CurrTime / t_0) : gamma;

            const ParentDomainT& surf_shape = shape.FacetShapeFunction(j);
            int nsi = surf_shape.NumIP();  /* 4 for 2x2 Gauss */

            shape.NodesOnFacet(j, face_nodes_index);
            face_nodes.Collect(face_nodes_index, element_card.NodesX());
            face_coords.SetLocal(face_nodes);

            CanonicalNodes3D(face_nodes_index, dof_indices);

            const double* w = surf_shape.Weight();

            for (int ip = 0; ip < nsi; ip++) {
                /* covariant basis vectors from DomainJacobian */
                surf_shape.DomainJacobian(face_coords, ip, jacobian);
                for (int d = 0; d < nsd; d++) {
                    g1[d] = jacobian(d, 0);
                    g2[d] = jacobian(d, 1);
                }

                /* covariant metric tensor components */
                double E = 0.0, G = 0.0, F_met = 0.0;
                for (int d = 0; d < nsd; d++) {
                    E     += g1[d]*g1[d];
                    G     += g2[d]*g2[d];
                    F_met += g1[d]*g2[d];
                }
                double H = sqrt(E*G - F_met*F_met);

                /* shape function derivatives at this IP */
                const double* dN1 = surf_shape.DShape(ip, 0);  /* dN/dxi1, length nfn */
                const double* dN2 = surf_shape.DShape(ip, 1);  /* dN/dxi2, length nfn */

                double coeff = -scaled_gamma * w[ip] / H;

                /* residual contribution: R_a = -(gamma/H)*w * [(G*n1a - F*n2a)*g1 + (E*n2a - F*n1a)*g2] */
                for (int a = 0; a < nfn; a++) {
                    double n1a = dN1[a], n2a = dN2[a];
                    double c1 = G*n1a - F_met*n2a;
                    double c2 = E*n2a - F_met*n1a;
                    for (int d = 0; d < nsd; d++) {
                        R_Total[dof_indices[3*a + d]] += coeff * (c1*g1[d] + c2*g2[d]);
                    }
                }
            }  /* end IP loop */
        }  /* end face loop */

        R = 0.0;
        R.CopyIn(0, R_Total);
        fRHS += R;
    }  /* end element loop */
}

/* ---- FormStiffness: surface tension tangent stiffness contribution ---- */
void SimoQ1P0_3D_Surface::FormStiffness(double constK)
{
    /* inherited bulk contribution */
    SimoQ1P0::FormStiffness(constK);

    const ShapeFunctionT& shape = ShapeFunction();
    int nsd = 3;
    int nfn = 4;
    int nen = NumElementNodes();
    int nme = nen * nsd;

    dMatrixT jacobian(nsd, nsd-1);
    iArrayT face_nodes(nfn), face_nodes_index(nfn);
    LocalArrayT face_coords(LocalArrayT::kCurrCoords, nfn, nsd);
    ElementSupport().RegisterCoordinates(face_coords);

    iArrayT dof_indices(nfn * nsd);

    double g1[3], g2[3];
    double Va[3], Vb[3];

    double CurrTime = ElementSupport().Time();

    for (int i = 0; i < fSurfaceElements.Length(); i++) {
        int element = fSurfaceElements[i];
        if (element != CurrElementNumber()) continue;

        const ElementCardT& element_card = ElementCard(element);
        fLocInitCoords.SetLocal(element_card.NodesX());
        fLocDisp.SetLocal(element_card.NodesU());
        fAmm_mat2 = 0.0;

        for (int j = 0; j < fSurfaceElementNeighbors.MinorDim(); j++) {
            if (fSurfaceElementNeighbors(i,j) != -1) continue;

            double gamma = fSurfaceGamma(i, j);
            if (gamma == 0.0) continue;

            double t_0 = fSurfaceT0(i, j);
            double scaled_gamma = (t_0 > 0.0) ? gamma * min(1.0, CurrTime / t_0) : gamma;

            const ParentDomainT& surf_shape = shape.FacetShapeFunction(j);
            int nsi = surf_shape.NumIP();

            shape.NodesOnFacet(j, face_nodes_index);
            face_nodes.Collect(face_nodes_index, element_card.NodesX());
            face_coords.SetLocal(face_nodes);

            CanonicalNodes3D(face_nodes_index, dof_indices);

            const double* w = surf_shape.Weight();

            for (int ip = 0; ip < nsi; ip++) {
                surf_shape.DomainJacobian(face_coords, ip, jacobian);
                for (int d = 0; d < nsd; d++) {
                    g1[d] = jacobian(d, 0);
                    g2[d] = jacobian(d, 1);
                }

                double E = 0.0, G = 0.0, F_met = 0.0;
                for (int d = 0; d < nsd; d++) {
                    E     += g1[d]*g1[d];
                    G     += g2[d]*g2[d];
                    F_met += g1[d]*g2[d];
                }
                double H = sqrt(E*G - F_met*F_met);
                double H2 = H*H;

                const double* dN1 = surf_shape.DShape(ip, 0);
                const double* dN2 = surf_shape.DShape(ip, 1);

                double pref = scaled_gamma * w[ip] / H;
                double pref_curv = scaled_gamma * w[ip] / (H2 * H);

                for (int a = 0; a < nfn; a++) {
                    double n1a = dN1[a], n2a = dN2[a];
                    double c1a = G*n1a - F_met*n2a;
                    double c2a = E*n2a - F_met*n1a;
                    for (int d = 0; d < nsd; d++)
                        Va[d] = c1a*g1[d] + c2a*g2[d];

                    for (int b = 0; b < nfn; b++) {
                        double n1b = dN1[b], n2b = dN2[b];
                        double c1b = G*n1b - F_met*n2b;
                        double c2b = E*n2b - F_met*n1b;
                        for (int d = 0; d < nsd; d++)
                            Vb[d] = c1b*g1[d] + c2b*g2[d];

                        /* Term 1 scalar: (G*n1a*n1b - F*(n1a*n2b+n2a*n1b) + E*n2a*n2b) */
                        double s_ab = G*n1a*n1b - F_met*(n1a*n2b + n2a*n1b) + E*n2a*n2b;

                        int row_base = 3*a;
                        int col_base = 3*b;

                        for (int d = 0; d < nsd; d++) {
                            for (int e = 0; e < nsd; e++) {
                                /* Term 1: identity block */
                                double K_de = pref * s_ab * (d == e ? 1.0 : 0.0);

                                /* Term 2: outer products from metric tensor variation
                                 * = n1a*g1[d] * (2*n2b*g2[e])
                                 *   - (n2a*g1[d] + n1a*g2[d]) * (n1b*g2[e] + n2b*g1[e])
                                 *   + n2a*g2[d] * (2*n1b*g1[e])                         */
                                K_de += pref * (
                                      n1a*g1[d] * 2.0*n2b*g2[e]
                                    - (n2a*g1[d] + n1a*g2[d]) * (n1b*g2[e] + n2b*g1[e])
                                    + n2a*g2[d] * 2.0*n1b*g1[e]
                                );

                                /* Term 3: curvature — Va ⊗ Vb / H^2 */
                                K_de -= pref_curv * Va[d] * Vb[e];

                                fAmm_mat2(dof_indices[row_base + d],
                                          dof_indices[col_base + e]) += K_de;
                            }
                        }
                    }  /* end b loop */
                }  /* end a loop */
            }  /* end IP loop */
        }  /* end face loop */

        fAmm_mat2 *= constK;
        fLHS.AddBlock(0, 0, fAmm_mat2);
    }  /* end element loop */
}

/* ---- CanonicalNodes3D: element DOF indices for a quad face ---- */
void SimoQ1P0_3D_Surface::CanonicalNodes3D(const iArrayT& face_nodes_index,
                                            iArrayT& dof_indices) const
{
    /* INTERLEAVED DOF ordering: node n, direction d → index 3*n + d */
    for (int a = 0; a < 4; a++) {
        int n = face_nodes_index[a];
        dof_indices[3*a + 0] = 3*n + 0;  /* x */
        dof_indices[3*a + 1] = 3*n + 1;  /* y */
        dof_indices[3*a + 2] = 3*n + 2;  /* z */
    }
}

/* ---- ComputeOutput ---- */
void SimoQ1P0_3D_Surface::ComputeOutput(const iArrayT& n_codes,
    dArray2DT& n_values, const iArrayT& e_codes, dArray2DT& e_values)
{
    SimoQ1P0::ComputeOutput(n_codes, n_values, e_codes, e_values);
}

} // namespace Tahoe
