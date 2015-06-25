/* Id: FSMicromorphic3DT.cpp,v 1.6 2006/10/10 19:55:23 regueiro Exp $ */
#include "FSMicromorphic3DT.h"

#include "OutputSetT.h"
#include "ParameterContainerT.h"
#include "CommunicatorT.h"
#include <cmath>

using namespace Tahoe;

/* constructor */
FSMicromorphic3DT::FSMicromorphic3DT(const ElementSupportT& support):
    ElementBaseT(support), //pass the solid displacement field to the base class
    u(LocalArrayT::kDisp),
    u_dot(LocalArrayT::kVel),
    u_dot_n(LocalArrayT::kLastVel),
    u_dotdot(LocalArrayT::kAcc),
    u_dotdot_n(LocalArrayT::kLastAcc),
    u_n(LocalArrayT::kLastDisp),
    Phi(LocalArrayT::kDisp),
    Phi_dot(LocalArrayT::kVel),
    Phi_dot_n(LocalArrayT::kLastVel),
    Phi_dotdot(LocalArrayT::kAcc),
    Phi_dotdot_n(LocalArrayT::kLastAcc),
    Phi_n(LocalArrayT::kLastDisp),
    fInitCoords_displ(LocalArrayT::kInitCoords),
    fCurrCoords_displ(LocalArrayT::kCurrCoords),
    fInitCoords_micro(LocalArrayT::kInitCoords),
    fCurrCoords_micro(LocalArrayT::kCurrCoords),
    fTractionBCSet(0),
    fDispl(NULL),
    fMicro(NULL),
    fShapes_displ(NULL),
    fShapes_micro(NULL),
    fKdd(ElementMatrixT::kNonSymmetric),
    fKdphi(ElementMatrixT::kNonSymmetric),
    fKphid(ElementMatrixT::kNonSymmetric),
    fKphiphi(ElementMatrixT::kNonSymmetric),
    bStep_Complete(0)
{
    SetName("micromorphic_FS_3D");
}

/* destructor */
FSMicromorphic3DT::~FSMicromorphic3DT(void)
{
    delete fShapes_displ;
    delete fShapes_micro;
}


void FSMicromorphic3DT::Echo_Input_Data(void)
{

    cout << "#######################################################" << endl;
    cout << "############### ECHO FSMicromorphic3D DATA #########################" << endl;
    cout << "#######################################################" << endl;

    //################## material parameters ##################
    cout << "iConstitutiveModelType "               << iConstitutiveModelType         << endl;

    //-- Elasticity parameters for solid
    cout << "fMaterial_Params[kMu] "                << fMaterial_Params[kMu]          << endl;
    cout << "fMaterial_Params[kLambda] "            << fMaterial_Params[kLambda] << endl;
    cout << "fMaterial_Params[kNu] "                << fMaterial_Params[kNu]          << endl;
    cout << "fMaterial_Params[kSigma_const] "       << fMaterial_Params[kSigma_const]  << endl;
    cout << "fMaterial_Params[kTau] "               << fMaterial_Params[kTau]          << endl;
    cout << "fMaterial_Params[kEta] "               << fMaterial_Params[kEta]           << endl;
    cout << "fMaterial_Params[kKappa] "             << fMaterial_Params[kKappa]  << endl;
    cout << "fMaterial_Params[kTau1] "                  << fMaterial_Params[kTau1]          << endl;
    cout << "fMaterial_Params[kTau2] "                  << fMaterial_Params[kTau2]          << endl;
    cout << "fMaterial_Params[kTau3] "                  << fMaterial_Params[kTau3]          << endl;
    cout << "fMaterial_Params[kTau4] "                  << fMaterial_Params[kTau4]          << endl;
    cout << "fMaterial_Params[kTau5] "                  << fMaterial_Params[kTau5]          << endl;
    cout << "fMaterial_Params[kTau6] "                  << fMaterial_Params[kTau6]          << endl;
    cout << "fMaterial_Params[kTau7] "                  << fMaterial_Params[kTau7]          << endl;
    cout << "fMaterial_Params[kTau8] "                  << fMaterial_Params[kTau8]          << endl;
    cout << "fMaterial_Params[kTau9] "                  << fMaterial_Params[kTau9]          << endl;
    cout << "fMaterial_Params[kTau10] "                 << fMaterial_Params[kTau10]          << endl;
    cout << "fMaterial_Params[kTau11] "                 << fMaterial_Params[kTau11]          << endl;

}


//---------------------------------------------------------------------

void FSMicromorphic3DT::RHSDriver(void)
{

    int curr_group = ElementSupport().CurrentGroup();

    /* traction boundary conditions acting on displacement equations */
    if (curr_group == fDispl->Group())
    ApplyTractionBC();

    /* choose solution method */
    if (fDispl->Group() == fMicro->Group())
        RHSDriver_monolithic();
    else
        RHSDriver_staggered();

}
//---------------------------------------------------------------------

void FSMicromorphic3DT::Equations(AutoArrayT<const iArray2DT*>& eq_d,
                 AutoArrayT<const RaggedArray2DT<int>*>& eq_phi)
{
    /* doing monolithic solution */
    if (fDispl->Group() == fMicro->Group())
    {
        int ndof_micro = fMicro->NumDOF();
        int ndof_displ = fDispl->NumDOF();

        /* loop over connectivity blocks */
        fEqnos_displ.Dimension(fEqnos.Length());
        fEqnos_micro.Dimension(fEqnos.Length());
        for (int i = 0; i < fEqnos.Length(); i++)
        {
            /* connectivities */
            const iArray2DT& connects_displ = *(fConnectivities_displ[i]);
            const iArray2DT& connects_micro = *(fConnectivities_micro[i]);
            int nel = connects_displ.MajorDim();

            /* dimension */
            fEqnos[i].Dimension(nel, n_en_displ*ndof_displ + n_en_micro*ndof_micro);
            iArray2DT& displ_eq = fEqnos_displ[i];
            iArray2DT& micro_eq = fEqnos_micro[i];
            displ_eq.Dimension(nel, n_en_displ*ndof_displ);
            micro_eq.Dimension(nel, n_en_micro*ndof_micro);

            /* get equation numbers */
            fDispl->SetLocalEqnos(connects_displ, displ_eq);
            fMicro->SetLocalEqnos(connects_micro, micro_eq);

            /* write into one array */
            fEqnos[i].BlockColumnCopyAt(displ_eq, 0);
            fEqnos[i].BlockColumnCopyAt(micro_eq, displ_eq.MinorDim());

            /* add to list of equation numbers */
            eq_d.Append(&fEqnos[i]);
        }

        /* reset pointers to element cards */
        SetElementCards(fBlockData, fConnectivities_displ, fEqnos_displ, fElementCards_displ);
        SetElementCards(fBlockData, fConnectivities_micro, fEqnos_micro, fElementCards_micro);
    }
    else
    /* doing staggered */
    {
#pragma message("initialization for staggered solution needs to be corrected")

        /* ElementBaseT handles equation array for displacements */
        if (ElementSupport().CurrentGroup() == fDispl->Group())
            ElementBaseT::Equations(eq_d, eq_phi);

        /* micro-displacement-gradient equation */
        if (ElementSupport().CurrentGroup() == fMicro->Group())
        {
            /* collect local equation numbers */
            //fPress.SetLocalEqnos(fConnectivities_micro, fEqnos_micro);

            //eq_d.Append(&fEqnos_micro);
        }
    }

    /* get the equation number for the nodes on the faces */
    /*
    for (int i = 0; i < fMicroFaceEqnos.Length(); i++)
    {
        iArray2DT& faces = fMicroFaces[i];
        iArray2DT& eqnos = fMicroFaceEqnos[i];
        eqnos.Dimension(faces.MajorDim(), faces.MajorDim()*fDispl->NumDOF());

        fDispl->SetLocalEqnos(faces, eqnos);
    }
    */

}


//---------------------------------------------------------------------

void FSMicromorphic3DT::LHSDriver(GlobalT::SystemTypeT)
{
/** Everything done in RHSDriver for efficiency */
//cout << "############### In LHS Driver ############### \n";
}

//---------------------------------------------------------------------


void FSMicromorphic3DT::Select_Equations (const int &iBalLinChoice, const int &iBalFirstMomMomChoice )
{
    /** Choices for Linear Momentum Balance Equation */

    switch ( iBalLinChoice )    {

    default :
    cout << "FSMicromorphic3DT::Select_Equations() .. currently only one linear momentum balance for micromorphic continuum \n";
    break;
    }

    /** Choices for First Moment of Momentum Balance Equation */

    switch ( iBalFirstMomMomChoice )    {

    default :
    cout << "FSMicromorphic3DT::Select_Equations() .. currently only one first moment of momentum balance equation for micromorphic continuum \n";
    break;
    }

}

//---------------------------------------------------------------------

/* return true if the element contributes to the solution of the
 * given group. */
bool FSMicromorphic3DT::InGroup(int group) const
{
    return group == fDispl->Group() || group == fMicro->Group();
}

//---------------------------------------------------------------------


/* initialize/finalize step */
void FSMicromorphic3DT::InitStep(void)
{
    /* inherited */
    ElementBaseT::InitStep();
}


/* close current time increment */
void FSMicromorphic3DT::CloseStep(void)
{
    /* inherited */
    ElementBaseT::CloseStep();

    //-- Store/Register initial values in classic tahoe manner
    if ( ElementSupport().Time()==0 )
    {
    bStep_Complete=1;
    RHSDriver();
    bStep_Complete=0;

    /* my output set */
    const OutputSetT& output_set = ElementSupport().OutputSet(fOutputID);

    /* my nodes used */
    const iArrayT& nodes_used = output_set.NodesUsed();

    /* smooth stresses to nodes */
    ElementSupport().ResetAverage(knumstrain+knumstress+knum_d_state+knumdispl);
    dArray2DT out_variable_all;
    dArrayT out_variable;
    dArray2DT nd_var(NumElementNodes(), knumstrain+knumstress+knum_d_state+knumdispl);
    Top();
    while (NextElement())
    {
        /* extrapolate */
        nd_var = 0.0;
        out_variable_all.Alias(fNumIP_displ, knumstrain+knumstress+knum_d_state+knumdispl, fIPVariable(CurrElementNumber()));
        fShapes_displ->TopIP();
        while (fShapes_displ->NextIP())
        {
        out_variable.Alias(knumstrain+knumstress+knum_d_state+knumdispl, out_variable_all(fShapes_displ->CurrIP()));
        fShapes_displ->Extrapolate(out_variable, nd_var);
        }

        /* accumulate - extrapolation done from ip's to corners => X nodes  */
        ElementSupport().AssembleAverage(CurrentElement().NodesX(), nd_var);
    }

    /* get nodally averaged values */
    dArray2DT extrap_values;
    ElementSupport().OutputUsedAverage(extrap_values);

    /* temp space for group displacements */
    int num_node_output = fDispl->NumDOF() + fMicro->NumDOF() + knumstrain + knumstress + knum_d_state + knumdispl;
    dArray2DT n_values(nodes_used.Length(), num_node_output);

    /* collect nodal values */
    const dArray2DT& fPhi = (*fMicro)[0];
    const dArray2DT& fU = (*fDispl)[0];
    for (int i = 0; i < nodes_used.Length(); i++)
    {
        int node = nodes_used[i];
        double* row = n_values(i);
        for (int j = 0; j < fPhi.MinorDim(); j++)
        *row++ = fPhi(node,j);

        for (int j = 0; j < fU.MinorDim(); j++)
        *row++ = fU(node,j);

        double* p_stress = extrap_values(i);
        for (int j = 0; j < (knumstrain+knumstress+knum_d_state+knumdispl); j++)
        *row++ = p_stress[j];
    }

    /* send */
    ElementSupport().WriteOutput(fOutputID, n_values, fIPVariable);

    }

    //will not need this for quasi-static micromorphic
    /* zero first derivative of fields which are created at time=0 during calculating geostatic equilibrium (Trapezoidal rule) */
    /*
    if ( ElementSupport().Time()==0 &&  kInitialConditionType==1 && kAnalysisType==1)
    {
        FieldT* fmicro = const_cast <FieldT*> (fMicro);
        FieldT* fdispl = const_cast <FieldT*> (fDispl);

        (*fdispl)[1] = 0;
        (*fmicro)[1] = 0;
    }
    */

    //will not need this for quasi-static micromorphic
    /* zero second derivative of fields which are created at time=0 during calculating geostatic equilibrium(Newmark method) */
    /*
    if ( ElementSupport().Time()==0 &&  kInitialConditionType==1 && kAnalysisType==2)
    {
        FieldT* fmicro = const_cast <FieldT*> (fMicro);
        FieldT* fdispl = const_cast <FieldT*> (fDispl);

        (*fdispl)[2] = 0;
        (*fmicro)[2] = 0;
    }
    */

    /* reassign initial 2nd time derivative of micro-displacement-gradient to 1st derivative */
    /*
    if ( ElementSupport().Time()==0 &&  kInitialConditionType==2)
    {
        FieldT* fmicro = const_cast <FieldT*> (fMicro);
        (*fmicro)[1] = (*fmicro)[2];
        (*fmicro)[2] = 0;
    }
    */

    /* store more recently updated values */
    fdState = fdState_new;
    fiState = fiState_new;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////this part keeps track of the parameters from the previous iteration "n" to be used in "n+1"/////////////////////////

    SigN_IPs_el_n      = SigN_IPs_el;
    GammaN_IPs_el_n    = GammaN_IPs_el;
    mn_IPs_el_n        = mn_IPs_el;
    sn_sigman_IPs_el_n = sn_sigman_IPs_el;

    Fn_ar_IPs_el=F_ar_IPs_el;
    FnInv_ar_IPs_el=FInv_ar_IPs_el;
    ChiN_ar_IPs_el_n=Chi_ar_IPs_el;
    GRAD_ChiN_ar_IPs_el_n=GRAD_Chi_ar_IPs_el;
    //Counter_IPs_el_n=Counter_IPs_el;
//Here is close step function


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//    fs_micromorph3D_out   << endl
//      << setw(outputFileWidth) << "time_step"
//      << endl;
    step_number = ElementSupport().StepNumber();
//    fs_micromorph3D_out   << setw(outputFileWidth) << step_number
//      << endl;
//    fs_micromorph3D_out   << endl << "**********************************************************************************************";
//    fs_micromorph3D_out   << endl << "**********************************************************************************************" << endl;

/*    fs_micromorph3D_out        << "xxxxxxxxxxxxxx " << endl;
    fs_micromorph3D_out<<"step number"<<":"<<step_number<< endl;
    fs_micromorph3D_out        << "xxxxxxxxxxxxxx " << endl;
    for (int i=0; i<3; i++)
    {

        for (int j=0; j<3; j++)
        {
            fs_micromorph3D_out<<"Sigma"<< "("<<i<<","<<j<<")"<< " :  " ;
            fs_micromorph3D_out << Sigma(i,j) << "\t";
        }
        fs_micromorph3D_out        << endl ;
    }

*/
/*    fs_micromorph3D_out        << "xxxxxxxxxxxxxx " << endl;
    fs_micromorph3D_out        << "xxxxxxxxxxxxxx " << endl;*/
/*    for (int i=0; i<3; i++)
    {

        for (int j=0; j<3; j++)
        {
           for(int k=0;k<3;k++)
           {
            fs_micromorph3D_out<<"Mnplus1"<< "("<<i<<","<<j<<","<<k<<")"<< " :  " ;
            fs_micromorph3D_out << Mnplus1[i][j][k] <<"\t";
           }
        fs_micromorph3D_out << endl;
        }
       // fs_micromorph3DMn_out        << endl ;
    }*/
/*    for (int i=0; i<3; i++)
       {

           for (int j=0; j<3; j++)
           {
              for(int k=0;k<3;k++)
              {
               fs_micromorph3D_out<<"fMKLM"<< "("<<i<<","<<j<<","<<k<<")"<< " :  " ;
               fs_micromorph3D_out << fMKLM(i,j,k) <<"\t";
              }
           fs_micromorph3D_out << endl;
           }
       }
*/


}


/* resets to the last converged solution */
/*
GlobalT::RelaxCodeT FSMicromorphic3DT::ResetStep(void)
{
    const char caller[] = "FSMicromorphic3DT::ResetStep";

    // inherited
    GlobalT::RelaxCodeT relax = ElementBaseT::ResetStep();

    // update material internal variables
    //needs to be implemented
#pragma message("reseting internal variables not implemented")
    //ExceptionT::GeneralFail(caller, "reseting internal variables not implemented");

    return relax;
}
*/

/* element level reconfiguration for the current time increment */
/*
GlobalT::RelaxCodeT FSMicromorphic3DT::RelaxSystem(void)
{
    const char caller[] = "FSMicromorphic3DT::RelaxSystem";

    // inherited
    GlobalT::RelaxCodeT relax = ElementBaseT::RelaxSystem();

    // loop over materials
    //needs to be implemented
#pragma message("relax step for materials not implemented")
    //ExceptionT::GeneralFail(caller, "relax step for materials not implemented");

    return relax;
}
*/


void FSMicromorphic3DT::SendOutput(int kincode)
{
#pragma unused(kincode)
//not implemented
}


/* return geometry and number of nodes on each facet */
void FSMicromorphic3DT::FacetGeometry(ArrayT<GeometryT::CodeT>& facet_geometry,
    iArrayT& num_facet_nodes) const
{
    /* from integration domain */
    ShapeFunctionDispl().FacetGeometry(facet_geometry, num_facet_nodes);
}


/* form of tangent matrix */
GlobalT::SystemTypeT FSMicromorphic3DT::TangentType(void) const
{
    return GlobalT::kNonSymmetric;
}

/*
void FSMicromorphic3DT::SetStatus(const ArrayT<ElementCardT::StatusT>& status)
{
    // loop over elements and initial state variables
    int elem_num = 0;
    Top();
    while (NextElement())
    {
        // current element
        ElementCardT::StatusT& flag = CurrentElement().Flag();
        flag = status[elem_num++];

        if (flag == ElementCardT::kMarkON)
            flag = ElementCardT::kON;
        else if (flag == ElementCardT::kMarkOFF)
            flag = ElementCardT::kOFF;
    }
}
*/

/* initial condition/restart functions (per time sequence) */
void FSMicromorphic3DT::InitialCondition(void)
{
    /* inherited */
    ElementBaseT::InitialCondition();
}



//#############################################################################
//#############################################################################
//#############################################################################
//#############################################################################
//############################### NODAL FORCE  ################################
//#############################################################################
//#############################################################################
//#############################################################################
//#############################################################################

/* accumulate the residual force on the specified node */
void FSMicromorphic3DT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
    const char caller[] = "FSMicromorphic3DT::AddNodalForce";

    /* displ, micro, or neither */
    bool is_displ = false;
    dArrayT* element_force = NULL;
    int num_force = 0;
    if (field.FieldName() == fDispl->FieldName())
    {
        is_displ = true;
        element_force = &fFd_int;
        num_force = fDispl->NumDOF();
    }
    else if (field.FieldName() == fMicro->FieldName())
    {
        is_displ = false;
        element_force = &fFphi_int;
        num_force = fMicro->NumDOF();
    }
    else
    return;

    /* time Step Increment */
    double delta_t = ElementSupport().TimeStep();
    time = ElementSupport().Time();
    step_number = ElementSupport().StepNumber();

    /* temp for nodal force */
    dArrayT nodalforce;

    dArray2DT fdstatenew_all, fdstate_all;

    /* loop over elements */
    int e;
    Top();
    while (NextElement())
    {
        int nodeposition;
        const iArrayT& nodes_u = CurrentElement().NodesU();
        if (nodes_u.HasValue(node, nodeposition))
        {
            e = CurrElementNumber();
            const iArrayT& nodes_displ = fElementCards_displ[e].NodesU();
            const iArrayT& nodes_micro = fElementCards_micro[e].NodesU();

            u.SetLocal(nodes_displ);
            u_n.SetLocal(nodes_displ);
            Phi.SetLocal(nodes_micro);
            Phi_n.SetLocal(nodes_micro);

            del_u.DiffOf (u, u_n);
            del_Phi.DiffOf (Phi, Phi_n);

            // calculate derivatives based on reference coordinates
            fInitCoords_displ.SetLocal(fElementCards_displ[e].NodesX());
            //fCurrCoords_displ.SetToCombination (1.0, fInitCoords_displ, 1.0, u);
            fCurrCoords_displ=fInitCoords_displ;
            fShapes_displ->SetDerivatives_DN_DDN();

            //
            fInitCoords_micro.SetLocal(fElementCards_micro[e].NodesX());
            fCurrCoords_micro=fInitCoords_micro;
            //fCurrCoords_micro.SetToCombination (1.0, fInitCoords_micro, 1.0, u);
            fShapes_micro->SetDerivatives();

            //update state variables
            fdstatenew_all.Alias(fNumIP_displ, knum_d_state, fdState_new(CurrElementNumber()));
            fdstate_all.Alias(fNumIP_displ, knum_d_state, fdState(CurrElementNumber()));

            const double* Det    = fShapes_displ->IPDets();
            const double* Weight = fShapes_displ->IPWeights();
            /* calculate displacement nodal force */
            if (is_displ)
            {
                /* residual for displacement field */
                //generate this vector fFd_int
                fShapes_displ->TopIP();
                while (fShapes_displ->NextIP())
                {
                    //nothing right now
                    fFd_int=0.0;
                }
            }
            else /* pressure nodal force */
            {
                /* residual for micro-displacement-gradient field */
                // generate this vector fFphi_int
                fShapes_displ->TopIP();
                while (fShapes_displ->NextIP())
                {
                    //nothing right now
                    fFphi_int=0.0;
                }
            }

            /* loop over nodes (double-noding OK) */
            int dex = 0;
            if (is_displ)
            {
                for (int i = 0; i < nodes_displ.Length(); i++)
                {
                    if (nodes_displ[i] == node)
                    {
                        /* components for node */
                        nodalforce.Set(num_force, element_force->Pointer(dex));

                        /* accumulate */
                        force += nodalforce;
                    }
                    dex += fDispl->NumDOF();
                }
            }
            else /* micro-displacement-gradient nodal dof */
            {
                for (int i = 0; i < nodes_micro.Length(); i++)
                {
                    if (nodes_micro[i] == node)
                    {
                        /* components for node */
                        nodalforce.Set(num_force, element_force->Pointer(dex));

                        /* accumulate */
                        force += nodalforce;
                    }
                    dex += fMicro->NumDOF();
                }
            }
        }
    }
//  cout << "F_int = \n" << fFd_int << endl;
}

//---------------------------------------------------------------------

double FSMicromorphic3DT::InternalEnergy ( void )
{
//not implemented
    return 0.0;
}

//---------------------------------------------------------------------

/* write restart data to the output stream */
void FSMicromorphic3DT::WriteRestart(ostream& out) const
{
    /* inherited */
    ElementBaseT::WriteRestart(out);

    /* write state variable data */
    out << fdState;
}

//---------------------------------------------------------------------

/* read restart data to the output stream */
void FSMicromorphic3DT::ReadRestart(istream& in)
{
    /* inherited */
    ElementBaseT::ReadRestart(in);

    /* write state variable data */
    in >> fdState;
}

//---------------------------------------------------------------------

void FSMicromorphic3DT::RegisterOutput(void)
{
    /* collect block ID's */
    ArrayT<StringT> block_ID(fBlockData.Length());
    for (int i = 0; i < block_ID.Length(); i++)
    block_ID[i] = fBlockData[i].ID();

    /* output per element - strain, stress, and ISVs at the integration points */
    ArrayT<StringT> e_labels(fNumIP_displ*(knumstrain+knumstress+knumdispl+knum_d_state));

    /* over integration points */
    // enter what values you need at integration points
    // stress and strain
  //  const char* slabels3D[] = {"s11", "s22", "s33","s23","s13","s12","e11","e22","e33","e23","e13","e12"};
    const char* slabels3D[] = {"s11","s22","s33","s12","s13","s21","s23","s31","s32","e11","e22","e33","e12","e13","e21","e23","e31","e32","u1","u2","u3"};

    // state variables; ?
    const char* svlabels3D[] = {"thing1","thing2","J"};
    int count = 0;
    for (int j = 0; j < fNumIP_micro; j++)
    {
        StringT ip_label;
        ip_label.Append("ip", j+1);

        /* over strain and stress components */
        for (int i = 0; i < knumstrain+knumstress+knumdispl; i++)
        {
            e_labels[count].Clear();
            e_labels[count].Append(ip_label, ".", slabels3D[i]);
            count++;
        }

        /* over state variables */
        for (int i = 0; i < knum_d_state; i++)
        {
            e_labels[count].Clear();
            e_labels[count].Append(ip_label, ".", svlabels3D[i]);
            count++;
        }
    }

    /* output per node */
    int num_node_output = fDispl->NumDOF() + fMicro->NumDOF() + knumstrain + knumstress + knum_d_state+knumdispl;
    ArrayT<StringT> n_labels(num_node_output);
    count = 0;

    /* labels from micro-displacement-gradient */
    const ArrayT<StringT>& micro_labels = fMicro->Labels();
    for (int i = 0; i < micro_labels.Length(); i++)
    n_labels[count++] = micro_labels[i];

    /* labels from displacement */
    const ArrayT<StringT>& displ_labels = fDispl->Labels();
    for (int i = 0; i < displ_labels.Length(); i++)
    n_labels[count++] = displ_labels[i];

    /* labels from strains and stresses at the nodes */
    for (int i = 0; i < knumstrain+knumstress+knumdispl; i++)
    n_labels[count++] = slabels3D[i];

    /* labels from state variables at the nodes */
    for (int i = 0; i < knum_d_state; i++)
    n_labels[count++] = svlabels3D[i];

    /* set output specifier */
#pragma message("FSMicromorphic3DT::RegisterOutput: is this right? ")
    OutputSetT output_set(fGeometryCode_displ, block_ID, fConnectivities, n_labels, e_labels, false);

    /* register and get output ID */
    fOutputID = ElementSupport().RegisterOutput(output_set);
}

//#############################################################################
//#############################################################################
//#############################################################################
//#############################################################################
//############################### WRITE OUTPUT ################################
//#############################################################################
//#############################################################################
//#############################################################################
//#############################################################################

void FSMicromorphic3DT::WriteOutput(void)
{
    bStep_Complete=1;
    RHSDriver();
    bStep_Complete=0;

    /* my output set */
    const OutputSetT& output_set = ElementSupport().OutputSet(fOutputID);

    /* my nodes used */
    const iArrayT& nodes_used = output_set.NodesUsed();

    /* smooth stresses to nodes */
    ElementSupport().ResetAverage(knumstrain+knumstress+knum_d_state+knumdispl);
    dArray2DT out_variable_all;
    dArrayT out_variable;
    dArray2DT nd_var(NumElementNodes(), knumstrain+knumstress+knum_d_state+knumdispl);
    Top();
    while (NextElement())
    {
        /* extrapolate */
        nd_var = 0.0;
        out_variable_all.Alias(fNumIP_displ, knumstrain+knumstress+knum_d_state+knumdispl, fIPVariable(CurrElementNumber()));
        fShapes_displ->TopIP();
        while (fShapes_displ->NextIP())
        {
            out_variable.Alias(knumstrain+knumstress+knum_d_state+knumdispl, out_variable_all(fShapes_displ->CurrIP()));
            fShapes_displ->Extrapolate(out_variable, nd_var);
        }

        /* accumulate - extrapolation done from ip's to corners => X nodes  */
        ElementSupport().AssembleAverage(CurrentElement().NodesX(), nd_var);
    }

    /* get nodally averaged values */
    dArray2DT extrap_values;
    ElementSupport().OutputUsedAverage(extrap_values);

    /* temp space for group displacements */
    int num_node_output = fDispl->NumDOF() + fMicro->NumDOF() + knumstrain + knumstress + knum_d_state+knumdispl;
    dArray2DT n_values(nodes_used.Length(), num_node_output);

    /* collect nodal values */
    const dArray2DT& fPhi = (*fMicro)[0];
    const dArray2DT& fU = (*fDispl)[0];
    for (int i = 0; i < nodes_used.Length(); i++)
    {
        int node = nodes_used[i];
        double* row = n_values(i);
        for (int j = 0; j < fPhi.MinorDim(); j++)
            *row++ = fPhi(node,j);

        for (int j = 0; j < fU.MinorDim(); j++)
            *row++ = fU(node,j);

        double* p_stress = extrap_values(i);
        for (int j = 0; j < (knumstrain+knumstress+knum_d_state+knumdispl); j++)
            *row++ = p_stress[j];
    }

    /* send */
    ElementSupport().WriteOutput(fOutputID, n_values, fIPVariable);
}


//#############################################################################
//#############################################################################
//#############################################################################
//#############################################################################
//###################### Actual Solver Routines Below  ########################
//#############################################################################
//#############################################################################
//#############################################################################
//#############################################################################

/*************************************************************************
 * Private
 *************************************************************************/

/* form group contribution to the stiffness matrix and RHS */
void FSMicromorphic3DT::RHSDriver_staggered(void)
{
    const char caller[] = "FSMicromorphic3DT::RHSDriver_staggered";
#pragma message("staggered solution not implemented")
}

//---------------------------------------------------------------------
/* form group contribution to the stiffness matrix and RHS */
void FSMicromorphic3DT::RHSDriver_monolithic(void)
{
    const char caller[] = "FSMicromorphic3DT::RHSDriver_monolithic";
    if (fDispl->Group() != fMicro->Group())
    ExceptionT::GeneralFail(caller, "displacement and micro-displacement-gradient groups must be the same: %d != %d",
                fDispl->Group(), fMicro->Group());

    int curr_group = ElementSupport().CurrentGroup();

    /* stress output work space */
    dArray2DT   out_variable_all, fdstatenew_all, fdstate_all;
    dArrayT     out_variable;

    /* time Step Increment */
    double delta_t = ElementSupport().TimeStep();
    time = ElementSupport().Time();
    step_number = ElementSupport().StepNumber();

    /* print time */
//    fs_micromorph3D_out   <<"delta_t "<<delta_t << endl ;
//    fs_micromorph3D_out   <<"time "<<time << endl ;

    /* loop over elements */
    int e,l;
    Top();

//   fs_micromorph3D_out    <<"kInitialConditionType "<<kInitialConditionType << endl ;
//   fs_micromorph3D_out    <<"kAnalysisType "<<kAnalysisType << endl ;

    //don't need now for micromorphic
    /* at time=0 when geostatic initial condition is calculated,
       trapezoidal integrator will calculate first time derivative of fields
       which by setting alpha_delta_t = 1 will be changed to
       displacement and pressure which should be assigned to them,
       note that at time=0, delta_t=0 and Trapezoidal scheme
       which is embeded in the integrator will do nothing by itself
       (in changing previous iteration values)*/
    /*
    if ( ElementSupport().Time()==0 &&  kInitialConditionType==1 && kAnalysisType==1)
    {
        FieldT* fmicro = const_cast <FieldT*> (fMicro);
        FieldT* fdispl = const_cast <FieldT*> (fDispl);
        (*fmicro)[0] = (*fmicro)[1];
        (*fdispl)[0] = (*fdispl)[1];
    }
    */

    /* at time=0 when geostatic initial condition is calculated,
       dynamic Newmark integrator will calculate second time derivative
       of fields which by setting beta_delta_t2 = 1 will be changed to
       displacement and pressure which should be assigned to them,
       note that at time=0, delta_t=0 and Newmark scheme which is embeded in
       dynamic integrator will do nothing by itself(in changing previous iteration value)*/
    /*
    if ( ElementSupport().Time()==0 &&  kInitialConditionType==1 && kAnalysisType==2)
    {
        FieldT* fmicro = const_cast <FieldT*> (fMicro);
        FieldT* fdispl = const_cast <FieldT*> (fDispl);
        (*fmicro)[0] = (*fmicro)[2];
        (*fdispl)[0] = (*fdispl)[2];
    }
    */


    while (NextElement())
    {
    /* initialize */
/*    fFd_int_N1_vector = 0.0;
    fK_dd_G3_1_matrix = 0.0;
    fK_dd_G3_2_matrix = 0.0;
    fK_dd_G3_3_matrix = 0.0;
    fK_dd_G3_4_matrix = 0.0;
    fK_dd_BTDB_matrix = 0.0;
    fFd_int_smallstrain_vector = 0.0;
    fFd_int_G4_vector = 0.0;
    fK_dd_G4_matrix = 0.0;*/
    Sigma=0.0;
    Uint_1=0.0;
    Uint_1_temp=0.0;
    fFd_int=0.0;
    fG1_1=0.0;
    fG1_2=0.0;
    fG1_3=0.0;
    fG1_4=0.0;
    fG1_5=0.0;
    fG1_6=0.0;
    fG1_7=0.0;
    fG1_8=0.0;
    fG1_9=0.0;
    fG1_10=0.0;
    fG1_11=0.0;
    fG1_12=0.0;
    fG1_13=0.0;
    fG1_14=0.0;

    Pint_1=0.0;
    Pint_2=0.0;
    Pint_3=0.0;
    Pint_1_temp=0.0;
    Pint_2_temp=0.0;
    Pint_3_temp=0.0;
    fFphi_int=0.0;
    fH1_Etagrad=0.0;
    fH1_1=0.0;
    fH1_2=0.0;
    fH1_3=0.0;
    fH1_4=0.0;
    fH1_5=0.0;
    fH1_6=0.0;
   // fH1_7=0.0;
    fH1_71=0.0;
    fH1_72=0.0;
    fH1_73=0.0;
    fH1_74=0.0;
    fH1_75=0.0;
    fH1_76=0.0;
    fH1_77=0.0;
    fH1_78=0.0;
/*    fH1_8=0.0;
    fH1_9=0.0;
    fH1_10=0.0;*/
    fH1_11=0.0;
    fH1_12=0.0;
    fH1_13=0.0;
    fH1_14=0.0;

    fH2_1=0.0;
    fH2_2=0.0;
    fH2_3=0.0;
    fH2_4=0.0;
    fH2_5=0.0;
    fH2_6=0.0;
    fH2_7=0.0;
    fH2_8=0.0;
    fH2_9=0.0;
    fH2_10=0.0;
    fH2_11=0.0;
    fH2_12=0.0;
    fH2_13=0.0;

    fH3_1=0.0;

    fKdd=0.0;
    fKdphi=0.0;
    fKphid=0.0;
    fKphiphi=0.0;


    fIota_w_temp_matrix=0.0;
    fIota_eta_temp_matrix=0.0;
    NCHI=0.0;
    fTemp_matrix_nudof_x_nchidof=0.0;
    fTemp_matrix_nchidof_x_nchidof=0.0;
    fTemp_matrix_nchidof_x_nudof=0.0;


    ////////////////////////////////////////////////////////////////
    //////////////FINITE STRAIN MATRICES INITIALIZE/////////////////
    ////////////////////////////////////////////////////////////////
    //FSF=0;
    SPK=0.0;
     KirchhoffST=0.0;// The second Piola-Kirchhoff Matrix
     Temp_SPK=0.0;//temporary Matrix used in calculation of SPK
  //  dMatrixT FSF;
  //  dMatrixT LST;//Lagrangian strain tensor used in some functions to get rid of long name
     LagrangianStn=0.0;
     MicroStnTensor=0.0;//Micro-strain tensor
     PSI=0.0;//deformation measure PSI=Transpose(F).chi
     ChiM=0.0; //Micro-deformation tensor Chi ( used a different tensor this time )
     I1_1=0.0;
     I1_2=0.0;
     I1_3=0.0;
     I1_4=0.0;
     I1_5=0.0;
     I1_6=0.0;
     I1_7=0.0;
     I2_1=0.0;
     I1_8=0.0;
     I2_2=0.0;
     I1_9=0.0;
     I2_3=0.0;
     fFJ=0.0;
     fJF=0.0;
     fJ1_1=0.0;
     fJ1_2=0.0;
     fJ1_3=0.0;
     fJ1_4=0.0;
     fJ2_1=0.0;
     fJ1_5=0.0;
     fJ2_2=0.0;
     fJ1_6=0.0;
     fJ2_3=0.0;
     Vint_1=0.0;
     Vint_1_temp=0.0;
     Vint_2=0.0;
     Vint_2_temp=0.0;
     Vint_3=0.0;
     Vint_3_temp=0.0;

     fMKLM=0.0;
     GAMMA=0.0;
     GRAD_CHIM=0.0;
     fV1=0.0;
     fV2=0.0;
     fV3=0.0;
     fKu_1=0.0;
     fKu_2=0.0;
     fKu_3=0.0;
     fKu_4=0.0;
     fKu_5=0.0;
     fKu_6=0.0;
     fKu_7=0.0;
     fKFJu=0.0;
     fKJFu=0.0;



     fKuphi_1=0.0;
     fKu_8=0.0;
     fKuphi_2=0.0;
     fKu_9=0.0;
     fKuphi_3=0.0;
     fKphiu_1=0.0;
     fKphiu_2=0.0;
     fKphiu_3=0.0;
     fKphiu_4=0.0;
     fKphiu_5=0.0;
     fKphiphi_1=0.0;
     fKphiphi_2=0.0;
     fKphiu_6=0.0;
     fKphiphi_3=0.0;

     fFM=0.0;
     fMF=0.0;
     fEtaM=0.0;
     fMpu_1=0.0;
     fMpp_1=0.0;;
     fMpu_2=0.0;
     fMpp_2=0.0;

     fKMFphiu=0.0;;
     fKMchiphiphi=0.0;
     fKMphiu_1=0.0;;
     fKMphiphi_1=0.0;
     fKMphiu_2=0.0;
     fKMphiphi_2=0.0;


     Jmat=0.0;
     KJmat=0.0;
    ////////////////////////////////////////////////////////////////
    //////////////FINITE STRAIN MATRICES INITIALIZE/////////////////
    ////////////////////////////////////////////////////////////////

    e = CurrElementNumber();
    const iArrayT& nodes_displ = fElementCards_displ[e].NodesU();
    const iArrayT& nodes_micro = fElementCards_micro[e].NodesU();

    u.SetLocal(nodes_displ);
    u_n.SetLocal(nodes_displ);
    if (u_dot.IsRegistered()) u_dot.SetLocal(nodes_displ);
    if (u_dot_n.IsRegistered()) u_dot_n.SetLocal(nodes_displ);
    if (u_dotdot.IsRegistered()) u_dotdot.SetLocal(nodes_displ);
    if (u_dotdot_n.IsRegistered())u_dotdot_n.SetLocal(nodes_displ);

    Phi.SetLocal(nodes_micro);
    Phi_n.SetLocal(nodes_micro);
    if (Phi_dot.IsRegistered()) Phi_dot.SetLocal(nodes_micro);
    if (Phi_dot_n.IsRegistered()) Phi_dot_n.SetLocal(nodes_micro);
    if (Phi_dotdot.IsRegistered()) Phi_dotdot.SetLocal(nodes_micro);
    if (Phi_dotdot_n.IsRegistered()) Phi_dotdot_n.SetLocal(nodes_micro);

/*
 // print solid displacement at current step (u)
/*  fs_micromorph3D_out        <<"nodal solid displacement at current step(u)"<< endl ;
    for (int i=0; i<n_en_displ; i++)
    {
        fs_micromorph3D_out        << "node number " << i+1 <<" :  " ;
        for (int j=0; j<n_sd; j++)
        fs_micromorph3D_out << u(i,j) << "\t";
        fs_micromorph3D_out        << endl ;
    }*/



    /* print micro-displacement-gradient at current step (micro)*/
/*  fs_micromorph3D_out        <<"nodal micro-displacement-gradient at current step(micro)"<< endl ;
    for (int i=0; i<n_en_micro; i++)
    {
        fs_micromorph3D_out        << "node number " << i+1 <<" :  " ;
        fs_micromorph3D_out        << Phi(i,0) << endl;
    }
*/

/*    for (int i=0; i<3; i++)
    {
       for(int j=0;j<3;j++)
       {
    	   for(int k=0;k<3;k++)
    	   {
    	    	fs_micromorph3D_out        << "fMKLM(i,j,k) " << i<<j<<k <<" :  " ;
    	        fs_micromorph3D_out        << fMKLM(i,j,k) << endl;
    	   }
       }

    }
*/
    /* populate solid displacement,solid velocity and
       solid accelration in vector form*/


    int index_u = 0;
    for (int i=0; i<n_en_displ; i++)
    {
        for (int j=0; j<n_sd; j++)
        {
            u_vec[index_u] = u(i,j);
            u_dot_vec[index_u] = u_dot(i,j);
            u_dotdot_vec[index_u] = u_dotdot(i,j);
            index_u += 1;
        }
    }

    /* [u_dot_column_matrix] will be formed */
    for (int i=0; i<n_en_displ_x_n_sd; i++)
        u_dot_column_matrix(i,0) = u_dot_vec[i];

    /* [u_dot_column_matrix_Transpose] will be formed */
    u_dot_column_matrix_Transpose.Transpose(u_dot_column_matrix);

    /* [u_dotdot_column_matrix] will be formed */
    for (int i=0; i<n_en_displ_x_n_sd; i++)
        u_dotdot_column_matrix(i,0) = u_dotdot_vec[i];

    /* populate micro-displacement-gradient, first and second time derivatives of
       micro-displacement-gradient in vector form*/
    int index_Phi = 0;
    for (int i=0; i<n_en_micro; i++)
    {
        for (int j=0; j<ndof_per_nd_micro; j++)
        {
            Phi_vec[index_Phi] = Phi(i,j);
            Phi_dot_vec[index_Phi] = Phi_dot(i,j);
            Phi_dotdot_vec[index_Phi] = Phi_dotdot(i,j);
            index_Phi += 1;
        }
    }

    /* [micro_dot_column_matrix] will be formed */
    /*
    for (int i=0; i<n_en_micro_x_ndof_per_nd_micro; i++)
        Phi_dot_column_matrix(i,0) = Phi_dot_vec[i];
    */

    del_u.DiffOf (u, u_n);
    del_Phi.DiffOf (Phi, Phi_n);

    // calculate derivatives based on reference coordinates
    fInitCoords_displ.SetLocal(fElementCards_displ[e].NodesX());
    fCurrCoords_displ=fInitCoords_displ;
    //fCurrCoords_displ.SetToCombination (1.0, fInitCoords_displ, 1.0, u);
    fShapes_displ->SetDerivatives_DN_DDN();
    //
    fInitCoords_micro.SetLocal(fElementCards_micro[e].NodesX());
    fCurrCoords_micro=fInitCoords_micro;
    //fCurrCoords_micro.SetToCombination (1.0, fInitCoords_micro, 1.0, u);
    fShapes_micro->SetDerivatives();

    //update state variables
    fdstatenew_all.Alias(fNumIP_displ, knum_d_state, fdState_new(CurrElementNumber()));
    fdstate_all.Alias(fNumIP_displ, knum_d_state, fdState(CurrElementNumber()));
    /*
    fdstatenew_all.Alias(fNumIP_micro, knum_d_state, fdState_new(CurrElementNumber()));
    fdstate_all.Alias(fNumIP_micro, knum_d_state, fdState(CurrElementNumber()));
    */

    if (bStep_Complete)
    {
        //-- Store/Register data in classic tahoe manner
        //out_variable_all.Alias(fNumIP_displ, knumstrain+knumstress+knum_d_state, fIPVariable(CurrElementNumber()));
        out_variable_all.Alias(fNumIP_micro, knumstrain+knumstress+knum_d_state+knumdispl, fIPVariable(CurrElementNumber()));
        for (l=0; l < fNumIP_displ; l++)
        //for (l=0; l < fNumIP_micro; l++)
        {
            out_variable.Alias(knumstrain+knumstress+knum_d_state+knumdispl, out_variable_all(l));

          Put_values_In_dArrayT_vector(fCauchy_stress_Elements_IPs, e,l,fTemp_nine_values);
//            Put_values_In_dArrayT_vector(fCauchy_stress_Elements_IPs, e,l,fTemp_six_values);

          out_variable.CopyIn(0,fTemp_nine_values);
//            out_variable.CopyIn(0,fTemp_six_values);

//        Put_values_In_dArrayT_vector(fEulerian_strain_Elements_IPs, e,l,fTemp_six_values);
          Put_values_In_dArrayT_vector(fEulerian_strain_Elements_IPs, e,l,fTemp_nine_values);

//        out_variable.CopyIn(6,fTemp_six_values);
          out_variable.CopyIn(9,fTemp_nine_values);//!!9->6?
            /*
            out_variable[13]=fState_variables_Elements_IPs(e,l*3+0);
            out_variable[14]=fState_variables_Elements_IPs(e,l*3+1);
            out_variable[15]=fState_variables_Elements_IPs(e,l*3+2);
            */

          /* Not sure if these are correct!!! */
      //    Put_values_In_dArrayT_vector(fDisplacement_Element_IPs,e,l,u_element);
          Put_values_In_Array(fDisplacement_Element_IPs,e,l,ftemp_u_element);
          out_variable.CopyIn(18,ftemp_u_element);


        }

/*
        int e;
        e=CurrElementNumber();
       cout<<e<<endl;
        el_num=0;
        if(e==8 || e== 9 || e==10 || e==11)
        {
         fs_micromorph3D_out<<"element number="<<e<< " :  " <<"\n";
         fShapeDispl.Multx(u_vec,u_element);
      	 for(int ii=0;ii<3;ii++)
      	 {
      		 u_el(el_num,ii)=u_element[ii];
      	 }


      	  for(int ii=0;ii<3;ii++)
      	  {
             fs_micromorph3D_out << u_el(el_num,ii) <<"\n";
      	  }
      	  el_num++;
        }


*/

    }
    else
    { //-- Still Iterating

            /* residual and tangent for displacements */

            const double* Det    = fShapes_displ->IPDets();
            const double* Weight = fShapes_displ->IPWeights();
            fShapes_displ->TopIP();
            fShapes_micro->TopIP();

            SigN_IPs_el_n.RowCopy(e,SigN_IPs_n);
            mn_IPs_el_n.RowCopy(e,mn_IPs_n);
            GammaN_IPs_el_n.RowCopy(e,GammaN_IPs_n);
            sn_sigman_IPs_el_n.RowCopy(e,sn_sigman_IPs_n);

            Fn_ar_IPs_el.RowCopy(e,Fn_ar_IPs);
            FnInv_ar_IPs_el.RowCopy(e,FnInv_ar_IPs);
            ChiN_ar_IPs_el_n.RowCopy(e,ChiN_ar_IPs_n);
            GRAD_ChiN_ar_IPs_el_n.RowCopy(e,GRAD_ChiN_ar_IPs);


            while (fShapes_displ->NextIP() && fShapes_micro->NextIP())
            {
                double scale_const = (*Weight++)*(*Det++);

                   const int IP = fShapes_displ->CurrIP();
                   dArrayT DisplIPCoordinate(n_sd), MicroIPCoordinate(n_sd);
                   fShapes_displ->IPCoords(DisplIPCoordinate);
                   fShapes_micro->IPCoords(MicroIPCoordinate);

                   const double* shapes_displ_X = fShapes_displ->IPShapeX();
                   /* [fShapeDispl] will be formed */
                   Form_solid_shape_functions(shapes_displ_X);//output:fShapeDispl

               //   fShapeDispl_Tr.Transpose(fShapeDispl);

                   fShapes_displ->GradNa(fShapeDisplGrad_temp);
                   /* [fShapeDisplGrad] will be formed *///in reference config
                   Form_Gradient_of_solid_shape_functions(fShapeDisplGrad_temp);//output:fShapeDisplGrad in Reference config.
                   Form_GRAD_Nuw_matrix(fShapeDisplGrad_temp) ;//output:GRAD_Nuw

                   const double* shapes_micro_X = fShapes_micro->IPShapeX();
                  //  {fShapeMicro} will be formed
                   Form_micro_shape_functions(shapes_micro_X);//output:fShapeMicro

                    //[fShapeMicro_row_matrix] will be formed
                   //need?
                   for (int i=0; i<n_en_micro ; i++)
                       fShapeMicro_row_matrix(0,i) = fShapeMicro[i];

                  //  [fShapeMicroGrad] will be formed
                   fShapes_micro->GradNa(fShapeMicroGrad_temp);

                   //the correct name should be NPHI, NCHI is not a proper name!
                   Form_NCHI_matrix(fShapeMicro_row_matrix); //output: NCHI matrix
                   NCHI_Tr.Transpose(NCHI);
                   Form_Gradient_of_micro_shape_eta_functions(fShapeMicroGrad_temp);//output: GRAD_NCHI

                   /* [fIdentity_matrix] will be formed */
                   fIdentity_matrix = 0.0;
                   for (int i=0; i<n_sd ; i++)
                       fIdentity_matrix(i,i) =1.0;

                   /* [fDeformation_Gradient] will be formed */
                   Form_deformation_gradient_tensor();//output: F (deform. grad. tensor)
                   /* [fDeformation_Gradient_Inverse] and [fDeformation_Gradient_Transpose] and [fDeformation_Gradient_Inverse_Transpose] will be formed */
                   if (fDeformation_Gradient.Det()==0)
                       fDeformation_Gradient = fIdentity_matrix;
                   fDeformation_Gradient_Inverse.Inverse(fDeformation_Gradient);
                  // Form Chi[i][j] and  dMatrixT ChiM(i,j)
                  Form_micro_deformation_tensor_Chi();//output: Chi[i][j]
                  //Form_ChiM();//It is also micro-deformation gradient tensor but defined as dMatrixT
                  Form_Chi_inv_matrix();//output: ChiInv

                   SigN_IPs_n.RowCopy(IP,SigN_ar);
                   sn_sigman_IPs_n.RowCopy(IP,sn_sigman);
                   GammaN_IPs_n.RowCopy(IP,GammaN_ar);
                   mn_IPs_n.RowCopy(IP,mn_ar);
                   Mapping_double_and_Array(-1);

                   Fn_ar_IPs.RowCopy(IP,Fn_ar);
                   FnInv_ar_IPs.RowCopy(IP,FnInv_ar);
                   ChiN_ar_IPs_n.RowCopy(IP,ChiN_ar);
                   GRAD_ChiN_ar_IPs.RowCopy(IP,GRAD_ChiN_ar);
                   Form_deformation_tensors_arrays(-1);


                   /* KroneckerDelta matrix is formed*/
                   Form_KroneckerDelta_matrix();//output: KrDelta
                   Form_CCof_tensor();//output: Coeff tensor


                   Form_double_Finv_from_Deformation_tensor_inverse();// output: Finv
                   Form_GRAD_Chi_matrix();////CHI=1+PHI ==> GRAD_CHI=GRAD_PHI output: GRAD_Chi[i][J][K] AND GRAD_CHIM which is the dTensor3DT form
                   Form_Gamma_tensor3D();

                   Form_Finv_w_matrix();//output: Finv_w
                   Form_Finv_eta_matrix();//output: Finv_eta



                   /* [fDefGradInv_Grad_grad] will be formed */
                   Form_Grad_grad_transformation_matrix();//output:fDefGradInv_Grad_grad


                   /* Calculating Jacobian */
                   double J = fDeformation_Gradient.Det();
                   /* Jacobian for the current IP will be saved */
                   fState_variables_IPs(IP,2)=J;

                   /*  fRho */
                   fRho_0 = fMaterial_Params[kRho_0];

                   /* Calculating fRho */
                   fRho = fRho_0/J;

                   /* [fRight_Cauchy_Green_tensor] will be formed */
                   fRight_Cauchy_Green_tensor.MultATB(fDeformation_Gradient, fDeformation_Gradient);
                   /* [fRight_Cauchy_Green_tensor_Inverse] will be formed */
                   if (fRight_Cauchy_Green_tensor.Det()==0)
                       fRight_Cauchy_Green_tensor = fIdentity_matrix;
                   fRight_Cauchy_Green_tensor_Inverse.Inverse(fRight_Cauchy_Green_tensor);

                   LagrangianStn=fIdentity_matrix;
                   LagrangianStn*=-1;
                   LagrangianStn+=fRight_Cauchy_Green_tensor;
                   LagrangianStn*=0.5;

                   // Micro-Strain tensor will be formed
                    MicroStnTensor  = fIdentity_matrix;
                    MicroStnTensor *= -1;
                    PSI.MultATB(fDeformation_Gradient,ChiM);
                    MicroStnTensor += PSI;


                    //GAMMA deformation measure will be formed
                //   GAMMA.ContractIndex(GRAD_CHIM,0,fDeformation_Gradient,1);

    /*               fs_micromorph3D_out<<"MicroStnTensor"<< endl ;
                    for (int i=0; i<3; i++)
                    {
                    	for(int j=0;j<3;j++)
                    	{
                    		fs_micromorph3D_out<< "MicroStnTensor(i,j)"<<i<<j<<endl;
                    		fs_micromorph3D_out<< MicroStnTensor(i,j)<<endl;
                    	}
                    }*/

                   /* [fLeft_Cauchy_Green_tensor] will be formed */
                   fLeft_Cauchy_Green_tensor.MultABT(fDeformation_Gradient, fDeformation_Gradient);
                   /* [fLeft_Cauchy_Green_tensor_Inverse] will be formed */
                   if (fLeft_Cauchy_Green_tensor.Det()==0)
                       fLeft_Cauchy_Green_tensor = fIdentity_matrix;
                   fLeft_Cauchy_Green_tensor_Inverse.Inverse(fLeft_Cauchy_Green_tensor);

                   /* [fEulerian_strain_tensor_current_IP] will be formed */
                   fEulerian_strain_tensor_current_IP = fLeft_Cauchy_Green_tensor_Inverse;
                   fEulerian_strain_tensor_current_IP *= -1;
                   fEulerian_strain_tensor_current_IP += fIdentity_matrix;
                   fEulerian_strain_tensor_current_IP *= 0.5;
               //    fEulerian_strain_tensor_current_IP = LagrangianStn;
                   Extract_six_values_from_symmetric_tensor(fEulerian_strain_tensor_current_IP,fTemp_nine_values);
 //                Extract_six_values_from_symmetric_tensor(fEulerian_strain_tensor_current_IP,fTemp_six_values);

                   /* Save Eulerian strain tensor of the current IP */
                   fEulerian_strain_IPs.SetRow(IP,fTemp_nine_values);
  //                 fEulerian_strain_IPs.SetRow(IP,fTemp_six_values);
                   fShapeDispl.Multx(u_vec,u_element);
                   fDisplacement_IPs.SetRow(IP,u_element);

                   /* [fIota_temp_matrix] will be formed */
                   fIota_temp_matrix.MultATB(fShapeDisplGrad,fDefGradInv_Grad_grad);

                   /* [fIota_w_temp_matrix] will be formed*/
                   fIota_w_temp_matrix.MultATBT(GRAD_Nuw,Finv_w);
                   /* [fIota_eta_temp_matrix] will be formed*/
                   fIota_eta_temp_matrix.MultATBT(GRAD_NCHI,Finv_eta);


                   //fShapeDisplGrad--> [GRAD(Ns,e)] so it in reference configuration


                   /* second derivatives of solid shape functions, [fShapeDisplGradGrad] will be formed */
                   //fShapes_displ->Grad_GradNa(fShapeDisplGradGrad);


                   double scale = scale_const;

                   if(iConstitutiveModelType==1)
                   {

                	 //  double invJ=1/J;

                       Form_Second_Piola_Kirchhoff_SPK();
                       KirchhoffST.MultABCT(fDeformation_Gradient,SPK,fDeformation_Gradient);
                       Form_fV1();
                      // fIota_temp_matrix.Multx(fV1,Vint_1_temp);
                       fShapeDisplGrad.MultTx(fV1,Vint_1_temp);
                     // fIota_w_temp_matrix.Multx(fV1,Vint_1_temp);
                       scale=scale_const;
                       Vint_1_temp*=scale;
                       Vint_1 +=Vint_1_temp;

                       // -eta(ml)s_sigma(ml)-eta(ml,k)m(klm)
                       // eta(ml)s_sima(ml)+eta(ml,k)m(klm)
                       Form_SIGMA_S();//in current configuration SIGMA_S=s_sigma, but what we use sigma_s, so it needs to be multiplied by "-1"
                       Form_fV2();//gives F.SIGMA_S.F^T = s_sigma
                       NCHI.MultTx(fV2,Vint_2_temp);
                       scale=scale_const;
                       Vint_2_temp*=scale;
                       Vint_2 +=Vint_2_temp;


                       Form_GAMMA();
					   Form_fMKLM();
					   Form_fV3();
					   //fIota_eta_temp_matrix.Multx(fV3,Vint_3_temp);
					   GRAD_NCHI.MultTx(fV3,Vint_3_temp);
					   scale=scale_const;
					   Vint_3_temp*=scale;
					   Vint_3+=Vint_3_temp;




                       //Sigma.SetToScaled(1/J,KirchhoffST);
                       //Sigma*=1.7;
                       Sigma=KirchhoffST;
                       Sigma.SetToScaled(1/J,KirchhoffST);
                     //  cout<< invJ<<endl;
                      fCauchy_stress_tensor_current_IP=Sigma;

                    //fCauchy_stress_tensor_current_IP=SPK;
                    //Extract_six_values_from_symmetric_tensor(fCauchy_stress_tensor_current_IP,fTemp_six_values);
                      Extract_six_values_from_symmetric_tensor(fCauchy_stress_tensor_current_IP,fTemp_nine_values);

                      // Save Cauchy effective stress tensor of the current IP
                    //fCauchy_stress_IPs.SetRow(IP,fTemp_six_values);
                      fCauchy_stress_IPs.SetRow(IP,fTemp_nine_values);


                       /*internal force is calculated from BLM */
/*                       Form_I1_1();
                       Form_I1_2();*/

                    //   Form_Jmat();

                       Form_I1_3();
                       Form_I1_4();
                       Form_I1_5();
                       Form_I1_6();
                       Form_I1_7();
                       Form_I2_1();
                       Form_I1_8();
                       Form_I2_2();
                       Form_I1_9();
                       Form_I2_3();
                       Form_fFJ();
                       Form_fJF();
                       Form_fJ1_1();
                       Form_fJ1_2();
                       Form_fJ1_3();
                       Form_fJ1_4();
                       Form_fJ2_1();
                       Form_fJ1_5();
                       Form_fJ2_2();
                       Form_fJ1_6();
                       Form_fJ2_3();

                     //  Form_fFM();
                       Form_fMF();
                       Form_fMchi();
                       Form_fMpu_1();
                       Form_fMpp_1();
                       Form_fMpu_2();
                       Form_fMpp_2();


/*
                       fTemp_matrix_nudof_x_nudof.MultABCT(fIota_temp_matrix,I1_1,fIota_temp_matrix);
                       scale =-scale_const;
                       fTemp_matrix_nudof_x_nudof *= scale;
                       // accumulate
                       fKu_1 += fTemp_matrix_nudof_x_nudof;

                       fTemp_matrix_nudof_x_nudof.MultABCT(fIota_temp_matrix,I1_2,fIota_temp_matrix);
                       scale = scale_const;
                       fTemp_matrix_nudof_x_nudof *= scale;
                       // accumulate
                       fKu_2 += fTemp_matrix_nudof_x_nudof;
*/

/*
                       fTemp_matrix_nudof_x_nudof.MultABCT(fIota_temp_matrix,I1_3,fIota_temp_matrix);
                       scale = scale_const;
                       fTemp_matrix_nudof_x_nudof *= scale;
                       // accumulate
                       fKu_3 += fTemp_matrix_nudof_x_nudof;
                       */
                       fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I1_3,fShapeDisplGrad);
                       scale = scale_const;
                       fTemp_matrix_nudof_x_nudof *= scale;
                       fKu_3 += fTemp_matrix_nudof_x_nudof;


                       //Matrices from variation of SPK
/*                       fTemp_matrix_nudof_x_nudof.MultABCT(fIota_temp_matrix,I1_4,fIota_temp_matrix);
                       scale = scale_const*(fMaterial_Params[kLambda]+fMaterial_Params[kTau]);
                       fTemp_matrix_nudof_x_nudof *= scale;
                       // accumulate
                       fKu_4 += fTemp_matrix_nudof_x_nudof;
*/
                       fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I1_4,fShapeDisplGrad);
                       scale = scale_const*(fMaterial_Params[kLambda]+fMaterial_Params[kTau]);
                       fTemp_matrix_nudof_x_nudof *= scale;
                       fKu_4 += fTemp_matrix_nudof_x_nudof;



/*
                       fTemp_matrix_nudof_x_nudof.MultABCT(fIota_temp_matrix,I1_5,fIota_temp_matrix);
                       scale =scale_const*(fMaterial_Params[kMu]+fMaterial_Params[kSigma_const]);
                       fTemp_matrix_nudof_x_nudof *= scale;
                       // accumulate
                       fKu_5 += fTemp_matrix_nudof_x_nudof;
*/


                       fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I1_5,fShapeDisplGrad);
                       scale = scale_const*(fMaterial_Params[kMu]+fMaterial_Params[kSigma_const]);
                       fTemp_matrix_nudof_x_nudof *= scale;
                       fKu_5 += fTemp_matrix_nudof_x_nudof;



/*
                       fTemp_matrix_nudof_x_nudof.MultABCT(fIota_temp_matrix,I1_6,fIota_temp_matrix);
                       scale = scale_const*(fMaterial_Params[kMu]+fMaterial_Params[kSigma_const]);
                       fTemp_matrix_nudof_x_nudof *= scale;
                       // accumulate
                       fKu_6 += fTemp_matrix_nudof_x_nudof;
*/


                       fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I1_6,fShapeDisplGrad);
                       scale = scale_const*(fMaterial_Params[kMu]+fMaterial_Params[kSigma_const]);
                       fTemp_matrix_nudof_x_nudof *= scale;
                       fKu_6 += fTemp_matrix_nudof_x_nudof;

                       //to be deleted
/*                     //  fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,Jmat,fShapeDisplGrad);
                       //scale = J*scale_const;
                       //fTemp_matrix_nudof_x_nudof *= scale;
                       //KJmat+=fTemp_matrix_nudof_x_nudof;*/
                       ///


                       fTemp_matrix_nudof_x_nudof.MultABCT(fIota_temp_matrix,I1_7,fIota_temp_matrix);
                       scale = scale_const*fMaterial_Params[kEta];
                       fTemp_matrix_nudof_x_nudof *= scale;
                       // accumulate
                       fKu_7 += fTemp_matrix_nudof_x_nudof;

                       fTemp_matrix_nudof_x_nchidof.MultABC(fIota_temp_matrix,I2_1,NCHI);//ABC not ABCT
                       scale = scale_const*fMaterial_Params[kEta];
                       fTemp_matrix_nudof_x_nchidof *= scale;
                       // accumulate
                       fKuphi_1 += fTemp_matrix_nudof_x_nchidof;

                       fTemp_matrix_nudof_x_nudof.MultABCT(fIota_temp_matrix,I1_8,fIota_temp_matrix);
                       scale = scale_const*fMaterial_Params[kKappa];
                       fTemp_matrix_nudof_x_nudof *= scale;
                       // accumulate
                       fKu_8 += fTemp_matrix_nudof_x_nudof;



                       fTemp_matrix_nudof_x_nchidof.MultABC(fIota_temp_matrix,I2_2,NCHI);//ABC not ABCT
                       scale = scale_const*fMaterial_Params[kKappa];
                       fTemp_matrix_nudof_x_nchidof *= scale;
                       // accumulate
                       fKuphi_2 += fTemp_matrix_nudof_x_nchidof;

                       fTemp_matrix_nudof_x_nudof.MultABCT(fIota_temp_matrix,I1_9,fIota_temp_matrix);
                       scale = scale_const*fMaterial_Params[kNu];
                       fTemp_matrix_nudof_x_nudof *= scale;
                       // accumulate
                       fKu_9 += fTemp_matrix_nudof_x_nudof;

                       fTemp_matrix_nudof_x_nchidof.MultABC(fIota_temp_matrix,I2_3,NCHI);//ABC not ABCT
                       scale = scale_const*fMaterial_Params[kNu];
                       fTemp_matrix_nudof_x_nchidof *= scale;
                       // accumulate
                       fKuphi_3 += fTemp_matrix_nudof_x_nchidof;

/*******************************************************************************************************/



                       //fTemp_matrix_nchidof_x_nudof.MultABCT(NCHI_Tr,fFJ,fIota_temp_matrix);
                       fTemp_matrix_nchidof_x_nudof.MultATBC(NCHI,fFJ,fShapeDisplGrad);
                       scale =scale_const;
                       fTemp_matrix_nchidof_x_nudof *= scale;
                       // accumulate
                       fKFJu += fTemp_matrix_nchidof_x_nudof;


                      // fTemp_matrix_nchidof_x_nudof.MultABCT(NCHI_Tr,fJF,fIota_temp_matrix);
                       fTemp_matrix_nchidof_x_nudof.MultATBC(NCHI,fJF,fShapeDisplGrad);
                       scale =scale_const;
                       fTemp_matrix_nchidof_x_nudof *= scale;
                       // accumulate
                       fKJFu += fTemp_matrix_nchidof_x_nudof;



                      // fTemp_matrix_nchidof_x_nudof.MultABCT(NCHI_Tr,fJ1_1 ,fIota_temp_matrix);
                       fTemp_matrix_nchidof_x_nudof.MultATBC(NCHI,fJ1_1,fShapeDisplGrad);
                       scale =scale_const*fMaterial_Params[kTau];
                       fTemp_matrix_nchidof_x_nudof *= scale;
                       // accumulate
                       fKphiu_1 += fTemp_matrix_nchidof_x_nudof;

                       fTemp_matrix_nchidof_x_nudof.MultABCT(NCHI_Tr,fJ1_2 ,fIota_temp_matrix);
                       scale =scale_const*fMaterial_Params[kSigma_const];
                       fTemp_matrix_nchidof_x_nudof *= scale;
                       // accumulate
                       fKphiu_2 += fTemp_matrix_nchidof_x_nudof;

                       fTemp_matrix_nchidof_x_nudof.MultABCT(NCHI_Tr,fJ1_3 ,fIota_temp_matrix);
                       scale =scale_const*fMaterial_Params[kSigma_const];
                       fTemp_matrix_nchidof_x_nudof *= scale;
                       // accumulate
                       fKphiu_3 += fTemp_matrix_nchidof_x_nudof;

                      // fTemp_matrix_nchidof_x_nudof.MultABCT(NCHI_Tr,fJ1_4,fIota_temp_matrix);
                       fTemp_matrix_nchidof_x_nudof.MultATBC(NCHI,fJ1_4,fShapeDisplGrad);
                       scale =scale_const*(fMaterial_Params[kEta]-fMaterial_Params[kTau]);
                       fTemp_matrix_nchidof_x_nudof *= scale;
                       // accumulate
                       fKphiu_4 += fTemp_matrix_nchidof_x_nudof;


                       fTemp_matrix_nchidof_x_nchidof.MultABC(NCHI_Tr,fJ2_1,NCHI);
                       scale =scale_const*(fMaterial_Params[kEta]-fMaterial_Params[kTau]);
                       fTemp_matrix_nchidof_x_nchidof *= scale;
                       // accumulate
                       fKphiphi_1 += fTemp_matrix_nchidof_x_nchidof;

                       fTemp_matrix_nchidof_x_nudof.MultABCT(NCHI_Tr,fJ1_5,fIota_temp_matrix);
                       scale =scale_const*(fMaterial_Params[kNu]-fMaterial_Params[kSigma_const]);
                       fTemp_matrix_nchidof_x_nudof *= scale;
                       // accumulate
                       fKphiu_5 += fTemp_matrix_nchidof_x_nudof;

                       fTemp_matrix_nchidof_x_nchidof.MultABC(NCHI_Tr,fJ2_2,NCHI);
                       scale =scale_const*(fMaterial_Params[kNu]-fMaterial_Params[kSigma_const]);
                       fTemp_matrix_nchidof_x_nchidof *= scale;
                       // accumulate
                       fKphiphi_2 += fTemp_matrix_nchidof_x_nchidof;

                       fTemp_matrix_nchidof_x_nudof.MultABCT(NCHI_Tr,fJ1_6,fIota_temp_matrix);
                       scale =scale_const*(fMaterial_Params[kKappa]-fMaterial_Params[kSigma_const]);
                       fTemp_matrix_nchidof_x_nudof *= scale;
                       // accumulate
                       fKphiu_5 += fTemp_matrix_nchidof_x_nudof;

                       fTemp_matrix_nchidof_x_nchidof.MultABC(NCHI_Tr,fJ2_3,NCHI);
                       scale =scale_const*(fMaterial_Params[kKappa]-fMaterial_Params[kSigma_const]);
                       fTemp_matrix_nchidof_x_nchidof *= scale;
                       // accumulate
                       fKphiphi_3 += fTemp_matrix_nchidof_x_nchidof;
/*************************************************************************************************************/

                      // fTemp_matrix_nchidof_x_nudof.MultABCT(fIota_eta_temp_matrix,fMF,fIota_temp_matrix);
                       fTemp_matrix_nchidof_x_nudof.MultATBC(GRAD_NCHI,fMF,fShapeDisplGrad);
                       scale =scale_const;
                       fTemp_matrix_nchidof_x_nudof *= scale;
                       // accumulate
                       fKMFphiu += fTemp_matrix_nchidof_x_nudof;

                       //fTemp_matrix_nchidof_x_nchidof.MultABC(fIota_eta_temp_matrix,fMchi,NCHI);
                       fTemp_matrix_nchidof_x_nchidof.MultATBC(GRAD_NCHI,fMchi,NCHI);
                       scale =scale_const;
                       fTemp_matrix_nchidof_x_nchidof *= scale;
                       // accumulate
                       fKMchiphiphi += fTemp_matrix_nchidof_x_nchidof;


                       //fTemp_matrix_nchidof_x_nudof.MultABCT(fIota_eta_temp_matrix,fMpu_1,fIota_temp_matrix);
                       fTemp_matrix_nchidof_x_nudof.MultATBC(GRAD_NCHI,fMpu_1,fShapeDisplGrad);
                       scale =scale_const*fMaterial_Params[kTau8];
                       fTemp_matrix_nchidof_x_nudof *= scale;
                       // accumulate
                       fKMphiu_1 += fTemp_matrix_nchidof_x_nudof;

                       //fTemp_matrix_nchidof_x_nchidof.MultABC(fIota_eta_temp_matrix,fMpp_1,GRAD_NCHI);
                       fTemp_matrix_nchidof_x_nchidof.MultATBC(GRAD_NCHI,fMpp_1,GRAD_NCHI);
                       scale =scale_const*fMaterial_Params[kTau8];
                       fTemp_matrix_nchidof_x_nchidof *= scale;
                       // accumulate
                       fKMphiphi_1 += fTemp_matrix_nchidof_x_nchidof;

                       //fTemp_matrix_nchidof_x_nudof.MultABCT(fIota_eta_temp_matrix,fMpu_2,fIota_temp_matrix);
                       fTemp_matrix_nchidof_x_nudof.MultATBC(GRAD_NCHI,fMpu_2,fShapeDisplGrad);
                       scale =scale_const*fMaterial_Params[kTau8];
                       fTemp_matrix_nchidof_x_nudof *= scale;
                       // accumulate
                       fKMphiu_2 += fTemp_matrix_nchidof_x_nudof;

                       //fTemp_matrix_nchidof_x_nchidof.MultABC(fIota_eta_temp_matrix,fMpp_2,GRAD_NCHI);
                       fTemp_matrix_nchidof_x_nchidof.MultATBC(GRAD_NCHI,fMpp_2,GRAD_NCHI);
                       scale =scale_const*fMaterial_Params[kTau8];
                       fTemp_matrix_nchidof_x_nchidof *= scale;
                       // accumulate
                       fKMphiphi_2 += fTemp_matrix_nchidof_x_nchidof;

                   }
                   else
                   {
   ////////////////////////////////////////////////////////////////////////////////////////
   /////////////////MicroMorphic Internal force vectors////////////////////////////////////
                   Form_G1_matrix();//output:G1 vector & Sigma matrix
                   fIota_w_temp_matrix.Multx(G1,Uint_1_temp);
                   scale=-1*scale_const*J;
                   Uint_1_temp*=scale;
                   Uint_1 +=Uint_1_temp;
   //              fShapeDispl.MultTx(fGravity_vector,Uext_1);
   //              Uext_1*=-fMaterial_Params[kRho_0];


                   Form_H1_matrix();//output:H1 vctor & Mnplus1 tensor
                   fIota_eta_temp_matrix.Multx(H1,Pint_1_temp);
                   scale=-1*scale_const*J;
                   Pint_1_temp*=scale;
                   Pint_1+=Pint_1_temp;

                   Form_H2_matrix();//output: H2 vector & s_sigma matrix
                   NCHI.MultTx(H2,Pint_2_temp);
                   scale=scale_const*J;
                   Pint_2_temp*=scale;//the sign is taken into account when forming H2 vector.
                   Pint_2+=Pint_2_temp;

                   Form_H3_matrix();//output H3 vector
                   NCHI.MultTx(H3,Pint_3_temp);
                   scale=scale_const*fMaterial_Params[kRho_0];
                   Pint_3_temp*=scale;
                   Pint_3+=Pint_3_temp;

                   SPK.MultABCT(fDeformation_Gradient_Inverse,Sigma,fDeformation_Gradient_Inverse);
                   SPK*=J;
                  // fCauchy_stress_tensor_current_IP=SPK;
                   fCauchy_stress_tensor_current_IP=Sigma;

                   // extract six values of stress from symmetric cauchy stress tensor
                  Extract_six_values_from_symmetric_tensor(fCauchy_stress_tensor_current_IP,fTemp_nine_values);
               // Extract_six_values_from_symmetric_tensor(fCauchy_stress_tensor_current_IP,fTemp_six_values);

                //Save Cauchy effective stress tensor of the current IP
                //fCauchy_stress_IPs.SetRow(IP,fTemp_six_values);
                  fCauchy_stress_IPs.SetRow(IP,fTemp_nine_values);

   ////////////////////////////////////////////////////////////////////////////////
   /////////////////MicroMorphic Internal force vectors finish here////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////
/*                   for(int iii=0;iii<3;iii++)
                   {
                       for(int jjj=0;jjj<3;jjj++)
                       {
                           fDeformation_Gradient(iii,jjj)=fIdentity_matrix(iii,jjj);
                           Chi[iii][jjj]=fIdentity_matrix(iii,jjj);
                           ChiN[iii][jjj]=fIdentity_matrix(iii,jjj);
                           Finv[iii][jjj]=fIdentity_matrix(iii,jjj);
                           Fn[iii][jjj]=fIdentity_matrix(iii,jjj);
                           for(int kkk=0;kkk<3;kkk++)
                           {
                               GRAD_Chi[iii][jjj][kkk]=0.0;
                               GRAD_ChiN[iii][jjj][kkk]=0.0;
                           }
                       }
                   }*/

   ////////////////Micromorphic 3-D Matrices are being formed coming from linearization process//////////////////////////
                   Form_CapitalLambda_matrix();//output:CapitalLambda
                   Form_Var_F_tensor();
                   Form_Tsigma_1_matrix();
                   Form_Tsigma_2_matrix();
                   Form_Tsigma_3_matrix();
                   Form_TFn_1_matrix();
                   Form_TFn_2_matrix();
                   Form_TFn_3_matrix();
                   Form_TFn_4_matrix();
                   Form_TFn_5_matrix();
                   Form_TFn_6_matrix();
                   Form_TChi_1_matrix();
                   Form_TChi_2_matrix();
                   Form_TChi_3_matrix();
                   Form_SigCurr_matrix();
   //////////////////////////////////////////////////////////
                   Form_Etagrad_matrix();
                   Form_Mm_1_matrix();
                   Form_Mm_2_matrix();
                   Form_Mm_3_matrix();
                   Form_Mm_4_matrix();
                   Form_Mm_5_matrix();
                   Form_Mm_6_matrix();
                   Form_Mm_7_matrix();
                   Form_Mm_71_matrix();
                   Form_Mm_72_matrix();
                   Form_Mm_73_matrix();
                   Form_Mm_74_matrix();
                   Form_Mm_75_matrix();
                   Form_Mm_76_matrix();
                   Form_Mm_77_matrix();
                   Form_Mm_78_matrix();
                   Form_Mm_8_matrix();
                   Form_Mm_9_matrix();
                   Form_Mm_10_matrix();
                   Form_Mm_11_matrix();
                   Form_Mm_12_matrix();
                   Form_Mm_13_matrix();
                   Form_Mm_14_matrix();
                   Form_Ru_1_matrix();
                   Form_Ru_2_matrix();
                   Form_Ru_3_matrix();
                   Form_RChi_1_matrix();
                   Form_Ru_4_matrix();
                   Form_RChi_2_matrix();
                   Form_Ru_5_matrix();
                   Form_Ru_6_matrix();
                   Form_Ru_7_matrix();
                   Form_RChi_3_matrix();
                   Form_Ru_8_matrix();
                   Form_Ru_9_matrix();

                   Form_Rs_sigma_matrix();//output:Rs_sigma
                   Form_R_Capital_Lambda_Chi_matrix();//output:R_Capital_Gamma_Chi


   ////////////////////////Finished here///////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////


   ////////////////////////////////////////////////////////////////////////////////
   ////////////////fG1_ matrices are constructed////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////

                   fTemp_matrix_nudof_x_nudof.MultABCT(fIota_w_temp_matrix,Tsigma_1,fIota_temp_matrix);
                   scale = scale_const*J;
                   fTemp_matrix_nudof_x_nudof *= scale;
                   // accumulate
                   fG1_1 += fTemp_matrix_nudof_x_nudof;

                   fTemp_matrix_nudof_x_nudof.MultABCT(fIota_w_temp_matrix,Tsigma_2,fIota_temp_matrix);
                   scale =-scale_const*J;
                   fTemp_matrix_nudof_x_nudof *= scale;
                   // accumulate
                   fG1_2 += fTemp_matrix_nudof_x_nudof;

                   fTemp_matrix_nudof_x_nudof.MultABCT(fIota_w_temp_matrix,Tsigma_3,fIota_temp_matrix);
                   scale = -scale_const*J;
                   fTemp_matrix_nudof_x_nudof *= scale;
                   // accumulate
                   fG1_3 += fTemp_matrix_nudof_x_nudof;

                   fTemp_matrix_nudof_x_nudof.MultABCT(fIota_temp_matrix,TFn_1,fIota_temp_matrix);
                   scale = -scale_const*J*(fMaterial_Params[kLambda]+fMaterial_Params[kTau]);
                   fTemp_matrix_nudof_x_nudof *= scale;
                   // accumulate
                   fG1_4 += fTemp_matrix_nudof_x_nudof;


                   fTemp_matrix_nudof_x_nudof.MultABCT(fIota_w_temp_matrix,TFn_2,fIota_temp_matrix);
                   scale = -scale_const*J*(fMaterial_Params[kMu]+fMaterial_Params[kSigma_const]);
                   fTemp_matrix_nudof_x_nudof *= scale;
                    //accumulate
                   fG1_5 += fTemp_matrix_nudof_x_nudof;


                   fTemp_matrix_nudof_x_nudof.MultABCT(fIota_w_temp_matrix,TFn_3,fIota_temp_matrix);
                   scale = -scale_const*J*(fMaterial_Params[kMu]+fMaterial_Params[kSigma_const]);
                   fTemp_matrix_nudof_x_nudof *= scale;
                   // accumulate
                   fG1_6 += fTemp_matrix_nudof_x_nudof;


                   fTemp_matrix_nudof_x_nchidof.MultABC(fIota_w_temp_matrix,TChi_1,NCHI);// ABC not ABCT
                   scale = -scale_const*J*(fMaterial_Params[kEta]);
                   fTemp_matrix_nudof_x_nchidof *= scale;
                   // accumulate
                   fG1_7 += fTemp_matrix_nudof_x_nchidof;


                   fTemp_matrix_nudof_x_nudof.MultABCT(fIota_w_temp_matrix,TFn_4,fIota_temp_matrix);
                   scale = -scale_const*J*(fMaterial_Params[kEta]);
                   fTemp_matrix_nudof_x_nudof *= scale;
                   // accumulate
                   fG1_8 += fTemp_matrix_nudof_x_nudof;


                   fTemp_matrix_nudof_x_nchidof.MultABC(fIota_w_temp_matrix,TChi_2,NCHI);//ABC not ABCT
                   scale = -scale_const*J*(fMaterial_Params[kKappa]);
                   fTemp_matrix_nudof_x_nchidof *= scale;
                   // accumulate
                   fG1_9 += fTemp_matrix_nudof_x_nchidof;

                   fTemp_matrix_nudof_x_nudof.MultABCT(fIota_w_temp_matrix,TFn_5,fIota_temp_matrix);
                   scale = -scale_const*J*(fMaterial_Params[kKappa]);
                   fTemp_matrix_nudof_x_nudof *= scale;
                   // accumulate
                   fG1_10 += fTemp_matrix_nudof_x_nudof;



                   fTemp_matrix_nudof_x_nchidof.MultABC(fIota_w_temp_matrix,TChi_3,NCHI);//ABC not ABCT
                   scale = -scale_const*J*(fMaterial_Params[kNu]);
                   fTemp_matrix_nudof_x_nchidof *= scale;
                   // accumulate
                   fG1_11 += fTemp_matrix_nudof_x_nchidof;

                   fTemp_matrix_nudof_x_nudof.MultABCT(fIota_w_temp_matrix,TFn_6,fIota_temp_matrix);
                   scale = -scale_const*J*(fMaterial_Params[kNu]);
                   fTemp_matrix_nudof_x_nudof *= scale;
                   // accumulate
                   fG1_12 += fTemp_matrix_nudof_x_nudof;

                  fTemp_matrix_nudof_x_nudof.MultABC(fIota_w_temp_matrix,SigCurr,fShapeDisplGrad);//ABC not ABCT
                  //fTemp_matrix_nudof_x_nudof.MultABCT(fIota_w_temp_matrix,SigCurr,fIota_temp_matrix);//ABC
                   scale = -scale_const*J;
                   fTemp_matrix_nudof_x_nudof *= scale;
                   // accumulate
                   fG1_13 += fTemp_matrix_nudof_x_nudof;

                   //TransShapeDisplGrad.Transpose(GRAD_Nuw);//
                   //fTemp_matrix_nudof_x_nudof.MultABCT(TransShapeDisplGrad,Var_F,fIota_temp_matrix);
                   fTemp_matrix_nudof_x_nudof.MultABCT(fIota_w_temp_matrix,Var_F,fIota_temp_matrix);
                   scale= scale_const*J;
                   fTemp_matrix_nudof_x_nudof *= scale;
                   fG1_14 +=fTemp_matrix_nudof_x_nudof;


   ////////////////////////////////////////////////////////////////////////////////
   /////////////////fG1_ matrices finish here//////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////////
   /////////////////fH_ matrices are constructed///////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////
                   fTemp_matrix_nchidof_x_nudof.MultABCT(fIota_eta_temp_matrix,Etagrad,fIota_temp_matrix);
                   scale =scale_const*J;
                   fTemp_matrix_nchidof_x_nudof *= scale;
                   // accumulate
                   fH1_Etagrad += fTemp_matrix_nchidof_x_nudof;



                   fTemp_matrix_nchidof_x_nudof.MultABCT(fIota_eta_temp_matrix,Mm_1,fIota_temp_matrix);
                   scale = scale_const*J;
                   fTemp_matrix_nchidof_x_nudof *= scale;
                   // accumulate
                   fH1_1 += fTemp_matrix_nchidof_x_nudof;


                   fTemp_matrix_nchidof_x_nudof.MultABCT(fIota_eta_temp_matrix,Mm_2,fIota_temp_matrix);
                   scale = -1*scale_const*J;
                   fTemp_matrix_nchidof_x_nudof *= scale;
                   // accumulate
                   fH1_2 += fTemp_matrix_nchidof_x_nudof;

                   fTemp_matrix_nchidof_x_nudof.MultABCT(fIota_eta_temp_matrix,Mm_3,fIota_temp_matrix);
                   scale = -1*scale_const*J;
                   fTemp_matrix_nchidof_x_nudof *= scale;
                   // accumulate
                   fH1_3 += fTemp_matrix_nchidof_x_nudof;


                   fTemp_matrix_nchidof_x_nchidof.MultABC(fIota_eta_temp_matrix,Mm_4,NCHI);
                   scale =-1*scale_const*J;
                   fTemp_matrix_nchidof_x_nchidof *= scale;
                   // accumulate
                   fH1_4 += fTemp_matrix_nchidof_x_nchidof;

                   fTemp_matrix_nchidof_x_nchidof.MultABC(fIota_eta_temp_matrix,Mm_5,NCHI);
                   scale = -1*scale_const*J;
                   fTemp_matrix_nchidof_x_nchidof *= scale;
                   // accumulate
                   fH1_5 += fTemp_matrix_nchidof_x_nchidof;


                   fTemp_matrix_nchidof_x_nchidof.MultABC(fIota_eta_temp_matrix,Mm_6,NCHI);
                   scale = scale_const*J;
                   fTemp_matrix_nchidof_x_nchidof *= scale;
                   // accumulate
                   fH1_6 += fTemp_matrix_nchidof_x_nchidof;

                   fTemp_matrix_nchidof_x_nchidof.MultABC(fIota_eta_temp_matrix,Mm_7,NCHI);
                   scale =-1*scale_const*J;
                   fTemp_matrix_nchidof_x_nchidof *= scale;
                   // accumulate
                   fH1_7 += fTemp_matrix_nchidof_x_nchidof;

                   fTemp_matrix_nchidof_x_nchidof.MultABC(fIota_eta_temp_matrix,Mm_71,GRAD_NCHI);
                   scale =-1*scale_const*J;
                   fTemp_matrix_nchidof_x_nchidof *= scale;
                   // accumulate
                   fH1_71 += fTemp_matrix_nchidof_x_nchidof;

                   fTemp_matrix_nchidof_x_nudof.MultABCT(fIota_eta_temp_matrix,Mm_72,fIota_temp_matrix);
                   scale =1*scale_const*J;
                   fTemp_matrix_nchidof_x_nudof *= scale;
                   // accumulate
                   fH1_72 += fTemp_matrix_nchidof_x_nudof;

                   fTemp_matrix_nchidof_x_nchidof.MultABC(fIota_eta_temp_matrix,Mm_73,NCHI);
                   scale =-1*scale_const*J;
                   fTemp_matrix_nchidof_x_nchidof *= scale;
                   // accumulate
                   fH1_73 += fTemp_matrix_nchidof_x_nchidof;

                   fTemp_matrix_nchidof_x_nchidof.MultABC(fIota_eta_temp_matrix,Mm_74,NCHI);
                   scale =-1*scale_const*J;
                   fTemp_matrix_nchidof_x_nchidof *= scale;
                   // accumulate
                   fH1_74 += fTemp_matrix_nchidof_x_nchidof;

                   fTemp_matrix_nchidof_x_nchidof.MultABC(fIota_eta_temp_matrix,Mm_75,NCHI);
                   scale =1*scale_const*J;
                   fTemp_matrix_nchidof_x_nchidof *= scale;
                   // accumulate
                   fH1_75 += fTemp_matrix_nchidof_x_nchidof;

                   fTemp_matrix_nchidof_x_nchidof.MultABC(fIota_eta_temp_matrix,Mm_76,GRAD_NCHI);
                   scale =-1*scale_const*J;
                   fTemp_matrix_nchidof_x_nchidof *= scale;
                   // accumulate
                   fH1_76 += fTemp_matrix_nchidof_x_nchidof;

                   fTemp_matrix_nchidof_x_nchidof.MultABC(fIota_eta_temp_matrix,Mm_77,NCHI);
                   scale =1*scale_const*J;
                   fTemp_matrix_nchidof_x_nchidof *= scale;
                   // accumulate
                   fH1_77 += fTemp_matrix_nchidof_x_nchidof;

                   fTemp_matrix_nchidof_x_nudof.MultABCT(fIota_eta_temp_matrix,Mm_78,fIota_temp_matrix);
                   scale =1*scale_const*J;
                   fTemp_matrix_nchidof_x_nudof *= scale;
                   // accumulate
                   fH1_78 += fTemp_matrix_nchidof_x_nudof;


                   fTemp_matrix_nchidof_x_nchidof.MultABC(fIota_eta_temp_matrix,Mm_8,NCHI);
                   scale =scale_const*J;
                   fTemp_matrix_nchidof_x_nchidof *= scale;
                   // accumulate
                   fH1_8 += fTemp_matrix_nchidof_x_nchidof;

                   fTemp_matrix_nchidof_x_nchidof.MultABC(fIota_eta_temp_matrix,Mm_9,GRAD_NCHI);
                   scale = -1*scale_const*J;
                   fTemp_matrix_nchidof_x_nchidof *= scale;
                   // accumulate
                   fH1_9 += fTemp_matrix_nchidof_x_nchidof;


                   fTemp_matrix_nchidof_x_nchidof.MultABC(fIota_eta_temp_matrix,Mm_10,NCHI);
                   scale = -1*scale_const*J;
                   fTemp_matrix_nchidof_x_nchidof *= scale;
                   // accumulate
                   fH1_10 += fTemp_matrix_nchidof_x_nchidof;



                   fTemp_matrix_nchidof_x_nudof.MultABCT(fIota_eta_temp_matrix,Mm_11,fIota_temp_matrix);
                   scale = -1*scale_const*J;
                   fTemp_matrix_nchidof_x_nudof *= scale;
                   // accumulate
                   fH1_11 += fTemp_matrix_nchidof_x_nudof;



                   fTemp_matrix_nchidof_x_nudof.MultABCT(fIota_eta_temp_matrix,Mm_12,fIota_temp_matrix);
                   scale = -1*scale_const*J;
                   fTemp_matrix_nchidof_x_nudof *= scale;
                   // accumulate
                   fH1_12 += fTemp_matrix_nchidof_x_nudof;


                   fTemp_matrix_nchidof_x_nchidof.MultABC(fIota_eta_temp_matrix,Mm_13,NCHI);
                   scale = -1*scale_const*J;
                   fTemp_matrix_nchidof_x_nchidof *= scale;
                   // accumulate
                   fH1_13 += fTemp_matrix_nchidof_x_nchidof;


                   fTemp_matrix_nchidof_x_nudof.MultABC(fIota_eta_temp_matrix,Mm_14,fShapeDisplGrad);
                   scale = -1*scale_const*J;
                   fTemp_matrix_nchidof_x_nudof *= scale;
                   // accumulate
                   fH1_14 += fTemp_matrix_nchidof_x_nudof;

                   fTemp_matrix_nchidof_x_nudof.MultABCT(NCHI_Tr,Ru_1,fIota_temp_matrix);
                   scale = scale_const*J;
                   fTemp_matrix_nchidof_x_nudof *= scale;
                   // accumulate
                   fH2_1 += fTemp_matrix_nchidof_x_nudof;

                   fTemp_matrix_nchidof_x_nudof.MultABCT(NCHI_Tr,Ru_2,fIota_temp_matrix);
                   scale = -scale_const*J;
                   fTemp_matrix_nchidof_x_nudof *= scale;
                   // accumulate
                   fH2_2 += fTemp_matrix_nchidof_x_nudof;

                   fTemp_matrix_nchidof_x_nudof.MultABCT(NCHI_Tr,Ru_3,fIota_temp_matrix);
                   scale = -scale_const*J;
                   fTemp_matrix_nchidof_x_nudof *= scale;
                   // accumulate
                   fH2_3 += fTemp_matrix_nchidof_x_nudof;

                   fTemp_matrix_nchidof_x_nchidof.MultABC(NCHI_Tr,RChi_1,NCHI);
                   scale = -scale_const*J*(fMaterial_Params[kKappa]-fMaterial_Params[kSigma_const]);
                   fTemp_matrix_nchidof_x_nudof *= scale;
                   // accumulate
                   fH2_4 += fTemp_matrix_nchidof_x_nchidof;

                   fTemp_matrix_nchidof_x_nudof.MultABCT(NCHI_Tr,Ru_4,fIota_temp_matrix);
                   scale = -scale_const*J*(fMaterial_Params[kKappa]-fMaterial_Params[kSigma_const]);
                   fTemp_matrix_nchidof_x_nudof *= scale;
                   // accumulate
                   fH2_5 += fTemp_matrix_nchidof_x_nudof;

                   fTemp_matrix_nchidof_x_nchidof.MultABC(NCHI_Tr,RChi_2,NCHI);
                   scale = -scale_const*J*(fMaterial_Params[kNu]-fMaterial_Params[kSigma_const]);
                   fTemp_matrix_nchidof_x_nudof *= scale;
                   // accumulate
                   fH2_6 += fTemp_matrix_nchidof_x_nchidof;

                   fTemp_matrix_nchidof_x_nudof.MultABCT(NCHI_Tr,Ru_5,fIota_temp_matrix);
                   scale = -scale_const*J*(fMaterial_Params[kNu]-fMaterial_Params[kSigma_const]);
                   fTemp_matrix_nchidof_x_nudof *= scale;
                   // accumulate
                   fH2_7 += fTemp_matrix_nchidof_x_nudof;


                   fTemp_matrix_nchidof_x_nudof.MultABCT(NCHI_Tr,Ru_6,fIota_temp_matrix);
                   scale = -1*scale_const*J*fMaterial_Params[kTau];
                   fTemp_matrix_nchidof_x_nudof *= scale;
                   // accumulate
                   fH2_8 += fTemp_matrix_nchidof_x_nudof;


                   fTemp_matrix_nchidof_x_nudof.MultABCT(NCHI_Tr,Ru_7,fIota_temp_matrix);
                   scale = -1*scale_const*J*fMaterial_Params[kSigma_const];
                   fTemp_matrix_nchidof_x_nudof *= scale;
                   // accumulate
                   fH2_9 += fTemp_matrix_nchidof_x_nudof;

                   fTemp_matrix_nchidof_x_nudof.MultABCT(NCHI_Tr,Ru_8,fIota_temp_matrix);
                   scale = -1*scale_const*J*fMaterial_Params[kSigma_const];
                   fTemp_matrix_nchidof_x_nudof *= scale;
                   // accumulate
                   fH2_10 += fTemp_matrix_nchidof_x_nudof;

                   fTemp_matrix_nchidof_x_nchidof.MultABC(NCHI_Tr,RChi_3,NCHI);
                   scale = -scale_const*J*(fMaterial_Params[kEta]-fMaterial_Params[kTau]);
                   fTemp_matrix_nchidof_x_nchidof *= scale;
                   // accumulate
                   fH2_11 += fTemp_matrix_nchidof_x_nchidof;

                   fTemp_matrix_nchidof_x_nudof.MultABCT(NCHI_Tr,Ru_9,fIota_temp_matrix);
                   scale = -scale_const*J*(fMaterial_Params[kEta]-fMaterial_Params[kTau]);
                   fTemp_matrix_nchidof_x_nudof *= scale;
                   // accumulate
                   fH2_12 += fTemp_matrix_nchidof_x_nudof;

                   fTemp_matrix_nchidof_x_nudof.MultABC(NCHI_Tr,Rs_sigma,fShapeDisplGrad);
                   scale = -1*scale_const*J;
                   fTemp_matrix_nchidof_x_nudof *= scale;
                   // accumulate
                   fH2_13 += fTemp_matrix_nchidof_x_nudof;

                   fTemp_matrix_nchidof_x_nchidof.MultABC(NCHI_Tr,R_Capital_Gamma_Chi,NCHI);
                   scale = 1*scale_const*fMaterial_Params[kRho_0];//*J?????
                   fTemp_matrix_nchidof_x_nudof *= scale;
                   // accumulate
                   fH3_1 += fTemp_matrix_nchidof_x_nchidof;


   /////////////////fH_ matrices finish here////////////////////////////////////////

                    //[fK_dd_G3_2_matrix] will be formed
/*
                   fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fIota_temp_matrix,fHbar_temp_matrix,fIota_temp_matrix);
                   scale = fMaterial_Params[kMu] * scale_const;
                   fTemp_matrix_ndof_se_x_ndof_se *= scale;
                   //accumulate
                   fK_dd_G3_2_matrix += fTemp_matrix_ndof_se_x_ndof_se;


                   //[fK_dd_G3_3_matrix] will be formed
                   fTemp_matrix_ndof_se_x_ndof_se.MultABCT(fIota_temp_matrix,fEll_temp_matrix,fIota_temp_matrix);
                   scale = fMaterial_Params[kMu] * scale_const;
                   fTemp_matrix_ndof_se_x_ndof_se *= scale;
                   // accumulate
                   fK_dd_G3_3_matrix += fTemp_matrix_ndof_se_x_ndof_se;


                   // [fK_dd_G3_4_matrix] will be formed
                   fTemp_matrix_ndof_se_x_ndof_se.MultABC(fIota_temp_matrix,fI_ij_column_matrix,fPi_temp_row_matrix);
                   scale = fMaterial_Params[kLambda] * scale_const;
                   fTemp_matrix_ndof_se_x_ndof_se *= scale;
                   // accumulate
                   fK_dd_G3_4_matrix += fTemp_matrix_ndof_se_x_ndof_se;

*/

                   //need?

                   // {fFd_int_G4_vector} will be formed
                   fShapeDispl.MultTx(fGravity_vector,fTemp_vector_ndof_se);
                   scale = -1*fRho_0*scale_const;
                   fTemp_vector_ndof_se *= scale;
                   // accumulate
                   fFd_int_G4_vector += fTemp_vector_ndof_se;


                   // [fK_dd_G4_matrix] will be formed
                   //pg57 Davoud's thesis

                   fTemp_matrix_ndof_se_x_ndof_se.MultATBC(fShapeDispl,fGravity_column_matrix,fPi_temp_row_matrix);
                   scale = -1*J*(fRho)*scale_const;
                   fTemp_matrix_ndof_se_x_ndof_se *= scale;
                   // accumulate
                   fK_dd_G4_matrix += fTemp_matrix_ndof_se_x_ndof_se;



                   /////////////////////saving matrices at Gauss Points////////////

                   Form_deformation_gradient_tensor();
                   Form_micro_deformation_tensor_Chi();
                   Form_GRAD_Chi_matrix();

                   Mapping_double_and_Array(1);
                   GammaN_IPs.SetRow(IP,GammaN_ar);
                   SigN_IPs.SetRow(IP,Sigma);
                   sn_sigman_IPs.SetRow(IP,s_sigma);
                   mn_IPs.SetRow(IP,mn_ar);

                   Form_deformation_tensors_arrays(1);
                   F_ar_IPs.SetRow(IP,fDeformation_Gradient);
                   FInv_ar_IPs.SetRow(IP,fDeformation_Gradient_Inverse);
                   Chi_ar_IPs.SetRow(IP,Chi_ar);
                   GRAD_Chi_ar_IPs.SetRow(IP,GRAD_Chi_ar);
            }//constitutive loop ends here

            } //end Gauss integration loop



            /* saving eulerian strain for each IPs of the current element */
            fEulerian_strain_Elements_IPs.SetRow(e,fEulerian_strain_IPs);

            /* saving cauchy stress for each IPs of the current element */
            fCauchy_stress_Elements_IPs.SetRow(e,fCauchy_stress_IPs);

            /* saving state variables for each IPs of the current element */
            fState_variables_Elements_IPs.SetRow(e,fState_variables_IPs);

            /*saving displacement ??? */
            fDisplacement_Element_IPs.SetRow(e,fDisplacement_IPs);

            GammaN_IPs_el.SetRow(e,GammaN_IPs);
            SigN_IPs_el.SetRow(e,SigN_IPs);
            sn_sigman_IPs_el.SetRow(e,sn_sigman_IPs);
            mn_IPs_el.SetRow(e,mn_IPs);

            F_ar_IPs_el.SetRow(e,F_ar_IPs);
            FInv_ar_IPs_el.SetRow(e,FInv_ar_IPs);
            Chi_ar_IPs_el.SetRow(e,Chi_ar_IPs);
            GRAD_Chi_ar_IPs_el.SetRow(e,GRAD_Chi_ar_IPs);

            if(iConstitutiveModelType==1)
            {

           	//{fFd_int} will be formed
//              fFd_int  = 0.0;
              fFd_int  = Vint_1;
              fFd_int *= -1;

            //Micromorphic case fKdd coming from bal. of linear momentum
              fKdd=0.0;
/*            fKdd  =  fKu_1;
              fKdd +=  fKu_2;*/
              fKdd +=  fKu_3;
              fKdd +=  fKu_4;
              fKdd +=  fKu_5;
              fKdd +=  fKu_6;
              fKdd +=  fKu_7;
              fKdd +=  fKu_8;
              fKdd +=  fKu_9;
          //    fKdd +=  KJmat;

           /* [fKdphi] will be formed */

//             fKdphi  = 0.0;
             fKdphi =fKuphi_1;
             fKdphi +=fKuphi_2;
             fKdphi +=fKuphi_3;


           /* [fKphid] will be formed */
            //fKphid = 0.0;
            fKphid =fKphiu_1;
            fKphid+=fKphiu_2;
            fKphid+=fKphiu_3;
            fKphid+=fKphiu_4;
            fKphid+=fKphiu_5;

            fKphid+=fKFJu;
            fKphid+=fKJFu;

            fKphid+=fKMFphiu;
            fKphid+=fKMphiu_1;
            fKphid+=fKMphiu_2;


           /* [fKphiphi] will be formed */
          //need to code
           //fKphiphi = 0.0;
           fKphiphi =fKphiphi_1;
           fKphiphi+=fKphiphi_2;
           fKphiphi+=fKphiphi_3;

           fKphiphi+=fKMchiphiphi;
           fKphiphi+=fKMphiphi_1;
           fKphiphi+=fKMphiphi_2;


            /* {fFphi_int} will be formed */
//           fFphi_int  = 0.0;
           fFphi_int  = Vint_2;
           fFphi_int +=Vint_3;
           fFphi_int *=-1;

            }
            else
            {

           //{fFd_int} will be formed
            fFd_int = fFd_int_N1_vector;
            fFd_int += fFd_int_G4_vector;
            fFd_int *= -1;

           //{fFd_int} will be formed
           fFd_int  = 0.0;
           fFd_int  = Uint_1;
           //fFd_ext =-Uext_1; //no external traction is assumed
           fFd_int *= -1;


            //[fKdd] will be formed

            fKdd  = fK_dd_G3_1_matrix;
            fKdd += fK_dd_G3_2_matrix;
            fKdd += fK_dd_G3_3_matrix;
            fKdd += fK_dd_G3_4_matrix;
            fKdd += fK_dd_G4_matrix;


  //Micromorphic case fKdd coming from bal. of linear momentum

            fKdd  =  fG1_1;
            fKdd +=  fG1_2;
            fKdd +=  fG1_3;
            fKdd +=  fG1_4;
            fKdd +=  fG1_5;
            fKdd +=  fG1_6;
            fKdd +=  fG1_8;
            fKdd +=  fG1_10;
            fKdd +=  fG1_12;
            fKdd +=  fG1_13;
            fKdd +=  fG1_14;





           //[fKdphi] will be formed
           //need to code
           // fKdphi = 0.0;
           // Micromorphic case fKdPhi  from coming from bal. of linear momentum

            fKdphi  = fG1_7 ;
            fKdphi += fG1_9 ;
            fKdphi += fG1_11;



            //[fKphid] will be formed
            //need to code
            //fKphid = 0.0;

            fKphid = fH1_Etagrad;
            fKphid += fH1_1;
            fKphid += fH1_2;
            fKphid += fH1_3;
            fKphid += fH1_11;
            fKphid += fH1_12;
            fKphid += fH1_14;
            fKphid += fH1_72;
            fKphid += fH1_78;

            fKphid += fH2_1;
            fKphid += fH2_2;
            fKphid += fH2_3;
            fKphid += fH2_5;
            fKphid += fH2_7;
            fKphid += fH2_8;
            fKphid += fH2_9;
            fKphid += fH2_10;
            fKphid += fH2_12;
            fKphid += fH2_13;


            //[fKphiphi] will be formed
            //need to code
            //fKphiphi = 0.0;

            fKphiphi  = fH1_4;
            fKphiphi += fH1_5;
            fKphiphi += fH1_6;
            fKphiphi += fH1_71;
            fKphiphi += fH1_73;
            fKphiphi += fH1_74;
            fKphiphi += fH1_75;
            fKphiphi += fH1_76;
            fKphiphi += fH1_77;
            fKphiphi += fH1_7;
            fKphiphi += fH1_8;
            fKphiphi += fH1_9;
            fKphiphi += fH1_10;
            fKphiphi += fH1_13;

            fKphiphi += fH2_4;
            fKphiphi += fH2_6;
            fKphiphi += fH2_11;

            fKphiphi += fH3_1;



            // {fFphi_int} will be formed
            //need to code
            //fFphi_int  = 0.0;
            fFphi_int  =Pint_1;
            fFphi_int +=Pint_2;
            fFphi_int +=Pint_3;//no external traction is assumed Pext=0
            fFphi_int *= -1;
            }

            /* equations numbers */
            const iArrayT& displ_eq = fElementCards_displ[e].Equations();
            const iArrayT& micro_eq = fElementCards_micro[e].Equations();

            /* assemble residuals */
            ElementSupport().AssembleRHS(curr_group, fFd_int, displ_eq);
            ElementSupport().AssembleRHS(curr_group, fFphi_int, micro_eq);

            /* assemble components of the tangent */
            ElementSupport().AssembleLHS(curr_group, fKdd, displ_eq);
            ElementSupport().AssembleLHS(curr_group, fKphiphi, micro_eq);
            ElementSupport().AssembleLHS(curr_group, fKdphi, displ_eq, micro_eq);
            ElementSupport().AssembleLHS(curr_group, fKphid, micro_eq, displ_eq);

    }
    }
}



/* form global shape function derivatives */
void FSMicromorphic3DT::SetGlobalShape(void)
{
    /* fetch (initial) coordinates */
    SetLocalX(fLocInitCoords);

    /* compute shape function derivatives */
    fShapes_displ->SetDerivatives_DN_DDN();
    fShapes_micro->SetDerivatives();

}



/* describe the parameters needed by the interface */
void FSMicromorphic3DT::DefineParameters(ParameterListT& list) const
{
    /* inherited */
    ElementBaseT::DefineParameters(list);

    /* displacement field */
    //already done in ElementBaseT
    //list.AddParameter(ParameterT::Word, "displ_field_name");

    /* micro-displacement-gradient field */
    list.AddParameter(ParameterT::Word, "micro_field_name");

    list.AddParameter(fGeometryCode_displ_int, "GeometryCode_displ");
 //   list.AddParameter(fGeometryCode_micro_int, "GeometryCode_micro");
    list.AddParameter(fNumIP_displ, "NumIP_displ");
    list.AddParameter(fGeometryCodeSurf_displ_int, "GeometryCodeSurf_displ");
    list.AddParameter(fNumIPSurf_displ, "NumIPSurf_displ");
    list.AddParameter(n_en_displ, "n_en_displ");
    list.AddParameter(n_en_micro, "n_en_micro");
    list.AddParameter(ndof_per_nd_micro, "ndof_per_nd_micro");

    list.AddParameter(iConstitutiveModelType, "constitutive_mod_type");

    double shearMu, sLambda, Rho_0, gravity_g, gravity_g1, gravity_g2, gravity_g3;
    double Kappa, Nu, Sigma_const, Tau, Eta;
    double Tau1,Tau2,Tau3,Tau4,Tau5,Tau6,Tau7,Tau8,Tau9,Tau10,Tau11;

    // solid elasticity
    list.AddParameter(shearMu, "mu");
    list.AddParameter(sLambda, "lambda");
    //Material Parameter
    list.AddParameter(Kappa, "Kappa");
    list.AddParameter(Nu, "Nu");
    list.AddParameter(Sigma_const, "Sigma_const");
    list.AddParameter(Tau, "Tau");
    list.AddParameter(Eta, "Eta");
    list.AddParameter(Tau1, "Tau1");
    list.AddParameter(Tau2, "Tau2");
    list.AddParameter(Tau3, "Tau3");
    list.AddParameter(Tau4, "Tau4");
    list.AddParameter(Tau5, "Tau5");
    list.AddParameter(Tau6, "Tau6");
    list.AddParameter(Tau7, "Tau7");
    list.AddParameter(Tau8, "Tau8");
    list.AddParameter(Tau9, "Tau9");
    list.AddParameter(Tau10, "Tau10");
    list.AddParameter(Tau11, "Tau11");


    // gravity
    list.AddParameter(gravity_g, "g");

    // gravity in each direction (depends on the coordinate system which we have chosen for the problem)
    list.AddParameter(gravity_g1, "g1");
    list.AddParameter(gravity_g2, "g2");
    list.AddParameter(gravity_g3, "g3");

    // reference mass density
    list.AddParameter(Rho_0, "rho_0");

    // Newmark time integration parameters
//    list.AddParameter(newBeta, "beta");
//    list.AddParameter(newGamma, "gamma");

}


/* accept parameter list */
void FSMicromorphic3DT::TakeParameterList(const ParameterListT& list)
{
    const char caller[] = "FSMicromorphic3DT::TakeParameterList";

    /* inherited */
    ElementBaseT::TakeParameterList(list);

    /* get form of tangent */
    GlobalT::SystemTypeT type = TangentType();

    /* set form of element stiffness matrix */
    if (type == GlobalT::kSymmetric)
    fLHS.SetFormat(ElementMatrixT::kSymmetricUpper);
    else if (type == GlobalT::kNonSymmetric)
    fLHS.SetFormat(ElementMatrixT::kNonSymmetric);
    else if (type == GlobalT::kDiagonal)
    fLHS.SetFormat(ElementMatrixT::kDiagonal);

    /* get displacement field */
    /*
      const StringT& displ_field_name = list.GetParameter("displ_field_name");
      fDispl = ElementSupport().Field(displ_field_name);
      if (!fDispl)
      ExceptionT::GeneralFail(caller, "could not resolve \"%s\" displ_field",
      displ_field_name.Pointer());
    */
    const StringT& displ_field_name = list.GetParameter("field_name");
    fDispl = ElementSupport().Field(displ_field_name);
    if (!fDispl)
    ExceptionT::GeneralFail(caller, "could not resolve \"%s\" displ_field",
                displ_field_name.Pointer());

    /* get micro-displacement-gradient field */

    const StringT& micro_field_name = list.GetParameter("micro_field_name");
    fMicro = ElementSupport().Field(micro_field_name);
    if (!fMicro)
    ExceptionT::GeneralFail(caller, "could not resolve \"%s\" micro_field",
                micro_field_name.Pointer());

    fGeometryCode_displ_int = list.GetParameter("GeometryCode_displ");
    fGeometryCode_displ = GeometryT::int2CodeT(fGeometryCode_displ_int);
 //   fGeometryCode_micro_int=list.GetParameter("GeometryCode_micro");
 //   fGeometryCode_micro= GeometryT::int2CodeT(fGeometryCode_micro_int);

    fNumIP_displ = list.GetParameter("NumIP_displ");
    fGeometryCodeSurf_displ_int = list.GetParameter("GeometryCodeSurf_displ");
    fGeometryCodeSurf_displ = GeometryT::int2CodeT(fGeometryCodeSurf_displ_int);
    fNumIPSurf_displ = list.GetParameter("NumIPSurf_displ");
    n_en_displ = list.GetParameter("n_en_displ");
    n_en_micro = list.GetParameter("n_en_micro");
    ndof_per_nd_micro = list.GetParameter("ndof_per_nd_micro");

    //kAnalysisType = list.GetParameter("type_of_analysis_1consolidation_2dynamic");
    //kInitialConditionType = list.GetParameter("initial_condition_1geostatic_2displ_vel_press");

    fGeometryCode_micro = fGeometryCode_displ;
    fNumIP_micro = fNumIP_displ;
    fGeometryCodeSurf_micro = fGeometryCodeSurf_displ;
    fNumIPSurf_micro = fNumIPSurf_displ;

    iConstitutiveModelType = list.GetParameter("constitutive_mod_type");

    fMaterial_Params.Dimension ( kNUM_FMATERIAL_TERMS );
//    fIntegration_Params.Dimension ( kNUM_FINTEGRATE_TERMS );

    fMaterial_Params[kMu] = list.GetParameter("mu");
    fMaterial_Params[kLambda] = list.GetParameter("lambda");
    fMaterial_Params[kKappa] = list.GetParameter("Kappa");
    fMaterial_Params[kNu] = list.GetParameter("Nu");
    fMaterial_Params[kSigma_const] = list.GetParameter("Sigma_const");
    fMaterial_Params[kTau] = list.GetParameter("Tau");
    fMaterial_Params[kEta] = list.GetParameter("Eta");
    fMaterial_Params[kTau1] = list.GetParameter("Tau1");
    fMaterial_Params[kTau2] = list.GetParameter("Tau2");
    fMaterial_Params[kTau3] = list.GetParameter("Tau3");
    fMaterial_Params[kTau4] = list.GetParameter("Tau4");
    fMaterial_Params[kTau5] = list.GetParameter("Tau5");
    fMaterial_Params[kTau6] = list.GetParameter("Tau6");
    fMaterial_Params[kTau7] = list.GetParameter("Tau7");
    fMaterial_Params[kTau8] = list.GetParameter("Tau8");
    fMaterial_Params[kTau9] = list.GetParameter("Tau9");
    fMaterial_Params[kTau10] = list.GetParameter("Tau10");
    fMaterial_Params[kTau11] = list.GetParameter("Tau11");

    fMaterial_Params[kg] = list.GetParameter("g");
    fMaterial_Params[kg1] = list.GetParameter("g1");
    fMaterial_Params[kg2] = list.GetParameter("g2");
    fMaterial_Params[kg3] = list.GetParameter("g3");
    fMaterial_Params[kRho_0] = list.GetParameter("rho_0");

//    fIntegration_Params[kBeta] = list.GetParameter("beta");
//    fIntegration_Params[kGamma] = list.GetParameter("gamma");

    Echo_Input_Data();

    //need to change for appropriate number of ISVs for micromorphic model
    knum_d_state = 3; // #? internal state variables
    knum_i_state = 0; // int's needed per ip, state variables

    //need to change these for non-symmetric stress, and higher order couple stress output
    knumstrain = 9; // number of strain outputs
    knumstress = 9; // number of stress outputs + higher order = ??

/*    knumstrain = 6; // number of strain outputs
    knumstress = 6; // number of stress outputs + higher order = ??
*/

    knumdispl=3;

    output = "out";

    /* dimensions (notation as per Hughes' Book) */
    int& n_ip_displ = fNumIP_displ;
    int& n_ip_micro = fNumIP_micro;
    n_sd = NumSD();
    int nen = NumElementNodes(); /* number of nodes/element in the mesh */

    /* initialize connectivities */
    fConnectivities_displ.Alias(fConnectivities);
    fConnectivities_micro.Alias(fConnectivities);

    /* pick element interpolations based on available number of element nodes
     * and the specified number of integration points */
    // only implemented for 3D, quadratic hexs
    //if (n_sd == 2 && n_en_micro != n_en_displ && fGeometryCode_displ == GeometryT::kQuadrilateral)
    if (n_sd == 3 && n_en_micro != n_en_displ && fGeometryCode_displ == GeometryT::kHexahedron)
    {
        // don't expect reduced integration for both fields
        // if (n_ip_displ == 4 && n_ip_micro == 4)
        //        ExceptionT::GeneralFail(caller, "not expecting 4 ips for both fields");
        //else if (n_ip_displ == 4 || n_ip_micro == 4) // create reduced connectivities
        //{
        // reduce the number of element nodes based on the number ip's
        int& nen_red = (n_ip_displ == 8) ? n_en_displ : n_en_micro;
        nen_red = 8;
        ArrayT<const iArray2DT*>& connects_red = (n_ip_displ == 8) ?
            fConnectivities_displ :
            fConnectivities_micro;

        //create reduced connectivities
        connects_red.Dimension(0);
        connects_red.Dimension(fConnectivities.Length());
        fConnectivities_reduced.Dimension(fConnectivities.Length());
        for (int i = 0; i < fConnectivities_reduced.Length(); i++) {
            iArray2DT& connects_red_store = fConnectivities_reduced[i];
            const iArray2DT& connects = *(fConnectivities[i]);
            connects_red_store.Dimension(connects.MajorDim(), nen_red);
            connects_red[i] = &connects_red_store;

            //take 1st eight element nodes (columns)
            for (int j = 0; j < nen_red; j++)
            connects_red_store.ColumnCopy(j, connects, j);
        }
        //}
    }


    n_el = NumElements();
    n_sd_surf = n_sd;


    /* set shape functions */
    // u
    fInitCoords_displ.Dimension(n_en_displ, n_sd);
    ElementSupport().RegisterCoordinates(fInitCoords_displ);
    fCurrCoords_displ.Dimension(n_en_displ, n_sd);
    fShapes_displ = new ShapeFunctionT(fGeometryCode_displ, fNumIP_displ, fCurrCoords_displ,1 );
    //fShapes_displ->Initialize();
    // phi
    fInitCoords_micro.Dimension(n_en_micro, n_sd);
    ElementSupport().RegisterCoordinates(fInitCoords_micro);
    fCurrCoords_micro.Dimension(n_en_micro, n_sd);
    fShapes_micro = new ShapeFunctionT(fGeometryCode_micro, fNumIP_micro, fCurrCoords_micro);
    //fShapes_micro = new ShapeFunctionT(fGeometryCode_micro, fNumIP_micro, fCurrCoords_displ);
    //fShapes_micro->Initialize();

    /* set local arrays for displacement field */
    u.Dimension (n_en_displ, n_sd);
    u_dot.Dimension (n_en_displ, n_sd);
    u_dot_n.Dimension (n_en_displ, n_sd);
    u_dotdot_n.Dimension (n_en_displ, n_sd);
    u_dotdot.Dimension (n_en_displ, n_sd);
    u_n.Dimension (n_en_displ, n_sd);
    del_u.Dimension (n_en_displ, n_sd);
    n_en_displ_x_n_sd = n_en_displ*n_sd;
    n_sd_x_n_sd_x_n_sd=n_sd*n_sd*n_sd;
    n_en_micro_x_n_sd=n_en_micro*n_sd;
    n_en_micro_x_n_sd_x_n_sd=n_en_micro*n_sd_x_n_sd;
    del_u_vec.Dimension (n_en_displ_x_n_sd);
    u_vec.Dimension (n_en_displ_x_n_sd);
    u_dot_vec.Dimension (n_en_displ_x_n_sd);
    u_dotdot_vec.Dimension (n_en_displ_x_n_sd);
    //ElementSupport().RegisterCoordinates(fInitCoords_displ);
    fDispl->RegisterLocal(u);
    fDispl->RegisterLocal(u_n);


    if (fIntegrator->Order() == 1)
    {
        fDispl->RegisterLocal(u_dot);
        fDispl->RegisterLocal(u_dot_n);
    }

    if (fIntegrator->Order() == 2)
    {
        fDispl->RegisterLocal(u_dot);
        fDispl->RegisterLocal(u_dot_n);
        fDispl->RegisterLocal(u_dotdot);
        fDispl->RegisterLocal(u_dotdot_n);
    }

    /* set local arrays for micro-displacement-gradient field */
    Phi.Dimension (n_en_micro, ndof_per_nd_micro);
    Phi_dot.Dimension (n_en_micro, ndof_per_nd_micro);
    Phi_dot_n.Dimension (n_en_micro, ndof_per_nd_micro);
    Phi_dotdot.Dimension (n_en_micro, ndof_per_nd_micro);
    Phi_dotdot_n.Dimension (n_en_micro, ndof_per_nd_micro);
    Phi_n.Dimension (n_en_micro, ndof_per_nd_micro);
    del_Phi.Dimension (n_en_micro, ndof_per_nd_micro);
    n_en_micro_x_ndof_per_nd_micro = n_en_micro*ndof_per_nd_micro;
    ndof_per_nd_micro_x_n_sd=ndof_per_nd_micro*n_sd;
    del_Phi_vec.Dimension (n_en_micro_x_ndof_per_nd_micro);
    Phi_vec.Dimension (n_en_micro_x_ndof_per_nd_micro);
    Phi_dot_vec.Dimension (n_en_micro_x_ndof_per_nd_micro);
    Phi_dotdot_vec.Dimension (n_en_micro_x_ndof_per_nd_micro);
    //ElementSupport().RegisterCoordinates(fInitCoords_micro);

    fMicro->RegisterLocal(Phi);
    fMicro->RegisterLocal(Phi_n);

    if (fIntegrator->Order() == 1)
    {
        fMicro->RegisterLocal(Phi_dot);
        fMicro->RegisterLocal(Phi_dot_n);
    }

    if (fIntegrator->Order() == 2)
    {
        fMicro->RegisterLocal(Phi_dot);
        fMicro->RegisterLocal(Phi_dot_n);
        fMicro->RegisterLocal(Phi_dotdot);
        fMicro->RegisterLocal(Phi_dotdot_n);
    }


    /* allocate state variable storage */
    // state variables are calculated at IPs for displacement field
    int num_ip = fNumIP_displ;
    //int num_ip = fNumIP_micro;
    fdState_new.Dimension(n_el, num_ip*knum_d_state);
    fdState.Dimension(n_el, num_ip*knum_d_state);
    fiState_new.Dimension(n_el, num_ip*knum_i_state);
    fiState.Dimension(n_el, num_ip*knum_i_state);

    /* initialize equations */
    fEqnos_displ.Alias(fEqnos_displ);
    fEqnos_micro.Dimension(fConnectivities_micro.Length());

    /* initialize state variables */
    fdState = 0;
    fdState_new = 0;
    fiState = 0;
    fiState_new = 0;

    /* initialize element cards */
    fElementCards_displ.Alias(fElementCards);
    fElementCards_micro.Dimension(fElementCards.Length());

    /* set cards to data in array - NOT NEEDED IF YOU'RE NOT
     * GOING TO USE THE ElementCardT ARRAY? */
    for (int i= 0; i < fElementCards.Length(); i++)
    fElementCards[i].Set(fiState.MinorDim(), fiState(i), fdState.MinorDim(), fdState(i));

    fKdd.Dimension          ( n_en_displ_x_n_sd, n_en_displ_x_n_sd );
    fKdphi.Dimension        ( n_en_displ_x_n_sd, n_en_micro_x_ndof_per_nd_micro );
    fKphid.Dimension        ( n_en_micro_x_ndof_per_nd_micro, n_en_displ_x_n_sd );
    fKphiphi.Dimension      ( n_en_micro_x_ndof_per_nd_micro, n_en_micro_x_ndof_per_nd_micro );

    fFd_int.Dimension       ( n_en_displ_x_n_sd );
    fFd_ext.Dimension       ( n_en_displ_x_n_sd );
    fFphi_int.Dimension     ( n_en_micro_x_ndof_per_nd_micro );
    fFphi_ext.Dimension     ( n_en_micro_x_ndof_per_nd_micro );

    /* workspace matricies */
    fShapeDispl.Dimension (n_sd, n_en_displ_x_n_sd);
    fShapeDispl_Tr.Dimension (n_en_displ_x_n_sd,n_sd);
    fShapeMicro.Dimension (ndof_per_nd_micro, n_en_micro_x_ndof_per_nd_micro);
    n_sd_x_n_sd = n_sd*n_sd;
    n_sd_x_n_sd_x_n_sd = n_sd*n_sd*n_sd;
    fShapeDisplGrad_temp.Dimension (n_sd, n_en_displ);
    fShapeDisplGrad.Dimension (n_sd_x_n_sd, n_en_displ_x_n_sd);
    fShapeDisplGrad_t.Dimension (n_sd_x_n_sd, n_en_displ_x_n_sd);
    fShapeDisplGradGrad.Dimension (n_sd*2 , n_en_displ);
   // ndof_per_nd_micro_x_n_sd_x_n_sd = ndof_per_nd_micro*n_sd*n_sd;
    fShapeMicroGrad_temp.Dimension(n_sd,n_en_micro);
    fShapeMicroGrad.Dimension (ndof_per_nd_micro*n_sd, n_en_micro*ndof_per_nd_micro);
    fDeformation_Gradient.Dimension (n_sd,n_sd);
    fGrad_disp_vector.Dimension (n_sd_x_n_sd);
    fGrad_disp_matrix.Dimension (n_sd,n_sd);
    fDeformation_Gradient_Inverse.Dimension (n_sd,n_sd);
    fDeformation_Gradient_Transpose.Dimension (n_sd,n_sd);
    fDeformation_Gradient_Inverse_Transpose.Dimension (n_sd,n_sd);
    fDefGradInv_Grad_grad.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fDefGradInv_Grad_grad_Transpose.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fDefGradT_9x9_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fDefGradInv_vector.Dimension (n_sd_x_n_sd);
    fRight_Cauchy_Green_tensor.Dimension (n_sd,n_sd);
    fRight_Cauchy_Green_tensor_Inverse.Dimension (n_sd,n_sd);
    fLeft_Cauchy_Green_tensor.Dimension (n_sd,n_sd);
    fIdentity_matrix.Dimension (n_sd,n_sd);
    fSecond_Piola_tensor.Dimension (n_sd,n_sd);
    fTemp_matrix_nsd_x_nsd.Dimension (n_sd,n_sd);
    fKirchhoff_tensor.Dimension (n_sd,n_sd);
    fKirchhoff_vector.Dimension (n_sd_x_n_sd);
    fIota_temp_matrix.Dimension (n_en_displ_x_n_sd,n_sd_x_n_sd);
    fVarpi_temp_matrix.Dimension (n_sd, n_en_displ_x_n_sd);
    ///////////////////////////////////////////////////////////////////////////
    /////////////DIMENSIONALIZE MICROMORPHIC MATRICES FOR 3D CASE//////////////
    ///////////////////////////////////////////////////////////////////////////
    fIota_w_temp_matrix.Dimension(n_en_displ_x_n_sd,n_sd_x_n_sd);
    fIota_eta_temp_matrix.Dimension(n_en_micro*n_sd_x_n_sd,n_sd_x_n_sd_x_n_sd);


    SigN_m.Dimension(n_sd,n_sd);
    Fn_m.Dimension(n_sd,n_sd);
    Finv_m.Dimension(n_sd,n_sd);
    deltaL.Dimension(n_sd,n_sd);
    deltaL_Tr.Dimension(n_sd,n_sd);
    deltad.Dimension(n_sd,n_sd);
    tempSig.Dimension(n_sd,n_sd);

    deltaEp.Dimension(n_sd,n_sd);
    deltaNu.Dimension(n_sd,n_sd);
    ChiN_m.Dimension(n_sd,n_sd);
    Chi_vec.Dimension(n_sd_x_n_sd);
    GRAD_Nuw.Dimension(n_sd_x_n_sd,n_en_displ_x_n_sd);
    GRAD_Chi_vec.Dimension(n_sd_x_n_sd_x_n_sd);
    NCHI.Dimension(n_sd_x_n_sd,n_en_micro*n_sd_x_n_sd);
    NCHI_Tr.Dimension(n_en_micro*n_sd_x_n_sd,n_sd_x_n_sd);
//    NCHI_eta.Dimension(n_sd_x_n_sd,n_en_micro_x_n_sd);same with the one above no need!

    fTemp_matrix_nudof_x_nudof.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fTemp_matrix_nudof_x_nchidof.Dimension (n_en_displ_x_n_sd,n_en_micro*n_sd_x_n_sd);
    fTemp_matrix_nchidof_x_nchidof.Dimension (n_en_micro*n_sd_x_n_sd,n_en_micro*n_sd_x_n_sd);
    fTemp_matrix_nchidof_x_nudof.Dimension (n_en_micro*n_sd_x_n_sd,n_en_displ_x_n_sd);


    Sigma.Dimension(n_sd,n_sd);
    SigN_ar.Dimension(n_sd,n_sd);
    SigN_ar=0.0;
    SigN_IPs.Dimension(fNumIP_displ,n_sd_x_n_sd);
    SigN_IPs_n.Dimension(fNumIP_displ,n_sd_x_n_sd);
    SigN_IPs_n=0.0;
    Sig_IPs.Dimension(fNumIP_displ,n_sd_x_n_sd);
    SigN_IPs_el.Dimension(n_el,n_sd_x_n_sd*fNumIP_displ);
    Sig_IPs_el.Dimension(n_el,n_sd_x_n_sd*fNumIP_displ);
    SigN_IPs_el_n.Dimension(n_el,n_sd_x_n_sd*fNumIP_displ);
    SigN_IPs_el=0.0;
    SigN_IPs_el_n=0.0;

 //   Counter.Dimension(1);
/*    Counter_IPs_el_n.Dimension(n_el,fNumIP_displ);
    Counter_IPs_el_n=0.0;
    Counter_IPs_el.Dimension(n_el,fNumIP_displ);
    Counter_IPs.Dimension(fNumIP_displ);*/

    TransShapeDisplGrad.Dimension(n_en_displ_x_n_sd,n_sd_x_n_sd);
    Var_F.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    Finv_w.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    Tsigma_1.Dimension (n_sd_x_n_sd,n_sd_x_n_sd);
    fG1_1.Dimension (n_en_displ_x_n_sd ,n_en_displ_x_n_sd );
    Tsigma_2.Dimension (n_sd_x_n_sd,n_sd_x_n_sd);
    fG1_2.Dimension (n_en_displ_x_n_sd ,n_en_displ_x_n_sd );
    Tsigma_3.Dimension (n_sd_x_n_sd,n_sd_x_n_sd);
    fG1_3.Dimension (n_en_displ_x_n_sd ,n_en_displ_x_n_sd );
    TFn_1.Dimension (n_sd_x_n_sd,n_sd_x_n_sd);
    fG1_4.Dimension (n_en_displ_x_n_sd ,n_en_displ_x_n_sd );
    TFn_2.Dimension (n_sd_x_n_sd,n_sd_x_n_sd);
    fG1_5.Dimension (n_en_displ_x_n_sd ,n_en_displ_x_n_sd );
    TFn_3.Dimension (n_sd_x_n_sd,n_sd_x_n_sd);
    fG1_6.Dimension (n_en_displ_x_n_sd ,n_en_displ_x_n_sd );
    TChi_1.Dimension (n_sd_x_n_sd,n_sd_x_n_sd);
    fG1_7.Dimension (n_en_displ_x_n_sd ,n_en_micro*n_sd_x_n_sd );
    TFn_4.Dimension (n_sd_x_n_sd,n_sd_x_n_sd);
    fG1_8.Dimension (n_en_displ_x_n_sd ,n_en_displ_x_n_sd);
    TChi_2.Dimension (n_sd_x_n_sd,n_sd_x_n_sd);
    fG1_9.Dimension (n_en_displ_x_n_sd ,n_en_micro*n_sd_x_n_sd);
    TFn_5.Dimension (n_sd_x_n_sd,n_sd_x_n_sd);
    fG1_10.Dimension (n_en_displ_x_n_sd ,n_en_displ_x_n_sd);
    TChi_3.Dimension (n_sd_x_n_sd,n_sd_x_n_sd);
    fG1_11.Dimension (n_en_displ_x_n_sd ,n_en_micro*n_sd_x_n_sd);
    TFn_6.Dimension (n_sd_x_n_sd,n_sd_x_n_sd);
    fG1_12.Dimension (n_en_displ_x_n_sd ,n_en_displ_x_n_sd);
    SigCurr.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    fG1_13.Dimension (n_en_displ_x_n_sd ,n_en_displ_x_n_sd);
    fG1_14.Dimension (n_en_displ_x_n_sd ,n_en_displ_x_n_sd);

    //Bal. of linear Mom of Momtm
    Etagrad.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd);
    Mm_1.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd);
    Mm_2.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd);
    Mm_3.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd);
    Mm_4.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd);
    Mm_5.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd);
    Mm_6.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd);
    Mm_7.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd);
    Mm_71.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd_x_n_sd);
    Mm_72.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd);
    Mm_73.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd);
    Mm_74.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd);
    Mm_75.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd);
    Mm_76.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd_x_n_sd);
    Mm_77.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd);
    Mm_78.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd);
    Mm_8.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd);
    Mm_9.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd_x_n_sd);
    Mm_10.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd);
    Mm_11.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd);
    Mm_12.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd);
    Mm_13.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd);
    Mm_14.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd);
    Ru_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    Ru_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    Ru_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    Ru_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    Ru_5.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    Ru_6.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    Ru_7.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    Ru_8.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    Ru_9.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    RChi_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    RChi_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    RChi_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    Rs_sigma.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    R_Capital_Gamma_Chi.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    CapitalLambda.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);

    mn_ar.Dimension(n_sd_x_n_sd_x_n_sd);
    mn_ar=0.0;
    mn_IPs.Dimension(fNumIP_displ,n_sd_x_n_sd_x_n_sd);
    mn_IPs=0.0;
    mn_IPs_n.Dimension(fNumIP_displ,n_sd_x_n_sd_x_n_sd);
    mn_IPs_n=0.0;
    mn_IPs_el.Dimension(n_el,n_sd_x_n_sd_x_n_sd*fNumIP_displ);
    mn_IPs_el=0.0;
    mn_IPs_el_n.Dimension(n_el,n_sd_x_n_sd_x_n_sd*fNumIP_displ);
    mn_IPs_el=0.0;
    mn_IPs_el_n=0.0;
    GammaN_ar.Dimension(n_sd_x_n_sd_x_n_sd);
    GammaN_ar=0.0;
    GammaN_IPs.Dimension(fNumIP_displ,n_sd_x_n_sd_x_n_sd);
    GammaN_IPs=0.0;
    GammaN_IPs_n.Dimension(fNumIP_displ,n_sd_x_n_sd_x_n_sd);
    GammaN_IPs_n=0.0;
    GammaN_IPs_el.Dimension(n_el,n_sd_x_n_sd_x_n_sd*fNumIP_displ);
    GammaN_IPs_el=0.0;
    GammaN_IPs_el_n.Dimension(n_el,n_sd_x_n_sd_x_n_sd*fNumIP_displ);
    GammaN_IPs_el=0.0;
    GammaN_IPs_el_n=0.0;

    sn_sigman.Dimension(n_sd,n_sd);
    sn_sigman=0.0;
    sn_sigman_IPs.Dimension(fNumIP_displ,n_sd_x_n_sd);
    sn_sigman_IPs=0.0;
    sn_sigman_IPs_n.Dimension(fNumIP_displ,n_sd_x_n_sd);
    SigN_IPs_n=0.0;
    sn_sigman_IPs_el.Dimension(n_el,fNumIP_displ*n_sd_x_n_sd);
    sn_sigman_IPs_el=0.0;
    sn_sigman_IPs_el_n.Dimension(n_el,fNumIP_displ*n_sd_x_n_sd);
    sn_sigman_IPs_el_n=0.0;
    s_sigma.Dimension(n_sd,n_sd);
    s_sigma=0.0;

    GRAD_NCHI.Dimension(n_sd_x_n_sd_x_n_sd,n_en_micro*n_sd_x_n_sd);
 //   GRAD_NCHI_Phi.Dimension(n_sd_x_n_sd_x_n_sd,n_en_micro_x_n_sd_x_n_sd);
    Finv_eta.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd_x_n_sd);
    fH1_Etagrad.Dimension(n_en_micro*n_sd_x_n_sd,n_en_displ_x_n_sd);
    fH1_1.Dimension(n_en_micro*n_sd_x_n_sd,n_en_displ_x_n_sd);
    fH1_2.Dimension(n_en_micro*n_sd_x_n_sd,n_en_displ_x_n_sd);
    fH1_3.Dimension(n_en_micro*n_sd_x_n_sd,n_en_displ_x_n_sd);
    fH1_4.Dimension(n_en_micro*n_sd_x_n_sd,n_en_micro*n_sd_x_n_sd);
    fH1_5.Dimension(n_en_micro*n_sd_x_n_sd,n_en_micro*n_sd_x_n_sd);
    fH1_6.Dimension(n_en_micro*n_sd_x_n_sd,n_en_micro*n_sd_x_n_sd);
    fH1_7.Dimension(n_en_micro*n_sd_x_n_sd,n_en_micro*n_sd_x_n_sd);
    fH1_71.Dimension(n_en_micro*n_sd_x_n_sd,n_en_micro*n_sd_x_n_sd);
    fH1_72.Dimension(n_en_micro*n_sd_x_n_sd,n_en_displ_x_n_sd);
    fH1_73.Dimension(n_en_micro*n_sd_x_n_sd,n_en_micro*n_sd_x_n_sd);
    fH1_74.Dimension(n_en_micro*n_sd_x_n_sd,n_en_micro*n_sd_x_n_sd);
    fH1_75.Dimension(n_en_micro*n_sd_x_n_sd,n_en_micro*n_sd_x_n_sd);
    fH1_76.Dimension(n_en_micro*n_sd_x_n_sd,n_en_micro*n_sd_x_n_sd);
    fH1_77.Dimension(n_en_micro*n_sd_x_n_sd,n_en_micro*n_sd_x_n_sd);
    fH1_78.Dimension(n_en_micro*n_sd_x_n_sd,n_en_displ_x_n_sd);
    fH1_8.Dimension(n_en_micro*n_sd_x_n_sd,n_en_micro*n_sd_x_n_sd);
    fH1_9.Dimension(n_en_micro*n_sd_x_n_sd,n_en_micro*n_sd_x_n_sd);
    fH1_10.Dimension(n_en_micro*n_sd_x_n_sd,n_en_micro*n_sd_x_n_sd);
    fH1_11.Dimension(n_en_micro*n_sd_x_n_sd,n_en_displ_x_n_sd);
    fH1_12.Dimension(n_en_micro*n_sd_x_n_sd,n_en_displ_x_n_sd);
    fH1_13.Dimension(n_en_micro*n_sd_x_n_sd,n_en_micro*n_sd_x_n_sd);
    fH1_14.Dimension(n_en_micro*n_sd_x_n_sd,n_en_displ_x_n_sd);

    fH2_1.Dimension(n_en_micro*n_sd_x_n_sd,n_en_displ_x_n_sd);
    fH2_2.Dimension(n_en_micro*n_sd_x_n_sd,n_en_displ_x_n_sd);
    fH2_3.Dimension(n_en_micro*n_sd_x_n_sd,n_en_displ_x_n_sd);
    fH2_4.Dimension(n_en_micro*n_sd_x_n_sd,n_en_micro*n_sd_x_n_sd);
    fH2_5.Dimension(n_en_micro*n_sd_x_n_sd,n_en_displ_x_n_sd);
    fH2_6.Dimension(n_en_micro*n_sd_x_n_sd,n_en_micro*n_sd_x_n_sd);
    fH2_7.Dimension(n_en_micro*n_sd_x_n_sd,n_en_displ_x_n_sd);
    fH2_8.Dimension(n_en_micro*n_sd_x_n_sd,n_en_displ_x_n_sd);
    fH2_9.Dimension(n_en_micro*n_sd_x_n_sd,n_en_displ_x_n_sd);
    fH2_10.Dimension(n_en_micro*n_sd_x_n_sd,n_en_displ_x_n_sd);
    fH2_11.Dimension(n_en_micro*n_sd_x_n_sd,n_en_micro*n_sd_x_n_sd);
    fH2_12.Dimension(n_en_micro*n_sd_x_n_sd,n_en_displ_x_n_sd);
    fH2_13.Dimension(n_en_micro*n_sd_x_n_sd,n_en_displ_x_n_sd);
    fH3_1.Dimension(n_en_micro*n_sd_x_n_sd,n_en_micro*n_sd_x_n_sd);

     Mat1.Dimension(9,9);
     Mat2.Dimension(9,9);
     Mat3.Dimension(9,9);
     Mat4.Dimension(9,9);
     Mat5.Dimension(9,9);
     Mat5_Inv.Dimension(9,9);
   //  RHS.Dimension(9);
     Sigma1.Dimension(3,3);
     Sigma2.Dimension(3,3);
     Sigma3.Dimension(3,3);
     Sigma4.Dimension(3,3);
     Sigma5.Dimension(3,3);
     Sigma6.Dimension(3,3);
     Sigma7.Dimension(3,3);
     Sigma8.Dimension(3,3);

     dMatrixT s_sigma;



    G1.Dimension(n_sd_x_n_sd);
    Uint_1.Dimension(n_en_displ_x_n_sd);
    Uint_1_temp.Dimension(n_en_displ_x_n_sd);
    Uext_1.Dimension(n_en_displ_x_n_sd);
    Gext.Dimension(n_en_displ_x_n_sd);

    H1.Dimension(n_sd_x_n_sd_x_n_sd);
    H2.Dimension(n_sd_x_n_sd);
    H3.Dimension(n_sd_x_n_sd);
    Pint_1.Dimension(n_en_micro*n_sd_x_n_sd);
    Pint_1_temp.Dimension(n_en_micro*n_sd_x_n_sd);
    Pint_2.Dimension(n_en_micro*n_sd_x_n_sd);
    Pint_2_temp.Dimension(n_en_micro*n_sd_x_n_sd);
    Pint_3.Dimension(n_en_micro*n_sd_x_n_sd);
    Pint_3_temp.Dimension(n_en_micro*n_sd_x_n_sd);
    Hext.Dimension(n_en_micro*n_sd_x_n_sd);

    Lambda.Dimension(n_sd,n_sd);//small lambda
    Lambda=0.0;//10^8
    Omega.Dimension(n_sd,n_sd);//small omega
    Omega=0.0;

    Fn_ar.Dimension(n_sd,n_sd);
    Fn_ar=0.0;
    Fn_ar_IPs.Dimension(fNumIP_displ,n_sd_x_n_sd);
    Fn_ar_IPs=0.0;
    F_ar_IPs.Dimension(fNumIP_displ,n_sd_x_n_sd);
    F_ar_IPs=0.0;
    Fn_ar_IPs_el.Dimension(n_el,fNumIP_displ*n_sd_x_n_sd);
    F_ar_IPs_el.Dimension(n_el,fNumIP_displ*n_sd_x_n_sd);

    FnInv_ar.Dimension(n_sd,n_sd);
    FnInv_ar=0.0;
    FnInv_ar_IPs.Dimension(fNumIP_displ,n_sd_x_n_sd);
    FnInv_ar_IPs=0.0;
    FInv_ar_IPs.Dimension(fNumIP_displ,n_sd_x_n_sd);
    FInv_ar_IPs=0.0;
  //  FInv_ar_IPs_n.Dimension(fNumIP_displ,n_sd_x_n_sd);
    FnInv_ar_IPs_el.Dimension(n_el,fNumIP_displ*n_sd_x_n_sd);
    FnInv_ar_IPs_el=0.0;
    FInv_ar_IPs_el.Dimension(n_el,fNumIP_displ*n_sd_x_n_sd);
    FInv_ar_IPs_el=0.0;

    Chi_m.Dimension(n_sd,n_sd);
    ChiInv_m.Dimension(n_sd,n_sd);
    ChiN_ar.Dimension(n_sd,n_sd);
    ChiN_ar=0.0;
    Chi_ar.Dimension(n_sd,n_sd);
    Chi_ar=0.0;
    ChiN_ar_IPs.Dimension(fNumIP_displ,n_sd_x_n_sd);
    ChiN_ar_IPs=0.0;
    Chi_ar_IPs.Dimension(fNumIP_displ,n_sd_x_n_sd);
    Chi_ar_IPs=0.0;
    ChiN_ar_IPs_n.Dimension(fNumIP_displ,n_sd_x_n_sd);
    ChiN_ar_IPs_n=0.0;
    ChiN_ar_IPs_el.Dimension(n_el,fNumIP_displ*n_sd_x_n_sd);
    ChiN_ar_IPs_el=0.0;
    Chi_ar_IPs_el.Dimension(n_el,fNumIP_displ*n_sd_x_n_sd);
    Chi_ar_IPs_el=0.0;
    ChiN_ar_IPs_el_n.Dimension(n_el,fNumIP_displ*n_sd_x_n_sd);
    ChiN_ar_IPs_el_n=0.0;

    GRAD_Chi_ar.Dimension(n_sd_x_n_sd_x_n_sd);
    GRAD_Chi_ar=0.0;
    GRAD_ChiN_ar.Dimension(n_sd_x_n_sd_x_n_sd);
    GRAD_Chi_ar.Dimension(n_sd_x_n_sd_x_n_sd);
    GRAD_Chi_ar=0.0;
    GRAD_ChiN_ar_IPs.Dimension(fNumIP_displ,n_sd_x_n_sd_x_n_sd);
    GRAD_Chi_ar_IPs.Dimension(fNumIP_displ,n_sd_x_n_sd_x_n_sd);
    GRAD_ChiN_ar_IPs=0.0;
    GRAD_ChiN_ar_IPs_n.Dimension(fNumIP_displ,n_sd_x_n_sd_x_n_sd);
    GRAD_ChiN_ar_IPs_el.Dimension(n_el,fNumIP_displ*n_sd_x_n_sd_x_n_sd);
    GRAD_ChiN_ar_IPs_el=0.0;
    GRAD_Chi_ar_IPs_el.Dimension(n_el,fNumIP_displ*n_sd_x_n_sd_x_n_sd);
    GRAD_Chi_ar_IPs_el=0.0;
    GRAD_ChiN_ar_IPs_el_n.Dimension(n_el,fNumIP_displ*n_sd_x_n_sd_x_n_sd);
    GRAD_ChiN_ar_IPs_el_n=0.0;
    fIdentity_matrix=0.0;
    Temp_Identity_array.Dimension(9);
    for (int i=0; i<9; i++) Temp_Identity_array[i] = 0.0;
    Temp_Identity_array[0]=1.0;
    Temp_Identity_array[4]=1.0;
    Temp_Identity_array[8]=1.0;


    /////////////////////////////////////////////////////////////
    /////////////FINITE STRAIN ELASTICITY MATRICES///////////////
    /////////////////////////////////////////////////////////////
    KirchhoffST.Dimension(n_sd,n_sd);
    SPK.Dimension(n_sd,n_sd);
    Temp_SPK.Dimension(n_sd,n_sd);
   // FSF.Dimension(n_sd,n_sd);
    LagrangianStn.Dimension(n_sd,n_sd);
    MicroStnTensor.Dimension(n_sd,n_sd);
    ChiM.Dimension(n_sd,n_sd);
    PSI.Dimension(n_sd,n_sd);
    I1_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    I1_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    I1_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    I1_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    I1_5.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    I1_6.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    I1_7.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    I2_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    I1_8.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    I2_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    I1_9.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    I2_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    fKu_1.Dimension(n_en_displ_x_n_sd ,n_en_displ_x_n_sd );
    fKu_2.Dimension(n_en_displ_x_n_sd ,n_en_displ_x_n_sd );
    fKu_3.Dimension(n_en_displ_x_n_sd ,n_en_displ_x_n_sd );
    fKu_4.Dimension(n_en_displ_x_n_sd ,n_en_displ_x_n_sd );
    fKu_5.Dimension(n_en_displ_x_n_sd ,n_en_displ_x_n_sd );
    fKu_6.Dimension(n_en_displ_x_n_sd ,n_en_displ_x_n_sd );
    fKu_7.Dimension(n_en_displ_x_n_sd ,n_en_displ_x_n_sd );
    fKuphi_1.Dimension (n_en_displ_x_n_sd,n_en_micro*n_sd_x_n_sd);
    fKu_8.Dimension(n_en_displ_x_n_sd ,n_en_displ_x_n_sd );
    fKuphi_2.Dimension (n_en_displ_x_n_sd,n_en_micro*n_sd_x_n_sd);
    fKu_9.Dimension(n_en_displ_x_n_sd ,n_en_displ_x_n_sd );
    fKuphi_3.Dimension (n_en_displ_x_n_sd,n_en_micro*n_sd_x_n_sd);
    Vint_1.Dimension(n_en_displ_x_n_sd);
    Vint_2.Dimension(n_en_micro*n_sd_x_n_sd);
    Vint_3.Dimension(n_en_micro*n_sd_x_n_sd);
    fV1.Dimension(n_sd_x_n_sd);
    fV2.Dimension(n_sd_x_n_sd);
    fV3.Dimension(n_sd_x_n_sd_x_n_sd);
    Vint_1_temp.Dimension(n_en_displ_x_n_sd);
    Vint_2_temp.Dimension(n_en_micro*n_sd_x_n_sd);
    Vint_3_temp.Dimension(n_en_micro*n_sd_x_n_sd);
    SIGMA_S.Dimension(n_sd,n_sd);
    fMKLM.Dimension(n_sd,n_sd,n_sd);
    GAMMA.Dimension(n_sd,n_sd,n_sd);
    GRAD_CHIM.Dimension(n_sd,n_sd,n_sd);
    fTemp_tensor_n_sd_x_n_sd_x_nsd.Dimension(n_sd,n_sd,n_sd);

    fKFJu.Dimension(n_en_micro*n_sd_x_n_sd,n_en_displ_x_n_sd);
    fKJFu.Dimension(n_en_micro*n_sd_x_n_sd,n_en_displ_x_n_sd);

    fKphiu_1.Dimension(n_en_micro*n_sd_x_n_sd,n_en_displ_x_n_sd);
    fKphiu_2.Dimension(n_en_micro*n_sd_x_n_sd,n_en_displ_x_n_sd);
    fKphiu_3.Dimension(n_en_micro*n_sd_x_n_sd,n_en_displ_x_n_sd);
    fKphiu_4.Dimension(n_en_micro*n_sd_x_n_sd,n_en_displ_x_n_sd);
    fKphiphi_1.Dimension(n_en_micro*n_sd_x_n_sd,n_en_micro*n_sd_x_n_sd);
    fKphiu_5.Dimension(n_en_micro*n_sd_x_n_sd,n_en_displ_x_n_sd);
    fKphiphi_2.Dimension(n_en_micro*n_sd_x_n_sd,n_en_micro*n_sd_x_n_sd);
    fKphiphi_3.Dimension(n_en_micro*n_sd_x_n_sd,n_en_micro*n_sd_x_n_sd);

    fKMFphiu.Dimension(n_en_micro*n_sd_x_n_sd,n_en_displ_x_n_sd);
    fKMchiphiphi.Dimension(n_en_micro*n_sd_x_n_sd,n_en_micro*n_sd_x_n_sd);
    fKMphiu_1.Dimension(n_en_micro*n_sd_x_n_sd,n_en_displ_x_n_sd);
    fKMphiphi_1.Dimension(n_en_micro*n_sd_x_n_sd,n_en_micro*n_sd_x_n_sd);
    fKMphiu_2.Dimension(n_en_micro*n_sd_x_n_sd,n_en_displ_x_n_sd);
    fKMphiphi_2.Dimension(n_en_micro*n_sd_x_n_sd,n_en_micro*n_sd_x_n_sd);
    fFJ.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    fJF.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    fJ1_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    fJ1_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    fJ1_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    fJ1_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    fJ2_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    fJ1_5.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    fJ1_6.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    fJ2_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    fJ2_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);

    fEtaM.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd);
    fFM.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd);
    fMF.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd);
    fMchi.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd);
    fMpu_1.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd);
    fMpp_1.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd_x_n_sd);
    fMpu_2.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd);
    fMpp_2.Dimension(n_sd_x_n_sd_x_n_sd,n_sd_x_n_sd_x_n_sd);

    Jmat.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
    KJmat.Dimension(n_en_displ_x_n_sd ,n_en_displ_x_n_sd );

    int element_number=4;

    u_el.Dimension(4,n_sd);
    u_element.Dimension(n_sd);
    ftemp_u_element.Dimension(n_sd);
    //////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////







/*    int row=0;
    for(int i=0;i<3;i++)
    {
        fIdentity_matrix(i,i)=1.0;
        for(int j=0;j<3;j++)
            {
                //Fn_ar[row]=fIdentity_matrix(i,j);
                //FnInv_ar[row]=fIdentity_matrix(i,j);
                //ChiN_ar[row]=fIdentity_matrix(i,j);
                row++;
            }
    }*/

    while (NextElement())
    {
        int e,l;
        e = CurrElementNumber();
        for (l=0; l < fNumIP_displ; l++)
        {
            Fn_ar_IPs.SetRow(l,Temp_Identity_array);
            FnInv_ar_IPs.SetRow(l,Temp_Identity_array);
            ChiN_ar_IPs.SetRow(l,Temp_Identity_array);
        }
        Fn_ar_IPs_el.SetRow(e,Fn_ar_IPs);
        FnInv_ar_IPs_el.SetRow(e,FnInv_ar_IPs);
        ChiN_ar_IPs_el_n.SetRow(e,ChiN_ar_IPs);
    }
/*    ChiN_ar_IPs_el= ChiN_ar_IPs_el_n;
    FInv_ar_IPs_el=FnInv_ar_IPs_el_n;
    F_ar_IPs_el=Fn_ar_IPs_el_n;*///no need for this because this is initializing and _el parts are cal_ed in gauss loop
//here is take parameter list function



    ///////////////////////////////////////////////////////////////////////////
    /////////////DIMENSIONALIZE MICROMORPHIC MATRICES FINISH HERE FOR 3D CASE//////////////
    ///////////////////////////////////////////////////////////////////////////


    fChi_temp_vector.Dimension (n_sd);
    fTemp_vector_9x1.Dimension (n_sd_x_n_sd);
    fFd_int_N1_vector.Dimension (n_en_displ_x_n_sd);
    fTemp_vector_ndof_se.Dimension (n_en_displ_x_n_sd);
    fTemp_vector_nen_micro.Dimension (n_en_micro);
    fIm_temp_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fHbar_temp_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fEll_temp_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fPi_temp_transpose_vector.Dimension (n_en_displ_x_n_sd);
    fPi_temp_row_matrix.Dimension (1,n_en_displ_x_n_sd);

    fK_dd_G3_1_matrix.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fK_dd_G3_2_matrix.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fK_dd_G3_3_matrix.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fK_dd_G3_4_matrix.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fI_ij_column_matrix.Dimension (n_sd_x_n_sd, 1);
    fShapeDisplGrad_t_Transpose.Dimension (n_en_displ_x_n_sd, n_sd_x_n_sd);
    fShapeMicro_row_matrix.Dimension (1,n_en_micro);
    fTemp_nsd_vector.Dimension (n_sd);
    fGrad_1_J_vector.Dimension (n_sd);
    fChi_temp_column_matrix.Dimension (n_sd, 1);
    fTemp_matrix_nsd_x_1.Dimension (n_sd,1);
    fTemp_matrix_ndof_se_x_ndof_se.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fTemp_matrix_ndof_se_x_nen_micro.Dimension (n_en_displ_x_n_sd,n_en_micro);
    fc_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fC_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fIm_Prim_temp_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fB_matrix.Dimension (6 ,n_en_displ_x_n_sd);
    fD_matrix.Dimension (6,6);
    fK_dd_BTDB_matrix.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fFd_int_smallstrain_vector.Dimension (n_en_displ_x_n_sd);
    fEulerian_strain_tensor_current_IP.Dimension (n_sd,n_sd);
    fCauchy_stress_tensor_current_IP.Dimension (n_sd,n_sd);
    //
    fDisplacements_current_IPs.Dimension(n_sd);
    fEulerian_strain_IPs.Dimension (fNumIP_displ,knumstrain);
    fCauchy_stress_IPs.Dimension (fNumIP_displ,knumstress);
    fState_variables_IPs.Dimension (fNumIP_displ,knum_d_state);
   //
    fDisplacement_IPs.Dimension(fNumIP_displ,knumdispl);
    fTemp_nine_values.Dimension(9);
    fTemp_six_values.Dimension(6);
    fEulerian_strain_Elements_IPs.Dimension (NumElements(),fNumIP_displ*knumstrain);
    fCauchy_stress_Elements_IPs.Dimension (NumElements(),fNumIP_displ*knumstress);
    fState_variables_Elements_IPs.Dimension (NumElements(),fNumIP_displ*knum_d_state);
    fDisplacement_Element_IPs.Dimension(NumElements(),fNumIP_displ*knumdispl);
    fGravity_vector.Dimension (n_sd);
    fFd_int_G4_vector.Dimension (n_en_displ_x_n_sd);
    fDefGradInv_column_matrix.Dimension (n_sd_x_n_sd,1);
    fDefGradInv_column_matrix_Transpose.Dimension (1,n_sd_x_n_sd);
    u_dotdot_column_matrix.Dimension (n_en_displ_x_n_sd,1);
    fGradv_vector.Dimension (n_sd_x_n_sd);
    fgradv_vector.Dimension (n_sd_x_n_sd);
    fXi_temp_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fVarsigma_temp_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    fI_ijkl_matrix.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);
    u_dot_column_matrix.Dimension (n_en_displ_x_n_sd,1);
    fTemp_matrix1_ndof_se_x_ndof_se.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fK_dd_G4_matrix.Dimension (n_en_displ_x_n_sd,n_en_displ_x_n_sd);
    fGravity_column_matrix.Dimension (n_sd, 1);
    fTemp_matrix_nsd_x_ndof_se.Dimension (n_sd,n_en_displ_x_n_sd);
    fTemp_matrix_nsd_x_nen_micro.Dimension (n_sd,n_en_micro);
    micro_dot_column_matrix.Dimension (n_en_micro,1);
    u_dot_column_matrix_Transpose.Dimension (1, n_en_displ_x_n_sd);
    fImath_temp_matrix.Dimension (n_sd,n_sd_x_n_sd);

    /* streams */
    ofstreamT& out = ElementSupport().Output();

    /* storage for integration point strain, stress, and ISVs*/
    fIPVariable.Dimension (n_el, fNumIP_displ*(knumstrain+knumstress+knum_d_state+knumdispl));
    //fIPVariable.Dimension (n_el, fNumIP_micro*(knumstrain+knumstress+knum_d_state));
    fIPVariable = 0.0;

    /* allocate storage for nodal forces */
    //fForces_at_Node.Dimension ( n_sd );

    /* extract natural boundary conditions */
    TakeNaturalBC(list);

    /* setup output file and format */
    outputPrecision = 10;
    outputFileWidth = outputPrecision + 8;
    fs_micromorph3D_out.open("fs_micromorph3D.info");
//   fs_micromorph3D_outMn.open("fs_micromorph3DMn.info");
}


/* information about subordinate parameter lists */
void FSMicromorphic3DT::DefineSubs(SubListT& sub_list) const
{
    /* inherited */
    ElementBaseT::DefineSubs(sub_list);

    /* element blocks */
    sub_list.AddSub("micromorphic_FS_3D_element_block");

    /* tractions */
    sub_list.AddSub("micromorphic_FS_3D_natural_bc", ParameterListT::Any);
}



/* return the description of the given inline subordinate parameter list */
void FSMicromorphic3DT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order,
                       SubListT& sub_lists) const
{
    ElementBaseT::DefineInlineSub(name, order, sub_lists);
}



/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FSMicromorphic3DT::NewSub(const StringT& name) const
{
    /* create non-const this */
    FSMicromorphic3DT* non_const_this = const_cast<FSMicromorphic3DT*>(this);

    if (name == "micromorphic_FS_3D_natural_bc") /* traction bc */
    {
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
    else if (name == "micromorphic_FS_3D_element_block")
    {
        ParameterContainerT* element_block = new ParameterContainerT(name);
        element_block->AddSub("block_ID_list");
        return element_block;
    }
    else /* inherited */
    return ElementBaseT::NewSub(name);
}



//##################################################################################
//###### Traction B.C. Methods (Cut and Paste from ContinuumElementT) ##############
//##################################################################################

//---------------------------------------------------------------------

//---------------------------------------------------------------------

/* update traction BC data */
void FSMicromorphic3DT::SetTractionBC(void)
{
//NOTE: With the possibility of variable global node numbers and
//      and equations, we assume as little as possible here with
//      regard to the validity of the node/equation numbers, requiring
//      only that NodesX in the element cards has the correct global
//      node numbers.

    /* dimensions */
    int ndof = NumDOF();

    /* echo values */
    iArray2DT nd_tmp, eq_tmp;
    for (int i = 0; i < fTractionList.Length(); i++)
    {
    Traction_CardT& BC_card = fTractionList[i];

    /* traction element/facet */
    int elem, facet;
    BC_card.Destination(elem, facet);

    /* set global node numbers */
    const iArrayT& loc_nodes = BC_card.LocalNodeNumbers();
    int nnd = loc_nodes.Length();

    iArrayT& nodes = BC_card.Nodes();
    nodes.Dimension(nnd);
    nodes.Collect(loc_nodes, fElementCards[elem].NodesX());

    /* set global equation numbers */
    iArrayT& eqnos = BC_card.Eqnos();
    eqnos.Dimension(ndof*nnd);

    /* get from node manager */
    nd_tmp.Set(1, nnd, nodes.Pointer());
    eq_tmp.Set(1, ndof*nnd, eqnos.Pointer());
    fDispl->SetLocalEqnos(nd_tmp, eq_tmp);
    }

    /* set flag */
    fTractionBCSet = 1;
}



/* extract natural boundary condition information */
void FSMicromorphic3DT::TakeNaturalBC(const ParameterListT& list)
{
    const char caller[] = "FSMicromorphic3DT::TakeTractionBC";

    int num_natural_bc = list.NumLists("natural_bc");
    if (num_natural_bc > 0)
    {
    /* model manager */
    ModelManagerT& model = ElementSupport().ModelManager();

    /* temp space */
    ArrayT<StringT> block_ID(num_natural_bc);
    ArrayT<iArray2DT> localsides(num_natural_bc);
    iArrayT LTf(num_natural_bc);
    ArrayT<Traction_CardT::CoordSystemT> coord_sys(num_natural_bc);
    ArrayT<dArray2DT> values(num_natural_bc);

    /* nodes on element facets */
    iArrayT num_facet_nodes;
    fShapes_displ->NumNodesOnFacets(num_facet_nodes);

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

        /* switch to elements numbering within the group */
        iArray2DT& side_set = localsides[i];
        iArrayT elems(num_sides);
        side_set.ColumnCopy(0, elems);
        BlockToGroupElementNumbers(elems, block_ID[i]);
        side_set.SetColumn(0, elems);

        /* all facets in set must have the same number of nodes */
        int num_nodes = num_facet_nodes[side_set(0,1)];
        for (int f = 0; f < num_sides; f++)
            if (num_facet_nodes[side_set(f,1)] != num_nodes)
            ExceptionT::BadInputValue(caller, "faces side set \"%s\" have different numbers of nodes",
                          ss_ID.Pointer());

        /* read traction nodal values */
        dArray2DT& nodal_values = values[i];
        nodal_values.Dimension(num_nodes, NumDOF());
        int num_traction_vectors = natural_bc.NumLists("DoubleList");
        if (num_traction_vectors != 1 && num_traction_vectors != num_nodes)
            ExceptionT::GeneralFail(caller, "expecting 1 or %d vectors not %d",
                        num_nodes, num_traction_vectors);

        /* constant over the face */
        if (num_traction_vectors == 1) {
            const ParameterListT& traction_vector = natural_bc.GetList("DoubleList");
            int dim = traction_vector.NumLists("Double");
            if (dim != NumDOF())
            ExceptionT::GeneralFail(caller, "expecting traction vector length %d not %d",
                        NumDOF(), dim);

            /* same for all face nodes */
            for (int f = 0; f < NumDOF(); f++) {
            double t = traction_vector.GetList("Double", f).GetParameter("value");
            nodal_values.SetColumn(f, t);
            }
        }
        else
        {
            /* read separate vector for each face node */
            dArrayT t;
            for (int f = 0; f < num_nodes; f++) {
            const ParameterListT& traction_vector = natural_bc.GetList("DoubleList", f);
            int dim = traction_vector.NumLists("Double");
            if (dim != NumDOF())
                ExceptionT::GeneralFail(caller, "expecting traction vector length %d not %d",
                            NumDOF(), dim);

            nodal_values.RowAlias(f, t);
            for (int j = 0; j < NumDOF(); j++)
                t[j] = traction_vector.GetList("Double", j).GetParameter("value");
            }
        }
        }
    }
#pragma message("OK with empty side sets?")

    /* allocate all traction BC cards */
    fTractionList.Dimension(tot_num_sides);

    /* correct numbering offset */
    LTf--;

    /* define traction cards */
    if (tot_num_sides > 0)
    {
        iArrayT loc_node_nums;
        int dex = 0;
        for (int i = 0; i < num_natural_bc; i++)
        {
        /* set traction BC cards */
        iArray2DT& side_set = localsides[i];
        int num_sides = side_set.MajorDim();
        for (int j = 0; j < num_sides; j++)
        {
            /* get facet local node numbers */
            fShapes_displ->NodesOnFacet(side_set(j, 1), loc_node_nums);

            /* set and echo */
            fTractionList[dex++].SetValues(ElementSupport(), side_set(j,0), side_set (j,1), LTf[i],
                           coord_sys[i], loc_node_nums, values[i]);
        }
        }
    }

    /* check coordinate system specifications */
    if (NumSD() != NumDOF())
        for (int i = 0; i < fTractionList.Length(); i++)
        if (fTractionList[i].CoordSystem() != Traction_CardT::kCartesian)
            ExceptionT::BadInputValue(caller, "coordinate system must be Cartesian if (nsd != ndof) for card %d", i+1);
    }
}


//---------------------------------------------------------------------

/* compute contribution to RHS from traction BC's */
void FSMicromorphic3DT::ApplyTractionBC(void)
{
    if (fTractionList.Length() > 0)
    {
    /* dimensions */
    int nsd = NumSD();
    int ndof = NumDOF();

    /* update equation numbers */
    if (!fTractionBCSet) SetTractionBC();

    /* force vector */
    dArrayT rhs;
    VariArrayT<double> rhs_man(25, rhs);

    /* local coordinates */
    LocalArrayT coords(LocalArrayT::kInitCoords);
    VariLocalArrayT coord_man(25, coords, nsd);
    ElementSupport().RegisterCoordinates(coords);

    /* nodal tractions */
    LocalArrayT tract(LocalArrayT::kUnspecified);
    VariLocalArrayT tract_man(25, tract, ndof);

    /* integration point tractions */
    dArray2DT ip_tract;
    nVariArray2DT<double> ip_tract_man(25, ip_tract, ndof);
    dArrayT tract_loc, tract_glb(ndof);
    dMatrixT Q(ndof);

    /* Jacobian of the surface mapping */
    dMatrixT jacobian(nsd, nsd-1);

    for (int i = 0; i < fTractionList.Length(); i++)
    {
        const Traction_CardT& BC_card = fTractionList[i];

        /* dimension */
        const iArrayT& nodes = BC_card.Nodes();
        int nnd = nodes.Length();
        rhs_man.SetLength(nnd*ndof, false);
        coord_man.SetNumberOfNodes(nnd);
        tract_man.SetNumberOfNodes(nnd);

        /* local coordinates */
        coords.SetLocal(nodes);

        /* nodal traction vectors: (ndof x nnd) */
        BC_card.CurrentValue(tract);

        /* BC destination */
        int elem, facet;
        BC_card.Destination(elem, facet);

        /* default thickness */
        double thick = 1.0;

        /* boundary shape functions */
        const ParentDomainT& surf_shape = ShapeFunctionDispl().FacetShapeFunction(facet);
        int nip = surf_shape.NumIP();

        /* all ip tractions: (nip x ndof) */
        ip_tract_man.SetMajorDimension(nip, false);
        surf_shape.Interpolate(tract, ip_tract);

        /* traction vector coordinate system */
        if (BC_card.CoordSystem() == Traction_CardT::kCartesian)
        {
        /* integrate */
        rhs = 0.0;
        const double* w = surf_shape.Weight();
        for (int j = 0; j < nip; j++)
        {
            /* coordinate mapping */
            surf_shape.DomainJacobian(coords, j, jacobian);
            double detj = surf_shape.SurfaceJacobian(jacobian);

            /* ip weight */
            double jwt = detj*w[j]*thick;

            /* ip traction */
            const double* tj = ip_tract(j);

            /* accumulate */
            for (int l = 0; l < ndof; l++)
            {
            /* nodal shape function */
            const double* Na = surf_shape.Shape(j);

            double* prhs = rhs.Pointer(l);
            double  fact = jwt*(*tj++);
            for (int k = 0; k < nnd; k++)
            {
                *prhs += fact*(*Na++);
                prhs += ndof;
            }
            }
        }
        }
        else if (BC_card.CoordSystem() == Traction_CardT::kLocal)
        {
        /* integrate */
        rhs = 0.0;
        const double* w = surf_shape.Weight();
        for (int j = 0; j < nip; j++)
        {
            /* coordinate mapping */
            surf_shape.DomainJacobian(coords, j, jacobian);
            double detj = surf_shape.SurfaceJacobian(jacobian, Q);

            /* ip weight */
            double jwt = detj*w[j]*thick;

            /* transform ip traction out of local frame */
            ip_tract.RowAlias(j, tract_loc);
            Q.Multx(tract_loc, tract_glb);

            /* ip traction */
            const double* tj = tract_glb.Pointer();

            /* accumulate */
            for (int l = 0; l < ndof; l++)
            {
            /* nodal shape function */
            const double* Na = surf_shape.Shape(j);

            double* prhs = rhs.Pointer(l);
            double  fact = jwt*(*tj++);
            for (int k = 0; k < nnd; k++)
            {
                *prhs += fact*(*Na++);
                prhs += ndof;
            }
            }
        }
        }
        else
        throw ExceptionT::kGeneralFail;

        /* assemble into displacement equations */
        ElementSupport().AssembleRHS(fDispl->Group(), rhs, BC_card.Eqnos());
    }
    }
}

void FSMicromorphic3DT::Form_solid_shape_functions(const double* &shapes_displ_X)
{
    fShapeDispl = 0.0;
    for (int i=0; i<27; i++)
    {
        fShapeDispl(0,i*3) = shapes_displ_X[i];
        fShapeDispl(1,1+i*3) = shapes_displ_X[i];
        fShapeDispl(2,2+i*3) = shapes_displ_X[i];
    }
}

void FSMicromorphic3DT::Form_Gradient_of_solid_shape_functions(const dMatrixT &fShapeDisplGrad_temp)
{
    fShapeDisplGrad = 0.0;
    for(int i=0; i<27; i++)
    {
        fShapeDisplGrad(0,i*3) = fShapeDisplGrad_temp(0,i);
        fShapeDisplGrad(1,1+i*3) = fShapeDisplGrad_temp(0,i);
        fShapeDisplGrad(2,2+i*3) = fShapeDisplGrad_temp(0,i);

        fShapeDisplGrad(3,i*3) = fShapeDisplGrad_temp(1,i);
        fShapeDisplGrad(4,1+i*3) = fShapeDisplGrad_temp(1,i);
        fShapeDisplGrad(5,2+i*3) = fShapeDisplGrad_temp(1,i);

        fShapeDisplGrad(6,i*3) = fShapeDisplGrad_temp(2,i);
        fShapeDisplGrad(7,1+i*3) = fShapeDisplGrad_temp(2,i);
        fShapeDisplGrad(8,2+i*3) = fShapeDisplGrad_temp(2,i);
    }

}

void FSMicromorphic3DT::Form_micro_shape_functions(const double* &shapes_micro_X)
{
    fShapeMicro = 0.0;
    //hard coded for n_en_micro=8; can change
    for (int i=0; i<n_en_micro; i++)
        fShapeMicro[i] = shapes_micro_X[i];
}

void FSMicromorphic3DT::Form_deformation_gradient_tensor(void)
{
    fShapeDisplGrad.Multx(u_vec,fGrad_disp_vector);
    fDeformation_Gradient(0,0) = fGrad_disp_vector[0]+1.0;
    fDeformation_Gradient(0,1) = fGrad_disp_vector[3];
    fDeformation_Gradient(0,2) = fGrad_disp_vector[6];
    fDeformation_Gradient(1,0) = fGrad_disp_vector[1];
    fDeformation_Gradient(1,1) = fGrad_disp_vector[4]+1.0;
    fDeformation_Gradient(1,2) = fGrad_disp_vector[7];
    fDeformation_Gradient(2,0) = fGrad_disp_vector[2];
    fDeformation_Gradient(2,1) = fGrad_disp_vector[5];
    fDeformation_Gradient(2,2) = fGrad_disp_vector[8]+1.0;

    fGrad_disp_matrix(0,0) = fGrad_disp_vector[0];
    fGrad_disp_matrix(0,1) = fGrad_disp_vector[3];
    fGrad_disp_matrix(0,2) = fGrad_disp_vector[6];
    fGrad_disp_matrix(1,0) = fGrad_disp_vector[1];
    fGrad_disp_matrix(1,1) = fGrad_disp_vector[4];
    fGrad_disp_matrix(1,2) = fGrad_disp_vector[7];
    fGrad_disp_matrix(2,0) = fGrad_disp_vector[2];
    fGrad_disp_matrix(2,1) = fGrad_disp_vector[5];
    fGrad_disp_matrix(2,2) = fGrad_disp_vector[8];

}

void FSMicromorphic3DT::Form_Grad_grad_transformation_matrix(void)
{
    fDefGradInv_Grad_grad = 0.0;
    fDefGradInv_Grad_grad(0,0) = fDeformation_Gradient_Inverse(0,0);
    fDefGradInv_Grad_grad(0,3) = fDeformation_Gradient_Inverse(0,1);
    fDefGradInv_Grad_grad(0,6) = fDeformation_Gradient_Inverse(0,2);
    fDefGradInv_Grad_grad(1,1) = fDeformation_Gradient_Inverse(0,0);
    fDefGradInv_Grad_grad(1,4) = fDeformation_Gradient_Inverse(0,1);
    fDefGradInv_Grad_grad(1,7) = fDeformation_Gradient_Inverse(0,2);
    fDefGradInv_Grad_grad(2,2) = fDeformation_Gradient_Inverse(0,0);
    fDefGradInv_Grad_grad(2,5) = fDeformation_Gradient_Inverse(0,1);
    fDefGradInv_Grad_grad(2,8) = fDeformation_Gradient_Inverse(0,2);

    fDefGradInv_Grad_grad(3,0) = fDeformation_Gradient_Inverse(1,0);
    fDefGradInv_Grad_grad(3,3) = fDeformation_Gradient_Inverse(1,1);
    fDefGradInv_Grad_grad(3,6) = fDeformation_Gradient_Inverse(1,2);
    fDefGradInv_Grad_grad(4,1) = fDeformation_Gradient_Inverse(1,0);
    fDefGradInv_Grad_grad(4,4) = fDeformation_Gradient_Inverse(1,1);
    fDefGradInv_Grad_grad(4,7) = fDeformation_Gradient_Inverse(1,2);
    fDefGradInv_Grad_grad(5,2) = fDeformation_Gradient_Inverse(1,0);
    fDefGradInv_Grad_grad(5,5) = fDeformation_Gradient_Inverse(1,1);
    fDefGradInv_Grad_grad(5,8) = fDeformation_Gradient_Inverse(1,2);

    fDefGradInv_Grad_grad(6,0) = fDeformation_Gradient_Inverse(2,0);
    fDefGradInv_Grad_grad(6,3) = fDeformation_Gradient_Inverse(2,1);
    fDefGradInv_Grad_grad(6,6) = fDeformation_Gradient_Inverse(2,2);
    fDefGradInv_Grad_grad(7,1) = fDeformation_Gradient_Inverse(2,0);
    fDefGradInv_Grad_grad(7,4) = fDeformation_Gradient_Inverse(2,1);
    fDefGradInv_Grad_grad(7,7) = fDeformation_Gradient_Inverse(2,2);
    fDefGradInv_Grad_grad(8,2) = fDeformation_Gradient_Inverse(2,0);
    fDefGradInv_Grad_grad(8,5) = fDeformation_Gradient_Inverse(2,1);
    fDefGradInv_Grad_grad(8,8) = fDeformation_Gradient_Inverse(2,2);
}

void FSMicromorphic3DT::Form_fDefGradT_9x9_matrix(void)
{
    fDefGradT_9x9_matrix = 0.0;
    fDefGradT_9x9_matrix(0,0) = fDeformation_Gradient(0,0);
    fDefGradT_9x9_matrix(0,1) = fDeformation_Gradient(1,0);
    fDefGradT_9x9_matrix(0,2) = fDeformation_Gradient(2,0);
    fDefGradT_9x9_matrix(1,3) = fDeformation_Gradient(0,0);
    fDefGradT_9x9_matrix(1,4) = fDeformation_Gradient(1,0);
    fDefGradT_9x9_matrix(1,5) = fDeformation_Gradient(2,0);
    fDefGradT_9x9_matrix(2,6) = fDeformation_Gradient(0,0);
    fDefGradT_9x9_matrix(2,7) = fDeformation_Gradient(1,0);
    fDefGradT_9x9_matrix(2,8) = fDeformation_Gradient(2,0);

    fDefGradT_9x9_matrix(3,0) = fDeformation_Gradient(0,1);
    fDefGradT_9x9_matrix(3,1) = fDeformation_Gradient(1,1);
    fDefGradT_9x9_matrix(3,2) = fDeformation_Gradient(2,1);
    fDefGradT_9x9_matrix(4,3) = fDeformation_Gradient(0,1);
    fDefGradT_9x9_matrix(4,4) = fDeformation_Gradient(1,1);
    fDefGradT_9x9_matrix(4,5) = fDeformation_Gradient(2,1);
    fDefGradT_9x9_matrix(5,6) = fDeformation_Gradient(0,1);
    fDefGradT_9x9_matrix(5,7) = fDeformation_Gradient(1,1);
    fDefGradT_9x9_matrix(5,8) = fDeformation_Gradient(2,1);

    fDefGradT_9x9_matrix(6,0) = fDeformation_Gradient(0,2);
    fDefGradT_9x9_matrix(6,1) = fDeformation_Gradient(1,2);
    fDefGradT_9x9_matrix(6,2) = fDeformation_Gradient(2,2);
    fDefGradT_9x9_matrix(7,3) = fDeformation_Gradient(0,2);
    fDefGradT_9x9_matrix(7,4) = fDeformation_Gradient(1,2);
    fDefGradT_9x9_matrix(7,5) = fDeformation_Gradient(2,2);
    fDefGradT_9x9_matrix(8,6) = fDeformation_Gradient(0,2);
    fDefGradT_9x9_matrix(8,7) = fDeformation_Gradient(1,2);
    fDefGradT_9x9_matrix(8,8) = fDeformation_Gradient(2,2);
}


void FSMicromorphic3DT::Form_deformation_gradient_inv_vector(void)
{
    fDefGradInv_vector[0] = fDeformation_Gradient_Inverse(0,0);
    fDefGradInv_vector[1] = fDeformation_Gradient_Inverse(0,1);
    fDefGradInv_vector[2] = fDeformation_Gradient_Inverse(0,2);
    fDefGradInv_vector[3] = fDeformation_Gradient_Inverse(1,0);
    fDefGradInv_vector[4] = fDeformation_Gradient_Inverse(1,1);
    fDefGradInv_vector[5] = fDeformation_Gradient_Inverse(1,2);
    fDefGradInv_vector[6] = fDeformation_Gradient_Inverse(2,0);
    fDefGradInv_vector[7] = fDeformation_Gradient_Inverse(2,1);
    fDefGradInv_vector[8] = fDeformation_Gradient_Inverse(2,2);

}

void FSMicromorphic3DT::Form_kirchhoff_stress_vector()
{
    fKirchhoff_vector[0] = fKirchhoff_tensor(0,0);
    fKirchhoff_vector[1] = fKirchhoff_tensor(1,0);
    fKirchhoff_vector[2] = fKirchhoff_tensor(2,0);
    fKirchhoff_vector[3] = fKirchhoff_tensor(0,1);
    fKirchhoff_vector[4] = fKirchhoff_tensor(1,1);
    fKirchhoff_vector[5] = fKirchhoff_tensor(2,1);
    fKirchhoff_vector[6] = fKirchhoff_tensor(0,2);
    fKirchhoff_vector[7] = fKirchhoff_tensor(1,2);
    fKirchhoff_vector[8] = fKirchhoff_tensor(2,2);
}
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
////////////////////// MATRICES FOR MICROMORPHIC 3-D CASE/////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

void FSMicromorphic3DT:: Form_CCof_tensor()
{

       for(int k=0;k<3;k++)
        {
            for(int l=0;l<3;l++)
            {
                for(int m=0;m<3;m++)
                {
                    for(int n=0;n<3;n++)
                    {
                        for(int p=0;p<3;p++)
                        {
                            for(int q=0;q<3;q++)
                            {
                            CCof[k][l][m][n][p][q]=0.0;
                            }}}}}}

    for(int k=0;k<3;k++)
    {
        for(int l=0;l<3;l++)
        {
            for(int m=0;m<3;m++)
            {
                for(int n=0;n<3;n++)
                {
                    for(int p=0;p<3;p++)
                    {
                        for(int q=0;q<3;q++)
                        {
                        CCof[k][l][m][n][p][q]=   fMaterial_Params[kTau1]*( KrDelta[k][l]*KrDelta[n][m]*KrDelta[q][p]+ KrDelta[l][m]*KrDelta[n][p]*KrDelta[k][q])
                                                + fMaterial_Params[kTau2]*( KrDelta[k][l]*KrDelta[m][p]*KrDelta[n][q]+ KrDelta[k][m]*KrDelta[l][q]*KrDelta[n][p])
                                                + fMaterial_Params[kTau3]*( KrDelta[k][l]*KrDelta[q][m]*KrDelta[n][p])
                                                + fMaterial_Params[kTau4]*( KrDelta[l][m]*KrDelta[k][n]*KrDelta[q][p])
                                                + fMaterial_Params[kTau5]*( KrDelta[k][m]*KrDelta[l][n]*KrDelta[q][p]+ KrDelta[l][m]*KrDelta[k][p]*KrDelta[n][q])
                                                + fMaterial_Params[kTau6]*( KrDelta[k][m]*KrDelta[l][p]*KrDelta[n][q])
                                                + fMaterial_Params[kTau7]*( KrDelta[k][n]*KrDelta[l][p]*KrDelta[q][m])
                                                + fMaterial_Params[kTau8]*( KrDelta[k][p]*KrDelta[l][q]*KrDelta[n][m]+ KrDelta[k][q]*KrDelta[l][n]*KrDelta[m][p])
                                                + fMaterial_Params[kTau9]*( KrDelta[k][n]*KrDelta[l][q]*KrDelta[m][p])
                                                +fMaterial_Params[kTau10]*( KrDelta[k][p]*KrDelta[l][n]*KrDelta[q][m])
                                                +fMaterial_Params[kTau11]*( KrDelta[k][q]*KrDelta[l][p]*KrDelta[n][m]);
                        }
                    }
                }
            }
        }
    }


}

void FSMicromorphic3DT::Form_micro_deformation_tensor_Chi()
{
    NCHI.Multx(Phi_vec,Chi_vec);
    Chi[0][0] = Chi_vec[0]+1.0;
    Chi[0][1] = Chi_vec[3];
    Chi[0][2] = Chi_vec[6];
    Chi[1][0] = Chi_vec[1];
    Chi[1][1] = Chi_vec[4]+1.0;
    Chi[1][2] = Chi_vec[7];
    Chi[2][0] = Chi_vec[2];
    Chi[2][1] = Chi_vec[5];
    Chi[2][2] = Chi_vec[8]+1.0;

    for(int i=0;i<3;i++)
    {
    	for(int j=0;j<3;j++)
    	{
    		ChiM(i,j)=Chi[i][j];
    	}
    }

}

void FSMicromorphic3DT:: Form_Chi_inv_matrix()
{
    ChiInv_m=0.0;
    Chi_m=0.0;
    for(int i=0;i<3;i++)
    {for(int j=0;j<3;j++)
    {Chi_m(i,j)=Chi[i][j];}}

    ChiInv_m.Inverse(Chi_m);
    for(int i=0;i<3;i++)
    {for(int j=0;j<3;j++)
    {ChiInv[i][j]=ChiInv_m(i,j);}}


}
void FSMicromorphic3DT:: Form_GRAD_Chi_matrix()
{
    int row;
    row=0;
    GRAD_NCHI.Multx(Phi_vec,GRAD_Chi_vec);
    for(int T=0;T<=2;T++)
    {
        for(int i=0;i<=2;i++)
        {
            for(int K=0;K<=2;K++)
            {
                GRAD_Chi[i][T][K]=GRAD_Chi_vec[row];//CHI=1+PHI ==> GRAD_CHI=GRAD_PHI
                row++;
            }
        }
    }

    row=0;
    for(int T=0;T<=2;T++)
    {
        for(int i=0;i<=2;i++)
        {
            for(int K=0;K<=2;K++)
            {
                GRAD_CHIM(i,T,K)=GRAD_Chi_vec[row];//CHI=1+PHI ==> GRAD_CHI=GRAD_PHI
                row++;
            }
        }
    }
}


//Forming the Matrices coming from the Bal. of Lin. Mom.
void FSMicromorphic3DT::Form_KroneckerDelta_matrix()
{
    for(int i=0;i<=2;i++)
    { for(int j=0;j<=2;j++)
            {KrDelta[i][j]=0.0;}}
    KrDelta[0][0]=1.0;
    KrDelta[1][1]=1.0;
    KrDelta[2][2]=1.0;
}

void FSMicromorphic3DT::Form_Gamma_tensor3D()
{

    for(int a=0;a<3;a++)
    {for(int p=0;p<3;p++)
       {for(int q=0;q<3;q++)
            {Gamma[a][p][q]=0.0;}}}

    for(int a=0;a<3;a++)
    {
        for(int p=0;p<3;p++)
        {
            for(int q=0;q<3;q++)
            {
                //summation over the same term starts here
                for(int A=0;A<3;A++)
                    {
                    for(int Q=0;Q<3;Q++)
                     {
                        Gamma[a][p][q]+=ChiInv[A][p]*GRAD_Chi[a][A][Q]*Finv[Q][q];
                     }
                    }
            }
        }
    }
}

void FSMicromorphic3DT::Form_G1_matrix()
{
    int row;
    int col;
    double scale;
    scale=0.0;
    row=0;
/*    Mat1=0.0;
    Mat2=0.0;
    Mat3=0.0;
    Mat4=0.0;
    Mat5=0.0;*/
    Sigma=0.0;
  //  RHS=0.0;
    //Sigma1=0.0;
    G1=0.0;
    deltaL=0.0;
    deltad=0.0;
    deltaL_Tr=0.0;
    Fn_m=0.0;
    Finv_m=0.0;
    SigN_m=0.0;
    ChiInv_m=0.0;
    ChiN_m=0.0;
    deltaEp=0.0;
    deltaNu=0.0;
    tempSig=0.0;
    double trdeltaEp=0.0;
    double trdeltad=0.0;


    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            Fn_m(i,j)=Fn[i][j];
            Finv_m(i,j)=Finv[i][j];
            SigN_m(i,j)=SigN[i][j];
            ChiInv_m(i,j)=ChiInv[i][j];
            ChiN_m(i,j)=ChiN[i][j];
          }
    }



deltaL.MultAB(Fn_m,Finv_m);
deltaL*=-1;
deltaL+=fIdentity_matrix;

deltaL_Tr.Transpose(deltaL);

deltad=deltaL;
deltad+=deltaL_Tr;
deltad*=0.5;
trdeltad=0.0;
for(int i=0;i<3;i++)
    trdeltad+=deltad(i,i);
//cout<<trdeltad<<endl;

deltaNu.MultAB(ChiN_m,ChiInv_m);
deltaNu*=-1;
deltaNu+=fIdentity_matrix;
deltaEp=deltaNu;
deltaEp+=deltaL_Tr;

for(int i=0;i<3;i++)
    trdeltaEp+=deltaEp(i,i);


Sigma=SigN_m;
//


tempSig=SigN_m;
scale=-1*trdeltad;
tempSig*=scale;
Sigma+=tempSig;
//

tempSig.MultAB(deltaL,SigN_m);
Sigma+=tempSig;


tempSig.MultABT(SigN_m,deltaL);
Sigma+=tempSig;
//


tempSig=fIdentity_matrix;
scale=trdeltad*(fMaterial_Params[kLambda]+fMaterial_Params[kTau]);
tempSig*=scale;
Sigma+=tempSig;
//

tempSig=deltad;
scale=2*(fMaterial_Params[kMu]+fMaterial_Params[kSigma_const]);
tempSig*=scale;
Sigma+=tempSig;

//

tempSig=fIdentity_matrix;
scale=trdeltaEp*fMaterial_Params[kEta];
tempSig*=scale;
Sigma+=tempSig;


tempSig=deltaEp;
scale=fMaterial_Params[kKappa];
tempSig*=scale;
Sigma+=tempSig;
//


tempSig.Transpose(deltaEp);
scale=fMaterial_Params[kNu];
tempSig*=scale;
Sigma+=tempSig;


///// to calculate fn+1/////
/*for(int l=0;l<9;l++)
{
     Mat1(l,l)=1.0;
}

for(int l=0;l<9;l++)
{

     Mat2(l,l)=trdeltad;

}

row=0;
col=0;
for(int k=0;k<3;k++)
    {
        Mat3(row,col)    =-deltaL(0,0);
        Mat3(row,col+1)  =-deltaL(0,1);
        Mat3(row,col+2)  =-deltaL(0,2);
        Mat3(row+1,col)  =-deltaL(1,0);
        Mat3(row+1,col+1)=-deltaL(1,1);
        Mat3(row+1,col+2)=-deltaL(1,2);
        Mat3(row+2,col)  =-deltaL(2,0);
        Mat3(row+2,col+1)=-deltaL(2,1);
        Mat3(row+2,col+2)=-deltaL(2,2);
        col=col+3;
        row=row+3;
    }

row=0;
col=0;

for(int l=0;l<3;l++)
{
    row=l*3;
    col=0;
    for(int k=0;k<3;k++)
    {
        Mat4(row  ,col)  =-deltaL(l,k);
        Mat4(row+1,col+1)=-deltaL(l,k);
        Mat4(row+2,col+2)=-deltaL(l,k);
        col=col+3;
    }
}

Mat5=Mat1;
Mat5+=Mat2;
Mat5+=Mat3;
Mat5+=Mat4;
Mat5_Inv.Inverse(Mat5);

RHS=0.0;
row=0;

for(int l=0;l<3;l++)
{
    for(int k=0;k<3;k++)
    {
        RHS[row]=SigN[k][l]+fMaterial_Params[kLambda]*trdeltad*fIdentity_matrix(l,k)+2*fMaterial_Params[kMu]*deltad(k,l);
        row++;
    }
}

Mat5_Inv.Multx(RHS,Sigma1);





row=0;
for(int k=0;k<=2;k++)
{
    for(int l=0;l<=2;l++)
    {
        Sigma(l,k)=Sigma1[row];
        G1[row]=Sigma1[row];
        row++;
    }
}*/


    for(int k=0;k<=2;k++)
    {
        for(int l=0;l<=2;l++)
        {
             G1[row]=Sigma(l,k);
            row++;
        }
    }

}

void FSMicromorphic3DT::Form_Finv_w_matrix()//checked correct
{

    int count, row, col;
    row=0;
    col=0;
    count=0;
    Finv_w=0.0;
    while(count<=2)
    {
        Finv_w(row,col)    =fDeformation_Gradient_Inverse(0,0);
        Finv_w(row,col+1)  =fDeformation_Gradient_Inverse(1,0);
        Finv_w(row,col+2)  =fDeformation_Gradient_Inverse(2,0);
        Finv_w(row+1,col)  =fDeformation_Gradient_Inverse(0,1);
        Finv_w(row+1,col+1)=fDeformation_Gradient_Inverse(1,1);
        Finv_w(row+1,col+2)=fDeformation_Gradient_Inverse(2,1);
        Finv_w(row+2,col)  =fDeformation_Gradient_Inverse(0,2);
        Finv_w(row+2,col+1)=fDeformation_Gradient_Inverse(1,2);
        Finv_w(row+2,col+2)=fDeformation_Gradient_Inverse(2,2);
        row=row+3;
        col=col+3;
        count++;
    }

}



void FSMicromorphic3DT::Form_NCHI_matrix(const dMatrixT &fShapeMicro_row_matrix)
{
    int row=0;
    int col=0;
  //  int counter;
    NCHI=0.0;
    for(int j=0;j<=8;j++)
    {
        col=j;
        for(int i=0;i<8;i++)
        {
            NCHI(row,col)=fShapeMicro_row_matrix(0,i);
            col=col+9;
        }
        row++;
    }

}


void FSMicromorphic3DT:: Form_GRAD_Nuw_matrix(const dMatrixT &fShapeDisplGrad_temp)
{
    int row,col;
    row=0;
    col=0;
    GRAD_Nuw=0.0;
    for(int j=0;j<3;j++)
    {
        col=j;
        for(int i=0;i<27;i++)
        {
            GRAD_Nuw(row,col)  =fShapeDisplGrad_temp(0,i);
            GRAD_Nuw(row+1,col)=fShapeDisplGrad_temp(1,i);
            GRAD_Nuw(row+2,col)=fShapeDisplGrad_temp(2,i);
            col=col+3;
        }
        row=row+3;
    }


}


void FSMicromorphic3DT::Form_Var_F_tensor()
{
    int row,col;
    row=0;
    col=0;
    Var_F=0.0;
/*    for (int l=0;l<3;l++)
    {
        for(int i=0;i<3;i++)
        {
            row=0;
            //row operations starts
            for (int k=0;k<3;k++)
            {
                    //summation on the same terms starts
                    for (int K=0;K<3;K++)
                    {
                        Var_F(row,col)=Finv[K][i]*Sigma(l,k);
                       row++;
                    }

                }
            col++;
        }
    }*/
  for (int l=0;l<3;l++)
    {
        for(int m=0;m<3;m++)
        {
            row=m;
            for (int k=0;k<3;k++)
                 {
                  Var_F(row,col)= Sigma(l,k);
                   row=row+3;
                }
            col++;
        }
    }



}

void FSMicromorphic3DT::Form_Tsigma_1_matrix()
{


        int row=0;
        int col=0;
        Tsigma_1=0.0;
        for(int i=0;i<3;i++)
        {
            for(int m=0;m<3;m++)
            {
                row=0;
                for(int k=0;k<3;k++)
                {
                    for(int l=0;l<3;l++)
                    {
                        //summation on the same term
                        for(int L=0;L<3;L++)
                        {
                            Tsigma_1(row,col)+=SigN[l][k]*Fn[i][L]*Finv[L][m];
                        }
                        row++;
                    }
                }
            col++;
            }
        }
/*    Tsigma_1=0.0;
    Tsigma_1(0,0)=(Finv[0][0]*Fn[0][0]*SigN[0][0] + Finv[1][0]*Fn[0][1]*SigN[0][0] + Finv[2][0]*Fn[0][2]*SigN[0][0]);
    Tsigma_1(1,0)=(Finv[0][0]*Fn[0][0]*SigN[1][0] + Finv[1][0]*Fn[0][1]*SigN[1][0] + Finv[2][0]*Fn[0][2]*SigN[1][0]);
    Tsigma_1(2,0)=(Finv[0][0]*Fn[0][0]*SigN[2][0] + Finv[1][0]*Fn[0][1]*SigN[2][0] + Finv[2][0]*Fn[0][2]*SigN[2][0]);
    Tsigma_1(3,0)=(Finv[0][0]*Fn[0][0]*SigN[0][1] + Finv[1][0]*Fn[0][1]*SigN[0][1] + Finv[2][0]*Fn[0][2]*SigN[0][1]);
    Tsigma_1(4,0)=(Finv[0][0]*Fn[0][0]*SigN[1][1] + Finv[1][0]*Fn[0][1]*SigN[1][1] + Finv[2][0]*Fn[0][2]*SigN[1][1]);
    Tsigma_1(5,0)=(Finv[0][0]*Fn[0][0]*SigN[2][1] + Finv[1][0]*Fn[0][1]*SigN[2][1] + Finv[2][0]*Fn[0][2]*SigN[2][1]);
    Tsigma_1(6,0)=(Finv[0][0]*Fn[0][0]*SigN[0][2] + Finv[1][0]*Fn[0][1]*SigN[0][2] + Finv[2][0]*Fn[0][2]*SigN[0][2]);
    Tsigma_1(7,0)=(Finv[0][0]*Fn[0][0]*SigN[1][2] + Finv[1][0]*Fn[0][1]*SigN[1][2] + Finv[2][0]*Fn[0][2]*SigN[1][2]);
    Tsigma_1(8,0)=(Finv[0][0]*Fn[0][0]*SigN[2][2] + Finv[1][0]*Fn[0][1]*SigN[2][2] + Finv[2][0]*Fn[0][2]*SigN[2][2]);

    Tsigma_1(0,1)=(Finv[0][1]*Fn[0][0]*SigN[0][0] + Finv[1][1]*Fn[0][1]*SigN[0][0] + Finv[2][1]*Fn[0][2]*SigN[0][0]);
    Tsigma_1(1,1)=(Finv[0][1]*Fn[0][0]*SigN[1][0] + Finv[1][1]*Fn[0][1]*SigN[1][0] + Finv[2][1]*Fn[0][2]*SigN[1][0]);
    Tsigma_1(2,1)=(Finv[0][1]*Fn[0][0]*SigN[2][0] + Finv[1][1]*Fn[0][1]*SigN[2][0] + Finv[2][1]*Fn[0][2]*SigN[2][0]);
    Tsigma_1(3,1)=(Finv[0][1]*Fn[0][0]*SigN[0][1] + Finv[1][1]*Fn[0][1]*SigN[0][1] + Finv[2][1]*Fn[0][2]*SigN[0][1]);
    Tsigma_1(4,1)=(Finv[0][1]*Fn[0][0]*SigN[1][1] + Finv[1][1]*Fn[0][1]*SigN[1][1] + Finv[2][1]*Fn[0][2]*SigN[1][1]);
    Tsigma_1(5,1)=(Finv[0][1]*Fn[0][0]*SigN[2][1] + Finv[1][1]*Fn[0][1]*SigN[2][1] + Finv[2][1]*Fn[0][2]*SigN[2][1]);
    Tsigma_1(6,1)=(Finv[0][1]*Fn[0][0]*SigN[0][2] + Finv[1][1]*Fn[0][1]*SigN[0][2] + Finv[2][1]*Fn[0][2]*SigN[0][2]);
    Tsigma_1(7,1)=(Finv[0][1]*Fn[0][0]*SigN[1][2] + Finv[1][1]*Fn[0][1]*SigN[1][2] + Finv[2][1]*Fn[0][2]*SigN[1][2]);
    Tsigma_1(8,1)=(Finv[0][1]*Fn[0][0]*SigN[2][2] + Finv[1][1]*Fn[0][1]*SigN[2][2] + Finv[2][1]*Fn[0][2]*SigN[2][2]);

    Tsigma_1(0,2)=(Finv[0][2]*Fn[0][0]*SigN[0][0] + Finv[1][2]*Fn[0][1]*SigN[0][0] + Finv[2][2]*Fn[0][2]*SigN[0][0]);
    Tsigma_1(1,2)=(Finv[0][2]*Fn[0][0]*SigN[1][0] + Finv[1][2]*Fn[0][1]*SigN[1][0] + Finv[2][2]*Fn[0][2]*SigN[1][0]);
    Tsigma_1(2,2)=(Finv[0][2]*Fn[0][0]*SigN[2][0] + Finv[1][2]*Fn[0][1]*SigN[2][0] + Finv[2][2]*Fn[0][2]*SigN[2][0]);
    Tsigma_1(3,2)=(Finv[0][2]*Fn[0][0]*SigN[0][1] + Finv[1][2]*Fn[0][1]*SigN[0][1] + Finv[2][2]*Fn[0][2]*SigN[0][1]);
    Tsigma_1(4,2)=(Finv[0][2]*Fn[0][0]*SigN[1][1] + Finv[1][2]*Fn[0][1]*SigN[1][1] + Finv[2][2]*Fn[0][2]*SigN[1][1]);
    Tsigma_1(5,2)=(Finv[0][2]*Fn[0][0]*SigN[2][1] + Finv[1][2]*Fn[0][1]*SigN[2][1] + Finv[2][2]*Fn[0][2]*SigN[2][1]);
    Tsigma_1(6,2)=(Finv[0][2]*Fn[0][0]*SigN[0][2] + Finv[1][2]*Fn[0][1]*SigN[0][2] + Finv[2][2]*Fn[0][2]*SigN[0][2]);
    Tsigma_1(7,2)=(Finv[0][2]*Fn[0][0]*SigN[1][2] + Finv[1][2]*Fn[0][1]*SigN[1][2] + Finv[2][2]*Fn[0][2]*SigN[1][2]);
    Tsigma_1(8,2)=(Finv[0][2]*Fn[0][0]*SigN[2][2] + Finv[1][2]*Fn[0][1]*SigN[2][2] + Finv[2][2]*Fn[0][2]*SigN[2][2]);

    Tsigma_1(0,3)=(Finv[0][0]*Fn[1][0]*SigN[0][0] + Finv[1][0]*Fn[1][1]*SigN[0][0] + Finv[2][0]*Fn[1][2]*SigN[0][0]);
    Tsigma_1(1,3)=(Finv[0][0]*Fn[1][0]*SigN[1][0] + Finv[1][0]*Fn[1][1]*SigN[1][0] + Finv[2][0]*Fn[1][2]*SigN[1][0]);
    Tsigma_1(2,3)=(Finv[0][0]*Fn[1][0]*SigN[2][0] + Finv[1][0]*Fn[1][1]*SigN[2][0] + Finv[2][0]*Fn[1][2]*SigN[2][0]);
    Tsigma_1(3,3)=(Finv[0][0]*Fn[1][0]*SigN[0][1] + Finv[1][0]*Fn[1][1]*SigN[0][1] + Finv[2][0]*Fn[1][2]*SigN[0][1]);
    Tsigma_1(4,3)=(Finv[0][0]*Fn[1][0]*SigN[1][1] + Finv[1][0]*Fn[1][1]*SigN[1][1] + Finv[2][0]*Fn[1][2]*SigN[1][1]);
    Tsigma_1(5,3)=(Finv[0][0]*Fn[1][0]*SigN[2][1] + Finv[1][0]*Fn[1][1]*SigN[2][1] + Finv[2][0]*Fn[1][2]*SigN[2][1]);
    Tsigma_1(6,3)=(Finv[0][0]*Fn[1][0]*SigN[0][2] + Finv[1][0]*Fn[1][1]*SigN[0][2] + Finv[2][0]*Fn[1][2]*SigN[0][2]);
    Tsigma_1(7,3)=(Finv[0][0]*Fn[1][0]*SigN[1][2] + Finv[1][0]*Fn[1][1]*SigN[1][2] + Finv[2][0]*Fn[1][2]*SigN[1][2]);
    Tsigma_1(8,3)=(Finv[0][0]*Fn[1][0]*SigN[2][2] + Finv[1][0]*Fn[1][1]*SigN[2][2] + Finv[2][0]*Fn[1][2]*SigN[2][2]);

    Tsigma_1(0,4)=(Finv[0][1]*Fn[1][0]*SigN[0][0] + Finv[1][1]*Fn[1][1]*SigN[0][0] + Finv[2][1]*Fn[1][2]*SigN[0][0]);
    Tsigma_1(1,4)=(Finv[0][1]*Fn[1][0]*SigN[1][0] + Finv[1][1]*Fn[1][1]*SigN[1][0] + Finv[2][1]*Fn[1][2]*SigN[1][0]);
    Tsigma_1(2,4)=(Finv[0][1]*Fn[1][0]*SigN[2][0] + Finv[1][1]*Fn[1][1]*SigN[2][0] + Finv[2][1]*Fn[1][2]*SigN[2][0]);
    Tsigma_1(3,4)=(Finv[0][1]*Fn[1][0]*SigN[0][1] + Finv[1][1]*Fn[1][1]*SigN[0][1] + Finv[2][1]*Fn[1][2]*SigN[0][1]);
    Tsigma_1(4,4)=(Finv[0][1]*Fn[1][0]*SigN[1][1] + Finv[1][1]*Fn[1][1]*SigN[1][1] + Finv[2][1]*Fn[1][2]*SigN[1][1]);
    Tsigma_1(5,4)=(Finv[0][1]*Fn[1][0]*SigN[2][1] + Finv[1][1]*Fn[1][1]*SigN[2][1] + Finv[2][1]*Fn[1][2]*SigN[2][1]);
    Tsigma_1(6,4)=(Finv[0][1]*Fn[1][0]*SigN[0][2] + Finv[1][1]*Fn[1][1]*SigN[0][2] + Finv[2][1]*Fn[1][2]*SigN[0][2]);
    Tsigma_1(7,4)=(Finv[0][1]*Fn[1][0]*SigN[1][2] + Finv[1][1]*Fn[1][1]*SigN[1][2] + Finv[2][1]*Fn[1][2]*SigN[1][2]);
    Tsigma_1(8,4)=(Finv[0][1]*Fn[1][0]*SigN[2][2] + Finv[1][1]*Fn[1][1]*SigN[2][2] + Finv[2][1]*Fn[1][2]*SigN[2][2]);

    Tsigma_1(0,5)=(Finv[0][2]*Fn[1][0]*SigN[0][0] + Finv[1][2]*Fn[1][1]*SigN[0][0] + Finv[2][2]*Fn[1][2]*SigN[0][0]);
    Tsigma_1(1,5)=(Finv[0][2]*Fn[1][0]*SigN[1][0] + Finv[1][2]*Fn[1][1]*SigN[1][0] + Finv[2][2]*Fn[1][2]*SigN[1][0]);
    Tsigma_1(2,5)=(Finv[0][2]*Fn[1][0]*SigN[2][0] + Finv[1][2]*Fn[1][1]*SigN[2][0] + Finv[2][2]*Fn[1][2]*SigN[2][0]);
    Tsigma_1(3,5)=(Finv[0][2]*Fn[1][0]*SigN[0][1] + Finv[1][2]*Fn[1][1]*SigN[0][1] + Finv[2][2]*Fn[1][2]*SigN[0][1]);
    Tsigma_1(4,5)=(Finv[0][2]*Fn[1][0]*SigN[1][1] + Finv[1][2]*Fn[1][1]*SigN[1][1] + Finv[2][2]*Fn[1][2]*SigN[1][1]);
    Tsigma_1(5,5)=(Finv[0][2]*Fn[1][0]*SigN[2][1] + Finv[1][2]*Fn[1][1]*SigN[2][1] + Finv[2][2]*Fn[1][2]*SigN[2][1]);
    Tsigma_1(6,5)=(Finv[0][2]*Fn[1][0]*SigN[0][2] + Finv[1][2]*Fn[1][1]*SigN[0][2] + Finv[2][2]*Fn[1][2]*SigN[0][2]);
    Tsigma_1(7,5)=(Finv[0][2]*Fn[1][0]*SigN[1][2] + Finv[1][2]*Fn[1][1]*SigN[1][2] + Finv[2][2]*Fn[1][2]*SigN[1][2]);
    Tsigma_1(8,5)=(Finv[0][2]*Fn[1][0]*SigN[2][2] + Finv[1][2]*Fn[1][1]*SigN[2][2] + Finv[2][2]*Fn[1][2]*SigN[2][2]);

    Tsigma_1(0,6)=(Finv[0][0]*Fn[2][0]*SigN[0][0] + Finv[1][0]*Fn[2][1]*SigN[0][0] + Finv[2][0]*Fn[2][2]*SigN[0][0]);
    Tsigma_1(1,6)=(Finv[0][0]*Fn[2][0]*SigN[1][0] + Finv[1][0]*Fn[2][1]*SigN[1][0] + Finv[2][0]*Fn[2][2]*SigN[1][0]);
    Tsigma_1(2,6)=(Finv[0][0]*Fn[2][0]*SigN[2][0] + Finv[1][0]*Fn[2][1]*SigN[2][0] + Finv[2][0]*Fn[2][2]*SigN[2][0]);
    Tsigma_1(3,6)=(Finv[0][0]*Fn[2][0]*SigN[0][1] + Finv[1][0]*Fn[2][1]*SigN[0][1] + Finv[2][0]*Fn[2][2]*SigN[0][1]);
    Tsigma_1(4,6)=(Finv[0][0]*Fn[2][0]*SigN[1][1] + Finv[1][0]*Fn[2][1]*SigN[1][1] + Finv[2][0]*Fn[2][2]*SigN[1][1]);
    Tsigma_1(5,6)=(Finv[0][0]*Fn[2][0]*SigN[2][1] + Finv[1][0]*Fn[2][1]*SigN[2][1] + Finv[2][0]*Fn[2][2]*SigN[2][1]);
    Tsigma_1(6,6)=(Finv[0][0]*Fn[2][0]*SigN[0][2] + Finv[1][0]*Fn[2][1]*SigN[0][2] + Finv[2][0]*Fn[2][2]*SigN[0][2]);
    Tsigma_1(7,6)=(Finv[0][0]*Fn[2][0]*SigN[1][2] + Finv[1][0]*Fn[2][1]*SigN[1][2] + Finv[2][0]*Fn[2][2]*SigN[1][2]);
    Tsigma_1(8,6)=(Finv[0][0]*Fn[2][0]*SigN[2][2] + Finv[1][0]*Fn[2][1]*SigN[2][2] + Finv[2][0]*Fn[2][2]*SigN[2][2]);


    Tsigma_1(0,7)=(Finv[0][1]*Fn[2][0]*SigN[0][0] + Finv[1][1]*Fn[2][1]*SigN[0][0] + Finv[2][1]*Fn[2][2]*SigN[0][0]);
    Tsigma_1(1,7)=(Finv[0][1]*Fn[2][0]*SigN[1][0] + Finv[1][1]*Fn[2][1]*SigN[1][0] + Finv[2][1]*Fn[2][2]*SigN[1][0]);
    Tsigma_1(2,7)=(Finv[0][1]*Fn[2][0]*SigN[2][0] + Finv[1][1]*Fn[2][1]*SigN[2][0] + Finv[2][1]*Fn[2][2]*SigN[2][0]);
    Tsigma_1(3,7)=(Finv[0][1]*Fn[2][0]*SigN[0][1] + Finv[1][1]*Fn[2][1]*SigN[0][1] + Finv[2][1]*Fn[2][2]*SigN[0][1]);
    Tsigma_1(4,7)=(Finv[0][1]*Fn[2][0]*SigN[1][1] + Finv[1][1]*Fn[2][1]*SigN[1][1] + Finv[2][1]*Fn[2][2]*SigN[1][1]);
    Tsigma_1(5,7)=(Finv[0][1]*Fn[2][0]*SigN[2][1] + Finv[1][1]*Fn[2][1]*SigN[2][1] + Finv[2][1]*Fn[2][2]*SigN[2][1]);
    Tsigma_1(6,7)=(Finv[0][1]*Fn[2][0]*SigN[0][2] + Finv[1][1]*Fn[2][1]*SigN[0][2] + Finv[2][1]*Fn[2][2]*SigN[0][2]);
    Tsigma_1(7,7)=(Finv[0][1]*Fn[2][0]*SigN[1][2] + Finv[1][1]*Fn[2][1]*SigN[1][2] + Finv[2][1]*Fn[2][2]*SigN[1][2]);
    Tsigma_1(8,7)=(Finv[0][1]*Fn[2][0]*SigN[2][2] + Finv[1][1]*Fn[2][1]*SigN[2][2] + Finv[2][1]*Fn[2][2]*SigN[2][2]);

    Tsigma_1(0,8)=(Finv[0][2]*Fn[2][0]*SigN[0][0] + Finv[1][2]*Fn[2][1]*SigN[0][0] + Finv[2][2]*Fn[2][2]*SigN[0][0]);
    Tsigma_1(1,8)=(Finv[0][2]*Fn[2][0]*SigN[1][0] + Finv[1][2]*Fn[2][1]*SigN[1][0] + Finv[2][2]*Fn[2][2]*SigN[1][0]);
    Tsigma_1(2,8)=(Finv[0][2]*Fn[2][0]*SigN[2][0] + Finv[1][2]*Fn[2][1]*SigN[2][0] + Finv[2][2]*Fn[2][2]*SigN[2][0]);
    Tsigma_1(3,8)=(Finv[0][2]*Fn[2][0]*SigN[0][1] + Finv[1][2]*Fn[2][1]*SigN[0][1] + Finv[2][2]*Fn[2][2]*SigN[0][1]);
    Tsigma_1(4,8)=(Finv[0][2]*Fn[2][0]*SigN[1][1] + Finv[1][2]*Fn[2][1]*SigN[1][1] + Finv[2][2]*Fn[2][2]*SigN[1][1]);
    Tsigma_1(5,8)=(Finv[0][2]*Fn[2][0]*SigN[2][1] + Finv[1][2]*Fn[2][1]*SigN[2][1] + Finv[2][2]*Fn[2][2]*SigN[2][1]);
    Tsigma_1(6,8)=(Finv[0][2]*Fn[2][0]*SigN[0][2] + Finv[1][2]*Fn[2][1]*SigN[0][2] + Finv[2][2]*Fn[2][2]*SigN[0][2]);
    Tsigma_1(7,8)=(Finv[0][2]*Fn[2][0]*SigN[1][2] + Finv[1][2]*Fn[2][1]*SigN[1][2] + Finv[2][2]*Fn[2][2]*SigN[1][2]);
    Tsigma_1(8,8)=(Finv[0][2]*Fn[2][0]*SigN[2][2] + Finv[1][2]*Fn[2][1]*SigN[2][2] + Finv[2][2]*Fn[2][2]*SigN[2][2]);*/




}

void FSMicromorphic3DT::Form_Tsigma_2_matrix()
{

    int row=0;
    int col=0;
    Tsigma_2=0.0;
    for(int i=0;i<3;i++)
    {
        for(int m=0;m<3;m++)
        {
            row=0;
            for(int k=0;k<3;k++)
            {
                for(int l=0;l<3;l++)
                {
                    //summation over the same term
                    for(int L=0;L<3;L++)
                    {
                        Tsigma_2(row,col)+=SigN[i][k]*Fn[l][L]*Finv[L][m];
                    }
                    row++;
                }
            }
        col++;
        }
    }

/*    Tsigma_2=0.0;
    Tsigma_2(0,0)=(Finv[0][0]*Fn[0][0]*SigN[0][0] + Finv[1][0]*Fn[0][1]*SigN[0][0] + Finv[2][0]*Fn[0][2]*SigN[0][0]);
    Tsigma_2(1,0)=(Finv[0][0]*Fn[1][0]*SigN[0][0] + Finv[1][0]*Fn[1][1]*SigN[0][0] + Finv[2][0]*Fn[1][2]*SigN[0][0]);
    Tsigma_2(2,0)=(Finv[0][0]*Fn[2][0]*SigN[0][0] + Finv[1][0]*Fn[2][1]*SigN[0][0] + Finv[2][0]*Fn[2][2]*SigN[0][0]);
    Tsigma_2(3,0)=(Finv[0][0]*Fn[0][0]*SigN[0][1] + Finv[1][0]*Fn[0][1]*SigN[0][1] + Finv[2][0]*Fn[0][2]*SigN[0][1]);
    Tsigma_2(4,0)=(Finv[0][0]*Fn[1][0]*SigN[0][1] + Finv[1][0]*Fn[1][1]*SigN[0][1] + Finv[2][0]*Fn[1][2]*SigN[0][1]);
    Tsigma_2(5,0)=(Finv[0][0]*Fn[2][0]*SigN[0][1] + Finv[1][0]*Fn[2][1]*SigN[0][1] + Finv[2][0]*Fn[2][2]*SigN[0][1]);
    Tsigma_2(6,0)=(Finv[0][0]*Fn[0][0]*SigN[0][2] + Finv[1][0]*Fn[0][1]*SigN[0][2] + Finv[2][0]*Fn[0][2]*SigN[0][2]);
    Tsigma_2(7,0)=(Finv[0][0]*Fn[1][0]*SigN[0][2] + Finv[1][0]*Fn[1][1]*SigN[0][2] + Finv[2][0]*Fn[1][2]*SigN[0][2]);
    Tsigma_2(8,0)=(Finv[0][0]*Fn[2][0]*SigN[0][2] + Finv[1][0]*Fn[2][1]*SigN[0][2] + Finv[2][0]*Fn[2][2]*SigN[0][2]);

    Tsigma_2(0,1)=(Finv[0][1]*Fn[0][0]*SigN[0][0] + Finv[1][1]*Fn[0][1]*SigN[0][0] + Finv[2][1]*Fn[0][2]*SigN[0][0]);
    Tsigma_2(1,1)=(Finv[0][1]*Fn[1][0]*SigN[0][0] + Finv[1][1]*Fn[1][1]*SigN[0][0] + Finv[2][1]*Fn[1][2]*SigN[0][0]);
    Tsigma_2(2,1)=(Finv[0][1]*Fn[2][0]*SigN[0][0] + Finv[1][1]*Fn[2][1]*SigN[0][0] + Finv[2][1]*Fn[2][2]*SigN[0][0]);
    Tsigma_2(3,1)=(Finv[0][1]*Fn[0][0]*SigN[0][1] + Finv[1][1]*Fn[0][1]*SigN[0][1] + Finv[2][1]*Fn[0][2]*SigN[0][1]);
    Tsigma_2(4,1)=(Finv[0][1]*Fn[1][0]*SigN[0][1] + Finv[1][1]*Fn[1][1]*SigN[0][1] + Finv[2][1]*Fn[1][2]*SigN[0][1]);
    Tsigma_2(5,1)=(Finv[0][1]*Fn[2][0]*SigN[0][1] + Finv[1][1]*Fn[2][1]*SigN[0][1] + Finv[2][1]*Fn[2][2]*SigN[0][1]);
    Tsigma_2(6,1)=(Finv[0][1]*Fn[0][0]*SigN[0][2] + Finv[1][1]*Fn[0][1]*SigN[0][2] + Finv[2][1]*Fn[0][2]*SigN[0][2]);
    Tsigma_2(7,1)=(Finv[0][1]*Fn[1][0]*SigN[0][2] + Finv[1][1]*Fn[1][1]*SigN[0][2] + Finv[2][1]*Fn[1][2]*SigN[0][2]);
    Tsigma_2(8,1)=(Finv[0][1]*Fn[2][0]*SigN[0][2] + Finv[1][1]*Fn[2][1]*SigN[0][2] + Finv[2][1]*Fn[2][2]*SigN[0][2]);

    Tsigma_2(0,2)=(Finv[0][2]*Fn[0][0]*SigN[0][0] + Finv[1][2]*Fn[0][1]*SigN[0][0] + Finv[2][2]*Fn[0][2]*SigN[0][0]);
    Tsigma_2(1,2)=(Finv[0][2]*Fn[1][0]*SigN[0][0] + Finv[1][2]*Fn[1][1]*SigN[0][0] + Finv[2][2]*Fn[1][2]*SigN[0][0]);
    Tsigma_2(2,2)=(Finv[0][2]*Fn[2][0]*SigN[0][0] + Finv[1][2]*Fn[2][1]*SigN[0][0] + Finv[2][2]*Fn[2][2]*SigN[0][0]);
    Tsigma_2(3,2)=(Finv[0][2]*Fn[0][0]*SigN[0][1] + Finv[1][2]*Fn[0][1]*SigN[0][1] + Finv[2][2]*Fn[0][2]*SigN[0][1]);
    Tsigma_2(4,2)=(Finv[0][2]*Fn[1][0]*SigN[0][1] + Finv[1][2]*Fn[1][1]*SigN[0][1] + Finv[2][2]*Fn[1][2]*SigN[0][1]);
    Tsigma_2(5,2)=(Finv[0][2]*Fn[2][0]*SigN[0][1] + Finv[1][2]*Fn[2][1]*SigN[0][1] + Finv[2][2]*Fn[2][2]*SigN[0][1]);
    Tsigma_2(6,2)=(Finv[0][2]*Fn[0][0]*SigN[0][2] + Finv[1][2]*Fn[0][1]*SigN[0][2] + Finv[2][2]*Fn[0][2]*SigN[0][2]);
    Tsigma_2(7,2)=(Finv[0][2]*Fn[1][0]*SigN[0][2] + Finv[1][2]*Fn[1][1]*SigN[0][2] + Finv[2][2]*Fn[1][2]*SigN[0][2]);
    Tsigma_2(8,2)=(Finv[0][2]*Fn[2][0]*SigN[0][2] + Finv[1][2]*Fn[2][1]*SigN[0][2] + Finv[2][2]*Fn[2][2]*SigN[0][2]);

    Tsigma_2(0,3)=(Finv[0][0]*Fn[0][0]*SigN[1][0] + Finv[1][0]*Fn[0][1]*SigN[1][0] + Finv[2][0]*Fn[0][2]*SigN[1][0]);
    Tsigma_2(1,3)=(Finv[0][0]*Fn[1][0]*SigN[1][0] + Finv[1][0]*Fn[1][1]*SigN[1][0] + Finv[2][0]*Fn[1][2]*SigN[1][0]);
    Tsigma_2(2,3)=(Finv[0][0]*Fn[2][0]*SigN[1][0] + Finv[1][0]*Fn[2][1]*SigN[1][0] + Finv[2][0]*Fn[2][2]*SigN[1][0]);
    Tsigma_2(3,3)=(Finv[0][0]*Fn[0][0]*SigN[1][1] + Finv[1][0]*Fn[0][1]*SigN[1][1] + Finv[2][0]*Fn[0][2]*SigN[1][1]);
    Tsigma_2(4,3)=(Finv[0][0]*Fn[1][0]*SigN[1][1] + Finv[1][0]*Fn[1][1]*SigN[1][1] + Finv[2][0]*Fn[1][2]*SigN[1][1]);
    Tsigma_2(5,3)=(Finv[0][0]*Fn[2][0]*SigN[1][1] + Finv[1][0]*Fn[2][1]*SigN[1][1] + Finv[2][0]*Fn[2][2]*SigN[1][1]);
    Tsigma_2(6,3)=(Finv[0][0]*Fn[0][0]*SigN[1][2] + Finv[1][0]*Fn[0][1]*SigN[1][2] + Finv[2][0]*Fn[0][2]*SigN[1][2]);
    Tsigma_2(7,3)=(Finv[0][0]*Fn[1][0]*SigN[1][2] + Finv[1][0]*Fn[1][1]*SigN[1][2] + Finv[2][0]*Fn[1][2]*SigN[1][2]);
    Tsigma_2(8,3)=(Finv[0][0]*Fn[2][0]*SigN[1][2] + Finv[1][0]*Fn[2][1]*SigN[1][2] + Finv[2][0]*Fn[2][2]*SigN[1][2]);

    Tsigma_2(0,4)=(Finv[0][1]*Fn[0][0]*SigN[1][0] + Finv[1][1]*Fn[0][1]*SigN[1][0] + Finv[2][1]*Fn[0][2]*SigN[1][0]);
    Tsigma_2(1,4)=(Finv[0][1]*Fn[1][0]*SigN[1][0] + Finv[1][1]*Fn[1][1]*SigN[1][0] + Finv[2][1]*Fn[1][2]*SigN[1][0]);
    Tsigma_2(2,4)=(Finv[0][1]*Fn[2][0]*SigN[1][0] + Finv[1][1]*Fn[2][1]*SigN[1][0] + Finv[2][1]*Fn[2][2]*SigN[1][0]);
    Tsigma_2(3,4)=(Finv[0][1]*Fn[0][0]*SigN[1][1] + Finv[1][1]*Fn[0][1]*SigN[1][1] + Finv[2][1]*Fn[0][2]*SigN[1][1]);
    Tsigma_2(4,4)=(Finv[0][1]*Fn[1][0]*SigN[1][1] + Finv[1][1]*Fn[1][1]*SigN[1][1] + Finv[2][1]*Fn[1][2]*SigN[1][1]);
    Tsigma_2(5,4)=(Finv[0][1]*Fn[2][0]*SigN[1][1] + Finv[1][1]*Fn[2][1]*SigN[1][1] + Finv[2][1]*Fn[2][2]*SigN[1][1]);
    Tsigma_2(6,4)=(Finv[0][1]*Fn[0][0]*SigN[1][2] + Finv[1][1]*Fn[0][1]*SigN[1][2] + Finv[2][1]*Fn[0][2]*SigN[1][2]);
    Tsigma_2(7,4)=(Finv[0][1]*Fn[1][0]*SigN[1][2] + Finv[1][1]*Fn[1][1]*SigN[1][2] + Finv[2][1]*Fn[1][2]*SigN[1][2]);
    Tsigma_2(8,4)=(Finv[0][1]*Fn[2][0]*SigN[1][2] + Finv[1][1]*Fn[2][1]*SigN[1][2] + Finv[2][1]*Fn[2][2]*SigN[1][2]);

    Tsigma_2(0,5)=(Finv[0][2]*Fn[0][0]*SigN[1][0] + Finv[1][2]*Fn[0][1]*SigN[1][0] + Finv[2][2]*Fn[0][2]*SigN[1][0]);
    Tsigma_2(1,5)=(Finv[0][2]*Fn[1][0]*SigN[1][0] + Finv[1][2]*Fn[1][1]*SigN[1][0] + Finv[2][2]*Fn[1][2]*SigN[1][0]);
    Tsigma_2(2,5)=(Finv[0][2]*Fn[2][0]*SigN[1][0] + Finv[1][2]*Fn[2][1]*SigN[1][0] + Finv[2][2]*Fn[2][2]*SigN[1][0]);
    Tsigma_2(3,5)=(Finv[0][2]*Fn[0][0]*SigN[1][1] + Finv[1][2]*Fn[0][1]*SigN[1][1] + Finv[2][2]*Fn[0][2]*SigN[1][1]);
    Tsigma_2(4,5)=(Finv[0][2]*Fn[1][0]*SigN[1][1] + Finv[1][2]*Fn[1][1]*SigN[1][1] + Finv[2][2]*Fn[1][2]*SigN[1][1]);
    Tsigma_2(5,5)=(Finv[0][2]*Fn[2][0]*SigN[1][1] + Finv[1][2]*Fn[2][1]*SigN[1][1] + Finv[2][2]*Fn[2][2]*SigN[1][1]);
    Tsigma_2(6,5)=(Finv[0][2]*Fn[0][0]*SigN[1][2] + Finv[1][2]*Fn[0][1]*SigN[1][2] + Finv[2][2]*Fn[0][2]*SigN[1][2]);
    Tsigma_2(7,5)=(Finv[0][2]*Fn[1][0]*SigN[1][2] + Finv[1][2]*Fn[1][1]*SigN[1][2] + Finv[2][2]*Fn[1][2]*SigN[1][2]);
    Tsigma_2(8,5)=(Finv[0][2]*Fn[2][0]*SigN[1][2] + Finv[1][2]*Fn[2][1]*SigN[1][2] + Finv[2][2]*Fn[2][2]*SigN[1][2]);

    Tsigma_2(0,6)=(Finv[0][0]*Fn[0][0]*SigN[2][0] + Finv[1][0]*Fn[0][1]*SigN[2][0] + Finv[2][0]*Fn[0][2]*SigN[2][0]);
    Tsigma_2(1,6)=(Finv[0][0]*Fn[1][0]*SigN[2][0] + Finv[1][0]*Fn[1][1]*SigN[2][0] + Finv[2][0]*Fn[1][2]*SigN[2][0]);
    Tsigma_2(2,6)=(Finv[0][0]*Fn[2][0]*SigN[2][0] + Finv[1][0]*Fn[2][1]*SigN[2][0] + Finv[2][0]*Fn[2][2]*SigN[2][0]);
    Tsigma_2(3,6)=(Finv[0][0]*Fn[0][0]*SigN[2][1] + Finv[1][0]*Fn[0][1]*SigN[2][1] + Finv[2][0]*Fn[0][2]*SigN[2][1]);
    Tsigma_2(4,6)=(Finv[0][0]*Fn[1][0]*SigN[2][1] + Finv[1][0]*Fn[1][1]*SigN[2][1] + Finv[2][0]*Fn[1][2]*SigN[2][1]);
    Tsigma_2(5,6)=(Finv[0][0]*Fn[2][0]*SigN[2][1] + Finv[1][0]*Fn[2][1]*SigN[2][1] + Finv[2][0]*Fn[2][2]*SigN[2][1]);
    Tsigma_2(6,6)=(Finv[0][0]*Fn[0][0]*SigN[2][2] + Finv[1][0]*Fn[0][1]*SigN[2][2] + Finv[2][0]*Fn[0][2]*SigN[2][2]);
    Tsigma_2(7,6)=(Finv[0][0]*Fn[1][0]*SigN[2][2] + Finv[1][0]*Fn[1][1]*SigN[2][2] + Finv[2][0]*Fn[1][2]*SigN[2][2]);
    Tsigma_2(8,6)=(Finv[0][0]*Fn[2][0]*SigN[2][2] + Finv[1][0]*Fn[2][1]*SigN[2][2] + Finv[2][0]*Fn[2][2]*SigN[2][2]);

    Tsigma_2(0,7)=(Finv[0][1]*Fn[0][0]*SigN[2][0] + Finv[1][1]*Fn[0][1]*SigN[2][0] + Finv[2][1]*Fn[0][2]*SigN[2][0]);
    Tsigma_2(1,7)=(Finv[0][1]*Fn[1][0]*SigN[2][0] + Finv[1][1]*Fn[1][1]*SigN[2][0] + Finv[2][1]*Fn[1][2]*SigN[2][0]);
    Tsigma_2(2,7)=(Finv[0][1]*Fn[2][0]*SigN[2][0] + Finv[1][1]*Fn[2][1]*SigN[2][0] + Finv[2][1]*Fn[2][2]*SigN[2][0]);
    Tsigma_2(3,7)=(Finv[0][1]*Fn[0][0]*SigN[2][1] + Finv[1][1]*Fn[0][1]*SigN[2][1] + Finv[2][1]*Fn[0][2]*SigN[2][1]);
    Tsigma_2(4,7)=(Finv[0][1]*Fn[1][0]*SigN[2][1] + Finv[1][1]*Fn[1][1]*SigN[2][1] + Finv[2][1]*Fn[1][2]*SigN[2][1]);
    Tsigma_2(5,7)=(Finv[0][1]*Fn[2][0]*SigN[2][1] + Finv[1][1]*Fn[2][1]*SigN[2][1] + Finv[2][1]*Fn[2][2]*SigN[2][1]);
    Tsigma_2(6,7)=(Finv[0][1]*Fn[0][0]*SigN[2][2] + Finv[1][1]*Fn[0][1]*SigN[2][2] + Finv[2][1]*Fn[0][2]*SigN[2][2]);
    Tsigma_2(7,7)=(Finv[0][1]*Fn[1][0]*SigN[2][2] + Finv[1][1]*Fn[1][1]*SigN[2][2] + Finv[2][1]*Fn[1][2]*SigN[2][2]);
    Tsigma_2(8,7)=(Finv[0][1]*Fn[2][0]*SigN[2][2] + Finv[1][1]*Fn[2][1]*SigN[2][2] + Finv[2][1]*Fn[2][2]*SigN[2][2]);

    Tsigma_2(0,8)=(Finv[0][2]*Fn[0][0]*SigN[2][0] + Finv[1][2]*Fn[0][1]*SigN[2][0] + Finv[2][2]*Fn[0][2]*SigN[2][0]);
    Tsigma_2(1,8)=(Finv[0][2]*Fn[1][0]*SigN[2][0] + Finv[1][2]*Fn[1][1]*SigN[2][0] + Finv[2][2]*Fn[1][2]*SigN[2][0]);
    Tsigma_2(2,8)=(Finv[0][2]*Fn[2][0]*SigN[2][0] + Finv[1][2]*Fn[2][1]*SigN[2][0] + Finv[2][2]*Fn[2][2]*SigN[2][0]);
    Tsigma_2(3,8)=(Finv[0][2]*Fn[0][0]*SigN[2][1] + Finv[1][2]*Fn[0][1]*SigN[2][1] + Finv[2][2]*Fn[0][2]*SigN[2][1]);
    Tsigma_2(4,8)=(Finv[0][2]*Fn[1][0]*SigN[2][1] + Finv[1][2]*Fn[1][1]*SigN[2][1] + Finv[2][2]*Fn[1][2]*SigN[2][1]);
    Tsigma_2(5,8)=(Finv[0][2]*Fn[2][0]*SigN[2][1] + Finv[1][2]*Fn[2][1]*SigN[2][1] + Finv[2][2]*Fn[2][2]*SigN[2][1]);
    Tsigma_2(6,8)=(Finv[0][2]*Fn[0][0]*SigN[2][2] + Finv[1][2]*Fn[0][1]*SigN[2][2] + Finv[2][2]*Fn[0][2]*SigN[2][2]);
    Tsigma_2(7,8)=(Finv[0][2]*Fn[1][0]*SigN[2][2] + Finv[1][2]*Fn[1][1]*SigN[2][2] + Finv[2][2]*Fn[1][2]*SigN[2][2]);
    Tsigma_2(8,8)=(Finv[0][2]*Fn[2][0]*SigN[2][2] + Finv[1][2]*Fn[2][1]*SigN[2][2] + Finv[2][2]*Fn[2][2]*SigN[2][2]);
*/

}

void FSMicromorphic3DT::Form_Tsigma_3_matrix()
{
    int row=0;
    int col=0;
    Tsigma_3=0.0;
    for(int i=0;i<3;i++)
    {
        for(int m=0;m<3;m++)
        {
            row=0;
            for(int k=0;k<3;k++)
            {
                for(int l=0;l<3;l++)
                {
                    //summation on the same term
                    for(int L=0;L<3;L++)
                    {
                        Tsigma_3(row,col)+=SigN[l][i]*Fn[k][L]*Finv[L][m];
                    }
                    row++;
                }
            }
        col++;
        }
    }
/*    int row=0;
    int col=0;
    Tsigma_3=0.0;
    for(int i=0;i<3;i++)
    {
        for(int m=0;m<3;m++)
        {
            row=0;
            for(int l=0;l<3;l++)
            {
                for(int k=0;k<3;k++)
                {
                    //summation on the same term
                    for(int L=0;L<3;L++)
                    {
                        Tsigma_3(row,col)=Tsigma_3(row,col)+SigN[l][i]*Fn[k][L]*Finv[L][m];
                    }
                    row++;
                }
            }
        col++;
        }
    }*/
}

void FSMicromorphic3DT::Form_TFn_1_matrix()
    {

       int row=0;
        int col=0;
        TFn_1=0.0;
        for(int i=0;i<3;i++)
        {
            for(int m=0;m<3;m++)
            {
                row=0;
                for(int k=0;k<3;k++)
                {
                        for(int l=0;l<3;l++)
                    {
                        //summation on the same term
                        for(int L=0;L<3;L++)
                        {
                            TFn_1(row,col)+=fIdentity_matrix(l,k)*Fn[i][L]*Finv[L][m];
                        }
                        row++;
                    }
                }
            col++;
            }
        }

/*

    int row=0;
    int col=0;
    TFn_1=0.0;
    for(int i=0;i<3;i++)
    {
        for(int m=0;m<3;m++)
        {
            row=0;
            for(int l=0;l<3;l++)
            {
                for(int k=0;k<3;k++)
                {
                    //summation on the same term
                    for(int L=0;L<3;L++)
                    {
                        TFn_1(row,col)=TFn_1(row,col)+KrDelta[l][k]*Fn[i][L]*Finv[L][m];
                    }
                    row++;
                }
            }
        col++;
        }
    }
*/

 }

void FSMicromorphic3DT::Form_TFn_2_matrix()
{
    int row=0;
    int col=0;
    TFn_2=0.0;
    for(int k=0;k<3;k++)
    {
        for(int m=0;m<3;m++)
        {
            //row operations
            row=3*k;
                for(int l=0;l<3;l++)
                {
                    //summation on the same term
                    for(int L=0;L<3;L++)
                    {
                        TFn_2(row,col)+=Fn[l][L]*Finv[L][m];
                    }
                    row++;
                }
            col++;
        }

    }

/*     int row=0;
        int col=0;
        TFn_2=0.0;
        for(int k=0;k<3;k++)
        {
            for(int m=0;m<3;m++)
            {
                //row operations
                row=k;
                    for(int l=0;l<3;l++)
                    {
                        //summation on the same term
                        for(int L=0;L<3;L++)
                        {
                            TFn_2(row,col)=TFn_2(row,col)+Fn[l][L]*Finv[L][m];
                        }
                        row+=3;
                    }

                col++;
            }
        }
*/

}

void FSMicromorphic3DT::Form_TFn_3_matrix()
    {

    int row=0;
    int col=0;
    TFn_3=0.0;
    for(int l=0;l<3;l++)
    {
        for(int m=0;m<3;m++)
        {
            row=l;
                for(int k=0;k<3;k++)
                {
                    //summation on the same term
                    for(int L=0;L<3;L++)
                    {
                        TFn_3(row,col)+=Fn[k][L]*Finv[L][m];
                    }
                    row=row+3;
                }
        col++;
        }
    }

    // mutliply with (Mu+sigma)
/*    TFn_3=0.0;
    TFn_3(0,0)=(Finv[0][0]*Fn[0][0] + Finv[1][0]*Fn[0][1] + Finv[2][0]*Fn[0][2]);//*w[0][0]
    TFn_3(3,0)=(Finv[0][0]*Fn[1][0] + Finv[1][0]*Fn[1][1] + Finv[2][0]*Fn[1][2]);//*w[1][0] +
    TFn_3(6,0)=(Finv[0][0]*Fn[2][0] + Finv[1][0]*Fn[2][1] + Finv[2][0]*Fn[2][2]);//*w[2][0]) +

    TFn_3(0,1)=(Finv[0][1]*Fn[0][0] + Finv[1][1]*Fn[0][1] + Finv[2][1]*Fn[0][2]);//*w[0][0] +
    TFn_3(3,1)=(Finv[0][1]*Fn[1][0] + Finv[1][1]*Fn[1][1] + Finv[2][1]*Fn[1][2]);//*w[1][0] +
    TFn_3(6,1)=(Finv[0][1]*Fn[2][0] + Finv[1][1]*Fn[2][1] + Finv[2][1]*Fn[2][2]);//*w[2][0]) +

    TFn_3(0,2)=(Finv[0][2]*Fn[0][0] + Finv[1][2]*Fn[0][1] + Finv[2][2]*Fn[0][2]);//*w[0][0] +
    TFn_3(3,2)=(Finv[0][2]*Fn[1][0] + Finv[1][2]*Fn[1][1] + Finv[2][2]*Fn[1][2]);//*w[1][0] +
    TFn_3(6,2)=(Finv[0][2]*Fn[2][0] + Finv[1][2]*Fn[2][1] + Finv[2][2]*Fn[2][2]);//*w[2][0]) +

    TFn_3(1,3)=(Finv[0][0]*Fn[0][0] + Finv[1][0]*Fn[0][1] + Finv[2][0]*Fn[0][2]);//*w[0][1] +
    TFn_3(4,3)=(Finv[0][0]*Fn[1][0] + Finv[1][0]*Fn[1][1] + Finv[2][0]*Fn[1][2]);//*w[1][1] +
    TFn_3(7,3)=(Finv[0][0]*Fn[2][0] + Finv[1][0]*Fn[2][1] + Finv[2][0]*Fn[2][2]);//*w[2][1]) +

    TFn_3(1,4)=(Finv[0][1]*Fn[0][0] + Finv[1][1]*Fn[0][1] + Finv[2][1]*Fn[0][2]);//*w[0][1] +
    TFn_3(4,4)=(Finv[0][1]*Fn[1][0] + Finv[1][1]*Fn[1][1] + Finv[2][1]*Fn[1][2]);//*w[1][1] +
    TFn_3(7,4)=(Finv[0][1]*Fn[2][0] + Finv[1][1]*Fn[2][1] + Finv[2][1]*Fn[2][2]);//*w[2][1]) +

    TFn_3(1,5)=(Finv[0][2]*Fn[0][0] + Finv[1][2]*Fn[0][1] + Finv[2][2]*Fn[0][2]);//*w[0][1] +
    TFn_3(4,5)=(Finv[0][2]*Fn[1][0] + Finv[1][2]*Fn[1][1] + Finv[2][2]*Fn[1][2]);//*w[1][1] +
    TFn_3(7,5)=(Finv[0][2]*Fn[2][0] + Finv[1][2]*Fn[2][1] + Finv[2][2]*Fn[2][2]);//*w[2][1]) +

    TFn_3(2,6)=(Finv[0][0]*Fn[0][0] + Finv[1][0]*Fn[0][1] + Finv[2][0]*Fn[0][2]);//*w[0][2] +
    TFn_3(5,6)=(Finv[0][0]*Fn[1][0] + Finv[1][0]*Fn[1][1] + Finv[2][0]*Fn[1][2]);//*w[1][2] +
    TFn_3(8,6)=(Finv[0][0]*Fn[2][0] + Finv[1][0]*Fn[2][1] + Finv[2][0]*Fn[2][2]);//*w[2][2]) +

    TFn_3(2,7)=(Finv[0][1]*Fn[0][0] + Finv[1][1]*Fn[0][1] + Finv[2][1]*Fn[0][2]);//*w[0][2] +
    TFn_3(5,7)=(Finv[0][1]*Fn[1][0] + Finv[1][1]*Fn[1][1] + Finv[2][1]*Fn[1][2]);//*w[1][2] +
    TFn_3(8,7)=(Finv[0][1]*Fn[2][0] + Finv[1][1]*Fn[2][1] + Finv[2][1]*Fn[2][2]);//*w[2][2]) +

    TFn_3(2,8)=(Finv[0][2]*Fn[0][0] + Finv[1][2]*Fn[0][1] + Finv[2][2]*Fn[0][2]);//*w[0][2] +
    TFn_3(5,8)=(Finv[0][2]*Fn[1][0] + Finv[1][2]*Fn[1][1] + Finv[2][2]*Fn[1][2]);//*w[1][2] +
    TFn_3(8,8)=(Finv[0][2]*Fn[2][0] + Finv[1][2]*Fn[2][1] + Finv[2][2]*Fn[2][2]);//*w[2][2])
*/
    }

void FSMicromorphic3DT::Form_TChi_1_matrix()
{
    int row=0;
    int col=0;
    TChi_1=0.0;
    for(int T=0;T<3;T++)
    {
        for(int m=0;m<3;m++)
        {
            row=0;
            for(int k=0;k<3;k++)
            {
                for(int l=0;l<3;l++)
                {
                    //summation over the same term starts here
                    for(int i=0;i<3;i++)
                    {
                        for(int L=0;L<3;L++)
                        {
                            TChi_1(row,col)+=ChiN[i][L]*ChiInv[L][m]*ChiInv[T][i]*fIdentity_matrix(l,k);
                        }
                    }
                    row++;
                }
            }
            col++;
        }
    }


/*TChi_1(0,0)=(ChiInv[0][0]*ChiInv[0][0]*ChiN[0][0]*KrDelta[0][0] + ChiInv[0][0]*ChiInv[1][0]*ChiN[0][1]*KrDelta[0][0] +
             ChiInv[0][0]*ChiInv[2][0]*ChiN[0][2]*KrDelta[0][0] + ChiInv[0][0]*ChiInv[0][1]*ChiN[1][0]*KrDelta[0][0] +
             ChiInv[0][1]*ChiInv[1][0]*ChiN[1][1]*KrDelta[0][0] + ChiInv[0][1]*ChiInv[2][0]*ChiN[1][2]*KrDelta[0][0] +
             ChiInv[0][0]*ChiInv[0][2]*ChiN[2][0]*KrDelta[0][0] + ChiInv[0][2]*ChiInv[1][0]*ChiN[2][1]*KrDelta[0][0] +
             ChiInv[0][2]*ChiInv[2][0]*ChiN[2][2]*KrDelta[0][0]);//        *w[0][0] +

TChi_1(1,0)=(ChiInv[0][0]*ChiInv[0][0]*ChiN[0][0]*KrDelta[1][0] + ChiInv[0][0]*ChiInv[1][0]*ChiN[0][1]*KrDelta[1][0] +
             ChiInv[0][0]*ChiInv[2][0]*ChiN[0][2]*KrDelta[1][0] + ChiInv[0][0]*ChiInv[0][1]*ChiN[1][0]*KrDelta[1][0] +
             ChiInv[0][1]*ChiInv[1][0]*ChiN[1][1]*KrDelta[1][0] + ChiInv[0][1]*ChiInv[2][0]*ChiN[1][2]*KrDelta[1][0] +
             ChiInv[0][0]*ChiInv[0][2]*ChiN[2][0]*KrDelta[1][0] + ChiInv[0][2]*ChiInv[1][0]*ChiN[2][1]*KrDelta[1][0] +
             ChiInv[0][2]*ChiInv[2][0]*ChiN[2][2]*KrDelta[1][0]);//*w[0][1] +

TChi_1(2,0)=(ChiInv[0][0]*ChiInv[0][0]*ChiN[0][0]*KrDelta[2][0] + ChiInv[0][0]*ChiInv[1][0]*ChiN[0][1]*KrDelta[2][0] +
             ChiInv[0][0]*ChiInv[2][0]*ChiN[0][2]*KrDelta[2][0] + ChiInv[0][0]*ChiInv[0][1]*ChiN[1][0]*KrDelta[2][0] +
             ChiInv[0][1]*ChiInv[1][0]*ChiN[1][1]*KrDelta[2][0] + ChiInv[0][1]*ChiInv[2][0]*ChiN[1][2]*KrDelta[2][0] +
             ChiInv[0][0]*ChiInv[0][2]*ChiN[2][0]*KrDelta[2][0] + ChiInv[0][2]*ChiInv[1][0]*ChiN[2][1]*KrDelta[2][0] +
             ChiInv[0][2]*ChiInv[2][0]*ChiN[2][2]*KrDelta[2][0]);//*w[0][2] +
TChi_1(3,0)=(ChiInv[0][0]*ChiInv[0][0]*ChiN[0][0]*KrDelta[0][1] + ChiInv[0][0]*ChiInv[1][0]*ChiN[0][1]*KrDelta[0][1] +
             ChiInv[0][0]*ChiInv[2][0]*ChiN[0][2]*KrDelta[0][1] + ChiInv[0][0]*ChiInv[0][1]*ChiN[1][0]*KrDelta[0][1] +
             ChiInv[0][1]*ChiInv[1][0]*ChiN[1][1]*KrDelta[0][1] + ChiInv[0][1]*ChiInv[2][0]*ChiN[1][2]*KrDelta[0][1] +
             ChiInv[0][0]*ChiInv[0][2]*ChiN[2][0]*KrDelta[0][1] + ChiInv[0][2]*ChiInv[1][0]*ChiN[2][1]*KrDelta[0][1] +
             ChiInv[0][2]*ChiInv[2][0]*ChiN[2][2]*KrDelta[0][1]);//*w[1][0] +
TChi_1(4,0)=(ChiInv[0][0]*ChiInv[0][0]*ChiN[0][0]*KrDelta[1][1] + ChiInv[0][0]*ChiInv[1][0]*ChiN[0][1]*KrDelta[1][1] +
             ChiInv[0][0]*ChiInv[2][0]*ChiN[0][2]*KrDelta[1][1] + ChiInv[0][0]*ChiInv[0][1]*ChiN[1][0]*KrDelta[1][1] +
             ChiInv[0][1]*ChiInv[1][0]*ChiN[1][1]*KrDelta[1][1] + ChiInv[0][1]*ChiInv[2][0]*ChiN[1][2]*KrDelta[1][1] +
             ChiInv[0][0]*ChiInv[0][2]*ChiN[2][0]*KrDelta[1][1] + ChiInv[0][2]*ChiInv[1][0]*ChiN[2][1]*KrDelta[1][1] +
             ChiInv[0][2]*ChiInv[2][0]*ChiN[2][2]*KrDelta[1][1]);//*w[1][1] +
TChi_1(5,0)=(ChiInv[0][0]*ChiInv[0][0]*ChiN[0][0]*KrDelta[2][1] + ChiInv[0][0]*ChiInv[1][0]*ChiN[0][1]*KrDelta[2][1] +
             ChiInv[0][0]*ChiInv[2][0]*ChiN[0][2]*KrDelta[2][1] + ChiInv[0][0]*ChiInv[0][1]*ChiN[1][0]*KrDelta[2][1] +
             ChiInv[0][1]*ChiInv[1][0]*ChiN[1][1]*KrDelta[2][1] + ChiInv[0][1]*ChiInv[2][0]*ChiN[1][2]*KrDelta[2][1] +
             ChiInv[0][0]*ChiInv[0][2]*ChiN[2][0]*KrDelta[2][1] + ChiInv[0][2]*ChiInv[1][0]*ChiN[2][1]*KrDelta[2][1] +
             ChiInv[0][2]*ChiInv[2][0]*ChiN[2][2]*KrDelta[2][1]);//*w[1][2] +
TChi_1(6,0)=(ChiInv[0][0]*ChiInv[0][0]*ChiN[0][0]*KrDelta[0][2] + ChiInv[0][0]*ChiInv[1][0]*ChiN[0][1]*KrDelta[0][2] +
             ChiInv[0][0]*ChiInv[2][0]*ChiN[0][2]*KrDelta[0][2] + ChiInv[0][0]*ChiInv[0][1]*ChiN[1][0]*KrDelta[0][2] +
             ChiInv[0][1]*ChiInv[1][0]*ChiN[1][1]*KrDelta[0][2] + ChiInv[0][1]*ChiInv[2][0]*ChiN[1][2]*KrDelta[0][2] +
             ChiInv[0][0]*ChiInv[0][2]*ChiN[2][0]*KrDelta[0][2] + ChiInv[0][2]*ChiInv[1][0]*ChiN[2][1]*KrDelta[0][2] +
             ChiInv[0][2]*ChiInv[2][0]*ChiN[2][2]*KrDelta[0][2]);//*w[2][0] +
TChi_1(7,0)=(ChiInv[0][0]*ChiInv[0][0]*ChiN[0][0]*KrDelta[1][2] + ChiInv[0][0]*ChiInv[1][0]*ChiN[0][1]*KrDelta[1][2] +
             ChiInv[0][0]*ChiInv[2][0]*ChiN[0][2]*KrDelta[1][2] + ChiInv[0][0]*ChiInv[0][1]*ChiN[1][0]*KrDelta[1][2] +
             ChiInv[0][1]*ChiInv[1][0]*ChiN[1][1]*KrDelta[1][2] + ChiInv[0][1]*ChiInv[2][0]*ChiN[1][2]*KrDelta[1][2] +
             ChiInv[0][0]*ChiInv[0][2]*ChiN[2][0]*KrDelta[1][2] + ChiInv[0][2]*ChiInv[1][0]*ChiN[2][1]*KrDelta[1][2] +
             ChiInv[0][2]*ChiInv[2][0]*ChiN[2][2]*KrDelta[1][2]);//*w[2][1] +
TChi_1(8,0)=(ChiInv[0][0]*ChiInv[0][0]*ChiN[0][0]*KrDelta[2][2] + ChiInv[0][0]*ChiInv[1][0]*ChiN[0][1]*KrDelta[2][2] +
             ChiInv[0][0]*ChiInv[2][0]*ChiN[0][2]*KrDelta[2][2] + ChiInv[0][0]*ChiInv[0][1]*ChiN[1][0]*KrDelta[2][2] +
             ChiInv[0][1]*ChiInv[1][0]*ChiN[1][1]*KrDelta[2][2] + ChiInv[0][1]*ChiInv[2][0]*ChiN[1][2]*KrDelta[2][2] +
             ChiInv[0][0]*ChiInv[0][2]*ChiN[2][0]*KrDelta[2][2] + ChiInv[0][2]*ChiInv[1][0]*ChiN[2][1]*KrDelta[2][2] +
             ChiInv[0][2]*ChiInv[2][0]*ChiN[2][2]*KrDelta[2][2]);//*w[2][2]

TChi_1(0,1)=(ChiInv[0][0]*ChiInv[0][1]*ChiN[0][0]*KrDelta[0][0] + ChiInv[0][0]*ChiInv[1][1]*ChiN[0][1]*KrDelta[0][0] +
            ChiInv[0][0]*ChiInv[2][1]*ChiN[0][2]*KrDelta[0][0] + ChiInv[0][1]*ChiInv[0][1]*ChiN[1][0]*KrDelta[0][0] +
            ChiInv[0][1]*ChiInv[1][1]*ChiN[1][1]*KrDelta[0][0] + ChiInv[0][1]*ChiInv[2][1]*ChiN[1][2]*KrDelta[0][0] +
            ChiInv[0][1]*ChiInv[0][2]*ChiN[2][0]*KrDelta[0][0] + ChiInv[0][2]*ChiInv[1][1]*ChiN[2][1]*KrDelta[0][0] +
            ChiInv[0][2]*ChiInv[2][1]*ChiN[2][2]*KrDelta[0][0]);//*w[0][0] +


TChi_1(1,1)=ChiInv[0][0]*ChiInv[0][1]*ChiN[0][0]*KrDelta[1][0] + ChiInv[0][0]*ChiInv[1][1]*ChiN[0][1]*KrDelta[1][0] +
         ChiInv[0][0]*ChiInv[2][1]*ChiN[0][2]*KrDelta[1][0] + ChiInv[0][1]*ChiInv[0][1]*ChiN[1][0]*KrDelta[1][0] +
         ChiInv[0][1]*ChiInv[1][1]*ChiN[1][1]*KrDelta[1][0] + ChiInv[0][1]*ChiInv[2][1]*ChiN[1][2]*KrDelta[1][0] +
         ChiInv[0][1]*ChiInv[0][2]*ChiN[2][0]*KrDelta[1][0] + ChiInv[0][2]*ChiInv[1][1]*ChiN[2][1]*KrDelta[1][0] +
         ChiInv[0][2]*ChiInv[2][1]*ChiN[2][2]*KrDelta[1][0];//*w[0][1] +
TChi_1(2,1)=(ChiInv[0][0]*ChiInv[0][1]*ChiN[0][0]*KrDelta[2][0] + ChiInv[0][0]*ChiInv[1][1]*ChiN[0][1]*KrDelta[2][0] +
         ChiInv[0][0]*ChiInv[2][1]*ChiN[0][2]*KrDelta[2][0] + ChiInv[0][1]*ChiInv[0][1]*ChiN[1][0]*KrDelta[2][0] +
         ChiInv[0][1]*ChiInv[1][1]*ChiN[1][1]*KrDelta[2][0] + ChiInv[0][1]*ChiInv[2][1]*ChiN[1][2]*KrDelta[2][0] +
         ChiInv[0][1]*ChiInv[0][2]*ChiN[2][0]*KrDelta[2][0] + ChiInv[0][2]*ChiInv[1][1]*ChiN[2][1]*KrDelta[2][0] +
         ChiInv[0][2]*ChiInv[2][1]*ChiN[2][2]*KrDelta[2][0]);//*w[0][2] +
TChi_1(3,1)=(ChiInv[0][0]*ChiInv[0][1]*ChiN[0][0]*KrDelta[0][1] + ChiInv[0][0]*ChiInv[1][1]*ChiN[0][1]*KrDelta[0][1] +
         ChiInv[0][0]*ChiInv[2][1]*ChiN[0][2]*KrDelta[0][1] + ChiInv[0][1]*ChiInv[0][1]*ChiN[1][0]*KrDelta[0][1] +
         ChiInv[0][1]*ChiInv[1][1]*ChiN[1][1]*KrDelta[0][1] + ChiInv[0][1]*ChiInv[2][1]*ChiN[1][2]*KrDelta[0][1] +
         ChiInv[0][1]*ChiInv[0][2]*ChiN[2][0]*KrDelta[0][1] + ChiInv[0][2]*ChiInv[1][1]*ChiN[2][1]*KrDelta[0][1] +
         ChiInv[0][2]*ChiInv[2][1]*ChiN[2][2]*KrDelta[0][1]);//*w[1][0] +
TChi_1(4,1)=(ChiInv[0][0]*ChiInv[0][1]*ChiN[0][0]*KrDelta[1][1] + ChiInv[0][0]*ChiInv[1][1]*ChiN[0][1]*KrDelta[1][1] +
         ChiInv[0][0]*ChiInv[2][1]*ChiN[0][2]*KrDelta[1][1] + ChiInv[0][1]*ChiInv[0][1]*ChiN[1][0]*KrDelta[1][1] +
         ChiInv[0][1]*ChiInv[1][1]*ChiN[1][1]*KrDelta[1][1] + ChiInv[0][1]*ChiInv[2][1]*ChiN[1][2]*KrDelta[1][1] +
         ChiInv[0][1]*ChiInv[0][2]*ChiN[2][0]*KrDelta[1][1] + ChiInv[0][2]*ChiInv[1][1]*ChiN[2][1]*KrDelta[1][1] +
         ChiInv[0][2]*ChiInv[2][1]*ChiN[2][2]*KrDelta[1][1]);//*w[1][1] +
TChi_1(5,1)=(ChiInv[0][0]*ChiInv[0][1]*ChiN[0][0]*KrDelta[2][1] + ChiInv[0][0]*ChiInv[1][1]*ChiN[0][1]*KrDelta[2][1] +
         ChiInv[0][0]*ChiInv[2][1]*ChiN[0][2]*KrDelta[2][1] + ChiInv[0][1]*ChiInv[0][1]*ChiN[1][0]*KrDelta[2][1] +
         ChiInv[0][1]*ChiInv[1][1]*ChiN[1][1]*KrDelta[2][1] + ChiInv[0][1]*ChiInv[2][1]*ChiN[1][2]*KrDelta[2][1] +
         ChiInv[0][1]*ChiInv[0][2]*ChiN[2][0]*KrDelta[2][1] + ChiInv[0][2]*ChiInv[1][1]*ChiN[2][1]*KrDelta[2][1] +
         ChiInv[0][2]*ChiInv[2][1]*ChiN[2][2]*KrDelta[2][1]);//*w[1][2] +
TChi_1(6,1)=(ChiInv[0][0]*ChiInv[0][1]*ChiN[0][0]*KrDelta[0][2] + ChiInv[0][0]*ChiInv[1][1]*ChiN[0][1]*KrDelta[0][2] +
         ChiInv[0][0]*ChiInv[2][1]*ChiN[0][2]*KrDelta[0][2] + ChiInv[0][1]*ChiInv[0][1]*ChiN[1][0]*KrDelta[0][2] +
         ChiInv[0][1]*ChiInv[1][1]*ChiN[1][1]*KrDelta[0][2] + ChiInv[0][1]*ChiInv[2][1]*ChiN[1][2]*KrDelta[0][2] +
         ChiInv[0][1]*ChiInv[0][2]*ChiN[2][0]*KrDelta[0][2] + ChiInv[0][2]*ChiInv[1][1]*ChiN[2][1]*KrDelta[0][2] +
         ChiInv[0][2]*ChiInv[2][1]*ChiN[2][2]*KrDelta[0][2]);//*w[2][0] +
TChi_1(7,1)=(ChiInv[0][0]*ChiInv[0][1]*ChiN[0][0]*KrDelta[1][2] + ChiInv[0][0]*ChiInv[1][1]*ChiN[0][1]*KrDelta[1][2] +
         ChiInv[0][0]*ChiInv[2][1]*ChiN[0][2]*KrDelta[1][2] + ChiInv[0][1]*ChiInv[0][1]*ChiN[1][0]*KrDelta[1][2] +
         ChiInv[0][1]*ChiInv[1][1]*ChiN[1][1]*KrDelta[1][2] + ChiInv[0][1]*ChiInv[2][1]*ChiN[1][2]*KrDelta[1][2] +
         ChiInv[0][1]*ChiInv[0][2]*ChiN[2][0]*KrDelta[1][2] + ChiInv[0][2]*ChiInv[1][1]*ChiN[2][1]*KrDelta[1][2] +
         ChiInv[0][2]*ChiInv[2][1]*ChiN[2][2]*KrDelta[1][2]);//*w[2][1] +

TChi_1(8,1)=(ChiInv[0][0]*ChiInv[0][1]*ChiN[0][0]*KrDelta[2][2] + ChiInv[0][0]*ChiInv[1][1]*ChiN[0][1]*KrDelta[2][2] +
         ChiInv[0][0]*ChiInv[2][1]*ChiN[0][2]*KrDelta[2][2] + ChiInv[0][1]*ChiInv[0][1]*ChiN[1][0]*KrDelta[2][2] +
         ChiInv[0][1]*ChiInv[1][1]*ChiN[1][1]*KrDelta[2][2] + ChiInv[0][1]*ChiInv[2][1]*ChiN[1][2]*KrDelta[2][2] +
         ChiInv[0][1]*ChiInv[0][2]*ChiN[2][0]*KrDelta[2][2] + ChiInv[0][2]*ChiInv[1][1]*ChiN[2][1]*KrDelta[2][2] +
         ChiInv[0][2]*ChiInv[2][1]*ChiN[2][2]*KrDelta[2][2]);//*w[2][2]

TChi_1(0,2)=(ChiInv[0][0]*ChiInv[0][2]*ChiN[0][0]*KrDelta[0][0] + ChiInv[0][0]*ChiInv[1][2]*ChiN[0][1]*KrDelta[0][0] +
         ChiInv[0][0]*ChiInv[2][2]*ChiN[0][2]*KrDelta[0][0] + ChiInv[0][1]*ChiInv[0][2]*ChiN[1][0]*KrDelta[0][0] +
         ChiInv[0][1]*ChiInv[1][2]*ChiN[1][1]*KrDelta[0][0] + ChiInv[0][1]*ChiInv[2][2]*ChiN[1][2]*KrDelta[0][0] +
         ChiInv[0][2]*ChiInv[0][2]*ChiN[2][0]*KrDelta[0][0] + ChiInv[0][2]*ChiInv[1][2]*ChiN[2][1]*KrDelta[0][0] +
         ChiInv[0][2]*ChiInv[2][2]*ChiN[2][2]*KrDelta[0][0]);//*w[0][0] +
TChi_1(1,2)=(ChiInv[0][0]*ChiInv[0][2]*ChiN[0][0]*KrDelta[1][0] + ChiInv[0][0]*ChiInv[1][2]*ChiN[0][1]*KrDelta[1][0] +
         ChiInv[0][0]*ChiInv[2][2]*ChiN[0][2]*KrDelta[1][0] + ChiInv[0][1]*ChiInv[0][2]*ChiN[1][0]*KrDelta[1][0] +
         ChiInv[0][1]*ChiInv[1][2]*ChiN[1][1]*KrDelta[1][0] + ChiInv[0][1]*ChiInv[2][2]*ChiN[1][2]*KrDelta[1][0] +
         ChiInv[0][2]*ChiInv[0][2]*ChiN[2][0]*KrDelta[1][0] + ChiInv[0][2]*ChiInv[1][2]*ChiN[2][1]*KrDelta[1][0] +
         ChiInv[0][2]*ChiInv[2][2]*ChiN[2][2]*KrDelta[1][0]);//*w[0][1] +
TChi_1(2,2)= (ChiInv[0][0]*ChiInv[0][2]*ChiN[0][0]*KrDelta[2][0] + ChiInv[0][0]*ChiInv[1][2]*ChiN[0][1]*KrDelta[2][0] +
         ChiInv[0][0]*ChiInv[2][2]*ChiN[0][2]*KrDelta[2][0] + ChiInv[0][1]*ChiInv[0][2]*ChiN[1][0]*KrDelta[2][0] +
         ChiInv[0][1]*ChiInv[1][2]*ChiN[1][1]*KrDelta[2][0] + ChiInv[0][1]*ChiInv[2][2]*ChiN[1][2]*KrDelta[2][0] +
         ChiInv[0][2]*ChiInv[0][2]*ChiN[2][0]*KrDelta[2][0] + ChiInv[0][2]*ChiInv[1][2]*ChiN[2][1]*KrDelta[2][0] +
         ChiInv[0][2]*ChiInv[2][2]*ChiN[2][2]*KrDelta[2][0]);//*w[0][2] +
TChi_1(3,2)=(ChiInv[0][0]*ChiInv[0][2]*ChiN[0][0]*KrDelta[0][1] + ChiInv[0][0]*ChiInv[1][2]*ChiN[0][1]*KrDelta[0][1] +
         ChiInv[0][0]*ChiInv[2][2]*ChiN[0][2]*KrDelta[0][1] + ChiInv[0][1]*ChiInv[0][2]*ChiN[1][0]*KrDelta[0][1] +
         ChiInv[0][1]*ChiInv[1][2]*ChiN[1][1]*KrDelta[0][1] + ChiInv[0][1]*ChiInv[2][2]*ChiN[1][2]*KrDelta[0][1] +
         ChiInv[0][2]*ChiInv[0][2]*ChiN[2][0]*KrDelta[0][1] + ChiInv[0][2]*ChiInv[1][2]*ChiN[2][1]*KrDelta[0][1] +
         ChiInv[0][2]*ChiInv[2][2]*ChiN[2][2]*KrDelta[0][1]);//*w[1][0] +
TChi_1(4,2)=(ChiInv[0][0]*ChiInv[0][2]*ChiN[0][0]*KrDelta[1][1] + ChiInv[0][0]*ChiInv[1][2]*ChiN[0][1]*KrDelta[1][1] +
         ChiInv[0][0]*ChiInv[2][2]*ChiN[0][2]*KrDelta[1][1] + ChiInv[0][1]*ChiInv[0][2]*ChiN[1][0]*KrDelta[1][1] +
         ChiInv[0][1]*ChiInv[1][2]*ChiN[1][1]*KrDelta[1][1] + ChiInv[0][1]*ChiInv[2][2]*ChiN[1][2]*KrDelta[1][1] +
         ChiInv[0][2]*ChiInv[0][2]*ChiN[2][0]*KrDelta[1][1] + ChiInv[0][2]*ChiInv[1][2]*ChiN[2][1]*KrDelta[1][1] +
         ChiInv[0][2]*ChiInv[2][2]*ChiN[2][2]*KrDelta[1][1]);//*w[1][1] +
TChi_1(5,2)=(ChiInv[0][0]*ChiInv[0][2]*ChiN[0][0]*KrDelta[2][1] + ChiInv[0][0]*ChiInv[1][2]*ChiN[0][1]*KrDelta[2][1] +
         ChiInv[0][0]*ChiInv[2][2]*ChiN[0][2]*KrDelta[2][1] + ChiInv[0][1]*ChiInv[0][2]*ChiN[1][0]*KrDelta[2][1] +
         ChiInv[0][1]*ChiInv[1][2]*ChiN[1][1]*KrDelta[2][1] + ChiInv[0][1]*ChiInv[2][2]*ChiN[1][2]*KrDelta[2][1] +
         ChiInv[0][2]*ChiInv[0][2]*ChiN[2][0]*KrDelta[2][1] + ChiInv[0][2]*ChiInv[1][2]*ChiN[2][1]*KrDelta[2][1] +
         ChiInv[0][2]*ChiInv[2][2]*ChiN[2][2]*KrDelta[2][1]);//*w[1][2] +
TChi_1(6,2)=(ChiInv[0][0]*ChiInv[0][2]*ChiN[0][0]*KrDelta[0][2] + ChiInv[0][0]*ChiInv[1][2]*ChiN[0][1]*KrDelta[0][2] +
         ChiInv[0][0]*ChiInv[2][2]*ChiN[0][2]*KrDelta[0][2] + ChiInv[0][1]*ChiInv[0][2]*ChiN[1][0]*KrDelta[0][2] +
         ChiInv[0][1]*ChiInv[1][2]*ChiN[1][1]*KrDelta[0][2] + ChiInv[0][1]*ChiInv[2][2]*ChiN[1][2]*KrDelta[0][2] +
         ChiInv[0][2]*ChiInv[0][2]*ChiN[2][0]*KrDelta[0][2] + ChiInv[0][2]*ChiInv[1][2]*ChiN[2][1]*KrDelta[0][2] +
         ChiInv[0][2]*ChiInv[2][2]*ChiN[2][2]*KrDelta[0][2]);//*w[2][0] +
TChi_1(7,2)=(ChiInv[0][0]*ChiInv[0][2]*ChiN[0][0]*KrDelta[1][2] + ChiInv[0][0]*ChiInv[1][2]*ChiN[0][1]*KrDelta[1][2] +
         ChiInv[0][0]*ChiInv[2][2]*ChiN[0][2]*KrDelta[1][2] + ChiInv[0][1]*ChiInv[0][2]*ChiN[1][0]*KrDelta[1][2] +
         ChiInv[0][1]*ChiInv[1][2]*ChiN[1][1]*KrDelta[1][2] + ChiInv[0][1]*ChiInv[2][2]*ChiN[1][2]*KrDelta[1][2] +
         ChiInv[0][2]*ChiInv[0][2]*ChiN[2][0]*KrDelta[1][2] + ChiInv[0][2]*ChiInv[1][2]*ChiN[2][1]*KrDelta[1][2] +
         ChiInv[0][2]*ChiInv[2][2]*ChiN[2][2]*KrDelta[1][2]);//*w[2][1] +
TChi_1(8,2)=(ChiInv[0][0]*ChiInv[0][2]*ChiN[0][0]*KrDelta[2][2] + ChiInv[0][0]*ChiInv[1][2]*ChiN[0][1]*KrDelta[2][2] +
         ChiInv[0][0]*ChiInv[2][2]*ChiN[0][2]*KrDelta[2][2] + ChiInv[0][1]*ChiInv[0][2]*ChiN[1][0]*KrDelta[2][2] +
         ChiInv[0][1]*ChiInv[1][2]*ChiN[1][1]*KrDelta[2][2] + ChiInv[0][1]*ChiInv[2][2]*ChiN[1][2]*KrDelta[2][2] +
         ChiInv[0][2]*ChiInv[0][2]*ChiN[2][0]*KrDelta[2][2] + ChiInv[0][2]*ChiInv[1][2]*ChiN[2][1]*KrDelta[2][2] +
         ChiInv[0][2]*ChiInv[2][2]*ChiN[2][2]*KrDelta[2][2]);//*w[2][2])


TChi_1(0,3)=(ChiInv[0][0]*ChiInv[1][0]*ChiN[0][0]*KrDelta[0][0] + ChiInv[1][0]*ChiInv[1][0]*ChiN[0][1]*KrDelta[0][0] +
            ChiInv[1][0]*ChiInv[2][0]*ChiN[0][2]*KrDelta[0][0] + ChiInv[0][0]*ChiInv[1][1]*ChiN[1][0]*KrDelta[0][0] +
            ChiInv[1][0]*ChiInv[1][1]*ChiN[1][1]*KrDelta[0][0] + ChiInv[1][1]*ChiInv[2][0]*ChiN[1][2]*KrDelta[0][0] +
            ChiInv[0][0]*ChiInv[1][2]*ChiN[2][0]*KrDelta[0][0] + ChiInv[1][0]*ChiInv[1][2]*ChiN[2][1]*KrDelta[0][0] +
            ChiInv[1][2]*ChiInv[2][0]*ChiN[2][2]*KrDelta[0][0]);//*w[0][0] +
    TChi_1(1,3)=(ChiInv[0][0]*ChiInv[1][0]*ChiN[0][0]*KrDelta[1][0] + ChiInv[1][0]*ChiInv[1][0]*ChiN[0][1]*KrDelta[1][0] +
         ChiInv[1][0]*ChiInv[2][0]*ChiN[0][2]*KrDelta[1][0] + ChiInv[0][0]*ChiInv[1][1]*ChiN[1][0]*KrDelta[1][0] +
         ChiInv[1][0]*ChiInv[1][1]*ChiN[1][1]*KrDelta[1][0] + ChiInv[1][1]*ChiInv[2][0]*ChiN[1][2]*KrDelta[1][0] +
         ChiInv[0][0]*ChiInv[1][2]*ChiN[2][0]*KrDelta[1][0] + ChiInv[1][0]*ChiInv[1][2]*ChiN[2][1]*KrDelta[1][0] +
         ChiInv[1][2]*ChiInv[2][0]*ChiN[2][2]*KrDelta[1][0]);//*w[0][1] +
    TChi_1(2,3)=(ChiInv[0][0]*ChiInv[1][0]*ChiN[0][0]*KrDelta[2][0] + ChiInv[1][0]*ChiInv[1][0]*ChiN[0][1]*KrDelta[2][0] +
         ChiInv[1][0]*ChiInv[2][0]*ChiN[0][2]*KrDelta[2][0] + ChiInv[0][0]*ChiInv[1][1]*ChiN[1][0]*KrDelta[2][0] +
         ChiInv[1][0]*ChiInv[1][1]*ChiN[1][1]*KrDelta[2][0] + ChiInv[1][1]*ChiInv[2][0]*ChiN[1][2]*KrDelta[2][0] +
         ChiInv[0][0]*ChiInv[1][2]*ChiN[2][0]*KrDelta[2][0] + ChiInv[1][0]*ChiInv[1][2]*ChiN[2][1]*KrDelta[2][0] +
         ChiInv[1][2]*ChiInv[2][0]*ChiN[2][2]*KrDelta[2][0]);//*w[0][2] +
    TChi_1(3,3)=(ChiInv[0][0]*ChiInv[1][0]*ChiN[0][0]*KrDelta[0][1] + ChiInv[1][0]*ChiInv[1][0]*ChiN[0][1]*KrDelta[0][1] +
         ChiInv[1][0]*ChiInv[2][0]*ChiN[0][2]*KrDelta[0][1] + ChiInv[0][0]*ChiInv[1][1]*ChiN[1][0]*KrDelta[0][1] +
         ChiInv[1][0]*ChiInv[1][1]*ChiN[1][1]*KrDelta[0][1] + ChiInv[1][1]*ChiInv[2][0]*ChiN[1][2]*KrDelta[0][1] +
         ChiInv[0][0]*ChiInv[1][2]*ChiN[2][0]*KrDelta[0][1] + ChiInv[1][0]*ChiInv[1][2]*ChiN[2][1]*KrDelta[0][1] +
         ChiInv[1][2]*ChiInv[2][0]*ChiN[2][2]*KrDelta[0][1]);//*w[1][0] +
    TChi_1(4,3)=(ChiInv[0][0]*ChiInv[1][0]*ChiN[0][0]*KrDelta[1][1] + ChiInv[1][0]*ChiInv[1][0]*ChiN[0][1]*KrDelta[1][1] +
         ChiInv[1][0]*ChiInv[2][0]*ChiN[0][2]*KrDelta[1][1] + ChiInv[0][0]*ChiInv[1][1]*ChiN[1][0]*KrDelta[1][1] +
         ChiInv[1][0]*ChiInv[1][1]*ChiN[1][1]*KrDelta[1][1] + ChiInv[1][1]*ChiInv[2][0]*ChiN[1][2]*KrDelta[1][1] +
         ChiInv[0][0]*ChiInv[1][2]*ChiN[2][0]*KrDelta[1][1] + ChiInv[1][0]*ChiInv[1][2]*ChiN[2][1]*KrDelta[1][1] +
         ChiInv[1][2]*ChiInv[2][0]*ChiN[2][2]*KrDelta[1][1]);//*w[1][1] +
    TChi_1(5,3)=(ChiInv[0][0]*ChiInv[1][0]*ChiN[0][0]*KrDelta[2][1] + ChiInv[1][0]*ChiInv[1][0]*ChiN[0][1]*KrDelta[2][1] +
         ChiInv[1][0]*ChiInv[2][0]*ChiN[0][2]*KrDelta[2][1] + ChiInv[0][0]*ChiInv[1][1]*ChiN[1][0]*KrDelta[2][1] +
         ChiInv[1][0]*ChiInv[1][1]*ChiN[1][1]*KrDelta[2][1] + ChiInv[1][1]*ChiInv[2][0]*ChiN[1][2]*KrDelta[2][1] +
         ChiInv[0][0]*ChiInv[1][2]*ChiN[2][0]*KrDelta[2][1] + ChiInv[1][0]*ChiInv[1][2]*ChiN[2][1]*KrDelta[2][1] +
         ChiInv[1][2]*ChiInv[2][0]*ChiN[2][2]*KrDelta[2][1]);//*w[1][2] +
    TChi_1(6,3)=(ChiInv[0][0]*ChiInv[1][0]*ChiN[0][0]*KrDelta[0][2] + ChiInv[1][0]*ChiInv[1][0]*ChiN[0][1]*KrDelta[0][2] +
         ChiInv[1][0]*ChiInv[2][0]*ChiN[0][2]*KrDelta[0][2] + ChiInv[0][0]*ChiInv[1][1]*ChiN[1][0]*KrDelta[0][2] +
         ChiInv[1][0]*ChiInv[1][1]*ChiN[1][1]*KrDelta[0][2] + ChiInv[1][1]*ChiInv[2][0]*ChiN[1][2]*KrDelta[0][2] +
         ChiInv[0][0]*ChiInv[1][2]*ChiN[2][0]*KrDelta[0][2] + ChiInv[1][0]*ChiInv[1][2]*ChiN[2][1]*KrDelta[0][2] +
         ChiInv[1][2]*ChiInv[2][0]*ChiN[2][2]*KrDelta[0][2]);//*w[2][0] +
    TChi_1(7,3)=(ChiInv[0][0]*ChiInv[1][0]*ChiN[0][0]*KrDelta[1][2] + ChiInv[1][0]*ChiInv[1][0]*ChiN[0][1]*KrDelta[1][2] +
         ChiInv[1][0]*ChiInv[2][0]*ChiN[0][2]*KrDelta[1][2] + ChiInv[0][0]*ChiInv[1][1]*ChiN[1][0]*KrDelta[1][2] +
         ChiInv[1][0]*ChiInv[1][1]*ChiN[1][1]*KrDelta[1][2] + ChiInv[1][1]*ChiInv[2][0]*ChiN[1][2]*KrDelta[1][2] +
         ChiInv[0][0]*ChiInv[1][2]*ChiN[2][0]*KrDelta[1][2] + ChiInv[1][0]*ChiInv[1][2]*ChiN[2][1]*KrDelta[1][2] +
         ChiInv[1][2]*ChiInv[2][0]*ChiN[2][2]*KrDelta[1][2]);//*w[2][1] +
    TChi_1(8,3)=(ChiInv[0][0]*ChiInv[1][0]*ChiN[0][0]*KrDelta[2][2] + ChiInv[1][0]*ChiInv[1][0]*ChiN[0][1]*KrDelta[2][2] +
         ChiInv[1][0]*ChiInv[2][0]*ChiN[0][2]*KrDelta[2][2] + ChiInv[0][0]*ChiInv[1][1]*ChiN[1][0]*KrDelta[2][2] +
         ChiInv[1][0]*ChiInv[1][1]*ChiN[1][1]*KrDelta[2][2] + ChiInv[1][1]*ChiInv[2][0]*ChiN[1][2]*KrDelta[2][2] +
         ChiInv[0][0]*ChiInv[1][2]*ChiN[2][0]*KrDelta[2][2] + ChiInv[1][0]*ChiInv[1][2]*ChiN[2][1]*KrDelta[2][2] +
         ChiInv[1][2]*ChiInv[2][0]*ChiN[2][2]*KrDelta[2][2]);//*w[2][2]

    TChi_1(0,4)=(ChiInv[0][1]*ChiInv[1][0]*ChiN[0][0]*KrDelta[0][0] + ChiInv[1][0]*ChiInv[1][1]*ChiN[0][1]*KrDelta[0][0] +
          ChiInv[1][0]*ChiInv[2][1]*ChiN[0][2]*KrDelta[0][0] + ChiInv[0][1]*ChiInv[1][1]*ChiN[1][0]*KrDelta[0][0] +
          ChiInv[1][1]*ChiInv[1][1]*ChiN[1][1]*KrDelta[0][0] + ChiInv[1][1]*ChiInv[2][1]*ChiN[1][2]*KrDelta[0][0] +
          ChiInv[0][1]*ChiInv[1][2]*ChiN[2][0]*KrDelta[0][0] + ChiInv[1][1]*ChiInv[1][2]*ChiN[2][1]*KrDelta[0][0] +
          ChiInv[1][2]*ChiInv[2][1]*ChiN[2][2]*KrDelta[0][0]);//*w[0][0] +

    TChi_1(1,4)=(ChiInv[0][1]*ChiInv[1][0]*ChiN[0][0]*KrDelta[1][0] + ChiInv[1][0]*ChiInv[1][1]*ChiN[0][1]*KrDelta[1][0] +
          ChiInv[1][0]*ChiInv[2][1]*ChiN[0][2]*KrDelta[1][0] + ChiInv[0][1]*ChiInv[1][1]*ChiN[1][0]*KrDelta[1][0] +
          ChiInv[1][1]*ChiInv[1][1]*ChiN[1][1]*KrDelta[1][0] + ChiInv[1][1]*ChiInv[2][1]*ChiN[1][2]*KrDelta[1][0] +
          ChiInv[0][1]*ChiInv[1][2]*ChiN[2][0]*KrDelta[1][0] + ChiInv[1][1]*ChiInv[1][2]*ChiN[2][1]*KrDelta[1][0] +
          ChiInv[1][2]*ChiInv[2][1]*ChiN[2][2]*KrDelta[1][0]);//*w[0][1] +
    TChi_1(2,4)=(ChiInv[0][1]*ChiInv[1][0]*ChiN[0][0]*KrDelta[2][0] + ChiInv[1][0]*ChiInv[1][1]*ChiN[0][1]*KrDelta[2][0] +
          ChiInv[1][0]*ChiInv[2][1]*ChiN[0][2]*KrDelta[2][0] + ChiInv[0][1]*ChiInv[1][1]*ChiN[1][0]*KrDelta[2][0] +
          ChiInv[1][1]*ChiInv[1][1]*ChiN[1][1]*KrDelta[2][0] + ChiInv[1][1]*ChiInv[2][1]*ChiN[1][2]*KrDelta[2][0] +
          ChiInv[0][1]*ChiInv[1][2]*ChiN[2][0]*KrDelta[2][0] + ChiInv[1][1]*ChiInv[1][2]*ChiN[2][1]*KrDelta[2][0] +
          ChiInv[1][2]*ChiInv[2][1]*ChiN[2][2]*KrDelta[2][0]);//*w[0][2] +
    TChi_1(3,4)=(ChiInv[0][1]*ChiInv[1][0]*ChiN[0][0]*KrDelta[0][1] + ChiInv[1][0]*ChiInv[1][1]*ChiN[0][1]*KrDelta[0][1] +
          ChiInv[1][0]*ChiInv[2][1]*ChiN[0][2]*KrDelta[0][1] + ChiInv[0][1]*ChiInv[1][1]*ChiN[1][0]*KrDelta[0][1] +
          ChiInv[1][1]*ChiInv[1][1]*ChiN[1][1]*KrDelta[0][1] + ChiInv[1][1]*ChiInv[2][1]*ChiN[1][2]*KrDelta[0][1] +
          ChiInv[0][1]*ChiInv[1][2]*ChiN[2][0]*KrDelta[0][1] + ChiInv[1][1]*ChiInv[1][2]*ChiN[2][1]*KrDelta[0][1] +
          ChiInv[1][2]*ChiInv[2][1]*ChiN[2][2]*KrDelta[0][1]);//*w[1][0] +
    TChi_1(4,4)=(ChiInv[0][1]*ChiInv[1][0]*ChiN[0][0]*KrDelta[1][1] + ChiInv[1][0]*ChiInv[1][1]*ChiN[0][1]*KrDelta[1][1] +
          ChiInv[1][0]*ChiInv[2][1]*ChiN[0][2]*KrDelta[1][1] + ChiInv[0][1]*ChiInv[1][1]*ChiN[1][0]*KrDelta[1][1] +
          ChiInv[1][1]*ChiInv[1][1]*ChiN[1][1]*KrDelta[1][1] + ChiInv[1][1]*ChiInv[2][1]*ChiN[1][2]*KrDelta[1][1] +
          ChiInv[0][1]*ChiInv[1][2]*ChiN[2][0]*KrDelta[1][1] + ChiInv[1][1]*ChiInv[1][2]*ChiN[2][1]*KrDelta[1][1] +
          ChiInv[1][2]*ChiInv[2][1]*ChiN[2][2]*KrDelta[1][1]);//*w[1][1] +
    TChi_1(5,4)=(ChiInv[0][1]*ChiInv[1][0]*ChiN[0][0]*KrDelta[2][1] + ChiInv[1][0]*ChiInv[1][1]*ChiN[0][1]*KrDelta[2][1] +
          ChiInv[1][0]*ChiInv[2][1]*ChiN[0][2]*KrDelta[2][1] + ChiInv[0][1]*ChiInv[1][1]*ChiN[1][0]*KrDelta[2][1] +
          ChiInv[1][1]*ChiInv[1][1]*ChiN[1][1]*KrDelta[2][1] + ChiInv[1][1]*ChiInv[2][1]*ChiN[1][2]*KrDelta[2][1] +
          ChiInv[0][1]*ChiInv[1][2]*ChiN[2][0]*KrDelta[2][1] + ChiInv[1][1]*ChiInv[1][2]*ChiN[2][1]*KrDelta[2][1] +
          ChiInv[1][2]*ChiInv[2][1]*ChiN[2][2]*KrDelta[2][1]);//*w[1][2] +
    TChi_1(6,4)=(ChiInv[0][1]*ChiInv[1][0]*ChiN[0][0]*KrDelta[0][2] + ChiInv[1][0]*ChiInv[1][1]*ChiN[0][1]*KrDelta[0][2] +
          ChiInv[1][0]*ChiInv[2][1]*ChiN[0][2]*KrDelta[0][2] + ChiInv[0][1]*ChiInv[1][1]*ChiN[1][0]*KrDelta[0][2] +
          ChiInv[1][1]*ChiInv[1][1]*ChiN[1][1]*KrDelta[0][2] + ChiInv[1][1]*ChiInv[2][1]*ChiN[1][2]*KrDelta[0][2] +
          ChiInv[0][1]*ChiInv[1][2]*ChiN[2][0]*KrDelta[0][2] + ChiInv[1][1]*ChiInv[1][2]*ChiN[2][1]*KrDelta[0][2] +
          ChiInv[1][2]*ChiInv[2][1]*ChiN[2][2]*KrDelta[0][2]);//*w[2][0] +
    TChi_1(7,4)=(ChiInv[0][1]*ChiInv[1][0]*ChiN[0][0]*KrDelta[1][2] + ChiInv[1][0]*ChiInv[1][1]*ChiN[0][1]*KrDelta[1][2] +
          ChiInv[1][0]*ChiInv[2][1]*ChiN[0][2]*KrDelta[1][2] + ChiInv[0][1]*ChiInv[1][1]*ChiN[1][0]*KrDelta[1][2] +
          ChiInv[1][1]*ChiInv[1][1]*ChiN[1][1]*KrDelta[1][2] + ChiInv[1][1]*ChiInv[2][1]*ChiN[1][2]*KrDelta[1][2] +
          ChiInv[0][1]*ChiInv[1][2]*ChiN[2][0]*KrDelta[1][2] + ChiInv[1][1]*ChiInv[1][2]*ChiN[2][1]*KrDelta[1][2] +
          ChiInv[1][2]*ChiInv[2][1]*ChiN[2][2]*KrDelta[1][2]);//*w[2][1] +
    TChi_1(8,4)=(ChiInv[0][1]*ChiInv[1][0]*ChiN[0][0]*KrDelta[2][2] + ChiInv[1][0]*ChiInv[1][1]*ChiN[0][1]*KrDelta[2][2] +
          ChiInv[1][0]*ChiInv[2][1]*ChiN[0][2]*KrDelta[2][2] + ChiInv[0][1]*ChiInv[1][1]*ChiN[1][0]*KrDelta[2][2] +
          ChiInv[1][1]*ChiInv[1][1]*ChiN[1][1]*KrDelta[2][2] + ChiInv[1][1]*ChiInv[2][1]*ChiN[1][2]*KrDelta[2][2] +
          ChiInv[0][1]*ChiInv[1][2]*ChiN[2][0]*KrDelta[2][2] + ChiInv[1][1]*ChiInv[1][2]*ChiN[2][1]*KrDelta[2][2] +
          ChiInv[1][2]*ChiInv[2][1]*ChiN[2][2]*KrDelta[2][2]);//*w[2][2] +

     TChi_1(0,5)=(ChiInv[0][2]*ChiInv[1][0]*ChiN[0][0]*KrDelta[0][0] + ChiInv[1][0]*ChiInv[1][2]*ChiN[0][1]*KrDelta[0][0] +
          ChiInv[1][0]*ChiInv[2][2]*ChiN[0][2]*KrDelta[0][0] + ChiInv[0][2]*ChiInv[1][1]*ChiN[1][0]*KrDelta[0][0] +
          ChiInv[1][1]*ChiInv[1][2]*ChiN[1][1]*KrDelta[0][0] + ChiInv[1][1]*ChiInv[2][2]*ChiN[1][2]*KrDelta[0][0] +
          ChiInv[0][2]*ChiInv[1][2]*ChiN[2][0]*KrDelta[0][0] + ChiInv[1][2]*ChiInv[1][2]*ChiN[2][1]*KrDelta[0][0] +
          ChiInv[1][2]*ChiInv[2][2]*ChiN[2][2]*KrDelta[0][0]);//*w[0][0] +
     TChi_1(1,5)=(ChiInv[0][2]*ChiInv[1][0]*ChiN[0][0]*KrDelta[1][0] + ChiInv[1][0]*ChiInv[1][2]*ChiN[0][1]*KrDelta[1][0] +
          ChiInv[1][0]*ChiInv[2][2]*ChiN[0][2]*KrDelta[1][0] + ChiInv[0][2]*ChiInv[1][1]*ChiN[1][0]*KrDelta[1][0] +
          ChiInv[1][1]*ChiInv[1][2]*ChiN[1][1]*KrDelta[1][0] + ChiInv[1][1]*ChiInv[2][2]*ChiN[1][2]*KrDelta[1][0] +
          ChiInv[0][2]*ChiInv[1][2]*ChiN[2][0]*KrDelta[1][0] + ChiInv[1][2]*ChiInv[1][2]*ChiN[2][1]*KrDelta[1][0] +
          ChiInv[1][2]*ChiInv[2][2]*ChiN[2][2]*KrDelta[1][0]);//*w[0][1] +
     TChi_1(2,5)=(ChiInv[0][2]*ChiInv[1][0]*ChiN[0][0]*KrDelta[2][0] + ChiInv[1][0]*ChiInv[1][2]*ChiN[0][1]*KrDelta[2][0] +
          ChiInv[1][0]*ChiInv[2][2]*ChiN[0][2]*KrDelta[2][0] + ChiInv[0][2]*ChiInv[1][1]*ChiN[1][0]*KrDelta[2][0] +
          ChiInv[1][1]*ChiInv[1][2]*ChiN[1][1]*KrDelta[2][0] + ChiInv[1][1]*ChiInv[2][2]*ChiN[1][2]*KrDelta[2][0] +
          ChiInv[0][2]*ChiInv[1][2]*ChiN[2][0]*KrDelta[2][0] + ChiInv[1][2]*ChiInv[1][2]*ChiN[2][1]*KrDelta[2][0] +
          ChiInv[1][2]*ChiInv[2][2]*ChiN[2][2]*KrDelta[2][0]);//*w[0][2] +
     TChi_1(3,5)=(ChiInv[0][2]*ChiInv[1][0]*ChiN[0][0]*KrDelta[0][1] + ChiInv[1][0]*ChiInv[1][2]*ChiN[0][1]*KrDelta[0][1] +
          ChiInv[1][0]*ChiInv[2][2]*ChiN[0][2]*KrDelta[0][1] + ChiInv[0][2]*ChiInv[1][1]*ChiN[1][0]*KrDelta[0][1] +
          ChiInv[1][1]*ChiInv[1][2]*ChiN[1][1]*KrDelta[0][1] + ChiInv[1][1]*ChiInv[2][2]*ChiN[1][2]*KrDelta[0][1] +
          ChiInv[0][2]*ChiInv[1][2]*ChiN[2][0]*KrDelta[0][1] + ChiInv[1][2]*ChiInv[1][2]*ChiN[2][1]*KrDelta[0][1] +
          ChiInv[1][2]*ChiInv[2][2]*ChiN[2][2]*KrDelta[0][1]);//*w[1][0] +
     TChi_1(4,5)=(ChiInv[0][2]*ChiInv[1][0]*ChiN[0][0]*KrDelta[1][1] + ChiInv[1][0]*ChiInv[1][2]*ChiN[0][1]*KrDelta[1][1] +
          ChiInv[1][0]*ChiInv[2][2]*ChiN[0][2]*KrDelta[1][1] + ChiInv[0][2]*ChiInv[1][1]*ChiN[1][0]*KrDelta[1][1] +
          ChiInv[1][1]*ChiInv[1][2]*ChiN[1][1]*KrDelta[1][1] + ChiInv[1][1]*ChiInv[2][2]*ChiN[1][2]*KrDelta[1][1] +
          ChiInv[0][2]*ChiInv[1][2]*ChiN[2][0]*KrDelta[1][1] + ChiInv[1][2]*ChiInv[1][2]*ChiN[2][1]*KrDelta[1][1] +
          ChiInv[1][2]*ChiInv[2][2]*ChiN[2][2]*KrDelta[1][1]);//*w[1][1] +
     TChi_1(5,5)=(ChiInv[0][2]*ChiInv[1][0]*ChiN[0][0]*KrDelta[2][1] + ChiInv[1][0]*ChiInv[1][2]*ChiN[0][1]*KrDelta[2][1] +
          ChiInv[1][0]*ChiInv[2][2]*ChiN[0][2]*KrDelta[2][1] + ChiInv[0][2]*ChiInv[1][1]*ChiN[1][0]*KrDelta[2][1] +
          ChiInv[1][1]*ChiInv[1][2]*ChiN[1][1]*KrDelta[2][1] + ChiInv[1][1]*ChiInv[2][2]*ChiN[1][2]*KrDelta[2][1] +
          ChiInv[0][2]*ChiInv[1][2]*ChiN[2][0]*KrDelta[2][1] + ChiInv[1][2]*ChiInv[1][2]*ChiN[2][1]*KrDelta[2][1] +
          ChiInv[1][2]*ChiInv[2][2]*ChiN[2][2]*KrDelta[2][1]);//*w[1][2] +
     TChi_1(6,5)=(ChiInv[0][2]*ChiInv[1][0]*ChiN[0][0]*KrDelta[0][2] + ChiInv[1][0]*ChiInv[1][2]*ChiN[0][1]*KrDelta[0][2] +
          ChiInv[1][0]*ChiInv[2][2]*ChiN[0][2]*KrDelta[0][2] + ChiInv[0][2]*ChiInv[1][1]*ChiN[1][0]*KrDelta[0][2] +
          ChiInv[1][1]*ChiInv[1][2]*ChiN[1][1]*KrDelta[0][2] + ChiInv[1][1]*ChiInv[2][2]*ChiN[1][2]*KrDelta[0][2] +
          ChiInv[0][2]*ChiInv[1][2]*ChiN[2][0]*KrDelta[0][2] + ChiInv[1][2]*ChiInv[1][2]*ChiN[2][1]*KrDelta[0][2] +
          ChiInv[1][2]*ChiInv[2][2]*ChiN[2][2]*KrDelta[0][2]);//*w[2][0] +
     TChi_1(7,5)=(ChiInv[0][2]*ChiInv[1][0]*ChiN[0][0]*KrDelta[1][2] + ChiInv[1][0]*ChiInv[1][2]*ChiN[0][1]*KrDelta[1][2] +
          ChiInv[1][0]*ChiInv[2][2]*ChiN[0][2]*KrDelta[1][2] + ChiInv[0][2]*ChiInv[1][1]*ChiN[1][0]*KrDelta[1][2] +
          ChiInv[1][1]*ChiInv[1][2]*ChiN[1][1]*KrDelta[1][2] + ChiInv[1][1]*ChiInv[2][2]*ChiN[1][2]*KrDelta[1][2] +
          ChiInv[0][2]*ChiInv[1][2]*ChiN[2][0]*KrDelta[1][2] + ChiInv[1][2]*ChiInv[1][2]*ChiN[2][1]*KrDelta[1][2] +
          ChiInv[1][2]*ChiInv[2][2]*ChiN[2][2]*KrDelta[1][2]);//*w[2][1] +
     TChi_1(8,5)=(ChiInv[0][2]*ChiInv[1][0]*ChiN[0][0]*KrDelta[2][2] + ChiInv[1][0]*ChiInv[1][2]*ChiN[0][1]*KrDelta[2][2] +
          ChiInv[1][0]*ChiInv[2][2]*ChiN[0][2]*KrDelta[2][2] + ChiInv[0][2]*ChiInv[1][1]*ChiN[1][0]*KrDelta[2][2] +
          ChiInv[1][1]*ChiInv[1][2]*ChiN[1][1]*KrDelta[2][2] + ChiInv[1][1]*ChiInv[2][2]*ChiN[1][2]*KrDelta[2][2] +
          ChiInv[0][2]*ChiInv[1][2]*ChiN[2][0]*KrDelta[2][2] + ChiInv[1][2]*ChiInv[1][2]*ChiN[2][1]*KrDelta[2][2] +
          ChiInv[1][2]*ChiInv[2][2]*ChiN[2][2]*KrDelta[2][2]);//*w[2][2]

     TChi_1(0,6)=(ChiInv[0][0]*ChiInv[2][0]*ChiN[0][0]*KrDelta[0][0] + ChiInv[1][0]*ChiInv[2][0]*ChiN[0][1]*KrDelta[0][0] +
          ChiInv[2][0]*ChiInv[2][0]*ChiN[0][2]*KrDelta[0][0] + ChiInv[0][0]*ChiInv[2][1]*ChiN[1][0]*KrDelta[0][0] +
          ChiInv[1][0]*ChiInv[2][1]*ChiN[1][1]*KrDelta[0][0] + ChiInv[2][0]*ChiInv[2][1]*ChiN[1][2]*KrDelta[0][0] +
          ChiInv[0][0]*ChiInv[2][2]*ChiN[2][0]*KrDelta[0][0] + ChiInv[1][0]*ChiInv[2][2]*ChiN[2][1]*KrDelta[0][0] +
          ChiInv[2][0]*ChiInv[2][2]*ChiN[2][2]*KrDelta[0][0]);//*w[0][0] +
     TChi_1(1,6)=(ChiInv[0][0]*ChiInv[2][0]*ChiN[0][0]*KrDelta[1][0] + ChiInv[1][0]*ChiInv[2][0]*ChiN[0][1]*KrDelta[1][0] +
          ChiInv[2][0]*ChiInv[2][0]*ChiN[0][2]*KrDelta[1][0] + ChiInv[0][0]*ChiInv[2][1]*ChiN[1][0]*KrDelta[1][0] +
          ChiInv[1][0]*ChiInv[2][1]*ChiN[1][1]*KrDelta[1][0] + ChiInv[2][0]*ChiInv[2][1]*ChiN[1][2]*KrDelta[1][0] +
          ChiInv[0][0]*ChiInv[2][2]*ChiN[2][0]*KrDelta[1][0] + ChiInv[1][0]*ChiInv[2][2]*ChiN[2][1]*KrDelta[1][0] +
          ChiInv[2][0]*ChiInv[2][2]*ChiN[2][2]*KrDelta[1][0]);//*w[0][1] +
     TChi_1(2,6)=(ChiInv[0][0]*ChiInv[2][0]*ChiN[0][0]*KrDelta[2][0] + ChiInv[1][0]*ChiInv[2][0]*ChiN[0][1]*KrDelta[2][0] +
          ChiInv[2][0]*ChiInv[2][0]*ChiN[0][2]*KrDelta[2][0] + ChiInv[0][0]*ChiInv[2][1]*ChiN[1][0]*KrDelta[2][0] +
          ChiInv[1][0]*ChiInv[2][1]*ChiN[1][1]*KrDelta[2][0] + ChiInv[2][0]*ChiInv[2][1]*ChiN[1][2]*KrDelta[2][0] +
          ChiInv[0][0]*ChiInv[2][2]*ChiN[2][0]*KrDelta[2][0] + ChiInv[1][0]*ChiInv[2][2]*ChiN[2][1]*KrDelta[2][0] +
          ChiInv[2][0]*ChiInv[2][2]*ChiN[2][2]*KrDelta[2][0]);//*w[0][2] +
     TChi_1(3,6)=(ChiInv[0][0]*ChiInv[2][0]*ChiN[0][0]*KrDelta[0][1] + ChiInv[1][0]*ChiInv[2][0]*ChiN[0][1]*KrDelta[0][1] +
          ChiInv[2][0]*ChiInv[2][0]*ChiN[0][2]*KrDelta[0][1] + ChiInv[0][0]*ChiInv[2][1]*ChiN[1][0]*KrDelta[0][1] +
          ChiInv[1][0]*ChiInv[2][1]*ChiN[1][1]*KrDelta[0][1] + ChiInv[2][0]*ChiInv[2][1]*ChiN[1][2]*KrDelta[0][1] +
          ChiInv[0][0]*ChiInv[2][2]*ChiN[2][0]*KrDelta[0][1] + ChiInv[1][0]*ChiInv[2][2]*ChiN[2][1]*KrDelta[0][1] +
          ChiInv[2][0]*ChiInv[2][2]*ChiN[2][2]*KrDelta[0][1]);//*w[1][0] +
     TChi_1(4,6)=(ChiInv[0][0]*ChiInv[2][0]*ChiN[0][0]*KrDelta[1][1] + ChiInv[1][0]*ChiInv[2][0]*ChiN[0][1]*KrDelta[1][1] +
          ChiInv[2][0]*ChiInv[2][0]*ChiN[0][2]*KrDelta[1][1] + ChiInv[0][0]*ChiInv[2][1]*ChiN[1][0]*KrDelta[1][1] +
          ChiInv[1][0]*ChiInv[2][1]*ChiN[1][1]*KrDelta[1][1] + ChiInv[2][0]*ChiInv[2][1]*ChiN[1][2]*KrDelta[1][1] +
          ChiInv[0][0]*ChiInv[2][2]*ChiN[2][0]*KrDelta[1][1] + ChiInv[1][0]*ChiInv[2][2]*ChiN[2][1]*KrDelta[1][1] +
          ChiInv[2][0]*ChiInv[2][2]*ChiN[2][2]*KrDelta[1][1]);//*w[1][1] +
     TChi_1(5,6)=(ChiInv[0][0]*ChiInv[2][0]*ChiN[0][0]*KrDelta[2][1] + ChiInv[1][0]*ChiInv[2][0]*ChiN[0][1]*KrDelta[2][1] +
          ChiInv[2][0]*ChiInv[2][0]*ChiN[0][2]*KrDelta[2][1] + ChiInv[0][0]*ChiInv[2][1]*ChiN[1][0]*KrDelta[2][1] +
          ChiInv[1][0]*ChiInv[2][1]*ChiN[1][1]*KrDelta[2][1] + ChiInv[2][0]*ChiInv[2][1]*ChiN[1][2]*KrDelta[2][1] +
          ChiInv[0][0]*ChiInv[2][2]*ChiN[2][0]*KrDelta[2][1] + ChiInv[1][0]*ChiInv[2][2]*ChiN[2][1]*KrDelta[2][1] +
          ChiInv[2][0]*ChiInv[2][2]*ChiN[2][2]*KrDelta[2][1]);//*w[1][2] +
     TChi_1(6,6)=(ChiInv[0][0]*ChiInv[2][0]*ChiN[0][0]*KrDelta[0][2] + ChiInv[1][0]*ChiInv[2][0]*ChiN[0][1]*KrDelta[0][2] +
          ChiInv[2][0]*ChiInv[2][0]*ChiN[0][2]*KrDelta[0][2] + ChiInv[0][0]*ChiInv[2][1]*ChiN[1][0]*KrDelta[0][2] +
          ChiInv[1][0]*ChiInv[2][1]*ChiN[1][1]*KrDelta[0][2] + ChiInv[2][0]*ChiInv[2][1]*ChiN[1][2]*KrDelta[0][2] +
          ChiInv[0][0]*ChiInv[2][2]*ChiN[2][0]*KrDelta[0][2] + ChiInv[1][0]*ChiInv[2][2]*ChiN[2][1]*KrDelta[0][2] +
          ChiInv[2][0]*ChiInv[2][2]*ChiN[2][2]*KrDelta[0][2]);//*w[2][0] +
     TChi_1(7,6)=(ChiInv[0][0]*ChiInv[2][0]*ChiN[0][0]*KrDelta[1][2] + ChiInv[1][0]*ChiInv[2][0]*ChiN[0][1]*KrDelta[1][2] +
          ChiInv[2][0]*ChiInv[2][0]*ChiN[0][2]*KrDelta[1][2] + ChiInv[0][0]*ChiInv[2][1]*ChiN[1][0]*KrDelta[1][2] +
          ChiInv[1][0]*ChiInv[2][1]*ChiN[1][1]*KrDelta[1][2] + ChiInv[2][0]*ChiInv[2][1]*ChiN[1][2]*KrDelta[1][2] +
          ChiInv[0][0]*ChiInv[2][2]*ChiN[2][0]*KrDelta[1][2] + ChiInv[1][0]*ChiInv[2][2]*ChiN[2][1]*KrDelta[1][2] +
          ChiInv[2][0]*ChiInv[2][2]*ChiN[2][2]*KrDelta[1][2]);//*w[2][1] +
     TChi_1(8,6)=(ChiInv[0][0]*ChiInv[2][0]*ChiN[0][0]*KrDelta[2][2] + ChiInv[1][0]*ChiInv[2][0]*ChiN[0][1]*KrDelta[2][2] +
          ChiInv[2][0]*ChiInv[2][0]*ChiN[0][2]*KrDelta[2][2] + ChiInv[0][0]*ChiInv[2][1]*ChiN[1][0]*KrDelta[2][2] +
          ChiInv[1][0]*ChiInv[2][1]*ChiN[1][1]*KrDelta[2][2] + ChiInv[2][0]*ChiInv[2][1]*ChiN[1][2]*KrDelta[2][2] +
          ChiInv[0][0]*ChiInv[2][2]*ChiN[2][0]*KrDelta[2][2] + ChiInv[1][0]*ChiInv[2][2]*ChiN[2][1]*KrDelta[2][2] +
          ChiInv[2][0]*ChiInv[2][2]*ChiN[2][2]*KrDelta[2][2]);//*w[2][2]) +

     TChi_1(0,7)=(ChiInv[0][1]*ChiInv[2][0]*ChiN[0][0]*KrDelta[0][0] + ChiInv[1][1]*ChiInv[2][0]*ChiN[0][1]*KrDelta[0][0] +
          ChiInv[2][0]*ChiInv[2][1]*ChiN[0][2]*KrDelta[0][0] + ChiInv[0][1]*ChiInv[2][1]*ChiN[1][0]*KrDelta[0][0] +
          ChiInv[1][1]*ChiInv[2][1]*ChiN[1][1]*KrDelta[0][0] + ChiInv[2][1]*ChiInv[2][1]*ChiN[1][2]*KrDelta[0][0] +
          ChiInv[0][1]*ChiInv[2][2]*ChiN[2][0]*KrDelta[0][0] + ChiInv[1][1]*ChiInv[2][2]*ChiN[2][1]*KrDelta[0][0] +
          ChiInv[2][1]*ChiInv[2][2]*ChiN[2][2]*KrDelta[0][0]);//*w[0][0] +
     TChi_1(1,7)=(ChiInv[0][1]*ChiInv[2][0]*ChiN[0][0]*KrDelta[1][0] + ChiInv[1][1]*ChiInv[2][0]*ChiN[0][1]*KrDelta[1][0] +
          ChiInv[2][0]*ChiInv[2][1]*ChiN[0][2]*KrDelta[1][0] + ChiInv[0][1]*ChiInv[2][1]*ChiN[1][0]*KrDelta[1][0] +
          ChiInv[1][1]*ChiInv[2][1]*ChiN[1][1]*KrDelta[1][0] + ChiInv[2][1]*ChiInv[2][1]*ChiN[1][2]*KrDelta[1][0] +
          ChiInv[0][1]*ChiInv[2][2]*ChiN[2][0]*KrDelta[1][0] + ChiInv[1][1]*ChiInv[2][2]*ChiN[2][1]*KrDelta[1][0] +
          ChiInv[2][1]*ChiInv[2][2]*ChiN[2][2]*KrDelta[1][0]);//*w[0][1] +
     TChi_1(2,7)=(ChiInv[0][1]*ChiInv[2][0]*ChiN[0][0]*KrDelta[2][0] + ChiInv[1][1]*ChiInv[2][0]*ChiN[0][1]*KrDelta[2][0] +
          ChiInv[2][0]*ChiInv[2][1]*ChiN[0][2]*KrDelta[2][0] + ChiInv[0][1]*ChiInv[2][1]*ChiN[1][0]*KrDelta[2][0] +
          ChiInv[1][1]*ChiInv[2][1]*ChiN[1][1]*KrDelta[2][0] + ChiInv[2][1]*ChiInv[2][1]*ChiN[1][2]*KrDelta[2][0] +
          ChiInv[0][1]*ChiInv[2][2]*ChiN[2][0]*KrDelta[2][0] + ChiInv[1][1]*ChiInv[2][2]*ChiN[2][1]*KrDelta[2][0] +
          ChiInv[2][1]*ChiInv[2][2]*ChiN[2][2]*KrDelta[2][0]);//*w[0][2] +
     TChi_1(3,7)=(ChiInv[0][1]*ChiInv[2][0]*ChiN[0][0]*KrDelta[0][1] + ChiInv[1][1]*ChiInv[2][0]*ChiN[0][1]*KrDelta[0][1] +
          ChiInv[2][0]*ChiInv[2][1]*ChiN[0][2]*KrDelta[0][1] + ChiInv[0][1]*ChiInv[2][1]*ChiN[1][0]*KrDelta[0][1] +
          ChiInv[1][1]*ChiInv[2][1]*ChiN[1][1]*KrDelta[0][1] + ChiInv[2][1]*ChiInv[2][1]*ChiN[1][2]*KrDelta[0][1] +
          ChiInv[0][1]*ChiInv[2][2]*ChiN[2][0]*KrDelta[0][1] + ChiInv[1][1]*ChiInv[2][2]*ChiN[2][1]*KrDelta[0][1] +
          ChiInv[2][1]*ChiInv[2][2]*ChiN[2][2]*KrDelta[0][1]);//*w[1][0] +
     TChi_1(4,7)=(ChiInv[0][1]*ChiInv[2][0]*ChiN[0][0]*KrDelta[1][1] + ChiInv[1][1]*ChiInv[2][0]*ChiN[0][1]*KrDelta[1][1] +
          ChiInv[2][0]*ChiInv[2][1]*ChiN[0][2]*KrDelta[1][1] + ChiInv[0][1]*ChiInv[2][1]*ChiN[1][0]*KrDelta[1][1] +
          ChiInv[1][1]*ChiInv[2][1]*ChiN[1][1]*KrDelta[1][1] + ChiInv[2][1]*ChiInv[2][1]*ChiN[1][2]*KrDelta[1][1] +
          ChiInv[0][1]*ChiInv[2][2]*ChiN[2][0]*KrDelta[1][1] + ChiInv[1][1]*ChiInv[2][2]*ChiN[2][1]*KrDelta[1][1] +
          ChiInv[2][1]*ChiInv[2][2]*ChiN[2][2]*KrDelta[1][1]);//*w[1][1] +
     TChi_1(5,7)=(ChiInv[0][1]*ChiInv[2][0]*ChiN[0][0]*KrDelta[2][1] + ChiInv[1][1]*ChiInv[2][0]*ChiN[0][1]*KrDelta[2][1] +
          ChiInv[2][0]*ChiInv[2][1]*ChiN[0][2]*KrDelta[2][1] + ChiInv[0][1]*ChiInv[2][1]*ChiN[1][0]*KrDelta[2][1] +
          ChiInv[1][1]*ChiInv[2][1]*ChiN[1][1]*KrDelta[2][1] + ChiInv[2][1]*ChiInv[2][1]*ChiN[1][2]*KrDelta[2][1] +
          ChiInv[0][1]*ChiInv[2][2]*ChiN[2][0]*KrDelta[2][1] + ChiInv[1][1]*ChiInv[2][2]*ChiN[2][1]*KrDelta[2][1] +
          ChiInv[2][1]*ChiInv[2][2]*ChiN[2][2]*KrDelta[2][1]);//*w[1][2] +
     TChi_1(6,7)=(ChiInv[0][1]*ChiInv[2][0]*ChiN[0][0]*KrDelta[0][2] + ChiInv[1][1]*ChiInv[2][0]*ChiN[0][1]*KrDelta[0][2] +
          ChiInv[2][0]*ChiInv[2][1]*ChiN[0][2]*KrDelta[0][2] + ChiInv[0][1]*ChiInv[2][1]*ChiN[1][0]*KrDelta[0][2] +
          ChiInv[1][1]*ChiInv[2][1]*ChiN[1][1]*KrDelta[0][2] + ChiInv[2][1]*ChiInv[2][1]*ChiN[1][2]*KrDelta[0][2] +
          ChiInv[0][1]*ChiInv[2][2]*ChiN[2][0]*KrDelta[0][2] + ChiInv[1][1]*ChiInv[2][2]*ChiN[2][1]*KrDelta[0][2] +
          ChiInv[2][1]*ChiInv[2][2]*ChiN[2][2]*KrDelta[0][2]);//*w[2][0] +
     TChi_1(7,7)=(ChiInv[0][1]*ChiInv[2][0]*ChiN[0][0]*KrDelta[1][2] + ChiInv[1][1]*ChiInv[2][0]*ChiN[0][1]*KrDelta[1][2] +
          ChiInv[2][0]*ChiInv[2][1]*ChiN[0][2]*KrDelta[1][2] + ChiInv[0][1]*ChiInv[2][1]*ChiN[1][0]*KrDelta[1][2] +
          ChiInv[1][1]*ChiInv[2][1]*ChiN[1][1]*KrDelta[1][2] + ChiInv[2][1]*ChiInv[2][1]*ChiN[1][2]*KrDelta[1][2] +
          ChiInv[0][1]*ChiInv[2][2]*ChiN[2][0]*KrDelta[1][2] + ChiInv[1][1]*ChiInv[2][2]*ChiN[2][1]*KrDelta[1][2] +
          ChiInv[2][1]*ChiInv[2][2]*ChiN[2][2]*KrDelta[1][2]);//*w[2][1] +
     TChi_1(8,7)=(ChiInv[0][1]*ChiInv[2][0]*ChiN[0][0]*KrDelta[2][2] + ChiInv[1][1]*ChiInv[2][0]*ChiN[0][1]*KrDelta[2][2] +
          ChiInv[2][0]*ChiInv[2][1]*ChiN[0][2]*KrDelta[2][2] + ChiInv[0][1]*ChiInv[2][1]*ChiN[1][0]*KrDelta[2][2] +
          ChiInv[1][1]*ChiInv[2][1]*ChiN[1][1]*KrDelta[2][2] + ChiInv[2][1]*ChiInv[2][1]*ChiN[1][2]*KrDelta[2][2] +
          ChiInv[0][1]*ChiInv[2][2]*ChiN[2][0]*KrDelta[2][2] + ChiInv[1][1]*ChiInv[2][2]*ChiN[2][1]*KrDelta[2][2] +
          ChiInv[2][1]*ChiInv[2][2]*ChiN[2][2]*KrDelta[2][2]);//*w[2][2]) +

     TChi_1(0,8)=(ChiInv[0][2]*ChiInv[2][0]*ChiN[0][0]*KrDelta[0][0] + ChiInv[1][2]*ChiInv[2][0]*ChiN[0][1]*KrDelta[0][0] +
          ChiInv[2][0]*ChiInv[2][2]*ChiN[0][2]*KrDelta[0][0] + ChiInv[0][2]*ChiInv[2][1]*ChiN[1][0]*KrDelta[0][0] +
          ChiInv[1][2]*ChiInv[2][1]*ChiN[1][1]*KrDelta[0][0] + ChiInv[2][1]*ChiInv[2][2]*ChiN[1][2]*KrDelta[0][0] +
          ChiInv[0][2]*ChiInv[2][2]*ChiN[2][0]*KrDelta[0][0] + ChiInv[1][2]*ChiInv[2][2]*ChiN[2][1]*KrDelta[0][0] +
          ChiInv[2][2]*ChiInv[2][2]*ChiN[2][2]*KrDelta[0][0]);//*w[0][0] +
     TChi_1(1,8)=(ChiInv[0][2]*ChiInv[2][0]*ChiN[0][0]*KrDelta[1][0] +
          ChiInv[1][2]*ChiInv[2][0]*ChiN[0][1]*KrDelta[1][0] + ChiInv[2][0]*ChiInv[2][2]*ChiN[0][2]*KrDelta[1][0] +
          ChiInv[0][2]*ChiInv[2][1]*ChiN[1][0]*KrDelta[1][0] + ChiInv[1][2]*ChiInv[2][1]*ChiN[1][1]*KrDelta[1][0] +
          ChiInv[2][1]*ChiInv[2][2]*ChiN[1][2]*KrDelta[1][0] + ChiInv[0][2]*ChiInv[2][2]*ChiN[2][0]*KrDelta[1][0] +
          ChiInv[1][2]*ChiInv[2][2]*ChiN[2][1]*KrDelta[1][0] + ChiInv[2][2]*ChiInv[2][2]*ChiN[2][2]*KrDelta[1][0]);//*w[0][1] +
     TChi_1(2,8)=(ChiInv[0][2]*ChiInv[2][0]*ChiN[0][0]*KrDelta[2][0] + ChiInv[1][2]*ChiInv[2][0]*ChiN[0][1]*KrDelta[2][0] +
          ChiInv[2][0]*ChiInv[2][2]*ChiN[0][2]*KrDelta[2][0] + ChiInv[0][2]*ChiInv[2][1]*ChiN[1][0]*KrDelta[2][0] +
          ChiInv[1][2]*ChiInv[2][1]*ChiN[1][1]*KrDelta[2][0] + ChiInv[2][1]*ChiInv[2][2]*ChiN[1][2]*KrDelta[2][0] +
          ChiInv[0][2]*ChiInv[2][2]*ChiN[2][0]*KrDelta[2][0] + ChiInv[1][2]*ChiInv[2][2]*ChiN[2][1]*KrDelta[2][0] +
          ChiInv[2][2]*ChiInv[2][2]*ChiN[2][2]*KrDelta[2][0]);//*w[0][2] +
     TChi_1(3,8)=(ChiInv[0][2]*ChiInv[2][0]*ChiN[0][0]*KrDelta[0][1] +
          ChiInv[1][2]*ChiInv[2][0]*ChiN[0][1]*KrDelta[0][1] + ChiInv[2][0]*ChiInv[2][2]*ChiN[0][2]*KrDelta[0][1] +
          ChiInv[0][2]*ChiInv[2][1]*ChiN[1][0]*KrDelta[0][1] + ChiInv[1][2]*ChiInv[2][1]*ChiN[1][1]*KrDelta[0][1] +
          ChiInv[2][1]*ChiInv[2][2]*ChiN[1][2]*KrDelta[0][1] + ChiInv[0][2]*ChiInv[2][2]*ChiN[2][0]*KrDelta[0][1] +
          ChiInv[1][2]*ChiInv[2][2]*ChiN[2][1]*KrDelta[0][1] + ChiInv[2][2]*ChiInv[2][2]*ChiN[2][2]*KrDelta[0][1]);//*w[1][0] +
     TChi_1(4,8)=(ChiInv[0][2]*ChiInv[2][0]*ChiN[0][0]*KrDelta[1][1] + ChiInv[1][2]*ChiInv[2][0]*ChiN[0][1]*KrDelta[1][1] +
          ChiInv[2][0]*ChiInv[2][2]*ChiN[0][2]*KrDelta[1][1] + ChiInv[0][2]*ChiInv[2][1]*ChiN[1][0]*KrDelta[1][1] +
          ChiInv[1][2]*ChiInv[2][1]*ChiN[1][1]*KrDelta[1][1] + ChiInv[2][1]*ChiInv[2][2]*ChiN[1][2]*KrDelta[1][1] +
          ChiInv[0][2]*ChiInv[2][2]*ChiN[2][0]*KrDelta[1][1] + ChiInv[1][2]*ChiInv[2][2]*ChiN[2][1]*KrDelta[1][1] +
          ChiInv[2][2]*ChiInv[2][2]*ChiN[2][2]*KrDelta[1][1]);//*w[1][1] +
     TChi_1(5,8)=(ChiInv[0][2]*ChiInv[2][0]*ChiN[0][0]*KrDelta[2][1] +
          ChiInv[1][2]*ChiInv[2][0]*ChiN[0][1]*KrDelta[2][1] + ChiInv[2][0]*ChiInv[2][2]*ChiN[0][2]*KrDelta[2][1] +
          ChiInv[0][2]*ChiInv[2][1]*ChiN[1][0]*KrDelta[2][1] + ChiInv[1][2]*ChiInv[2][1]*ChiN[1][1]*KrDelta[2][1] +
          ChiInv[2][1]*ChiInv[2][2]*ChiN[1][2]*KrDelta[2][1] + ChiInv[0][2]*ChiInv[2][2]*ChiN[2][0]*KrDelta[2][1] +
          ChiInv[1][2]*ChiInv[2][2]*ChiN[2][1]*KrDelta[2][1] + ChiInv[2][2]*ChiInv[2][2]*ChiN[2][2]*KrDelta[2][1]);//*w[1][2] +
     TChi_1(6,8)=(ChiInv[0][2]*ChiInv[2][0]*ChiN[0][0]*KrDelta[0][2] + ChiInv[1][2]*ChiInv[2][0]*ChiN[0][1]*KrDelta[0][2] +
          ChiInv[2][0]*ChiInv[2][2]*ChiN[0][2]*KrDelta[0][2] + ChiInv[0][2]*ChiInv[2][1]*ChiN[1][0]*KrDelta[0][2] +
          ChiInv[1][2]*ChiInv[2][1]*ChiN[1][1]*KrDelta[0][2] + ChiInv[2][1]*ChiInv[2][2]*ChiN[1][2]*KrDelta[0][2] +
          ChiInv[0][2]*ChiInv[2][2]*ChiN[2][0]*KrDelta[0][2] + ChiInv[1][2]*ChiInv[2][2]*ChiN[2][1]*KrDelta[0][2] +
          ChiInv[2][2]*ChiInv[2][2]*ChiN[2][2]*KrDelta[0][2]);//*w[2][0] +
     TChi_1(7,8)=(ChiInv[0][2]*ChiInv[2][0]*ChiN[0][0]*KrDelta[1][2] +
          ChiInv[1][2]*ChiInv[2][0]*ChiN[0][1]*KrDelta[1][2] + ChiInv[2][0]*ChiInv[2][2]*ChiN[0][2]*KrDelta[1][2] +
          ChiInv[0][2]*ChiInv[2][1]*ChiN[1][0]*KrDelta[1][2] + ChiInv[1][2]*ChiInv[2][1]*ChiN[1][1]*KrDelta[1][2] +
          ChiInv[2][1]*ChiInv[2][2]*ChiN[1][2]*KrDelta[1][2] + ChiInv[0][2]*ChiInv[2][2]*ChiN[2][0]*KrDelta[1][2] +
          ChiInv[1][2]*ChiInv[2][2]*ChiN[2][1]*KrDelta[1][2] + ChiInv[2][2]*ChiInv[2][2]*ChiN[2][2]*KrDelta[1][2]);//*w[2][1] +
     TChi_1(8,8)=(ChiInv[0][2]*ChiInv[2][0]*ChiN[0][0]*KrDelta[2][2] + ChiInv[1][2]*ChiInv[2][0]*ChiN[0][1]*KrDelta[2][2] +
          ChiInv[2][0]*ChiInv[2][2]*ChiN[0][2]*KrDelta[2][2] + ChiInv[0][2]*ChiInv[2][1]*ChiN[1][0]*KrDelta[2][2] +
          ChiInv[1][2]*ChiInv[2][1]*ChiN[1][1]*KrDelta[2][2] + ChiInv[2][1]*ChiInv[2][2]*ChiN[1][2]*KrDelta[2][2] +
          ChiInv[0][2]*ChiInv[2][2]*ChiN[2][0]*KrDelta[2][2] + ChiInv[1][2]*ChiInv[2][2]*ChiN[2][1]*KrDelta[2][2] +
          ChiInv[2][2]*ChiInv[2][2]*ChiN[2][2]*KrDelta[2][2]);//*w[2][2])


*/

}

void FSMicromorphic3DT::Form_TFn_4_matrix()
{
    int row;
    int col;
    row=0;
    col=0;
    TFn_4=0.0;
    for(int i=0;i<3;i++)
    {
        for(int m=0;m<3;m++)
        {
            row=0;
            for(int k=0;k<3;k++)
            {
                for(int l=0;l<3;l++)
                {
                    //summation over the same term starts here
                    for(int L=0;L<3;L++)
                    {
                        TFn_4(row,col)+=Fn[i][L]*Finv[L][m]*fIdentity_matrix(k,l);
                    }
                    row++;
                }
            }
            col++;
        }
    }

/*     TFn_4(0,0)=(Finv[0][0]*Fn[0][0]*KrDelta[0][0] + Finv[1][0]*Fn[0][1]*KrDelta[0][0] + Finv[2][0]*Fn[0][2]*KrDelta[0][0]);//*w[0][0] +
     TFn_4(1,0)=(Finv[0][0]*Fn[0][0]*KrDelta[1][0] + Finv[1][0]*Fn[0][1]*KrDelta[1][0] + Finv[2][0]*Fn[0][2]*KrDelta[1][0]);//*w[0][1] +
     TFn_4(2,0)=(Finv[0][0]*Fn[0][0]*KrDelta[2][0] + Finv[1][0]*Fn[0][1]*KrDelta[2][0] + Finv[2][0]*Fn[0][2]*KrDelta[2][0]);//*w[0][2] +
     TFn_4(3,0)=(Finv[0][0]*Fn[0][0]*KrDelta[0][1] + Finv[1][0]*Fn[0][1]*KrDelta[0][1] + Finv[2][0]*Fn[0][2]*KrDelta[0][1]);//*w[1][0] +
     TFn_4(4,0)=(Finv[0][0]*Fn[0][0]*KrDelta[1][1] + Finv[1][0]*Fn[0][1]*KrDelta[1][1] + Finv[2][0]*Fn[0][2]*KrDelta[1][1]);//*w[1][1] +
     TFn_4(5,0)=(Finv[0][0]*Fn[0][0]*KrDelta[2][1] + Finv[1][0]*Fn[0][1]*KrDelta[2][1] + Finv[2][0]*Fn[0][2]*KrDelta[2][1]);//*w[1][2] +
     TFn_4(6,0)=(Finv[0][0]*Fn[0][0]*KrDelta[0][2] + Finv[1][0]*Fn[0][1]*KrDelta[0][2] + Finv[2][0]*Fn[0][2]*KrDelta[0][2]);//*w[2][0] +
     TFn_4(7,0)=(Finv[0][0]*Fn[0][0]*KrDelta[1][2] + Finv[1][0]*Fn[0][1]*KrDelta[1][2] + Finv[2][0]*Fn[0][2]*KrDelta[1][2]);//*w[2][1] +
     TFn_4(8,0)=(Finv[0][0]*Fn[0][0]*KrDelta[2][2] + Finv[1][0]*Fn[0][1]*KrDelta[2][2] + Finv[2][0]*Fn[0][2]*KrDelta[2][2]);//*w[2][2]) +

     TFn_4(0,1)=(Finv[0][1]*Fn[0][0]*KrDelta[0][0] + Finv[1][1]*Fn[0][1]*KrDelta[0][0] + Finv[2][1]*Fn[0][2]*KrDelta[0][0]);//*w[0][0] +
     TFn_4(1,1)=(Finv[0][1]*Fn[0][0]*KrDelta[1][0] + Finv[1][1]*Fn[0][1]*KrDelta[1][0] + Finv[2][1]*Fn[0][2]*KrDelta[1][0]);//*w[0][1] +
     TFn_4(2,1)=(Finv[0][1]*Fn[0][0]*KrDelta[2][0] + Finv[1][1]*Fn[0][1]*KrDelta[2][0] + Finv[2][1]*Fn[0][2]*KrDelta[2][0]);//*w[0][2] +
     TFn_4(3,1)=(Finv[0][1]*Fn[0][0]*KrDelta[0][1] + Finv[1][1]*Fn[0][1]*KrDelta[0][1] + Finv[2][1]*Fn[0][2]*KrDelta[0][1]);//*w[1][0] +
     TFn_4(4,1)=(Finv[0][1]*Fn[0][0]*KrDelta[1][1] + Finv[1][1]*Fn[0][1]*KrDelta[1][1] + Finv[2][1]*Fn[0][2]*KrDelta[1][1]);//*w[1][1] +
     TFn_4(5,1)=(Finv[0][1]*Fn[0][0]*KrDelta[2][1] + Finv[1][1]*Fn[0][1]*KrDelta[2][1] + Finv[2][1]*Fn[0][2]*KrDelta[2][1]);//*w[1][2] +
     TFn_4(6,1)=(Finv[0][1]*Fn[0][0]*KrDelta[0][2] + Finv[1][1]*Fn[0][1]*KrDelta[0][2] + Finv[2][1]*Fn[0][2]*KrDelta[0][2]);//*w[2][0] +
     TFn_4(7,1)=(Finv[0][1]*Fn[0][0]*KrDelta[1][2] + Finv[1][1]*Fn[0][1]*KrDelta[1][2] + Finv[2][1]*Fn[0][2]*KrDelta[1][2]);//*w[2][1] +
     TFn_4(8,1)=(Finv[0][1]*Fn[0][0]*KrDelta[2][2] + Finv[1][1]*Fn[0][1]*KrDelta[2][2] + Finv[2][1]*Fn[0][2]*KrDelta[2][2]);//*w[2][2]) +

     TFn_4(0,2)=(Finv[0][2]*Fn[0][0]*KrDelta[0][0] + Finv[1][2]*Fn[0][1]*KrDelta[0][0] + Finv[2][2]*Fn[0][2]*KrDelta[0][0]);//*w[0][0] +
     TFn_4(1,2)=(Finv[0][2]*Fn[0][0]*KrDelta[1][0] + Finv[1][2]*Fn[0][1]*KrDelta[1][0] + Finv[2][2]*Fn[0][2]*KrDelta[1][0]);//*w[0][1] +
     TFn_4(2,2)=(Finv[0][2]*Fn[0][0]*KrDelta[2][0] + Finv[1][2]*Fn[0][1]*KrDelta[2][0] + Finv[2][2]*Fn[0][2]*KrDelta[2][0]);//*w[0][2] +
     TFn_4(3,2)=(Finv[0][2]*Fn[0][0]*KrDelta[0][1] + Finv[1][2]*Fn[0][1]*KrDelta[0][1] + Finv[2][2]*Fn[0][2]*KrDelta[0][1]);//*w[1][0] +
     TFn_4(4,2)=(Finv[0][2]*Fn[0][0]*KrDelta[1][1] + Finv[1][2]*Fn[0][1]*KrDelta[1][1] + Finv[2][2]*Fn[0][2]*KrDelta[1][1]);//*w[1][1] +
     TFn_4(5,2)=(Finv[0][2]*Fn[0][0]*KrDelta[2][1] + Finv[1][2]*Fn[0][1]*KrDelta[2][1] + Finv[2][2]*Fn[0][2]*KrDelta[2][1]);//*w[1][2] +
     TFn_4(6,2)=(Finv[0][2]*Fn[0][0]*KrDelta[0][2] + Finv[1][2]*Fn[0][1]*KrDelta[0][2] + Finv[2][2]*Fn[0][2]*KrDelta[0][2]);//*w[2][0] +
     TFn_4(7,2)=(Finv[0][2]*Fn[0][0]*KrDelta[1][2] + Finv[1][2]*Fn[0][1]*KrDelta[1][2] + Finv[2][2]*Fn[0][2]*KrDelta[1][2]);//*w[2][1] +
     TFn_4(8,2)=(Finv[0][2]*Fn[0][0]*KrDelta[2][2] + Finv[1][2]*Fn[0][1]*KrDelta[2][2] + Finv[2][2]*Fn[0][2]*KrDelta[2][2]);//*w[2][2]) +

     TFn_4(0,3)=(Finv[0][0]*Fn[1][0]*KrDelta[0][0] + Finv[1][0]*Fn[1][1]*KrDelta[0][0] + Finv[2][0]*Fn[1][2]*KrDelta[0][0]);//*w[0][0] +
     TFn_4(1,3)=(Finv[0][0]*Fn[1][0]*KrDelta[1][0] + Finv[1][0]*Fn[1][1]*KrDelta[1][0] + Finv[2][0]*Fn[1][2]*KrDelta[1][0]);//*w[0][1] +
     TFn_4(3,3)=(Finv[0][0]*Fn[1][0]*KrDelta[2][0] + Finv[1][0]*Fn[1][1]*KrDelta[2][0] + Finv[2][0]*Fn[1][2]*KrDelta[2][0]);//*w[0][2] +
     TFn_4(3,3)=(Finv[0][0]*Fn[1][0]*KrDelta[0][1] + Finv[1][0]*Fn[1][1]*KrDelta[0][1] + Finv[2][0]*Fn[1][2]*KrDelta[0][1]);//*w[1][0] +
     TFn_4(4,3)=(Finv[0][0]*Fn[1][0]*KrDelta[1][1] + Finv[1][0]*Fn[1][1]*KrDelta[1][1] + Finv[2][0]*Fn[1][2]*KrDelta[1][1]);//*w[1][1] +
     TFn_4(5,3)=(Finv[0][0]*Fn[1][0]*KrDelta[2][1] + Finv[1][0]*Fn[1][1]*KrDelta[2][1] + Finv[2][0]*Fn[1][2]*KrDelta[2][1]);//*w[1][2] +
     TFn_4(6,3)=(Finv[0][0]*Fn[1][0]*KrDelta[0][2] + Finv[1][0]*Fn[1][1]*KrDelta[0][2] + Finv[2][0]*Fn[1][2]*KrDelta[0][2]);//*w[2][0] +
     TFn_4(7,3)=(Finv[0][0]*Fn[1][0]*KrDelta[1][2] + Finv[1][0]*Fn[1][1]*KrDelta[1][2] + Finv[2][0]*Fn[1][2]*KrDelta[1][2]);//*w[2][1] +
     TFn_4(8,3)=(Finv[0][0]*Fn[1][0]*KrDelta[2][2] + Finv[1][0]*Fn[1][1]*KrDelta[2][2] + Finv[2][0]*Fn[1][2]*KrDelta[2][2]);//*w[2][2]) +

     TFn_4(0,4)=(Finv[0][1]*Fn[1][0]*KrDelta[0][0] + Finv[1][1]*Fn[1][1]*KrDelta[0][0] + Finv[2][1]*Fn[1][2]*KrDelta[0][0]);//*w[0][0] +
     TFn_4(1,4)=(Finv[0][1]*Fn[1][0]*KrDelta[1][0] + Finv[1][1]*Fn[1][1]*KrDelta[1][0] + Finv[2][1]*Fn[1][2]*KrDelta[1][0]);//*w[0][1] +
     TFn_4(2,4)=(Finv[0][1]*Fn[1][0]*KrDelta[2][0] + Finv[1][1]*Fn[1][1]*KrDelta[2][0] + Finv[2][1]*Fn[1][2]*KrDelta[2][0]);//*w[0][2] +
     TFn_4(3,4)=(Finv[0][1]*Fn[1][0]*KrDelta[0][1] + Finv[1][1]*Fn[1][1]*KrDelta[0][1] + Finv[2][1]*Fn[1][2]*KrDelta[0][1]);//*w[1][0] +
     TFn_4(4,4)=(Finv[0][1]*Fn[1][0]*KrDelta[1][1] + Finv[1][1]*Fn[1][1]*KrDelta[1][1] + Finv[2][1]*Fn[1][2]*KrDelta[1][1]);//*w[1][1] +
     TFn_4(5,4)=(Finv[0][1]*Fn[1][0]*KrDelta[2][1] + Finv[1][1]*Fn[1][1]*KrDelta[2][1] + Finv[2][1]*Fn[1][2]*KrDelta[2][1]);//*w[1][2] +
     TFn_4(6,4)=(Finv[0][1]*Fn[1][0]*KrDelta[0][2] + Finv[1][1]*Fn[1][1]*KrDelta[0][2] + Finv[2][1]*Fn[1][2]*KrDelta[0][2]);//*w[2][0] +
     TFn_4(7,4)=(Finv[0][1]*Fn[1][0]*KrDelta[1][2] + Finv[1][1]*Fn[1][1]*KrDelta[1][2] + Finv[2][1]*Fn[1][2]*KrDelta[1][2]);//*w[2][1] +
     TFn_4(8,4)=(Finv[0][1]*Fn[1][0]*KrDelta[2][2] + Finv[1][1]*Fn[1][1]*KrDelta[2][2] + Finv[2][1]*Fn[1][2]*KrDelta[2][2]);//*w[2][2]) +

     TFn_4(0,5)=(Finv[0][2]*Fn[1][0]*KrDelta[0][0] + Finv[1][2]*Fn[1][1]*KrDelta[0][0] + Finv[2][2]*Fn[1][2]*KrDelta[0][0]);//*w[0][0] +
     TFn_4(1,5)=(Finv[0][2]*Fn[1][0]*KrDelta[1][0] + Finv[1][2]*Fn[1][1]*KrDelta[1][0] + Finv[2][2]*Fn[1][2]*KrDelta[1][0]);//*w[0][1] +
     TFn_4(2,5)=(Finv[0][2]*Fn[1][0]*KrDelta[2][0] + Finv[1][2]*Fn[1][1]*KrDelta[2][0] + Finv[2][2]*Fn[1][2]*KrDelta[2][0]);//*w[0][2] +
     TFn_4(3,5)=(Finv[0][2]*Fn[1][0]*KrDelta[0][1] + Finv[1][2]*Fn[1][1]*KrDelta[0][1] + Finv[2][2]*Fn[1][2]*KrDelta[0][1]);//*w[1][0] +
     TFn_4(4,5)=(Finv[0][2]*Fn[1][0]*KrDelta[1][1] + Finv[1][2]*Fn[1][1]*KrDelta[1][1] + Finv[2][2]*Fn[1][2]*KrDelta[1][1]);//*w[1][1] +
     TFn_4(5,5)=(Finv[0][2]*Fn[1][0]*KrDelta[2][1] + Finv[1][2]*Fn[1][1]*KrDelta[2][1] + Finv[2][2]*Fn[1][2]*KrDelta[2][1]);//*w[1][2] +
     TFn_4(6,5)=(Finv[0][2]*Fn[1][0]*KrDelta[0][2] + Finv[1][2]*Fn[1][1]*KrDelta[0][2] + Finv[2][2]*Fn[1][2]*KrDelta[0][2]);//*w[2][0] +
     TFn_4(7,5)=(Finv[0][2]*Fn[1][0]*KrDelta[1][2] + Finv[1][2]*Fn[1][1]*KrDelta[1][2] + Finv[2][2]*Fn[1][2]*KrDelta[1][2]);//*w[2][1] +
     TFn_4(8,5)=(Finv[0][2]*Fn[1][0]*KrDelta[2][2] + Finv[1][2]*Fn[1][1]*KrDelta[2][2] + Finv[2][2]*Fn[1][2]*KrDelta[2][2]);//*w[2][2]) +

     TFn_4(0,6)=(Finv[0][0]*Fn[2][0]*KrDelta[0][0] + Finv[1][0]*Fn[2][1]*KrDelta[0][0] + Finv[2][0]*Fn[2][2]*KrDelta[0][0]);//*w[0][0] +
     TFn_4(1,6)=(Finv[0][0]*Fn[2][0]*KrDelta[1][0] + Finv[1][0]*Fn[2][1]*KrDelta[1][0] + Finv[2][0]*Fn[2][2]*KrDelta[1][0]);//*w[0][1] +
     TFn_4(2,6)=(Finv[0][0]*Fn[2][0]*KrDelta[2][0] + Finv[1][0]*Fn[2][1]*KrDelta[2][0] + Finv[2][0]*Fn[2][2]*KrDelta[2][0]);//*w[0][2] +
     TFn_4(3,6)=(Finv[0][0]*Fn[2][0]*KrDelta[0][1] + Finv[1][0]*Fn[2][1]*KrDelta[0][1] + Finv[2][0]*Fn[2][2]*KrDelta[0][1]);//*w[1][0] +
     TFn_4(4,6)=(Finv[0][0]*Fn[2][0]*KrDelta[1][1] + Finv[1][0]*Fn[2][1]*KrDelta[1][1] + Finv[2][0]*Fn[2][2]*KrDelta[1][1]);//*w[1][1] +
     TFn_4(5,6)=(Finv[0][0]*Fn[2][0]*KrDelta[2][1] + Finv[1][0]*Fn[2][1]*KrDelta[2][1] + Finv[2][0]*Fn[2][2]*KrDelta[2][1]);//*w[1][2] +
     TFn_4(6,6)=(Finv[0][0]*Fn[2][0]*KrDelta[0][2] + Finv[1][0]*Fn[2][1]*KrDelta[0][2] + Finv[2][0]*Fn[2][2]*KrDelta[0][2]);//*w[2][0] +
     TFn_4(7,6)=(Finv[0][0]*Fn[2][0]*KrDelta[1][2] + Finv[1][0]*Fn[2][1]*KrDelta[1][2] + Finv[2][0]*Fn[2][2]*KrDelta[1][2]);//*w[2][1] +
     TFn_4(8,6)=(Finv[0][0]*Fn[2][0]*KrDelta[2][2] + Finv[1][0]*Fn[2][1]*KrDelta[2][2] + Finv[2][0]*Fn[2][2]*KrDelta[2][2]);//*w[2][2]) +

     TFn_4(0,7)=(Finv[0][1]*Fn[2][0]*KrDelta[0][0] + Finv[1][1]*Fn[2][1]*KrDelta[0][0] + Finv[2][1]*Fn[2][2]*KrDelta[0][0]);//*w[0][0] +
     TFn_4(1,7)=(Finv[0][1]*Fn[2][0]*KrDelta[1][0] + Finv[1][1]*Fn[2][1]*KrDelta[1][0] + Finv[2][1]*Fn[2][2]*KrDelta[1][0]);//*w[0][1] +
     TFn_4(2,7)=(Finv[0][1]*Fn[2][0]*KrDelta[2][0] + Finv[1][1]*Fn[2][1]*KrDelta[2][0] + Finv[2][1]*Fn[2][2]*KrDelta[2][0]);//*w[0][2] +
     TFn_4(3,7)=(Finv[0][1]*Fn[2][0]*KrDelta[0][1] + Finv[1][1]*Fn[2][1]*KrDelta[0][1] + Finv[2][1]*Fn[2][2]*KrDelta[0][1]);//*w[1][0] +
     TFn_4(4,7)=(Finv[0][1]*Fn[2][0]*KrDelta[1][1] + Finv[1][1]*Fn[2][1]*KrDelta[1][1] + Finv[2][1]*Fn[2][2]*KrDelta[1][1]);//*w[1][1] +
     TFn_4(5,7)=(Finv[0][1]*Fn[2][0]*KrDelta[2][1] + Finv[1][1]*Fn[2][1]*KrDelta[2][1] + Finv[2][1]*Fn[2][2]*KrDelta[2][1]);//*w[1][2] +
     TFn_4(6,7)=(Finv[0][1]*Fn[2][0]*KrDelta[0][2] + Finv[1][1]*Fn[2][1]*KrDelta[0][2] + Finv[2][1]*Fn[2][2]*KrDelta[0][2]);//*w[2][0] +
     TFn_4(7,7)=(Finv[0][1]*Fn[2][0]*KrDelta[1][2] + Finv[1][1]*Fn[2][1]*KrDelta[1][2] + Finv[2][1]*Fn[2][2]*KrDelta[1][2]);//*w[2][1] +
     TFn_4(8,7)=(Finv[0][1]*Fn[2][0]*KrDelta[2][2] + Finv[1][1]*Fn[2][1]*KrDelta[2][2] + Finv[2][1]*Fn[2][2]*KrDelta[2][2]);//*w[2][2]) +

     TFn_4(0,8)=(Finv[0][2]*Fn[2][0]*KrDelta[0][0] + Finv[1][2]*Fn[2][1]*KrDelta[0][0] + Finv[2][2]*Fn[2][2]*KrDelta[0][0]);
     TFn_4(1,8)=(Finv[0][2]*Fn[2][0]*KrDelta[1][0] + Finv[1][2]*Fn[2][1]*KrDelta[1][0] + Finv[2][2]*Fn[2][2]*KrDelta[1][0]);
     TFn_4(2,8)=(Finv[0][2]*Fn[2][0]*KrDelta[2][0] + Finv[1][2]*Fn[2][1]*KrDelta[2][0] + Finv[2][2]*Fn[2][2]*KrDelta[2][0]);
     TFn_4(3,8)=(Finv[0][2]*Fn[2][0]*KrDelta[0][1] + Finv[1][2]*Fn[2][1]*KrDelta[0][1] + Finv[2][2]*Fn[2][2]*KrDelta[0][1]);
     TFn_4(4,8)=(Finv[0][2]*Fn[2][0]*KrDelta[1][1] + Finv[1][2]*Fn[2][1]*KrDelta[1][1] + Finv[2][2]*Fn[2][2]*KrDelta[1][1]);
     TFn_4(5,8)=(Finv[0][2]*Fn[2][0]*KrDelta[2][1] + Finv[1][2]*Fn[2][1]*KrDelta[2][1] + Finv[2][2]*Fn[2][2]*KrDelta[2][1]);
     TFn_4(6,8)=(Finv[0][2]*Fn[2][0]*KrDelta[0][2] + Finv[1][2]*Fn[2][1]*KrDelta[0][2] + Finv[2][2]*Fn[2][2]*KrDelta[0][2]);
     TFn_4(7,8)=(Finv[0][2]*Fn[2][0]*KrDelta[1][2] + Finv[1][2]*Fn[2][1]*KrDelta[1][2] + Finv[2][2]*Fn[2][2]*KrDelta[1][2]);
     TFn_4(8,8)=(Finv[0][2]*Fn[2][0]*KrDelta[2][2] + Finv[1][2]*Fn[2][1]*KrDelta[2][2] + Finv[2][2]*Fn[2][2]*KrDelta[2][2]);*/


}

void FSMicromorphic3DT::Form_TChi_2_matrix()
{
   int row,col;
   row=0;
   col=0;
   TChi_2=0.0;

   for(int T=0;T<3;T++)
   {
       for(int m=0;m<3;m++)
       {
           row=0;
           for(int k=0;k<3;k++)
           {
               for(int l=0;l<3;l++)
               {
                   //summation
                   for(int L=0;L<3;L++)
                   {
                       TChi_2(row,col)+=ChiN[l][L]*ChiInv[L][m]*ChiInv[T][k];
                   }
                   row++;
               }
           }
        col++;
       }
   }


/*  //should be mutliplied by kappa
   TChi_2(0,0)=(ChiInv[0][0]*ChiInv[0][0]*ChiN[0][0] + ChiInv[0][0]*ChiInv[1][0]*ChiN[0][1] + ChiInv[0][0]*ChiInv[2][0]*ChiN[0][2]);
   TChi_2(1,0)=(ChiInv[0][0]*ChiInv[0][0]*ChiN[1][0] + ChiInv[0][0]*ChiInv[1][0]*ChiN[1][1] + ChiInv[0][0]*ChiInv[2][0]*ChiN[1][2]);
   TChi_2(2,0)=(ChiInv[0][0]*ChiInv[0][0]*ChiN[2][0] + ChiInv[0][0]*ChiInv[1][0]*ChiN[2][1] + ChiInv[0][0]*ChiInv[2][0]*ChiN[2][2]);
   TChi_2(3,0)=(ChiInv[0][0]*ChiInv[0][1]*ChiN[0][0] + ChiInv[0][1]*ChiInv[1][0]*ChiN[0][1] + ChiInv[0][1]*ChiInv[2][0]*ChiN[0][2]);
   TChi_2(4,0)=(ChiInv[0][0]*ChiInv[0][1]*ChiN[1][0] + ChiInv[0][1]*ChiInv[1][0]*ChiN[1][1] + ChiInv[0][1]*ChiInv[2][0]*ChiN[1][2]);
   TChi_2(5,0)=(ChiInv[0][0]*ChiInv[0][1]*ChiN[2][0] + ChiInv[0][1]*ChiInv[1][0]*ChiN[2][1] + ChiInv[0][1]*ChiInv[2][0]*ChiN[2][2]);
   TChi_2(6,0)=(ChiInv[0][0]*ChiInv[0][2]*ChiN[0][0] + ChiInv[0][2]*ChiInv[1][0]*ChiN[0][1] + ChiInv[0][2]*ChiInv[2][0]*ChiN[0][2]);
   TChi_2(7,0)=(ChiInv[0][0]*ChiInv[0][2]*ChiN[1][0] + ChiInv[0][2]*ChiInv[1][0]*ChiN[1][1] + ChiInv[0][2]*ChiInv[2][0]*ChiN[1][2]);
   TChi_2(8,0)=(ChiInv[0][0]*ChiInv[0][2]*ChiN[2][0] + ChiInv[0][2]*ChiInv[1][0]*ChiN[2][1] + ChiInv[0][2]*ChiInv[2][0]*ChiN[2][2]);


    TChi_2(0,1)=(ChiInv[0][0]*ChiInv[0][1]*ChiN[0][0] + ChiInv[0][0]*ChiInv[1][1]*ChiN[0][1] + ChiInv[0][0]*ChiInv[2][1]*ChiN[0][2]);
    TChi_2(1,1)=(ChiInv[0][0]*ChiInv[0][1]*ChiN[1][0] + ChiInv[0][0]*ChiInv[1][1]*ChiN[1][1] + ChiInv[0][0]*ChiInv[2][1]*ChiN[1][2]);
    TChi_2(2,1)=(ChiInv[0][0]*ChiInv[0][1]*ChiN[2][0] + ChiInv[0][0]*ChiInv[1][1]*ChiN[2][1] + ChiInv[0][0]*ChiInv[2][1]*ChiN[2][2]);
    TChi_2(3,1)=(ChiInv[0][1]*ChiInv[0][1]*ChiN[0][0] + ChiInv[0][1]*ChiInv[1][1]*ChiN[0][1] + ChiInv[0][1]*ChiInv[2][1]*ChiN[0][2]);
    TChi_2(4,1)=(ChiInv[0][1]*ChiInv[0][1]*ChiN[1][0] + ChiInv[0][1]*ChiInv[1][1]*ChiN[1][1] + ChiInv[0][1]*ChiInv[2][1]*ChiN[1][2]);
    TChi_2(5,1)=(ChiInv[0][1]*ChiInv[0][1]*ChiN[2][0] + ChiInv[0][1]*ChiInv[1][1]*ChiN[2][1] + ChiInv[0][1]*ChiInv[2][1]*ChiN[2][2]);
    TChi_2(6,1)=(ChiInv[0][1]*ChiInv[0][2]*ChiN[0][0] + ChiInv[0][2]*ChiInv[1][1]*ChiN[0][1] + ChiInv[0][2]*ChiInv[2][1]*ChiN[0][2]);
    TChi_2(7,1)=(ChiInv[0][1]*ChiInv[0][2]*ChiN[1][0] + ChiInv[0][2]*ChiInv[1][1]*ChiN[1][1] + ChiInv[0][2]*ChiInv[2][1]*ChiN[1][2]);
    TChi_2(8,1)=(ChiInv[0][1]*ChiInv[0][2]*ChiN[2][0] + ChiInv[0][2]*ChiInv[1][1]*ChiN[2][1] + ChiInv[0][2]*ChiInv[2][1]*ChiN[2][2]);

    TChi_2(0,2)=(ChiInv[0][0]*ChiInv[0][2]*ChiN[0][0] + ChiInv[0][0]*ChiInv[1][2]*ChiN[0][1] + ChiInv[0][0]*ChiInv[2][2]*ChiN[0][2]);
    TChi_2(1,2)=(ChiInv[0][0]*ChiInv[0][2]*ChiN[1][0] + ChiInv[0][0]*ChiInv[1][2]*ChiN[1][1] + ChiInv[0][0]*ChiInv[2][2]*ChiN[1][2]);
    TChi_2(2,2)=(ChiInv[0][0]*ChiInv[0][2]*ChiN[2][0] + ChiInv[0][0]*ChiInv[1][2]*ChiN[2][1] + ChiInv[0][0]*ChiInv[2][2]*ChiN[2][2]);
    TChi_2(3,2)=(ChiInv[0][1]*ChiInv[0][2]*ChiN[0][0] + ChiInv[0][1]*ChiInv[1][2]*ChiN[0][1] + ChiInv[0][1]*ChiInv[2][2]*ChiN[0][2]);
    TChi_2(4,2)=(ChiInv[0][1]*ChiInv[0][2]*ChiN[1][0] + ChiInv[0][1]*ChiInv[1][2]*ChiN[1][1] + ChiInv[0][1]*ChiInv[2][2]*ChiN[1][2]);
    TChi_2(5,2)=(ChiInv[0][1]*ChiInv[0][2]*ChiN[2][0] + ChiInv[0][1]*ChiInv[1][2]*ChiN[2][1] + ChiInv[0][1]*ChiInv[2][2]*ChiN[2][2]);
    TChi_2(6,2)=(ChiInv[0][2]*ChiInv[0][2]*ChiN[0][0] + ChiInv[0][2]*ChiInv[1][2]*ChiN[0][1] + ChiInv[0][2]*ChiInv[2][2]*ChiN[0][2]);
    TChi_2(7,2)=(ChiInv[0][2]*ChiInv[0][2]*ChiN[1][0] + ChiInv[0][2]*ChiInv[1][2]*ChiN[1][1] + ChiInv[0][2]*ChiInv[2][2]*ChiN[1][2]);
    TChi_2(8,2)=(ChiInv[0][2]*ChiInv[0][2]*ChiN[2][0] + ChiInv[0][2]*ChiInv[1][2]*ChiN[2][1] + ChiInv[0][2]*ChiInv[2][2]*ChiN[2][2]);


    TChi_2(0,3)=(ChiInv[0][0]*ChiInv[1][0]*ChiN[0][0] + ChiInv[1][0]*ChiInv[1][0]*ChiN[0][1] + ChiInv[1][0]*ChiInv[2][0]*ChiN[0][2]);
    TChi_2(1,3)=(ChiInv[0][0]*ChiInv[1][0]*ChiN[1][0] + ChiInv[1][0]*ChiInv[1][0]*ChiN[1][1] + ChiInv[1][0]*ChiInv[2][0]*ChiN[1][2]);
    TChi_2(2,3)=(ChiInv[0][0]*ChiInv[1][0]*ChiN[2][0] + ChiInv[1][0]*ChiInv[1][0]*ChiN[2][1] + ChiInv[1][0]*ChiInv[2][0]*ChiN[2][2]);
    TChi_2(3,3)=(ChiInv[0][0]*ChiInv[1][1]*ChiN[0][0] + ChiInv[1][0]*ChiInv[1][1]*ChiN[0][1] + ChiInv[1][1]*ChiInv[2][0]*ChiN[0][2]);
    TChi_2(4,3)=(ChiInv[0][0]*ChiInv[1][1]*ChiN[1][0] + ChiInv[1][0]*ChiInv[1][1]*ChiN[1][1] + ChiInv[1][1]*ChiInv[2][0]*ChiN[1][2]);
    TChi_2(5,3)=(ChiInv[0][0]*ChiInv[1][1]*ChiN[2][0] + ChiInv[1][0]*ChiInv[1][1]*ChiN[2][1] + ChiInv[1][1]*ChiInv[2][0]*ChiN[2][2]);
    TChi_2(6,3)=(ChiInv[0][0]*ChiInv[1][2]*ChiN[0][0] + ChiInv[1][0]*ChiInv[1][2]*ChiN[0][1] + ChiInv[1][2]*ChiInv[2][0]*ChiN[0][2]);
    TChi_2(7,3)=(ChiInv[0][0]*ChiInv[1][2]*ChiN[1][0] + ChiInv[1][0]*ChiInv[1][2]*ChiN[1][1] + ChiInv[1][2]*ChiInv[2][0]*ChiN[1][2]);
    TChi_2(8,3)=(ChiInv[0][0]*ChiInv[1][2]*ChiN[2][0] + ChiInv[1][0]*ChiInv[1][2]*ChiN[2][1] + ChiInv[1][2]*ChiInv[2][0]*ChiN[2][2]);

    TChi_2(0,4)=(ChiInv[0][1]*ChiInv[1][0]*ChiN[0][0] + ChiInv[1][0]*ChiInv[1][1]*ChiN[0][1] + ChiInv[1][0]*ChiInv[2][1]*ChiN[0][2]);
    TChi_2(1,4)=(ChiInv[0][1]*ChiInv[1][0]*ChiN[1][0] + ChiInv[1][0]*ChiInv[1][1]*ChiN[1][1] + ChiInv[1][0]*ChiInv[2][1]*ChiN[1][2]);
    TChi_2(2,4)=(ChiInv[0][1]*ChiInv[1][0]*ChiN[2][0] + ChiInv[1][0]*ChiInv[1][1]*ChiN[2][1] + ChiInv[1][0]*ChiInv[2][1]*ChiN[2][2]);
    TChi_2(3,4)=(ChiInv[0][1]*ChiInv[1][1]*ChiN[0][0] + ChiInv[1][1]*ChiInv[1][1]*ChiN[0][1] + ChiInv[1][1]*ChiInv[2][1]*ChiN[0][2]);
    TChi_2(4,4)=(ChiInv[0][1]*ChiInv[1][1]*ChiN[1][0] + ChiInv[1][1]*ChiInv[1][1]*ChiN[1][1] + ChiInv[1][1]*ChiInv[2][1]*ChiN[1][2]);
    TChi_2(5,4)=(ChiInv[0][1]*ChiInv[1][1]*ChiN[2][0] + ChiInv[1][1]*ChiInv[1][1]*ChiN[2][1] + ChiInv[1][1]*ChiInv[2][1]*ChiN[2][2]);
    TChi_2(6,4)=(ChiInv[0][1]*ChiInv[1][2]*ChiN[0][0] + ChiInv[1][1]*ChiInv[1][2]*ChiN[0][1] + ChiInv[1][2]*ChiInv[2][1]*ChiN[0][2]);
    TChi_2(7,4)=(ChiInv[0][1]*ChiInv[1][2]*ChiN[1][0] + ChiInv[1][1]*ChiInv[1][2]*ChiN[1][1] + ChiInv[1][2]*ChiInv[2][1]*ChiN[1][2]);
    TChi_2(8,4)=(ChiInv[0][1]*ChiInv[1][2]*ChiN[2][0] + ChiInv[1][1]*ChiInv[1][2]*ChiN[2][1] + ChiInv[1][2]*ChiInv[2][1]*ChiN[2][2]);

    TChi_2(0,5)=(ChiInv[0][2]*ChiInv[1][0]*ChiN[0][0] + ChiInv[1][0]*ChiInv[1][2]*ChiN[0][1] + ChiInv[1][0]*ChiInv[2][2]*ChiN[0][2]);
    TChi_2(1,5)=(ChiInv[0][2]*ChiInv[1][0]*ChiN[1][0] + ChiInv[1][0]*ChiInv[1][2]*ChiN[1][1] + ChiInv[1][0]*ChiInv[2][2]*ChiN[1][2]);
    TChi_2(2,5)=(ChiInv[0][2]*ChiInv[1][0]*ChiN[2][0] + ChiInv[1][0]*ChiInv[1][2]*ChiN[2][1] + ChiInv[1][0]*ChiInv[2][2]*ChiN[2][2]);
    TChi_2(3,5)=(ChiInv[0][2]*ChiInv[1][1]*ChiN[0][0] + ChiInv[1][1]*ChiInv[1][2]*ChiN[0][1] + ChiInv[1][1]*ChiInv[2][2]*ChiN[0][2]);
    TChi_2(4,5)=(ChiInv[0][2]*ChiInv[1][1]*ChiN[1][0] + ChiInv[1][1]*ChiInv[1][2]*ChiN[1][1] + ChiInv[1][1]*ChiInv[2][2]*ChiN[1][2]);
    TChi_2(5,5)=(ChiInv[0][2]*ChiInv[1][1]*ChiN[2][0] + ChiInv[1][1]*ChiInv[1][2]*ChiN[2][1] + ChiInv[1][1]*ChiInv[2][2]*ChiN[2][2]);
    TChi_2(6,5)=(ChiInv[0][2]*ChiInv[1][2]*ChiN[0][0] + ChiInv[1][2]*ChiInv[1][2]*ChiN[0][1] + ChiInv[1][2]*ChiInv[2][2]*ChiN[0][2]);
    TChi_2(7,5)=(ChiInv[0][2]*ChiInv[1][2]*ChiN[1][0] + ChiInv[1][2]*ChiInv[1][2]*ChiN[1][1] + ChiInv[1][2]*ChiInv[2][2]*ChiN[1][2]);
    TChi_2(8,5)=(ChiInv[0][2]*ChiInv[1][2]*ChiN[2][0] + ChiInv[1][2]*ChiInv[1][2]*ChiN[2][1] + ChiInv[1][2]*ChiInv[2][2]*ChiN[2][2]);

    TChi_2(0,6)=(ChiInv[0][0]*ChiInv[2][0]*ChiN[0][0] + ChiInv[1][0]*ChiInv[2][0]*ChiN[0][1] + ChiInv[2][0]*ChiInv[2][0]*ChiN[0][2]);
    TChi_2(1,6)=(ChiInv[0][0]*ChiInv[2][0]*ChiN[1][0] + ChiInv[1][0]*ChiInv[2][0]*ChiN[1][1] + ChiInv[2][0]*ChiInv[2][0]*ChiN[1][2]);
    TChi_2(2,6)=(ChiInv[0][0]*ChiInv[2][0]*ChiN[2][0] + ChiInv[1][0]*ChiInv[2][0]*ChiN[2][1] + ChiInv[2][0]*ChiInv[2][0]*ChiN[2][2]);
    TChi_2(3,6)=(ChiInv[0][0]*ChiInv[2][1]*ChiN[0][0] + ChiInv[1][0]*ChiInv[2][1]*ChiN[0][1] + ChiInv[2][0]*ChiInv[2][1]*ChiN[0][2]);
    TChi_2(4,6)=(ChiInv[0][0]*ChiInv[2][1]*ChiN[1][0] + ChiInv[1][0]*ChiInv[2][1]*ChiN[1][1] + ChiInv[2][0]*ChiInv[2][1]*ChiN[1][2]);
    TChi_2(5,6)=(ChiInv[0][0]*ChiInv[2][1]*ChiN[2][0] + ChiInv[1][0]*ChiInv[2][1]*ChiN[2][1] + ChiInv[2][0]*ChiInv[2][1]*ChiN[2][2]);
    TChi_2(6,6)=(ChiInv[0][0]*ChiInv[2][2]*ChiN[0][0] + ChiInv[1][0]*ChiInv[2][2]*ChiN[0][1] + ChiInv[2][0]*ChiInv[2][2]*ChiN[0][2]);
    TChi_2(7,6)=(ChiInv[0][0]*ChiInv[2][2]*ChiN[1][0] + ChiInv[1][0]*ChiInv[2][2]*ChiN[1][1] + ChiInv[2][0]*ChiInv[2][2]*ChiN[1][2]);
    TChi_2(8,6)=(ChiInv[0][0]*ChiInv[2][2]*ChiN[2][0] + ChiInv[1][0]*ChiInv[2][2]*ChiN[2][1] + ChiInv[2][0]*ChiInv[2][2]*ChiN[2][2]);

    TChi_2(0,7)=(ChiInv[0][1]*ChiInv[2][0]*ChiN[0][0] + ChiInv[1][1]*ChiInv[2][0]*ChiN[0][1] + ChiInv[2][0]*ChiInv[2][1]*ChiN[0][2]);
    TChi_2(1,7)=(ChiInv[0][1]*ChiInv[2][0]*ChiN[1][0] + ChiInv[1][1]*ChiInv[2][0]*ChiN[1][1] + ChiInv[2][0]*ChiInv[2][1]*ChiN[1][2]);
    TChi_2(2,7)=(ChiInv[0][1]*ChiInv[2][0]*ChiN[2][0] + ChiInv[1][1]*ChiInv[2][0]*ChiN[2][1] + ChiInv[2][0]*ChiInv[2][1]*ChiN[2][2]);
    TChi_2(3,7)=(ChiInv[0][1]*ChiInv[2][1]*ChiN[0][0] + ChiInv[1][1]*ChiInv[2][1]*ChiN[0][1] + ChiInv[2][1]*ChiInv[2][1]*ChiN[0][2]);
    TChi_2(4,7)=(ChiInv[0][1]*ChiInv[2][1]*ChiN[1][0] + ChiInv[1][1]*ChiInv[2][1]*ChiN[1][1] + ChiInv[2][1]*ChiInv[2][1]*ChiN[1][2]);
    TChi_2(5,7)=(ChiInv[0][1]*ChiInv[2][1]*ChiN[2][0] + ChiInv[1][1]*ChiInv[2][1]*ChiN[2][1] + ChiInv[2][1]*ChiInv[2][1]*ChiN[2][2]);
    TChi_2(6,7)=(ChiInv[0][1]*ChiInv[2][2]*ChiN[0][0] + ChiInv[1][1]*ChiInv[2][2]*ChiN[0][1] + ChiInv[2][1]*ChiInv[2][2]*ChiN[0][2]);
    TChi_2(7,7)=(ChiInv[0][1]*ChiInv[2][2]*ChiN[1][0] + ChiInv[1][1]*ChiInv[2][2]*ChiN[1][1] + ChiInv[2][1]*ChiInv[2][2]*ChiN[1][2]);
    TChi_2(8,7)=(ChiInv[0][1]*ChiInv[2][2]*ChiN[2][0] + ChiInv[1][1]*ChiInv[2][2]*ChiN[2][1] + ChiInv[2][1]*ChiInv[2][2]*ChiN[2][2]);

    TChi_2(0,8)=(ChiInv[0][2]*ChiInv[2][0]*ChiN[0][0] + ChiInv[1][2]*ChiInv[2][0]*ChiN[0][1] + ChiInv[2][0]*ChiInv[2][2]*ChiN[0][2]);
    TChi_2(1,8)=(ChiInv[0][2]*ChiInv[2][0]*ChiN[1][0] + ChiInv[1][2]*ChiInv[2][0]*ChiN[1][1] + ChiInv[2][0]*ChiInv[2][2]*ChiN[1][2]);
    TChi_2(2,8)=(ChiInv[0][2]*ChiInv[2][0]*ChiN[2][0] + ChiInv[1][2]*ChiInv[2][0]*ChiN[2][1] + ChiInv[2][0]*ChiInv[2][2]*ChiN[2][2]);
    TChi_2(3,8)=(ChiInv[0][2]*ChiInv[2][1]*ChiN[0][0] + ChiInv[1][2]*ChiInv[2][1]*ChiN[0][1] + ChiInv[2][1]*ChiInv[2][2]*ChiN[0][2]);
    TChi_2(4,8)=(ChiInv[0][2]*ChiInv[2][1]*ChiN[1][0] + ChiInv[1][2]*ChiInv[2][1]*ChiN[1][1] + ChiInv[2][1]*ChiInv[2][2]*ChiN[1][2]);
    TChi_2(5,8)=(ChiInv[0][2]*ChiInv[2][1]*ChiN[2][0] + ChiInv[1][2]*ChiInv[2][1]*ChiN[2][1] + ChiInv[2][1]*ChiInv[2][2]*ChiN[2][2]);
    TChi_2(6,8)=(ChiInv[0][2]*ChiInv[2][2]*ChiN[0][0] + ChiInv[1][2]*ChiInv[2][2]*ChiN[0][1] + ChiInv[2][2]*ChiInv[2][2]*ChiN[0][2]);
    TChi_2(7,8)=(ChiInv[0][2]*ChiInv[2][2]*ChiN[1][0] + ChiInv[1][2]*ChiInv[2][2]*ChiN[1][1] + ChiInv[2][2]*ChiInv[2][2]*ChiN[1][2]);
    TChi_2(8,8)=(ChiInv[0][2]*ChiInv[2][2]*ChiN[2][0] + ChiInv[1][2]*ChiInv[2][2]*ChiN[2][1] + ChiInv[2][2]*ChiInv[2][2]*ChiN[2][2]);
*/
}

void FSMicromorphic3DT::Form_TFn_5_matrix()
{
    int row,col;
    row=0;
    col=0;
    TFn_5=0.0;
    for(int l=0;l<3;l++)
    {
        for(int m=0;m<3;m++)
        {
            //
            row=l;
            for(int k=0;k<3;k++)
            {
                //summation
                for(int L=0;L<3;L++)
                {
                    TFn_5(row,col)+=Fn[k][L]*Finv[L][m];
                }
                row=row+3;
            }
            col++;
        }
    }
/*    TFn_5=0.0;
    TFn_5(0,0)=(Finv[0][0]*Fn[0][0] + Finv[1][0]*Fn[0][1] + Finv[2][0]*Fn[0][2]);
    TFn_5(3,0)=(Finv[0][0]*Fn[1][0] + Finv[1][0]*Fn[1][1] + Finv[2][0]*Fn[1][2]);
    TFn_5(6,0)=(Finv[0][0]*Fn[2][0] + Finv[1][0]*Fn[2][1] + Finv[2][0]*Fn[2][2]);

    TFn_5(0,1)=(Finv[0][1]*Fn[0][0] + Finv[1][1]*Fn[0][1] + Finv[2][1]*Fn[0][2]);
    TFn_5(3,1)=(Finv[0][1]*Fn[1][0] + Finv[1][1]*Fn[1][1] + Finv[2][1]*Fn[1][2]);
    TFn_5(6,1)=(Finv[0][1]*Fn[2][0] + Finv[1][1]*Fn[2][1] + Finv[2][1]*Fn[2][2]);

    TFn_5(0,2)=(Finv[0][2]*Fn[0][0] + Finv[1][2]*Fn[0][1] + Finv[2][2]*Fn[0][2]);
    TFn_5(3,2)=(Finv[0][2]*Fn[1][0] + Finv[1][2]*Fn[1][1] + Finv[2][2]*Fn[1][2]);
    TFn_5(6,2)=(Finv[0][2]*Fn[2][0] + Finv[1][2]*Fn[2][1] + Finv[2][2]*Fn[2][2]);

    TFn_5(1,3)=(Finv[0][0]*Fn[0][0] + Finv[1][0]*Fn[0][1] + Finv[2][0]*Fn[0][2]);
    TFn_5(4,3)=(Finv[0][0]*Fn[1][0] + Finv[1][0]*Fn[1][1] + Finv[2][0]*Fn[1][2]);
    TFn_5(7,3)=(Finv[0][0]*Fn[2][0] + Finv[1][0]*Fn[2][1] + Finv[2][0]*Fn[2][2]);

    TFn_5(1,4)=(Finv[0][1]*Fn[0][0] + Finv[1][1]*Fn[0][1] + Finv[2][1]*Fn[0][2]);
    TFn_5(4,4)=(Finv[0][1]*Fn[1][0] + Finv[1][1]*Fn[1][1] + Finv[2][1]*Fn[1][2]);
    TFn_5(7,4)=(Finv[0][1]*Fn[2][0] + Finv[1][1]*Fn[2][1] + Finv[2][1]*Fn[2][2]);

    TFn_5(1,5)=(Finv[0][2]*Fn[0][0] + Finv[1][2]*Fn[0][1] + Finv[2][2]*Fn[0][2]);
    TFn_5(4,5)=(Finv[0][2]*Fn[1][0] + Finv[1][2]*Fn[1][1] + Finv[2][2]*Fn[1][2]);
    TFn_5(7,5)=(Finv[0][2]*Fn[2][0] + Finv[1][2]*Fn[2][1] + Finv[2][2]*Fn[2][2]);

    TFn_5(2,6)=(Finv[0][0]*Fn[0][0] + Finv[1][0]*Fn[0][1] + Finv[2][0]*Fn[0][2]);
    TFn_5(5,6)=(Finv[0][0]*Fn[1][0] + Finv[1][0]*Fn[1][1] + Finv[2][0]*Fn[1][2]);
    TFn_5(8,6)=(Finv[0][0]*Fn[2][0] + Finv[1][0]*Fn[2][1] + Finv[2][0]*Fn[2][2]);


    TFn_5(2,7)=(Finv[0][1]*Fn[0][0] + Finv[1][1]*Fn[0][1] + Finv[2][1]*Fn[0][2]);
    TFn_5(5,7)=(Finv[0][1]*Fn[1][0] + Finv[1][1]*Fn[1][1] + Finv[2][1]*Fn[1][2]);
    TFn_5(8,7)=(Finv[0][1]*Fn[2][0] + Finv[1][1]*Fn[2][1] + Finv[2][1]*Fn[2][2]);

    TFn_5(2,8)=(Finv[0][2]*Fn[0][0] + Finv[1][2]*Fn[0][1] + Finv[2][2]*Fn[0][2]);
    TFn_5(5,8)=(Finv[0][2]*Fn[1][0] + Finv[1][2]*Fn[1][1] + Finv[2][2]*Fn[1][2]);
    TFn_5(8,8)=(Finv[0][2]*Fn[2][0] + Finv[1][2]*Fn[2][1] + Finv[2][2]*Fn[2][2]);*/

}

void FSMicromorphic3DT::Form_TChi_3_matrix()
{
    int row,col;
    row=0;
    col=0;
    TChi_3=0.0;

    for(int T=0;T<3;T++)
    {
        for(int m=0;m<3;m++)
        {
            //
            row=0;
            for(int k=0;k<3;k++)
            {
                for(int l=0;l<3;l++)
                {
                    //summation over the same term
                    for(int L=0;L<3;L++)
                    {
                        TChi_3(row,col)+=ChiN[k][L]*ChiInv[L][m]*ChiInv[T][l];
                    }
                    row++;
                }
            }
            col++;
        }
    }


/*    TChi_3(0,0)=(ChiInv[0][0]*ChiInv[0][0]*ChiN[0][0] + ChiInv[0][0]*ChiInv[1][0]*ChiN[0][1] + ChiInv[0][0]*ChiInv[2][0]*ChiN[0][2]);
    TChi_3(1,0)=(ChiInv[0][0]*ChiInv[0][1]*ChiN[0][0] + ChiInv[0][1]*ChiInv[1][0]*ChiN[0][1] + ChiInv[0][1]*ChiInv[2][0]*ChiN[0][2]);
    TChi_3(2,0)=(ChiInv[0][0]*ChiInv[0][2]*ChiN[0][0] + ChiInv[0][2]*ChiInv[1][0]*ChiN[0][1] + ChiInv[0][2]*ChiInv[2][0]*ChiN[0][2]);
    TChi_3(3,0)=(ChiInv[0][0]*ChiInv[0][0]*ChiN[1][0] + ChiInv[0][0]*ChiInv[1][0]*ChiN[1][1] + ChiInv[0][0]*ChiInv[2][0]*ChiN[1][2]);
    TChi_3(4,0)=(ChiInv[0][0]*ChiInv[0][1]*ChiN[1][0] + ChiInv[0][1]*ChiInv[1][0]*ChiN[1][1] + ChiInv[0][1]*ChiInv[2][0]*ChiN[1][2]);
    TChi_3(5,0)=(ChiInv[0][0]*ChiInv[0][2]*ChiN[1][0] + ChiInv[0][2]*ChiInv[1][0]*ChiN[1][1] + ChiInv[0][2]*ChiInv[2][0]*ChiN[1][2]);
    TChi_3(6,0)=(ChiInv[0][0]*ChiInv[0][0]*ChiN[2][0] + ChiInv[0][0]*ChiInv[1][0]*ChiN[2][1] + ChiInv[0][0]*ChiInv[2][0]*ChiN[2][2]);
    TChi_3(7,0)=(ChiInv[0][0]*ChiInv[0][1]*ChiN[2][0] + ChiInv[0][1]*ChiInv[1][0]*ChiN[2][1] + ChiInv[0][1]*ChiInv[2][0]*ChiN[2][2]);
    TChi_3(8,0)=(ChiInv[0][0]*ChiInv[0][2]*ChiN[2][0] + ChiInv[0][2]*ChiInv[1][0]*ChiN[2][1] + ChiInv[0][2]*ChiInv[2][0]*ChiN[2][2]);

    TChi_3(0,1)=(ChiInv[0][0]*ChiInv[0][1]*ChiN[0][0] + ChiInv[0][0]*ChiInv[1][1]*ChiN[0][1] + ChiInv[0][0]*ChiInv[2][1]*ChiN[0][2]);
    TChi_3(1,1)=(ChiInv[0][1]*ChiInv[0][1]*ChiN[0][0] + ChiInv[0][1]*ChiInv[1][1]*ChiN[0][1] + ChiInv[0][1]*ChiInv[2][1]*ChiN[0][2]);
    TChi_3(2,1)=(ChiInv[0][1]*ChiInv[0][2]*ChiN[0][0] + ChiInv[0][2]*ChiInv[1][1]*ChiN[0][1] + ChiInv[0][2]*ChiInv[2][1]*ChiN[0][2]);
    TChi_3(3,1)=(ChiInv[0][0]*ChiInv[0][1]*ChiN[1][0] + ChiInv[0][0]*ChiInv[1][1]*ChiN[1][1] + ChiInv[0][0]*ChiInv[2][1]*ChiN[1][2]);
    TChi_3(4,1)=(ChiInv[0][1]*ChiInv[0][1]*ChiN[1][0] + ChiInv[0][1]*ChiInv[1][1]*ChiN[1][1] + ChiInv[0][1]*ChiInv[2][1]*ChiN[1][2]);
    TChi_3(5,1)=(ChiInv[0][1]*ChiInv[0][2]*ChiN[1][0] + ChiInv[0][2]*ChiInv[1][1]*ChiN[1][1] + ChiInv[0][2]*ChiInv[2][1]*ChiN[1][2]);
    TChi_3(6,1)=(ChiInv[0][0]*ChiInv[0][1]*ChiN[2][0] + ChiInv[0][0]*ChiInv[1][1]*ChiN[2][1] + ChiInv[0][0]*ChiInv[2][1]*ChiN[2][2]);
    TChi_3(7,1)=(ChiInv[0][1]*ChiInv[0][1]*ChiN[2][0] + ChiInv[0][1]*ChiInv[1][1]*ChiN[2][1] + ChiInv[0][1]*ChiInv[2][1]*ChiN[2][2]);
    TChi_3(8,1)=(ChiInv[0][1]*ChiInv[0][2]*ChiN[2][0] + ChiInv[0][2]*ChiInv[1][1]*ChiN[2][1] + ChiInv[0][2]*ChiInv[2][1]*ChiN[2][2]);

    TChi_3(0,2)=(ChiInv[0][0]*ChiInv[0][2]*ChiN[0][0] + ChiInv[0][0]*ChiInv[1][2]*ChiN[0][1] + ChiInv[0][0]*ChiInv[2][2]*ChiN[0][2]);
    TChi_3(1,2)=(ChiInv[0][1]*ChiInv[0][2]*ChiN[0][0] + ChiInv[0][1]*ChiInv[1][2]*ChiN[0][1] + ChiInv[0][1]*ChiInv[2][2]*ChiN[0][2]);
    TChi_3(2,2)=(ChiInv[0][2]*ChiInv[0][2]*ChiN[0][0] + ChiInv[0][2]*ChiInv[1][2]*ChiN[0][1] + ChiInv[0][2]*ChiInv[2][2]*ChiN[0][2]);
    TChi_3(3,2)=(ChiInv[0][0]*ChiInv[0][2]*ChiN[1][0] + ChiInv[0][0]*ChiInv[1][2]*ChiN[1][1] + ChiInv[0][0]*ChiInv[2][2]*ChiN[1][2]);
    TChi_3(4,2)=(ChiInv[0][1]*ChiInv[0][2]*ChiN[1][0] + ChiInv[0][1]*ChiInv[1][2]*ChiN[1][1] + ChiInv[0][1]*ChiInv[2][2]*ChiN[1][2]);
    TChi_3(5,2)=(ChiInv[0][2]*ChiInv[0][2]*ChiN[1][0] + ChiInv[0][2]*ChiInv[1][2]*ChiN[1][1] + ChiInv[0][2]*ChiInv[2][2]*ChiN[1][2]);
    TChi_3(6,2)=(ChiInv[0][0]*ChiInv[0][2]*ChiN[2][0] + ChiInv[0][0]*ChiInv[1][2]*ChiN[2][1] + ChiInv[0][0]*ChiInv[2][2]*ChiN[2][2]);
    TChi_3(7,2)=(ChiInv[0][1]*ChiInv[0][2]*ChiN[2][0] + ChiInv[0][1]*ChiInv[1][2]*ChiN[2][1] + ChiInv[0][1]*ChiInv[2][2]*ChiN[2][2]);
    TChi_3(8,2)=(ChiInv[0][2]*ChiInv[0][2]*ChiN[2][0] + ChiInv[0][2]*ChiInv[1][2]*ChiN[2][1] + ChiInv[0][2]*ChiInv[2][2]*ChiN[2][2]);

    TChi_3(0,3)=(ChiInv[0][0]*ChiInv[1][0]*ChiN[0][0] + ChiInv[1][0]*ChiInv[1][0]*ChiN[0][1] + ChiInv[1][0]*ChiInv[2][0]*ChiN[0][2]);
    TChi_3(1,3)=(ChiInv[0][0]*ChiInv[1][1]*ChiN[0][0] + ChiInv[1][0]*ChiInv[1][1]*ChiN[0][1] + ChiInv[1][1]*ChiInv[2][0]*ChiN[0][2]);
    TChi_3(2,3)=(ChiInv[0][0]*ChiInv[1][2]*ChiN[0][0] + ChiInv[1][0]*ChiInv[1][2]*ChiN[0][1] + ChiInv[1][2]*ChiInv[2][0]*ChiN[0][2]);
    TChi_3(3,3)=(ChiInv[0][0]*ChiInv[1][0]*ChiN[1][0] + ChiInv[1][0]*ChiInv[1][0]*ChiN[1][1] + ChiInv[1][0]*ChiInv[2][0]*ChiN[1][2]);
    TChi_3(4,3)=(ChiInv[0][0]*ChiInv[1][1]*ChiN[1][0] + ChiInv[1][0]*ChiInv[1][1]*ChiN[1][1] + ChiInv[1][1]*ChiInv[2][0]*ChiN[1][2]);
    TChi_3(5,3)=(ChiInv[0][0]*ChiInv[1][2]*ChiN[1][0] + ChiInv[1][0]*ChiInv[1][2]*ChiN[1][1] + ChiInv[1][2]*ChiInv[2][0]*ChiN[1][2]);
    TChi_3(6,3)=(ChiInv[0][0]*ChiInv[1][0]*ChiN[2][0] + ChiInv[1][0]*ChiInv[1][0]*ChiN[2][1] + ChiInv[1][0]*ChiInv[2][0]*ChiN[2][2]);
    TChi_3(7,3)=(ChiInv[0][0]*ChiInv[1][1]*ChiN[2][0] + ChiInv[1][0]*ChiInv[1][1]*ChiN[2][1] + ChiInv[1][1]*ChiInv[2][0]*ChiN[2][2]);
    TChi_3(8,3)=(ChiInv[0][0]*ChiInv[1][2]*ChiN[2][0] + ChiInv[1][0]*ChiInv[1][2]*ChiN[2][1] + ChiInv[1][2]*ChiInv[2][0]*ChiN[2][2]);

    TChi_3(0,4)=(ChiInv[0][1]*ChiInv[1][0]*ChiN[0][0] + ChiInv[1][0]*ChiInv[1][1]*ChiN[0][1] + ChiInv[1][0]*ChiInv[2][1]*ChiN[0][2]);
    TChi_3(1,4)=(ChiInv[0][1]*ChiInv[1][1]*ChiN[0][0] + ChiInv[1][1]*ChiInv[1][1]*ChiN[0][1] + ChiInv[1][1]*ChiInv[2][1]*ChiN[0][2]);
    TChi_3(2,4)=(ChiInv[0][1]*ChiInv[1][2]*ChiN[0][0] + ChiInv[1][1]*ChiInv[1][2]*ChiN[0][1] + ChiInv[1][2]*ChiInv[2][1]*ChiN[0][2]);
    TChi_3(3,4)=(ChiInv[0][1]*ChiInv[1][0]*ChiN[1][0] + ChiInv[1][0]*ChiInv[1][1]*ChiN[1][1] + ChiInv[1][0]*ChiInv[2][1]*ChiN[1][2]);
    TChi_3(4,4)=(ChiInv[0][1]*ChiInv[1][1]*ChiN[1][0] + ChiInv[1][1]*ChiInv[1][1]*ChiN[1][1] + ChiInv[1][1]*ChiInv[2][1]*ChiN[1][2]);
    TChi_3(5,4)=(ChiInv[0][1]*ChiInv[1][2]*ChiN[1][0] + ChiInv[1][1]*ChiInv[1][2]*ChiN[1][1] + ChiInv[1][2]*ChiInv[2][1]*ChiN[1][2]);
    TChi_3(6,4)=(ChiInv[0][1]*ChiInv[1][0]*ChiN[2][0] + ChiInv[1][0]*ChiInv[1][1]*ChiN[2][1] + ChiInv[1][0]*ChiInv[2][1]*ChiN[2][2]);
    TChi_3(7,4)=(ChiInv[0][1]*ChiInv[1][1]*ChiN[2][0] + ChiInv[1][1]*ChiInv[1][1]*ChiN[2][1] + ChiInv[1][1]*ChiInv[2][1]*ChiN[2][2]);
    TChi_3(8,4)=(ChiInv[0][1]*ChiInv[1][2]*ChiN[2][0] + ChiInv[1][1]*ChiInv[1][2]*ChiN[2][1] + ChiInv[1][2]*ChiInv[2][1]*ChiN[2][2]);

    TChi_3(0,5)=(ChiInv[0][2]*ChiInv[1][0]*ChiN[0][0] + ChiInv[1][0]*ChiInv[1][2]*ChiN[0][1] + ChiInv[1][0]*ChiInv[2][2]*ChiN[0][2]);
    TChi_3(1,5)=(ChiInv[0][2]*ChiInv[1][1]*ChiN[0][0] + ChiInv[1][1]*ChiInv[1][2]*ChiN[0][1] + ChiInv[1][1]*ChiInv[2][2]*ChiN[0][2]);
    TChi_3(2,5)=(ChiInv[0][2]*ChiInv[1][2]*ChiN[0][0] + ChiInv[1][2]*ChiInv[1][2]*ChiN[0][1] + ChiInv[1][2]*ChiInv[2][2]*ChiN[0][2]);
    TChi_3(3,5)=(ChiInv[0][2]*ChiInv[1][0]*ChiN[1][0] + ChiInv[1][0]*ChiInv[1][2]*ChiN[1][1] + ChiInv[1][0]*ChiInv[2][2]*ChiN[1][2]);
    TChi_3(4,5)=(ChiInv[0][2]*ChiInv[1][1]*ChiN[1][0] + ChiInv[1][1]*ChiInv[1][2]*ChiN[1][1] + ChiInv[1][1]*ChiInv[2][2]*ChiN[1][2]);
    TChi_3(5,5)=(ChiInv[0][2]*ChiInv[1][2]*ChiN[1][0] + ChiInv[1][2]*ChiInv[1][2]*ChiN[1][1] + ChiInv[1][2]*ChiInv[2][2]*ChiN[1][2]);
    TChi_3(6,5)=(ChiInv[0][2]*ChiInv[1][0]*ChiN[2][0] + ChiInv[1][0]*ChiInv[1][2]*ChiN[2][1] + ChiInv[1][0]*ChiInv[2][2]*ChiN[2][2]);
    TChi_3(7,5)=(ChiInv[0][2]*ChiInv[1][1]*ChiN[2][0] + ChiInv[1][1]*ChiInv[1][2]*ChiN[2][1] + ChiInv[1][1]*ChiInv[2][2]*ChiN[2][2]);
    TChi_3(8,5)=(ChiInv[0][2]*ChiInv[1][2]*ChiN[2][0] + ChiInv[1][2]*ChiInv[1][2]*ChiN[2][1] + ChiInv[1][2]*ChiInv[2][2]*ChiN[2][2]);

    TChi_3(0,6)=(ChiInv[0][0]*ChiInv[2][0]*ChiN[0][0] + ChiInv[1][0]*ChiInv[2][0]*ChiN[0][1] + ChiInv[2][0]*ChiInv[2][0]*ChiN[0][2]);
    TChi_3(1,6)=(ChiInv[0][0]*ChiInv[2][1]*ChiN[0][0] + ChiInv[1][0]*ChiInv[2][1]*ChiN[0][1] + ChiInv[2][0]*ChiInv[2][1]*ChiN[0][2]);
    TChi_3(2,6)=(ChiInv[0][0]*ChiInv[2][2]*ChiN[0][0] + ChiInv[1][0]*ChiInv[2][2]*ChiN[0][1] + ChiInv[2][0]*ChiInv[2][2]*ChiN[0][2]);
    TChi_3(3,6)=(ChiInv[0][0]*ChiInv[2][0]*ChiN[1][0] + ChiInv[1][0]*ChiInv[2][0]*ChiN[1][1] + ChiInv[2][0]*ChiInv[2][0]*ChiN[1][2]);
    TChi_3(4,6)=(ChiInv[0][0]*ChiInv[2][1]*ChiN[1][0] + ChiInv[1][0]*ChiInv[2][1]*ChiN[1][1] + ChiInv[2][0]*ChiInv[2][1]*ChiN[1][2]);
    TChi_3(5,6)=(ChiInv[0][0]*ChiInv[2][2]*ChiN[1][0] + ChiInv[1][0]*ChiInv[2][2]*ChiN[1][1] + ChiInv[2][0]*ChiInv[2][2]*ChiN[1][2]);
    TChi_3(6,6)=(ChiInv[0][0]*ChiInv[2][0]*ChiN[2][0] + ChiInv[1][0]*ChiInv[2][0]*ChiN[2][1] + ChiInv[2][0]*ChiInv[2][0]*ChiN[2][2]);
    TChi_3(7,6)=(ChiInv[0][0]*ChiInv[2][1]*ChiN[2][0] + ChiInv[1][0]*ChiInv[2][1]*ChiN[2][1] + ChiInv[2][0]*ChiInv[2][1]*ChiN[2][2]);
    TChi_3(8,6)=(ChiInv[0][0]*ChiInv[2][2]*ChiN[2][0] + ChiInv[1][0]*ChiInv[2][2]*ChiN[2][1] + ChiInv[2][0]*ChiInv[2][2]*ChiN[2][2]);

    TChi_3(0,7)=(ChiInv[0][1]*ChiInv[2][0]*ChiN[0][0] + ChiInv[1][1]*ChiInv[2][0]*ChiN[0][1] + ChiInv[2][0]*ChiInv[2][1]*ChiN[0][2]);
    TChi_3(1,7)=(ChiInv[0][1]*ChiInv[2][1]*ChiN[0][0] + ChiInv[1][1]*ChiInv[2][1]*ChiN[0][1] + ChiInv[2][1]*ChiInv[2][1]*ChiN[0][2]);
    TChi_3(2,7)=(ChiInv[0][1]*ChiInv[2][2]*ChiN[0][0] + ChiInv[1][1]*ChiInv[2][2]*ChiN[0][1] + ChiInv[2][1]*ChiInv[2][2]*ChiN[0][2]);
    TChi_3(3,7)=(ChiInv[0][1]*ChiInv[2][0]*ChiN[1][0] + ChiInv[1][1]*ChiInv[2][0]*ChiN[1][1] + ChiInv[2][0]*ChiInv[2][1]*ChiN[1][2]);
    TChi_3(4,7)=(ChiInv[0][1]*ChiInv[2][1]*ChiN[1][0] + ChiInv[1][1]*ChiInv[2][1]*ChiN[1][1] + ChiInv[2][1]*ChiInv[2][1]*ChiN[1][2]);
    TChi_3(5,7)=(ChiInv[0][1]*ChiInv[2][2]*ChiN[1][0] + ChiInv[1][1]*ChiInv[2][2]*ChiN[1][1] + ChiInv[2][1]*ChiInv[2][2]*ChiN[1][2]);
    TChi_3(6,7)=(ChiInv[0][1]*ChiInv[2][0]*ChiN[2][0] + ChiInv[1][1]*ChiInv[2][0]*ChiN[2][1] + ChiInv[2][0]*ChiInv[2][1]*ChiN[2][2]);
    TChi_3(7,7)=(ChiInv[0][1]*ChiInv[2][1]*ChiN[2][0] + ChiInv[1][1]*ChiInv[2][1]*ChiN[2][1] + ChiInv[2][1]*ChiInv[2][1]*ChiN[2][2]);
    TChi_3(8,7)=(ChiInv[0][1]*ChiInv[2][2]*ChiN[2][0] + ChiInv[1][1]*ChiInv[2][2]*ChiN[2][1] + ChiInv[2][1]*ChiInv[2][2]*ChiN[2][2]);

    TChi_3(0,8)=(ChiInv[0][2]*ChiInv[2][0]*ChiN[0][0] + ChiInv[1][2]*ChiInv[2][0]*ChiN[0][1] + ChiInv[2][0]*ChiInv[2][2]*ChiN[0][2]);
    TChi_3(1,8)=(ChiInv[0][2]*ChiInv[2][1]*ChiN[0][0] + ChiInv[1][2]*ChiInv[2][1]*ChiN[0][1] + ChiInv[2][1]*ChiInv[2][2]*ChiN[0][2]);
    TChi_3(2,8)=(ChiInv[0][2]*ChiInv[2][2]*ChiN[0][0] + ChiInv[1][2]*ChiInv[2][2]*ChiN[0][1] + ChiInv[2][2]*ChiInv[2][2]*ChiN[0][2]);
    TChi_3(3,8)=(ChiInv[0][2]*ChiInv[2][0]*ChiN[1][0] + ChiInv[1][2]*ChiInv[2][0]*ChiN[1][1] + ChiInv[2][0]*ChiInv[2][2]*ChiN[1][2]);
    TChi_3(4,8)=(ChiInv[0][2]*ChiInv[2][1]*ChiN[1][0] + ChiInv[1][2]*ChiInv[2][1]*ChiN[1][1] + ChiInv[2][1]*ChiInv[2][2]*ChiN[1][2]);
    TChi_3(5,8)=(ChiInv[0][2]*ChiInv[2][2]*ChiN[1][0] + ChiInv[1][2]*ChiInv[2][2]*ChiN[1][1] + ChiInv[2][2]*ChiInv[2][2]*ChiN[1][2]);
    TChi_3(6,8)=(ChiInv[0][2]*ChiInv[2][0]*ChiN[2][0] + ChiInv[1][2]*ChiInv[2][0]*ChiN[2][1] + ChiInv[2][0]*ChiInv[2][2]*ChiN[2][2]);
    TChi_3(7,8)=(ChiInv[0][2]*ChiInv[2][1]*ChiN[2][0] + ChiInv[1][2]*ChiInv[2][1]*ChiN[2][1] + ChiInv[2][1]*ChiInv[2][2]*ChiN[2][2]);
    TChi_3(8,8)=(ChiInv[0][2]*ChiInv[2][2]*ChiN[2][0] + ChiInv[1][2]*ChiInv[2][2]*ChiN[2][1] + ChiInv[2][2]*ChiInv[2][2]*ChiN[2][2]);
*/


}

void FSMicromorphic3DT::Form_TFn_6_matrix()

{
    int row,col;
    TFn_6=0.0;
    row=0;
    col=0;
    for(int k=0;k<3;k++)
    {
        for(int m=0;m<3;m++)
        {
            row=k*3;
            for(int l=0;l<3;l++)
            {
                //summation over the same term
                for(int L=0;L<3;L++)
                {
                    TFn_6(row,col)+=Fn[l][L]*Finv[L][m];
                }
                row++;
            }
            col++;
        }
    }

/*    TFn_6(0,0)=(Finv[0][0]*Fn[0][0] + Finv[1][0]*Fn[0][1] + Finv[2][0]*Fn[0][2]);
    TFn_6(3,0)=(Finv[0][0]*Fn[1][0] + Finv[1][0]*Fn[1][1] + Finv[2][0]*Fn[1][2]);
    TFn_6(6,0)=(Finv[0][0]*Fn[2][0] + Finv[1][0]*Fn[2][1] + Finv[2][0]*Fn[2][2]);

    TFn_6(0,1)=(Finv[0][1]*Fn[0][0] + Finv[1][1]*Fn[0][1] + Finv[2][1]*Fn[0][2]);
    TFn_6(3,1)=(Finv[0][1]*Fn[1][0] + Finv[1][1]*Fn[1][1] + Finv[2][1]*Fn[1][2]);
    TFn_6(6,1)=(Finv[0][1]*Fn[2][0] + Finv[1][1]*Fn[2][1] + Finv[2][1]*Fn[2][2]);

    TFn_6(0,2)=(Finv[0][2]*Fn[0][0] + Finv[1][2]*Fn[0][1] + Finv[2][2]*Fn[0][2]);
    TFn_6(3,2)=(Finv[0][2]*Fn[1][0] + Finv[1][2]*Fn[1][1] + Finv[2][2]*Fn[1][2]);
    TFn_6(6,2)=(Finv[0][2]*Fn[2][0] + Finv[1][2]*Fn[2][1] + Finv[2][2]*Fn[2][2]);

    TFn_6(1,3)=(Finv[0][0]*Fn[0][0] + Finv[1][0]*Fn[0][1] + Finv[2][0]*Fn[0][2]);
    TFn_6(4,3)=(Finv[0][0]*Fn[1][0] + Finv[1][0]*Fn[1][1] + Finv[2][0]*Fn[1][2]);
    TFn_6(7,3)=(Finv[0][0]*Fn[2][0] + Finv[1][0]*Fn[2][1] + Finv[2][0]*Fn[2][2]);

    TFn_6(1,4)=(Finv[0][1]*Fn[0][0] + Finv[1][1]*Fn[0][1] + Finv[2][1]*Fn[0][2]);
    TFn_6(4,4)=(Finv[0][1]*Fn[1][0] + Finv[1][1]*Fn[1][1] + Finv[2][1]*Fn[1][2]);
    TFn_6(7,4)=(Finv[0][1]*Fn[2][0] + Finv[1][1]*Fn[2][1] + Finv[2][1]*Fn[2][2]);

    TFn_6(1,5)=(Finv[0][2]*Fn[0][0] + Finv[1][2]*Fn[0][1] + Finv[2][2]*Fn[0][2]);
    TFn_6(4,5)=(Finv[0][2]*Fn[1][0] + Finv[1][2]*Fn[1][1] + Finv[2][2]*Fn[1][2]);
    TFn_6(7,5)=(Finv[0][2]*Fn[2][0] + Finv[1][2]*Fn[2][1] + Finv[2][2]*Fn[2][2]);

    TFn_6(2,6)=(Finv[0][0]*Fn[0][0] + Finv[1][0]*Fn[0][1] + Finv[2][0]*Fn[0][2]);
    TFn_6(5,6)=(Finv[0][0]*Fn[1][0] + Finv[1][0]*Fn[1][1] + Finv[2][0]*Fn[1][2]);
    TFn_6(8,6)=(Finv[0][0]*Fn[2][0] + Finv[1][0]*Fn[2][1] + Finv[2][0]*Fn[2][2]);

    TFn_6(2,7)=(Finv[0][1]*Fn[0][0] + Finv[1][1]*Fn[0][1] + Finv[2][1]*Fn[0][2]);
    TFn_6(5,7)=(Finv[0][1]*Fn[1][0] + Finv[1][1]*Fn[1][1] + Finv[2][1]*Fn[1][2]);
    TFn_6(8,7)=(Finv[0][1]*Fn[2][0] + Finv[1][1]*Fn[2][1] + Finv[2][1]*Fn[2][2]);

    TFn_6(2,8)=(Finv[0][2]*Fn[0][0] + Finv[1][2]*Fn[0][1] + Finv[2][2]*Fn[0][2]);
    TFn_6(5,8)=(Finv[0][2]*Fn[1][0] + Finv[1][2]*Fn[1][1] + Finv[2][2]*Fn[1][2]);
    TFn_6(8,8)=(Finv[0][2]*Fn[2][0] + Finv[1][2]*Fn[2][1] + Finv[2][2]*Fn[2][2]);*/

}

void FSMicromorphic3DT::Form_double_Finv_from_Deformation_tensor_inverse()
{
    for(int i=0;i<=2;i++)
        {
            for(int j=0;j<=2;j++)
                {
                Finv[i][j]=fDeformation_Gradient_Inverse(i,j);}}
}

void FSMicromorphic3DT::Form_SigCurr_matrix()
{

        int row=0;
        int col=0;
        SigCurr=0.0;
/*        for(int j=0;j<3;j++)
        {
            for(int n=0;n<3;n++)
            {
                //col starts
                row=0;
                for(int k=0;k<3;k++)
                {
                    for(int l=0;l<3;l++)
                    {
                        SigCurr(row,col)=Sigma(l,k)*fIdentity_matrix(j,n);
                        row++;
                    }
                }
                col++;
            }
        }
*/

        for(int N=0;N<3;N++)
           {
               for(int n=0;n<3;n++)
               {
                   //col starts
                   row=0;
                   for(int k=0;k<3;k++)
                   {
                       for(int l=0;l<3;l++)
                       {
                           SigCurr(row,col)=Sigma(l,k)*fDeformation_Gradient_Inverse(N,n);
                           row++;
                       }
                   }
                   col++;
               }
           }


}




// Forming the matrices coming from the Bal. of First Mom. of Momtm
void FSMicromorphic3DT::Form_Finv_eta_matrix()
{
    Finv_eta=0.0;
    int row,col,count;
    row=0;
    col=0;
    count=1;
    while (count<=9)
    {
        Finv_eta(row,col)    =Finv[0][0];
        Finv_eta(row,col+1)  =Finv[1][0];
        Finv_eta(row,col+2)  =Finv[2][0];
        Finv_eta(row+1,col)  =Finv[0][1];
        Finv_eta(row+1,col+1)=Finv[1][1];
        Finv_eta(row+1,col+2)=Finv[2][1];
        Finv_eta(row+2,col)  =Finv[0][2];
        Finv_eta(row+2,col+1)=Finv[1][2];
        Finv_eta(row+2,col+2)=Finv[2][2];
        row=row+3;
        col=col+3;
        count++;
    }

}

void FSMicromorphic3DT::Form_Gradient_of_micro_shape_eta_functions(const dMatrixT &fShapeMicroGrad_temp)
{

    int row=0;
    int col=0;
    int i=0;
    GRAD_NCHI=0.0;

    for(int i=0;i<9;i++)
    {
        col=i;
        for(int j=0; j<8; j++)
        {
            GRAD_NCHI(row  ,col)  =fShapeMicroGrad_temp(0,j);
            GRAD_NCHI(row+1,col)  =fShapeMicroGrad_temp(1,j);
            GRAD_NCHI(row+2,col)  =fShapeMicroGrad_temp(2,j);
            col=col+9;
        }
        row=row+3;
    }


}


void FSMicromorphic3DT:: Form_Etagrad_matrix()
{
    int row=0;
    int col=0;
    Etagrad=0.0;
    for(int k=0;k<3;k++)
    {
        for(int i=0;i<3;i++)
        {
            row=i;
            for(int m=0;m<3;m++)
            {
                for(int l=0;l<3;l++)
                {
                    Etagrad(row,col)=Mnplus1[k][l][m];
                    row=row+3;
                }
            }
            col++;
        }
    }
}

void FSMicromorphic3DT::Form_Mm_1_matrix()
{
    int row;
    int col;
    Mm_1=0.0;
    row=0;
    col=0;
    for(int i=0;i<3;i++)
    {
        for(int p=0;p<3;p++)
        {
            //
            row=0;
            for(int m=0;m<3;m++)
            {
                for(int l=0;l<3;l++)
                {
                    for(int k=0;k<3;k++)
                    {
                        //summation
                        for(int L=0;L<3;L++)
                        {
                            Mm_1(row,col)+=Fn[i][L]*Finv[L][p]*mn[k][l][m];
                        }
                        row++;
                    }
                }
            }
            col++;
        }
    }
}


void FSMicromorphic3DT::Form_Mm_2_matrix()
{
    int row,col;
    Mm_2=0.0;
    row=0;
    col=0;
    for(int i=0;i<3;i++)
    {
        for(int p=0;p<3;p++)
        {
            //
            row=0;
            for(int m=0;m<3;m++)
            {
                for(int l=0;l<3;l++)
                {
                    for(int k=0;k<3;k++)
                    {
                        //summation
                        for(int L=0;L<3;L++)
                        {
                            Mm_2(row,col)+=Fn[k][L]*Finv[L][p]*mn[i][l][m];
                        }
                        row++;
                    }
                }
            }
            col++;
        }
    }

}

void FSMicromorphic3DT::Form_Mm_3_matrix()
{
    int row,col;
    Mm_3=0.0;
    row=0;
    col=0;
    for(int i=0;i<3;i++)
    {
        for(int p=0;p<3;p++)
        {
            //
            row=0;
            for(int m=0;m<3;m++)
            {
                for(int l=0;l<3;l++)
                {
                    for(int k=0;k<3;k++)
                    {
                        //summation
                        for(int L=0;L<3;L++)
                        {
                            Mm_3(row,col)+=Fn[l][L]*Finv[L][p]*mn[k][i][m];
                        }
                        row++;
                    }
                }
            }
            col++;
        }
    }


}

void FSMicromorphic3DT:: Form_Mm_4_matrix()
{
    Mm_4=0.0;
    int col;
    int row;
    col=0;
    for(int T = 0;T<= 2;T++)
        {
        for(int n=0;n<=2;n++)
            {
            row = 0;//row calculations start here
            for(int m = 0;m <= 2; m++)
            {
                for(int l = 0; l <= 2; l++)
                {
                    for(int k = 0; k <= 2; k++)
                    {
                    //summation on the same term starts here
                        for(int L = 0; L <= 2; L++)
                            {
                            for(int i = 0; i <= 2; i++)
                                {
                                    Mm_4(row, col) +=mn[k][l][i]*ChiN[m][L]*ChiInv[L][n]*ChiInv[T][i];
                                }
                            }
                        row++;
                    }
                }
            }
            col++;
            }
        }


}

void FSMicromorphic3DT:: Form_Mm_5_matrix()
{
Mm_5=0.0;
int col;
int row;
col=0;
for(int T = 0;T<= 2;T++)
    {
    for(int n=0;n<=2;n++)
        {
        row = 0;//row calculations start here
        for(int m = 0;m <= 2; m++)
        {
            for(int l = 0; l <= 2; l++)
            {
                for(int k = 0; k <= 2; k++)
                {
                    //summation on the same term starts here
                    for(int p = 0; p <= 2; p++)
                    {
                        for(int r = 0; r <= 2; r++)
                        {
                            for(int s = 0; s <= 2; s++)
                            {
                                for(int L = 0; L <= 2; L++)
                                {
                                    for(int i = 0; i <= 2; i++)
                                    {
                                        Mm_5(row,col) +=CCof[k][l][m][p][r][s]*ChiN[p][L]*ChiInv[L][n]*
                                                        ChiInv[T][i]*GammaN[i][r][s];}}}}}
                    row++;}
                }
            }
        col++;
        }
    }

}

void FSMicromorphic3DT:: Form_Mm_6_matrix()
{
    Mm_6=0.0;
    int col;
    int row;
    col=0;
    for(int T = 0;T<= 2;T++)
        {
        for(int n=0;n<=2;n++)
            {
            row = 0;//row calculations start here
            for(int m = 0;m <= 2; m++)
            {
                for(int l = 0; l <= 2; l++)
                {
                    for(int k = 0; k <= 2; k++)
                    {
                        //summation on the same term starts here
                        for(int p = 0; p <= 2; p++)
                        {
                            for(int r = 0; r <= 2; r++)
                            {
                                for(int s = 0; s <= 2; s++)
                                {
                                    for(int L = 0; L <= 2; L++)
                                    {
                                        for(int i = 0; i <= 2; i++)
                                        {
                                            Mm_6(row,col)+= CCof[k][l][m][p][r][s]*ChiN[i][L]*ChiInv[L][n]*
                                                            ChiInv[T][r]*GammaN[p][i][s];}}}}}
                        row++;}
                    }
                }
            col++;
            }
        }


}


void FSMicromorphic3DT:: Form_Mm_7_matrix()
{
    Mm_7=0.0;
    int col;
    int row;
    col=0;
    for(int T = 0;T<= 2;T++)
        {
        for(int n=0;n<=2;n++)
            {
            row = 0;//row calculations start here
            for(int m = 0;m <= 2; m++)
            {
                for(int l = 0; l <= 2; l++)
                {
                    for(int k = 0; k <= 2; k++)
                    {
                        // summation on the same term starts
                        for(int p = 0; p <= 2; p++)
                        {
                            for(int r = 0; r <= 2; r++)
                            {
                                for(int s = 0; s <= 2; s++)
                                {
                                    for(int L = 0; L <= 2; L++)
                                    {
                                        for(int R=0; R<=2;R++)
                                        {
                                            Mm_7(row, col) +=CCof[k][l][m][p][r][s]*GRAD_ChiN[p][L][R]*Finv[R][r]*ChiInv[L][n]*
                                                            ChiInv[T][s];}}}}}
                        //summation on the same term ends
                        row++;}
                    }
                }
            col++;
            }
        }
}


void FSMicromorphic3DT:: Form_Mm_71_matrix()
{
    Mm_71=0.0;
    int row=0;
    int col=0;
    for(int A=0; A<3;A++)
    {
        for(int p=0;p<3;p++)
        {
            for(int T=0;T<3;T++)
            {
                //
                row=0;
                for(int m=0;m<3;m++)
                {
                    for(int l=0;l<3;l++)
                    {
                        for(int k=0;k<3;k++)
                        {
                            //summation starts here
                            for(int r=0;r<3;r++)
                            {
                                for(int s=0;s<3;s++)
                                {
                                                            Mm_71(row,col)+=CCof[k][l][m][p][r][s]*Finv[T][s]*ChiInv[A][r];
                                }
                            }
                            row++;
                        }
                    }
                }
                col++;
            }
        }
    }


}

void FSMicromorphic3DT:: Form_Mm_72_matrix()
{
    Mm_72=0.0;
    int row=0;
    int col=0;
    for(int s=0;s<3;s++)
    {
        for(int a=0;a<3;a++)
        {
            row=0;
            for(int m=0;m<3;m++)
            {
                for(int l=0;l<3;l++)
                {
                    for(int k=0;k<3;k++)
                    {
                        //summation
                        for(int p=0;p<3;p++)
                            {
                            for( int r=0;r<3;r++)
                            {
                            for(int A=0;A<3;A++)
                                {
                                    for(int T=0;T<3;T++)
                                    {
                                Mm_72(row,col)+=CCof[k][l][m][p][r][s]*(GRAD_Chi[p][A][T]-GRAD_ChiN[p][A][T])*Finv[T][a]*ChiInv[A][r];
                                    }
                                }
                            }
                            }
                    row++;
                    }
                }
            }
            col++;
        }
    }
}

void FSMicromorphic3DT:: Form_Mm_73_matrix()
{
    Mm_73=0.0;
    int row=0;
    int col=0;
    for(int B=0;B<3;B++)
    {
        for(int a=0;a<3;a++)
        {
            row=0;
            for(int m=0;m<3;m++)
            {
                for(int l=0;l<3;l++)
                {
                    for(int k=0;k<3;k++)
                    {
                        //summation
                        for(int A=0;A<3;A++)
                        {
                            for(int T=0;T<3;T++)
                            {
                                for(int p=0;p<3;p++)
                                {
                                    for(int r=0;r<3;r++)
                                    {
                                        for(int s=0;s<3;s++)
                                        {
                                            Mm_73(row,col)+=CCof[k][l][m][p][r][s]*(GRAD_Chi[p][A][T]-GRAD_ChiN[p][A][T])*Finv[T][s]*ChiInv[A][a]*ChiInv[B][r];

                                        }
                                    }
                                }
                            }
                        }
                        row++;
                    }
                }
            }
            col++;
        }
    }
}


void FSMicromorphic3DT:: Form_Mm_74_matrix()
{
    Mm_74=0.0;
    int row=0;
    int col=0;
    for(int A=0;A<3;A++)
    {
        for(int p=0;p<3;p++)
        {
            row=0;
            for(int m=0;m<3;m++)
            {
                for(int l=0;l<3;l++)
                {
                    for(int k=0;k<3;k++)
                    {
                        //summation
                        for(int B=0;B<3;B++)
                        {
                            for(int T=0;T<3;T++)
                            {
                                for(int r=0; r<3;r++)
                                {
                                    for(int s=0;s<3;s++)
                                    {
                                        for(int a=0;a<3;a++)
                                        {
                                            Mm_74(row,col)+=CCof[k][l][m][p][r][s]*ChiInv[A][a]*GRAD_Chi[a][B][T]*ChiInv[B][r]*Finv[T][s];
                                        }

                                    }
                                }
                            }
                        }
                        row++;
                    }
                }
            }
            col++;
        }
    }
}


void FSMicromorphic3DT:: Form_Mm_75_matrix()
{
    Mm_75=0.0;
    int row=0;
    int col=0;
    for(int L=0;L<3;L++)
    {
        for(int i=0;i<3;i++)
        {
            row=0;
            for(int m=0;m<3;m++)
            {
                for(int l=0;l<3;l++)
                {
                    for(int k=0;k<3;k++)
                    {
                        //summation
                        for(int p=0;p<3;p++)
                        {
                            for(int A=0;A<3;A++)
                            {
                                for(int r=0;r<3;r++)
                                {
                                    for(int s=0;s<3;s++)
                                    {
                                        for(int B=0;B<3;B++)
                                        {
                                            for(int a=0;a<3;a++)
                                            {
                                                for(int T=0;T<3;T++)
                                                {
                                                    Mm_75(row,col)+=CCof[k][l][m][p][r][s]*(Chi[p][A]-ChiN[p][A])*ChiInv[A][i]*ChiInv[L][a]*GRAD_Chi[a][B][T]
                                                                  *ChiInv[B][r]*Finv[T][s];
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        row++;

                    }
                }
            }
            col++;
        }
    }
}


void FSMicromorphic3DT:: Form_Mm_76_matrix()
{
    Mm_76=0.0;
    int row=0;
    int col=0;
    for(int B=0;B<3;B++)
    {
        for(int a=0;a<3;a++)
        {
            for(int T=0;T<3;T++)
            {
                row=0;
                for(int m=0;m<3;m++)
                {
                    for(int l=0;l<3;l++)
                    {
                        for(int k=0;k<3;k++)
                        {
                            //summation

                            for(int p=0;p<3;p++)
                            {
                                for(int r=0;r<3;r++)
                                {
                                    for(int s=0;s<3;s++)
                                    {
                                        for(int A=0;A<3;A++)
                                        {
                                            Mm_76(row,col)+=CCof[k][l][m][p][r][s]*(Chi[p][A]-ChiN[p][A])*ChiInv[A][a]*ChiInv[B][r]*Finv[T][s];
                                        }
                                    }
                                }
                            }
                            row++;
                        }
                    }
                }
                col++;
            }

        }
    }
}

void FSMicromorphic3DT:: Form_Mm_77_matrix()
{
    Mm_77=0.0;
    int row=0;
    int col=0;
    for(int L=0;L<3;L++)
    {
        for(int i=0;i<3;i++)
        {
            row=0;
            for(int m=0;m<3;m++)
            {
                for(int l=0;l<3;l++)
                {
                    for(int k=0;k<3;k++)
                    {
                        //summation
                        for(int p=0;p<3;p++)
                        {
                            for(int r=0;r<3;r++)
                            {
                                for(int s=0;s<3;s++)
                                {
                                    for(int A=0;A<3;A++)
                                    {
                                        for(int a=0;a<3;a++)
                                        {
                                            for(int B=0;B<3;B++)
                                            {
                                                for(int T=0;T<3;T++)
                                                {
                                                    Mm_77(row,col)+=CCof[k][l][m][p][r][s]*(Chi[p][A]-ChiN[p][A])*ChiInv[A][a]*GRAD_Chi[a][B][T]*ChiInv[B][i]
                                                            *ChiInv[L][r]*Finv[T][s];
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        row++;
                    }
                }
            }
            col++;
        }
    }
}

void FSMicromorphic3DT:: Form_Mm_78_matrix()
{
    Mm_78=0.0;
    int row=0;
    int col=0;
    for(int s=0;s<3;s++)
    {
        for(int i=0;i<3;i++)
        {
            row=0;
            for(int m=0;m<3;m++)
            {
                for(int l=0;l<3;l++)
                {
                    for(int k=0;k<3;k++)
                    {
                        //summation
                        for(int p=0;p<3;p++)
                        {
                            for(int r=0;r<3;r++)
                            {
                                for(int A=0;A<3;A++)
                                {
                                    for(int a=0;a<3;a++)
                                    {
                                        for(int B=0;B<3;B++)
                                        {
                                            for(int T=0;T<3;T++)
                                            {
                                                Mm_78(row,col)+=CCof[k][l][m][p][r][s]*(Chi[p][A]-ChiN[p][A])*ChiInv[A][a]*GRAD_Chi[a][B][T]*ChiInv[B][r]*Finv[T][i];
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        row++;
                    }
                }
            }
            col++;
        }
    }
}

void FSMicromorphic3DT:: Form_Mm_8_matrix()
{

    Mm_8=0.0;
    int col;
    int row;
    col=0;
    for(int T = 0;T<= 2;T++)
        {
        for(int n=0;n<=2;n++)
            {
            row = 0;//row calculations start here
            for(int m = 0;m <= 2; m++)
            {
                for(int l = 0; l <= 2; l++)
                {
                    for(int k = 0; k <= 2; k++)
                    {
                        //summation on the same term starts here
                        for(int p = 0; p <= 2; p++)
                        {
                            for(int r = 0; r <= 2; r++)
                            {
                                for(int s = 0; s <= 2; s++)
                                {
                                    for(int L = 0; L <= 2; L++)
                                    {
                                        for(int i = 0; i <= 2; i++)
                                        {
                                            for(int H=0;H<=2;H++)
                                            {
                                                for(int R=0;R<=2;R++)
                                                {
                                                    Mm_8(row, col) +=CCof[k][l][m][p][r][s]*ChiN[p][L]*ChiInv[L][n]*
                                                                 ChiInv[T][i]*GRAD_Chi[i][H][R]*Finv[R][r]*ChiInv[H][s];}}}}}}}
                        row++;}
                    }
                }
            col++;
            }
        }
}

void FSMicromorphic3DT:: Form_Mm_9_matrix()
{

    Mm_9=0.0;
    int col;
    int row;
    col=0;
    for(int T = 0;T<= 2;T++)
        {
            for(int n=0;n<=2;n++)
            {
                for(int K=0;K<=2;K++)
                {
                    row = 0;//row calculations start here
                    for(int m = 0;m <= 2; m++)
                    {
                        for(int l = 0; l <= 2; l++)
                        {
                            for(int k = 0; k <= 2; k++)
                            {
                                //summation on the same term starts here
                                for(int p = 0; p <= 2; p++)
                                {
                                    for(int r = 0; r <= 2; r++)
                                    {
                                        for(int s = 0; s <= 2; s++)
                                        {
                                            for(int L = 0; L <= 2; L++)
                                            {
                                                Mm_9(row, col) +=CCof[k][l][m][p][r][s]*ChiN[p][L]*ChiInv[L][n]
                                                               *Finv[K][r]*ChiInv[T][s];}}}}
                            row++;}
                        }
                    }
                    col++;
                }
            }
        }




}

void FSMicromorphic3DT:: Form_Mm_10_matrix()
{

    Mm_10=0.0;
    int col;
    int row;
    col=0;
    for(int T = 0;T<= 2;T++)
        {
            for(int n=0;n<=2;n++)
            {
                    row = 0;//row calculations start here
                    for(int m = 0;m <= 2; m++)
                    {
                        for(int l = 0; l <= 2; l++)
                        {
                            for(int k = 0; k <= 2; k++)
                            {
                                //summation on the same term starts here
                                for(int p = 0; p <= 2; p++)
                                {
                                    for(int r = 0; r <= 2; r++)
                                    {
                                        for(int s = 0; s <= 2; s++)
                                        {
                                            for(int L = 0; L <= 2; L++)
                                            {
                                                for(int a = 0; a <= 2; a++)
                                                {
                                                    for(int K = 0; K <= 2; K++)
                                                    {
                                                        for(int R = 0; R <= 2; R++)
                                                        {
                                                            Mm_10(row,col) += CCof[k][l][m][p][r][s]*ChiN[p][L]*ChiInv[L][a]*GRAD_Chi[a][K][R]
                                                                            *Finv[R][r]*ChiInv[K][n]*ChiInv[T][s];}}}}}}}
                            row++;}
                        }
                    }
                    col++;
                }
            }
}

void FSMicromorphic3DT:: Form_Mm_11_matrix()
{

    Mm_11=0.0;
    int col;
    int row;
    col=0;
    row=0;
/*    for(int k = 0;k<= 2;k++)
        {
            for(int n=0;n<=2;n++)
            {
                    row = k;//row calculations start here//attention here is different
                    for(int m = 0;m <= 2; m++)
                    {
                        for(int l = 0; l <= 2; l++)
                        {
                                //summation on the same term starts here
                                for(int i=0;i<=2;i++)
                                {
                                    for(int L=0;L<=2;L++)
                                    {

                                        Mm_11(row, col) +=Fn[i][L]*Finv[L][n]*GammaN[i][l][m];}}
                            row=row+3;
                        }
                    }
                    col++;
                }
            }
*/
    for(int p = 0;p<= 2;p++)
        {
            for(int n=0;n<=2;n++)
            {
                   //row calculations start here//attention here is different
                row=0;
                for(int m=0;m<3;m++)
                        {
                        for(int l = 0;l <= 2; l++)
                            {
                            for(int k = 0; k <= 2; k++)
                                {
                                //summation on the same term starts here

                                    for(int r=0;r<=2;r++)
                                        {
                                        for(int s=0;s<3;s++)
                                            {
                                                for(int L=0;L<3;L++)
                                                    {
                                                        for(int i=0;i<3;i++)
                                                            {
                                                                for(int n=0;n<3;n++)
                                                                    {
                                                                    Mm_11(row, col) +=CCof[k][l][m][p][r][s]*Fn[i][L]*Finv[L][n]*GammaN[i][r][s];
                                                                    }
                                                                }
                                                        }
                                                }
                                    }
                                row++;
                                }
                            }
                        }
                    col++;
            }
        }

}


void FSMicromorphic3DT:: Form_Mm_12_matrix()
{

    Mm_12=0.0;
    int col;
    int row;
    col=0;
/*    for(int m = 0;m<= 2;m++)
        {
            for(int n=0;n<=2;n++)
            {
                    row = m*9;//row calculations start here, attention here  is different too
                    for(int l = 0; l <= 2; l++)
                        {
                            for(int k = 0; k <= 2; k++)
                            {
                                //summation on the same term starts here
                                for(int i=0;i<=2;i++)
                                {
                                    for(int L=0;L<=2;L++)
                                    {
                                        Mm_12(row, col) +=GammaN[k][l][i]*Fn[i][L]*Finv[L][n];}}
                            row++;}
                        }
                    col++;
                }
            }*/

    for(int s = 0;s<= 2;s++)
        {
            for(int n=0;n<=2;n++)
            {
                    row = 0;//row calculations start here, attention here  is different too
                    for(int m = 0; m <= 2; m++)
                        {
                            for(int l = 0; l <= 2; l++)
                            {
                               for(int k=0;k<3;k++)
                               {
                                //summation on the same term starts here
                                for(int p=0;p<=2;p++)
                                {
                                    for(int r=0;r<=2;r++)
                                    {
                                      for(int i=0;i<3;i++)
                                      {
                                          for(int L=0;L<3;L++)
                                          {
                                              Mm_12(row, col) +=CCof[k][l][m][p][r][s]*GammaN[p][r][i]*Fn[i][L]*Finv[L][n];
                                          }
                                      }

                                    }
                                }
                            row++;
                               }
                            }
                        }
                    col++;
            }
        }

}

void FSMicromorphic3DT:: Form_Mm_13_matrix()
{
    Mm_13=0.0;
    int col;
    int row;
    col=0;

/*
    for(int T = 0;T<= 2;T++)
        {
            for(int n=0;n<=2;n++)
            {
                    row = 0;//row calculations start here
                    for(int m = 0;m <= 2; m++)
                    {
                        for(int l = 0; l <= 2; l++)
                        {
                            for(int k = 0; k <= 2; k++)
                            {
                                //summation on the same term starts here
                                for(int i=0;i<=2;i++)
                                {
                                    for(int L=0;L<=2;L++)
                                    {
                                        Mm_13(row, col) +=GammaN[k][i][m]*ChiN[i][L]*ChiInv[L][n]*ChiInv[T][l];}}
                            row++;}
                        }
                    }
                    col++;
                }
            }
*/




    for(int T = 0;T<= 2;T++)
        {
            for(int t=0;t<=2;t++)
            {
                    row = 0;//row calculations start here
                    for(int m = 0;m <= 2; m++)
                    {
                        for(int l = 0; l <= 2; l++)
                        {
                            for(int k = 0; k <= 2; k++)
                            {
                                //summation on the same term starts here
                                for(int p=0;p<3;p++)
                                {
                                    for(int r=0;r<3;r++)
                                    {
                                        for(int s=0;s<3;s++)
                                        {

                                            for(int i=0;i<=2;i++)
                                            {
                                                for(int L=0;L<=2;L++)
                                                {
                                                    Mm_13(row, col) +=CCof[k][l][m][p][r][s]*GammaN[p][i][s]*ChiN[i][L]*ChiInv[L][t]*ChiInv[T][r];}}

                                        }
                                    }
                                }

                            row++;}
                        }
                    }
                    col++;
                }
            }

}


void FSMicromorphic3DT:: Form_Mm_14_matrix()
{

    Mm_14=0.0;
    int col;
    int row;
    col=0;
    for(int T = 0;T<= 2;T++)
        {
        for(int n=0;n<=2;n++)
            {
            row = 0;//row calculations start here
            for(int m = 0;m <= 2; m++)
            {
                for(int l = 0; l <= 2; l++)
                {
                    for(int k = 0; k <= 2; k++)
                    {
                      Mm_14(row, col) +=Mnplus1[k][l][m]*Finv[T][n];
                      row++;
                     }
                }
            }
            col++;
            }
        }


}

void FSMicromorphic3DT:: Form_Ru_1_matrix()
{
    Ru_1=0.0;
    int col;
    int row;
    col=0;
    for(int i = 0;i<= 2;i++)
        {
            for(int k=0;k<=2;k++)
            {
                    row = 0;//row calculations start here
                    for(int m = 0; m <= 2; m++)
                        {
                            for(int l = 0; l <= 2; l++)
                            {
                                //summation on the same term starts here
                                    for(int L=0;L<=2;L++)
                                    {
                                        Ru_1(row, col) +=Fn[i][L]*Finv[L][k]*sn_sigman(m,l);
                                    }
                            row++;
                            }
                        }
                    col++;
                }
            }
}

void FSMicromorphic3DT:: Form_Ru_2_matrix()
{

    Ru_2=0.0;
    int col;
    int row;
    col=0;
    for(int i = 0;i<= 2;i++)
        {
            for(int k=0;k<=2;k++)
            {
                    row = 0;//row calculations start here
                    for(int m = 0; m <= 2; m++)
                        {
                            for(int l = 0; l <= 2; l++)
                            {
                                //summation on the same term starts here
                                    for(int L=0;L<=2;L++)
                                    {
                                        Ru_2(row, col)+=Fn[m][L]*Finv[L][k]*sn_sigman(i,l);
                                     }

                            row++;}
                        }
                    col++;
                }
            }


}

void FSMicromorphic3DT:: Form_Ru_3_matrix()
{

    Ru_3=0.0;
    int col;
    int row;
    col=0;
    for(int i = 0;i<= 2;i++)
        {
            for(int k=0;k<=2;k++)
            {
                    row = 0;//row calculations start here
                    for(int m = 0; m <= 2; m++)
                        {
                            for(int l = 0; l <= 2; l++)
                            {
                                //summation on the same term starts here
                                    for(int L=0;L<=2;L++)
                                    {
                                        Ru_3(row, col) +=sn_sigman(m,i)*Fn[l][L]*Finv[L][k];}
                            row++;}
                        }
                    col++;
                }
            }
}

void FSMicromorphic3DT:: Form_RChi_1_matrix()
{
    RChi_1=0.0;
    int col;
    int row;
    col=0;
    for(int K = 0; K<= 2; K++)
        {
            for(int p=0; p<=2; p++)
            {
                    row = 0;//row calculations start here
                    for(int m = 0; m <= 2; m++)
                        {
                            for(int l = 0; l <= 2; l++)
                            {
                                //summation on the same term starts here
                                    for(int T=0; T<=2; T++)
                                    {
                                        RChi_1(row, col) +=ChiN[l][T]*ChiInv[T][p]*ChiInv[K][m];
                                    }
                            row++;
                           }
                        }
                    col++;
                }
            }


}

void FSMicromorphic3DT::Form_Ru_4_matrix()
{
    Ru_4=0.0;
    int col;
    int row;
    col=0;
    for(int l = 0; l<= 2; l++)
        {
            for(int k=0;k<=2;k++)
            {
                row = l;//row calculations start here
                for(int m = 0; m <= 2; m++)
                    {
                            //summation on the same term starts here
                                for(int K=0;K<=2;K++)
                                {
                                    Ru_4(row, col) +=Fn[m][K]*Finv[K][k];
                                }
                                row=row+3;
                    }
                    col++;
                }
            }


}

void FSMicromorphic3DT:: Form_RChi_2_matrix()
{

    RChi_2=0.0;
    int col;
    int row;
    col=0;
    for(int K = 0; K<= 2; K++)
        {
            for(int p=0; p<=2; p++)
            {
                row = 0;//row calculations start here
                for(int m = 0; m <= 2; m++)
                    {
                        for(int l = 0; l <= 2; l++)
                        {
                            //summation on the same term starts here
                            for(int T=0; T<=2; T++)
                            {
                                RChi_2(row, col) +=ChiN[m][T]*ChiInv[T][p]*ChiInv[K][l];
                             }
                            row++;
                        }
                    }
                col++;
                }
            }
}


void FSMicromorphic3DT:: Form_Ru_5_matrix()
{

    Ru_5=0.0;
    int col;
    int row;
    col=0;
    for(int m = 0; m<= 2; m++)
        {
            for(int k=0;k<=2;k++)
            {
                row =m*3;//row calculations start here , here is also different
                for(int l = 0; l <= 2; l++)
                    {
                            //summation on the same term starts here
                                for(int K=0;K<=2;K++)
                                {
                                    Ru_5(row, col) +=Fn[l][K]*Finv[K][k];
                                }
                        row++;
                      }
                    col++;
                }
            }



}

void FSMicromorphic3DT:: Form_Ru_6_matrix()
{
    Ru_6=0.0;
    int col=0;
    int row=0;
    for(int i=0; i<3;i++)
    {
        for(int k=0;k<3;k++)
        {
            row=0;
            for(int m=0;m<3;m++)
            {
                for(int l=0;l<3;l++)
                {
                    for(int K=0;K<3;K++)
                    {
                        Ru_6(row,col)+=Fn[i][K]*Finv[K][k]*fIdentity_matrix(m,l);
                    }
                    row++;
                }
            }
            col++;
        }
    }

}
void FSMicromorphic3DT:: Form_Ru_7_matrix()
{
    Ru_7=0.0;
    int col=0;
    int row=0;

 for(int l = 0; l<= 2; l++)
     {
         for(int k=0;k<=2;k++)
         {
             row = l;//row calculations start here
             for(int m = 0; m <= 2; m++)
                 {
                         //summation on the same term starts here
                             for(int K=0;K<=2;K++)
                             {
                                 Ru_7(row, col) =(Ru_7(row, col) +Fn[m][K]*Finv[K][k]);
                             }
                             row=row+3;
                 }
                 col++;
             }
         }


}

void FSMicromorphic3DT:: Form_RChi_3_matrix()
{
    RChi_3=0.0;
    int row=0;
    int col=0;
    for(int K=0;K<3;K++)
    {
        for(int p=0;p<3;p++)
        {
            row=0;
            for(int m=0;m<3;m++)
            {
                for(int l=0;l<3;l++)
                {
                //summation
                    for(int T=0;T<3;T++)
                    {
                        for(int i=0;i<3;i++)
                        {
                            RChi_3(row,col)+=ChiN[i][T]*ChiInv[T][p]*ChiInv[K][i]*fIdentity_matrix(m,l);

                        }
                    }
                    row++;
                }
            }
            col++;
        }
    }
}

void FSMicromorphic3DT:: Form_Ru_8_matrix()
{
    Ru_8=0.0;
    int row=0;
    int col=0;
    for(int m=0;m<3;m++)
    {
        for(int k=0;k<3;k++)
        {
            row=m*3;
                for(int l=0;l<3;l++)
                {
                    //summation over the same term
                    for(int L=0;L<3;L++)
                    {
                        Ru_8(row,col)+=Fn[l][L]*Finv[L][k];
                    }
                    row++;
                }

            col++;
        }
    }

}

void FSMicromorphic3DT:: Form_Ru_9_matrix()
{
    Ru_9=0.0;
    int row=0;
    int col=0;
    for(int i=0;i<3;i++)
    {
        for(int p=0;p<3;p++)
        {
            row=0;
            for(int m=0;m<3;m++)
            {
                for(int l=0;l<3;l++)
                {
                   //sum
                    for(int K=0;K<3;K++)
                    {
                        Ru_9(row,col)+=Fn[i][K]*Finv[K][p]*fIdentity_matrix(m,l);
                    }
                    row++;
                }
            }
            col++;
        }
    }

}



void FSMicromorphic3DT:: Form_Rs_sigma_matrix()
{

    Rs_sigma=0.0;
    int col;
    int row;
    col=0;
    for(int T = 0;T<= 2;T++)
        {
            for(int n=0;n<=2;n++)
            {
                    row = 0;//row calculations start here
                    for(int m = 0; m <= 2; m++)
                        {
                            for(int l = 0; l <= 2; l++)
                            {

                                Rs_sigma(row, col)=s_sigma(m,l)*Finv[T][n];
                            row++;
                            }
                        }
                    col++;
                }
            }


}

void FSMicromorphic3DT:: Form_R_Capital_Lambda_Chi_matrix()
{


    R_Capital_Gamma_Chi=0.0;
    int col;
    int row;
    col=0;
            for(int K=0;K<=2;K++)
            {
                for( int m=0; m<=2; m++)
                {
                    row = m*3;//row calculations start here
                    for(int l = 0; l <= 2; l++)
                    {
                        //summation on the same term starts here
                        R_Capital_Gamma_Chi(row, col)=CapitalLambda(l,K);
                        row++;
                    }
                }
            col=col++;
            }


}

void FSMicromorphic3DT:: Form_CapitalLambda_matrix()
{
    for(int i=0;i<=2;i++)
    {
        for(int j=0;j<=2;j++)
        {
            CapitalLambda(i,j)=0.0;//10^8
        }
    }
}

void FSMicromorphic3DT::Form_H1_matrix()
{
    int row=0;
    double trdeltad=0.0;
    double dtgd[3][3][3];
    double grad_Nu[3][3][3];
    double Cgamma[3][3][3];
    double dgcir[3][3][3];
    H1=0.0;
/*    double DtDnu[3][3][3];
    double Dtnu[3][3];
    double DtGdot[3][3][3];
    double DTL[3][3];
    double grad_Chi[3][3][3];
    double grad_ChiN[3][3][3];
    double DtGC[3][3][3];
    double CklmprsDtGC[3][3][3];
    double trd;
    trd=0.0;*/
 /********************************************************************************
 *********************************************************************************
 *********************************************************************************
 *********************************************************************************/
    deltaEp=0.0;
    deltaNu=0.0;
    deltad=0.0;
    deltaL=0.0;
    tempSig=0.0;
    Fn_m=0.0;
    Finv_m=0.0;
    SigN_m=0.0;
    ChiInv_m=0.0;
    ChiN_m=0.0;


    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {

            Fn_m(i,j)=Fn[i][j];
            Finv_m(i,j)=Finv[i][j];
            SigN_m(i,j)=SigN[i][j];
            ChiInv_m(i,j)=ChiInv[i][j];
            ChiN_m(i,j)=ChiN[i][j];
        }
    }

    deltaL.MultAB(Fn_m,Finv_m);
    deltaL*=-1;
    deltaL+=fIdentity_matrix;

    deltaL_Tr.Transpose(deltaL);

    deltad=deltaL;
    deltad+=deltaL_Tr;
    deltad*=0.5;
    trdeltad=0.0;
    for(int i=0;i<3;i++)
        trdeltad+=deltad(i,i);

    deltaNu.MultAB(ChiN_m,ChiInv_m);
    deltaNu*=-1;
    deltaNu+=fIdentity_matrix;
    deltaEp+=deltaNu;
    deltaEp+=deltaL_Tr;
 /********************************************************************************
 *********************************************************************************
 *********************************************************************************
 *********************************************************************************/
//constructing Mnplus1
//initiliazting  the tensors
for(int m=0;m<3;m++)
{
    for(int l=0;l<3;l++)
    {
        for(int k=0;k<3;k++)
        {
            Mnplus1[k][l][m]=0.0;
            dtgd[k][l][m]=0.0;
            grad_Nu[k][l][m]=0.0;
            Cgamma[k][l][m]=0.0;
            dgcir[k][l][m]=0.0;}}}
//calculating the dChiInvdX appearing in grad_Nu(pr,s) in equation 101

for(int p=0;p<3;p++)
{
    for(int r=0;r<3;r++)
    {
        for(int s=0;s<3;s++)
        {
            //summation
            for(int L=0;L<3;L++)
            {
                for(int l=0;l<3;l++)
                {
                    for(int K=0;K<3;K++)
                    {
                        for(int T=0;T<3;T++)
                        {
                            grad_Nu[p][r][s]+=-(Chi[p][L]-ChiN[p][L])*ChiInv[L][l]*GRAD_Chi[l][K][T]*ChiInv[K][r]*fDeformation_Gradient_Inverse(T,s);
                        }
                    }
                }
            }
        }
    }
}

for(int p=0;p<3;p++)
{
    for(int r=0;r<3;r++)
    {
        for(int s=0;s<3;s++)
        {
            //summation
            for(int K=0;K<3;K++)
            {
                for(int L=0;L<3;L++)
                {
                    grad_Nu[p][r][s]+=(GRAD_Chi[p][L][K]-GRAD_ChiN[p][L][K])*fDeformation_Gradient_Inverse(K,s)*ChiInv[L][r];
                }
            }
        }
    }
}



for(int p=0;p<3;p++)
{
    for(int r=0;r<3;r++)
    {
        for(int s=0;s<3;s++)
        {
            //
             dtgd[p][r][s]+=grad_Nu[p][r][s];
            for(int i=0;i<3;i++)
            {
               dtgd[p][r][s]+=deltaNu(p,i)*GammaN[i][r][s]-deltaNu(i,r)*GammaN[p][i][s];
            }
        }
    }
}

for(int p=0;p<3;p++)
{
    for(int r=0;r<3;r++)
    {
        for(int s=0;s<3;s++)
        {        //
            dgcir[p][r][s]=dtgd[p][r][s];
            for(int i=0;i<3;i++)
            {
                dgcir[p][r][s]+=deltaL(i,p)*GammaN[i][r][s]+GammaN[p][r][i]*deltaL(i,s)+GammaN[p][i][s]*deltaNu(i,r);
            }
        }
    }
}


for(int k=0;k<3;k++)
{
    for(int l=0;l<3;l++)
    {
        for(int m=0;m<3;m++)
        {
            for(int p=0;p<3;p++)
            {
                for(int r=0;r<3;r++)
                {
                    for(int s=0;s<3;s++)
                    {
                        Cgamma[k][l][m]+=CCof[k][l][m][p][r][s]*dgcir[p][r][s];
                    }
                }
            }
        }
    }
}

for(int m=0;m<3;m++)
    {
        for(int l=0;l<3;l++)
        {
            for(int k=0;k<3;k++)
            {
                   //
                Mnplus1[k][l][m]+=(1-trdeltad)*mn[k][l][m];
                for(int p=0;p<3;p++)
                {
                    //Mnplus1[k][l][m]+=DTL[k][p]*mn[p][l][m]+mn[k][p][m]*DTL[l][p]+mn[k][l][p]*Dtnu[m][p];
                    Mnplus1[k][l][m]+=deltaL(k,p)*mn[p][l][m]+mn[k][p][m]*deltaL(l,p)+mn[k][l][p]*deltaNu(m,p);
                }
            }
        }
    }

for(int m=0;m<3;m++)
    {
        for(int l=0;l<3;l++)
        {
            for(int k=0;k<3;k++)
            {
                Mnplus1[k][l][m]+=Cgamma[k][l][m];
            }
        }
    }

    for(int m=0;m<3;m++)
    {
        for(int l=0;l<3;l++)
        {
            for(int k=0;k<3;k++)
            {
                H1[row]=Mnplus1[k][l][m];
                row++;
            }
        }
    }


}

void FSMicromorphic3DT::Form_H2_matrix()
{
    int row;
    double trdeltad=0.0;
    double trdeltaEp=0.0;
    double scale=0.0;
    row=0;
    H2=0.0;
    s_sigma=0.0;
    deltaEp=0.0;
    deltaNu=0.0;
    deltad=0.0;
    deltaL=0.0;
    tempSig=0.0;
    Fn_m=0.0;
    Finv_m=0.0;
    SigN_m=0.0;
    ChiInv_m=0.0;
    ChiN_m=0.0;
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            Fn_m(i,j)=Fn[i][j];
            Finv_m(i,j)=Finv[i][j];
            SigN_m(i,j)=SigN[i][j];
            ChiInv_m(i,j)=ChiInv[i][j];
            ChiN_m(i,j)=ChiN[i][j];
        }
    }

    deltaL.MultAB(Fn_m,Finv_m);
    deltaL*=-1;
    deltaL+=fIdentity_matrix;

    deltaL_Tr.Transpose(deltaL);

    deltad=deltaL;
    deltad+=deltaL_Tr;
    deltad*=0.5;
    trdeltad=0.0;
    for(int i=0;i<3;i++)
        trdeltad+=deltad(i,i);

    deltaNu.MultAB(ChiN_m,ChiInv_m);
    deltaNu*=-1;
    deltaNu+=fIdentity_matrix;
    deltaEp+=deltaNu;
    deltaEp+=deltaL_Tr;

    for(int i=0;i<3;i++)
        trdeltaEp+=deltaEp(i,i);

    s_sigma+=sn_sigman;

    tempSig=sn_sigman;
    tempSig*=-trdeltad;
    s_sigma+=tempSig;

    tempSig.MultAB(deltaL,sn_sigman);
    s_sigma+=tempSig;

    tempSig.MultABT(sn_sigman,deltaL);
    s_sigma+=tempSig;

    tempSig.Transpose(deltaEp);
    scale=(fMaterial_Params[kKappa]-fMaterial_Params[kKappa]);
    tempSig*=scale;
    s_sigma+=tempSig;

    tempSig=deltaEp;
    scale=(fMaterial_Params[kNu]-fMaterial_Params[kSigma_const]);
    tempSig*=scale;
    s_sigma+=tempSig;

    tempSig=fIdentity_matrix;
    scale=trdeltad*fMaterial_Params[kTau];
    tempSig*=scale;
    s_sigma+=tempSig;

    tempSig=deltad;
    scale=2*fMaterial_Params[kSigma_const];
    tempSig*=scale;
    s_sigma+=tempSig;

    tempSig=fIdentity_matrix;
    scale=trdeltaEp*(fMaterial_Params[kEta]-fMaterial_Params[kTau]);
    tempSig*=scale;
    s_sigma+=tempSig;



 //  s_sigma*=-1;//because in formulation it is sigma-s!!

    for(int m=0;m<=2;m++)
    {
        for(int l=0;l<=2;l++)
        {
            H2[row]=-s_sigma(m,l);
            row++;
        }
    }
}

void FSMicromorphic3DT:: Form_H3_matrix()
{
    int row;
    row=0;
    H3=0.0;
    for(int m=0;m<=2;m++)
    {
        for(int l=0;l<=2;l++)
        {
            H3[row]=Lambda(l,m)-Omega(l,m);
            row++;
        }
    }
}

void FSMicromorphic3DT:: Mapping_double_and_Array(const int condition)
{
    int row;
    row=0;
  //
    if(condition==1)
    {
        for(int i=0;i<=2;i++)
        {
            for(int j=0;j<=2;j++)
            {
                for(int k=0;k<=2;k++)
                {
                    mn_ar[row]=Mnplus1[i][j][k];
                    GammaN_ar[row]=Gamma[i][j][k];
                    row++;
                }
            }
        }
    }
   if(condition==-1)
    {
        for(int i=0;i<=2;i++)
        {
            for(int j=0;j<=2;j++)
            {
                SigN[i][j]=SigN_ar(i,j);
              //  row1++;
                for(int k=0;k<=2;k++)
                {
                    mn[i][j][k]=mn_ar[row];
                    GammaN[i][j][k]=GammaN_ar[row];
                    row++;
                }
            }
        }
    }

}

void FSMicromorphic3DT:: Form_deformation_tensors_arrays(const int condition) //
{
    int row,row1;
    row=0;

    if(condition==1)
    {
        for (int i=0;i<=2;i++)
        {
            for(int j=0;j<=2;j++)
            {

                Chi_ar(i,j)=Chi[i][j];

                for(int k=0;k<=2;k++)
                {
                    GRAD_Chi_ar[row]=GRAD_Chi[i][j][k];
                    row++;
                }
            }
        }
    }
    if(condition==-1)
    {
        for (int i=0;i<=2;i++)
        {
            for(int j=0;j<=2;j++)
            {
                Fn[i][j]=Fn_ar(i,j);
                FnInv[i][j]=FnInv_ar(i,j);
                ChiN[i][j]=ChiN_ar(i,j);

                for(int k=0;k<=2;k++)
                {
                    GRAD_ChiN[i][j][k]=GRAD_ChiN_ar[row];
                    row++;
                }
            }
        }
    }

}


////////////////////////////////////////////////////////////////
//////////////FINITE STRAIN MATRICES FUNCTIONS//////////////////
////////////////////////////////////////////////////////////////
void FSMicromorphic3DT::Form_ChiM()
{
    NCHI.Multx(Phi_vec,Chi_vec);
    ChiM(0,0) = Chi_vec[0]+1.0;
    ChiM(0,1) = Chi_vec[3];
    ChiM(0,2) = Chi_vec[6];
    ChiM(1,0) = Chi_vec[1];
    ChiM(1,1) = Chi_vec[4]+1.0;
    ChiM(1,2) = Chi_vec[7];
    ChiM(2,0) = Chi_vec[2];
    ChiM(2,1) = Chi_vec[5];
    ChiM(2,2) = Chi_vec[8]+1.0;

}

void FSMicromorphic3DT::Form_Second_Piola_Kirchhoff_SPK()
{
    double trLST=0.0;
    double  trcE=0.0;
    double scale=0.0;
/*    double trcE=0.0;
    double scale=0.0;*/
    SPK=0.0;
    Temp_SPK=0.0;
    KirchhoffST=0.0;

  /*  trLST=LagrangianStn(0,0)+LagrangianStn(1,1)+LagrangianStn(2,2);
    SPK.SetToScaled(trLST*fMaterial_Params[kLambda],fIdentity_matrix);
    Temp_SPK.SetToScaled(2*fMaterial_Params[kMu],LagrangianStn);
    SPK+=Temp_SPK;
    KirchhoffST.MultABCT(fDeformation_Gradient,SPK,fDeformation_Gradient);*/

    trLST= LagrangianStn(0,0)+LagrangianStn(1,1)+LagrangianStn(2,2);
    trcE = MicroStnTensor(0,0)+MicroStnTensor(1,1)+MicroStnTensor(2,2);

   Temp_SPK=fIdentity_matrix;
   scale=(fMaterial_Params[kLambda]+fMaterial_Params[kTau]);
   scale*=trLST;
   Temp_SPK*=scale;
   SPK+=Temp_SPK;

   scale=2*(fMaterial_Params[kMu]+fMaterial_Params[kSigma_const]);
   Temp_SPK=LagrangianStn;
   Temp_SPK*=scale;
   SPK+=Temp_SPK;

   scale=fMaterial_Params[kEta]*trcE;
   Temp_SPK=fIdentity_matrix;
   Temp_SPK*=scale;
   SPK+=Temp_SPK;

   Temp_SPK=MicroStnTensor;
   Temp_SPK*=fMaterial_Params[kKappa];
   SPK+=Temp_SPK;

   Temp_SPK.Transpose(MicroStnTensor);
   Temp_SPK*=fMaterial_Params[kNu];
   SPK+=Temp_SPK;


}

void FSMicromorphic3DT:: Form_I1_1()
{
    int row=0;
    int col=0;
    I1_1=0.0;
  for(int k=0;k<3;k++)
    {
        for(int n=0;n<3;n++)
        {
            row=3*n;
            for(int l=0;l<3;l++)
            {
                //summation over the dummy indices
                I1_1(row,col)+=KirchhoffST(k,l);
              // row=row+3;
                row++;
            }
            col++;
        }

    }
}


void FSMicromorphic3DT:: Form_I1_2()
{
    int row=0;
    int col=0;
    I1_2=0.0;
     for(int n=0;n<3;n++)
       {
        for(int k=0;k<3;k++)
        {
            row=3*k;
            for(int l=0;l<3;l++)
            {
            	I1_2(row,col)+=KirchhoffST(n,l);
                row++;
            }
            col++;
        }
       }


}

void FSMicromorphic3DT:: Form_I1_3()
{
    int row=0;
    int col=0;
    I1_3=0.0;

/*    for(int m=0;m<3;m++)
    {
        for(int l=0;l<3;l++)
        {
            //row=3*l;
        	row=l;
            for(int k=0;k<3;k++)
            {
            	I1_3(row,col)+= KirchhoffST(k,m);
            	row=row+3;
            	//row=row+1;
            }
            col++;
        }
     }*/

    for(int L=0;L<3;L++)
    {
 	   for(int l=0;l<3;l++)
 	   {
 		   row=l;
 		   for(int K=0;K<3;K++)
 		   {
 			  I1_3(row,col)=SPK(K,L);
 			   row=row+3;
 		   }
 		   col++;
 	   }
    }


}

void FSMicromorphic3DT:: Form_I1_4()
{
    int row=0;
    int col=0;
    I1_4=0.0;
/*    for(int n=0;n<3;n++)
    {
        for(int i=0;i<3;i++)
        {
            //
            row=0;
            for(int k=0;k<3;k++)
            {
                for(int l=0;l<3;l++)
                {
                    //summation over the same term
                    for(int K=0;K<3;K++)
                    {
                    	for(int M=0;M<3;M++)
                    	{
                    		I1_4(row,col)+=fDeformation_Gradient(k,K)*fDeformation_Gradient(n,M)*fDeformation_Gradient(i,M)
                                          *fDeformation_Gradient(l,K);
                    	}
                    }
                    row++;
                }
            }
            col++;
        }
    }
*/


	for(int M=0;M<3;M++)
	{
		for(int i=0;i<3;i++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int l=0;l<3;l++)
				{
					I1_4(row,col)=fDeformation_Gradient(l,L)*fDeformation_Gradient(i,M);
					row++;
				}
			}
			col++;
		}
	}


}

void FSMicromorphic3DT:: Form_I1_5()
{
    int row;
    int col=0;
    I1_5=0.0;
/*    for(int n=0;n<3;n++)
    {
        for(int i=0;i<3;i++)
        {
            row=0;
            //
            for(int k=0;k<3;k++)
            {
                for(int l=0;l<3;l++)
                {
                    //summation
                    for(int K=0;K<3;K++)
                    {
                        for(int L=0;L<3;L++)
                        {
                            I1_5(row,col)+=fDeformation_Gradient(k,K)*fDeformation_Gradient(n,K)
                                         *fDeformation_Gradient(i,L)*fDeformation_Gradient(l,L);
                        }
                    }
                    row++;
                }
            }
            col++;
        }
    }
    */

	for(int K=0;K<3;K++)
	{
		for(int i=0;i<3;i++)
		{

			row=K*3;
			for(int l=0;l<3;l++)
			{
				for(int L=0;L<3;L++)
				{
					I1_5(row,col)+=fDeformation_Gradient(l,L)
											   *fDeformation_Gradient(i,L);
				}
				row++;
			}
			col++;
		}
	}



}

void FSMicromorphic3DT:: Form_I1_6()
{
    int row;
    int col=0;
    I1_6=0.0;
/*    for(int n=0;n<3;n++)
    {
        for(int i=0;i<3;i++)
        {
            row=0;
            //
            for(int k=0;k<3;k++)
            {
                for(int l=0;l<3;l++)
                {
                	//summation
                    for(int K=0;K<3;K++)
                    {
                        for(int L=0;L<3;L++)
                        {
                            I1_6(row,col)+=fDeformation_Gradient(k,K)*fDeformation_Gradient(i,K)
                                          *fDeformation_Gradient(n,L)*fDeformation_Gradient(l,L);
                        }
                    }
                row++;
                }
            }
            col++;
        }
    }*/

	for(int L=0;L<3;L++)
	{
		for(int i=0;i<3;i++)
		{
			row=0;
			for(int K=0;K<3;K++)
			{
				for(int l=0;l<3;l++)
				{
					I1_6(row,col)=fDeformation_Gradient(l,L)*fDeformation_Gradient(i,K);
					row++;
				}
			}
			col++;
		}
	}

}

void FSMicromorphic3DT:: Form_I1_7()
{
    int row;
    int col=0;
    I1_7=0.0;
    for(int n=0;n<3;n++)
    {
        for(int i=0;i<3;i++)
        {
            row=0;
            for(int k=0;k<3;k++)
            {
                for(int l=0;l<3;l++)
                {
                    //summation
                    for(int K=0;K<3;K++)
                    {
                        for(int M=0;M<3;M++)
                        {
                            I1_7(row,col)+=fDeformation_Gradient(k,K)*fDeformation_Gradient(n,M)
                                          *ChiM(i,M)*fDeformation_Gradient(l,K);
                        }
                    }
                    row++;
                }
            }
            col++;
        }
    }
}

void FSMicromorphic3DT:: Form_I2_1()
{
    int row;
    int col=0;
    I2_1=0.0;
    for(int M=0;M<3;M++)
    {
        for(int i=0;i<3;i++)
        {
            row=0;
            for(int k=0;k<3;k++)
            {
                for(int l=0;l<3;l++)
                {
                    //summation
                    for(int K=0;K<3;K++)
                    {
                        I2_1(row,col)+=fDeformation_Gradient(k,K)*fDeformation_Gradient(i,M)
                                      *fDeformation_Gradient(l,K);
                    }
                row++;
                }
            }
            col++;
        }
    }
}


void FSMicromorphic3DT:: Form_I1_8()
{
    int row;
    int col=0;
    I1_8=0.0;
    for(int n=0;n<3;n++)
    {
        for(int i=0;i<3;i++)
        {
            //
            row=0;
            for(int k=0;k<3;k++)
            {
                for (int l=0;l<3;l++)
                {
                    //summation
                    for(int K=0;K<3;K++)
                    {
                        for(int L=0;L<3;L++)
                        {
                            I1_8(row,col)+=fDeformation_Gradient(k,K)*fDeformation_Gradient(n,K)
                                          *ChiM(i,L)*fDeformation_Gradient(l,L);
                        }
                    }
                    row++;
                }
            }
            col++;
        }
    }

}

void FSMicromorphic3DT:: Form_I2_2()
{
    int row;
    int col=0;
    I2_2=0.0;
    for(int L=0;L<3;L++)
    {
        for(int i=0;i<3;i++)
        {
            //
            row=0;
            for(int k=0;k<3;k++)
            {
                for(int l=0;l<3;l++)
                {
                    //summation
                	for(int K=0;K<3;K++)
                     {
                         I2_2(row,col)+=fDeformation_Gradient(k,K)*fDeformation_Gradient(i,K)
                                       *fDeformation_Gradient(l,L);
                     }
                     row++;
                }
            }
            col++;
        }
    }
}

void FSMicromorphic3DT:: Form_I1_9()
{
    int row;
    int col=0;
    I1_9=0.0;
    for(int n=0;n<3;n++)
    {
        for(int i=0;i<3;i++)
        {
            //
            row=0;
            for(int k=0;k<3;k++)
            {
                for(int l=0;l<3;l++)
                {
                    //summation
                    for(int K=0;K<3;K++)
                    {
                        for(int L=0;L<3;L++)
                        {
                            I1_9(row,col)+=fDeformation_Gradient(k,K)*fDeformation_Gradient(n,L)
                                          *ChiM(i,K)*fDeformation_Gradient(l,L);
                        }
                    }
                    row++;
                }
            }
            col++;
        }
    }
}


void FSMicromorphic3DT:: Form_I2_3()
{
    int row;
    int col=0;
    I2_3=0.0;
    for(int K=0;K<3;K++)
    {
        for(int i=0;i<3;i++)
        {
            //
            row=0;
            for(int k=0;k<3;k++)
            {
                for(int l=0;l<3;l++)
                {
                    //summation
                    for(int L=0;L<3;L++)
                    {
                        I2_3(row,col)+=fDeformation_Gradient(k,K)*fDeformation_Gradient(i,L)*fDeformation_Gradient(l,L);
                    }
                    row++;
                }
            }
            col++;
        }
    }
}

void FSMicromorphic3DT:: Form_fV1()
{
    int row=0;
    fV1=0.0;
    //fTemp_matrix_nsd_x_nsd=0.0;
    //fTemp_matrix_nsd_x_nsd.MultABCT(fDeformation_Gradient,SPK,fDeformation_Gradient);
/*    for(int l=0;l<3;l++)
    {
        for(int k=0;k<3;k++)
        {
            fV1[row]=KirchhoffST(l,k);
            row++;
        }
    }*/

    fTemp_matrix_nsd_x_nsd.MultABT(SPK,fDeformation_Gradient);
    for (int K=0;K<3;K++)
    {
    	for(int l=0;l<3;l++)
    	{
    		fV1[row]=fTemp_matrix_nsd_x_nsd(K,l);
    		row++;
    	}
    }


}

void FSMicromorphic3DT:: Form_fV2()
{
	int row=0;
	fV2=0.0;
	Temp_SPK=0.0;
	Temp_SPK.MultABCT(fDeformation_Gradient,SIGMA_S,fDeformation_Gradient);
	//Temp_SPK*=-1;
	for(int m=0;m<3;m++)
	{
		for(int l=0;l<3;l++)
		{
			fV2[row]=Temp_SPK(m,l);//this is s_sigma
			row++;
		}
	}
}


void FSMicromorphic3DT:: Form_fV3()
{
	int row=0;
	fV3=0.0;
	fTemp_tensor_n_sd_x_n_sd_x_nsd=0.0;
/*	for(int k=0;k<3;k++)
	{
		for(int l=0;l<3;l++)
		{
			for(int m=0;m<3;m++)
			{
				//summation
				for(int K=0;K<3;K++)
				{
					for(int L=0;L<3;L++)
					{
						for(int M=0;M<3;M++)
						{
							fTemp_tensor_n_sd_x_n_sd_x_nsd(k,l,m)+=fDeformation_Gradient(k,K)
																  *fDeformation_Gradient(l,L)
																  *fMKLM(K,L,M)
																  *ChiM(m,M);
						}
					}
				}
			}
		}
	}*/

	for(int K=0;K<3;K++)
	{
		for(int l=0;l<3;l++)
		{
			for(int m=0;m<3;m++)
			{
				//summation
 					for(int L=0;L<3;L++)
					{
						for(int M=0;M<3;M++)
						{
							fTemp_tensor_n_sd_x_n_sd_x_nsd(K,l,m)+=fDeformation_Gradient(l,L)
																  *fMKLM(K,L,M)
																  *ChiM(m,M);
						}
					}
 			}
		}
	}


	for(int l=0;l<3;l++)
	{
		for(int m=0;m<3;m++)
		{
			for(int K=0;K<3;K++)
			{
				fV3[row]+=fTemp_tensor_n_sd_x_n_sd_x_nsd(K,l,m);
				row++;
			}
		}
	}

}

void FSMicromorphic3DT:: Form_SIGMA_S()
{
	double trLST=0.0;
	double trcE =0.0;
	double scale=0.0;
	SIGMA_S=0.0;
	fTemp_matrix_nsd_x_nsd=0.0;

	trLST= LagrangianStn(0,0)+LagrangianStn(1,1)+LagrangianStn(2,2);
	trcE = MicroStnTensor(0,0)+MicroStnTensor(1,1)+MicroStnTensor(2,2);

	SIGMA_S =fIdentity_matrix;
	scale   =(LagrangianStn(0,0)+LagrangianStn(1,1)+LagrangianStn(2,2))*fMaterial_Params[kTau];//or trLST
	SIGMA_S*=scale;

	fTemp_matrix_nsd_x_nsd=LagrangianStn;
	scale=2*fMaterial_Params[kSigma_const];
	fTemp_matrix_nsd_x_nsd*=scale;
	SIGMA_S+=fTemp_matrix_nsd_x_nsd;

	fTemp_matrix_nsd_x_nsd=fIdentity_matrix;
	scale=trcE*(fMaterial_Params[kEta]-fMaterial_Params[kTau]);
	fTemp_matrix_nsd_x_nsd*=scale;
	SIGMA_S+=fTemp_matrix_nsd_x_nsd;

	fTemp_matrix_nsd_x_nsd=MicroStnTensor;
	scale=(fMaterial_Params[kNu]-fMaterial_Params[kSigma_const]);
	fTemp_matrix_nsd_x_nsd*=scale;
	SIGMA_S+=fTemp_matrix_nsd_x_nsd;

	fTemp_matrix_nsd_x_nsd.Transpose(MicroStnTensor);
	scale=(fMaterial_Params[kKappa]-fMaterial_Params[kSigma_const]);
	fTemp_matrix_nsd_x_nsd*=scale;
	SIGMA_S+=fTemp_matrix_nsd_x_nsd;


}


void FSMicromorphic3DT:: Form_fFJ()
{
	int row=0;
	int col=0;
	fFJ=0.0;
	Temp_SPK=0.0;
/*	fTemp_matrix_nsd_x_nsd=0.0;
	fTemp_matrix_nsd_x_nsd=SIGMA_S;
	Temp_SPK.MultABCT(fDeformation_Gradient,fTemp_matrix_nsd_x_nsd,fDeformation_Gradient);*/
/*
	for(int n=0;n<3;n++)
	{
		for(int m=0;m<3;m++)
		{
			row=m;
			for(int l=0;l<3;l++)
			{
				fFJ(row,col)+=Temp_SPK(n,l);
				row=row+3;
			}
			col++;
		}
	}*/

	for(int M=0;M<3;M++)
	{
		for(int m=0;m<3;m++)
		{
			row=m;
			for(int l=0;l<3;l++)
			{
				//summation
				for(int L=0;L<3;L++)
				{
					fFJ(row,col)+=SIGMA_S(M,L)*fDeformation_Gradient(l,L);
				}
				row=row+3;
			}
			col++;
		}
	}

}


void FSMicromorphic3DT:: Form_fJF()
{
	int row=0;
	int col=0;
	fJF=0.0;
/*	Temp_SPK=0.0;
	fTemp_matrix_nsd_x_nsd=0.0;
	fTemp_matrix_nsd_x_nsd=SIGMA_S;
	Temp_SPK.MultABCT(fDeformation_Gradient,fTemp_matrix_nsd_x_nsd,fDeformation_Gradient);
	for(int n=0;n<3;n++)
	{
		for(int l=0;l<3;l++)
		{
			row=3*l;
			for(int m=0;m<3;m++)
			{
				fJF(row,col)+=Temp_SPK(m,n);
				row++;
			}
			col++;
		}
	}*/

	for(int L=0;L<3;L++)
	{
		for(int l=0;l<3;l++)
		{
			row=3*l;
			for(int m=0;m<3;m++)
			{
				//summation
				for(int M=0;M<3;M++)
				{
					fJF(row,col)+=fDeformation_Gradient(m,M)*SIGMA_S(M,L);
				}
				row++;
			}
			col++;
		}
	}


}

void FSMicromorphic3DT:: Form_fJ1_1()
{
	int row=0;
	int col=0;
	fJ1_1=0.0;
/*	for(int n=0;n<3;n++)
	{
		for(int i=0;i<3;i++)
		{
			row=0;
			for(int l=0;l<3;l++)
			{
				for(int m=0;m<3;m++)
				{
					//summation
					for(int M=0;M<3;M++)
					{
						for(int N=0;N<3;N++)
						{
							fJ1_1(row,col)+=fDeformation_Gradient(m,M)*fDeformation_Gradient(n,N)
										 *fDeformation_Gradient(i,N)*fDeformation_Gradient(l,M);
						}
					}
					row++;
				}
			}
			col++;
		}
	}*/
	for(int M=0;M<3;M++)
	{
		for(int i=0;i<3;i++)
		{
			row=0;
			for(int l=0;l<3;l++)
			{
				for(int m=0;m<3;m++)
				{
					//summation
						for(int K=0;K<3;K++)
						{
							fJ1_1(row,col)+=fDeformation_Gradient(m,K)*fDeformation_Gradient(i,M)
										 *fDeformation_Gradient(l,K);
						}

					row++;
				}
			}
			col++;
		}
	}


}
void FSMicromorphic3DT::Form_fJ1_2()
{
	int row=0;
	int col=0;
	fJ1_2=0.0;
	for(int n=0;n<3;n++)
	{
		for(int i=0;i<3;i++)
		{
			row=0;
			for(int l=0;l<3;l++)
			{
				for(int m=0;m<3;m++)
				{
					//summation
					for(int M=0;M<3;M++)
					{
						for(int L=0;L<3;L++)
						{
							fJ1_2(row,col)+=fDeformation_Gradient(m,M)*fDeformation_Gradient(n,M)
										 *fDeformation_Gradient(i,L)*fDeformation_Gradient(l,L);
						}
					}
					row++;
				}
			}
			col++;
		}
	}

}


void FSMicromorphic3DT::Form_fJ1_3()
{
	int row=0;
	int col=0;
	fJ1_3=0.0;
	for(int n=0;n<3;n++)
	{
		for(int i=0;i<3;i++)
		{
			row=0;
			for(int l=0;l<3;l++)
			{
				for(int m=0;m<3;m++)
				{
					//summation
					for(int M=0;M<3;M++)
					{
						for(int L=0;L<3;L++)
						{
							fJ1_3(row,col)+=fDeformation_Gradient(m,M)*fDeformation_Gradient(i,M)
										 *fDeformation_Gradient(n,L)*fDeformation_Gradient(l,L);
						}
					}
					row++;
				}
			}
			col++;
		}
	}

}

void FSMicromorphic3DT:: Form_fJ1_4()
{
	int row=0;
	int col=0;
	fJ1_4=0.0;
/*		for(int n=0;n<3;n++)
		{
			for(int i=0;i<3;i++)
			{
				row=0;
				for(int l=0;l<3;l++)
				{
					for(int m=0;m<3;m++)
					{
						//summation
						for(int M=0;M<3;M++)
						{
							for(int N=0;N<3;N++)
							{
								fJ1_4(row,col)+=fDeformation_Gradient(m,M)*fDeformation_Gradient(n,N)
											 *ChiM(i,N)*fDeformation_Gradient(l,M);
							}
						}
						row++;
					}
				}
				col++;
			}
		}*/

	for(int M=0;M<3;M++)
	{
		for(int i=0;i<3;i++)
		{
			row=0;
			for(int l=0;l<3;l++)
			{
				for(int m=0;m<3;m++)
				{
					//summation

						for(int K=0;K<3;K++)
						{
							fJ1_4(row,col)+=fDeformation_Gradient(m,K)*ChiM(i,M)*fDeformation_Gradient(l,K);
						}

					row++;
				}
			}
			col++;
		}
	}


}

void FSMicromorphic3DT::Form_fJ2_1()
{
	int row=0;
	int col=0;
	fJ2_1=0.0;
	for(int N=0;N<3;N++)
	{
		for(int i=0;i<3;i++)
		{
			row=0;
			for(int l=0;l<3;l++)
			{
				for(int m=0;m<3;m++)
				{
					//summation
					for(int M=0;M<3;M++)
					{
						fJ2_1(row,col)+=fDeformation_Gradient(m,M)*fDeformation_Gradient(i,N)*fDeformation_Gradient(l,M);
					}
					row++;
				}
			}
			col++;
		}
	}
}


void FSMicromorphic3DT:: Form_fJ1_5()
{
	int row=0;
	int col=0;
	fJ1_5=0.0;
	for(int n=0;n<3;n++)
	{
		for(int i=0;i<3;i++)
		{
			row=0;
			for(int l=0;l<3;l++)
			{
				for(int m=0;m<3;m++)
				{
					//summation
					for(int M=0;M<3;M++)
					{
						for(int L=0;L<3;L++)
						{
							fJ1_5(row,col)+=fDeformation_Gradient(m,M)*fDeformation_Gradient(n,M)
										 *ChiM(i,L)*fDeformation_Gradient(l,L);
						}
					}
					row++;
				}
			}
			col++;
		}
	}

}

void FSMicromorphic3DT:: Form_fJ2_2()
{
	int row=0;
	int col=0;
	fJ2_2=0.0;
	for(int L=0;L<3;L++)
	{
		for(int i=0;i<3;i++)
		{
			//
			row=0;
			for(int l=0;l<3;l++)
			{
				for(int m=0;m<3;m++)
				{
					//summation
					for(int M=0;M<3;M++)
					{
						fJ2_2(row,col)+=fDeformation_Gradient(m,M)*fDeformation_Gradient(i,M)*fDeformation_Gradient(l,L);
					}
				row++;
				}
			}
			col++;
		}
	}
}

void FSMicromorphic3DT:: Form_fJ1_6()
{
	int row=0;
	int col=0;
	fJ1_6=0.0;
	for(int n=0;n<3;n++)
	{
		for(int i=0;i<3;i++)
		{
			row=0;
			for(int l=0;l<3;l++)
			{
				for(int m=0;m<3;m++)
				{
					//summation
					for(int M=0;M<3;M++)
					{
						for(int L=0;L<3;L++)
						{
							fJ1_6(row,col)+=fDeformation_Gradient(m,M)*fDeformation_Gradient(n,L)
										  *ChiM(i,M)*fDeformation_Gradient(l,L);
						}
					}
					row++;
				}
			}
			col++;
		}
	}
}

void FSMicromorphic3DT::Form_fJ2_3()
{
	int row=0;
	int col=0;
	fJ2_3=0.0;
	for(int M=0;M<3;M++)
	{
		for(int i=0;i<3;i++)
		{
			row=0;
			for(int l=0;l<3;l++)
			{
				for(int m=0;m<3;m++)
				{
					//summation
					for(int L=0;L<3;L++)
					{
						fJ2_3(row,col)+=fDeformation_Gradient(m,M)*fDeformation_Gradient(i,L)*fDeformation_Gradient(l,L);
					}
					row++;
				}
			}
			col++;
		}
	}
}

void FSMicromorphic3DT:: Form_GAMMA()
{
	GAMMA=0.0;
	for(int K=0;K<3;K++)
	{
		for(int L=0;L<3;L++)
		{
			for(int M=0;M<3;M++)
			{
				//summation
				for(int i=0;i<3;i++)
				{
					GAMMA(K,L,M)+=fDeformation_Gradient(i,K)*GRAD_CHIM(i,L,M);
				}
			}
		}
	}
}


void FSMicromorphic3DT:: Form_fMKLM()
{
	fMKLM=0.0;
	for(int K=0;K<3;K++ )
	{
		for(int L=0;L<3;L++)
		{
			for(int M=0;M<3;M++)
			{
				fMKLM(K,L,M)=fMaterial_Params[kTau8]*(GAMMA(K,L,M)+GAMMA(M,L,K));
			}
		}
	}
}

void FSMicromorphic3DT:: Form_fEtaM()
{
	int row=0;
	int col=0;
	fEtaM=0.0;
	for(int k=0;k<3;k++)
	{
		for(int i=0;i<3;i++)
		{
			//
			row=i;
			for(int l=0;l<3;l++)
			{
				for(int m=0;m<3;m++)
				{
					//summation
					for(int K=0;K<3;K++)
					{
						for(int L=0;L<3;L++)
						{
							for(int M=0;M<3;M++)
							{
								fEtaM(row,col)+=fDeformation_Gradient(k,K)
											  *fDeformation_Gradient(l,L)
											  *fMKLM(K,L,M)
											  *ChiM(m,M);
							}
						}
					}
					row=row+3;
				}
			}
			col++;
		}
	}
}


void FSMicromorphic3DT:: Form_fFM()
{
	int row=0;
	int col=0;
	fFM=0.0;
	for(int n=0;n<3;n++)
	{
		for(int k=0;k<3;k++)
		{
			//
			row=k;
			for(int l=0;l<3;l++)
			{
				for(int m=0;m<3;m++)
				{
					//summation
					for(int K=0;K<3;K++)
					{
						for(int L=0;L<3;L++)
						{
							for(int M=0;M<3;M++)
							{
								fEtaM(row,col)+=fDeformation_Gradient(n,K)
											  *fDeformation_Gradient(l,L)
											  *fMKLM(K,L,M)
											  *ChiM(m,M);
							}
						}
					}
					row=row+3;
				}
			}
			col++;
		}
	}
}

void FSMicromorphic3DT:: Form_fMF()
{
	int row=0;
	int col=0;
	fMF=0.0;
/*	for(int n=0;n<3;n++)
	{
		for(int l=0;l<3;l++)
		{
			//
			row=l*9;
			for(int m=0;m<3;m++)
			{
				for(int k=0;k<3;k++)
				{
					//summation
					for(int K=0;K<3;K++)
					{
						for(int L=0;L<3;L++)
						{
							for(int M=0;M<3;M++)
							{
								for(int R=0;R<3;R++)
								{
								fMF(row,col)+=fDeformation_Gradient(k,K)
											 *fDeformation_Gradient(n,L)
											 *fMKLM(K,L,M)
											 *ChiM(m,M)
											 *fDeformation_Gradient_Inverse(R,k);
								}
							}
						}
					}
					row++;
				}
			}
			col++;
		}
	}*/

	for(int L=0;L<3;L++)
	{
		for(int l=0;l<3;l++)
		{
			//
			row=l*9;
			for(int m=0;m<3;m++)
			{
				for(int K=0;K<3;K++)
				{
					//summation
					for(int M=0;M<3;M++)
					{
					fMF(row,col)+=ChiM(m,M)*fMKLM(K,L,M);
					}
					row++;
				}
			}
			col++;

		}
	}



}

void FSMicromorphic3DT:: Form_fMchi()
{
	int row=0;
	int col=0;
	fMchi=0.0;
/*	for(int M=0;M<3;M++)
	{
		for(int m=0;m<3;m++)
		{
			row=m*3;
			for(int l=0;l<3;l++)
			{
				for(int k=0;k<3;k++)
				{
					//summation
					for(int K=0;K<3;K++)
					{
						for(int L=0;L<3;L++)
						{
							for(int R=0;R<3;R++)
							{
							fMchi(row,col)+=fDeformation_Gradient(k,K)*fDeformation_Gradient(l,L)*fMKLM(K,L,M)
										   *fDeformation_Gradient_Inverse(R,k);
							}
						}
					}
					row++;
				}
				row=row+6;
			}
			col++;
		}
	}*/

	for(int M=0;M<3;M++)
	{
		for(int m=0;m<3;m++)
		{
			row=m*3;
			for(int l=0;l<3;l++)
			{
				for(int K=0;K<3;K++)
				{
					//summation
					for(int L=0;L<3;L++)
					{
						fMchi(row,col)+=fDeformation_Gradient(l,L)*fMKLM(K,L,M);
					}
					row++;
				}
				row=row+6;
			}
			col++;
		}
	}




}

void FSMicromorphic3DT:: Form_fMpu_1()
{
	int row=0;
	int col=0;
	fMpu_1=0.0;
/*	for(int n=0;n<3;n++)
	{
		for(int i=0;i<3;i++)
		{
			//
			row=0;
			for(int l=0;l<3;l++)
			{
				for(int m=0;m<3;m++)
				{
					for(int k=0;k<3;k++)
					{
						//summation
						for(int K=0;K<3;K++)
						{
							for(int L=0;L<3;L++)
							{
								for(int M=0;M<3;M++)
								{
									for(int R=0;R<3;R++)
									{
									fMpu_1(row,col)+=fDeformation_Gradient(k,K)
													*fDeformation_Gradient(l,L)
													*fDeformation_Gradient(n,K)
													*GRAD_CHIM(i,L,M)
													*ChiM(m,M)
													*fDeformation_Gradient_Inverse(K,k)!this is worn K is repeated 3 times but already commented out);
									}
								}
							}
						}
						row++;
					}
				}
			}
			col++;
		}
	}*/
	for(int K=0;K<3;K++)
	{
		for(int i=0;i<3;i++)
		{
			row=K;
			for(int l=0;l<3;l++)
			{
				for(int m=0;m<3;m++)
				{
					//summation
					for(int L=0;L<3;L++)
					{
						for(int M=0;M<3;M++)
						{
							fMpu_1(row,col)+=fDeformation_Gradient(l,L)
											*GRAD_CHIM(i,L,M)
											*ChiM(m,M);
						}
					}
					row=row+3;
				}
			}
			col++;
		}
	}

}

void FSMicromorphic3DT:: Form_fMpp_1()
{
	int row=0;
	int col=0;
	fMpp_1=0.0;
/*	for(int L=0;L<3;L++)
	{
		for(int i=0;i<3;i++)
		{
			for(int M=0;M<3;M++)
			{
				//
				row=0;
				for(int l=0;l<3;l++)
				{
					for(int m=0;m<3;m++)
					{
						for(int k=0;k<3;k++)
						{
							//summation
							for(int K=0;K<3;K++)
							{
								fMpp_1(row,col)+=fDeformation_Gradient(k,K)
												*fDeformation_Gradient(l,L)
												*fDeformation_Gradient(i,K)
												*ChiM(m,M);
							}
						row++;
						}
					}
				}
				col++;
			}
		}
	}*/

	for(int L=0;L<3;L++)
	{
		for(int i=0;i<3;i++)
		{
			for(int M=0;M<3;M++)
			{
				//
				row=0;
				for(int l=0;l<3;l++)
				{
					for(int m=0;m<3;m++)
					{
						for(int K=0;K<3;K++)
						{
							//summation

						 fMpp_1(row,col)=fDeformation_Gradient(l,L)
										*fDeformation_Gradient(i,K)
										*ChiM(m,M);

						row++;
						}
					}
				}
				col++;
			}
		}
	}


}

void FSMicromorphic3DT:: Form_fMpu_2()
{
	int row=0;
	int col=0;
	fMpu_2=0.0;
/*	for(int n=0;n<3;n++)
	{
		for(int i=0;i<3;i++)
		{
			//
			row=0;
			for(int l=0;l<3;l++)
			{
				for(int m=0;m<3;m++)
				{
					for(int k=0;k<3;k++)
					{
						//summation
						for(int K=0;K<3;K++)
						{
							for(int L=0;L<3;L++)
							{
								for(int M=0;M<3;M++)
								{
									fMpu_2(row,col)+=fDeformation_Gradient(k,K)
													*fDeformation_Gradient(l,L)
													*fDeformation_Gradient(n,M)
													*GRAD_CHIM(i,L,K)
													*ChiM(m,M);
								}
							}
						}
						row++;
					}
				}
			}
			col++;
		}
	}*/
	for(int M=0;M<3;M++)
		{
			for(int i=0;i<3;i++)
			{
				row=0;
				for(int l=0;l<3;l++)
				{
					for(int m=0;m<3;m++)
					{
						for(int K=0;K<3;K++)
						{
							//summation
							for(int L=0;L<3;L++)
							{
								fMpu_2(row,col)+=fDeformation_Gradient(l,L)
												*ChiM(m,M)
												*GRAD_CHIM(i,L,K);
							}
							row++;
						}
					}
				}
				col++;
			}
		}


}

void FSMicromorphic3DT:: Form_fMpp_2()
{
	int row=0;
	int col=0;
	fMpp_2=0.0;
/*	for(int L=0;L<3;L++)
	{
		for(int i=0;i<3;i++)
		{
			for(int K=0;K<3;K++)
			{
				//
				row=0;
				for(int l=0;l<3;l++)
				{
					for(int m=0;m<3;m++)
					{
						for(int k=0;k<3;k++)
						{
							//summation
							for(int M=0;M<3;M++)
							{
								fMpp_2(row,col)+=fDeformation_Gradient(k,K)
												*fDeformation_Gradient(l,L)
												*fDeformation_Gradient(i,M)
												*ChiM(m,M);
							}
							row++;
						}
					}
				}
				col++;
			}
		}
	}*/

	for(int L=0;L<3;L++)
		{
			for(int i=0;i<3;i++)
			{
				for(int K=0;K<3;K++)
				{
					row=K;
					for(int l=0;l<3;l++)
					{
						for(int m=0;m<3;m++)
						{
							//summation
							for(int M=0;M<3;M++)
							{
								fMpp_2(row,col)+=fDeformation_Gradient(l,L)
												*ChiM(m,M)
												*fDeformation_Gradient(i,M);
							}

						row=row+3;
						}
					}
					col++;
				}
			}
		}
}

void FSMicromorphic3DT:: Form_Jmat()
{
	int row=0;
	int col=0;
	Jmat=0.0;

	for(int A=0;A<3;A++)
	{
		for(int i=0;i<3;i++)
		{
			row=0;
			for(int K=0;K<3;K++)
			{
				for(int l=0;l<3;l++)
				{
					for(int L=0;L<3;L++)
					{
						Jmat(row,col)+=fDeformation_Gradient(l,L)*SPK(K,L)*fIdentity_matrix(A,i);
					}
					row++;
				}
			}
			col++;
		}
	}


}


////////////////////////////////////////////////////////////////
//////////////FINITE STRAIN MATRICES ENDS//////////////////
////////////////////////////////////////////////////////////////

/*void FSMicromorphic3DT:: Extract_six_values_from_symmetric_tensor(const dMatrixT &fTensor,dArrayT& fTemp_six_values)*/
void FSMicromorphic3DT:: Extract_six_values_from_symmetric_tensor(const dMatrixT &fTensor,dArrayT& fTemp_nine_values)
{

/*    fTemp_six_values[0]=fTensor(0,0);
    fTemp_six_values[1]=fTensor(1,1);
    fTemp_six_values[2]=fTensor(2,2);
    fTemp_six_values[3]=fTensor(1,2);
    fTemp_six_values[4]=fTensor(2,0);
    fTemp_six_values[5]=fTensor(0,1);*/

    fTemp_nine_values[0]=fTensor(0,0);//sigma11
    fTemp_nine_values[1]=fTensor(1,1);//sigma22
    fTemp_nine_values[2]=fTensor(2,2);//sigma33
    fTemp_nine_values[3]=fTensor(0,1);//sigma12
    fTemp_nine_values[4]=fTensor(0,2);//sigma13
    fTemp_nine_values[5]=fTensor(1,0);//sigma21
    fTemp_nine_values[6]=fTensor(1,2);//sigma23
    fTemp_nine_values[7]=fTensor(2,0);//sigma31
    fTemp_nine_values[8]=fTensor(2,1);//sigma32

}

void FSMicromorphic3DT::Put_values_In_dArrayT_vector(const dArray2DT &f2DArrayT,const int& e,const int& IP,dArrayT& fArrayT)
{
/*    fArrayT[0]=f2DArrayT(e,IP*6+0);
    fArrayT[1]=f2DArrayT(e,IP*6+1);
    fArrayT[2]=f2DArrayT(e,IP*6+2);
    fArrayT[3]=f2DArrayT(e,IP*6+3);
    fArrayT[4]=f2DArrayT(e,IP*6+4);
    fArrayT[5]=f2DArrayT(e,IP*6+5);*/

    fArrayT[0]=f2DArrayT(e,IP*9+0);
    fArrayT[1]=f2DArrayT(e,IP*9+1);
    fArrayT[2]=f2DArrayT(e,IP*9+2);
    fArrayT[3]=f2DArrayT(e,IP*9+3);
    fArrayT[4]=f2DArrayT(e,IP*9+4);
    fArrayT[5]=f2DArrayT(e,IP*9+5);
    fArrayT[6]=f2DArrayT(e,IP*9+6);
    fArrayT[7]=f2DArrayT(e,IP*9+7);
    fArrayT[8]=f2DArrayT(e,IP*9+8);

/*    fArrayT[0]=f2DArrayT(e,IP*12+0);
    fArrayT[1]=f2DArrayT(e,IP*12+1);
    fArrayT[2]=f2DArrayT(e,IP*12+2);
    fArrayT[3]=f2DArrayT(e,IP*12+3);
    fArrayT[4]=f2DArrayT(e,IP*12+4);
    fArrayT[5]=f2DArrayT(e,IP*12+5);
    fArrayT[6]=f2DArrayT(e,IP*12+6);
    fArrayT[7]=f2DArrayT(e,IP*12+7);
    fArrayT[8]=f2DArrayT(e,IP*12+8);
    fArrayT[9]=f2DArrayT(e,IP*12+9);
    fArrayT[10]=f2DArrayT(e,IP*12+10);
    fArrayT[11]=f2DArrayT(e,IP*12+11);*/


}


void FSMicromorphic3DT::Put_values_In_Array(const dArray2DT &f2DArrayT,const int& e,const int& IP,dArrayT& fArrayT)
{


    fArrayT[0]=f2DArrayT(e,IP*3+0);
    fArrayT[1]=f2DArrayT(e,IP*3+1);
    fArrayT[2]=f2DArrayT(e,IP*3+2);


}



//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
////////////////////// FINISH HERE FOR MICROMORPHIC 3-D CASE/////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
/*
void FSMicromorphic3DT::Form_Varpi_temp_matrix()
{
    double N_A_1I, N_A_2I, N_A_3I;
    int j, temp_j;
    j = 0 ;
    for (int A=1; A <= n_en_displ; A++)
    {
    temp_j = j;
    for (int i=0; i<3; i++)
    {
        j = temp_j;
        for (int n=1; n<=3; n++)
        {
        switch (i+1)
        {
        case 1:
        {
            N_A_1I = fShapeDisplGradGrad(0,A-1);
            N_A_2I = fShapeDisplGradGrad(5,A-1);
            N_A_3I = fShapeDisplGradGrad(4,A-1);
        }break;
        case 2:
        {
            N_A_1I = fShapeDisplGradGrad(5,A-1);
            N_A_2I = fShapeDisplGradGrad(1,A-1);
            N_A_3I = fShapeDisplGradGrad(3,A-1);
        }break;
        case 3:
        {
            N_A_1I = fShapeDisplGradGrad(4,A-1);
            N_A_2I = fShapeDisplGradGrad(3,A-1);
            N_A_3I = fShapeDisplGradGrad(2,A-1);
        }break;
        }
        fVarpi_temp_matrix(i,j) = N_A_1I * fDeformation_Gradient_Inverse(0, n-1)
            + N_A_2I * fDeformation_Gradient_Inverse(1, n-1)
            + N_A_3I * fDeformation_Gradient_Inverse(2, n-1);
        j += 1;
        }
    }
    }
}

*/
/*
void FSMicromorphic3DT::Form_Gradient_t_of_solid_shape_functions(const dMatrixT &fShapeDisplGrad_temp)
{
    fShapeDisplGrad_t = 0.0;
    for (int i=0; i<27; i++)
    {
    fShapeDisplGrad_t(0,i*3) = fShapeDisplGrad_temp(0,i);
    fShapeDisplGrad_t(1,i*3) = fShapeDisplGrad_temp(1,i);
    fShapeDisplGrad_t(2,i*3) = fShapeDisplGrad_temp(2,i);

    fShapeDisplGrad_t(3,1+i*3) = fShapeDisplGrad_temp(0,i);
    fShapeDisplGrad_t(4,1+i*3) = fShapeDisplGrad_temp(1,i);
    fShapeDisplGrad_t(5,1+i*3) = fShapeDisplGrad_temp(2,i);

    fShapeDisplGrad_t(6,2+i*3) = fShapeDisplGrad_temp(0,i);
    fShapeDisplGrad_t(7,2+i*3) = fShapeDisplGrad_temp(1,i);
    fShapeDisplGrad_t(8,2+i*3) = fShapeDisplGrad_temp(2,i);
    }

}

void FSMicromorphic3DT::Form_Im_temp_matrix()
{
    fIm_temp_matrix = 0.0;
    fIm_temp_matrix(0,0) = fKirchhoff_tensor(0,0);
    fIm_temp_matrix(1,0) = fKirchhoff_tensor(1,0);
    fIm_temp_matrix(2,0) = fKirchhoff_tensor(2,0);

    fIm_temp_matrix(3,1) = fKirchhoff_tensor(0,0);
    fIm_temp_matrix(4,1) = fKirchhoff_tensor(1,0);
    fIm_temp_matrix(5,1) = fKirchhoff_tensor(2,0);

    fIm_temp_matrix(6,2) = fKirchhoff_tensor(0,0);
    fIm_temp_matrix(7,2) = fKirchhoff_tensor(1,0);
    fIm_temp_matrix(8,2) = fKirchhoff_tensor(2,0);

    fIm_temp_matrix(0,3) = fKirchhoff_tensor(0,1);
    fIm_temp_matrix(1,3) = fKirchhoff_tensor(1,1);
    fIm_temp_matrix(2,3) = fKirchhoff_tensor(2,1);

    fIm_temp_matrix(3,4) = fKirchhoff_tensor(0,1);
    fIm_temp_matrix(4,4) = fKirchhoff_tensor(1,1);
    fIm_temp_matrix(5,4) = fKirchhoff_tensor(2,1);

    fIm_temp_matrix(6,5) = fKirchhoff_tensor(0,1);
    fIm_temp_matrix(7,5) = fKirchhoff_tensor(1,1);
    fIm_temp_matrix(8,5) = fKirchhoff_tensor(2,1);

    fIm_temp_matrix(0,6) = fKirchhoff_tensor(0,2);
    fIm_temp_matrix(1,6) = fKirchhoff_tensor(1,2);
    fIm_temp_matrix(2,6) = fKirchhoff_tensor(2,2);

    fIm_temp_matrix(3,7) = fKirchhoff_tensor(0,2);
    fIm_temp_matrix(4,7) = fKirchhoff_tensor(1,2);
    fIm_temp_matrix(5,7) = fKirchhoff_tensor(2,2);

    fIm_temp_matrix(6,8) = fKirchhoff_tensor(0,2);
    fIm_temp_matrix(7,8) = fKirchhoff_tensor(1,2);
    fIm_temp_matrix(8,8) = fKirchhoff_tensor(2,2);
}

void FSMicromorphic3DT::Form_Hbar_temp_matrix()
{
    fHbar_temp_matrix =0.0;
    fHbar_temp_matrix(0,0) = fLeft_Cauchy_Green_tensor(0,0);
    fHbar_temp_matrix(3,0) = fLeft_Cauchy_Green_tensor(0,1);
    fHbar_temp_matrix(6,0) = fLeft_Cauchy_Green_tensor(0,2);

    fHbar_temp_matrix(1,1) = fLeft_Cauchy_Green_tensor(0,0);
    fHbar_temp_matrix(4,1) = fLeft_Cauchy_Green_tensor(0,1);
    fHbar_temp_matrix(7,1) = fLeft_Cauchy_Green_tensor(0,2);

    fHbar_temp_matrix(2,2) = fLeft_Cauchy_Green_tensor(0,0);
    fHbar_temp_matrix(5,2) = fLeft_Cauchy_Green_tensor(0,1);
    fHbar_temp_matrix(8,2) = fLeft_Cauchy_Green_tensor(0,2);

    fHbar_temp_matrix(0,3) = fLeft_Cauchy_Green_tensor(1,0);
    fHbar_temp_matrix(3,3) = fLeft_Cauchy_Green_tensor(1,1);
    fHbar_temp_matrix(6,3) = fLeft_Cauchy_Green_tensor(1,2);

    fHbar_temp_matrix(1,4) = fLeft_Cauchy_Green_tensor(1,0);
    fHbar_temp_matrix(4,4) = fLeft_Cauchy_Green_tensor(1,1);
    fHbar_temp_matrix(7,4) = fLeft_Cauchy_Green_tensor(1,2);

    fHbar_temp_matrix(2,5) = fLeft_Cauchy_Green_tensor(1,0);
    fHbar_temp_matrix(5,5) = fLeft_Cauchy_Green_tensor(1,1);
    fHbar_temp_matrix(8,5) = fLeft_Cauchy_Green_tensor(1,2);

    fHbar_temp_matrix(0,6) = fLeft_Cauchy_Green_tensor(2,0);
    fHbar_temp_matrix(3,6) = fLeft_Cauchy_Green_tensor(2,1);
    fHbar_temp_matrix(6,6) = fLeft_Cauchy_Green_tensor(2,2);

    fHbar_temp_matrix(1,7) = fLeft_Cauchy_Green_tensor(2,0);
    fHbar_temp_matrix(4,7) = fLeft_Cauchy_Green_tensor(2,1);
    fHbar_temp_matrix(7,7) = fLeft_Cauchy_Green_tensor(2,2);

    fHbar_temp_matrix(2,8) = fLeft_Cauchy_Green_tensor(2,0);
    fHbar_temp_matrix(5,8) = fLeft_Cauchy_Green_tensor(2,1);
    fHbar_temp_matrix(8,8) = fLeft_Cauchy_Green_tensor(2,2);
}

void FSMicromorphic3DT::Form_Ell_temp_matrix()
{
    fEll_temp_matrix(0,0) = 0.0;
    fEll_temp_matrix(0,0) = fLeft_Cauchy_Green_tensor(0,0);
    fEll_temp_matrix(1,0) = fLeft_Cauchy_Green_tensor(1,0);
    fEll_temp_matrix(2,0) = fLeft_Cauchy_Green_tensor(2,0);

    fEll_temp_matrix(3,1) = fLeft_Cauchy_Green_tensor(0,0);
    fEll_temp_matrix(4,1) = fLeft_Cauchy_Green_tensor(1,0);
    fEll_temp_matrix(5,1) = fLeft_Cauchy_Green_tensor(2,0);

    fEll_temp_matrix(6,2) = fLeft_Cauchy_Green_tensor(0,0);
    fEll_temp_matrix(7,2) = fLeft_Cauchy_Green_tensor(1,0);
    fEll_temp_matrix(8,2) = fLeft_Cauchy_Green_tensor(2,0);

    fEll_temp_matrix(0,3) = fLeft_Cauchy_Green_tensor(0,1);
    fEll_temp_matrix(1,3) = fLeft_Cauchy_Green_tensor(1,1);
    fEll_temp_matrix(2,3) = fLeft_Cauchy_Green_tensor(2,1);

    fEll_temp_matrix(3,4) = fLeft_Cauchy_Green_tensor(0,1);
    fEll_temp_matrix(4,4) = fLeft_Cauchy_Green_tensor(1,1);
    fEll_temp_matrix(5,4) = fLeft_Cauchy_Green_tensor(2,1);

    fEll_temp_matrix(6,5) = fLeft_Cauchy_Green_tensor(0,1);
    fEll_temp_matrix(7,5) = fLeft_Cauchy_Green_tensor(1,1);
    fEll_temp_matrix(8,5) = fLeft_Cauchy_Green_tensor(2,1);

    fEll_temp_matrix(0,6) = fLeft_Cauchy_Green_tensor(0,2);
    fEll_temp_matrix(1,6) = fLeft_Cauchy_Green_tensor(1,2);
    fEll_temp_matrix(2,6) = fLeft_Cauchy_Green_tensor(2,2);

    fEll_temp_matrix(3,7) = fLeft_Cauchy_Green_tensor(0,2);
    fEll_temp_matrix(4,7) = fLeft_Cauchy_Green_tensor(1,2);
    fEll_temp_matrix(5,7) = fLeft_Cauchy_Green_tensor(2,2);

    fEll_temp_matrix(6,8) = fLeft_Cauchy_Green_tensor(0,2);
    fEll_temp_matrix(7,8) = fLeft_Cauchy_Green_tensor(1,2);
    fEll_temp_matrix(8,8) = fLeft_Cauchy_Green_tensor(2,2);
}



void FSMicromorphic3DT::Form_Im_Prim_temp_matrix()
{
    fIm_Prim_temp_matrix = 0.0;
    fIm_Prim_temp_matrix(0,0) = fKirchhoff_tensor(0,0);
    fIm_Prim_temp_matrix(3,0) = fKirchhoff_tensor(1,0);
    fIm_Prim_temp_matrix(6,0) = fKirchhoff_tensor(2,0);

    fIm_Prim_temp_matrix(1,1) = fKirchhoff_tensor(0,0);
    fIm_Prim_temp_matrix(4,1) = fKirchhoff_tensor(1,0);
    fIm_Prim_temp_matrix(7,1) = fKirchhoff_tensor(2,0);

    fIm_Prim_temp_matrix(2,2) = fKirchhoff_tensor(0,0);
    fIm_Prim_temp_matrix(5,2) = fKirchhoff_tensor(1,0);
    fIm_Prim_temp_matrix(8,2) = fKirchhoff_tensor(2,0);

    fIm_Prim_temp_matrix(0,3) = fKirchhoff_tensor(0,1);
    fIm_Prim_temp_matrix(3,3) = fKirchhoff_tensor(1,1);
    fIm_Prim_temp_matrix(6,3) = fKirchhoff_tensor(2,1);

    fIm_Prim_temp_matrix(1,4) = fKirchhoff_tensor(0,1);
    fIm_Prim_temp_matrix(4,4) = fKirchhoff_tensor(1,1);
    fIm_Prim_temp_matrix(7,4) = fKirchhoff_tensor(2,1);

    fIm_Prim_temp_matrix(2,5) = fKirchhoff_tensor(0,1);
    fIm_Prim_temp_matrix(5,5) = fKirchhoff_tensor(1,1);
    fIm_Prim_temp_matrix(8,5) = fKirchhoff_tensor(2,1);

    fIm_Prim_temp_matrix(0,6) = fKirchhoff_tensor(0,2);
    fIm_Prim_temp_matrix(3,6) = fKirchhoff_tensor(1,2);
    fIm_Prim_temp_matrix(6,6) = fKirchhoff_tensor(2,2);

    fIm_Prim_temp_matrix(1,7) = fKirchhoff_tensor(0,2);
    fIm_Prim_temp_matrix(4,7) = fKirchhoff_tensor(1,2);
    fIm_Prim_temp_matrix(7,7) = fKirchhoff_tensor(2,2);

    fIm_Prim_temp_matrix(2,8) = fKirchhoff_tensor(0,2);
    fIm_Prim_temp_matrix(5,8) = fKirchhoff_tensor(1,2);
    fIm_Prim_temp_matrix(8,8) = fKirchhoff_tensor(2,2);
}

void FSMicromorphic3DT::Form_D_matrix(void)
{
    fD_matrix = 0.0;
    fD_matrix(0,0) = 2*fMaterial_Params[kMu]+ fMaterial_Params[kLambda];
    fD_matrix(1,1) = 2*fMaterial_Params[kMu]+ fMaterial_Params[kLambda];
    fD_matrix(2,2) = 2*fMaterial_Params[kMu]+ fMaterial_Params[kLambda];
    fD_matrix(3,3) = fMaterial_Params[kMu];
    fD_matrix(4,4) = fMaterial_Params[kMu];
    fD_matrix(5,5) = fMaterial_Params[kMu];
    fD_matrix(0,1) = fMaterial_Params[kLambda];
    fD_matrix(1,0) = fMaterial_Params[kLambda];
    fD_matrix(0,2) = fMaterial_Params[kLambda];
    fD_matrix(2,0) = fMaterial_Params[kLambda];
    fD_matrix(1,2) = fMaterial_Params[kLambda];
    fD_matrix(2,1) = fMaterial_Params[kLambda];
}

void FSMicromorphic3DT::Form_B_matrix(void)
{
    fB_matrix = 0.0;
    for(int i=0; i<27; i++)
    {
        fB_matrix(0,i*3)=fShapeDisplGrad_temp(0,i);
        fB_matrix(1,i*3+1)=fShapeDisplGrad_temp(1,i);
        fB_matrix(2,i*3+2)=fShapeDisplGrad_temp(2,i);
        fB_matrix(3,i*3+1)=fShapeDisplGrad_temp(2,i);
        fB_matrix(3,i*3+2)=fShapeDisplGrad_temp(1,i);
        fB_matrix(4,i*3)=fShapeDisplGrad_temp(2,i);
        fB_matrix(4,i*3+2)=fShapeDisplGrad_temp(0,i);
        fB_matrix(5,i*3)=fShapeDisplGrad_temp(1,i);
        fB_matrix(5,i*3+1)=fShapeDisplGrad_temp(0,i);
    }
}
*/
//void FSMicromorphic3DT::Extract_six_values_from_symmetric_tensor(const dMatrixT &fTensor,dArrayT& fTemp_nine_values)


/*void FSMicromorphic3DT::Form_gradv_vector(void)
{
    fShapeDisplGrad.Multx(u_dot_vec,fGradv_vector);
    fDefGradInv_Grad_grad.MultTx(fGradv_vector,fgradv_vector);
}

void FSMicromorphic3DT::Form_Xi_temp_matrix(void)
{
    fXi_temp_matrix = 0.0;

    double c11kl_vkl, c21kl_vkl, c31kl_vkl;
    double c12kl_vkl, c22kl_vkl, c32kl_vkl;
    double c13kl_vkl, c23kl_vkl, c33kl_vkl;
    int col,index;
    c11kl_vkl=c21kl_vkl=c31kl_vkl=0.0;
    c12kl_vkl=c22kl_vkl=c32kl_vkl=0.0;
    c13kl_vkl=c23kl_vkl=c33kl_vkl=0.0;

    for (int k=0; k<3; k++)
    for (int l=0; l<3; l++)
    {
        if (k==0 && l==0)
        col=index=0;
        else if (k==1 && l==0)
        col=index=1;
        else if (k==2 && l==0)
        col=index=2;
        else if (k==0 && l==1)
        col=index=3;
        else if (k==1 && l==1)
        col=index=4;
        else if (k==2 && l==1)
        col=index=5;
        else if (k==0 && l==2)
        col=index=6;
        else if (k==1 && l==2)
        col=index=7;
        else
        col=index=8;

        c11kl_vkl+=fc_matrix(0,col)*fgradv_vector[index];
        c21kl_vkl+=fc_matrix(1,col)*fgradv_vector[index];
        c31kl_vkl+=fc_matrix(2,col)*fgradv_vector[index];
        c12kl_vkl+=fc_matrix(3,col)*fgradv_vector[index];
        c22kl_vkl+=fc_matrix(4,col)*fgradv_vector[index];
        c32kl_vkl+=fc_matrix(5,col)*fgradv_vector[index];
        c13kl_vkl+=fc_matrix(6,col)*fgradv_vector[index];
        c23kl_vkl+=fc_matrix(7,col)*fgradv_vector[index];
        c33kl_vkl+=fc_matrix(8,col)*fgradv_vector[index];
    }
    for (int i=0; i<3; i++)
    {
    fXi_temp_matrix(i*3,i)=c11kl_vkl;
    fXi_temp_matrix(i*3+1,i)=c21kl_vkl;
    fXi_temp_matrix(i*3+2,i)=c31kl_vkl;

    fXi_temp_matrix(i*3,i+3)=c12kl_vkl;
    fXi_temp_matrix(i*3+1,i+3)=c22kl_vkl;
    fXi_temp_matrix(i*3+2,i+3)=c32kl_vkl;

    fXi_temp_matrix(i*3,i+6)=c13kl_vkl;
    fXi_temp_matrix(i*3+1,i+6)=c23kl_vkl;
    fXi_temp_matrix(i*3+2,i+6)=c33kl_vkl;

    }

}

void FSMicromorphic3DT::Form_Varsigma_temp_matrix(void)
{
    fVarsigma_temp_matrix = 0.0;
    int col,index;
    for (int i=0; i<9; i++)
    {
    for (int k=0; k<3; k++)
    {
        if (k==0)
        {
        col=0;
        index=0;
        }
        else if (k==1)
        {
        col=3;
        index=1;
        }
        if (k==2)
        {
        col=6;
        index=2;
        }

        fVarsigma_temp_matrix(i,0) += fc_matrix(i,col)*fgradv_vector[index];
    }

    for (int k=0; k<3; k++)
    {
        if (k==0)
        {
        col=0;
        index=3;
        }
        else if (k==1)
        {
        col=3;
        index=4;
        }
        if (k==2)
        {
        col=6;
        index=5;
        }

        fVarsigma_temp_matrix(i,1) += fc_matrix(i,col)*fgradv_vector[index];
    }

    for (int k=0; k<3; k++)
    {
        if (k==0)
        {
        col=0;
        index=6;
        }
        else if (k==1)
        {
        col=3;
        index=7;
        }
        if (k==2)
        {
        col=6;
        index=8;
        }

        fVarsigma_temp_matrix(i,2) += fc_matrix(i,col)*fgradv_vector[index];
    }

    for (int k=0; k<3; k++)
    {
        if (k==0)
        {
        col=1;
        index=0;
        }
        else if (k==1)
        {
        col=4;
        index=1;
        }
        if (k==2)
        {
        col=7;
        index=2;
        }

        fVarsigma_temp_matrix(i,3) += fc_matrix(i,col)*fgradv_vector[index];
    }

    for (int k=0; k<3; k++)
    {
        if (k==0)
        {
        col=1;
        index=3;
        }
        else if (k==1)
        {
        col=4;
        index=4;
        }
        if (k==2)
        {
        col=7;
        index=5;
        }

        fVarsigma_temp_matrix(i,4) += fc_matrix(i,col)*fgradv_vector[index];
    }

    for (int k=0; k<3; k++)
    {
        if (k==0)
        {
        col=1;
        index=6;
        }
        else if (k==1)
        {
        col=4;
        index=7;
        }
        if (k==2)
        {
        col=7;
        index=8;
        }

        fVarsigma_temp_matrix(i,5) += fc_matrix(i,col)*fgradv_vector[index];
    }

    for (int k=0; k<3; k++)
    {
        if (k==0)
        {
        col=2;
        index=0;
        }
        else if (k==1)
        {
        col=5;
        index=1;
        }
        if (k==2)
        {
        col=8;
        index=2;
        }

        fVarsigma_temp_matrix(i,6) += fc_matrix(i,col)*fgradv_vector[index];
    }

    for (int k=0; k<3; k++)
    {
        if (k==0)
        {
        col=2;
        index=3;
        }
        else if (k==1)
        {
        col=5;
        index=4;
        }
        if (k==2)
        {
        col=8;
        index=5;
        }

        fVarsigma_temp_matrix(i,7) += fc_matrix(i,col)*fgradv_vector[index];
    }

    for (int k=0; k<3; k++)
    {
        if (k==0)
        {
        col=2;
        index=6;
        }
        else if (k==1)
        {
        col=5;
        index=7;
        }
        if (k==2)
        {
        col=8;
        index=8;
        }

        fVarsigma_temp_matrix(i,8) += fc_matrix(i,col)*fgradv_vector[index];
    }
    }

}

void FSMicromorphic3DT::Form_I_ijkl_matrix(void)
{
    double delta_ik,delta_jl,delta_il,delta_jk;
    int row,col;
    for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
        for (int k=0; k<3; k++)
        for (int l=0; l<3; l++)
        {
            if (i==k)
            delta_ik =1;
            else
            delta_ik=0;

            if (j==l)
            delta_jl=1;
            else
            delta_jl=0;

            if (i==l)
            delta_il=1;
            else
            delta_il=0;

            if (j==k)
            delta_jk=1;
            else
            delta_jk=0;

            if (i==0 && j==0)
            row=0;
            else if (i==1 && j==0)
            row=1;
            else if (i==2 && j==0)
            row=2;
            else if (i==0 && j==1)
            row=3;
            else if (i==1 && j==1)
            row=4;
            else if (i==2 && j==1)
            row=5;
            else if (i==0 && j==2)
            row=6;
            else if (i==1 && j==2)
            row=7;
            else
            row=8;

            if (k==0 && l==0)
            col=0;
            else if (k==1 && l==0)
            col=1;
            else if (k==2 && l==0)
            col=2;
            else if (k==0 && l==1)
            col=3;
            else if (k==1 && l==1)
            col=4;
            else if (k==2 && l==1)
            col=5;
            else if (k==0 && l==2)
            col=6;
            else if (k==1 && l==2)
            col=7;
            else
            col=8;

            fI_ijkl_matrix(row,col)= 0.5*(delta_ik*delta_jl+delta_il*delta_jk);

        }
}
*/


/*
void FSMicromorphic3DT::Compute_norm_of_array(double& norm,const LocalArrayT& B)
{
    int index = 0;
    double sum = 0;
    for (int i=0; i<n_en_displ; i++)
    {
    for (int j=0; j<n_sd; j++)
        sum = sum + u_dotdot(i,j)*u_dotdot(i,j);
    }
    norm = sqrt(sum);
}
*/

/*
<kinematic_BC dof="3" node_ID="6" type="fixed"/>
<kinematic_BC dof="1" node_ID="6" type="fixed"/>
<kinematic_BC dof="2" node_ID="1" type="fixed"/>
<kinematic_BC dof="2" node_ID="4" schedule="1" type="u" value="-0.5"/>
*/
