#include "FSMicromorphic3DCurrConfigT.h"

#include "OutputSetT.h"
#include "ParameterContainerT.h"
#include "CommunicatorT.h"
#include <cmath>

using namespace Tahoe;

/* constructor */

FSMicromorphic3DCurrConfigT::FSMicromorphic3DCurrConfigT(const ElementSupportT& support):
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
    SetName("micromorphic_FS_3D_inc_curr_config");
}

/* destructor */
FSMicromorphic3DCurrConfigT::~FSMicromorphic3DCurrConfigT(void)
{
    delete fShapes_displ;
    delete fShapes_micro;
}


void FSMicromorphic3DCurrConfigT::Echo_Input_Data(void)
{
    cout << "#######################################################" << endl;
    cout << "############### ECHO FSMicromorphic3D DATA #########################" << endl;
    cout << "#######################################################" << endl;

    //################## material parameters ##################
    cout << "iConstitutiveModelType "               << iConstitutiveModelType         << endl;
    //--Plasticity parameters
    cout << "fMaterial_Params[kc0] "               << fMaterial_Params[kc0]          << endl;
    cout << "fMaterial_Params[kHc] "               << fMaterial_Params[kHc]          << endl;
    cout << "fMaterial_Params[kc0_chi] "               << fMaterial_Params[kc0_chi]          << endl;
    cout << "fMaterial_Params[kHc_chi] "               << fMaterial_Params[kHc_chi]          << endl;
//    cout << "fMaterial_Params[kGc0_chi1] "               << fMaterial_Params[kGc0_chi1]          << endl;
//    cout << "fMaterial_Params[kHGc_chi] "               << fMaterial_Params[kHGc_chi]          << endl;
//    cout << "fMaterial_Params[kGc0_chi2] "               << fMaterial_Params[kGc0_chi2]          << endl;
//    cout << "fMaterial_Params[kGc0_chi3] "               << fMaterial_Params[kGc0_chi3]          << endl;
    //cout << "fMaterial_Params[kZ0c] "              << fMaterial_Params[kZ0c]       << endl;
    cout << "fMaterial_Params[kFphi] "             << fMaterial_Params[kFphi]        << endl;
    cout << "fMaterial_Params[kDpsi] "             << fMaterial_Params[kDpsi]        << endl;
    cout << "fMaterial_Params[kFphi_chi] "             << fMaterial_Params[kFphi_chi]        << endl;
    cout << "fMaterial_Params[kDpsi_chi] "             << fMaterial_Params[kDpsi_chi]        << endl;
//    cout << "fMaterial_Params[kFGphi_chi] "             << fMaterial_Params[kFphi_chi]        << endl;
//    cout << "fMaterial_Params[kDGpsi_chi] "             << fMaterial_Params[kDpsi_chi]        << endl;
    //-- Elasticity parameters for solid
    cout << "fMaterial_Params[kMu] "                << fMaterial_Params[kMu]          << endl;
    cout << "fMaterial_Params[kLambda] "            << fMaterial_Params[kLambda] << endl;
    cout << "fMaterial_Params[kNu] "                << fMaterial_Params[kNu]          << endl;
    cout << "fMaterial_Params[kSigma_const] "       << fMaterial_Params[kSigma_const]  << endl;
    cout << "fMaterial_Params[kTau] "               << fMaterial_Params[kTau]          << endl;
    cout << "fMaterial_Params[kEta] "               << fMaterial_Params[kEta]           << endl;
    cout << "fMaterial_Params[kKappa] "             << fMaterial_Params[kKappa]  << endl;
    cout << "fMaterial_Params[kTau1] "              << fMaterial_Params[kTau1]          << endl;
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

void FSMicromorphic3DCurrConfigT::RHSDriver(void)
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

void FSMicromorphic3DCurrConfigT::Equations(AutoArrayT<const iArray2DT*>& eq_d,
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

void FSMicromorphic3DCurrConfigT::LHSDriver(GlobalT::SystemTypeT)
{
/** Everything done in RHSDriver for efficiency */
//cout << "############### In LHS Driver ############### \n";
}

//---------------------------------------------------------------------


void FSMicromorphic3DCurrConfigT::Select_Equations (const int &iBalLinChoice, const int &iBalFirstMomMomChoice )
{
    /** Choices for Linear Momentum Balance Equation */

    switch ( iBalLinChoice )    {

    default :
    cout << "FSMicromorphic3DCurrConfigT::Select_Equations() .. currently only one linear momentum balance for micromorphic continuum \n";
    break;
    }

    /** Choices for First Moment of Momentum Balance Equation */

    switch ( iBalFirstMomMomChoice )    {

    default :
    cout << "FSMicromorphic3DCurrConfigT::Select_Equations() .. currently only one first moment of momentum balance equation for micromorphic continuum \n";
    break;
    }

}

//---------------------------------------------------------------------

/* return true if the element contributes to the solution of the
 * given group. */
bool FSMicromorphic3DCurrConfigT::InGroup(int group) const
{

    return group == fDispl->Group() || group == fMicro->Group();
}

//---------------------------------------------------------------------


/* initialize/finalize step */
void FSMicromorphic3DCurrConfigT::InitStep(void)
{

    /* inherited */
    ElementBaseT::InitStep();
}


/* close current time increment */
void FSMicromorphic3DCurrConfigT::CloseStep(void)
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
   // ElementSupport().ResetAverage(knumstrain+knumstress+knum_d_state+knumdispl);
    ElementSupport().ResetAverage(knumstrain+knumstress+knum_d_state);
    dArray2DT out_variable_all;
    dArrayT out_variable;
    dArray2DT nd_var(NumElementNodes(), knumstrain+knumstress+knum_d_state);
    Top();
    while (NextElement())
    {
        /* extrapolate */
        nd_var = 0.0;
        out_variable_all.Alias(fNumIP_displ, knumstrain+knumstress+knum_d_state, fIPVariable(CurrElementNumber()));
        fShapes_displ->TopIP();
        while (fShapes_displ->NextIP())
        {
        out_variable.Alias(knumstrain+knumstress+knum_d_state, out_variable_all(fShapes_displ->CurrIP()));
        fShapes_displ->Extrapolate(out_variable, nd_var);
        }

        /* accumulate - extrapolation done from ip's to corners => X nodes  */
        ElementSupport().AssembleAverage(CurrentElement().NodesX(), nd_var);
    }

    /* get nodally averaged values */
    dArray2DT extrap_values;
    ElementSupport().OutputUsedAverage(extrap_values);

    /* temp space for group displacements */
    int num_node_output = fDispl->NumDOF() + fMicro->NumDOF() + knumstrain + knumstress + knum_d_state ;
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
        for (int j = 0; j < (knumstrain+knumstress+knum_d_state); j++)
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

  //if(iConstitutiveModelType==2)
  //{
    SigN_IPs_el_n      = SigN_IPs_el;
    GammaN_IPs_el_n    = GammaN_IPs_el;
    mn_IPs_el_n        = mn_IPs_el;
    sn_sigman_IPs_el_n = sn_sigman_IPs_el;

    Fn_ar_IPs_el=F_ar_IPs_el;
    FnInv_ar_IPs_el=FInv_ar_IPs_el;
    ChiN_ar_IPs_el_n=Chi_ar_IPs_el;
    GRAD_ChiN_ar_IPs_el_n=GRAD_Chi_ar_IPs_el;
   // }
  if(iConstitutiveModelType==3)
  {
	    fState_variables_n_Elements_IPs = fState_variables_Elements_IPs;
	    fdGdCauchy_Stress_n_Elements_IPs = fdGdCauchy_Stress_Elements_IPs;
	    fCauchy_stress_Elements_n_IPs = fCauchy_stress_Elements_IPs;
	    fDeformation_Gradient_n_Elements_IPs = fDeformation_Gradient_Elements_IPs;
  }


    step_number = ElementSupport().StepNumber();
    fs_micromorph3D_out << "xxxxxxxxxxxxxx " << endl;
    fs_micromorph3D_out << "step number" << ":" << step_number << endl;
    fs_micromorph3D_out << "xxxxxxxxxxxxxx " << endl;

}








/* resets to the last converged solution */
/*
GlobalT::RelaxCodeT FSMicromorphic3DCurrConfigT::ResetStep(void)
{
    const char caller[] = "FSMicromorphic3DCurrConfigT::ResetStep";

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
GlobalT::RelaxCodeT FSMicromorphic3DCurrConfigT::RelaxSystem(void)
{
    const char caller[] = "FSMicromorphic3DCurrConfigT::RelaxSystem";

    // inherited
    GlobalT::RelaxCodeT relax = ElementBaseT::RelaxSystem();

    // loop over materials
    //needs to be implemented
#pragma message("relax step for materials not implemented")
    //ExceptionT::GeneralFail(caller, "relax step for materials not implemented");

    return relax;
}
*/


void FSMicromorphic3DCurrConfigT::SendOutput(int kincode)
{
#pragma unused(kincode)
//not implemented
}


/* return geometry and number of nodes on each facet */
void FSMicromorphic3DCurrConfigT::FacetGeometry(ArrayT<GeometryT::CodeT>& facet_geometry,
    iArrayT& num_facet_nodes) const
{
    /* from integration domain */
    ShapeFunctionDispl().FacetGeometry(facet_geometry, num_facet_nodes);
}


/* form of tangent matrix */
GlobalT::SystemTypeT FSMicromorphic3DCurrConfigT::TangentType(void) const
{
    return GlobalT::kNonSymmetric;
}

/*
void FSMicromorphic3DCurrConfigT::SetStatus(const ArrayT<ElementCardT::StatusT>& status)
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
void FSMicromorphic3DCurrConfigT::InitialCondition(void)
{
    //           cout<<"CHECK POINT-12"<<endl;
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
void FSMicromorphic3DCurrConfigT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
    const char caller[] = "FSMicromorphic3DCurrConfigT::AddNodalForce";
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
    global_iteration = IterationNumber();

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
            //fShapes_displ->SetDerivatives_DN_DDN(); Commented out for Q8P8
            fShapes_displ->SetDerivatives();

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
                 fShapes_displ->TopIP();
                // residual for displacement field
                //generate this vector fFd_int
                //fFd_int=0.0;
                //Vint_1_temp=0.0;
                //Vint_1=0.0;

                while (fShapes_displ->NextIP())
                {
                    fFd_int=0.0;

                }

            }
            else /* pressure nodal force */
            {
                /* residual for micro-displacement-gradient field */
                // generate this vector fFphi_int
                fShapes_displ->TopIP();
             //   fFphi_int=0.0;
              //  Vint_2_temp=0.0;
               // Vint_2=0.0;
               // Vint_3_temp=0.0;
               // Vint_3=0.0;
                while (fShapes_displ->NextIP())
                {
                    //nothing right now
                    fFphi_int =0.0;
                /*  double scale;
                    double scale_const = (*Weight++)*(*Det++);
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
                Vint_3+=Vint_3_temp;*/

                }
               //  fFphi_int  = Vint_2;
               // fFphi_int +=Vint_3;
               // fFphi_int *=-1;

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

double FSMicromorphic3DCurrConfigT::InternalEnergy ( void )
{
//not implemented

    return 0.0;
}

//---------------------------------------------------------------------

/* write restart data to the output stream */
void FSMicromorphic3DCurrConfigT::WriteRestart(ostream& out) const
{

    /* inherited */
    ElementBaseT::WriteRestart(out);

    /* write state variable data */
    out << fdState;
}

//---------------------------------------------------------------------

/* read restart data to the output stream */
void FSMicromorphic3DCurrConfigT::ReadRestart(istream& in)
{

    /* inherited */
    ElementBaseT::ReadRestart(in);

    /* write state variable data */
    in >> fdState;
}

//---------------------------------------------------------------------

void FSMicromorphic3DCurrConfigT::RegisterOutput(void)
{

    /* collect block ID's */
    ArrayT<StringT> block_ID(fBlockData.Length());
    for (int i = 0; i < block_ID.Length(); i++)
    block_ID[i] = fBlockData[i].ID();

    /* output per element - strain, stress, and ISVs at the integration points */
    ArrayT<StringT> e_labels(fNumIP_displ*(knumstrain+knumstress+knum_d_state));

    /* over integration points */
    // enter what values you need at integration points
    // stress and strain
  //  const char* slabels3D[] = {"s11", "s22", "s33","s23","s13","s12","e11","e22","e33","e23","e13","e12"};
    const char* slabels3D[] = {"s11","s12","s13","s21","s22","s23","s31","s32","s33","e11","e22","e33","e12","e13","e21","e23","e31","e32"};


  //  if(iConstitutiveModelType==3)
   //  {
    const char* svlabels3D[] = {"kc","khc","kc_chi","khc_chi","kDelgamma","kDelgammachi","trCauchy_Stress","||dev(Cauchy_Stress)||","trRel","||dev(Rel)||","||trm||","||dev(m)||","trSPK","||dev(SPK)||","trSIGMA_S","||dev(SIGMA_S)||","||trM||","||dev(M)||","||PHI||","||GPHI||","treps","deveps","invtrgammastn","invdevgammastn"};
//          ,"kF11","kF12","kF13","kF21","kF22","kF23","kF31","kF32","kF33","kFe11","kFe12","kFe13","kFe21","kFe22","kFe23","kFe31","kFe32","kFe33",
//          "kX11","kX12","kX13","kX21","kX22","kX23","kX31","kX32","kX33","kXe11","kXe12","kXe13","kXe21","kXe22","kXe23","kXe31","kXe32","kXe33"};
  //  }
  // else
   // {

    // state variables; ?
 //   const char* svlabels3D[] = {"||devs||","||devrel||","||devmklm||","tr(sigma)","tr(s_sigma)","trmklm","E11","E22","E33","E12","E13","E21","E23","E31","E32","VE11","VE22","VE33","VE12","VE13","VE21","VE23","VE31","VE32"
 //   ,"F11","F22","F33","F12","F13","F21","F23","F31","F32","GAMMA(1,1,1)","GAMMA(2,1,2)","GAMMA(3,1,3)","GAMMA(1,2,1)","GAMMA(2,2,2)","GAMMA(3,2,3)",
  //  "GAMMA(1,3,1)","GAMMA(2,3,2)","GAMMA(3,3,3)"};
    //}
    int count = 0;
    for (int j = 0; j < fNumIP_micro; j++)
    {
        StringT ip_label;
        ip_label.Append("ip", j+1);

        /* over strain and stress components */
        for (int i = 0; i < knumstrain+knumstress; i++)
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
    int num_node_output = fDispl->NumDOF() + fMicro->NumDOF() + knumstrain + knumstress + knum_d_state;
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
    for (int i = 0; i < knumstrain+knumstress; i++)
    n_labels[count++] = slabels3D[i];

    /* labels from state variables at the nodes */
    for (int i = 0; i < knum_d_state; i++)
    n_labels[count++] = svlabels3D[i];

    /* set output specifier */
#pragma message("FSMicromorphic3DCurrConfigT::RegisterOutput: is this right? ")
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

void FSMicromorphic3DCurrConfigT::WriteOutput(void)
{

    bStep_Complete=1;
    RHSDriver();
    bStep_Complete=0;

    /* my output set */
    const OutputSetT& output_set = ElementSupport().OutputSet(fOutputID);

    /* my nodes used */
    const iArrayT& nodes_used = output_set.NodesUsed();

    /* smooth stresses to nodes */
    ElementSupport().ResetAverage(knumstrain+knumstress+knum_d_state);
    dArray2DT out_variable_all;
    dArrayT out_variable;
    dArray2DT nd_var(NumElementNodes(), knumstrain+knumstress+knum_d_state);
    Top();
    while (NextElement())
    {
        /* extrapolate */
        nd_var = 0.0;
        out_variable_all.Alias(fNumIP_displ, knumstrain+knumstress+knum_d_state, fIPVariable(CurrElementNumber()));
        fShapes_displ->TopIP();
        while (fShapes_displ->NextIP())
        {
            out_variable.Alias(knumstrain+knumstress+knum_d_state, out_variable_all(fShapes_displ->CurrIP()));
            fShapes_displ->Extrapolate(out_variable, nd_var);
        }

        /* accumulate - extrapolation done from ip's to corners => X nodes  */
        ElementSupport().AssembleAverage(CurrentElement().NodesX(), nd_var);
    }

    /* get nodally averaged values */
    dArray2DT extrap_values;
    ElementSupport().OutputUsedAverage(extrap_values);

    /* temp space for group displacements */
    int num_node_output = fDispl->NumDOF() + fMicro->NumDOF() + knumstrain + knumstress + knum_d_state;
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
        for (int j = 0; j < (knumstrain+knumstress+knum_d_state); j++)
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
void FSMicromorphic3DCurrConfigT::RHSDriver_staggered(void)
{

    const char caller[] = "FSMicromorphic3DCurrConfigT::RHSDriver_staggered";
#pragma message("staggered solution not implemented")
}

//---------------------------------------------------------------------
/* form group contribution to the stiffness matrix and RHS */
void FSMicromorphic3DCurrConfigT::RHSDriver_monolithic(void)
{

    const char caller[] = "FSMicromorphic3DCurrConfigT::RHSDriver_monolithic";
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
    global_iteration = IterationNumber();

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

//////////////////////// Initialization of Current configuration plasticity matrices////////////////////////////
    fV1p = 0.0;
    fV1e = 0.0;
    fdGdCauchy_Stress = 0.0;
    fdGdCauchy_Stress_tr = 0.0;
    fdGdCauchy_Stress_tr_Transpose = 0.0;
    fdGdCauchy_Stress_tr_trace = 0.0;
    fCauchy_stress_tensor_current_IP_tr = 0.0;
    mean_Cauchy_stress_tr = 0.0;
    dev_Cauchy_stress_tr = 0.0;
    fNorm_dev_Cauchy_stress_tensor_current_IP_tr = 0.0;
    Vintp_1_temp = 0.0;
    Vintp_1 = 0.0;
    Vinte_1_temp = 0.0;
    Vinte_1 = 0.0;
    fTemp_matrix_one_x_one = 0.0;
    fTemp_matrix_nsd_x_nsd = 0.0;
    fTemp_matrix_nsd_x_nsd2 = 0.0;
    fTemp_matrix_nsd_x_nsd3 = 0.0;
    fVelocity_Gradient_current_IP = 0.0;
    Comp11 = 0.0;
    Comp22 = 0.0;
    Comp33 = 0.0;
    Comp44 = 0.0;
    fCauchy_stress_tensor_current_IP = 0.0;
    fNorm_dev_Cauchy_stress_tensor_current_IP = 0.0;
    fdev_Cauchy_stress_tensor_current_IP = 0.0;
    fElastic_Velocity_Gradient_current_IP = 0.0;


    IJe_1 = 0.0;
    IJe_2 = 0.0;
    IJe_3 = 0.0;
    IJe_4 = 0.0;
    IJe_5 = 0.0;
    IJe_6 = 0.0;
    IJe_7 = 0.0;
    IJe_8 = 0.0;


    I1p_1 = 0.0;
    I1p_2 = 0.0;
    I1p_3 = 0.0;
    I1p_4 = 0.0;
    I1p_5 = 0.0;
    I1p_6 = 0.0;


    I2p_1 = 0.0;
    I2p_2 = 0.0;
    I2p_3 = 0.0;
    I2p_4 = 0.0;
    I2p_5 = 0.0;
    I2p_6 = 0.0;


    I3p_1 = 0.0;
    I3p_2 = 0.0;
    I3p_3 = 0.0;
    I3p_4 = 0.0;
    I3p_5 = 0.0;
    I3p_6 = 0.0;


    I4p_1 = 0.0;
    I4p_2 = 0.0;
    I4p_3 = 0.0;
    I4p_4 = 0.0;
    I4p_5 = 0.0;
    I4p_6 = 0.0;


    I5p_1 = 0.0;
    I5p_2 = 0.0;
    I5p_3 = 0.0;
    I5p_4 = 0.0;
    I5p_5 = 0.0;
    I5p_6 = 0.0;



    I6p_1 = 0.0;
    I6p_2 = 0.0;
    I6p_3 = 0.0;
    I6p_4 = 0.0;
    I6p_5 = 0.0;
    I6p_6 = 0.0;

    I7p_1 = 0.0;
    I7p_2 = 0.0;
    I7p_3 = 0.0;
    I7p_4 = 0.0;
    I7p_5 = 0.0;
    I7p_6 = 0.0;
    I7p_7 = 0.0;
    I7p_8 = 0.0;
    I7p_9 = 0.0;
    I7p_10 = 0.0;
    I7p_11 = 0.0;
    I7p_12 = 0.0;
    I7p_13 = 0.0;


    I8p_1 = 0.0;
    I8p_2 = 0.0;
    I8p_3 = 0.0;
    I8p_4 = 0.0;
    I8p_5 = 0.0;
    I8p_6 = 0.0;
    I8p_7 = 0.0;
    I8p_8 = 0.0;
    I8p_9 = 0.0;
    I8p_10 = 0.0;
    I8p_11 = 0.0;
    I8p_12 = 0.0;
    I8p_13 = 0.0;


    I9p_1 = 0.0;
    I9p_2 = 0.0;
    I9p_3 = 0.0;
    I9p_4 = 0.0;
    I9p_5 = 0.0;
    I9p_6 = 0.0;
    I9p_7 = 0.0;
    I9p_8 = 0.0;
    I9p_9 = 0.0;
    I9p_10 = 0.0;
    I9p_11 = 0.0;
    I9p_12 = 0.0;
    I9p_13 = 0.0;


    I10p_1 = 0.0;
    I10p_2 = 0.0;
    I10p_3 = 0.0;
    I10p_4 = 0.0;
    I10p_5 = 0.0;
    I10p_6 = 0.0;
    I10p_7 = 0.0;
    I10p_8 = 0.0;
    I10p_9 = 0.0;
    I10p_10 = 0.0;
    I10p_11 = 0.0;
    I10p_12 = 0.0;
    I10p_13 = 0.0;


    I_temp_11p_1 = 0.0;
    I_temp_11p_2 = 0.0;
    I_temp_11p_3 = 0.0;
    I_temp_11p_4 = 0.0;
    I_temp_11p_5 = 0.0;
    I_temp_11p_6 = 0.0;
    I_temp_11p_7 = 0.0;
    I_temp_11p_8 = 0.0;
    I_temp_11p_9 = 0.0;
    I_temp_11p_10 = 0.0;
    I_temp_11p_11 = 0.0;
    I_temp_11p_12 = 0.0;
    I_temp_11p_13 = 0.0;
    I_temp_11p_4_Transpose = 0.0;


    I11p_1 = 0.0;
    I11p_2 = 0.0;

    I12p_1 = 0.0;
    I12p_2 = 0.0;
    I12p_3 = 0.0;
    I12p_4 = 0.0;
    I12p_5 = 0.0;
    I12p_6 = 0.0;

    I13p_1 = 0.0;
    I13p_2 = 0.0;
    I13p_3 = 0.0;
    I13p_4 = 0.0;
    I13p_5 = 0.0;
    I13p_6 = 0.0;
////////////////////////////Test///////////////////////
    I12p_test_2 = 0.0;
    I12p_test_3 = 0.0;
    I12p_test_4 = 0.0;
    I_temp_11p_test_1 = 0.0;
    I_temp_11p_test_2 = 0.0;
    I_temp_11p_test_3 = 0.0;
    I_temp_11p_test_4 = 0.0;
    I_temp_11p_20 = 0.0;
    I_temp_11p_21 = 0.0;
    I_temp_11p_22 = 0.0;
    I_temp_11p_23 = 0.0;
    I_temp_11p_24 = 0.0;

   	II_temp_11p_1_1 = 0.0;
    II_temp_11p_1_2 = 0.0;
    II_temp_11p_1_3 = 0.0;
    II_temp_11p_1_4 = 0.0;
    II_temp_11p_1_5 = 0.0;
    II_temp_11p_1_6 = 0.0;

   	II_temp_11p_2_1 = 0.0;
    II_temp_11p_2_2 = 0.0;
    II_temp_11p_2_3 = 0.0;
    II_temp_11p_2_4 = 0.0;
    II_temp_11p_2_5 = 0.0;
    II_temp_11p_2_6 = 0.0;

   	II_temp_11p_3_1 = 0.0;
    II_temp_11p_3_2 = 0.0;
    II_temp_11p_3_3 = 0.0;
    II_temp_11p_3_4 = 0.0;
    II_temp_11p_3_5 = 0.0;
    II_temp_11p_3_6 = 0.0;

   	II_temp_11p_4_1 = 0.0;
    II_temp_11p_4_2 = 0.0;
    II_temp_11p_4_3 = 0.0;
    II_temp_11p_4_4 = 0.0;
    II_temp_11p_4_5 = 0.0;
    II_temp_11p_4_6 = 0.0;

   	II_temp_11p_5_1 = 0.0;
    II_temp_11p_5_2 = 0.0;
    II_temp_11p_5_3 = 0.0;
    II_temp_11p_5_4 = 0.0;
    II_temp_11p_5_5 = 0.0;
    II_temp_11p_5_6 = 0.0;

   	II_temp_11p_6_1 = 0.0;
    II_temp_11p_6_2 = 0.0;
    II_temp_11p_6_3 = 0.0;
    II_temp_11p_6_4 = 0.0;
    II_temp_11p_6_5 = 0.0;
    II_temp_11p_6_6 = 0.0;

   	II_temp_11p_7_1 = 0.0;
    II_temp_11p_7_2 = 0.0;
    II_temp_11p_7_3 = 0.0;
    II_temp_11p_7_4 = 0.0;
    II_temp_11p_7_5 = 0.0;
    II_temp_11p_7_6 = 0.0;

   	II_temp_11p_8_1 = 0.0;
    II_temp_11p_8_2 = 0.0;
    II_temp_11p_8_3 = 0.0;
    II_temp_11p_8_4 = 0.0;
    II_temp_11p_8_5 = 0.0;
    II_temp_11p_8_6 = 0.0;

   	II_temp_11p_9_1 = 0.0;
    II_temp_11p_9_2 = 0.0;
    II_temp_11p_9_3 = 0.0;
    II_temp_11p_9_4 = 0.0;
    II_temp_11p_9_5 = 0.0;
    II_temp_11p_9_6 = 0.0;

   	II_temp_11p_10_1 = 0.0;
    II_temp_11p_10_2 = 0.0;
    II_temp_11p_10_3 = 0.0;
    II_temp_11p_10_4 = 0.0;
    II_temp_11p_10_5 = 0.0;
    II_temp_11p_10_6 = 0.0;

   	II_temp_11p_11_1 = 0.0;
    II_temp_11p_11_2 = 0.0;
    II_temp_11p_11_3 = 0.0;
    II_temp_11p_11_4 = 0.0;
    II_temp_11p_11_5 = 0.0;
    II_temp_11p_11_6 = 0.0;

   	II_temp_11p_12_1 = 0.0;
    II_temp_11p_12_2 = 0.0;
    II_temp_11p_12_3 = 0.0;
    II_temp_11p_12_4 = 0.0;
    II_temp_11p_12_5 = 0.0;
    II_temp_11p_12_6 = 0.0;

   	II_temp_11p_13_1 = 0.0;
    II_temp_11p_13_2 = 0.0;
    II_temp_11p_13_3 = 0.0;
    II_temp_11p_13_4 = 0.0;
    II_temp_11p_13_5 = 0.0;
    II_temp_11p_13_6 = 0.0;



    fKu_I_temp_11p_1_1 = 0.0;
    fKu_I_temp_11p_1_2 = 0.0;
    fKu_I_temp_11p_1_3 = 0.0;
    fKu_I_temp_11p_1_4 = 0.0;
    fKu_I_temp_11p_1_5 = 0.0;
    fKu_I_temp_11p_1_6 = 0.0;

    fKu_I_temp_11p_2_1 = 0.0;
    fKu_I_temp_11p_2_2 = 0.0;
    fKu_I_temp_11p_2_3 = 0.0;
    fKu_I_temp_11p_2_4 = 0.0;
    fKu_I_temp_11p_2_5 = 0.0;
    fKu_I_temp_11p_2_6 = 0.0;

    fKu_I_temp_11p_3_1 = 0.0;
    fKu_I_temp_11p_3_2 = 0.0;
    fKu_I_temp_11p_3_3 = 0.0;
    fKu_I_temp_11p_3_4 = 0.0;
    fKu_I_temp_11p_3_5 = 0.0;
    fKu_I_temp_11p_3_6 = 0.0;

    fKu_I_temp_11p_4_1 = 0.0;
    fKu_I_temp_11p_4_2 = 0.0;
    fKu_I_temp_11p_4_3 = 0.0;
    fKu_I_temp_11p_4_4 = 0.0;
    fKu_I_temp_11p_4_5 = 0.0;
    fKu_I_temp_11p_4_6 = 0.0;

    fKu_I_temp_11p_5_1 = 0.0;
    fKu_I_temp_11p_5_2 = 0.0;
    fKu_I_temp_11p_5_3 = 0.0;
    fKu_I_temp_11p_5_4 = 0.0;
    fKu_I_temp_11p_5_5 = 0.0;
    fKu_I_temp_11p_5_6 = 0.0;

    fKu_I_temp_11p_6_1 = 0.0;
    fKu_I_temp_11p_6_2 = 0.0;
    fKu_I_temp_11p_6_3 = 0.0;
    fKu_I_temp_11p_6_4 = 0.0;
    fKu_I_temp_11p_6_5 = 0.0;
    fKu_I_temp_11p_6_6 = 0.0;

    fKu_I_temp_11p_7_1 = 0.0;
    fKu_I_temp_11p_7_2 = 0.0;
    fKu_I_temp_11p_7_3 = 0.0;
    fKu_I_temp_11p_7_4 = 0.0;
    fKu_I_temp_11p_7_5 = 0.0;
    fKu_I_temp_11p_7_6 = 0.0;

    fKu_I_temp_11p_8_1 = 0.0;
    fKu_I_temp_11p_8_2 = 0.0;
    fKu_I_temp_11p_8_3 = 0.0;
    fKu_I_temp_11p_8_4 = 0.0;
    fKu_I_temp_11p_8_5 = 0.0;
    fKu_I_temp_11p_8_6 = 0.0;

    fKu_I_temp_11p_9_1 = 0.0;
    fKu_I_temp_11p_9_2 = 0.0;
    fKu_I_temp_11p_9_3 = 0.0;
    fKu_I_temp_11p_9_4 = 0.0;
    fKu_I_temp_11p_9_5 = 0.0;
    fKu_I_temp_11p_9_6 = 0.0;

    fKu_I_temp_11p_10_1 = 0.0;
    fKu_I_temp_11p_10_2 = 0.0;
    fKu_I_temp_11p_10_3 = 0.0;
    fKu_I_temp_11p_10_4 = 0.0;
    fKu_I_temp_11p_10_5 = 0.0;
    fKu_I_temp_11p_10_6 = 0.0;

    fKu_I_temp_11p_11_1 = 0.0;
    fKu_I_temp_11p_11_2 = 0.0;
    fKu_I_temp_11p_11_3 = 0.0;
    fKu_I_temp_11p_11_4 = 0.0;
    fKu_I_temp_11p_11_5 = 0.0;
    fKu_I_temp_11p_11_6 = 0.0;

    fKu_I_temp_11p_12_1 = 0.0;
    fKu_I_temp_11p_12_2 = 0.0;
    fKu_I_temp_11p_12_3 = 0.0;
    fKu_I_temp_11p_12_4 = 0.0;
    fKu_I_temp_11p_12_5 = 0.0;
    fKu_I_temp_11p_12_6 = 0.0;

    fKu_I_temp_11p_13_1 = 0.0;
    fKu_I_temp_11p_13_2 = 0.0;
    fKu_I_temp_11p_13_3 = 0.0;
    fKu_I_temp_11p_13_4 = 0.0;
    fKu_I_temp_11p_13_5 = 0.0;
    fKu_I_temp_11p_13_6 = 0.0;


    //////////////////////some other comparison which is made/////////////////
	I_temp_11p_test_2_1_1 = 0.0;
	I_temp_11p_test_2_1_2 = 0.0;
	I_temp_11p_test_2_1_3 = 0.0;
	I_temp_11p_test_2_1_4 = 0.0;

	I_temp_11p_test_2_2_1 = 0.0;
	I_temp_11p_test_2_2_2 = 0.0;
   	I_temp_11p_test_2_2_3 = 0.0;
	I_temp_11p_test_2_2_4 = 0.0;

	I_temp_11p_test_2_3_1 = 0.0;
	I_temp_11p_test_2_3_2 = 0.0;


	I_temp_11p_test_2_4_1 = 0.0;
	I_temp_11p_test_2_4_2 = 0.0;
	I_temp_11p_test_2_4_3 = 0.0;
	I_temp_11p_test_2_4_4 = 0.0;

	I_temp_11p_test_2_5_1 = 0.0;
	I_temp_11p_test_2_5_2 = 0.0;


	I_temp_11p_test_2_6_1 = 0.0;
	I_temp_11p_test_2_6_2 = 0.0;
	I_temp_11p_test_2_6_3 = 0.0;
	I_temp_11p_test_2_6_4 = 0.0;

	I_temp_11p_test_2_7_1 = 0.0;
	I_temp_11p_test_2_7_2 = 0.0;


	I_temp_11p_test_2_8_1 = 0.0;
	I_temp_11p_test_2_8_2 = 0.0;
	I_temp_11p_test_2_8_3 = 0.0;
	I_temp_11p_test_2_8_4 = 0.0;

	I_temp_11p_test_2_9_1 = 0.0;
	I_temp_11p_test_2_9_2 = 0.0;
	I_temp_11p_test_2_9_3 = 0.0;
	I_temp_11p_test_2_9_4 = 0.0;

	I_temp_11p_test_2_10_1 = 0.0;
	I_temp_11p_test_2_10_2 = 0.0;
	I_temp_11p_test_2_10_3 = 0.0;
	I_temp_11p_test_2_10_4 = 0.0;

	I_temp_11p_test_2_11_1 = 0.0;
	I_temp_11p_test_2_11_2 = 0.0;
	I_temp_11p_test_2_11_3 = 0.0;
	I_temp_11p_test_2_11_4 = 0.0;

	I_temp_11p_test_2_12_1 = 0.0;
	I_temp_11p_test_2_12_2 = 0.0;
	I_temp_11p_test_2_12_3 = 0.0;
	I_temp_11p_test_2_12_4 = 0.0;


	I_temp_11p_test_2_13_1 = 0.0;
	I_temp_11p_test_2_13_2 = 0.0;
	I_temp_11p_test_2_13_3 = 0.0;
	I_temp_11p_test_2_13_4 = 0.0;

    fKu_I_temp_11p_test_2_1_1 = 0.0;
    fKu_I_temp_11p_test_2_1_2 = 0.0;
    fKu_I_temp_11p_test_2_1_3 = 0.0;
    fKu_I_temp_11p_test_2_1_4 = 0.0;

    fKu_I_temp_11p_test_2_2_1 = 0.0;
    fKu_I_temp_11p_test_2_2_2 = 0.0;
    fKu_I_temp_11p_test_2_2_3 = 0.0;
    fKu_I_temp_11p_test_2_2_4 = 0.0;

    fKu_I_temp_11p_test_2_3_1 = 0.0;
    fKu_I_temp_11p_test_2_3_2 = 0.0;


    fKu_I_temp_11p_test_2_4_1 = 0.0;
    fKu_I_temp_11p_test_2_4_2 = 0.0;
    fKu_I_temp_11p_test_2_4_3 = 0.0;
    fKu_I_temp_11p_test_2_4_4 = 0.0;

    fKu_I_temp_11p_test_2_5_1 = 0.0;
    fKu_I_temp_11p_test_2_5_2 = 0.0;


    fKu_I_temp_11p_test_2_6_1 = 0.0;
    fKu_I_temp_11p_test_2_6_2 = 0.0;
    fKu_I_temp_11p_test_2_6_3 = 0.0;
    fKu_I_temp_11p_test_2_6_4 = 0.0;

    fKu_I_temp_11p_test_2_7_1 = 0.0;
    fKu_I_temp_11p_test_2_7_2 = 0.0;


    fKu_I_temp_11p_test_2_8_1 = 0.0;
    fKu_I_temp_11p_test_2_8_2 = 0.0;
    fKu_I_temp_11p_test_2_8_3 = 0.0;
    fKu_I_temp_11p_test_2_8_4 = 0.0;


    fKu_I_temp_11p_test_2_9_1 = 0.0;
    fKu_I_temp_11p_test_2_9_2 = 0.0;
    fKu_I_temp_11p_test_2_9_3 = 0.0;
    fKu_I_temp_11p_test_2_9_4 = 0.0;

    fKu_I_temp_11p_test_2_10_1 = 0.0;
    fKu_I_temp_11p_test_2_10_2 = 0.0;
    fKu_I_temp_11p_test_2_10_3 = 0.0;
    fKu_I_temp_11p_test_2_10_4 = 0.0;

    fKu_I_temp_11p_test_2_11_1 = 0.0;
    fKu_I_temp_11p_test_2_11_2 = 0.0;
    fKu_I_temp_11p_test_2_11_3 = 0.0;
    fKu_I_temp_11p_test_2_11_4 = 0.0;

    fKu_I_temp_11p_test_2_12_1 = 0.0;
    fKu_I_temp_11p_test_2_12_2 = 0.0;
    fKu_I_temp_11p_test_2_12_3 = 0.0;
    fKu_I_temp_11p_test_2_12_4 = 0.0;

    fKu_I_temp_11p_test_2_13_1 = 0.0;
    fKu_I_temp_11p_test_2_13_2 = 0.0;
    fKu_I_temp_11p_test_2_13_3 = 0.0;
    fKu_I_temp_11p_test_2_13_4 = 0.0;


    fKdd_previous = 0.0;
    fKdd_full_implemet = 0.0;
    difference = 0.0;
    /////////////////////////////////////////////////////////


    fKu_IJe_1 = 0.0;
    fKu_IJe_2 = 0.0;
    fKu_IJe_3 = 0.0;
    fKu_IJe_4 = 0.0;
    fKu_IJe_5 = 0.0;
    fKu_IJe_6 = 0.0;
    fKu_IJe_7 = 0.0;
    fKu_IJe_8 = 0.0;


    fKu_I1p_1 = 0.0;
    fKu_I1p_2 = 0.0;
    fKu_I1p_3 = 0.0;
    fKu_I1p_4 = 0.0;
    fKu_I1p_5 = 0.0;
    fKu_I1p_6 = 0.0;



    fKu_I2p_1 = 0.0;
    fKu_I2p_2 = 0.0;
    fKu_I2p_3 = 0.0;
    fKu_I2p_4 = 0.0;
    fKu_I2p_5 = 0.0;
    fKu_I2p_6 = 0.0;



    fKu_I3p_1 = 0.0;
    fKu_I3p_2 = 0.0;
    fKu_I3p_3 = 0.0;
    fKu_I3p_4 = 0.0;
    fKu_I3p_5 = 0.0;
    fKu_I3p_6 = 0.0;




    fKu_I4p_1 = 0.0;
    fKu_I4p_2 = 0.0;
    fKu_I4p_3 = 0.0;
    fKu_I4p_4 = 0.0;
    fKu_I4p_5 = 0.0;
    fKu_I4p_6 = 0.0;



    fKu_I5p_1 = 0.0;
    fKu_I5p_2 = 0.0;
    fKu_I5p_3 = 0.0;
    fKu_I5p_4 = 0.0;
    fKu_I5p_5 = 0.0;
    fKu_I5p_6 = 0.0;



    fKu_I6p_1 = 0.0;
    fKu_I6p_2 = 0.0;
    fKu_I6p_3 = 0.0;
    fKu_I6p_4 = 0.0;
    fKu_I6p_5 = 0.0;
    fKu_I6p_6 = 0.0;

    fKu_I7p_1 = 0.0;
    fKu_I7p_2 = 0.0;
    fKu_I7p_3 = 0.0;
    fKu_I7p_4 = 0.0;
    fKu_I7p_5 = 0.0;
    fKu_I7p_6 = 0.0;
    fKu_I7p_7 = 0.0;
    fKu_I7p_8 = 0.0;
    fKu_I7p_9 = 0.0;
    fKu_I7p_10 = 0.0;
    fKu_I7p_11 = 0.0;
    fKu_I7p_12 = 0.0;
    fKu_I7p_13 = 0.0;


    fKu_I8p_1 = 0.0;
    fKu_I8p_2 = 0.0;
    fKu_I8p_3 = 0.0;
    fKu_I8p_4 = 0.0;
    fKu_I8p_5 = 0.0;
    fKu_I8p_6 = 0.0;
    fKu_I8p_7 = 0.0;
    fKu_I8p_8 = 0.0;
    fKu_I8p_9 = 0.0;
    fKu_I8p_10 = 0.0;
    fKu_I8p_11 = 0.0;
    fKu_I8p_12 = 0.0;
    fKu_I8p_13 = 0.0;


    fKu_I9p_1 = 0.0;
    fKu_I9p_2 = 0.0;
    fKu_I9p_3 = 0.0;
    fKu_I9p_4 = 0.0;
    fKu_I9p_5 = 0.0;
    fKu_I9p_6 = 0.0;
    fKu_I9p_7 = 0.0;
    fKu_I9p_8 = 0.0;
    fKu_I9p_9 = 0.0;
    fKu_I9p_10 = 0.0;
    fKu_I9p_11 = 0.0;
    fKu_I9p_12 = 0.0;
    fKu_I9p_13 = 0.0;


    fKu_I10p_1 = 0.0;
    fKu_I10p_2 = 0.0;
    fKu_I10p_3 = 0.0;
    fKu_I10p_4 = 0.0;
    fKu_I10p_5 = 0.0;
    fKu_I10p_6 = 0.0;
    fKu_I10p_7 = 0.0;
    fKu_I10p_8 = 0.0;
    fKu_I10p_9 = 0.0;
    fKu_I10p_10 = 0.0;
    fKu_I10p_11 = 0.0;
    fKu_I10p_12 = 0.0;
    fKu_I10p_13 = 0.0;


    fKu_I12p_1 = 0.0;
    fKu_I12p_2 = 0.0;
    fKu_I12p_3 = 0.0;
    fKu_I12p_4 = 0.0;
    fKu_I12p_5 = 0.0;
    fKu_I12p_6 = 0.0;

    fKu_I13p_1 = 0.0;
    fKu_I13p_2 = 0.0;
    fKu_I13p_3 = 0.0;
    fKu_I13p_4 = 0.0;
    fKu_I13p_5 = 0.0;
    fKu_I13p_6 = 0.0;


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



     /* Initializing Plasticity Matrices */
     /* All the matrices used in summation have to be initialized to 0.0; */



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


    del_u.DiffOf (u, u_n);
    del_Phi.DiffOf (Phi, Phi_n);

    // calculate derivatives based on reference coordinates
    fInitCoords_displ.SetLocal(fElementCards_displ[e].NodesX());
    fCurrCoords_displ=fInitCoords_displ;
    //fCurrCoords_displ.SetToCombination (1.0, fInitCoords_displ, 1.0, u);
    // fShapes_displ->SetDerivatives_DN_DDN(); Commented out for Q8P8
    fShapes_displ->SetDerivatives();
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
        out_variable_all.Alias(fNumIP_micro, knumstrain+knumstress+knum_d_state, fIPVariable(CurrElementNumber()));
        for (l=0; l < fNumIP_displ; l++)
        //for (l=0; l < fNumIP_micro; l++)
        {
          out_variable.Alias(knumstrain+knumstress+knum_d_state, out_variable_all(l));

          Put_values_In_dArrayT_vector(fCauchy_stress_Elements_IPs, e,l,fTemp_nine_values);
//            Put_values_In_dArrayT_vector(fCauchy_stress_Elements_IPs, e,l,fTemp_six_values);

          out_variable.CopyIn(0,fTemp_nine_values);
//            out_variable.CopyIn(0,fTemp_six_values);

//        Put_values_In_dArrayT_vector(fEulerian_strain_Elements_IPs, e,l,fTemp_six_values);
          Put_values_In_dArrayT_vector(fEulerian_strain_Elements_IPs, e,l,fTemp_nine_values);

//        out_variable.CopyIn(6,fTemp_six_values);
          out_variable.CopyIn(9,fTemp_nine_values);//!!9->6?
          out_variable[18]=fState_variables_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kc);
          out_variable[19]=fState_variables_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+khc);
          out_variable[20]=fState_variables_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kc_chi);
          out_variable[21]=fState_variables_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+khc_chi);
          out_variable[22]=fState_variables_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kDelgamma);
          out_variable[23]=fState_variables_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kDelgammachi);
          out_variable[24]=fState_variables_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+ktrCauchy_Stress);
          out_variable[25]=fState_variables_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kNorm_Dev_Cauchy_Stress);
          out_variable[26]=fState_variables_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+ktrRel);
          out_variable[27]=fState_variables_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kRel_inv);
          out_variable[28]=fState_variables_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+ktrm);
          out_variable[29]=fState_variables_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+km_inv);
          out_variable[30]=fState_variables_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+ktrS);
          out_variable[31]=fState_variables_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kinvdevS);
          out_variable[32]=fState_variables_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+ktrSIGMA_S);
          out_variable[33]=fState_variables_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kinvdevSIGMA_S);
          out_variable[34]=fState_variables_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kinvtrM);
          out_variable[35]=fState_variables_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kinvdevM);
          out_variable[36]=fState_variables_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kinvPhi);
          out_variable[37]=fState_variables_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kinvGPhi);
          out_variable[38]=fState_variables_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+ktreps);
          out_variable[39]=fState_variables_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kdeveps);
          out_variable[40]=fState_variables_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kinvtrgammastn);
          out_variable[41]=fState_variables_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kinvdevgammastn);

        }
    }
    else
    { //-- Still Iterating

        /* residual and tangent for displacements */
        const double* Det    = fShapes_displ->IPDets();
        const double* Weight = fShapes_displ->IPWeights();
        fShapes_displ->TopIP();
        fShapes_micro->TopIP();
        if(iConstitutiveModelType==2)
        {
            SigN_IPs_el_n.RowCopy(e,SigN_IPs_n);
            mn_IPs_el_n.RowCopy(e,mn_IPs_n);
            GammaN_IPs_el_n.RowCopy(e,GammaN_IPs_n);
            sn_sigman_IPs_el_n.RowCopy(e,sn_sigman_IPs_n);

            Fn_ar_IPs_el.RowCopy(e,Fn_ar_IPs);
            FnInv_ar_IPs_el.RowCopy(e,FnInv_ar_IPs);
            ChiN_ar_IPs_el_n.RowCopy(e,ChiN_ar_IPs_n);
            GRAD_ChiN_ar_IPs_el_n.RowCopy(e,GRAD_ChiN_ar_IPs);
        }
        if(iConstitutiveModelType==3)
        {

        	         /* retrieve dGdCauchy_Stress and dGdCauchy_Stress_n in element */
        	         fdGdCauchy_Stress_n_Elements_IPs.RowCopy(e,fdGdCauchy_Stress_n_IPs);
        	         fdGdCauchy_Stress_Elements_IPs.RowCopy(e,fdGdCauchy_Stress_IPs);

        	         /* retrieve Cauchy_Stress and Cauchy_Stress_n in element */
        	         fCauchy_stress_Elements_n_IPs.RowCopy(e,fCauchy_stress_n_IPs);
        	         fCauchy_stress_Elements_IPs.RowCopy(e,fCauchy_stress_IPs);

        	         /* retrieve Deformation Gradient and Deformation Gradient_n in element */
        	         fDeformation_Gradient_n_Elements_IPs.RowCopy(e,fDeformation_Gradient_n_IPs);
        	         fDeformation_Gradient_Elements_IPs.RowCopy(e,fDeformation_Gradient_IPs);

        	         /* retrieve ISVs and ISVs_n in element */
        	         fState_variables_n_Elements_IPs.RowCopy(e,fState_variables_n_IPs);
        	         fState_variables_Elements_IPs.RowCopy(e,fState_variables_IPs);

        }

        while (fShapes_displ->NextIP() && fShapes_micro->NextIP())
        {

                double scale_const = (*Weight++)*(*Det++);

                const int  IP = fShapes_displ->CurrIP();
                dArrayT DisplIPCoordinate(n_sd), MicroIPCoordinate(n_sd);
                fShapes_displ->IPCoords(DisplIPCoordinate);
                fShapes_micro->IPCoords(MicroIPCoordinate);


                const double* shapes_displ_X = fShapes_displ->IPShapeX();
                /* [fShapeDispl] will be formed */
                Form_solid_shape_functions(shapes_displ_X);//output:fShapeDispl
                Form_Gradient_of_solid_shape_functions(fShapeDisplGrad_temp);//output:fShapeDisplGrad in Reference config.
                Form_GRAD_Nuw_matrix(fShapeDisplGrad_temp) ;//output:GRAD_Nuw

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
                //fs_micromorph3D_out << "Gauss Point = "<< IP << endl;
                //fs_micromorph3D_out << "Current Deformation Gradient = " << fDeformation_Gradient << endl;
                /* [fDeformation_Gradient_Inverse] and [fDeformation_Gradient_Transpose] will be formed */
                if (fDeformation_Gradient.Det()==0)
                fDeformation_Gradient = fIdentity_matrix;
                fDeformation_Gradient_Inverse.Inverse(fDeformation_Gradient);
                // Form Chi[i][j] and  dMatrixT ChiM(i,j)
                Form_micro_deformation_tensor_Chi();//output: Chi[i][j]
                ChiM_Inverse=0.0;
                ChiM_Inverse.Inverse(ChiM);
                Form_ChiM();//It is also micro-deformation gradient tensor but defined as dMatrixT
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
                double invJ=1/J;
                /* Jacobian for the current IP will be saved */


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

                /* [fIota_temp_matrix] will be formed */
                fIota_temp_matrix.MultATB(fShapeDisplGrad,fDefGradInv_Grad_grad);

                /* [fIota_w_temp_matrix] will be formed*/
                fIota_w_temp_matrix.MultATBT(GRAD_Nuw,Finv_w);
                /* [fIota_eta_temp_matrix] will be formed*/
                fIota_eta_temp_matrix.MultATBT(GRAD_NCHI,Finv_eta);

                double	scale = scale_const;

                if(iConstitutiveModelType==3)
                {

                    fdGdCauchy_Stress_n_IPs.RowCopy(IP,fdGdCauchy_Stress_n);
                    fCauchy_stress_n_IPs.RowCopy(IP,fCauchy_stress_tensor_current_n_IP);
                    fDeformation_Gradient_n_IPs.RowCopy(IP,fDeformation_Gradient_n);

					fTemp_matrix_nsd_x_nsd = fDeformation_Gradient;
					fTemp_matrix_nsd_x_nsd-=fDeformation_Gradient_n;
					fDeformation_Gradient_Inverse.Inverse(fDeformation_Gradient);
					fVelocity_Gradient_current_IP.MultAB(fTemp_matrix_nsd_x_nsd,fDeformation_Gradient_Inverse);
					fVelocity_Gradient_current_IP_transpose.Transpose(fVelocity_Gradient_current_IP);
					fVelocity_Gradient_current_IP_trace = fVelocity_Gradient_current_IP.Trace();
					fdGdCauchy_Stress_n_Transpose.Transpose(fdGdCauchy_Stress_n);
					fdGdCauchy_Stress_n_trace = fdGdCauchy_Stress_n.Trace();

					fTemp_matrix_nsd_x_nsd = fVelocity_Gradient_current_IP;
					fTemp_matrix_nsd_x_nsd+= fVelocity_Gradient_current_IP_transpose;
					fSymmetric_part_Velocity_Gradient_current_IP = fTemp_matrix_nsd_x_nsd;
					fSymmetric_part_Velocity_Gradient_current_IP*= 0.5;
					fSymmetric_part_Velocity_Gradient_current_IP_trace = fSymmetric_part_Velocity_Gradient_current_IP.Trace();

					fs_micromorph3D_out << "Element = " << e << endl;
					fs_micromorph3D_out << "Integration point = " << IP << endl;
					fs_micromorph3D_out << "Global Iteration = " << global_iteration << endl;
					fs_micromorph3D_out << "Current Deformation Gradient = " << fDeformation_Gradient << endl;
					fs_micromorph3D_out << "Previous time step Deformation Gradient = " << fDeformation_Gradient_n << endl;
					fs_micromorph3D_out << "fCauchy_stress_tensor_current_n_IP = " << fCauchy_stress_tensor_current_n_IP << endl;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// Forming Cauchy Stress Trial///////////////////////////////////////
					/////////////////////////////////////////////////////////////////////////
					/////////////////////////////////////////////////////////////////////////
					fCauchy_stress_tensor_current_IP_tr = fCauchy_stress_tensor_current_n_IP;
					fTemp_matrix_nsd_x_nsd.SetToScaled(fSymmetric_part_Velocity_Gradient_current_IP_trace,fCauchy_stress_tensor_current_n_IP);
					fCauchy_stress_tensor_current_IP_tr-= fTemp_matrix_nsd_x_nsd;

					fTemp_matrix_nsd_x_nsd.MultAB(fVelocity_Gradient_current_IP,fCauchy_stress_tensor_current_n_IP);
					fCauchy_stress_tensor_current_IP_tr+= fTemp_matrix_nsd_x_nsd;

					fTemp_matrix_nsd_x_nsd.MultABT(fCauchy_stress_tensor_current_n_IP,fVelocity_Gradient_current_IP);
					fCauchy_stress_tensor_current_IP_tr+= fTemp_matrix_nsd_x_nsd;


					fTemp_matrix_one_x_one = fMaterial_Params[kLambda];
					fTemp_matrix_one_x_one*= fVelocity_Gradient_current_IP_trace;
					fTemp_matrix_nsd_x_nsd.SetToScaled(fTemp_matrix_one_x_one,fIdentity_matrix);
					fCauchy_stress_tensor_current_IP_tr+= fTemp_matrix_nsd_x_nsd;

					fTemp_matrix_nsd_x_nsd.SetToScaled(2*(fMaterial_Params[kMu]),fSymmetric_part_Velocity_Gradient_current_IP);
					fCauchy_stress_tensor_current_IP_tr+= fTemp_matrix_nsd_x_nsd;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
					//////////////////////////Mean Cauchy Stress/////////////////////////////
					////////////////////////////////////////////////////////////////////////
					mean_Cauchy_stress_tr = fCauchy_stress_tensor_current_IP_tr.Trace()/3;
					fTemp_matrix_nsd_x_nsd.SetToScaled(mean_Cauchy_stress_tr,fIdentity_matrix);
					dev_Cauchy_stress_tr = fCauchy_stress_tensor_current_IP_tr;
					dev_Cauchy_stress_tr-= fTemp_matrix_nsd_x_nsd;
					fTemp_matrix_one_x_one= dev_Cauchy_stress_tr.ScalarProduct();
					fNorm_dev_Cauchy_stress_tensor_current_IP_tr = sqrt(fTemp_matrix_one_x_one);
					////////////////////////////////////////////////////////////////////////////
					////////////////////////////////////////////////////////////////////////////
					fdGdCauchy_Stress_tr = 0.0;
					fdGdCauchy_Stress_tr.SetToScaled(Bpsi*1.0/3.0,fIdentity_matrix);
					fTemp_matrix_nsd_x_nsd.SetToScaled(1/fNorm_dev_Cauchy_stress_tensor_current_IP_tr,dev_Cauchy_stress_tr);
					fdGdCauchy_Stress_tr+=fTemp_matrix_nsd_x_nsd;
					fdGdCauchy_Stress_tr_Transpose.Transpose(fdGdCauchy_Stress_tr);
                    fdGdCauchy_Stress_tr_trace = fdGdCauchy_Stress_tr.Trace();
					////////////////////////////////////////////////////////////////////////////

                    /* Form terms related the cohesion and friction angle  in D-P yield function */
                    /* Initially Aphi is already assigined to fState_variables_n_IPs(IP,khc)=Apsi in TakeParameterList function*/
                    double Beta=-1.0;
                    double Aphi=2*sqrt(6)*cos(fMaterial_Params[kFphi])/(3+Beta*sin(fMaterial_Params[kFphi]));
                    double Bphi=2*sqrt(6)*sin(fMaterial_Params[kFphi])/(3+Beta*sin(fMaterial_Params[kFphi]));
                    // Form the cohesion and dilation angle related terms in Plastic potential function
                    double Apsi=2*sqrt(6)*cos(fMaterial_Params[kDpsi])/(3+Beta*sin(fMaterial_Params[kDpsi] ));
                    double Bpsi=2*sqrt(6)*sin(fMaterial_Params[kDpsi])/(3+Beta*sin(fMaterial_Params[kDpsi] ));


                    /* Form terms related the cohesion and friction angle  in Micro D-P yield function */
                    /* Initially Aphi_chi is already assigined to fState_variables_n_IPs(IP,khc_chi)=Apsi_chi in TakeParameterList function*/
                    //double Beta=-1.0;
                    double Aphi_chi=2*sqrt(6)*cos(fMaterial_Params[kFphi_chi])/(3+Beta*sin(fMaterial_Params[kFphi_chi]));
                    double Bphi_chi=2*sqrt(6)*sin(fMaterial_Params[kFphi_chi])/(3+Beta*sin(fMaterial_Params[kFphi_chi]));
                    // Form the cohesion and dilation angle related terms in Plastic potential function
                    double Apsi_chi=2*sqrt(6)*cos(fMaterial_Params[kDpsi_chi])/(3+Beta*sin(fMaterial_Params[kDpsi_chi] ));
                    double Bpsi_chi=2*sqrt(6)*sin(fMaterial_Params[kDpsi_chi])/(3+Beta*sin(fMaterial_Params[kDpsi_chi] ));




                    	fYield_function_tr=0.0;
                    	fYield_function_tr=fNorm_dev_Cauchy_stress_tensor_current_IP_tr-(Aphi*fState_variables_n_IPs(IP,kc)-Bphi*mean_Cauchy_stress_tr);
                        fF_tr_fact = fYield_function_tr/(fMaterial_Params[kMu]);



					if (global_iteration < 0) fF_tr_fact = -1.0;





			if(fF_tr_fact>dYieldTrialTol)
			{

				fs_micromorph3D_out<<"Macro-Plasticity at the Current Configuration" <<endl;
				// initialize before iteration
				fYield_function=fYield_function_tr;
				fCauchy_stress_tensor_current_IP=fCauchy_stress_tensor_current_IP_tr;
				fdev_Cauchy_stress_tensor_current_IP=dev_Cauchy_stress_tr;
				fNorm_dev_Cauchy_stress_tensor_current_IP=fNorm_dev_Cauchy_stress_tensor_current_IP_tr;

				fdelDelgamma = 0.0;
				fDelgamma = 0.0;

				// iterate using Newton-Raphson to solve for fDelgamma
				iter_count = 0;
				fs_micromorph3D_out << "Gauss Point = " << IP << endl;
				fs_micromorph3D_out << "Current Macro Yield function = " << fYield_function << endl;

				if (iPlasticityCheck==0)
				{
					while (fabs(fYield_function) > dAbsTol && fabs(fYield_function/fYield_function_tr) > dRelTol
							&& iter_count < iIterationMax)
					{
						iter_count += 1;

						fdCauchy_stressdDelgamma.SetToScaled(fdGdCauchy_Stress_tr_trace,fCauchy_stress_tensor_current_n_IP);

						fTemp_matrix_nsd_x_nsd.MultAB(fdGdCauchy_Stress_tr,fCauchy_stress_tensor_current_n_IP);
						fdCauchy_stressdDelgamma-= fTemp_matrix_nsd_x_nsd;

						fTemp_matrix_nsd_x_nsd.MultABT(fCauchy_stress_tensor_current_n_IP,fdGdCauchy_Stress_tr);
						fdCauchy_stressdDelgamma-= fTemp_matrix_nsd_x_nsd;

						fTemp_matrix_nsd_x_nsd.SetToScaled((fMaterial_Params[kLambda]*fdGdCauchy_Stress_tr_trace),fIdentity_matrix);
						fdCauchy_stressdDelgamma-= fTemp_matrix_nsd_x_nsd;

						fTemp_matrix_nsd_x_nsd.SetToScaled(fMaterial_Params[kMu],fdGdCauchy_Stress_tr);
						fdCauchy_stressdDelgamma-= fTemp_matrix_nsd_x_nsd;

						fTemp_matrix_nsd_x_nsd.SetToScaled(fMaterial_Params[kMu],fdGdCauchy_Stress_tr_Transpose);
						fdCauchy_stressdDelgamma-= fTemp_matrix_nsd_x_nsd;

						fdmean_Cauchy_stressdDelgamma = fdCauchy_stressdDelgamma.Trace()/3;
						fTemp_matrix_nsd_x_nsd.SetToScaled(fdmean_Cauchy_stressdDelgamma,fIdentity_matrix);
						fdDev_Cauchy_stressdDelgamma = fTemp_matrix_nsd_x_nsd;
						fdDev_Cauchy_stressdDelgamma*= -1;
						fdDev_Cauchy_stressdDelgamma+= fdCauchy_stressdDelgamma;
						fTemp_matrix_nsd_x_nsd.SetToScaled(1/fNorm_dev_Cauchy_stress_tensor_current_IP,fdev_Cauchy_stress_tensor_current_IP);
						fNorm_dDev_Cauchy_stressdDelgamma = dMatrixT::Dot(fdDev_Cauchy_stressdDelgamma,fTemp_matrix_nsd_x_nsd);


						dcdDelgamma=fState_variables_n_IPs(IP,khc)*fMaterial_Params[kHc];

						dFYdDelgamma=fNorm_dDev_Cauchy_stressdDelgamma-(Aphi*dcdDelgamma-Bphi*fdmean_Cauchy_stressdDelgamma);

						if (fabs(dFYdDelgamma) >= 1e-12) fdelDelgamma = -fYield_function/dFYdDelgamma;
						else fdelDelgamma = 0.0;

						fDelgamma+= fdelDelgamma;
						if (fDelgamma < 0.0) fDelgamma = 0.0;

						fTemp_matrix_nsd_x_nsd.SetToScaled(fDelgamma,fdGdCauchy_Stress_tr);
						fElastic_Velocity_Gradient_current_IP = fVelocity_Gradient_current_IP;
						fElastic_Velocity_Gradient_current_IP-= fTemp_matrix_nsd_x_nsd;
						fElastic_Velocity_Gradient_current_IP_trace = fElastic_Velocity_Gradient_current_IP.Trace();
						fElastic_Velocity_Gradient_current_IP_transpose.Transpose(fElastic_Velocity_Gradient_current_IP);

						fCauchy_stress_tensor_current_IP = fCauchy_stress_tensor_current_n_IP;
						fTemp_matrix_nsd_x_nsd.SetToScaled(fElastic_Velocity_Gradient_current_IP_trace,fCauchy_stress_tensor_current_n_IP);
						fCauchy_stress_tensor_current_IP-= fTemp_matrix_nsd_x_nsd;

						fTemp_matrix_nsd_x_nsd.MultAB(fElastic_Velocity_Gradient_current_IP,fCauchy_stress_tensor_current_n_IP);
						fCauchy_stress_tensor_current_IP+= fTemp_matrix_nsd_x_nsd;

						fTemp_matrix_nsd_x_nsd.MultAB(fCauchy_stress_tensor_current_n_IP,fElastic_Velocity_Gradient_current_IP_transpose);
						fCauchy_stress_tensor_current_IP+= fTemp_matrix_nsd_x_nsd;

						fTemp_matrix_one_x_one = fMaterial_Params[kLambda];
						fTemp_matrix_one_x_one*= fElastic_Velocity_Gradient_current_IP_trace;
						fTemp_matrix_nsd_x_nsd.SetToScaled(fTemp_matrix_one_x_one,fIdentity_matrix);
						fCauchy_stress_tensor_current_IP+= fTemp_matrix_nsd_x_nsd;

						fTemp_matrix_nsd_x_nsd.SetToScaled(fMaterial_Params[kMu],fElastic_Velocity_Gradient_current_IP);
						fCauchy_stress_tensor_current_IP+= fTemp_matrix_nsd_x_nsd;

						fTemp_matrix_nsd_x_nsd.SetToScaled(fMaterial_Params[kMu],fElastic_Velocity_Gradient_current_IP_transpose);
						fCauchy_stress_tensor_current_IP+= fTemp_matrix_nsd_x_nsd;

						fCauchy_stress_tensor_current_IP_trace = fCauchy_stress_tensor_current_IP.Trace();
						mean_Cauchy_stress = fCauchy_stress_tensor_current_IP_trace;
						mean_Cauchy_stress*= 1.0/3.0;

						fdev_Cauchy_stress_tensor_current_IP = fCauchy_stress_tensor_current_IP;
						fTemp_matrix_nsd_x_nsd.SetToScaled(mean_Cauchy_stress,fIdentity_matrix);
						fdev_Cauchy_stress_tensor_current_IP-= fTemp_matrix_nsd_x_nsd;


						fTemp_matrix_one_x_one = fdev_Cauchy_stress_tensor_current_IP.ScalarProduct();
						fNorm_dev_Cauchy_stress_tensor_current_IP = sqrt(fTemp_matrix_one_x_one);

						fs_micromorph3D_out  << "fCauchy_stress_tensor_current_IP= " << fCauchy_stress_tensor_current_IP << endl;
						fs_micromorph3D_out  << "mean_Cauchy_stress= " << mean_Cauchy_stress << endl;
						fs_micromorph3D_out  << "fDelgamma= " << fDelgamma << endl;

						cohesion=fState_variables_n_IPs(IP,kc)+fDelgamma*fState_variables_n_IPs(IP,khc)*fMaterial_Params[kHc];
						if (cohesion < 0.0)
						{
							cohesion = 0.0;
							fState_variables_IPs(IP,kc)= fState_variables_n_IPs(IP,kc);
						}

						fYield_function=fNorm_dev_Cauchy_stress_tensor_current_IP-(Aphi*cohesion-Bphi*mean_Cauchy_stress);

						fs_micromorph3D_out  << "Current relative residual = " << fabs(fYield_function/fYield_function_tr) << endl;
					}
				}
				if (iPlasticityCheck==1)
				{

					mean_Cauchy_stress_n=fCauchy_stress_tensor_current_n_IP.Trace()/3;
					dev_Cauchy_stress_n.SetToScaled(mean_Cauchy_stress_n,fIdentity_matrix);
					dev_Cauchy_stress_n*=-1;
					dev_Cauchy_stress_n+=fCauchy_stress_tensor_current_n_IP;

					Predictor_stress_terms.SetToScaled(fVelocity_Gradient_current_IP_trace,dev_Cauchy_stress_n);
					Predictor_stress_terms*=-1;
					Predictor_stress_terms+=dev_Cauchy_stress_n;

					fTemp_matrix_nsd_x_nsd.MultAB(fVelocity_Gradient_current_IP,fCauchy_stress_tensor_current_n_IP);
					Predictor_stress_terms+=fTemp_matrix_nsd_x_nsd;

					fTemp_matrix_nsd_x_nsd.MultABT(fCauchy_stress_tensor_current_n_IP,fVelocity_Gradient_current_IP);
					Predictor_stress_terms+=fTemp_matrix_nsd_x_nsd;

					fTemp_matrix_nsd_x_nsd.SetToScaled(fMaterial_Params[kMu],fVelocity_Gradient_current_IP);
					Predictor_stress_terms+= fTemp_matrix_nsd_x_nsd;

					fTemp_matrix_nsd_x_nsd.SetToScaled(fMaterial_Params[kMu],fVelocity_Gradient_current_IP_transpose);
					Predictor_stress_terms+= fTemp_matrix_nsd_x_nsd;

					fTemp_matrix_nsd_x_nsd.MultAB(fVelocity_Gradient_current_IP,fCauchy_stress_tensor_current_n_IP);
					Temp_trace_value = fTemp_matrix_nsd_x_nsd.Trace()/3;
					fTemp_matrix_nsd_x_nsd.SetToScaled(Temp_trace_value,fIdentity_matrix);
					Predictor_stress_terms-= fTemp_matrix_nsd_x_nsd;

					fTemp_matrix_nsd_x_nsd.MultABT(fCauchy_stress_tensor_current_n_IP,fVelocity_Gradient_current_IP);
					Temp_trace_value = fTemp_matrix_nsd_x_nsd.Trace()/3;
					fTemp_matrix_nsd_x_nsd.SetToScaled(Temp_trace_value,fIdentity_matrix);
					Predictor_stress_terms-= fTemp_matrix_nsd_x_nsd;

					fTemp_matrix_nsd_x_nsd.SetToScaled(fMaterial_Params[kMu]*fVelocity_Gradient_current_IP_trace*2.0/3.0,fIdentity_matrix);
					Predictor_stress_terms-= fTemp_matrix_nsd_x_nsd;

					////////////////////////////////////////// Corrector//////////////////////////////////////////////////////////

					Corrector_stress_terms.SetToScaled(fdGdCauchy_Stress_tr_trace,dev_Cauchy_stress_n);

					fTemp_matrix_nsd_x_nsd.MultAB(fdGdCauchy_Stress_tr,fCauchy_stress_tensor_current_n_IP);
					Corrector_stress_terms-=fTemp_matrix_nsd_x_nsd;

					fTemp_matrix_nsd_x_nsd.MultAB(fdGdCauchy_Stress_tr,fCauchy_stress_tensor_current_n_IP);
					Temp_trace_value = fTemp_matrix_nsd_x_nsd.Trace()/3;
					fTemp_matrix_nsd_x_nsd.SetToScaled(Temp_trace_value,fIdentity_matrix);
					Corrector_stress_terms+= fTemp_matrix_nsd_x_nsd;

					fTemp_matrix_nsd_x_nsd.MultABT(fCauchy_stress_tensor_current_n_IP,fdGdCauchy_Stress_tr);
					Corrector_stress_terms-=fTemp_matrix_nsd_x_nsd;

					fTemp_matrix_nsd_x_nsd.MultABT(fCauchy_stress_tensor_current_n_IP,fdGdCauchy_Stress_tr);
					Temp_trace_value = fTemp_matrix_nsd_x_nsd.Trace()/3;
					fTemp_matrix_nsd_x_nsd.SetToScaled(Temp_trace_value,fIdentity_matrix);
					Corrector_stress_terms+=fTemp_matrix_nsd_x_nsd;

					fTemp_matrix_nsd_x_nsd.SetToScaled(fMaterial_Params[kMu],fdGdCauchy_Stress_tr);
					Corrector_stress_terms-= fTemp_matrix_nsd_x_nsd;

					fTemp_matrix_nsd_x_nsd2.SetToScaled(fMaterial_Params[kMu],fdGdCauchy_Stress_tr_Transpose);
					Corrector_stress_terms-= fTemp_matrix_nsd_x_nsd;

					fTemp_matrix_nsd_x_nsd.SetToScaled(fMaterial_Params[kMu]*fdGdCauchy_Stress_tr_trace*2.0/3.0,fIdentity_matrix);
					Corrector_stress_terms+= fTemp_matrix_nsd_x_nsd;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////// Mean Cauchy stress = predictor_Meanstress_terms + (delgamma at current configuration)* corrector_Mean stress_terms//////////

					fTemp_matrix_one_x_one = mean_Cauchy_stress_n;
					fTemp_matrix_one_x_one*= fVelocity_Gradient_current_IP_trace;
					fTemp_matrix_one_x_one*= -1.0;
					Predictor_mean_stress_terms = fTemp_matrix_one_x_one;
					Predictor_mean_stress_terms+= mean_Cauchy_stress_n;

					fTemp_matrix_nsd_x_nsd.MultAB(fVelocity_Gradient_current_IP,fCauchy_stress_tensor_current_n_IP);
					Temp_trace_value = fTemp_matrix_nsd_x_nsd.Trace()/3;
					Predictor_mean_stress_terms+= Temp_trace_value;

					fTemp_matrix_nsd_x_nsd.MultABT(fCauchy_stress_tensor_current_n_IP,fVelocity_Gradient_current_IP);
					Temp_trace_value = fTemp_matrix_nsd_x_nsd.Trace()/3;
					Predictor_mean_stress_terms+= Temp_trace_value;

					Temp_trace_value = fVelocity_Gradient_current_IP_trace;
					Temp_trace_value*= fMaterial_Params[kLambda];
					Predictor_mean_stress_terms+= Temp_trace_value;

					Temp_trace_value = fVelocity_Gradient_current_IP_trace;
					Temp_trace_value*= fMaterial_Params[kMu];
					Temp_trace_value*= 2.0/3.0;
					Predictor_mean_stress_terms+= Temp_trace_value;


					fTemp_matrix_one_x_one = fdGdCauchy_Stress_tr.Trace();
					fTemp_matrix_one_x_one*= mean_Cauchy_stress_n;
					Corrector_mean_stress_terms+= fTemp_matrix_one_x_one;

					fTemp_matrix_nsd_x_nsd.MultAB(fdGdCauchy_Stress_tr,fCauchy_stress_tensor_current_n_IP);
					Temp_trace_value = fTemp_matrix_nsd_x_nsd.Trace()/3;
					Corrector_mean_stress_terms-= Temp_trace_value;

					fTemp_matrix_nsd_x_nsd.MultABT(fCauchy_stress_tensor_current_n_IP,fdGdCauchy_Stress_tr);
					Temp_trace_value = fTemp_matrix_nsd_x_nsd.Trace()/3;
					Corrector_mean_stress_terms-= Temp_trace_value;

					Temp_trace_value = fMaterial_Params[kLambda]*fdGdCauchy_Stress_tr_trace;
					Corrector_mean_stress_terms-= Temp_trace_value;

					Temp_trace_value = fdGdCauchy_Stress_tr.Trace();
					Temp_trace_value*= fMaterial_Params[kMu];
					Temp_trace_value*= 2.0/3.0;
					Corrector_mean_stress_terms-= Temp_trace_value;


									///////////////////////////////////////////////////////////////////////
//////////////////////////////////// Forming A(delgama_current_config)^2+B(delgama_current_config)+D=0///////////////////////////////////////////
									///////////////////////////////////////////////////////////////////////



					fTemp_matrix_one_x_one = Aphi;
					fTemp_matrix_one_x_one*= fState_variables_n_IPs(IP,kc);
					ftemp_D_coeff = Predictor_mean_stress_terms;
					ftemp_D_coeff*= Bphi;
					ftemp_D_coeff*= -1.0;
					ftemp_D_coeff+= fTemp_matrix_one_x_one;
					fD_coeff = ftemp_D_coeff;
					fD_coeff*= fD_coeff;
					fD_coeff*= -1.0;


					fTemp_matrix_one_x_one = Predictor_stress_terms.ScalarProduct();
					fD_coeff+= fTemp_matrix_one_x_one;

					fTemp_matrix_one_x_one = Aphi;
					fTemp_matrix_one_x_one*= fState_variables_n_IPs(IP,khc);
					fTemp_matrix_one_x_one*= fMaterial_Params[kHc];


					ftemp_B_coeff = Bphi;
					ftemp_B_coeff*= Corrector_mean_stress_terms;
					ftemp_B_coeff*= -1.0;
					ftemp_B_coeff+= fTemp_matrix_one_x_one;
					fB_coeff = ftemp_B_coeff;
					fB_coeff*= 2.0;
					fB_coeff*= ftemp_D_coeff;
					fB_coeff*= -1.0;


					fTemp_matrix_one_x_one = dMatrixT::Dot(Predictor_stress_terms,Corrector_stress_terms);
					fB_coeff+= fTemp_matrix_one_x_one;

					fTemp_matrix_one_x_one = dMatrixT::Dot(Corrector_stress_terms,Predictor_stress_terms);
					fB_coeff+= fTemp_matrix_one_x_one;


					fA_coeff = ftemp_B_coeff;
					fA_coeff*= fA_coeff;
					fA_coeff*= -1.0;

					fTemp_matrix_one_x_one = Corrector_stress_terms.ScalarProduct();
					fA_coeff+= fTemp_matrix_one_x_one;



					fDelta_function = fB_coeff;
					fDelta_function*= fDelta_function;
					fTemp_matrix_one_x_one = fA_coeff;
					fTemp_matrix_one_x_one*= 4.0;
					fTemp_matrix_one_x_one*= fD_coeff;
					fDelta_function-= fTemp_matrix_one_x_one;

					if (fDelta_function < 0.0 )
							{
							fs_micromorph3D_out << "The fDelta_function is negative So imaginary value for delgamma will be obtained " << endl;
							ExceptionT::GeneralFail(caller,"The fDelta_function is negative So imaginary value for delgamma will be obtained for %d element", e);
							}
					fDelta_function = sqrt(fDelta_function);



					fDelgamma_current_configuration_root_one = 0.0;
					fDelgamma_current_configuration_root_one = fB_coeff;
					fDelgamma_current_configuration_root_one*= -1.0;

					fDelgamma_current_configuration_root_one+= fDelta_function;
					fDelgamma_current_configuration_root_one/= fA_coeff;
					fDelgamma_current_configuration_root_one/= 2.0;

					fDelgamma_current_configuration_root_two = 0.0;
					fDelgamma_current_configuration_root_two = fB_coeff;
					fDelgamma_current_configuration_root_two*= -1.0;

					fDelgamma_current_configuration_root_two-= fDelta_function;
					fDelgamma_current_configuration_root_two/= fA_coeff;
					fDelgamma_current_configuration_root_two/= 2.0;


					if (fDelgamma_current_configuration_root_two > 0.0 && fDelgamma_current_configuration_root_one > 0.0 )
						{
						if (fDelgamma_current_configuration_root_two > fDelgamma_current_configuration_root_one)
							{
							fDelgamma = fDelgamma_current_configuration_root_one;
							}
						else
							{
							fDelgamma = fDelgamma_current_configuration_root_two;
							}
						}
					if (fDelgamma_current_configuration_root_two > 0.0 && fDelgamma_current_configuration_root_one < 0.0 )
					{
						fDelgamma = fDelgamma_current_configuration_root_two;
					}

					if (fDelgamma_current_configuration_root_two < 0.0 && fDelgamma_current_configuration_root_one > 0.0 )
					{
						fDelgamma = fDelgamma_current_configuration_root_one;
					}

					fs_micromorph3D_out<< "fDelgamma_current_configuration_root_one="<< fDelgamma_current_configuration_root_one<<endl;
					fs_micromorph3D_out<< "fDelgamma_current_configuration_root_two="<< fDelgamma_current_configuration_root_two<<endl;



					fTemp_matrix_nsd_x_nsd.SetToScaled(fDelgamma,fdGdCauchy_Stress_tr);
					fElastic_Velocity_Gradient_current_IP = fVelocity_Gradient_current_IP;
					fElastic_Velocity_Gradient_current_IP-= fTemp_matrix_nsd_x_nsd;
					fElastic_Velocity_Gradient_current_IP_trace = fElastic_Velocity_Gradient_current_IP.Trace();
					fElastic_Velocity_Gradient_current_IP_transpose.Transpose(fElastic_Velocity_Gradient_current_IP);

					fCauchy_stress_tensor_current_IP = fCauchy_stress_tensor_current_n_IP;
					fTemp_matrix_nsd_x_nsd.SetToScaled(fElastic_Velocity_Gradient_current_IP_trace,fCauchy_stress_tensor_current_n_IP);
					fCauchy_stress_tensor_current_IP-= fTemp_matrix_nsd_x_nsd;

					fTemp_matrix_nsd_x_nsd.MultAB(fElastic_Velocity_Gradient_current_IP,fCauchy_stress_tensor_current_n_IP);
					fCauchy_stress_tensor_current_IP+= fTemp_matrix_nsd_x_nsd;

					fTemp_matrix_nsd_x_nsd.MultABT(fCauchy_stress_tensor_current_n_IP,fElastic_Velocity_Gradient_current_IP);
					fCauchy_stress_tensor_current_IP+= fTemp_matrix_nsd_x_nsd;

					fTemp_matrix_one_x_one = fMaterial_Params[kLambda];
					fTemp_matrix_one_x_one*= fElastic_Velocity_Gradient_current_IP_trace;
					fTemp_matrix_nsd_x_nsd.SetToScaled(fTemp_matrix_one_x_one,fIdentity_matrix);
					fCauchy_stress_tensor_current_IP+= fTemp_matrix_nsd_x_nsd;

					fTemp_matrix_nsd_x_nsd.SetToScaled(fMaterial_Params[kMu],fElastic_Velocity_Gradient_current_IP);
					fCauchy_stress_tensor_current_IP+= fTemp_matrix_nsd_x_nsd;

					fTemp_matrix_nsd_x_nsd.SetToScaled(fMaterial_Params[kMu],fElastic_Velocity_Gradient_current_IP_transpose);
					fCauchy_stress_tensor_current_IP+= fTemp_matrix_nsd_x_nsd;

					fs_micromorph3D_out  << "fCauchy_stress_tensor_current_IP_its_formulation = " << fCauchy_stress_tensor_current_IP << endl;

					fCauchy_stress_tensor_current_IP_trace = fCauchy_stress_tensor_current_IP.Trace();
					mean_Cauchy_stress = fCauchy_stress_tensor_current_IP_trace;
					mean_Cauchy_stress*= 1.0/3.0;
					fs_micromorph3D_out  << "mean_Cauchy_stress formulation = " << mean_Cauchy_stress << endl;

					/*
					mean_Cauchy_stress = Corrector_mean_stress_terms;
					mean_Cauchy_stress*= fDelgamma;
					mean_Cauchy_stress+= Predictor_mean_stress_terms;

					fCauchy_stress_tensor_current_IP = 0.0;
					fCauchy_stress_tensor_current_IP = Predictor_stress_terms;
					fTemp_matrix_nsd_x_nsd.SetToScaled(fDelgamma,Corrector_stress_terms);
					fCauchy_stress_tensor_current_IP+= fTemp_matrix_nsd_x_nsd;
					fTemp_matrix_nsd_x_nsd2.SetToScaled(mean_Cauchy_stress,fIdentity_matrix);
					fCauchy_stress_tensor_current_IP+= fTemp_matrix_nsd_x_nsd2;
*/

					fdev_Cauchy_stress_tensor_current_IP.SetToScaled(fDelgamma,Corrector_stress_terms);
					fdev_Cauchy_stress_tensor_current_IP+= Predictor_stress_terms;
					fNorm_dev_Cauchy_stress_tensor_current_IP = fdev_Cauchy_stress_tensor_current_IP.ScalarProduct();
					fNorm_dev_Cauchy_stress_tensor_current_IP = sqrt(fNorm_dev_Cauchy_stress_tensor_current_IP);

					cohesion=fState_variables_n_IPs(IP,kc)+fDelgamma*fState_variables_n_IPs(IP,khc)*fMaterial_Params[kHc];

					if (cohesion < 0.0)
						{
						cohesion = 0.0;
						fState_variables_IPs(IP,kc)= fState_variables_n_IPs(IP,kc);
						}

					fYield_function=fNorm_dev_Cauchy_stress_tensor_current_IP-(Aphi*cohesion-Bphi*mean_Cauchy_stress);
					fs_micromorph3D_out  << "Current relative residual = " << fabs(fYield_function/fYield_function_tr) << endl;
					fs_micromorph3D_out  << "Current residual = " << fabs(fYield_function) << endl;
				}
                fs_micromorph3D_out  << "Current relative residual = " << fabs(fYield_function/fYield_function_tr) << endl;
                fs_micromorph3D_out  << "Yield Function= " << fabs(fYield_function) << endl;


		//	if (fabs(fYield_function/fYield_function_tr) > dRelTol)
			if (abs(fYield_function) > 1e-5)
			{
				fs_micromorph3D_out << "Local Delgamma Newton-Raphson algorithm did not converge" << endl;
				ExceptionT::GeneralFail(caller,"The value of Yield function is %d .", fYield_function);
			}
			else
			{
				fState_variables_IPs(IP,kDelgamma) = fDelgamma;
				fState_variables_IPs(IP,kc) = cohesion;
				fs_micromorph3D_out << "fDelgamma = " << fDelgamma << endl;
			}




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

                            /* Total Lagraingian strain tensor is saved as Eulerian strain tensor for plotting purposes */

                            Extract_six_values_from_symmetric_tensor(fEulerian_strain_tensor_current_IP,fTemp_nine_values);
                             /* Save Eulerian strain tensor of the current IP */
                            fEulerian_strain_IPs.SetRow(IP,fTemp_nine_values);


                                if(fNorm_dev_Cauchy_stress_tensor_current_IP==0.0)
                                {
                                	fdGdCauchy_Stress=0.0;
                                }
                                else
                                {
                                    /* calculate stress derivative of yield function */
                                	fdFYdCauchy_Stress = 0.0;
                                	fdFYdCauchy_Stress.SetToScaled(Bphi*1/3,fIdentity_matrix);
                                    fTemp_matrix_nsd_x_nsd.SetToScaled(1/fNorm_dev_Cauchy_stress_tensor_current_IP,fdev_Cauchy_stress_tensor_current_IP);
                                    fdFYdCauchy_Stress+=fTemp_matrix_nsd_x_nsd;
                                    /* calculate stress derivative of plastic potential function */
                                    fdGdCauchy_Stress = 0.0;
                                    fdGdCauchy_Stress.SetToScaled(Bpsi*1/3,fIdentity_matrix);
                                    fTemp_matrix_nsd_x_nsd.SetToScaled(1/fNorm_dev_Cauchy_stress_tensor_current_IP,fdev_Cauchy_stress_tensor_current_IP);
                                    fdGdCauchy_Stress+=fTemp_matrix_nsd_x_nsd;
                                }

                            fdFYdc=-Aphi;
                            fTemp_matrix_nsd_x_nsd.SetToScaled(fdGdCauchy_Stress_tr_trace,fCauchy_stress_tensor_current_n_IP);
                            fTemp_matrix_nsd_x_nsd*= -1;

                            fTemp_matrix_nsd_x_nsd2.MultAB(fdGdCauchy_Stress_tr,fCauchy_stress_tensor_current_n_IP);
                            fTemp_matrix_nsd_x_nsd+= fTemp_matrix_nsd_x_nsd2;

                            fTemp_matrix_nsd_x_nsd2.MultAB(fCauchy_stress_tensor_current_n_IP,fdGdCauchy_Stress_tr_Transpose);
                            fTemp_matrix_nsd_x_nsd+= fTemp_matrix_nsd_x_nsd2;

                            fTemp_matrix_nsd_x_nsd2.SetToScaled((fMaterial_Params[kLambda]*fdGdCauchy_Stress_tr_trace),fIdentity_matrix);
                            fTemp_matrix_nsd_x_nsd+= fTemp_matrix_nsd_x_nsd2;

                            fTemp_matrix_nsd_x_nsd2.SetToScaled(fMaterial_Params[kMu],fdGdCauchy_Stress_tr);
                            fTemp_matrix_nsd_x_nsd+= fTemp_matrix_nsd_x_nsd2;

                            fTemp_matrix_nsd_x_nsd2.SetToScaled(fMaterial_Params[kMu],fdGdCauchy_Stress_tr_Transpose);
                            fTemp_matrix_nsd_x_nsd+= fTemp_matrix_nsd_x_nsd2;

                            fTemp_matrix_one_x_one=dMatrixT::Dot(fdFYdCauchy_Stress,fTemp_matrix_nsd_x_nsd);


                            Comp11 = -1*fdFYdc*fMaterial_Params[kHc]*fState_variables_n_IPs(IP,khc);
                            Comp11+= fTemp_matrix_one_x_one;
                            Comp11 = 1.0/Comp11;

          ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				/////////////////////////////////dev Cauchy Stress_n////////////////////////////////////////////////
                            mean_Cauchy_stress_n = fCauchy_stress_tensor_current_n_IP.Trace()/3;
                            fTemp_matrix_nsd_x_nsd.SetToScaled(mean_Cauchy_stress_n,fIdentity_matrix);
                            fdev_Cauchy_stress_tensor_current_n_IP = fCauchy_stress_tensor_current_n_IP;
                            fdev_Cauchy_stress_tensor_current_n_IP-= fTemp_matrix_nsd_x_nsd;
/////////////////////////////////////////////////////////////////////////////////////////////////
                             //////////////////////////////////////////////////////////////////////////
               ////////////////////////////// dev(sigma_n):dev(sigma_tr)///////////////////////////////////////////////////
                            Comp22=dMatrixT::Dot(fdev_Cauchy_stress_tensor_current_n_IP,dev_Cauchy_stress_tr);
				///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// (dev(sigma_tr):dev(sigma_tr))^3/2/////////////////////////////////////////////////
                            fTemp_matrix_one_x_one = fNorm_dev_Cauchy_stress_tensor_current_IP_tr;
                            fTemp_matrix_one_x_one*= fNorm_dev_Cauchy_stress_tensor_current_IP_tr;
                            fTemp_matrix_one_x_one*= fNorm_dev_Cauchy_stress_tensor_current_IP_tr;
                            Comp33 = fTemp_matrix_one_x_one;
                            Comp33 = 1.0/Comp33;
                            Comp44 = fNorm_dev_Cauchy_stress_tensor_current_IP_tr;
                            Comp44 = 1.0/Comp44;


                            // saving Fp, Ce, dG/dS and dF/dS, etc. at each IP of the current element
                            fdGdCauchy_Stress_IPs.SetRow(IP,fdGdCauchy_Stress);
                            fDeformation_Gradient_IPs.SetRow(IP,fDeformation_Gradient);
                            //Extract_six_values_from_symmetric_tensor(fCauchy_stress_tensor_current_IP,fTemp_nine_values);
                            fCauchy_stress_IPs.SetRow(IP,fCauchy_stress_tensor_current_IP);
                            //fCauchy_stress_IPs.SetRow(IP,fTemp_nine_values);
                            fdFYdCauchy_Stress_IPs.SetRow(IP,fdFYdCauchy_Stress);


                            Form_fV1p();
                            fShapeDisplGrad.MultTx(fV1p,Vintp_1_temp);
                            scale=scale_const*J;
                            Vintp_1_temp*= scale;
                            Vintp_1+= Vintp_1_temp;


                            Form_IJe_1();
                            Form_IJe_2();
                            Form_IJe_3();
                            Form_IJe_4();
                            Form_IJe_5();
                            Form_IJe_6();
                            Form_IJe_7();
                            Form_IJe_8();



                            Form_I1p_1();
                            Form_I1p_2();
                            Form_I1p_3();
                            Form_I1p_4();
                            Form_I1p_5();
                            Form_I1p_6();


                            Form_I2p_1();
                            Form_I2p_2();
                            Form_I2p_3();
                            Form_I2p_4();
                            Form_I2p_5();
                            Form_I2p_6();


                            Form_I3p_1();
                            Form_I3p_2();
                            Form_I3p_3();
                            Form_I3p_4();
                            Form_I3p_5();
                            Form_I3p_6();


                            Form_I4p_1();
                            Form_I4p_2();
                            Form_I4p_3();
                            Form_I4p_4();
                            Form_I4p_5();
                            Form_I4p_6();


                            Form_I5p_1();
                            Form_I5p_2();
                            Form_I5p_3();
                            Form_I5p_4();
                            Form_I5p_5();
                            Form_I5p_6();

                            Form_I6p_1();
                            Form_I6p_2();
                            Form_I6p_3();
                            Form_I6p_4();
                            Form_I6p_5();
                            Form_I6p_6();

                            Form_I7p_1();
                            Form_I7p_2();
                            Form_I7p_3();
                            Form_I7p_4();
                            Form_I7p_5();
                            Form_I7p_6();
                            Form_I7p_7();
                            Form_I7p_8();
                            Form_I7p_9();
                            Form_I7p_10();
                            Form_I7p_11();
                            Form_I7p_12();
                            Form_I7p_13();

                            Form_I8p_1();
                            Form_I8p_2();
                            Form_I8p_3();
                            Form_I8p_4();
                            Form_I8p_5();
                            Form_I8p_6();
                            Form_I8p_7();
                            Form_I8p_8();
                            Form_I8p_9();
                            Form_I8p_10();
                            Form_I8p_11();
                            Form_I8p_12();
                            Form_I8p_13();

                            Form_I9p_1();
                            Form_I9p_2();
                            Form_I9p_3();
                            Form_I9p_4();
                            Form_I9p_5();
                            Form_I9p_6();
                            Form_I9p_7();
                            Form_I9p_8();
                            Form_I9p_9();
                            Form_I9p_10();
                            Form_I9p_11();
                            Form_I9p_12();
                            Form_I9p_13();

                            Form_I10p_1();
                            Form_I10p_2();
                            Form_I10p_3();
                            Form_I10p_4();
                            Form_I10p_5();
                            Form_I10p_6();
                            Form_I10p_7();
                            Form_I10p_8();
                            Form_I10p_9();
                            Form_I10p_10();
                            Form_I10p_11();
                            Form_I10p_12();
                            Form_I10p_13();

/*
                            fs_micromorph3D_out << "I9p_1 = " << I9p_1 << endl;
                            fs_micromorph3D_out << "I9p_2 = " << I9p_2 << endl;
                            fs_micromorph3D_out << "I9p_3 = " << I9p_3 << endl;
                            fs_micromorph3D_out << "I9p_4 = " << I9p_4 << endl;
                            fs_micromorph3D_out << "I9p_5 = " << I9p_5 << endl;
                            fs_micromorph3D_out << "I9p_6 = " << I9p_6 << endl;
                            fs_micromorph3D_out << "I9p_7 = " << I9p_7 << endl;
                            fs_micromorph3D_out << "I9p_8 = " << I9p_8 << endl;
                            fs_micromorph3D_out << "I9p_9 = " << I9p_9 << endl;
                            fs_micromorph3D_out << "I9p_10 = " << I9p_10 << endl;
                            fs_micromorph3D_out << "I9p_11 = " << I9p_11 << endl;
                            fs_micromorph3D_out << "I9p_12 = " << I9p_12 << endl;
                            fs_micromorph3D_out << "I9p_13 = " << I9p_13 << endl;

                            fs_micromorph3D_out << "I10p_1 = " << I10p_1 << endl;
                            fs_micromorph3D_out << "I10p_2 = " << I10p_2 << endl;
                            fs_micromorph3D_out << "I10p_3 = " << I10p_3 << endl;
                            fs_micromorph3D_out << "I10p_4 = " << I10p_4 << endl;
                            fs_micromorph3D_out << "I10p_5 = " << I10p_5 << endl;
                            fs_micromorph3D_out << "I10p_6 = " << I10p_6 << endl;
                            fs_micromorph3D_out << "I10p_7 = " << I10p_7 << endl;
                            fs_micromorph3D_out << "I10p_8 = " << I10p_8 << endl;
                            fs_micromorph3D_out << "I10p_9 = " << I10p_9 << endl;
                            fs_micromorph3D_out << "I10p_10 = " << I10p_10 << endl;
                            fs_micromorph3D_out << "I10p_11 = " << I10p_11 << endl;
                            fs_micromorph3D_out << "I10p_12 = " << I10p_12 << endl;
                            fs_micromorph3D_out << "I10p_13 = " << I10p_13 << endl;

*/


                            Form_I_temp_11p_1();
                            Form_I_temp_11p_2();
                            Form_I_temp_11p_3();
                            Form_I_temp_11p_4();
                            Form_I_temp_11p_5();
                            Form_I_temp_11p_6();
                            Form_I_temp_11p_7();
                            Form_I_temp_11p_8();
                            Form_I_temp_11p_9();
                            Form_I_temp_11p_10();
                            Form_I_temp_11p_11();
                            Form_I_temp_11p_12();
                            Form_I_temp_11p_13();

     /*
                            I_temp_11p_1*= -1;
                            I_temp_11p_3*= -1.0/3.0;
                            I_temp_11p_5*= -1.0/3.0;
                            I_temp_11p_6*= fMaterial_Params[kMu];
                            I_temp_11p_7*= -2.0/3.0;
                            I_temp_11p_7*= fMaterial_Params[kMu];
                            I_temp_11p_8*= fMaterial_Params[kMu];
                            I_temp_11p_10*= -1;
                            I_temp_11p_11*= -1;
                            I_temp_11p_12*= -1;
                            I_temp_11p_12*= fMaterial_Params[kMu];
                            I_temp_11p_13*= -1;
                            I_temp_11p_13*= fMaterial_Params[kMu];


                            I11p_1 = I_temp_11p_1;
                            I11p_1+= I_temp_11p_2;
                            I11p_1+= I_temp_11p_3;
                            I11p_1+= I_temp_11p_4;
                            I11p_1+= I_temp_11p_5;
                            I11p_1+= I_temp_11p_6;
                            I11p_1+= I_temp_11p_7;
                            I11p_1+= I_temp_11p_8;



                            I11p_1*= Comp44;
                            fs_micromorph3D_out << "I11p_1 = " << I11p_1 << endl;
                            fs_micromorph3D_out << "Comp44 = " << Comp44 << endl;
                            I11p_2 = I_temp_11p_9;
                            I11p_2+= I_temp_11p_10;
                            I11p_2+= I_temp_11p_11;
                            I11p_2+= I_temp_11p_12;
                            I11p_2+= I_temp_11p_13;
                            I11p_2*= Comp33;
                            I11p_1+= I11p_2;



                            Form_I12p_1();
                            Form_I12p_2();
                            Form_I12p_3();
                            Form_I12p_4();
                            Form_I12p_5();
                            Form_I12p_6();
*/


                            Form_I13p_1();
                            Form_I13p_2();
                            Form_I13p_3();
                            Form_I13p_4();
                            Form_I13p_5();
                            Form_I13p_6();

                            Form_II_temp_11p_1_1();
                            Form_II_temp_11p_1_2();
                            Form_II_temp_11p_1_3();
                            Form_II_temp_11p_1_4();
                            Form_II_temp_11p_1_5();
                            Form_II_temp_11p_1_6();

                            Form_II_temp_11p_2_1();
                            Form_II_temp_11p_2_2();
                            Form_II_temp_11p_2_3();
                            Form_II_temp_11p_2_4();
                            Form_II_temp_11p_2_5();
                            Form_II_temp_11p_2_6();

                            Form_II_temp_11p_3_1();
                            Form_II_temp_11p_3_2();
                            Form_II_temp_11p_3_3();
                            Form_II_temp_11p_3_4();
                            Form_II_temp_11p_3_5();
                            Form_II_temp_11p_3_6();

                            Form_II_temp_11p_4_1();
                            Form_II_temp_11p_4_2();
                            Form_II_temp_11p_4_3();
                            Form_II_temp_11p_4_4();
                            Form_II_temp_11p_4_5();
                            Form_II_temp_11p_4_6();

                            Form_II_temp_11p_5_1();
                            Form_II_temp_11p_5_2();
                            Form_II_temp_11p_5_3();
                            Form_II_temp_11p_5_4();
                            Form_II_temp_11p_5_5();
                            Form_II_temp_11p_5_6();

                            Form_II_temp_11p_6_1();
                            Form_II_temp_11p_6_2();
                            Form_II_temp_11p_6_3();
                            Form_II_temp_11p_6_4();
                            Form_II_temp_11p_6_5();
                            Form_II_temp_11p_6_6();

                            Form_II_temp_11p_7_1();
                            Form_II_temp_11p_7_2();
                            Form_II_temp_11p_7_3();
                            Form_II_temp_11p_7_4();
                            Form_II_temp_11p_7_5();
                            Form_II_temp_11p_7_6();

                            Form_II_temp_11p_8_1();
                            Form_II_temp_11p_8_2();
                            Form_II_temp_11p_8_3();
                            Form_II_temp_11p_8_4();
                            Form_II_temp_11p_8_5();
                            Form_II_temp_11p_8_6();

                            Form_II_temp_11p_9_1();
                            Form_II_temp_11p_9_2();
                            Form_II_temp_11p_9_3();
                            Form_II_temp_11p_9_4();
                            Form_II_temp_11p_9_5();
                            Form_II_temp_11p_9_6();

                            Form_II_temp_11p_10_1();
                            Form_II_temp_11p_10_2();
                            Form_II_temp_11p_10_3();
                            Form_II_temp_11p_10_4();
                            Form_II_temp_11p_10_5();
                            Form_II_temp_11p_10_6();

                            Form_II_temp_11p_11_1();
                            Form_II_temp_11p_11_2();
                            Form_II_temp_11p_11_3();
                            Form_II_temp_11p_11_4();
                            Form_II_temp_11p_11_5();
                            Form_II_temp_11p_11_6();

                            Form_II_temp_11p_12_1();
                            Form_II_temp_11p_12_2();
                            Form_II_temp_11p_12_3();
                            Form_II_temp_11p_12_4();
                            Form_II_temp_11p_12_5();
                            Form_II_temp_11p_12_6();

                            Form_II_temp_11p_13_1();
                            Form_II_temp_11p_13_2();
                            Form_II_temp_11p_13_3();
                            Form_II_temp_11p_13_4();
                            Form_II_temp_11p_13_5();
                            Form_II_temp_11p_13_6();


/*
                            ////////////////////////////comparison///////////////////////////////
                            Form_I_temp_11p_test_2_1_1();
                            Form_I_temp_11p_test_2_1_2();
                            Form_I_temp_11p_test_2_1_3();
                            Form_I_temp_11p_test_2_1_4();

                            Form_I_temp_11p_test_2_2_1();
                            Form_I_temp_11p_test_2_2_2();
                            Form_I_temp_11p_test_2_2_3();
                            Form_I_temp_11p_test_2_2_4();


                            Form_I_temp_11p_test_2_3_1();
                            Form_I_temp_11p_test_2_3_2();


                            Form_I_temp_11p_test_2_4_1();
                            Form_I_temp_11p_test_2_4_2();
                            Form_I_temp_11p_test_2_4_3();
                            Form_I_temp_11p_test_2_4_4();


                            Form_I_temp_11p_test_2_5_1();
                            Form_I_temp_11p_test_2_5_2();


                            Form_I_temp_11p_test_2_6_1();
                            Form_I_temp_11p_test_2_6_2();
                            Form_I_temp_11p_test_2_6_3();
                            Form_I_temp_11p_test_2_6_4();


                            Form_I_temp_11p_test_2_7_1();
                            Form_I_temp_11p_test_2_7_2();


                            Form_I_temp_11p_test_2_8_1();
                            Form_I_temp_11p_test_2_8_2();
                            Form_I_temp_11p_test_2_8_3();
                            Form_I_temp_11p_test_2_8_4();

                            Form_I_temp_11p_test_2_9_1();
                            Form_I_temp_11p_test_2_9_2();
                            Form_I_temp_11p_test_2_9_3();
                            Form_I_temp_11p_test_2_9_4();

                            Form_I_temp_11p_test_2_10_1();
                            Form_I_temp_11p_test_2_10_2();
                            Form_I_temp_11p_test_2_10_3();
                            Form_I_temp_11p_test_2_10_4();

                            Form_I_temp_11p_test_2_11_1();
                            Form_I_temp_11p_test_2_11_2();
                            Form_I_temp_11p_test_2_11_3();
                            Form_I_temp_11p_test_2_11_4();

                            Form_I_temp_11p_test_2_12_1();
                            Form_I_temp_11p_test_2_12_2();
                            Form_I_temp_11p_test_2_12_3();
                            Form_I_temp_11p_test_2_12_4();

                            Form_I_temp_11p_test_2_13_1();
                            Form_I_temp_11p_test_2_13_2();
                            Form_I_temp_11p_test_2_13_3();
                            Form_I_temp_11p_test_2_13_4();
*/
/*
                            I_temp_11p_test_2_1_1+= I_temp_11p_test_2_1_2;
                            I_temp_11p_test_2_1_1+= I_temp_11p_test_2_1_3;
                            I_temp_11p_test_2_1_1+= I_temp_11p_test_2_1_4;

                            I_temp_11p_test_2_2_1+= I_temp_11p_test_2_2_2;
                            I_temp_11p_test_2_2_1+= I_temp_11p_test_2_2_3;
                            I_temp_11p_test_2_2_1+= I_temp_11p_test_2_2_4;

                            I_temp_11p_test_2_3_1+= I_temp_11p_test_2_3_2;

                            I_temp_11p_test_2_4_1+= I_temp_11p_test_2_4_2;
                            I_temp_11p_test_2_4_1+= I_temp_11p_test_2_4_3;
                            I_temp_11p_test_2_4_1+= I_temp_11p_test_2_4_4;

                            I_temp_11p_test_2_5_1+= I_temp_11p_test_2_5_2;

                            I_temp_11p_test_2_6_1+= I_temp_11p_test_2_6_2;
                            I_temp_11p_test_2_6_1+= I_temp_11p_test_2_6_3;
                            I_temp_11p_test_2_6_1+= I_temp_11p_test_2_6_4;

                            I_temp_11p_test_2_7_1+= I_temp_11p_test_2_7_2;

                            I_temp_11p_test_2_8_1+= I_temp_11p_test_2_8_2;
                            I_temp_11p_test_2_8_1+= I_temp_11p_test_2_8_3;
                            I_temp_11p_test_2_8_1+= I_temp_11p_test_2_8_4;

                            fs_micromorph3D_out << "II_temp_11p_1_2 = " << II_temp_11p_1_2 << endl;
                            fs_micromorph3D_out << "I_temp_11p_test_2_1_1 = " << I_temp_11p_test_2_1_1 << endl;

                            I_temp_11p_test_2_1_1+= I_temp_11p_test_2_2_1;
                            I_temp_11p_test_2_1_1+= I_temp_11p_test_2_3_1;
                            I_temp_11p_test_2_1_1+= I_temp_11p_test_2_4_1;
                            I_temp_11p_test_2_1_1+= I_temp_11p_test_2_5_1;
                            I_temp_11p_test_2_1_1+= I_temp_11p_test_2_6_1;
                            I_temp_11p_test_2_1_1+= I_temp_11p_test_2_7_1;
                            I_temp_11p_test_2_1_1+= I_temp_11p_test_2_8_1;






       //                     I12p_2+= I13p_2;





                            fs_micromorph3D_out << "II_temp_11p_2_2 = " << II_temp_11p_2_2 << endl;
                            fs_micromorph3D_out << "I_temp_11p_test_2_2_1 = " << I_temp_11p_test_2_2_1 << endl;

                            fs_micromorph3D_out << "II_temp_11p_3_2 = " << II_temp_11p_3_2 << endl;
                            fs_micromorph3D_out << "I_temp_11p_test_2_3_1 = " << I_temp_11p_test_2_3_1 << endl;

                            fs_micromorph3D_out << "II_temp_11p_4_2 = " << II_temp_11p_4_2 << endl;
                            fs_micromorph3D_out << "I_temp_11p_test_2_4_1 = " << I_temp_11p_test_2_4_1 << endl;

                            fs_micromorph3D_out << "II_temp_11p_5_2 = " << II_temp_11p_5_2 << endl;
                            fs_micromorph3D_out << "I_temp_11p_test_2_5_1 = " << I_temp_11p_test_2_5_1 << endl;

                            fs_micromorph3D_out << "II_temp_11p_6_2 = " << II_temp_11p_6_2 << endl;
                            fs_micromorph3D_out << "I_temp_11p_test_2_6_1 = " << I_temp_11p_test_2_6_1 << endl;

                            fs_micromorph3D_out << "II_temp_11p_7_2 = " << II_temp_11p_7_2 << endl;
                            fs_micromorph3D_out << "I_temp_11p_test_2_7_1 = " << I_temp_11p_test_2_7_1 << endl;

                            fs_micromorph3D_out << "II_temp_11p_8_2 = " << II_temp_11p_8_2 << endl;
                            fs_micromorph3D_out << "I_temp_11p_test_2_8_1 = " << I_temp_11p_test_2_8_1 << endl;


                            fs_micromorph3D_out << "I12p_2 = " << I12p_2 << endl;
                            fs_micromorph3D_out << "I_temp_11p_test_2_1_1 = " << I_temp_11p_test_2_1_1 << endl;

*/
                            //////////////////////////////////////////Test///////////////////////////////////////////
/*                            Form_I_temp_11p_test_1();
                            Form_I_temp_11p_test_2();
                            Form_I_temp_11p_test_3();
                            Form_I_temp_11p_test_4();
                            Form_I12p_test_2();

                            I_temp_11p_test_4*= fMaterial_Params[kMu];
                            I_temp_11p_test_3*= fMaterial_Params[kMu];
                            I_temp_11p_test_4+= I_temp_11p_test_3;
                            I_temp_11p_test_4+= I_temp_11p_test_2;
                            I_temp_11p_test_4+= I_temp_11p_test_1;


                            I12p_test_4 = I12p_test_3;
                            I12p_test_4-= I_temp_11p_test_4;


                            fs_micromorph3D_out << "I_temp_11p_test_4 = " << I_temp_11p_test_4 << endl;
                            fs_micromorph3D_out << "I12p_test_3 = " << I12p_test_3 << endl;
                            fs_micromorph3D_out << "I12p_test_4 = " << I12p_test_4 << endl;
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////Elastic Part of Consistant Tangent////////////////////////////////////////////////////
                            /////////////////////////////////////////
                            /////////////////////////////////////////

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,IJe_1,fShapeDisplGrad);
                            scale = scale_const*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_IJe_1+= fTemp_matrix_nudof_x_nudof;


                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,IJe_2,fShapeDisplGrad);
                            scale = -1*scale_const*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_IJe_2+= fTemp_matrix_nudof_x_nudof;


                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,IJe_3,fShapeDisplGrad);
                            scale = -1*scale_const*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_IJe_3+= fTemp_matrix_nudof_x_nudof;


                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,IJe_4,fShapeDisplGrad);
                            scale = scale_const*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_IJe_4+= fTemp_matrix_nudof_x_nudof;


                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,IJe_5,fShapeDisplGrad);
                            scale = scale_const*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_IJe_5+= fTemp_matrix_nudof_x_nudof;


                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,IJe_6,fShapeDisplGrad);
                            scale = scale_const*fMaterial_Params[kLambda]*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_IJe_6+= fTemp_matrix_nudof_x_nudof;


                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,IJe_7,fShapeDisplGrad);
                            scale = scale_const*fMaterial_Params[kMu]*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_IJe_7+= fTemp_matrix_nudof_x_nudof;


                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,IJe_8,fShapeDisplGrad);
                            scale = scale_const*fMaterial_Params[kMu]*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_IJe_8+= fTemp_matrix_nudof_x_nudof;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                            ///////////////////////Plastic part of Global Consistent Tangent////////////////////////
                            ////////////////////////////////////////////////////////////////////////////////////////



                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I1p_1,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I1p_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I1p_2,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I1p_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I1p_3,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I1p_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I1p_4,fShapeDisplGrad);
                            scale = scale_const*fMaterial_Params[kLambda]*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I1p_4+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I1p_5,fShapeDisplGrad);
                            scale = scale_const*fMaterial_Params[kMu]*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I1p_5+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I1p_6,fShapeDisplGrad);
                            scale = scale_const*fMaterial_Params[kMu]*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I1p_6+= fTemp_matrix_nudof_x_nudof;



                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I2p_1,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I2p_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I2p_2,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I2p_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I2p_3,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I2p_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I2p_4,fShapeDisplGrad);
                            scale = -1*scale_const*fMaterial_Params[kLambda]*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I2p_4+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I2p_5,fShapeDisplGrad);
                            scale = -1*scale_const*fMaterial_Params[kMu]*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I2p_5+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I2p_6,fShapeDisplGrad);
                            scale = -1*scale_const*fMaterial_Params[kMu]*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I2p_6+= fTemp_matrix_nudof_x_nudof;



                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I3p_1,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I3p_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I3p_2,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I3p_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I3p_3,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I3p_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I3p_4,fShapeDisplGrad);
                            scale = -1*scale_const*fMaterial_Params[kLambda]*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I3p_4+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I3p_5,fShapeDisplGrad);
                            scale = -1*scale_const*fMaterial_Params[kMu]*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I3p_5+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I3p_6,fShapeDisplGrad);
                            scale = -1*scale_const*fMaterial_Params[kMu]*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I3p_6+= fTemp_matrix_nudof_x_nudof;



                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I4p_1,fShapeDisplGrad);
                            scale = scale_const*fMaterial_Params[kLambda]*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I4p_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I4p_2,fShapeDisplGrad);
                            scale = -1*scale_const*fMaterial_Params[kLambda]*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I4p_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I4p_3,fShapeDisplGrad);
                            scale = -1*scale_const*fMaterial_Params[kLambda]*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I4p_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I4p_4,fShapeDisplGrad);
                            scale = -1*scale_const*fMaterial_Params[kLambda]*fMaterial_Params[kLambda]*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I4p_4+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I4p_5,fShapeDisplGrad);
                            scale = -1*scale_const*fMaterial_Params[kMu]*fMaterial_Params[kLambda]*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I4p_5+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I4p_6,fShapeDisplGrad);
                            scale = -1*scale_const*fMaterial_Params[kMu]*fMaterial_Params[kLambda]*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I4p_6+= fTemp_matrix_nudof_x_nudof;



                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I5p_1,fShapeDisplGrad);
                            scale = scale_const*fMaterial_Params[kMu]*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I5p_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I5p_2,fShapeDisplGrad);
                            scale = -1*scale_const*fMaterial_Params[kMu]*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I5p_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I5p_3,fShapeDisplGrad);
                            scale = -1*scale_const*fMaterial_Params[kMu]*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I5p_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I5p_4,fShapeDisplGrad);
                            scale = -1*scale_const*fMaterial_Params[kLambda]*fMaterial_Params[kMu]*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I5p_4+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I5p_5,fShapeDisplGrad);
                            scale = -1*scale_const*fMaterial_Params[kMu]*fMaterial_Params[kMu]*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I5p_5+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I5p_6,fShapeDisplGrad);
                            scale = -1*scale_const*fMaterial_Params[kMu]*fMaterial_Params[kMu]*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I5p_6+= fTemp_matrix_nudof_x_nudof;



                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I6p_1,fShapeDisplGrad);
                            scale = scale_const*fMaterial_Params[kMu]*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I6p_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I6p_2,fShapeDisplGrad);
                            scale = -1*scale_const*fMaterial_Params[kMu]*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I6p_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I6p_3,fShapeDisplGrad);
                            scale = -1*scale_const*fMaterial_Params[kMu]*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I6p_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I6p_4,fShapeDisplGrad);
                            scale = -1*scale_const*fMaterial_Params[kLambda]*fMaterial_Params[kMu]*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I6p_4+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I6p_5,fShapeDisplGrad);
                            scale = -1*scale_const*fMaterial_Params[kMu]*fMaterial_Params[kMu]*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I6p_5+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I6p_6,fShapeDisplGrad);
                            scale = -1*scale_const*fMaterial_Params[kMu]*fMaterial_Params[kMu]*(Comp11)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I6p_6+= fTemp_matrix_nudof_x_nudof;

////////////////////////////////////// Terms related to the variation of (dG/dSigmaTr)/////////////////////////////////////////////
                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I7p_1,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*(Comp44)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I7p_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I7p_2,fShapeDisplGrad);
                            scale = -1*scale_const*fDelgamma*(Comp44)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I7p_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I7p_3,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*(Comp44)*J*(1.0/3.0);
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I7p_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I7p_4,fShapeDisplGrad);
                            scale = -1*scale_const*fDelgamma*(Comp44)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I7p_4+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I7p_5,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*(Comp44)*J*(1.0/3.0);
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I7p_5+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I7p_6,fShapeDisplGrad);
                            scale = -1*scale_const*fDelgamma*fMaterial_Params[kMu]*(Comp44)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I7p_6+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I7p_7,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*fMaterial_Params[kMu]*(Comp44)*J*(2.0/3.0);
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I7p_7+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I7p_8,fShapeDisplGrad);
                            scale = -1*scale_const*fDelgamma*fMaterial_Params[kMu]*(Comp44)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I7p_8+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I7p_9,fShapeDisplGrad);
                            scale = -1*scale_const*fDelgamma*(Comp33)*J*(Comp22);
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I7p_9+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I7p_10,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*(Comp33)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I7p_10+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I7p_11,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*(Comp33)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I7p_11+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I7p_12,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*fMaterial_Params[kMu]*(Comp33)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I7p_12+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I7p_13,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*fMaterial_Params[kMu]*(Comp33)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I7p_13+= fTemp_matrix_nudof_x_nudof;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I8p_1,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*(Comp44)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I8p_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I8p_2,fShapeDisplGrad);
                            scale = -1*scale_const*fDelgamma*(Comp44)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I8p_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I8p_3,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*(Comp44)*J*(1.0/3.0);
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I8p_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I8p_4,fShapeDisplGrad);
                            scale = -1*scale_const*fDelgamma*(Comp44)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I8p_4+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I8p_5,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*(Comp44)*J*(1.0/3.0);
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I8p_5+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I8p_6,fShapeDisplGrad);
                            scale = -1*scale_const*fDelgamma*fMaterial_Params[kMu]*(Comp44)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I8p_6+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I8p_7,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*fMaterial_Params[kMu]*(Comp44)*J*(2.0/3.0);
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I8p_7+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I8p_8,fShapeDisplGrad);
                            scale = -1*scale_const*fDelgamma*fMaterial_Params[kMu]*(Comp44)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I8p_8+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I8p_9,fShapeDisplGrad);
                            scale = -1*scale_const*fDelgamma*(Comp33)*J*(Comp22);
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I8p_9+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I8p_10,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*(Comp33)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I8p_10+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I8p_11,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*(Comp33)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I8p_11+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I8p_12,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*fMaterial_Params[kMu]*(Comp33)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I8p_12+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I8p_13,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*fMaterial_Params[kMu]*(Comp33)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I8p_13+= fTemp_matrix_nudof_x_nudof;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                            ////////////////////////////////////////////////////////////////////////////
                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I9p_1,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*(Comp44)*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I9p_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I9p_2,fShapeDisplGrad);
                            scale = -1*scale_const*fDelgamma*(Comp44)*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I9p_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I9p_3,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*(Comp44)*J*(1.0/3.0)*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I9p_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I9p_4,fShapeDisplGrad);
                            scale = -1*scale_const*fDelgamma*(Comp44)*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I9p_4+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I9p_5,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*(Comp44)*J*(1.0/3.0)*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I9p_5+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I9p_6,fShapeDisplGrad);
                            scale = -1*scale_const*fDelgamma*fMaterial_Params[kMu]*(Comp44)*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I9p_6+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I9p_7,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*fMaterial_Params[kMu]*(Comp44)*J*(2.0/3.0)*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I9p_7+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I9p_8,fShapeDisplGrad);
                            scale = -1*scale_const*fDelgamma*fMaterial_Params[kMu]*(Comp44)*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I9p_8+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I9p_9,fShapeDisplGrad);
                            scale = -1*scale_const*fDelgamma*(Comp33)*J*(Comp22)*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I9p_9+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I9p_10,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*(Comp33)*fMaterial_Params[kMu]*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I9p_10+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I9p_11,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*(Comp33)*fMaterial_Params[kMu]*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I9p_11+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I9p_12,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*fMaterial_Params[kMu]*fMaterial_Params[kMu]*(Comp33)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I9p_12+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I9p_13,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*fMaterial_Params[kMu]*fMaterial_Params[kMu]*(Comp33)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I9p_13+= fTemp_matrix_nudof_x_nudof;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I10p_1,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*(Comp44)*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I10p_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I10p_2,fShapeDisplGrad);
                            scale = -1*scale_const*fDelgamma*(Comp44)*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I10p_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I10p_3,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*(Comp44)*J*(1.0/3.0)*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I10p_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I10p_4,fShapeDisplGrad);
                            scale = -1*scale_const*fDelgamma*(Comp44)*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I10p_4+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I10p_5,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*(Comp44)*J*(1.0/3.0)*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I10p_5+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I10p_6,fShapeDisplGrad);
                            scale = -1*scale_const*fDelgamma*fMaterial_Params[kMu]*(Comp44)*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I10p_6+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I10p_7,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*fMaterial_Params[kMu]*(Comp44)*J*(2.0/3.0)*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I10p_7+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I10p_8,fShapeDisplGrad);
                            scale = -1*scale_const*fDelgamma*fMaterial_Params[kMu]*(Comp44)*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I10p_8+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I10p_9,fShapeDisplGrad);
                            scale = -1*scale_const*fDelgamma*(Comp33)*J*(Comp22)*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I10p_9+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I10p_10,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*(Comp33)*fMaterial_Params[kMu]*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I10p_10+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I10p_11,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*(Comp33)*fMaterial_Params[kMu]*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I10p_11+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I10p_12,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*fMaterial_Params[kMu]*fMaterial_Params[kMu]*(Comp33)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I10p_12+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I10p_13,fShapeDisplGrad);
                            scale = scale_const*fDelgamma*fMaterial_Params[kMu]*fMaterial_Params[kMu]*(Comp33)*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I10p_13+= fTemp_matrix_nudof_x_nudof;


							////////////////////////////////////////////////////////////////////////////
                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I12p_1,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I12p_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I12p_2,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I12p_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I12p_3,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I12p_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I12p_4,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*fMaterial_Params[kLambda]*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I12p_4+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I12p_5,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*fMaterial_Params[kMu]*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I12p_5+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I12p_6,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*fMaterial_Params[kMu]*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I12p_6+= fTemp_matrix_nudof_x_nudof;


                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I13p_1,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I13p_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I13p_2,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I13p_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I13p_3,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I13p_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I13p_4,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fMaterial_Params[kLambda]*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I13p_4+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I13p_5,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fMaterial_Params[kMu]*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I13p_5+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I13p_6,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fMaterial_Params[kMu]*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I13p_6+= fTemp_matrix_nudof_x_nudof;



                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_1_1,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_1_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_1_2,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_1_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_1_3,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_1_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_1_4,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fMaterial_Params[kLambda]*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_1_4+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_1_5,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fMaterial_Params[kMu]*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_1_5+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_1_6,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fMaterial_Params[kMu]*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_1_6+= fTemp_matrix_nudof_x_nudof;


                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_2_1,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_2_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_2_2,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_2_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_2_3,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_2_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_2_4,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fMaterial_Params[kLambda]*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_2_4+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_2_5,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fMaterial_Params[kMu]*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_2_5+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_2_6,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fMaterial_Params[kMu]*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_2_6+= fTemp_matrix_nudof_x_nudof;


                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_3_1,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fDelgamma*J*1.0/3.0;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_3_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_3_2,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fDelgamma*J*1.0/3.0;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_3_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_3_3,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fDelgamma*J*1.0/3.0;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_3_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_3_4,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fMaterial_Params[kLambda]*fDelgamma*J*1.0/3.0;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_3_4+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_3_5,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fMaterial_Params[kMu]*fDelgamma*J*1.0/3.0;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_3_5+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_3_6,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fMaterial_Params[kMu]*fDelgamma*J*1.0/3.0;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_3_6+= fTemp_matrix_nudof_x_nudof;



                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_4_1,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_4_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_4_2,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_4_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_4_3,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_4_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_4_4,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fMaterial_Params[kLambda]*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_4_4+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_4_5,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fMaterial_Params[kMu]*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_4_5+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_4_6,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fMaterial_Params[kMu]*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_4_6+= fTemp_matrix_nudof_x_nudof;



                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_5_1,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fDelgamma*J*1.0/3.0;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_5_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_5_2,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fDelgamma*J*1.0/3.0;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_5_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_5_3,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fDelgamma*J*1.0/3.0;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_5_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_5_4,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fMaterial_Params[kLambda]*fDelgamma*J*1.0/3.0;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_5_4+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_5_5,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fMaterial_Params[kMu]*fDelgamma*J*1.0/3.0;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_5_5+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_5_6,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fMaterial_Params[kMu]*fDelgamma*J*1.0/3.0;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_5_6+= fTemp_matrix_nudof_x_nudof;



                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_6_1,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_6_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_6_2,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_6_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_6_3,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_6_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_6_4,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fMaterial_Params[kLambda]*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_6_4+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_6_5,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fMaterial_Params[kMu]*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_6_5+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_6_6,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fMaterial_Params[kMu]*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_6_6+= fTemp_matrix_nudof_x_nudof;



                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_7_1,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fDelgamma*J*2.0/3.0*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_7_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_7_2,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fDelgamma*J*2.0/3.0*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_7_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_7_3,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fDelgamma*J*2.0/3.0*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_7_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_7_4,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fMaterial_Params[kLambda]*fDelgamma*J*2.0/3.0*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_7_4+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_7_5,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fMaterial_Params[kMu]*fDelgamma*J*2.0/3.0*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_7_5+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_7_6,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fMaterial_Params[kMu]*fDelgamma*J*2.0/3.0*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_7_6+= fTemp_matrix_nudof_x_nudof;



                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_8_1,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_8_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_8_2,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_8_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_8_3,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_8_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_8_4,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fMaterial_Params[kLambda]*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_8_4+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_8_5,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fMaterial_Params[kMu]*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_8_5+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_8_6,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fMaterial_Params[kMu]*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_8_6+= fTemp_matrix_nudof_x_nudof;



                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_9_1,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_9_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_9_2,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp33)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_9_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_9_3,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp33)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_9_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_9_4,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp33)*fMaterial_Params[kLambda]*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_9_4+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_9_5,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp33)*fMaterial_Params[kMu]*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_9_5+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_9_6,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp33)*fMaterial_Params[kMu]*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_9_6+= fTemp_matrix_nudof_x_nudof;



                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_10_1,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp33)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_10_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_10_2,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_10_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_10_3,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_10_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_10_4,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fMaterial_Params[kLambda]*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_10_4+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_10_5,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fMaterial_Params[kMu]*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_10_5+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_10_6,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fMaterial_Params[kMu]*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_10_6+= fTemp_matrix_nudof_x_nudof;



                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_11_1,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp33)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_11_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_11_2,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_11_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_11_3,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_11_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_11_4,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fMaterial_Params[kLambda]*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_11_4+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_11_5,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fMaterial_Params[kMu]*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_11_5+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_11_6,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fMaterial_Params[kMu]*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_11_6+= fTemp_matrix_nudof_x_nudof;



                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_12_1,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp33)*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_12_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_12_2,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_12_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_12_3,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_12_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_12_4,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fMaterial_Params[kLambda]*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_12_4+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_12_5,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fMaterial_Params[kMu]*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_12_5+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_12_6,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fMaterial_Params[kMu]*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_12_6+= fTemp_matrix_nudof_x_nudof;



                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_13_1,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp33)*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_13_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_13_2,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_13_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_13_3,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_13_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_13_4,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fMaterial_Params[kLambda]*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_13_4+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_13_5,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fMaterial_Params[kMu]*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_13_5+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,II_temp_11p_13_6,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fMaterial_Params[kMu]*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_13_6+= fTemp_matrix_nudof_x_nudof;

/*
//////////////////////////////////////constructing full implementation//////////////////////////////////////////////////////
                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_1_1,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_1_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_1_2,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_1_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_1_3,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_1_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_1_4,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_1_4+= fTemp_matrix_nudof_x_nudof;


                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_2_1,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_2_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_2_2,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_2_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_2_3,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_2_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_2_4,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_2_4+= fTemp_matrix_nudof_x_nudof;


                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_3_1,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fDelgamma*J*(1.0/3.0);
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_3_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_3_2,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fDelgamma*J*(1.0/3.0);
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_3_2+= fTemp_matrix_nudof_x_nudof;


                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_4_1,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_4_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_4_2,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_4_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_4_3,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_4_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_4_4,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_4_4+= fTemp_matrix_nudof_x_nudof;


                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_5_1,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fDelgamma*J*(1.0/3.0);
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_5_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_5_2,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fDelgamma*J*(1.0/3.0);
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_5_2+= fTemp_matrix_nudof_x_nudof;


                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_6_1,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_6_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_6_2,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_6_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_6_3,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_6_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_6_4,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_6_4+= fTemp_matrix_nudof_x_nudof;


                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_7_1,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fDelgamma*J*(2.0/3.0)*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_7_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_7_2,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp44)*fDelgamma*J*(2.0/3.0)*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_7_2+= fTemp_matrix_nudof_x_nudof;


                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_8_1,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_8_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_8_2,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_8_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_8_3,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_8_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_8_4,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp44)*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_8_4+= fTemp_matrix_nudof_x_nudof;


                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_9_1,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp33)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_9_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_9_2,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp33)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_9_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_9_3,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp33)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_9_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_9_4,fShapeDisplGrad);
                            scale = scale_const*(Comp11)*(Comp33)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_9_4+= fTemp_matrix_nudof_x_nudof;


                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_10_1,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_10_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_10_2,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_10_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_10_3,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_10_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_10_4,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_10_4+= fTemp_matrix_nudof_x_nudof;


                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_11_1,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_11_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_11_2,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_11_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_11_3,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_11_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_11_4,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fDelgamma*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_11_4+= fTemp_matrix_nudof_x_nudof;


                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_12_1,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_12_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_12_2,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_12_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_12_3,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_12_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_12_4,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_12_4+= fTemp_matrix_nudof_x_nudof;


                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_13_1,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_13_1+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_13_2,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_13_2+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_13_3,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_13_3+= fTemp_matrix_nudof_x_nudof;

                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,I_temp_11p_test_2_13_4,fShapeDisplGrad);
                            scale = -1*scale_const*(Comp11)*(Comp33)*fDelgamma*J*fMaterial_Params[kMu];
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_I_temp_11p_test_2_13_4+= fTemp_matrix_nudof_x_nudof;
*/
                      /////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                            /////////////////////////////////////////////////////////////////
                            /////////////////////////////////////////////////////////////////
                            }

                    else //(yielding did not occur / elastic step only/
                    {
                    		//fdelDelgamma = 0.0;


                    		fDelgamma = 0.0;

                    		fDel_Deformation_Gradient = fDeformation_Gradient;
                    		fDel_Deformation_Gradient-=fDeformation_Gradient_n;
        					fDeformation_Gradient_Inverse.Inverse(fDeformation_Gradient);
        					fVelocity_Gradient_current_IP.MultAB(fDel_Deformation_Gradient,fDeformation_Gradient_Inverse);

        					fCauchy_stress_tensor_current_IP = fCauchy_stress_tensor_current_IP_tr;
        					mean_Cauchy_stress = fCauchy_stress_tensor_current_IP.Trace()/3;
        					fTemp_matrix_nsd_x_nsd2.SetToScaled(mean_Cauchy_stress,fIdentity_matrix);
        					fdev_Cauchy_stress_tensor_current_IP = fCauchy_stress_tensor_current_IP;
        					fdev_Cauchy_stress_tensor_current_IP-= fTemp_matrix_nsd_x_nsd2;
        					fNorm_dev_Cauchy_stress_tensor_current_IP = fdev_Cauchy_stress_tensor_current_IP.ScalarProduct();
        					fNorm_dev_Cauchy_stress_tensor_current_IP = sqrt(fNorm_dev_Cauchy_stress_tensor_current_IP);

                            fDeformation_Gradient_IPs.SetRow(IP,fDeformation_Gradient);

                            //Extract_six_values_from_symmetric_tensor(fCauchy_stress_tensor_current_IP,fTemp_nine_values);
                           // fCauchy_stress_IPs.SetRow(IP,fTemp_nine_values);
                            fCauchy_stress_IPs.SetRow(IP,fCauchy_stress_tensor_current_IP);


                            fs_micromorph3D_out  << "Element = " << e << endl;
                            fs_micromorph3D_out  << "Gauss Point = " << IP << endl;
                            fs_micromorph3D_out  << "fDeformation_Gradient_n = " << fDeformation_Gradient_n << endl;
                            fs_micromorph3D_out  << "fDeformation_Gradient = " << fDeformation_Gradient << endl;
                            fs_micromorph3D_out  << "fCauchy_stress_tensor_current_n_IP = " << fCauchy_stress_tensor_current_n_IP << endl;
                            fs_micromorph3D_out  << "fCauchy_stress_tensor_current_IP = " << fCauchy_stress_tensor_current_IP << endl;


                            if(fYield_function_tr < dAbsTol && time > 0.0 )
								{
									/* calculate stress derivative of plastic potential function */
									fdGdCauchy_Stress = 0.0;
									fdGdCauchy_Stress.SetToScaled(Bpsi*1/3,fIdentity_matrix);
									fTemp_matrix_nsd_x_nsd.SetToScaled(1/fNorm_dev_Cauchy_stress_tensor_current_IP,fdev_Cauchy_stress_tensor_current_IP);
									fdGdCauchy_Stress+=fTemp_matrix_nsd_x_nsd;
								}

                            fdGdCauchy_Stress_IPs.SetRow(IP,fdGdCauchy_Stress);
                            Form_fV1e();
                            fShapeDisplGrad.MultTx(fV1e,Vinte_1_temp);
                            scale=scale_const*J;
                            Vinte_1_temp*=scale;
                            Vinte_1 +=Vinte_1_temp;

                            //cout << "Vintp_1 = " << Vintp_1 << endl;

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

                            /* [fRight_Cauchy_Green_tensor] will be formed */
                            fRight_Cauchy_Green_tensor.MultATB(fDeformation_Gradient,fDeformation_Gradient);
                            if (fRight_Cauchy_Green_tensor.Det()==0)
                                    fRight_Cauchy_Green_tensor = fIdentity_matrix;


                            //Extract_nine_values (fEulerian_strain_tensor_current_IP,fTemp_six_values);
                            Extract_six_values_from_symmetric_tensor(fEulerian_strain_tensor_current_IP,fTemp_nine_values);

                            /* Save Eulerian strain tensor of the current IP */
                            fEulerian_strain_IPs.SetRow(IP,fTemp_nine_values);


                            Form_IJe_1();
                            Form_IJe_2();
                            Form_IJe_3();
                            Form_IJe_4();
                            Form_IJe_5();
                            Form_IJe_6();
                            Form_IJe_7();
                            Form_IJe_8();



							fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,IJe_1,fShapeDisplGrad);
                            scale = scale_const*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_IJe_1+= fTemp_matrix_nudof_x_nudof;


                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,IJe_2,fShapeDisplGrad);
                            scale = -1*scale_const*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_IJe_2+= fTemp_matrix_nudof_x_nudof;


                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,IJe_3,fShapeDisplGrad);
                            scale = -1*scale_const*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_IJe_3+= fTemp_matrix_nudof_x_nudof;


                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,IJe_4,fShapeDisplGrad);
                            scale = scale_const*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_IJe_4+= fTemp_matrix_nudof_x_nudof;


                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,IJe_5,fShapeDisplGrad);
                            scale = scale_const*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_IJe_5+= fTemp_matrix_nudof_x_nudof;


                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,IJe_6,fShapeDisplGrad);
                            scale = scale_const*fMaterial_Params[kLambda]*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_IJe_6+= fTemp_matrix_nudof_x_nudof;


                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,IJe_7,fShapeDisplGrad);
                            scale = scale_const*fMaterial_Params[kMu]*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_IJe_7+= fTemp_matrix_nudof_x_nudof;


                            fTemp_matrix_nudof_x_nudof.MultATBC(fShapeDisplGrad,IJe_8,fShapeDisplGrad);
                            scale = scale_const*fMaterial_Params[kMu]*J;
                            fTemp_matrix_nudof_x_nudof*= scale;
                            fKu_IJe_8+= fTemp_matrix_nudof_x_nudof;




/*************************************************************************************************************************/

                    }//elastic part ends
                }// constitutive model =3 ends
                if(iConstitutiveModelType==2)
                {


                	////////////////////////////////////////////////////////////////////////////////////////
                	/////////////////MicroMorphic Internal force vectors////////////////////////////////////
                		Form_G1_matrix();//output:G1 vector & Sigma matrix
                        fIota_w_temp_matrix.Multx(G1,Uint_1_temp);
                        scale=-1*scale_const*J;
                        Uint_1_temp*=scale;
                        Uint_1 +=Uint_1_temp;
                        fs_micromorph3D_out << "Uint_1 = " << Uint_1 << endl;

                        //fShapeDispl.MultTx(fGravity_vector,Uext_1);
                        //Uext_1*=-fMaterial_Params[kRho_0];


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
                        fs_micromorph3D_out << "fCauchy_stress_tensor_current_IP = " << fCauchy_stress_tensor_current_IP << endl;



                        // extract six values of stress from symmetric cauchy stress tensor
                        Extract_six_values_from_symmetric_tensor(fCauchy_stress_tensor_current_IP,fTemp_nine_values);
                        // Extract_six_values_from_symmetric_tensor(fCauchy_stress_tensor_current_IP,fTemp_six_values);


                        //Save Cauchy effective stress tensor of the current IP
                        //fCauchy_stress_IPs.SetRow(IP,fTemp_six_values);
                        fCauchy_stress_IPs.SetRow(IP,fTemp_nine_values);



   ////////////////////////////////////////////////////////////////////////////////
   /////////////////Micromorphic Internal force vectors finish here////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////
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
                        fs_micromorph3D_out << "fG1_1 = " << fG1_1 << endl;

                        fTemp_matrix_nudof_x_nudof.MultABCT(fIota_w_temp_matrix,Tsigma_2,fIota_temp_matrix);
                        scale =-scale_const*J;
                        fTemp_matrix_nudof_x_nudof *= scale;
                        // accumulate
                        fG1_2 += fTemp_matrix_nudof_x_nudof;
                        fs_micromorph3D_out << "fG1_2 = " << fG1_2 << endl;

                        fTemp_matrix_nudof_x_nudof.MultABCT(fIota_w_temp_matrix,Tsigma_3,fIota_temp_matrix);
                        scale = -scale_const*J;
                        fTemp_matrix_nudof_x_nudof *= scale;
                        // accumulate
                        fG1_3 += fTemp_matrix_nudof_x_nudof;
                        fs_micromorph3D_out << "fG1_3 = " << fG1_3 << endl;

                        fTemp_matrix_nudof_x_nudof.MultABCT(fIota_temp_matrix,TFn_1,fIota_temp_matrix);
                        scale = -scale_const*J*(fMaterial_Params[kLambda]+fMaterial_Params[kTau]);
                        fTemp_matrix_nudof_x_nudof *= scale;
                        // accumulate
                        fG1_4 += fTemp_matrix_nudof_x_nudof;
                        fs_micromorph3D_out << "fG1_4 = " << fG1_4 << endl;

                        fTemp_matrix_nudof_x_nudof.MultABCT(fIota_w_temp_matrix,TFn_2,fIota_temp_matrix);
                        scale = -scale_const*J*(fMaterial_Params[kMu]+fMaterial_Params[kSigma_const]);
                        fTemp_matrix_nudof_x_nudof *= scale;
                        //accumulate
                        fG1_5 += fTemp_matrix_nudof_x_nudof;
                        fs_micromorph3D_out << "fG1_5 = " << fG1_5 << endl;

                        fTemp_matrix_nudof_x_nudof.MultABCT(fIota_w_temp_matrix,TFn_3,fIota_temp_matrix);
                        scale = -scale_const*J*(fMaterial_Params[kMu]+fMaterial_Params[kSigma_const]);
                        fTemp_matrix_nudof_x_nudof *= scale;
                        // accumulate
                        fG1_6 += fTemp_matrix_nudof_x_nudof;
                        fs_micromorph3D_out << "fG1_6 = " << fG1_6 << endl;

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
                        fs_micromorph3D_out << "fG1_13 = " << fG1_13 << endl;

                        //TransShapeDisplGrad.Transpose(GRAD_Nuw);//
                        //fTemp_matrix_nudof_x_nudof.MultABCT(TransShapeDisplGrad,Var_F,fIota_temp_matrix);
                        fTemp_matrix_nudof_x_nudof.MultABCT(fIota_w_temp_matrix,Var_F,fIota_temp_matrix);
                        scale= scale_const*J;
                        fTemp_matrix_nudof_x_nudof *= scale;
                        fG1_14 +=fTemp_matrix_nudof_x_nudof;
                        fs_micromorph3D_out << "fG1_14 = " << fG1_14 << endl;

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


                fState_variables_IPs(IP,ktrCauchy_Stress)=fCauchy_stress_tensor_current_IP_trace;
                fState_variables_IPs(IP,kNorm_Dev_Cauchy_Stress)=fNorm_dev_Cauchy_stress_tensor_current_IP;
                fState_variables_IPs(IP,ktrRel)=trs_sigma;
                fState_variables_IPs(IP,kRel_inv)=Rel_strs_inv;
                fState_variables_IPs(IP,ktrm)=trmklm;
                fState_variables_IPs(IP,km_inv)=Higher_orderT_inv;
                fState_variables_IPs(IP,ktrS)=trS;
                fState_variables_IPs(IP,kinvdevS)=invdevS;
                fState_variables_IPs(IP,ktrSIGMA_S)=trSIGMA_S;
                fState_variables_IPs(IP,kinvdevSIGMA_S)=invdevSIGMA_S;
                fState_variables_IPs(IP,kinvtrM)=invtrMKLM;
                fState_variables_IPs(IP,kinvdevM)=invdevMKLM;
                fState_variables_IPs(IP,kinvPhi)=invPhi;
                fState_variables_IPs(IP,kinvGPhi)=invGPhi;
                fState_variables_IPs(IP,ktreps)=treps;
                fState_variables_IPs(IP,kdeveps)=deveps;
                fState_variables_IPs(IP,kinvtrgammastn)=invtrgammastn;
                fState_variables_IPs(IP,kinvdevgammastn)=invdevgammastn;


//
        } //end Gauss integration loop




        /* saving eulerian strain for each IPs of the current element */
        fEulerian_strain_Elements_IPs.SetRow(e,fEulerian_strain_IPs);

        /* saving cauchy stress for each IPs of the current element */
        fCauchy_stress_Elements_IPs.SetRow(e,fCauchy_stress_IPs);

        // saving dGdCauchy Stress for each IP of the current element //
        fdGdCauchy_Stress_Elements_IPs.SetRow(e,fdGdCauchy_Stress_IPs);

        // saving Cauchy Stress for each IP of the current element //
        fDeformation_Gradient_Elements_IPs.SetRow(e,fDeformation_Gradient_IPs);

        // saving state variables for each IPs of the current element //
        fState_variables_Elements_IPs.SetRow(e,fState_variables_IPs);




        if(iConstitutiveModelType==2)
        {
            GammaN_IPs_el.SetRow(e,GammaN_IPs);
            SigN_IPs_el.SetRow(e,SigN_IPs);
            sn_sigman_IPs_el.SetRow(e,sn_sigman_IPs);
            mn_IPs_el.SetRow(e,mn_IPs);

            F_ar_IPs_el.SetRow(e,F_ar_IPs);
            FInv_ar_IPs_el.SetRow(e,FInv_ar_IPs);
            Chi_ar_IPs_el.SetRow(e,Chi_ar_IPs);
            GRAD_Chi_ar_IPs_el.SetRow(e,GRAD_Chi_ar_IPs);
        }

        if(iConstitutiveModelType==1)
        {

        }
        if(iConstitutiveModelType==3)
        {
        	//{fFd_int} will be formed
            fFd_int = 0.0;
            fFd_int = Vinte_1;
            fFd_int+= Vintp_1;
            fFd_int*= -1;
            //cout << "fFd_int = " << fFd_int << endl;


        //Micromorphic case fKdd coming from bal. of linear momentum
        //fKdd=0.0;
// The elastic part of global consistent tangent


            fKdd=	fKu_IJe_1;
            fKdd+=  fKu_IJe_2;
            fKdd+=  fKu_IJe_3;
            fKdd+=  fKu_IJe_4;
            fKdd+=  fKu_IJe_5;
            fKdd+=  fKu_IJe_6;
            fKdd+=  fKu_IJe_7;
            fKdd+=  fKu_IJe_8;


            /* [fKdphi] will be formed */
            fKdphi  = 0.0;
            /* [fKphid] will be formed */
            fKphid = 0.0;
            /* [fKphiphi] will be formed */
            fKphiphi = 0.0;
            /* {fFphi_int} will be formed */
            fFphi_int = 0.0;
            /* Plasticity contribution*/
            fKdd+= 	fKu_I1p_1;
            fKdd+=  fKu_I1p_2;
            fKdd+=  fKu_I1p_3;
            fKdd+=  fKu_I1p_4;
            fKdd+=  fKu_I1p_5;
            fKdd+=  fKu_I1p_6;


            fKdd+= 	fKu_I2p_1;
            fKdd+=  fKu_I2p_2;
            fKdd+=  fKu_I2p_3;
            fKdd+=  fKu_I2p_4;
            fKdd+=  fKu_I2p_5;
            fKdd+=  fKu_I2p_6;


            fKdd+= 	fKu_I3p_1;
            fKdd+=  fKu_I3p_2;
            fKdd+=  fKu_I3p_3;
            fKdd+=  fKu_I3p_4;
            fKdd+=  fKu_I3p_5;
            fKdd+=  fKu_I3p_6;


            fKdd+= 	fKu_I4p_1;
            fKdd+=  fKu_I4p_2;
            fKdd+=  fKu_I4p_3;
            fKdd+=  fKu_I4p_4;
            fKdd+=  fKu_I4p_5;
            fKdd+=  fKu_I4p_6;


            fKdd+= 	fKu_I5p_1;
            fKdd+=  fKu_I5p_2;
            fKdd+=  fKu_I5p_3;
            fKdd+=  fKu_I5p_4;
            fKdd+=  fKu_I5p_5;
            fKdd+=  fKu_I5p_6;


            fKdd+= 	fKu_I6p_1;
            fKdd+=  fKu_I6p_2;
            fKdd+=  fKu_I6p_3;
            fKdd+=  fKu_I6p_4;
            fKdd+=  fKu_I6p_5;
            fKdd+=  fKu_I6p_6;


            fKdd+= 	fKu_I7p_1;
            fKdd+=  fKu_I7p_2;
            fKdd+=  fKu_I7p_3;
            fKdd+=  fKu_I7p_4;
            fKdd+=  fKu_I7p_5;
            fKdd+=  fKu_I7p_6;
            fKdd+= 	fKu_I7p_7;
            fKdd+=  fKu_I7p_8;
            fKdd+=  fKu_I7p_9;
            fKdd+=  fKu_I7p_10;
            fKdd+=  fKu_I7p_11;
            fKdd+=  fKu_I7p_12;
            fKdd+=  fKu_I7p_13;

            fKdd+= 	fKu_I8p_1;
            fKdd+=  fKu_I8p_2;
            fKdd+=  fKu_I8p_3;
            fKdd+=  fKu_I8p_4;
            fKdd+=  fKu_I8p_5;
            fKdd+=  fKu_I8p_6;
            fKdd+= 	fKu_I8p_7;
            fKdd+=  fKu_I8p_8;
            fKdd+=  fKu_I8p_9;
            fKdd+=  fKu_I8p_10;
            fKdd+=  fKu_I8p_11;
            fKdd+=  fKu_I8p_12;
            fKdd+=  fKu_I8p_13;

            fKdd+= 	fKu_I9p_1;
            fKdd+=  fKu_I9p_2;
            fKdd+=  fKu_I9p_3;
            fKdd+=  fKu_I9p_4;
            fKdd+=  fKu_I9p_5;
            fKdd+=  fKu_I9p_6;
            fKdd+= 	fKu_I9p_7;
            fKdd+=  fKu_I9p_8;
            fKdd+=  fKu_I9p_9;
            fKdd+=  fKu_I9p_10;
            fKdd+=  fKu_I9p_11;
            fKdd+=  fKu_I9p_12;
            fKdd+=  fKu_I9p_13;

            fKdd+= 	fKu_I10p_1;
            fKdd+=  fKu_I10p_2;
            fKdd+=  fKu_I10p_3;
            fKdd+=  fKu_I10p_4;
            fKdd+=  fKu_I10p_5;
            fKdd+=  fKu_I10p_6;
            fKdd+= 	fKu_I10p_7;
            fKdd+=  fKu_I10p_8;
            fKdd+=  fKu_I10p_9;
            fKdd+=  fKu_I10p_10;
            fKdd+=  fKu_I10p_11;
            fKdd+=  fKu_I10p_12;
            fKdd+=  fKu_I10p_13;

/*
            fKdd+= 	fKu_I12p_1;
            fKdd+=  fKu_I12p_2;
            fKdd+=  fKu_I12p_3;
            fKdd+=  fKu_I12p_4;
            fKdd+=  fKu_I12p_5;
            fKdd+=  fKu_I12p_6;
*/

            fKdd+= 	fKu_I13p_1;
			fKdd+=  fKu_I13p_2;
			fKdd+=  fKu_I13p_3;
			fKdd+=  fKu_I13p_4;
			fKdd+=  fKu_I13p_5;
			fKdd+=  fKu_I13p_6;


			fKdd+=  fKu_I_temp_11p_1_1;
			fKdd+=  fKu_I_temp_11p_1_2;
			fKdd+=  fKu_I_temp_11p_1_3;
			fKdd+=  fKu_I_temp_11p_1_4;
			fKdd+=  fKu_I_temp_11p_1_5;
			fKdd+=  fKu_I_temp_11p_1_6;


			fKdd+=  fKu_I_temp_11p_2_1;
			fKdd+=  fKu_I_temp_11p_2_2;
			fKdd+=  fKu_I_temp_11p_2_3;
			fKdd+=  fKu_I_temp_11p_2_4;
			fKdd+=  fKu_I_temp_11p_2_5;
			fKdd+=  fKu_I_temp_11p_2_6;

			fKdd+=  fKu_I_temp_11p_3_1;
			fKdd+=  fKu_I_temp_11p_3_2;
			fKdd+=  fKu_I_temp_11p_3_3;
			fKdd+=  fKu_I_temp_11p_3_4;
			fKdd+=  fKu_I_temp_11p_3_5;
			fKdd+=  fKu_I_temp_11p_3_6;

			fKdd+=  fKu_I_temp_11p_4_1;
			fKdd+=  fKu_I_temp_11p_4_2;
			fKdd+=  fKu_I_temp_11p_4_3;
			fKdd+=  fKu_I_temp_11p_4_4;
			fKdd+=  fKu_I_temp_11p_4_5;
			fKdd+=  fKu_I_temp_11p_4_6;

			fKdd+=  fKu_I_temp_11p_5_1;
			fKdd+=  fKu_I_temp_11p_5_2;
			fKdd+=  fKu_I_temp_11p_5_3;
			fKdd+=  fKu_I_temp_11p_5_4;
			fKdd+=  fKu_I_temp_11p_5_5;
			fKdd+=  fKu_I_temp_11p_5_6;

			fKdd+=  fKu_I_temp_11p_6_1;
			fKdd+=  fKu_I_temp_11p_6_2;
			fKdd+=  fKu_I_temp_11p_6_3;
			fKdd+=  fKu_I_temp_11p_6_4;
			fKdd+=  fKu_I_temp_11p_6_5;
			fKdd+=  fKu_I_temp_11p_6_6;

			fKdd+=  fKu_I_temp_11p_7_1;
			fKdd+=  fKu_I_temp_11p_7_2;
			fKdd+=  fKu_I_temp_11p_7_3;
			fKdd+=  fKu_I_temp_11p_7_4;
			fKdd+=  fKu_I_temp_11p_7_5;
			fKdd+=  fKu_I_temp_11p_7_6;

            fKdd+= 	fKu_I_temp_11p_8_1;
            fKdd+=  fKu_I_temp_11p_8_2;
            fKdd+=  fKu_I_temp_11p_8_3;
            fKdd+=  fKu_I_temp_11p_8_4;
            fKdd+=  fKu_I_temp_11p_8_5;
            fKdd+=  fKu_I_temp_11p_8_6;

            fKdd+= 	fKu_I_temp_11p_9_1;
            fKdd+=  fKu_I_temp_11p_9_2;
            fKdd+=  fKu_I_temp_11p_9_3;
            fKdd+=  fKu_I_temp_11p_9_4;
            fKdd+=  fKu_I_temp_11p_9_5;
            fKdd+=  fKu_I_temp_11p_9_6;

            fKdd+= 	fKu_I_temp_11p_10_1;
            fKdd+=  fKu_I_temp_11p_10_2;
            fKdd+=  fKu_I_temp_11p_10_3;
            fKdd+=  fKu_I_temp_11p_10_4;
            fKdd+=  fKu_I_temp_11p_10_5;
            fKdd+=  fKu_I_temp_11p_10_6;

            fKdd+= 	fKu_I_temp_11p_11_1;
            fKdd+=  fKu_I_temp_11p_11_2;
            fKdd+=  fKu_I_temp_11p_11_3;
            fKdd+=  fKu_I_temp_11p_11_4;
            fKdd+=  fKu_I_temp_11p_11_5;
            fKdd+=  fKu_I_temp_11p_11_6;

            fKdd+= 	fKu_I_temp_11p_12_1;
            fKdd+=  fKu_I_temp_11p_12_2;
            fKdd+=  fKu_I_temp_11p_12_3;
            fKdd+=  fKu_I_temp_11p_12_4;
            fKdd+=  fKu_I_temp_11p_12_5;
            fKdd+=  fKu_I_temp_11p_12_6;

            fKdd+= 	fKu_I_temp_11p_13_1;
            fKdd+=  fKu_I_temp_11p_13_2;
            fKdd+=  fKu_I_temp_11p_13_3;
            fKdd+=  fKu_I_temp_11p_13_4;
            fKdd+=  fKu_I_temp_11p_13_5;
            fKdd+=  fKu_I_temp_11p_13_6;


/*

            fKdd_previous =  fKu_I13p_2;

            fKdd_previous+=  fKu_I_temp_11p_1_2;

            fKdd_previous+=  fKu_I_temp_11p_2_2;

            fKdd_previous+=  fKu_I_temp_11p_3_2;

            fKdd_previous+=  fKu_I_temp_11p_4_2;

            fKdd_previous+=  fKu_I_temp_11p_5_2;

            fKdd_previous+=  fKu_I_temp_11p_6_2;

            fKdd_previous+=  fKu_I_temp_11p_7_2;

            fKdd_previous+=  fKu_I_temp_11p_8_2;

            fKdd_previous+=  fKu_I_temp_11p_9_2;

            fKdd_previous+=  fKu_I_temp_11p_10_2;

            fKdd_previous+=  fKu_I_temp_11p_11_2;

            fKdd_previous+=  fKu_I_temp_11p_12_2;

            fKdd_previous+=  fKu_I_temp_11p_13_2;





            fKdd_full_implemet = fKu_I_temp_11p_test_2_1_1;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_1_2;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_1_3;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_1_4;

            fKdd_full_implemet+= fKu_I_temp_11p_test_2_2_1;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_2_2;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_2_3;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_2_4;

            fKdd_full_implemet+= fKu_I_temp_11p_test_2_3_1;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_3_2;

            fKdd_full_implemet+= fKu_I_temp_11p_test_2_4_1;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_4_2;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_4_3;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_4_4;

            fKdd_full_implemet+= fKu_I_temp_11p_test_2_5_1;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_5_2;

            fKdd_full_implemet+= fKu_I_temp_11p_test_2_6_1;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_6_2;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_6_3;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_6_4;

            fKdd_full_implemet+= fKu_I_temp_11p_test_2_7_1;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_7_2;


            fKdd_full_implemet+= fKu_I_temp_11p_test_2_8_1;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_8_2;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_8_3;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_8_4;


            fKdd_full_implemet+= fKu_I_temp_11p_test_2_9_1;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_9_2;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_9_3;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_9_4;

            fKdd_full_implemet+= fKu_I_temp_11p_test_2_10_1;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_10_2;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_10_3;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_10_4;

            fKdd_full_implemet+= fKu_I_temp_11p_test_2_11_1;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_11_2;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_11_3;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_11_4;

            fKdd_full_implemet+= fKu_I_temp_11p_test_2_12_1;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_12_2;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_12_3;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_12_4;

            fKdd_full_implemet+= fKu_I_temp_11p_test_2_13_1;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_13_2;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_13_3;
            fKdd_full_implemet+= fKu_I_temp_11p_test_2_13_4;

            difference = fKdd_full_implemet;
            difference-= fKdd_previous;

            fs_micromorph3D_out << "fKdd_full_implemet = " << fKdd_full_implemet << endl;
            fs_micromorph3D_out << "fKdd_previous = " << fKdd_previous << endl;
            fs_micromorph3D_out << "difference = " << difference << endl;

            //fs_micromorph3D_out << "Global Tangent = " << fKdd << endl;
*/
//******************************************************************
//******************************************************************
//******************************************************************

//******************************************************************
//*************Plastic matrices ************************************
//******************************************************************

        }
        if(iConstitutiveModelType==2)
        {

            //{fFd_int} will be formed
              fFd_int  = 0.0;
              fFd_int  = Uint_1;
              //fFd_ext =-Uext_1; //no external traction is assumed
              fFd_int *= -1;
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

          //fKdphi = 0.0;
          //Micromorphic case fKdPhi  from coming from bal. of linear momentum

          fKdphi  = fG1_7 ;
          fKdphi += fG1_9 ;
          fKdphi += fG1_11;

          //[fKphid] will be formed
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
          //fFphi_int  = 0.0;
          fFphi_int  =Pint_1;
          fFphi_int +=Pint_2;
          fFphi_int +=Pint_3;//no external traction is assumed Pext=0
          fFphi_int *= -1;
           // fFphi_int.ScalarProduct();
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
void FSMicromorphic3DCurrConfigT::SetGlobalShape(void)
{
    //          cout<<"CHECK POINT-21"<<endl;
    /* fetch (initial) coordinates */
    SetLocalX(fLocInitCoords);

    /* compute shape function derivatives */
    //fShapes_displ->SetDerivatives_DN_DDN(); Commented out for Q8P8
    fShapes_displ->SetDerivatives();
    fShapes_micro->SetDerivatives();

}



/* describe the parameters needed by the interface */
void FSMicromorphic3DCurrConfigT::DefineParameters(ParameterListT& list) const
{
     //           cout<<"CHECK POINT-22"<<endl;
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
    list.AddParameter(iPlasticityCheck, "plasticity_check");

  //  list.AddParameter(iplasticity, "plasticity");
    list.AddParameter(iIterationMax, "max_local_iterations");
    double shearMu, sLambda, Rho_0, gravity_g, gravity_g1, gravity_g2, gravity_g3;
    double Kappa, Nu, Sigma_const, Tau, Eta;
    double Tau1,Tau2,Tau3,Tau4,Tau5,Tau6,Tau7,Tau8,Tau9,Tau10,Tau11;
    double c0,Hc,Fphi,Dpsi;
    double c0_chi,Hc_chi,Fphi_chi,Dpsi_chi;
    //double Gc0_chi1,Gc0_chi2,Gc0_chi3,HGc_chi,FGphi_chi,DGpsi_chi;

    // solid elasticity
    list.AddParameter(shearMu, "mu");
    list.AddParameter(sLambda, "lambda");
    // plasticity
    list.AddParameter(c0,"c0");
    list.AddParameter(Hc,"Hc");
    list.AddParameter(Fphi,"Fphi");
    list.AddParameter(Dpsi,"Dpsi");
    // Micro-scale plasticity
    list.AddParameter(c0_chi,"c0_chi");
    list.AddParameter(Hc_chi,"Hc_chi");
    list.AddParameter(Fphi_chi,"Fphi_chi");
    list.AddParameter(Dpsi_chi,"Dpsi_chi");

    // Micro Scale Gradient plasticity
    //list.AddParameter(Gc0_chi1,"Gc0_chi1");
    //list.AddParameter(Gc0_chi2,"Gc0_chi2");
    //list.AddParameter(Gc0_chi3,"Gc0_chi3");
    //list.AddParameter(HGc_chi,"HGc_chi");
    //list.AddParameter(FGphi_chi,"FGphi_chi");
    //list.AddParameter(DGpsi_chi,"DGpsi_chi");


    // tolerance for yield check
    list.AddParameter(dYieldTrialTol, "local_yield_tr_tol");

    // convergence tolerances for local Newton-Raphson iteration
    list.AddParameter(dAbsTol, "local_tol_absolute");
    list.AddParameter(dRelTol, "local_tol_relative");
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
void FSMicromorphic3DCurrConfigT::TakeParameterList(const ParameterListT& list)
{
  //              cout<<"CHECK POINT-23"<<endl;
    const char caller[] = "FSMicromorphic3DCurrConfigT::TakeParameterList";

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
    iPlasticityCheck=list.GetParameter("plasticity_check");
    //  iplasticity = list.GetParameter("plasticity");
    iIterationMax = list.GetParameter("max_local_iterations");

    dAbsTol = list.GetParameter("local_tol_absolute");
    dRelTol = list.GetParameter("local_tol_relative");
    dYieldTrialTol = list.GetParameter("local_yield_tr_tol");

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
    fMaterial_Params[kc0] = list.GetParameter("c0");
    fMaterial_Params[kHc] = list.GetParameter("Hc");
    fMaterial_Params[kc0_chi] = list.GetParameter("c0_chi");
    fMaterial_Params[kHc_chi] = list.GetParameter("Hc_chi");
/*    fMaterial_Params[kGc0_chi1] = list.GetParameter("Gc0_chi1");
    fMaterial_Params[kGc0_chi2] = list.GetParameter("Gc0_chi2");
    fMaterial_Params[kGc0_chi3] = list.GetParameter("Gc0_chi3");
    fMaterial_Params[kHGc_chi] = list.GetParameter("HGc_chi");
*/
    //fMaterial_Params[kZ0c] =0.0;
    fMaterial_Params[kFphi] = list.GetParameter("Fphi");
    fMaterial_Params[kDpsi] = list.GetParameter("Dpsi");
    fMaterial_Params[kFphi_chi] = list.GetParameter("Fphi_chi");
    fMaterial_Params[kDpsi_chi] = list.GetParameter("Dpsi_chi");
//    fMaterial_Params[kFGphi_chi] = list.GetParameter("FGphi_chi");
//    fMaterial_Params[kDGpsi_chi] = list.GetParameter("DGpsi_chi");
    fMaterial_Params[kg] = list.GetParameter("g");
    fMaterial_Params[kg1] = list.GetParameter("g1");
    fMaterial_Params[kg2] = list.GetParameter("g2");
    fMaterial_Params[kg3] = list.GetParameter("g3");
    fMaterial_Params[kRho_0] = list.GetParameter("rho_0");

//    fIntegration_Params[kBeta] = list.GetParameter("beta");
//    fIntegration_Params[kGamma] = list.GetParameter("gamma");

    Echo_Input_Data();

    knum_d_state=24;//with extra strain measures

    knum_i_state = 0; // int's needed per ip, state variables

    //need to change these for non-symmetric stress, and higher order couple stress output
    knumstrain = 9; // number of strain outputs
    knumstress = 9; // number of stress outputs + higher order = ??

/*    knumstrain = 6; // number of strain outputs
    knumstress = 6; // number of stress outputs + higher order = ??
*/

  //  knumdispl=3;

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
    n_sd_x_n_sd = n_sd*n_sd;
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


    fShapeMicro_row_matrix.Dimension (1,n_en_micro);

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

    fDefGradInv_Grad_grad.Dimension (n_sd_x_n_sd, n_sd_x_n_sd);

    fRight_Cauchy_Green_tensor.Dimension (n_sd,n_sd);
    fRight_Cauchy_Green_tensor_Inverse.Dimension (n_sd,n_sd);
    fRight_Elastic_Cauchy_Green_tensor.Dimension (n_sd,n_sd);
    fRight_Elastic_Cauchy_Green_tensor_tr.Dimension (n_sd,n_sd);

    fMicroRight_Cauchy_Green_tensor.Dimension (n_sd,n_sd);
    fMicroRight_Cauchy_Green_tensor_tr.Dimension (n_sd,n_sd);



    fEulerian_strain_tensor_current_IP.Dimension (n_sd,n_sd);

    fTemp_nine_values.Dimension(9);
    fCauchy_stress_IPs.Dimension (fNumIP_displ,knumstress);



    fLeft_Cauchy_Green_tensor.Dimension (n_sd,n_sd);
    fIdentity_matrix.Dimension (n_sd,n_sd);
    fTemp_matrix_nsd_x_nsd.Dimension (n_sd,n_sd);
    fTemp_matrix_nsd_x_nsd2.Dimension (n_sd,n_sd);
    fTemp_matrix_nsd_x_nsd3.Dimension (n_sd,n_sd);

    fIota_temp_matrix.Dimension (n_en_displ_x_n_sd,n_sd_x_n_sd);

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


///////////////////////Used For the current Configuration Plasticity/////////////////////////
       fdGdCauchy_Stress.Dimension(n_sd,n_sd);
       fdGdCauchy_Stress_Transpose.Dimension(n_sd,n_sd);
       fdGdCauchy_Stress_IPs.Dimension(fNumIP_displ,n_sd_x_n_sd);
       fdGdCauchy_Stress_Elements_IPs.Dimension(NumElements(),fNumIP_displ*n_sd_x_n_sd);
       fdGdCauchy_Stress_Elements_IPs = 0.0;
       fdGdCauchy_Stress = 0.0;
       fdGdCauchy_Stress_Transpose = 0.0;


       fdGdCauchy_Stress_n_Transpose.Dimension(n_sd,n_sd);
       fdGdCauchy_Stress_n.Dimension(n_sd,n_sd);
       fdGdCauchy_Stress_n_IPs.Dimension(fNumIP_displ,n_sd_x_n_sd);
       fdGdCauchy_Stress_n_Elements_IPs.Dimension(NumElements(),fNumIP_displ*n_sd_x_n_sd);
       fdGdCauchy_Stress_n_Transpose = 0.0;
       fdGdCauchy_Stress_n_Elements_IPs = 0.0;
       fdGdCauchy_Stress_n = 0.0;

       fdGdCauchy_Stress_tr.Dimension(n_sd,n_sd);
       fdGdCauchy_Stress_tr_Transpose.Dimension(n_sd,n_sd);
       fdGdCauchy_Stress_tr_Transpose = 0.0;
       fdGdCauchy_Stress_tr = 0.0;


       fdFYdCauchy_Stress.Dimension(n_sd,n_sd);
       fdFYdCauchy_Stress_IPs.Dimension(fNumIP_displ,n_sd_x_n_sd);
       fdFYdCauchy_Stress_Elements_IPs.Dimension(NumElements(),fNumIP_displ*n_sd_x_n_sd);
       fdFYdCauchy_Stress_Elements_IPs = 0.0;
       fdFYdCauchy_Stress = 0.0;


       fdFYdCauchy_Stress_n_Transpose.Dimension(n_sd,n_sd);
       fdFYdCauchy_Stress_n.Dimension(n_sd,n_sd);
       fdFYdCauchy_Stress_n_IPs.Dimension(fNumIP_displ,n_sd_x_n_sd);
       fdFYdCauchy_Stress_n_Elements_IPs.Dimension(NumElements(),fNumIP_displ*n_sd_x_n_sd);
       fdFYdCauchy_Stress_n_Transpose = 0.0;
       fdFYdCauchy_Stress_n_Elements_IPs = 0.0;
       fdFYdCauchy_Stress_n = 0.0;


       fCauchy_stress_tensor_current_IP_Predictor.Dimension(n_sd,n_sd);
       fCauchy_stress_tensor_current_IP_Corrector.Dimension(n_sd,n_sd);
       fCauchy_stress_tensor_current_IP.Dimension(n_sd,n_sd);
       fCauchy_stress_IPs.Dimension(fNumIP_displ,n_sd_x_n_sd);
       fCauchy_stress_Elements_IPs.Dimension(NumElements(),fNumIP_displ*n_sd_x_n_sd);
       fdCauchy_stressdDelgamma.Dimension(n_sd,n_sd);
       fdDev_Cauchy_stressdDelgamma.Dimension(n_sd,n_sd);
       fCauchy_stress_Elements_IPs = 0.0;
       fCauchy_stress_tensor_current_IP = 0.0;
       fdCauchy_stressdDelgamma = 0.0;
       fdDev_Cauchy_stressdDelgamma = 0.0;
       fCauchy_stress_tensor_current_IP_Corrector = 0.0;
       fCauchy_stress_tensor_current_IP_Predictor = 0.0;

       fCauchy_stress_tensor_current_IP_tr.Dimension(n_sd,n_sd);
       fCauchy_stress_tensor_current_IP_tr = 0.0;


       dev_Cauchy_stress_tr.Dimension(n_sd,n_sd);
       fdGdS_tr.Dimension(n_sd,n_sd);


       dev_Cauchy_stress_n.Dimension(n_sd,n_sd);
       Predictor_stress_terms.Dimension(n_sd,n_sd);
       Corrector_stress_terms.Dimension(n_sd,n_sd);
       mean_Cauchy_stress_tr = 0.0;
       dev_Cauchy_stress_tr = 0.0;
       Cauchy_Stress_Norm_tr = 0.0;
       fdGdS_tr = 0.0;
       mean_Cauchy_stress_n = 0.0;
       dev_Cauchy_stress_n = 0.0;
       Predictor_stress_terms = 0.0;
       Corrector_stress_terms = 0.0;
       Predictor_mean_stress_terms = 0.0;
       Corrector_mean_stress_terms = 0.0;
       fTemp_matrix_one_x_one = 0.0;

       fdev_Cauchy_stress_tensor_current_n_IP.Dimension(n_sd,n_sd);
       fCauchy_stress_tensor_current_n_IP.Dimension(n_sd,n_sd);
       fdev_Cauchy_stress_tensor_current_IP.Dimension(n_sd,n_sd);
       fCauchy_stress_n_IPs.Dimension(fNumIP_displ,n_sd_x_n_sd);
       fCauchy_stress_Elements_n_IPs.Dimension(NumElements(),fNumIP_displ*n_sd_x_n_sd);
       fCauchy_stress_Elements_n_IPs = 0.0;
       fCauchy_stress_tensor_current_n_IP = 0.0;
       fdev_Cauchy_stress_tensor_current_IP = 0.0;

       fDeformation_Gradient.Dimension(n_sd,n_sd);
       fDeformation_Gradient_IPs.Dimension(fNumIP_displ,n_sd_x_n_sd);
       fDeformation_Gradient_Elements_IPs.Dimension(NumElements(),fNumIP_displ*n_sd_x_n_sd);
       fDel_Deformation_Gradient.Dimension(n_sd,n_sd);
       fDel_Deformation_Gradient = 0.0;


       fDeformation_Gradient_n.Dimension(n_sd,n_sd);
       fDeformation_Gradient_n_IPs.Dimension(fNumIP_displ,n_sd_x_n_sd);
       fDeformation_Gradient_n_Elements_IPs.Dimension(NumElements(),fNumIP_displ*n_sd_x_n_sd);
       fDeformation_Gradient_n_Elements_IPs = 0.0;

       fVelocity_Gradient_current_IP.Dimension(n_sd,n_sd);
       fSymmetric_part_Velocity_Gradient_current_IP.Dimension(n_sd,n_sd);
       fVelocity_Gradient_current_IP_transpose.Dimension(n_sd,n_sd);
       fElastic_Velocity_Gradient_current_IP.Dimension(n_sd,n_sd);
       fElastic_Velocity_Gradient_current_IP_transpose.Dimension(n_sd,n_sd);
       fVelocity_Gradient_current_IP = 0.0;
       fSymmetric_part_Velocity_Gradient_current_IP = 0.0;
       fElastic_Velocity_Gradient_current_IP = 0.0;
       fElastic_Velocity_Gradient_current_IP_transpose = 0.0;



       fTemp_Runge_Kutta_K1_nsd_x_nsd.Dimension(n_sd,n_sd);
       fTemp_Runge_Kutta_K2_nsd_x_nsd.Dimension(n_sd,n_sd);
       fTemp_Runge_Kutta_K3_nsd_x_nsd.Dimension(n_sd,n_sd);
       fTemp_Runge_Kutta_K4_nsd_x_nsd.Dimension(n_sd,n_sd);

       fTemp_Runge_Kutta_K1_nsd_x_nsd = 0.0;
       fTemp_Runge_Kutta_K2_nsd_x_nsd = 0.0;
       fTemp_Runge_Kutta_K3_nsd_x_nsd = 0.0;
       fTemp_Runge_Kutta_K4_nsd_x_nsd = 0.0;

       fCauchy_stress_tensor_current_IP_from_piola_stress.Dimension(n_sd,n_sd);
       fCauchy_stress_tensor_current_IP_from_piola_stress = 0.0;

       fLeft_Cauchy_Green_tensor_current_IP.Dimension(n_sd,n_sd);
       fLeft_Cauchy_Green_tensor_current_IP_transpose.Dimension(n_sd,n_sd);
       fLeft_Cauchy_Green_tensor_current_IP_Inverse.Dimension(n_sd,n_sd);


       fLeft_Cauchy_Green_tensor_current_IP = 0.0;
       fLeft_Cauchy_Green_tensor_current_IP_transpose = 0.0;
       fLeft_Cauchy_Green_tensor_current_IP_Trace = 0.0;
       fLeft_Cauchy_Green_tensor_current_IP_Inverse = 0.0;

       fRight_Elastic_Cauchy_Green_tensor_Inverse.Dimension(n_sd,n_sd);
       fRight_Elastic_Cauchy_Green_tensor_Inverse = 0.0;
///////////////////////////////////////////////////////////////////////////////////////////////
       ////////////////////// Global Consistant tangent Matrices///////////////////////////////
       fV1p.Dimension(n_sd_x_n_sd);
       Vintp_1.Dimension(n_en_displ_x_n_sd);
       Vintp_1_temp.Dimension(n_en_displ_x_n_sd);

       fV1e.Dimension(n_sd_x_n_sd);
       Vinte_1.Dimension(n_en_displ_x_n_sd);
       Vinte_1_temp.Dimension(n_en_displ_x_n_sd);

       /////////////////////////////// Elastic matrices////////////////////////////////////////
       IJe_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       IJe_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       IJe_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       IJe_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       IJe_5.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       IJe_6.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       IJe_7.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       IJe_8.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);

       ////////////////////////////////////////////////////////////////////////////////////////
       fKu_IJe_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_IJe_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_IJe_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_IJe_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_IJe_5.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_IJe_6.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_IJe_7.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_IJe_8.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);

       //////////////////////////// Plastic matrices/////////////////////////////////////////
       I1p_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I1p_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I1p_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I1p_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I1p_5.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I1p_6.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);


       I2p_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I2p_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I2p_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I2p_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I2p_5.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I2p_6.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);


       I3p_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I3p_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I3p_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I3p_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I3p_5.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I3p_6.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);


       I4p_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I4p_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I4p_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I4p_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I4p_5.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I4p_6.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);


       I5p_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I5p_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I5p_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I5p_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I5p_5.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I5p_6.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);


       I6p_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I6p_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I6p_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I6p_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I6p_5.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I6p_6.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);

       I7p_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I7p_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I7p_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I7p_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I7p_5.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I7p_6.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I7p_7.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I7p_8.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I7p_9.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I7p_10.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I7p_11.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I7p_12.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I7p_13.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);


       I8p_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I8p_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I8p_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I8p_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I8p_5.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I8p_6.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I8p_7.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I8p_8.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I8p_9.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I8p_10.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I8p_11.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I8p_12.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I8p_13.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);


       I9p_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I9p_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I9p_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I9p_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I9p_5.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I9p_6.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I9p_7.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I9p_8.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I9p_9.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I9p_10.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I9p_11.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I9p_12.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I9p_13.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);


       I10p_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I10p_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I10p_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I10p_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I10p_5.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I10p_6.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I10p_7.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I10p_8.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I10p_9.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I10p_10.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I10p_11.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I10p_12.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I10p_13.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);


       ////////////////////////////Test///////////////////////////
       /////////////////comparison////////////////////////////////
       I_temp_11p_test_2_1_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_1_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_1_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_1_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);

       I_temp_11p_test_2_2_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_2_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_2_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_2_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);

       I_temp_11p_test_2_3_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_3_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);


       I_temp_11p_test_2_4_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_4_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_4_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_4_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);

       I_temp_11p_test_2_5_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_5_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);


       I_temp_11p_test_2_6_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_6_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_6_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_6_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);

       I_temp_11p_test_2_7_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_7_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);


       I_temp_11p_test_2_8_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_8_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_8_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_8_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);

       I_temp_11p_test_2_9_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_9_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_9_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_9_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);

       I_temp_11p_test_2_10_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_10_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_10_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_10_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);

       I_temp_11p_test_2_11_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_11_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_11_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_11_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);

       I_temp_11p_test_2_12_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_12_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_12_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_12_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);

       I_temp_11p_test_2_13_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_13_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_13_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2_13_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);


       fKu_I_temp_11p_test_2_1_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_1_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_1_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_1_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);

       fKu_I_temp_11p_test_2_2_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_2_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_2_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_2_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);

       fKu_I_temp_11p_test_2_3_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_3_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);


       fKu_I_temp_11p_test_2_4_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_4_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_4_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_4_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);

       fKu_I_temp_11p_test_2_5_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_5_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);


       fKu_I_temp_11p_test_2_6_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_6_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_6_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_6_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);

       fKu_I_temp_11p_test_2_7_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_7_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);


       fKu_I_temp_11p_test_2_8_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_8_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_8_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_8_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);

       fKu_I_temp_11p_test_2_9_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_9_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_9_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_9_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);

       fKu_I_temp_11p_test_2_10_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_10_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_10_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_10_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);

       fKu_I_temp_11p_test_2_11_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_11_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_11_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_11_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);

       fKu_I_temp_11p_test_2_12_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_12_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_12_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_12_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);

       fKu_I_temp_11p_test_2_13_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_13_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_13_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_test_2_13_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);

       ///////////////////////////////////////////////////////////
       I_temp_11p_test_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_test_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I12p_test_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I12p_test_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I12p_test_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I_temp_11p_20.Dimension(n_sd,n_sd);
       I_temp_11p_21.Dimension(n_sd,n_sd);
       I_temp_11p_22.Dimension(n_sd,n_sd);
       I_temp_11p_23.Dimension(n_sd,n_sd);
       I_temp_11p_24.Dimension(n_sd,n_sd);
       I_temp_11p_4_Transpose.Dimension(n_sd,n_sd);

       II_temp_11p_1_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_1_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_1_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_1_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_1_5.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_1_6.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);

       II_temp_11p_2_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_2_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_2_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_2_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_2_5.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_2_6.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);

       II_temp_11p_3_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_3_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_3_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_3_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_3_5.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_3_6.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);

       II_temp_11p_4_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_4_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_4_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_4_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_4_5.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_4_6.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);

       II_temp_11p_5_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_5_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_5_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_5_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_5_5.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_5_6.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);

       II_temp_11p_6_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_6_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_6_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_6_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_6_5.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_6_6.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);

       II_temp_11p_7_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_7_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_7_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_7_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_7_5.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_7_6.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);

       II_temp_11p_8_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_8_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_8_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_8_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_8_5.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_8_6.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);

       II_temp_11p_9_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_9_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_9_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_9_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_9_5.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_9_6.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);

       II_temp_11p_10_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_10_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_10_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_10_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_10_5.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_10_6.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);

       II_temp_11p_11_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_11_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_11_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_11_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_11_5.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_11_6.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);

       II_temp_11p_12_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_12_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_12_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_12_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_12_5.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_12_6.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);

       II_temp_11p_13_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_13_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_13_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_13_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_13_5.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       II_temp_11p_13_6.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);



       fKu_I_temp_11p_1_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_1_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_1_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_1_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_1_5.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_1_6.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);

       fKu_I_temp_11p_2_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_2_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_2_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_2_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_2_5.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_2_6.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);

       fKu_I_temp_11p_3_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_3_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_3_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_3_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_3_5.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_3_6.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);

       fKu_I_temp_11p_4_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_4_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_4_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_4_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_4_5.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_4_6.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);

       fKu_I_temp_11p_5_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_5_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_5_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_5_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_5_5.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_5_6.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);

       fKu_I_temp_11p_6_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_6_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_6_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_6_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_6_5.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_6_6.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);

       fKu_I_temp_11p_7_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_7_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_7_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_7_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_7_5.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_7_6.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);

       fKu_I_temp_11p_8_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_8_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_8_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_8_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_8_5.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_8_6.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);

       fKu_I_temp_11p_9_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_9_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_9_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_9_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_9_5.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_9_6.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);

       fKu_I_temp_11p_10_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_10_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_10_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_10_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_10_5.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_10_6.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);

       fKu_I_temp_11p_11_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_11_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_11_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_11_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_11_5.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_11_6.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);

       fKu_I_temp_11p_12_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_12_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_12_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_12_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_12_5.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_12_6.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);

       fKu_I_temp_11p_13_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_13_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_13_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_13_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_13_5.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I_temp_11p_13_6.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);

       fKdd_previous.Dimension( n_en_displ_x_n_sd, n_en_displ_x_n_sd );
       fKdd_full_implemet.Dimension( n_en_displ_x_n_sd, n_en_displ_x_n_sd );
       difference.Dimension( n_en_displ_x_n_sd, n_en_displ_x_n_sd );

       ///////////////////////////////////////////////////////////


       I_temp_11p_1.Dimension(n_sd,n_sd);
       I_temp_11p_2.Dimension(n_sd,n_sd);
       I_temp_11p_3.Dimension(n_sd,n_sd);
       I_temp_11p_4.Dimension(n_sd,n_sd);
       I_temp_11p_5.Dimension(n_sd,n_sd);
       I_temp_11p_6.Dimension(n_sd,n_sd);
       I_temp_11p_7.Dimension(n_sd,n_sd);
       I_temp_11p_8.Dimension(n_sd,n_sd);
       I_temp_11p_9.Dimension(n_sd,n_sd);
       I_temp_11p_10.Dimension(n_sd,n_sd);
       I_temp_11p_11.Dimension(n_sd,n_sd);
       I_temp_11p_12.Dimension(n_sd,n_sd);
       I_temp_11p_13.Dimension(n_sd,n_sd);


       I11p_1.Dimension(n_sd,n_sd);
       I11p_2.Dimension(n_sd,n_sd);

       I12p_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I12p_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I12p_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I12p_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I12p_5.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I12p_6.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);

       I13p_1.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I13p_2.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I13p_3.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I13p_4.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I13p_5.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);
       I13p_6.Dimension(n_sd_x_n_sd,n_sd_x_n_sd);

       fKu_I1p_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I1p_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I1p_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I1p_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I1p_5.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I1p_6.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);




       fKu_I2p_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I2p_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I2p_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I2p_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I2p_5.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I2p_6.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);




       fKu_I3p_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I3p_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I3p_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I3p_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I3p_5.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I3p_6.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);




       fKu_I4p_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I4p_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I4p_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I4p_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I4p_5.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I4p_6.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);




       fKu_I5p_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I5p_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I5p_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I5p_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I5p_5.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I5p_6.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);




       fKu_I6p_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I6p_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I6p_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I6p_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I6p_5.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I6p_6.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);



       fKu_I7p_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I7p_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I7p_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I7p_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I7p_5.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I7p_6.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I7p_7.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I7p_8.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I7p_9.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I7p_10.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I7p_11.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I7p_12.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I7p_13.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);


       fKu_I8p_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I8p_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I8p_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I8p_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I8p_5.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I8p_6.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I8p_7.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I8p_8.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I8p_9.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I8p_10.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I8p_11.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I8p_12.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I8p_13.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);


       fKu_I9p_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I9p_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I9p_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I9p_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I9p_5.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I9p_6.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I9p_7.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I9p_8.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I9p_9.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I9p_10.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I9p_11.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I9p_12.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I9p_13.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);


       fKu_I10p_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I10p_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I10p_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I10p_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I10p_5.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I10p_6.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I10p_7.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I10p_8.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I10p_9.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I10p_10.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I10p_11.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I10p_12.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I10p_13.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);

       fKu_I12p_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I12p_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I12p_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I12p_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I12p_5.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I12p_6.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);

       fKu_I13p_1.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I13p_2.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I13p_3.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I13p_4.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I13p_5.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);
       fKu_I13p_6.Dimension(n_en_displ_x_n_sd,n_en_displ_x_n_sd);


///////////////////////////////////////////////////////////////////////////////////////////////


    /////////////////////////////////////////////////////////////
    /////////////FINITE STRAIN ELASTICITY MATRICES///////////////
    /////////////////////////////////////////////////////////////
    KirchhoffST.Dimension(n_sd,n_sd);
    SPK.Dimension(n_sd,n_sd);
    Temp_SPK.Dimension(n_sd,n_sd);
   // FSF.Dimension(n_sd,n_sd);


    ChiM.Dimension(n_sd,n_sd);
    ChiM_Inverse.Dimension(n_sd,n_sd);


    gammastn.Dimension(n_sd,n_sd,n_sd);
    trgammastn.Dimension(n_sd);
    devgammastn.Dimension(n_sd,n_sd,n_sd);

    fEulerian_strain_tensor_current_IP.Dimension (n_sd,n_sd);
    fCauchy_stress_tensor_current_IP.Dimension (n_sd,n_sd);
    //
    fDisplacements_current_IPs.Dimension(n_sd);
    fEulerian_strain_IPs.Dimension (fNumIP_displ,knumstrain);
    fCauchy_stress_IPs.Dimension (fNumIP_displ,knumstress);

   //
 //   fDisplacement_IPs.Dimension(fNumIP_displ,knumdispl);
    fTemp_nine_values.Dimension(9);
    fTemp_six_values.Dimension(6);
    fEulerian_strain_Elements_IPs.Dimension (NumElements(),fNumIP_displ*knumstrain);
    fCauchy_stress_Elements_IPs.Dimension (NumElements(),fNumIP_displ*knumstress);


    if(iConstitutiveModelType==2)
    {
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
   }


//    GAMMA.Dimension(n_sd,n_sd,n_sd);
    GRAD_CHIM.Dimension(n_sd,n_sd,n_sd);
    fState_variables_IPs.Dimension (fNumIP_displ,kNUM_FMATERIAL_STATE_TERMS);
    fState_variables_IPs=0.0;
    fState_variables_Elements_IPs.Dimension (NumElements(),fNumIP_displ*kNUM_FMATERIAL_STATE_TERMS);
    fState_variables_Elements_IPs=0.0;
    fState_variables_n_IPs.Dimension (fNumIP_displ,kNUM_FMATERIAL_STATE_TERMS);
    fState_variables_n_IPs=0.0;
    fState_variables_n_Elements_IPs.Dimension (NumElements(),fNumIP_displ*kNUM_FMATERIAL_STATE_TERMS);
    fState_variables_n_Elements_IPs=0.0;


    ////////////////stress measures/////////////////

    devsigma.Dimension(n_sd,n_sd);
    devRelsts.Dimension(n_sd,n_sd);
    s_sigma_temp.Dimension(n_sd,n_sd);
    fmklm.Dimension(n_sd,n_sd,n_sd);
    devmklm.Dimension(n_sd,n_sd,n_sd);
    trvecmklm.Dimension(n_sd);


    // if(iConstitutiveModelType==3)
     // {
      fState_variables_IPs.Dimension (fNumIP_displ,kNUM_FMATERIAL_STATE_TERMS);
      fState_variables_IPs=0.0;
      fState_variables_Elements_IPs.Dimension (NumElements(),fNumIP_displ*kNUM_FMATERIAL_STATE_TERMS);
      fState_variables_Elements_IPs=0.0;
      fState_variables_n_IPs.Dimension (fNumIP_displ,kNUM_FMATERIAL_STATE_TERMS);
      fState_variables_n_IPs=0.0;
      fState_variables_n_Elements_IPs.Dimension (NumElements(),fNumIP_displ*kNUM_FMATERIAL_STATE_TERMS);
      fState_variables_n_Elements_IPs=0.0;
     //}

     fIdentity_matrix=0.0;
     fIdentity_matrix(0,0)=1.0;
     fIdentity_matrix(1,1)=1.0;
     fIdentity_matrix(2,2)=1.0;
    // fCe_n=fIdentity_matrix;
      Beta=-1.0;
      Aphi=2*sqrt(6)*cos(fMaterial_Params[kFphi])/(3+Beta*sin(fMaterial_Params[kFphi]));
      Bphi=2*sqrt(6)*sin(fMaterial_Params[kFphi])/(3+Beta*sin(fMaterial_Params[kFphi]));
      Apsi=2*sqrt(6)*cos(fMaterial_Params[kDpsi])/(3+Beta*sin(fMaterial_Params[kDpsi]));
      Bpsi=2*sqrt(6)*sin(fMaterial_Params[kDpsi])/(3+Beta*sin(fMaterial_Params[kDpsi]));



   Top();



     while (NextElement())
     {

         int e,l;
         e = CurrElementNumber();

         for (l=0; l < fNumIP_displ; l++)
         {

             fState_variables_n_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kc)
             =fMaterial_Params[kc0];
             fState_variables_n_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kc_chi)
             =fMaterial_Params[kc0_chi];
             /*fState_variables_n_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kGc_chi1)
             =fMaterial_Params[kGc0_chi1];// Micro Gradient Pl.
             fState_variables_n_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kGc_chi2)
             =fMaterial_Params[kGc0_chi2];// Micro Gradient Pl.
             fState_variables_n_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kGc_chi3)
             =fMaterial_Params[kGc0_chi3];// Micro Gradient Pl.*/

             //fState_variables_n_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+khkappa)=0.0;
             fState_variables_n_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+khc)=Apsi;
             fState_variables_n_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+khc_chi)=Apsi_chi;
             //fState_variables_n_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+khGc_chi)=AGpsi_chi;// Micro Gradient Pl.

             fState_variables_n_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kDelgamma)=0.0;
             fState_variables_n_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kDelgammachi)=0.0;
             //fState_variables_n_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kDelgammaGchi)=0.0;// Micro Gradient Pl.
             fState_variables_n_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+ktrCauchy_Stress)=0.0;
             fState_variables_n_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kNorm_Dev_Cauchy_Stress)=0.0;
             fState_variables_n_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+ktrRel)=0.0;
             fState_variables_n_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kRel_inv)=0.0;
             fState_variables_n_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kinvtrM)=0.0;
             fState_variables_n_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kinvdevM)=0.0;
             fState_variables_n_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kinvPhi)=0.0;
             fState_variables_n_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kinvGPhi)=0.0;
             fState_variables_n_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+ktreps)=0.0;
             fState_variables_n_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kdeveps)=0.0;
             fState_variables_n_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kinvtrgammastn)=0.0;
             fState_variables_n_Elements_IPs(e,l*kNUM_FMATERIAL_STATE_TERMS+kinvdevgammastn)=0.0;


             fdGdCauchy_Stress_n_IPs.SetRow(l,fdGdCauchy_Stress_n);//
             fCauchy_stress_n_IPs.SetRow(l,fCauchy_stress_tensor_current_n_IP);//
             fDeformation_Gradient_n_IPs.SetRow(l,fIdentity_matrix);//
         }
         fdGdCauchy_Stress_n_Elements_IPs.SetRow(e,fdGdCauchy_Stress_n_IPs);
         fCauchy_stress_Elements_n_IPs.SetRow(e,fCauchy_stress_n_IPs);
         fDeformation_Gradient_n_Elements_IPs.SetRow(e,fDeformation_Gradient_n_IPs);
     }

     fdGdCauchy_Stress_Elements_IPs = fdGdCauchy_Stress_n_Elements_IPs;
     fCauchy_stress_Elements_IPs = fCauchy_stress_Elements_n_IPs;
     fDeformation_Gradient_Elements_IPs = fDeformation_Gradient_n_Elements_IPs;
     fState_variables_Elements_IPs=fState_variables_n_Elements_IPs;



    ///////////////////////////////////////////////////////////////////////////
    /////////////DIMENSIONALIZE MICROMORPHIC MATRICES FINISH HERE FOR 3D CASE//////////////
    ///////////////////////////////////////////////////////////////////////////




    u_dotdot_column_matrix.Dimension (n_en_displ_x_n_sd,1);
    u_dot_column_matrix.Dimension (n_en_displ_x_n_sd,1);
    u_dot_column_matrix_Transpose.Dimension (1, n_en_displ_x_n_sd);


    /* streams */
    ofstreamT& out = ElementSupport().Output();

    /* storage for integration point strain, stress, and ISVs*/
    fIPVariable.Dimension (n_el, fNumIP_displ*(knumstrain+knumstress+knum_d_state));
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
void FSMicromorphic3DCurrConfigT::DefineSubs(SubListT& sub_list) const
{
  //            cout<<"CHECK POINT-24"<<endl;
    /* inherited */
    ElementBaseT::DefineSubs(sub_list);

    /* element blocks */
    sub_list.AddSub("micromorphic_FS_3D_inc_curr_config_element_block");

    /* tractions */
    sub_list.AddSub("micromorphic_FS_3D_inc_curr_config_natural_bc", ParameterListT::Any);
}



/* return the description of the given inline subordinate parameter list */
void FSMicromorphic3DCurrConfigT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order,
                       SubListT& sub_lists) const
{
   //              cout<<"CHECK POINT-25"<<endl;
    ElementBaseT::DefineInlineSub(name, order, sub_lists);
}



/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FSMicromorphic3DCurrConfigT::NewSub(const StringT& name) const
{
   //           cout<<"CHECK POINT-26"<<endl;
    /* create non-const this */
    FSMicromorphic3DCurrConfigT* non_const_this = const_cast<FSMicromorphic3DCurrConfigT*>(this);

    if (name == "micromorphic_FS_3D_inc_curr_config_natural_bc") /* traction bc */
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
    else if (name == "micromorphic_FS_3D_inc_curr_config_element_block")
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
void FSMicromorphic3DCurrConfigT::SetTractionBC(void)
{
//NOTE: With the possibility of variable global node numbers and
//      and equations, we assume as little as possible here with
//      regard to the validity of the node/equation numbers, requiring
//      only that NodesX in the element cards has the correct global
//      node numbers.
      //        cout<<"CHECK POINT-27"<<endl;
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
void FSMicromorphic3DCurrConfigT::TakeNaturalBC(const ParameterListT& list)
{
  //             cout<<"CHECK POINT-28"<<endl;
    const char caller[] = "FSMicromorphic3DCurrConfigT::TakeTractionBC";

    int num_natural_bc = list.NumLists("micromorphic_FS_3D_inc_curr_config_natural_bc");
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
        const ParameterListT& natural_bc = list.GetList("micromorphic_FS_3D_inc_curr_config_natural_bc", i);

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
void FSMicromorphic3DCurrConfigT::ApplyTractionBC(void)
{
    //           cout<<"CHECK POINT-29"<<endl;
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

void FSMicromorphic3DCurrConfigT::Form_solid_shape_functions(const double* &shapes_displ_X)
{
    fShapeDispl = 0.0;
    for (int i=0; i<8; i++)
    {
        fShapeDispl(0,i*3) = shapes_displ_X[i];
        fShapeDispl(1,1+i*3) = shapes_displ_X[i];
        fShapeDispl(2,2+i*3) = shapes_displ_X[i];
    }
}

void FSMicromorphic3DCurrConfigT::Form_Gradient_of_solid_shape_functions(const dMatrixT &fShapeDisplGrad_temp)
{
    fShapeDisplGrad = 0.0;
    for(int i=0; i<8; i++)
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

void FSMicromorphic3DCurrConfigT::Form_micro_shape_functions(const double* &shapes_micro_X)
{
    fShapeMicro = 0.0;
    //hard coded for n_en_micro=8; can change
    for (int i=0; i<n_en_micro; i++)
        fShapeMicro[i] = shapes_micro_X[i];
}

void FSMicromorphic3DCurrConfigT::Form_deformation_gradient_tensor(void)
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

void FSMicromorphic3DCurrConfigT::Form_Grad_grad_transformation_matrix(void)
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




//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
////////////////////// MATRICES FOR MICROMORPHIC 3-D CASE/////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

void FSMicromorphic3DCurrConfigT:: Form_CCof_tensor()
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

void FSMicromorphic3DCurrConfigT::Form_micro_deformation_tensor_Chi()
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

void FSMicromorphic3DCurrConfigT:: Form_Chi_inv_matrix()
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
void FSMicromorphic3DCurrConfigT:: Form_GRAD_Chi_matrix()
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
void FSMicromorphic3DCurrConfigT::Form_KroneckerDelta_matrix()
{
    for(int i=0;i<=2;i++)
    { for(int j=0;j<=2;j++)
            {KrDelta[i][j]=0.0;}}
    KrDelta[0][0]=1.0;
    KrDelta[1][1]=1.0;
    KrDelta[2][2]=1.0;
}

void FSMicromorphic3DCurrConfigT::Form_Gamma_tensor3D()
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

void FSMicromorphic3DCurrConfigT::Form_G1_matrix()
{
    int row=0;
 //   int col;
    double scale;
    scale=0.0;
  //  row=0;
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



    for(int k=0;k<=2;k++)
    {
        for(int l=0;l<=2;l++)
        {
             G1[row]=Sigma(l,k);
            row++;
        }
    }

}

void FSMicromorphic3DCurrConfigT::Form_Finv_w_matrix()//checked correct
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



void FSMicromorphic3DCurrConfigT::Form_NCHI_matrix(const dMatrixT &fShapeMicro_row_matrix)
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


void FSMicromorphic3DCurrConfigT:: Form_GRAD_Nuw_matrix(const dMatrixT &fShapeDisplGrad_temp)
{
    int row,col;
    row=0;
    col=0;
    GRAD_Nuw=0.0;
    for(int j=0;j<3;j++)
    {
        col=j;
        for(int i=0;i<8;i++)
        {
            GRAD_Nuw(row,col)  =fShapeDisplGrad_temp(0,i);
            GRAD_Nuw(row+1,col)=fShapeDisplGrad_temp(1,i);
            GRAD_Nuw(row+2,col)=fShapeDisplGrad_temp(2,i);
            col=col+3;
        }
        row=row+3;
    }


}


void FSMicromorphic3DCurrConfigT::Form_Var_F_tensor()
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
                        Var_F(row,col)+=Finv[K][i]*Sigma(l,k);
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

void FSMicromorphic3DCurrConfigT::Form_Tsigma_1_matrix()
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





}

void FSMicromorphic3DCurrConfigT::Form_Tsigma_2_matrix()
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


}

void FSMicromorphic3DCurrConfigT::Form_Tsigma_3_matrix()
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

void FSMicromorphic3DCurrConfigT::Form_TFn_1_matrix()
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

void FSMicromorphic3DCurrConfigT::Form_TFn_2_matrix()
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

void FSMicromorphic3DCurrConfigT::Form_TFn_3_matrix()
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
    }

void FSMicromorphic3DCurrConfigT::Form_TChi_1_matrix()
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



}

void FSMicromorphic3DCurrConfigT::Form_TFn_4_matrix()
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

}

void FSMicromorphic3DCurrConfigT::Form_TChi_2_matrix()
{
   int row,col;
   row=0;
   col=0;
   TChi_2=0.0;
/*
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
*/

   for(int R=0;R<3;R++)
    {
        for(int m=0;m<3;m++)
        {
            row=0;
            for(int k=0;k<3;k++)
            {
                for(int l=0;l<3;l++)
                {
                    //summation
                    for(int K=0;K<3;K++)
                    {
            //TChi_2(row,col)+=ChiN[l][K]*ChiInv[K][m]*ChiInv[R,k];
                        TChi_2(row,col)+=ChiN[l][K]*ChiInv_m(K,m)*ChiInv_m(R,k);

                    }
                    row++;
                }
            }
         col++;
        }
    }


}

void FSMicromorphic3DCurrConfigT::Form_TFn_5_matrix()
{
    int row,col;
    row=0;
    col=0;
    TFn_5=0.0;
 /*   for(int l=0;l<3;l++)
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
*/

for(int i=0;i<3;i++)
{
    for(int m=0;m<3;m++)
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
                    TFn_5(row,col)+=Fn[k][K]*Finv[K][m]*fIdentity_matrix(l,i);
                }
                row++;
            }
        }
        col++;
    }
}

}

void FSMicromorphic3DCurrConfigT::Form_TChi_3_matrix()
{
    int row;
    int col;
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

}

void FSMicromorphic3DCurrConfigT::Form_TFn_6_matrix()

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


}

void FSMicromorphic3DCurrConfigT::Form_double_Finv_from_Deformation_tensor_inverse()
{
    for(int i=0;i<=2;i++)
        {
            for(int j=0;j<=2;j++)
                {
                Finv[i][j]=fDeformation_Gradient_Inverse(i,j);}}
}

void FSMicromorphic3DCurrConfigT::Form_SigCurr_matrix()
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
void FSMicromorphic3DCurrConfigT::Form_Finv_eta_matrix()
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

void FSMicromorphic3DCurrConfigT::Form_Gradient_of_micro_shape_eta_functions(const dMatrixT &fShapeMicroGrad_temp)
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


void FSMicromorphic3DCurrConfigT:: Form_Etagrad_matrix()
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

void FSMicromorphic3DCurrConfigT::Form_Mm_1_matrix()
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


void FSMicromorphic3DCurrConfigT::Form_Mm_2_matrix()
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

void FSMicromorphic3DCurrConfigT::Form_Mm_3_matrix()
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

void FSMicromorphic3DCurrConfigT:: Form_Mm_4_matrix()
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

void FSMicromorphic3DCurrConfigT:: Form_Mm_5_matrix()
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

void FSMicromorphic3DCurrConfigT:: Form_Mm_6_matrix()
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


void FSMicromorphic3DCurrConfigT:: Form_Mm_7_matrix()
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


void FSMicromorphic3DCurrConfigT:: Form_Mm_71_matrix()
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

void FSMicromorphic3DCurrConfigT:: Form_Mm_72_matrix()
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

void FSMicromorphic3DCurrConfigT:: Form_Mm_73_matrix()
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


void FSMicromorphic3DCurrConfigT:: Form_Mm_74_matrix()
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


void FSMicromorphic3DCurrConfigT:: Form_Mm_75_matrix()
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


void FSMicromorphic3DCurrConfigT:: Form_Mm_76_matrix()
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

void FSMicromorphic3DCurrConfigT:: Form_Mm_77_matrix()
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

void FSMicromorphic3DCurrConfigT:: Form_Mm_78_matrix()
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

void FSMicromorphic3DCurrConfigT:: Form_Mm_8_matrix()
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

void FSMicromorphic3DCurrConfigT:: Form_Mm_9_matrix()
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

void FSMicromorphic3DCurrConfigT:: Form_Mm_10_matrix()
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

void FSMicromorphic3DCurrConfigT:: Form_Mm_11_matrix()
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


void FSMicromorphic3DCurrConfigT:: Form_Mm_12_matrix()
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

void FSMicromorphic3DCurrConfigT:: Form_Mm_13_matrix()
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


void FSMicromorphic3DCurrConfigT:: Form_Mm_14_matrix()
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

void FSMicromorphic3DCurrConfigT:: Form_Ru_1_matrix()
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

void FSMicromorphic3DCurrConfigT:: Form_Ru_2_matrix()
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

void FSMicromorphic3DCurrConfigT:: Form_Ru_3_matrix()
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
                                        Ru_3(row, col) +=sn_sigman(m,i)*Fn[l][L]*Finv[L][k];
                                     }
                            row++;
                            }
                        }
                    col++;
                }
            }
}

void FSMicromorphic3DCurrConfigT:: Form_RChi_1_matrix()
{

    int col;
    int row;
    RChi_1=0.0;
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
                                        RChi_1(row,col) +=ChiN[l][T]*ChiInv[T][p]*ChiInv[K][m];
                                    }
                            row++;
                           }
                        }
                    col++;
                }
            }


}

void FSMicromorphic3DCurrConfigT::Form_Ru_4_matrix()
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
                                    Ru_4(row,col)+=Fn[m][K]*Finv[K][k];
                                }
                                row=row+3;
                    }
                    col++;
                }
            }

  /*  for(int p=0;p<3;p++)
    {
        for(int k=0;k<3;k++)
        {
            row=0;
            for(int m=0;m<3;m++)
            {
                for(int l=0;l<3;l++)
                {
                    //summation
                    for(int K=0;K<3;K++)
                    {
                        Ru_4(row, col) +=Fn[m][K]*Finv[K][k]*fIdentity_matrix(p,l);
                    }
                    row++;
                }
            }
        col++;
        }
    }

*/
}

void FSMicromorphic3DCurrConfigT:: Form_RChi_2_matrix()
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


void FSMicromorphic3DCurrConfigT:: Form_Ru_5_matrix()
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

void FSMicromorphic3DCurrConfigT:: Form_Ru_6_matrix()
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
void FSMicromorphic3DCurrConfigT:: Form_Ru_7_matrix()
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

void FSMicromorphic3DCurrConfigT:: Form_RChi_3_matrix()
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

void FSMicromorphic3DCurrConfigT:: Form_Ru_8_matrix()
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

void FSMicromorphic3DCurrConfigT:: Form_Ru_9_matrix()
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



void FSMicromorphic3DCurrConfigT:: Form_Rs_sigma_matrix()
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

void FSMicromorphic3DCurrConfigT:: Form_R_Capital_Lambda_Chi_matrix()
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

void FSMicromorphic3DCurrConfigT:: Form_CapitalLambda_matrix()
{
    for(int i=0;i<=2;i++)
    {
        for(int j=0;j<=2;j++)
        {
            CapitalLambda(i,j)=0.0;//10^8
        }
    }
}

void FSMicromorphic3DCurrConfigT::Form_H1_matrix()
{
    int row=0;
    double trdeltad=0.0;
    double dtgd[3][3][3];
    double grad_Nu[3][3][3];
    double Cgamma[3][3][3];
    double dgcir[3][3][3];
    H1=0.0;

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

void FSMicromorphic3DCurrConfigT::Form_H2_matrix()
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
    scale=(fMaterial_Params[kKappa]-fMaterial_Params[kSigma_const]);
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

void FSMicromorphic3DCurrConfigT:: Form_H3_matrix()
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

void FSMicromorphic3DCurrConfigT:: Mapping_double_and_Array(const int condition)
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

void FSMicromorphic3DCurrConfigT:: Form_deformation_tensors_arrays(const int condition) //
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
void FSMicromorphic3DCurrConfigT::Form_ChiM()
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


void FSMicromorphic3DCurrConfigT::Form_Second_Piola_Kirchhoff_SPK(const dMatrixT& LagStn, const dMatrixT& MicroStn)
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

    trLST= LagStn(0,0)+LagStn(1,1)+LagStn(2,2);
    trcE = MicroStn(0,0)+MicroStn(1,1)+MicroStn(2,2);

   Temp_SPK=fIdentity_matrix;
   scale=(fMaterial_Params[kLambda]+fMaterial_Params[kTau]);
   scale*=trLST;
   Temp_SPK*=scale;
   SPK+=Temp_SPK;

   scale=2*(fMaterial_Params[kMu]+fMaterial_Params[kSigma_const]);
   Temp_SPK=LagStn;
   Temp_SPK*=scale;
   SPK+=Temp_SPK;

   scale=fMaterial_Params[kEta]*trcE;
   Temp_SPK=fIdentity_matrix;
   Temp_SPK*=scale;
   SPK+=Temp_SPK;

   Temp_SPK=MicroStn;
   Temp_SPK*=fMaterial_Params[kKappa];
   SPK+=Temp_SPK;

   Temp_SPK.Transpose(MicroStn);
   Temp_SPK*=fMaterial_Params[kNu];
   SPK+=Temp_SPK;



}


void FSMicromorphic3DCurrConfigT:: Calculate_relative_strain_INV()
{
     eps=0.0;
     psi=0.0;
     deveps=0.0;
     treps=0.0;
     temp_inv=0.0;
     psi.MultATB(fDeformation_Gradient_Inverse,ChiM);
     eps  = psi;
     eps *= -1;
     eps += fIdentity_matrix;
     //trace of epsilon
     treps=eps(0,0)+eps(1,1)+eps(2,2);
     //||deveps||
     temp_inv=eps.ScalarProduct();
     deveps=sqrt(temp_inv);


}


void FSMicromorphic3DCurrConfigT:: Calculate_HOST_INV()
{
	gammastn=0.0;
	invtrgammastn=0.0;
        invdevgammastn=0.0;
        devgammastn=0.0;
        trgammastn=0.0;
	for(int k=0;k<3;k++)
	   {
	    for(int l=0;l<3;l++)
	       {
               for(int m=0;m<3;m++)
                  {
                  //summation
		  for(int K=0;K<3;K++)
		     {
 		     for(int a=0;a<3;a++)
		        {
                        for(int M=0;M<3;M++)
		           {
                           for(int R=0;R<3;R++)
   	       	              {
                                gammastn(k,l,m)+=fDeformation_Gradient_Inverse(K,k)*ChiM_Inverse(K,a)*GRAD_CHIM(a,R,M)*ChiM_Inverse(R,l)
						*fDeformation_Gradient_Inverse(M,m);
 			      }
 			    }
                         }
	 	     }
                   }
               }
           }

    temp_inv=0;
    for(int i =0;i<3;i++)
    {
        for(int a=0;a<3;a++)
        {
            trgammastn[i]+=gammastn(a,i,a);
        }
        temp_inv+=trgammastn[i]*trgammastn[i];
    }

    invtrgammastn=sqrt(temp_inv);


    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            for( int k=0;k<3;k++)
            {
                //for(a=0;a<3;a++)
                //{
                devgammastn(i,j,k)=gammastn(i,j,k)-(1/3)*fIdentity_matrix(i,j)*gammastn(0,k,0)
                                                  -(1/3)*fIdentity_matrix(i,j)*gammastn(1,k,1)
                                                  -(1/3)*fIdentity_matrix(i,j)*gammastn(2,k,2);
                //}
            }
        }
    }


    temp_inv=0;
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            for(int k=0;k<3;k++)
            {
                temp_inv+=devgammastn(i,j,k)*devgammastn(i,j,k);
            }
        }
    }

    invdevgammastn=sqrt(temp_inv);

}








////////////////////////////////////////////////////////////////
//////////////FINITE STRAIN MATRICES ENDS//////////////////
////////////////////////////////////////////////////////////////

/*void FSMicromorphic3DCurrConfigT:: Extract_six_values_from_symmetric_tensor(const dMatrixT &fTensor,dArrayT& fTemp_six_values)*/
void FSMicromorphic3DCurrConfigT:: Extract_six_values_from_symmetric_tensor(const dMatrixT &fTensor,dArrayT& fTemp_nine_values)
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

void FSMicromorphic3DCurrConfigT::Put_values_In_dArrayT_vector(const dArray2DT &f2DArrayT,const int& e,const int& IP,dArrayT& fArrayT)
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


void FSMicromorphic3DCurrConfigT::Put_values_In_Array(const dArray2DT &f2DArrayT,const int& e,const int& IP,dArrayT& fArrayT)
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


//////////////////////// Plasticity Matrices for the Global Consistent Tangent//////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FSMicromorphic3DCurrConfigT:: Form_fV1p()
{
    int row=0;
    fTemp_matrix_nsd_x_nsd=0.0;
    fV1p=0.0;
    fTemp_matrix_nsd_x_nsd.MultAB(fDeformation_Gradient_Inverse,fCauchy_stress_tensor_current_IP);


    for (int Kbar=0;Kbar<3;Kbar++)
    {
        for(int l=0;l<3;l++)
        {
            fV1p[row]=fTemp_matrix_nsd_x_nsd(Kbar,l);
            row++;
        }
    }


}


void FSMicromorphic3DCurrConfigT:: Form_fV1e()
{
    int row=0;
    fTemp_matrix_nsd_x_nsd=0.0;
    fV1e=0.0;
    fTemp_matrix_nsd_x_nsd.MultAB(fDeformation_Gradient_Inverse,fCauchy_stress_tensor_current_IP);


    for (int Kbar=0;Kbar<3;Kbar++)
    {
        for(int l=0;l<3;l++)
        {
            fV1e[row]=fTemp_matrix_nsd_x_nsd(Kbar,l);
            row++;
        }
    }


}



void FSMicromorphic3DCurrConfigT:: Form_IJe_1()
{
    int row=0;
    int col=0;
    IJe_1=0.0;
    for(int I=0;I<3;I++)
    {
        for(int l=0;l<3;l++)
        {
            row=0;
            for(int L=0;L<3;L++)
            {
                for(int k=0;k<3;k++)
                {
                    //summation over the same term starts
                    for(int j=0;j<3;j++)
                    {
                        IJe_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_IP(j,k)*fDeformation_Gradient_Inverse(I,l);
                    }
                    row++;
                }

            }
            col++;
        }
    }
}

void FSMicromorphic3DCurrConfigT:: Form_IJe_2()
{
    int row=0;
    int col=0;
    IJe_2=0.0;
    for(int I=0;I<3;I++)
    {
        for(int l=0;l<3;l++)
        {
            row=0;
            for(int L=0;L<3;L++)
            {
                for(int k=0;k<3;k++)
                {
                    //summation over the same term starts
                    for(int j=0;j<3;j++)
                    {
                        IJe_2(row,col)+=fDeformation_Gradient_Inverse(L,l)*fDeformation_Gradient_Inverse(I,j)*fCauchy_stress_tensor_current_IP(j,k);
                    }
                    row++;
                }

            }
            col++;
        }
    }
}



void FSMicromorphic3DCurrConfigT:: Form_IJe_3()
{
    int row=0;
    int col=0;
    IJe_3=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							//summation over the same term starts
							for(int i=0;i<3;i++)
							{
								IJe_3(row,col)+=fDeformation_Gradient_Inverse(L,i)*fDeformation_Gradient_n(j,K)*fDeformation_Gradient_Inverse(K,l)
								*fDeformation_Gradient_Inverse(I,j)*fCauchy_stress_tensor_current_n_IP(i,k);
							}

						}

					}
				 row++;
				}
			}
		col++;
		}

	}
}


void FSMicromorphic3DCurrConfigT:: Form_IJe_4()
{

    int row=0;
    int col=0;
    IJe_4=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
						{
							//summation over the same term starts
							for(int i=0;i<3;i++)
							{
								IJe_4(row,col)+=fDeformation_Gradient_Inverse(L,i)*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(K,l)
								*fDeformation_Gradient_Inverse(I,p)*fCauchy_stress_tensor_current_n_IP(p,k);
							}

						}

					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_IJe_5()
{

    int row=0;
    int col=0;
    IJe_5=0.0;
    for (int J=0;J<3;J++)
	{
		for(int i=0;i<3;i++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int I=0;I<3;I++)
					{
						for(int p=0;p<3;p++)
						{
							//summation over the same term starts
							for(int j=0;j<3;j++)
							{
								IJe_5(row,col)+=fDeformation_Gradient_Inverse(L,j)*fDeformation_Gradient_n(k,I)*fDeformation_Gradient_Inverse(I,i)
								*fDeformation_Gradient_Inverse(J,p)*fCauchy_stress_tensor_current_n_IP(j,p);
							}

						}

					}
				 row++;
				}
			}
		col++;
		}
	}
}




void FSMicromorphic3DCurrConfigT:: Form_IJe_6()
{

    int row=0;
    int col=0;
    IJe_6=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int i=0;i<3;i++)
						{
							//summation over the same term starts
							for(int j=0;j<3;j++)
							{
								IJe_6(row,col)+=fDeformation_Gradient_Inverse(L,j)*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(K,l)
								*fDeformation_Gradient_Inverse(I,i)*fIdentity_matrix(j,k);
							}

						}

					}
				 row++;
				}
			}
		col++;
		}
	}
}


void FSMicromorphic3DCurrConfigT:: Form_IJe_7()
{

    int row=0;
    int col=0;
    IJe_7=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
							for(int j=0;j<3;j++)
							{
								IJe_7(row,col)+=fDeformation_Gradient_Inverse(L,j)*fDeformation_Gradient_n(j,K)*fDeformation_Gradient_Inverse(K,l)
								*fDeformation_Gradient_Inverse(I,k);
							}

					}
				 row++;
				}
			}
		col++;
		}
	}
}


void FSMicromorphic3DCurrConfigT:: Form_IJe_8()
{

    int row=0;
    int col=0;
    IJe_8=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
							for(int j=0;j<3;j++)
							{
								IJe_8(row,col)+=fDeformation_Gradient_Inverse(L,j)*fDeformation_Gradient_n(k,K)*fDeformation_Gradient_Inverse(K,l)
								*fDeformation_Gradient_Inverse(I,j);
							}

					}
				 row++;
				}
			}
		col++;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////Plastic Terms related to the variation of del gamma/////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////



void FSMicromorphic3DCurrConfigT:: Form_I1p_1()
{

    int row=0;
    int col=0;
    I1p_1=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
						{
							for(int i=0;i<3;i++)
							{
								for(int n=0;n<3;n++)
								{
									for(int m=0;m<3;m++)
									{
											for(int j=0;j<3;j++)
											{
												I1p_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)*fCauchy_stress_tensor_current_n_IP(j,k)
												*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,i)*fCauchy_stress_tensor_current_n_IP(m,n);
											}

									}
								}
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_I1p_2()
{

    int row=0;
    int col=0;
    I1p_2=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
						{
							for(int i=0;i<3;i++)
							{
								for(int n=0;n<3;n++)
								{
									for(int m=0;m<3;m++)
									{
											for(int j=0;j<3;j++)
											{
												I1p_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(p,p)*fCauchy_stress_tensor_current_n_IP(j,k)
												*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(m,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,i)*fCauchy_stress_tensor_current_n_IP(i,n);
											}

									}
								}
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I1p_3()
{

    int row=0;
    int col=0;
    I1p_3=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int i=0;i<3;i++)
						{
							for(int n=0;n<3;n++)
							{
								for(int m=0;m<3;m++)
								{
									for(int q=0;q<3;q++)
									{
										for(int j=0;j<3;j++)
										{
												I1p_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(i,i)*fCauchy_stress_tensor_current_n_IP(j,k)
												*fdFYdCauchy_Stress(m,q)*fCauchy_stress_tensor_current_n_IP(m,n)*fDeformation_Gradient_n(q,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,n);
										}
									}
								}
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}


void FSMicromorphic3DCurrConfigT:: Form_I1p_4()
{

    int row=0;
    int col=0;
    I1p_4=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int p=0;p<3;p++)
					{
						for(int K=0;K<3;K++)
						{
							for(int i=0;i<3;i++)
							{
								for(int n=0;n<3;n++)
								{
									for(int m=0;m<3;m++)
									{
											for(int j=0;j<3;j++)
											{
												I1p_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(i,i)*fCauchy_stress_tensor_current_n_IP(j,k)
												*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(p,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,p)*fIdentity_matrix(m,n);
											}

									}
								}
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_I1p_5()
{

    int row=0;
    int col=0;
    I1p_5=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
							for(int i=0;i<3;i++)
							{
								for(int n=0;n<3;n++)
								{
									for(int m=0;m<3;m++)
									{
											for(int j=0;j<3;j++)
											{
												I1p_5(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(i,i)*fCauchy_stress_tensor_current_n_IP(j,k)
												*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(m,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,n);
											}

									}
								}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}


void FSMicromorphic3DCurrConfigT:: Form_I1p_6()
{

    int row=0;
    int col=0;
    I1p_6=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int i=0;i<3;i++)
						{
							for(int n=0;n<3;n++)
							{
								for(int m=0;m<3;m++)
								{
									for(int j=0;j<3;j++)
									{
											I1p_6(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(i,i)*fCauchy_stress_tensor_current_n_IP(j,k)
											*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(n,K)*fDeformation_Gradient_Inverse(K,l)
											*fDeformation_Gradient_Inverse(I,m);
									}
								}
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}




void FSMicromorphic3DCurrConfigT:: Form_I2p_1()
{

    int row=0;
    int col=0;
    I2p_1=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
						{
							for(int i=0;i<3;i++)
							{
								for(int n=0;n<3;n++)
								{
									for(int m=0;m<3;m++)
									{
											for(int j=0;j<3;j++)
											{
												I2p_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,p)*fCauchy_stress_tensor_current_n_IP(p,k)
												*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,i)*fCauchy_stress_tensor_current_n_IP(m,n);
											}

									}
								}
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_I2p_2()
{

    int row=0;
    int col=0;
    I2p_2=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
						{
							for(int i=0;i<3;i++)
							{
								for(int n=0;n<3;n++)
								{
									for(int m=0;m<3;m++)
									{
										for(int j=0;j<3;j++)
										{
											I2p_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,p)*fCauchy_stress_tensor_current_n_IP(p,k)
											*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(m,K)*fDeformation_Gradient_Inverse(K,l)
											*fDeformation_Gradient_Inverse(I,i)*fCauchy_stress_tensor_current_n_IP(i,n);
										}
									}
								}
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_I2p_3()
{

    int row=0;
    int col=0;
    I2p_3=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
						{
							for(int n=0;n<3;n++)
							{
								for(int m=0;m<3;m++)
								{
									for(int q=0;q<3;q++)
									{
										for(int j=0;j<3;j++)
										{
												I2p_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,p)*fCauchy_stress_tensor_current_n_IP(p,k)
												*fdFYdCauchy_Stress(m,q)*fCauchy_stress_tensor_current_n_IP(m,n)*fDeformation_Gradient_n(q,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,n);
										}
									}
								}
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_I2p_4()
{

    int row=0;
    int col=0;
    I2p_4=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int p=0;p<3;p++)
					{
						for(int K=0;K<3;K++)
						{
							for(int q=0;q<3;q++)
							{
								for(int n=0;n<3;n++)
								{
									for(int m=0;m<3;m++)
									{
											for(int j=0;j<3;j++)
											{
												I2p_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)*fCauchy_stress_tensor_current_n_IP(q,k)
												*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(p,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,p)*fIdentity_matrix(m,n);
											}

									}
								}
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}




void FSMicromorphic3DCurrConfigT:: Form_I2p_5()
{

    int row=0;
    int col=0;
    I2p_5=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
							for(int p=0;p<3;p++)
							{
								for(int n=0;n<3;n++)
								{
									for(int m=0;m<3;m++)
									{
											for(int j=0;j<3;j++)
											{
												I2p_5(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,p)*fCauchy_stress_tensor_current_n_IP(p,k)
												*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(m,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,n);
											}

									}
								}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}


void FSMicromorphic3DCurrConfigT:: Form_I2p_6()
{

    int row=0;
    int col=0;
    I2p_6=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
						{
							for(int n=0;n<3;n++)
							{
								for(int m=0;m<3;m++)
								{
									for(int j=0;j<3;j++)
									{
											I2p_6(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,p)*fCauchy_stress_tensor_current_n_IP(p,k)
											*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(n,K)*fDeformation_Gradient_Inverse(K,l)
											*fDeformation_Gradient_Inverse(I,m);
									}
								}
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I3p_1()
{

    int row=0;
    int col=0;
    I3p_1=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
						{
							for(int i=0;i<3;i++)
							{
								for(int n=0;n<3;n++)
								{
									for(int m=0;m<3;m++)
									{
											for(int j=0;j<3;j++)
											{
												I3p_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,p)*fdGdCauchy_Stress_tr(k,p)
												*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,i)*fCauchy_stress_tensor_current_n_IP(m,n);
											}

									}
								}
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}


void FSMicromorphic3DCurrConfigT:: Form_I3p_2()
{

    int row=0;
    int col=0;
    I3p_2=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
						{
							for(int i=0;i<3;i++)
							{
								for(int n=0;n<3;n++)
								{
									for(int m=0;m<3;m++)
									{
											for(int j=0;j<3;j++)
											{
												I3p_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,p)*fdGdCauchy_Stress_tr(k,p)
												*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(m,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,i)*fCauchy_stress_tensor_current_n_IP(i,n);
											}

									}
								}
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}


void FSMicromorphic3DCurrConfigT:: Form_I3p_3()
{

    int row=0;
    int col=0;
    I3p_3=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
						{
							for(int n=0;n<3;n++)
							{
								for(int m=0;m<3;m++)
								{
									for(int q=0;q<3;q++)
									{
										for(int j=0;j<3;j++)
										{
												I3p_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,p)*fdGdCauchy_Stress_tr(k,p)
												*fdFYdCauchy_Stress(m,q)*fCauchy_stress_tensor_current_n_IP(m,n)*fDeformation_Gradient_n(q,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,n);
										}
									}
								}
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}


void FSMicromorphic3DCurrConfigT:: Form_I3p_4()
{

    int row=0;
    int col=0;
    I3p_4=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int p=0;p<3;p++)
					{
						for(int K=0;K<3;K++)
						{
							for(int q=0;q<3;q++)
							{
								for(int n=0;n<3;n++)
								{
									for(int m=0;m<3;m++)
									{
											for(int j=0;j<3;j++)
											{
												I3p_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,q)*fdGdCauchy_Stress_tr(k,q)
												*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(p,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,p)*fIdentity_matrix(m,n);
											}

									}
								}
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}




void FSMicromorphic3DCurrConfigT:: Form_I3p_5()
{

    int row=0;
    int col=0;
    I3p_5=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
							for(int p=0;p<3;p++)
							{
								for(int n=0;n<3;n++)
								{
									for(int m=0;m<3;m++)
									{
											for(int j=0;j<3;j++)
											{
												I3p_5(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,p)*fdGdCauchy_Stress_tr(k,p)
												*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(m,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,n);
											}

									}
								}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_I3p_6()
{

    int row=0;
    int col=0;
    I3p_6=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
						{
							for(int n=0;n<3;n++)
							{
								for(int m=0;m<3;m++)
								{
									for(int j=0;j<3;j++)
									{
											I3p_6(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,p)*fdGdCauchy_Stress_tr(k,p)
											*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(n,K)*fDeformation_Gradient_Inverse(K,l)
											*fDeformation_Gradient_Inverse(I,m);
									}
								}
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_I4p_1()
{

    int row=0;
    int col=0;
    I4p_1=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
						{
							for(int i=0;i<3;i++)
							{
								for(int n=0;n<3;n++)
								{
									for(int m=0;m<3;m++)
									{
											for(int j=0;j<3;j++)
											{
												I4p_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)*fIdentity_matrix(j,k)
												*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,i)*fCauchy_stress_tensor_current_n_IP(m,n);
											}

									}
								}
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}


void FSMicromorphic3DCurrConfigT:: Form_I4p_2()
{

    int row=0;
    int col=0;
    I4p_2=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
						{
							for(int i=0;i<3;i++)
							{
								for(int n=0;n<3;n++)
								{
									for(int m=0;m<3;m++)
									{
											for(int j=0;j<3;j++)
											{
												I4p_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(p,p)*fIdentity_matrix(j,k)
												*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(m,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,i)*fCauchy_stress_tensor_current_n_IP(i,n);
											}

									}
								}
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_I4p_3()
{

    int row=0;
    int col=0;
    I4p_3=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int i=0;i<3;i++)
						{
							for(int n=0;n<3;n++)
							{
								for(int m=0;m<3;m++)
								{
									for(int q=0;q<3;q++)
									{
										for(int j=0;j<3;j++)
										{
												I4p_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(i,i)*fIdentity_matrix(j,k)
												*fdFYdCauchy_Stress(m,q)*fCauchy_stress_tensor_current_n_IP(m,n)*fDeformation_Gradient_n(q,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,n);
										}
									}
								}
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}


void FSMicromorphic3DCurrConfigT:: Form_I4p_4()
{

    int row=0;
    int col=0;
    I4p_4=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int p=0;p<3;p++)
					{
						for(int K=0;K<3;K++)
						{
							for(int i=0;i<3;i++)
							{
								for(int n=0;n<3;n++)
								{
									for(int m=0;m<3;m++)
									{
											for(int j=0;j<3;j++)
											{
												I4p_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(i,i)*fIdentity_matrix(j,k)
												*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(p,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,p)*fIdentity_matrix(m,n);
											}

									}
								}
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}




void FSMicromorphic3DCurrConfigT:: Form_I4p_5()
{

    int row=0;
    int col=0;
    I4p_5=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
							for(int i=0;i<3;i++)
							{
								for(int n=0;n<3;n++)
								{
									for(int m=0;m<3;m++)
									{
											for(int j=0;j<3;j++)
											{
												I4p_5(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(i,i)*fIdentity_matrix(j,k)
												*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(m,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,n);
											}

									}
								}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_I4p_6()
{

    int row=0;
    int col=0;
    I4p_6=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int i=0;i<3;i++)
						{
							for(int n=0;n<3;n++)
							{
								for(int m=0;m<3;m++)
								{
									for(int j=0;j<3;j++)
									{
											I4p_6(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(i,i)*fIdentity_matrix(j,k)
											*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(n,K)*fDeformation_Gradient_Inverse(K,l)
											*fDeformation_Gradient_Inverse(I,m);
									}
								}
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_I5p_1()
{

    int row=0;
    int col=0;
    I5p_1=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int i=0;i<3;i++)
							{
								for(int n=0;n<3;n++)
								{
									for(int m=0;m<3;m++)
									{
											for(int j=0;j<3;j++)
											{
												I5p_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,k)
												*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,i)*fCauchy_stress_tensor_current_n_IP(m,n);
											}

									}
								}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_I5p_2()
{

    int row=0;
    int col=0;
    I5p_2=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int i=0;i<3;i++)
							{
								for(int n=0;n<3;n++)
								{
									for(int m=0;m<3;m++)
									{
											for(int j=0;j<3;j++)
											{
												I5p_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,k)
												*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(m,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,i)*fCauchy_stress_tensor_current_n_IP(i,n);
											}

									}
								}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_I5p_3()
{

    int row=0;
    int col=0;
    I5p_3=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int n=0;n<3;n++)
							{
								for(int m=0;m<3;m++)
								{
									for(int q=0;q<3;q++)
									{
										for(int j=0;j<3;j++)
										{
												I5p_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,k)
												*fdFYdCauchy_Stress(m,q)*fCauchy_stress_tensor_current_n_IP(m,n)*fDeformation_Gradient_n(q,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,n);
										}
									}
								}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_I5p_4()
{

    int row=0;
    int col=0;
    I5p_4=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int p=0;p<3;p++)
					{
						for(int K=0;K<3;K++)
						{
								for(int n=0;n<3;n++)
								{
									for(int m=0;m<3;m++)
									{
											for(int j=0;j<3;j++)
											{
												I5p_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,k)
												*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(p,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,p)*fIdentity_matrix(m,n);
											}

									}
								}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}




void FSMicromorphic3DCurrConfigT:: Form_I5p_5()
{

    int row=0;
    int col=0;
    I5p_5=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
								for(int n=0;n<3;n++)
								{
									for(int m=0;m<3;m++)
									{
											for(int j=0;j<3;j++)
											{
												I5p_5(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,k)
												*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(m,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,n);
											}

									}
								}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_I5p_6()
{

    int row=0;
    int col=0;
    I5p_6=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
							for(int n=0;n<3;n++)
							{
								for(int m=0;m<3;m++)
								{
									for(int j=0;j<3;j++)
									{
											I5p_6(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,k)
											*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(n,K)*fDeformation_Gradient_Inverse(K,l)
											*fDeformation_Gradient_Inverse(I,m);
									}
								}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_I6p_1()
{

    int row=0;
    int col=0;
    I6p_1=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int i=0;i<3;i++)
							{
								for(int n=0;n<3;n++)
								{
									for(int m=0;m<3;m++)
									{
											for(int j=0;j<3;j++)
											{
												I6p_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(k,j)
												*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,i)*fCauchy_stress_tensor_current_n_IP(m,n);
											}

									}
								}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_I6p_2()
{

    int row=0;
    int col=0;
    I6p_2=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int i=0;i<3;i++)
							{
								for(int n=0;n<3;n++)
								{
									for(int m=0;m<3;m++)
									{
											for(int j=0;j<3;j++)
											{
												I6p_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(k,j)
												*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(m,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,i)*fCauchy_stress_tensor_current_n_IP(i,n);
											}

									}
								}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_I6p_3()
{

    int row=0;
    int col=0;
    I6p_3=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int n=0;n<3;n++)
							{
								for(int m=0;m<3;m++)
								{
									for(int q=0;q<3;q++)
									{
										for(int j=0;j<3;j++)
										{
												I6p_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(k,j)
												*fdFYdCauchy_Stress(m,q)*fCauchy_stress_tensor_current_n_IP(m,n)*fDeformation_Gradient_n(q,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,n);
										}
									}
								}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_I6p_4()
{

    int row=0;
    int col=0;
    I6p_4=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int p=0;p<3;p++)
					{
						for(int K=0;K<3;K++)
						{
								for(int n=0;n<3;n++)
								{
									for(int m=0;m<3;m++)
									{
											for(int j=0;j<3;j++)
											{
												I6p_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(k,j)
												*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(p,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,p)*fIdentity_matrix(m,n);
											}

									}
								}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}




void FSMicromorphic3DCurrConfigT:: Form_I6p_5()
{

    int row=0;
    int col=0;
    I6p_5=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
								for(int n=0;n<3;n++)
								{
									for(int m=0;m<3;m++)
									{
											for(int j=0;j<3;j++)
											{
												I6p_5(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(k,j)
												*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(m,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,n);
											}

									}
								}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_I6p_6()
{

    int row=0;
    int col=0;
    I6p_6=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
							for(int n=0;n<3;n++)
							{
								for(int m=0;m<3;m++)
								{
									for(int j=0;j<3;j++)
									{
											I6p_6(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(k,j)
											*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(n,K)*fDeformation_Gradient_Inverse(K,l)
											*fDeformation_Gradient_Inverse(I,m);
									}
								}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FSMicromorphic3DCurrConfigT:: Form_I7p_1()
{

    int row=0;
    int col=0;
    I7p_1=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int i=0;i<3;i++)
						{
							for(int n=0;n<3;n++)
							{
								for(int j=0;j<3;j++)
								{
									I7p_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(K,l)
									*fDeformation_Gradient_Inverse(I,i)*fdev_Cauchy_stress_tensor_current_n_IP(j,n)
									*fCauchy_stress_tensor_current_n_IP(n,k);
								}
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I7p_2()
{

    int row=0;
    int col=0;
    I7p_2=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int i=0;i<3;i++)
							{
								for(int n=0;n<3;n++)
									{
										for(int j=0;j<3;j++)
											{
												I7p_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fDeformation_Gradient_n(j,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,i)*fCauchy_stress_tensor_current_n_IP(i,n)*fCauchy_stress_tensor_current_n_IP(n,k);
											}

									}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I7p_3()
{

    int row=0;
    int col=0;
    I7p_3=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int i=0;i<3;i++)
							{
								for(int p=0;p<3;p++)
									{
										for(int j=0;j<3;j++)
											{
												for(int n=0;n<3;n++)
													{
														I7p_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(K,l)
														*fDeformation_Gradient_Inverse(I,p)*fCauchy_stress_tensor_current_n_IP(p,i)
														*fIdentity_matrix(j,n)*fCauchy_stress_tensor_current_n_IP(n,k);
													}
											}

									}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I7p_4()
{

    int row=0;
    int col=0;
    I7p_4=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int i=0;i<3;i++)
							{
							for(int j=0;j<3;j++)
								{
								for(int n=0;n<3;n++)
									{
										I7p_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,i)
										*fDeformation_Gradient_n(n,K)*fDeformation_Gradient_Inverse(K,l)
										*fDeformation_Gradient_Inverse(I,i)*fCauchy_stress_tensor_current_n_IP(n,k);
									}
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I7p_5()
{

    int row=0;
    int col=0;
    I7p_5=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int i=0;i<3;i++)
							{
								for(int p=0;p<3;p++)
									{
										for(int j=0;j<3;j++)
											{
												for(int n=0;n<3;n++)
													{
														I7p_5(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(p,i)
														*fDeformation_Gradient_n(p,K)*fDeformation_Gradient_Inverse(K,l)
														*fDeformation_Gradient_Inverse(I,i)*fIdentity_matrix(j,n)*fCauchy_stress_tensor_current_n_IP(n,k);
													}
											}

									}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I7p_6()
{

    int row=0;
    int col=0;
    I7p_6=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							for(int n=0;n<3;n++)
							{
								I7p_6(row,col)+=fDeformation_Gradient_Inverse(L,j)*fDeformation_Gradient_n(j,K)
								*fDeformation_Gradient_Inverse(K,l)*fDeformation_Gradient_Inverse(I,n)*fCauchy_stress_tensor_current_n_IP(n,k);
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I7p_7()
{

    int row=0;
    int col=0;
    I7p_7=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
							for(int n=0;n<3;n++)
								{
								for(int j=0;j<3;j++)
									{
										I7p_7(row,col)+=fDeformation_Gradient_Inverse(L,j)*fDeformation_Gradient_n(p,K)
										*fDeformation_Gradient_Inverse(K,l)*fDeformation_Gradient_Inverse(I,p)*fIdentity_matrix(j,n)
										*fCauchy_stress_tensor_current_n_IP(n,k);
									}
								}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I7p_8()
{

    int row=0;
    int col=0;
    I7p_8=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int n=0;n<3;n++)
							{
								for(int j=0;j<3;j++)
									{
										I7p_8(row,col)+=fDeformation_Gradient_Inverse(L,j)*fDeformation_Gradient_n(k,K)
										*fDeformation_Gradient_Inverse(K,l)*fDeformation_Gradient_Inverse(I,n)
										*fCauchy_stress_tensor_current_n_IP(n,j);
									}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I7p_9()
{

    int row=0;
    int col=0;
    I7p_9=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int n=0;n<3;n++)
							{
								for(int j=0;j<3;j++)
									{
									for(int i=0;i<3;i++)
										{
											I7p_9(row,col)+=fDeformation_Gradient_Inverse(L,j)*fDeformation_Gradient_n(i,K)
											*fDeformation_Gradient_Inverse(K,l)*fDeformation_Gradient_Inverse(I,i)*dev_Cauchy_stress_tr(j,n)
											*fCauchy_stress_tensor_current_n_IP(n,k);
										}
									}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I7p_10()
{

    int row=0;
    int col=0;
    I7p_10=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int n=0;n<3;n++)
							{
								for(int j=0;j<3;j++)
									{
									for(int i=0;i<3;i++)
										{
										for(int a=0;a<3;a++)
											{
											for(int b=0;b<3;b++)
												{
													I7p_10(row,col)+=fDeformation_Gradient_Inverse(L,j)*fDeformation_Gradient_n(a,K)
													*fDeformation_Gradient_Inverse(K,l)*fDeformation_Gradient_Inverse(I,i)*fCauchy_stress_tensor_current_n_IP(i,b)
													*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(j,n)
													*fCauchy_stress_tensor_current_n_IP(n,k);
												}
											}
										}
									}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I7p_11()
{

    int row=0;
    int col=0;
    I7p_11=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int n=0;n<3;n++)
							{
								for(int j=0;j<3;j++)
									{
									for(int i=0;i<3;i++)
										{
										for(int a=0;a<3;a++)
											{
											for(int b=0;b<3;b++)
												{
													I7p_11(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(a,i)
													*fDeformation_Gradient_n(b,K)*fDeformation_Gradient_Inverse(K,l)*fDeformation_Gradient_Inverse(I,i)
													*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(j,n)
													*fCauchy_stress_tensor_current_n_IP(n,k);
												}
											}
										}
									}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I7p_12()
{

    int row=0;
    int col=0;
    I7p_12=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int n=0;n<3;n++)
							{
								for(int j=0;j<3;j++)
									{
										for(int a=0;a<3;a++)
											{
											for(int b=0;b<3;b++)
												{
													I7p_12(row,col)+=fDeformation_Gradient_Inverse(L,j)*fDeformation_Gradient_n(a,K)
													*fDeformation_Gradient_Inverse(K,l)*fDeformation_Gradient_Inverse(I,b)
													*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(j,n)
													*fCauchy_stress_tensor_current_n_IP(n,k);
												}
											}
									}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I7p_13()
{

    int row=0;
    int col=0;
    I7p_13=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int n=0;n<3;n++)
						{
							for(int j=0;j<3;j++)
							{
								for(int a=0;a<3;a++)
								{
									for(int b=0;b<3;b++)
									{
										I7p_13(row,col)+=fDeformation_Gradient_Inverse(L,j)*fDeformation_Gradient_n(b,K)
										*fDeformation_Gradient_Inverse(K,l)*fDeformation_Gradient_Inverse(I,a)
										*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(j,n)
										*fCauchy_stress_tensor_current_n_IP(n,k);
									}
								}
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
void FSMicromorphic3DCurrConfigT:: Form_I8p_1()
{

    int row=0;
    int col=0;
    I8p_1=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int i=0;i<3;i++)
							{
								for(int n=0;n<3;n++)
									{
										for(int j=0;j<3;j++)
											{
												I8p_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,n)
												*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,i)*fdev_Cauchy_stress_tensor_current_n_IP(k,n);
											}
									}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I8p_2()
{

    int row=0;
    int col=0;
    I8p_2=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int i=0;i<3;i++)
							{
								for(int n=0;n<3;n++)
									{
										for(int j=0;j<3;j++)
											{
												I8p_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,n)
												*fDeformation_Gradient_n(k,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,i)*fCauchy_stress_tensor_current_n_IP(i,n);
											}

									}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I8p_3()
{

    int row=0;
    int col=0;
    I8p_3=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int i=0;i<3;i++)
							{
								for(int p=0;p<3;p++)
									{
										for(int j=0;j<3;j++)
											{
												for(int n=0;n<3;n++)
													{
														I8p_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,n)
														*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(K,l)
														*fDeformation_Gradient_Inverse(I,p)*fCauchy_stress_tensor_current_n_IP(p,i)
														*fIdentity_matrix(k,n);
													}
											}

									}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I8p_4()
{

    int row=0;
    int col=0;
    I8p_4=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int i=0;i<3;i++)
							{
							for(int j=0;j<3;j++)
								{
								for(int n=0;n<3;n++)
									{
										I8p_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,n)
										*fCauchy_stress_tensor_current_n_IP(k,i)*fDeformation_Gradient_n(n,K)*fDeformation_Gradient_Inverse(K,l)
										*fDeformation_Gradient_Inverse(I,i);
									}
								}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I8p_5()
{

    int row=0;
    int col=0;
    I8p_5=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int i=0;i<3;i++)
							{
								for(int p=0;p<3;p++)
									{
										for(int j=0;j<3;j++)
											{
												for(int n=0;n<3;n++)
													{
														I8p_5(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,n)
														*fCauchy_stress_tensor_current_n_IP(p,i)*fDeformation_Gradient_n(p,K)*fDeformation_Gradient_Inverse(K,l)
														*fDeformation_Gradient_Inverse(I,i)*fIdentity_matrix(k,n);
													}
											}

									}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I8p_6()
{

    int row=0;
    int col=0;
    I8p_6=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
							{
								for(int n=0;n<3;n++)
								{
									I8p_6(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,n)*fDeformation_Gradient_n(k,K)
									*fDeformation_Gradient_Inverse(K,l)*fDeformation_Gradient_Inverse(I,n);
								}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I8p_7()
{

    int row=0;
    int col=0;
    I8p_7=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
							for(int n=0;n<3;n++)
								{
								for(int j=0;j<3;j++)
									{
										I8p_7(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,n)*fDeformation_Gradient_n(p,K)
										*fDeformation_Gradient_Inverse(K,l)*fDeformation_Gradient_Inverse(I,p)*fIdentity_matrix(k,n)
										;
									}
								}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I8p_8()
{

    int row=0;
    int col=0;
    I8p_8=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int n=0;n<3;n++)
							{
								for(int j=0;j<3;j++)
									{
										I8p_8(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,n)*fDeformation_Gradient_n(n,K)
										*fDeformation_Gradient_Inverse(K,l)*fDeformation_Gradient_Inverse(I,k);
									}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I8p_9()
{

    int row=0;
    int col=0;
    I8p_9=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int n=0;n<3;n++)
							{
								for(int j=0;j<3;j++)
									{
									for(int i=0;i<3;i++)
										{
											I8p_9(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,n)
											*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(K,l)
											*fDeformation_Gradient_Inverse(I,i)*dev_Cauchy_stress_tr(k,n);
										}
									}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I8p_10()
{

    int row=0;
    int col=0;
    I8p_10=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int n=0;n<3;n++)
							{
								for(int j=0;j<3;j++)
									{
									for(int i=0;i<3;i++)
										{
										for(int a=0;a<3;a++)
											{
											for(int b=0;b<3;b++)
												{
													I8p_10(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,n)*fDeformation_Gradient_n(a,K)
													*fDeformation_Gradient_Inverse(K,l)*fDeformation_Gradient_Inverse(I,i)*fCauchy_stress_tensor_current_n_IP(i,b)
													*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(k,n);
												}
											}
										}
									}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I8p_11()
{

    int row=0;
    int col=0;
    I8p_11=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int n=0;n<3;n++)
							{
								for(int j=0;j<3;j++)
									{
									for(int i=0;i<3;i++)
										{
										for(int a=0;a<3;a++)
											{
											for(int b=0;b<3;b++)
												{
													I8p_11(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,n)
													*fCauchy_stress_tensor_current_n_IP(a,i)*fDeformation_Gradient_n(b,K)
													*fDeformation_Gradient_Inverse(K,l)*fDeformation_Gradient_Inverse(I,i)
													*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(k,n);
												}
											}
										}
									}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I8p_12()
{

    int row=0;
    int col=0;
    I8p_12=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int n=0;n<3;n++)
							{
								for(int j=0;j<3;j++)
									{
										for(int a=0;a<3;a++)
											{
											for(int b=0;b<3;b++)
												{
													I8p_12(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,n)
													*fDeformation_Gradient_n(a,K)
													*fDeformation_Gradient_Inverse(K,l)*fDeformation_Gradient_Inverse(I,b)
													*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(k,n);
												}
											}
									}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I8p_13()
{

    int row=0;
    int col=0;
    I8p_13=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int n=0;n<3;n++)
							{
								for(int j=0;j<3;j++)
									{
										for(int a=0;a<3;a++)
											{
											for(int b=0;b<3;b++)
												{
													I8p_13(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,n)
													*fDeformation_Gradient_n(b,K)
													*fDeformation_Gradient_Inverse(K,l)*fDeformation_Gradient_Inverse(I,a)
													*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(k,n);
												}
											}
									}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
void FSMicromorphic3DCurrConfigT:: Form_I9p_1()
{

    int row=0;
    int col=0;
    I9p_1=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int i=0;i<3;i++)
						{
							for(int j=0;j<3;j++)
							{
								I9p_1(row,col)+=fDeformation_Gradient_Inverse(L,j)
								*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(K,l)
								*fDeformation_Gradient_Inverse(I,i)*fdev_Cauchy_stress_tensor_current_n_IP(j,k);
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I9p_2()
{

    int row=0;
    int col=0;
    I9p_2=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int i=0;i<3;i++)
						{
								for(int j=0;j<3;j++)
								{
									I9p_2(row,col)+=fDeformation_Gradient_Inverse(L,j)
									*fDeformation_Gradient_n(j,K)*fDeformation_Gradient_Inverse(K,l)
									*fDeformation_Gradient_Inverse(I,i)*fCauchy_stress_tensor_current_n_IP(i,k);
								}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I9p_3()
{

    int row=0;
    int col=0;
    I9p_3=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int i=0;i<3;i++)
							{
								for(int p=0;p<3;p++)
									{
										for(int j=0;j<3;j++)
											{
												I9p_3(row,col)+=fDeformation_Gradient_Inverse(L,j)
												*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,p)*fCauchy_stress_tensor_current_n_IP(p,i)
												*fIdentity_matrix(j,k);
											}
									}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I9p_4()
{

    int row=0;
    int col=0;
    I9p_4=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int i=0;i<3;i++)
							{
							for(int j=0;j<3;j++)
								{
									I9p_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,i)
									*fDeformation_Gradient_n(k,K)*fDeformation_Gradient_Inverse(K,l)
									*fDeformation_Gradient_Inverse(I,i);
								}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I9p_5()
{

    int row=0;
    int col=0;
    I9p_5=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int i=0;i<3;i++)
							{
								for(int p=0;p<3;p++)
									{
										for(int j=0;j<3;j++)
											{
												I9p_5(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(p,i)
												*fDeformation_Gradient_n(p,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,i)*fIdentity_matrix(j,k);
											}
									}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I9p_6()
{

    int row=0;
    int col=0;
    I9p_6=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
							{
								I9p_6(row,col)+=fDeformation_Gradient_Inverse(L,j)*fDeformation_Gradient_n(j,K)
								*fDeformation_Gradient_Inverse(K,l)*fDeformation_Gradient_Inverse(I,k);
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I9p_7()
{

    int row=0;
    int col=0;
    I9p_7=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
							for(int j=0;j<3;j++)
								{
									I9p_7(row,col)+=fDeformation_Gradient_Inverse(L,j)*fDeformation_Gradient_n(p,K)
									*fDeformation_Gradient_Inverse(K,l)*fDeformation_Gradient_Inverse(I,p)*fIdentity_matrix(j,k);
								}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I9p_8()
{

    int row=0;
    int col=0;
    I9p_8=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
							{
								I9p_8(row,col)+=fDeformation_Gradient_Inverse(L,j)*fDeformation_Gradient_n(k,K)
								*fDeformation_Gradient_Inverse(K,l)*fDeformation_Gradient_Inverse(I,j);
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I9p_9()
{

    int row=0;
    int col=0;
    I9p_9=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							for(int i=0;i<3;i++)
								{
										I9p_9(row,col)+=fDeformation_Gradient_Inverse(L,j)
										*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(K,l)
										*fDeformation_Gradient_Inverse(I,i)*dev_Cauchy_stress_tr(j,k);
								}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I9p_10()
{

    int row=0;
    int col=0;
    I9p_10=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
							{
							for(int i=0;i<3;i++)
								{
									for(int a=0;a<3;a++)
										{
											for(int b=0;b<3;b++)
												{
													I9p_10(row,col)+=fDeformation_Gradient_Inverse(L,j)*fDeformation_Gradient_n(a,K)
													*fDeformation_Gradient_Inverse(K,l)*fDeformation_Gradient_Inverse(I,i)*fCauchy_stress_tensor_current_n_IP(i,b)
													*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(j,k);
												}
										}
								}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I9p_11()
{

    int row=0;
    int col=0;
    I9p_11=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							for(int i=0;i<3;i++)
							{
								for(int a=0;a<3;a++)
								{
									for(int b=0;b<3;b++)
									{
										I9p_11(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(a,i)
										*fDeformation_Gradient_n(b,K)
										*fDeformation_Gradient_Inverse(K,l)*fDeformation_Gradient_Inverse(I,i)
										*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(j,k);
									}
								}
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I9p_12()
{

    int row=0;
    int col=0;
    I9p_12=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
							for(int j=0;j<3;j++)
							{
									for(int a=0;a<3;a++)
									{
										for(int b=0;b<3;b++)
										{
											I9p_12(row,col)+=fDeformation_Gradient_Inverse(L,j)*fDeformation_Gradient_n(a,K)
											*fDeformation_Gradient_Inverse(K,l)*fDeformation_Gradient_Inverse(I,b)
											*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(j,k);
										}
									}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I9p_13()
{

    int row=0;
    int col=0;
    I9p_13=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
							for(int j=0;j<3;j++)
							{
									for(int a=0;a<3;a++)
									{
										for(int b=0;b<3;b++)
										{
											I9p_13(row,col)+=fDeformation_Gradient_Inverse(L,j)*fDeformation_Gradient_n(b,K)
											*fDeformation_Gradient_Inverse(K,l)*fDeformation_Gradient_Inverse(I,a)
											*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(j,k);
										}
									}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
void FSMicromorphic3DCurrConfigT:: Form_I10p_1()
{

    int row=0;
    int col=0;
    I10p_1=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int i=0;i<3;i++)
						{
							for(int j=0;j<3;j++)
							{
								I10p_1(row,col)+=fDeformation_Gradient_Inverse(L,j)
								*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(K,l)
								*fDeformation_Gradient_Inverse(I,i)*fdev_Cauchy_stress_tensor_current_n_IP(k,j);
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I10p_2()
{

    int row=0;
    int col=0;
    I10p_2=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int i=0;i<3;i++)
						{
								for(int j=0;j<3;j++)
								{
									I10p_2(row,col)+=fDeformation_Gradient_Inverse(L,j)
									*fDeformation_Gradient_n(k,K)*fDeformation_Gradient_Inverse(K,l)
									*fDeformation_Gradient_Inverse(I,i)*fCauchy_stress_tensor_current_n_IP(i,j);
								}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I10p_3()
{

    int row=0;
    int col=0;
    I10p_3=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int i=0;i<3;i++)
							{
								for(int p=0;p<3;p++)
									{
										for(int j=0;j<3;j++)
											{
												I10p_3(row,col)+=fDeformation_Gradient_Inverse(L,j)
												*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,p)*fCauchy_stress_tensor_current_n_IP(p,i)
												*fIdentity_matrix(k,j);
											}
									}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I10p_4()
{

    int row=0;
    int col=0;
    I10p_4=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
					for(int i=0;i<3;i++)
						{
							for(int j=0;j<3;j++)
								{
										I10p_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(k,i)
										*fDeformation_Gradient_n(j,K)*fDeformation_Gradient_Inverse(K,l)
										*fDeformation_Gradient_Inverse(I,i);
								}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I10p_5()
{

    int row=0;
    int col=0;
    I10p_5=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int i=0;i<3;i++)
							{
								for(int p=0;p<3;p++)
									{
										for(int j=0;j<3;j++)
											{
														I10p_5(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(p,i)
														*fDeformation_Gradient_n(p,K)*fDeformation_Gradient_Inverse(K,l)
														*fDeformation_Gradient_Inverse(I,i)*fIdentity_matrix(k,j);
											}
									}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I10p_6()
{

    int row=0;
    int col=0;
    I10p_6=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							I10p_6(row,col)+=fDeformation_Gradient_Inverse(L,j)*fDeformation_Gradient_n(k,K)
							*fDeformation_Gradient_Inverse(K,l)*fDeformation_Gradient_Inverse(I,j);
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I10p_7()
{

    int row=0;
    int col=0;
    I10p_7=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
						{
							for(int j=0;j<3;j++)
							{
									I10p_7(row,col)+=fDeformation_Gradient_Inverse(L,j)*fDeformation_Gradient_n(p,K)
									*fDeformation_Gradient_Inverse(K,l)*fDeformation_Gradient_Inverse(I,p)*fIdentity_matrix(k,j);
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I10p_8()
{

    int row=0;
    int col=0;
    I10p_8=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
							{
								I10p_8(row,col)+=fDeformation_Gradient_Inverse(L,j)*fDeformation_Gradient_n(j,K)
								*fDeformation_Gradient_Inverse(K,l)*fDeformation_Gradient_Inverse(I,k);
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I10p_9()
{

    int row=0;
    int col=0;
    I10p_9=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							for(int i=0;i<3;i++)
								{
										I10p_9(row,col)+=fDeformation_Gradient_Inverse(L,j)
										*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(K,l)
										*fDeformation_Gradient_Inverse(I,i)*dev_Cauchy_stress_tr(k,j);
								}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I10p_10()
{

    int row=0;
    int col=0;
    I10p_10=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
							{
							for(int i=0;i<3;i++)
								{
									for(int a=0;a<3;a++)
										{
											for(int b=0;b<3;b++)
												{
													I10p_10(row,col)+=fDeformation_Gradient_Inverse(L,j)*fDeformation_Gradient_n(a,K)
													*fDeformation_Gradient_Inverse(K,l)*fDeformation_Gradient_Inverse(I,i)*fCauchy_stress_tensor_current_n_IP(i,b)
													*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(k,j);
												}
										}
								}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I10p_11()
{

    int row=0;
    int col=0;
    I10p_11=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							for(int i=0;i<3;i++)
							{
								for(int a=0;a<3;a++)
								{
									for(int b=0;b<3;b++)
									{
										I10p_11(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(a,i)
										*fDeformation_Gradient_n(b,K)
										*fDeformation_Gradient_Inverse(K,l)*fDeformation_Gradient_Inverse(I,i)
										*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(k,j);
									}
								}
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I10p_12()
{

    int row=0;
    int col=0;
    I10p_12=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
							for(int j=0;j<3;j++)
							{
									for(int a=0;a<3;a++)
									{
										for(int b=0;b<3;b++)
										{
											I10p_12(row,col)+=fDeformation_Gradient_Inverse(L,j)*fDeformation_Gradient_n(a,K)
											*fDeformation_Gradient_Inverse(K,l)*fDeformation_Gradient_Inverse(I,b)
											*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(k,j);
										}
									}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I10p_13()
{

    int row=0;
    int col=0;
    I10p_13=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							for(int a=0;a<3;a++)
							{
								for(int b=0;b<3;b++)
								{
									I10p_13(row,col)+=fDeformation_Gradient_Inverse(L,j)*fDeformation_Gradient_n(b,K)
									*fDeformation_Gradient_Inverse(K,l)*fDeformation_Gradient_Inverse(I,a)
									*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(k,j);
								}
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_1()
{

    int row=0;
    int col=0;
    I_temp_11p_1=0.0;
    fTemp_matrix_nsd_x_nsd = 0.0;
    fTemp_matrix_nsd_x_nsd2 = 0.0;
    fTemp_matrix_nsd_x_nsd3 = 0.0;


		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int i=0;i<3;i++)
				{
					for(int p=0;p<3;p++)
					{
						for(int n=0;n<3;n++)
						{
							for(int m=0;m<3;m++)
							{
								fTemp_matrix_nsd_x_nsd(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(i,K)
								*fDeformation_Gradient_Inverse(I,i)*fdev_Cauchy_stress_tensor_current_n_IP(m,p)
								*fCauchy_stress_tensor_current_n_IP(p,n);
							}
						}
					}
				}
			row++;
			}
		col++;
		}

		row=0;
		col=0;

		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int i=0;i<3;i++)
				{
					for(int p=0;p<3;p++)
					{
						for(int n=0;n<3;n++)
						{
							for(int m=0;m<3;m++)
							{
								fTemp_matrix_nsd_x_nsd2(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(i,K)
								*fDeformation_Gradient_Inverse(I,i)*fCauchy_stress_tensor_current_n_IP(m,p)*fdev_Cauchy_stress_tensor_current_n_IP(n,p);
							}
						}
					}
				}
			row++;
			}
		col++;
		}

		row=0;
		col=0;

		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int i=0;i<3;i++)
				{
					for(int n=0;n<3;n++)
					{
						for(int m=0;m<3;m++)
						{
							fTemp_matrix_nsd_x_nsd3(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(i,K)
							*fDeformation_Gradient_Inverse(I,i)*fdev_Cauchy_stress_tensor_current_n_IP(m,n);
						}
					}
				}
			row++;
			}
		col++;
		}
		fTemp_matrix_nsd_x_nsd3*= fMaterial_Params[kMu];
		row=0;
		col=0;

		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int i=0;i<3;i++)
				{
					for(int n=0;n<3;n++)
					{
						for(int m=0;m<3;m++)
						{
							I_temp_11p_1(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(i,K)
							*fDeformation_Gradient_Inverse(I,i)*fdev_Cauchy_stress_tensor_current_n_IP(n,m);
						}
					}
				}
			row++;
			}
		col++;
		}
		I_temp_11p_1*= fMaterial_Params[kMu];

		I_temp_11p_1+= fTemp_matrix_nsd_x_nsd;
		I_temp_11p_1+= fTemp_matrix_nsd_x_nsd2;
		I_temp_11p_1+= fTemp_matrix_nsd_x_nsd3;
}

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_2()
{

    int row=0;
    int col=0;
    I_temp_11p_2=0.0;
    fTemp_matrix_nsd_x_nsd = 0.0;
    fTemp_matrix_nsd_x_nsd2 = 0.0;
    fTemp_matrix_nsd_x_nsd3 = 0.0;


		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int i=0;i<3;i++)
				{
					for(int p=0;p<3;p++)
					{
						for(int n=0;n<3;n++)
						{
							for(int m=0;m<3;m++)
							{
								fTemp_matrix_nsd_x_nsd(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_Inverse(I,i)
								*fDeformation_Gradient_n(m,K)*fCauchy_stress_tensor_current_n_IP(i,p)
								*fCauchy_stress_tensor_current_n_IP(p,n);
							}
						}
					}
				}
			row++;
			}
		col++;
		}



		row=0;
		col=0;

		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int i=0;i<3;i++)
				{
					for(int p=0;p<3;p++)
					{
						for(int n=0;n<3;n++)
						{
							for(int m=0;m<3;m++)
							{
								fTemp_matrix_nsd_x_nsd2(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_Inverse(I,i)
								*fCauchy_stress_tensor_current_n_IP(m,p)*fDeformation_Gradient_n(n,K)*fCauchy_stress_tensor_current_n_IP(i,p);
							}
						}
					}
				}
			row++;
			}
		col++;
		}

		row=0;
		col=0;



		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int i=0;i<3;i++)
				{
					for(int n=0;n<3;n++)
					{
						for(int m=0;m<3;m++)
						{
							fTemp_matrix_nsd_x_nsd3(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_Inverse(I,i)*fDeformation_Gradient_n(m,K)
							*fCauchy_stress_tensor_current_n_IP(i,n);
						}
					}
				}
			row++;
			}
		col++;
		}
		fTemp_matrix_nsd_x_nsd3*= fMaterial_Params[kMu];
		row=0;
		col=0;


		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int i=0;i<3;i++)
				{
					for(int n=0;n<3;n++)
					{
						for(int m=0;m<3;m++)
						{
							I_temp_11p_2(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_Inverse(I,i)*fDeformation_Gradient_n(n,K)
							*fCauchy_stress_tensor_current_n_IP(i,m);
						}
					}
				}
			row++;
			}
		col++;
		}
		I_temp_11p_2*= fMaterial_Params[kMu];

		I_temp_11p_2+= fTemp_matrix_nsd_x_nsd;
		I_temp_11p_2+= fTemp_matrix_nsd_x_nsd2;
		I_temp_11p_2+= fTemp_matrix_nsd_x_nsd3;
}

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_3()
{

    int row=0;
    int col=0;
    I_temp_11p_3=0.0;
    fTemp_matrix_nsd_x_nsd = 0.0;


		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int i=0;i<3;i++)
				{
					for(int p=0;p<3;p++)
					{
						for(int n=0;n<3;n++)
						{
							for(int m=0;m<3;m++)
							{
								fTemp_matrix_nsd_x_nsd(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(I,p)
								*fCauchy_stress_tensor_current_n_IP(p,i)*fCauchy_stress_tensor_current_n_IP(m,n);
							}
						}
					}
				}
			row++;
			}
		col++;
		}
		fTemp_matrix_nsd_x_nsd*= 2;

		row=0;
		col=0;

		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int i=0;i<3;i++)
				{
					for(int p=0;p<3;p++)
					{
						for(int n=0;n<3;n++)
						{
							for(int m=0;m<3;m++)
							{
								I_temp_11p_3(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(I,p)
								*fCauchy_stress_tensor_current_n_IP(p,i)*fIdentity_matrix(m,n);
							}
						}
					}
				}
			row++;
			}
		col++;
		}
		I_temp_11p_3*= 2;
		I_temp_11p_3*= fMaterial_Params[kMu];
		I_temp_11p_3+= fTemp_matrix_nsd_x_nsd;
}

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_4()
{

    int row=0;
    int col=0;
    I_temp_11p_4=0.0;
    I_temp_11p_4_Transpose = 0.0;
    fTemp_matrix_nsd_x_nsd2 = 0.0;
    fTemp_matrix_nsd_x_nsd3 = 0.0;


		for(int I=0;I<3;I++)
		{
			row=0;
			for(int K=0;K<3;K++)
			{
				for(int i=0;i<3;i++)
				{
					for(int p=0;p<3;p++)
					{
						for(int n=0;n<3;n++)
						{
							for(int m=0;m<3;m++)
							{
								I_temp_11p_4_Transpose(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_Inverse(I,i)
								*fCauchy_stress_tensor_current_n_IP(m,i)*fDeformation_Gradient_n(p,K)*fCauchy_stress_tensor_current_n_IP(p,n);
							}
						}
					}
				}
			row++;
			}
		col++;
		}

		row=0;
		col=0;


		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int i=0;i<3;i++)
				{
					for(int p=0;p<3;p++)
					{
						for(int n=0;n<3;n++)
						{
							for(int m=0;m<3;m++)
							{
								fTemp_matrix_nsd_x_nsd2(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_Inverse(I,i)
								*fCauchy_stress_tensor_current_n_IP(m,p)*fCauchy_stress_tensor_current_n_IP(n,i)*fDeformation_Gradient_n(p,K);
							}
						}
					}
				}
			row++;
			}
		col++;
		}

		row=0;
		col=0;

		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int i=0;i<3;i++)
				{
					for(int n=0;n<3;n++)
					{
						for(int m=0;m<3;m++)
						{
							fTemp_matrix_nsd_x_nsd3(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_Inverse(I,i)*fCauchy_stress_tensor_current_n_IP(m,i)
							*fDeformation_Gradient_n(n,K);
						}
					}
				}
			row++;
			}
		col++;
		}
		fTemp_matrix_nsd_x_nsd3*= fMaterial_Params[kMu];
		row=0;
		col=0;

		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int i=0;i<3;i++)
				{
					for(int n=0;n<3;n++)
					{
						for(int m=0;m<3;m++)
						{
							I_temp_11p_4(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_Inverse(I,i)*fCauchy_stress_tensor_current_n_IP(n,i)
							*fDeformation_Gradient_n(m,K);
						}
					}
				}
			row++;
			}
		col++;
		}
		I_temp_11p_4*= fMaterial_Params[kMu];

		I_temp_11p_4+= fTemp_matrix_nsd_x_nsd2;
		I_temp_11p_4+= fTemp_matrix_nsd_x_nsd3;
}

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_5()
{

    int row=0;
    int col=0;
    I_temp_11p_5=0.0;
    fTemp_matrix_nsd_x_nsd = 0.0;


		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int i=0;i<3;i++)
				{
					for(int p=0;p<3;p++)
					{
						for(int n=0;n<3;n++)
						{
							for(int m=0;m<3;m++)
							{
								fTemp_matrix_nsd_x_nsd(row,col)+=fdFYdCauchy_Stress(m,n)*fCauchy_stress_tensor_current_n_IP(p,i)
								*fDeformation_Gradient_n(p,K)*fDeformation_Gradient_Inverse(I,i)*fCauchy_stress_tensor_current_n_IP(m,n);
							}
						}
					}
				}
			row++;
			}
		col++;
		}
		fTemp_matrix_nsd_x_nsd*= 2;

		row=0;
		col=0;

		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int i=0;i<3;i++)
				{
					for(int p=0;p<3;p++)
					{
						for(int n=0;n<3;n++)
						{
							for(int m=0;m<3;m++)
							{
								I_temp_11p_5(row,col)+=fdFYdCauchy_Stress(m,n)*fCauchy_stress_tensor_current_n_IP(p,i)
								*fDeformation_Gradient_n(p,K)*fDeformation_Gradient_Inverse(I,i)*fIdentity_matrix(m,n);
							}
						}
					}
				}
			row++;
			}
		col++;
		}
		I_temp_11p_5*= fMaterial_Params[kMu];
		I_temp_11p_5*= 2;
		I_temp_11p_5+= fTemp_matrix_nsd_x_nsd;
}

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_6()
{

    int row=0;
    int col=0;
    I_temp_11p_6=0.0;
    fTemp_matrix_nsd_x_nsd = 0.0;
    fTemp_matrix_nsd_x_nsd2 = 0.0;
    fTemp_matrix_nsd_x_nsd3 = 0.0;


		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int p=0;p<3;p++)
				{
					for(int n=0;n<3;n++)
					{
						for(int m=0;m<3;m++)
						{
							fTemp_matrix_nsd_x_nsd(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(m,K)*fDeformation_Gradient_Inverse(I,p)
							*fCauchy_stress_tensor_current_n_IP(p,n);
						}
					}
				}
			row++;
			}
		col++;
		}

		row=0;
		col=0;

		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int p=0;p<3;p++)
				{
					for(int n=0;n<3;n++)
					{
						for(int m=0;m<3;m++)
						{
							fTemp_matrix_nsd_x_nsd2(row,col)+=fdFYdCauchy_Stress(m,n)*fCauchy_stress_tensor_current_n_IP(m,p)
							*fDeformation_Gradient_n(n,K)*fDeformation_Gradient_Inverse(I,p);
						}
					}
				}
			row++;
			}
		col++;
		}

		row=0;
		col=0;

		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int n=0;n<3;n++)
				{
					for(int m=0;m<3;m++)
					{
						fTemp_matrix_nsd_x_nsd3(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(m,K)*fDeformation_Gradient_Inverse(I,n);
					}
				}
			row++;
			}
		col++;
		}
		fTemp_matrix_nsd_x_nsd3*= fMaterial_Params[kMu];
		row=0;
		col=0;

		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int n=0;n<3;n++)
				{
					for(int m=0;m<3;m++)
					{
						I_temp_11p_6(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(n,K)*fDeformation_Gradient_Inverse(I,m);
					}
				}
			row++;
			}
		col++;
		}
		I_temp_11p_6*= fMaterial_Params[kMu];

		I_temp_11p_6+= fTemp_matrix_nsd_x_nsd;
		I_temp_11p_6+= fTemp_matrix_nsd_x_nsd2;
		I_temp_11p_6+= fTemp_matrix_nsd_x_nsd3;
}

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_7()
{

    int row=0;
    int col=0;
    I_temp_11p_7=0.0;
    fTemp_matrix_nsd_x_nsd = 0.0;


    for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int p=0;p<3;p++)
				{
					for(int n=0;n<3;n++)
					{
						for(int m=0;m<3;m++)
						{
							fTemp_matrix_nsd_x_nsd(row,col)+=fdFYdCauchy_Stress(m,n)
							*fDeformation_Gradient_n(p,K)*fDeformation_Gradient_Inverse(I,p)*fCauchy_stress_tensor_current_n_IP(m,n);
						}
					}
				}
			row++;
			}
		col++;
		}
		fTemp_matrix_nsd_x_nsd*= 2;

		row=0;
		col=0;

		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int p=0;p<3;p++)
				{
					for(int n=0;n<3;n++)
					{
						for(int m=0;m<3;m++)
						{
							I_temp_11p_7(row,col)+=fdFYdCauchy_Stress(m,n)
							*fDeformation_Gradient_n(p,K)*fDeformation_Gradient_Inverse(I,p)*fIdentity_matrix(m,n);
						}
					}
				}
			row++;
			}
		col++;
		}
		I_temp_11p_7*= fMaterial_Params[kMu];
		I_temp_11p_7*= 2;
		I_temp_11p_7+= fTemp_matrix_nsd_x_nsd;
}


void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_8()
{

    int row=0;
    int col=0;
    I_temp_11p_8=0.0;
    fTemp_matrix_nsd_x_nsd = 0.0;
    fTemp_matrix_nsd_x_nsd2 = 0.0;
    fTemp_matrix_nsd_x_nsd3 = 0.0;


    for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int p=0;p<3;p++)
				{
					for(int n=0;n<3;n++)
					{
						for(int m=0;m<3;m++)
						{
							fTemp_matrix_nsd_x_nsd(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(n,K)*fDeformation_Gradient_Inverse(I,p)
							*fCauchy_stress_tensor_current_n_IP(p,m);
						}
					}
				}
			row++;
			}
		col++;
		}

		row=0;
		col=0;

		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int p=0;p<3;p++)
				{
					for(int n=0;n<3;n++)
					{
						for(int m=0;m<3;m++)
						{
							fTemp_matrix_nsd_x_nsd2(row,col)+=fdFYdCauchy_Stress(m,n)*fCauchy_stress_tensor_current_n_IP(n,p)
							*fDeformation_Gradient_n(m,K)*fDeformation_Gradient_Inverse(I,p);
						}
					}
				}
			row++;
			}
		col++;
		}

		row=0;
		col=0;

		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int n=0;n<3;n++)
				{
					for(int m=0;m<3;m++)
					{
						fTemp_matrix_nsd_x_nsd3(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(n,K)*fDeformation_Gradient_Inverse(I,m);
					}
				}
			row++;
			}
		col++;
		}
		fTemp_matrix_nsd_x_nsd3*= fMaterial_Params[kMu];
		row=0;
		col=0;

		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int n=0;n<3;n++)
				{
					for(int m=0;m<3;m++)
					{
						I_temp_11p_8(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(m,K)*fDeformation_Gradient_Inverse(I,n);
					}
				}
			row++;
			}
		col++;
		}
		I_temp_11p_8*= fMaterial_Params[kMu];
		I_temp_11p_8+= fTemp_matrix_nsd_x_nsd;
		I_temp_11p_8+= fTemp_matrix_nsd_x_nsd2;
		I_temp_11p_8+= fTemp_matrix_nsd_x_nsd3;
}


void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_9()
{

    int row=0;
    int col=0;
    I_temp_11p_9=0.0;
    fTemp_matrix_nsd_x_nsd = 0.0;
    fTemp_matrix_nsd_x_nsd2 = 0.0;
    fTemp_matrix_nsd_x_nsd3 = 0.0;


    for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int i=0;i<3;i++)
				{
					for(int p=0;p<3;p++)
					{
						for(int n=0;n<3;n++)
						{
							for(int m=0;m<3;m++)
							{
								for(int a=0;a<3;a++)
								{
									for(int b=0;b<3;b++)
									{
										fTemp_matrix_nsd_x_nsd(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(I,i)
										*fdev_Cauchy_stress_tensor_current_n_IP(a,b)*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(m,p)
										*fCauchy_stress_tensor_current_n_IP(p,n);
									}
								}
							}
						}
					}
				}
			row++;
			}
		col++;
		}

		row=0;
		col=0;

		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int i=0;i<3;i++)
				{
					for(int p=0;p<3;p++)
					{
						for(int n=0;n<3;n++)
						{
							for(int m=0;m<3;m++)
							{
								for(int a=0;a<3;a++)
								{
									for(int b=0;b<3;b++)
									{
										fTemp_matrix_nsd_x_nsd2(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(I,i)
										*fdev_Cauchy_stress_tensor_current_n_IP(a,b)*dev_Cauchy_stress_tr(a,b)*fCauchy_stress_tensor_current_n_IP(m,p)
										*dev_Cauchy_stress_tr(n,p);
									}
								}
							}
						}
					}
				}
			row++;
			}
		col++;
		}

		row=0;
		col=0;

		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int i=0;i<3;i++)
				{
					for(int n=0;n<3;n++)
					{
						for(int m=0;m<3;m++)
						{
							for(int a=0;a<3;a++)
							{
								for(int b=0;b<3;b++)
								{
									fTemp_matrix_nsd_x_nsd3(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(I,i)
									*fdev_Cauchy_stress_tensor_current_n_IP(a,b)*dev_Cauchy_stress_tr(a,b)
									*dev_Cauchy_stress_tr(m,n);
								}
							}
						}
					}
				}
			row++;
			}
		col++;
		}

		fTemp_matrix_nsd_x_nsd3*= fMaterial_Params[kMu];
		row=0;
		col=0;

		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int i=0;i<3;i++)
				{
					for(int n=0;n<3;n++)
					{
						for(int m=0;m<3;m++)
						{
							for(int a=0;a<3;a++)
							{
								for(int b=0;b<3;b++)
								{
									I_temp_11p_9(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(I,i)
									*fdev_Cauchy_stress_tensor_current_n_IP(a,b)*dev_Cauchy_stress_tr(a,b)
									*dev_Cauchy_stress_tr(n,m);
								}
							}
						}
					}
				}
			row++;
			}
		col++;
		}
		I_temp_11p_9*= fMaterial_Params[kMu];
		I_temp_11p_9+= fTemp_matrix_nsd_x_nsd;
		I_temp_11p_9+= fTemp_matrix_nsd_x_nsd2;
		I_temp_11p_9+= fTemp_matrix_nsd_x_nsd3;
}

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_10()
{

    int row=0;
    int col=0;
    I_temp_11p_10=0.0;
    fTemp_matrix_nsd_x_nsd = 0.0;
    fTemp_matrix_nsd_x_nsd2 = 0.0;
    fTemp_matrix_nsd_x_nsd3 = 0.0;


    for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int i=0;i<3;i++)
				{
					for(int p=0;p<3;p++)
					{
						for(int n=0;n<3;n++)
						{
							for(int m=0;m<3;m++)
							{
								for(int a=0;a<3;a++)
								{
									for(int b=0;b<3;b++)
									{
										fTemp_matrix_nsd_x_nsd(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(a,K)*fDeformation_Gradient_Inverse(I,i)
										*fCauchy_stress_tensor_current_n_IP(i,b)*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(m,p)
										*fCauchy_stress_tensor_current_n_IP(p,n);
									}
								}
							}
						}
					}
				}
			row++;
			}
		col++;
		}

		row=0;
		col=0;

		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int i=0;i<3;i++)
				{
					for(int p=0;p<3;p++)
					{
						for(int n=0;n<3;n++)
						{
							for(int m=0;m<3;m++)
							{
								for(int a=0;a<3;a++)
								{
									for(int b=0;b<3;b++)
									{
										fTemp_matrix_nsd_x_nsd2(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(a,K)*fDeformation_Gradient_Inverse(I,i)
										*fCauchy_stress_tensor_current_n_IP(i,b)*dev_Cauchy_stress_tr(a,b)*fCauchy_stress_tensor_current_n_IP(m,p)
										*dev_Cauchy_stress_tr(n,p);
									}
								}
							}
						}
					}
				}
			row++;
			}
		col++;
		}

		row=0;
		col=0;

		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int i=0;i<3;i++)
				{
					for(int n=0;n<3;n++)
					{
						for(int m=0;m<3;m++)
						{
							for(int a=0;a<3;a++)
							{
								for(int b=0;b<3;b++)
								{
									fTemp_matrix_nsd_x_nsd3(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(a,K)*fDeformation_Gradient_Inverse(I,i)
									*fCauchy_stress_tensor_current_n_IP(i,b)*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(m,n);
								}
							}
						}
					}
				}
			row++;
			}
		col++;
		}
		fTemp_matrix_nsd_x_nsd3*= fMaterial_Params[kMu];
		row=0;
		col=0;

		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int i=0;i<3;i++)
				{
					for(int n=0;n<3;n++)
					{
						for(int m=0;m<3;m++)
						{
							for(int a=0;a<3;a++)
							{
								for(int b=0;b<3;b++)
								{
									I_temp_11p_10(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(a,K)*fDeformation_Gradient_Inverse(I,i)
									*fCauchy_stress_tensor_current_n_IP(i,b)*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(n,m);
								}
							}
						}
					}
				}
			row++;
			}
		col++;
		}
		I_temp_11p_10*= fMaterial_Params[kMu];
		I_temp_11p_10+= fTemp_matrix_nsd_x_nsd;
		I_temp_11p_10+= fTemp_matrix_nsd_x_nsd2;
		I_temp_11p_10+= fTemp_matrix_nsd_x_nsd3;
}

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_11()
{

    int row=0;
    int col=0;
    I_temp_11p_11=0.0;
    fTemp_matrix_nsd_x_nsd = 0.0;
    fTemp_matrix_nsd_x_nsd2 = 0.0;
    fTemp_matrix_nsd_x_nsd3 = 0.0;


    for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int i=0;i<3;i++)
				{
					for(int p=0;p<3;p++)
					{
						for(int n=0;n<3;n++)
						{
							for(int m=0;m<3;m++)
							{
								for(int a=0;a<3;a++)
								{
									for(int b=0;b<3;b++)
									{
										fTemp_matrix_nsd_x_nsd(row,col)+=fdFYdCauchy_Stress(m,n)*fCauchy_stress_tensor_current_n_IP(a,i)
										*fDeformation_Gradient_n(b,K)*fDeformation_Gradient_Inverse(I,i)
										*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(m,p)*fCauchy_stress_tensor_current_n_IP(p,n);
									}
								}
							}
						}
					}
				}
			row++;
			}
		col++;
		}

		row=0;
		col=0;

		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int i=0;i<3;i++)
				{
					for(int p=0;p<3;p++)
					{
						for(int n=0;n<3;n++)
						{
							for(int m=0;m<3;m++)
							{
								for(int a=0;a<3;a++)
								{
									for(int b=0;b<3;b++)
									{
										fTemp_matrix_nsd_x_nsd2(row,col)+=fdFYdCauchy_Stress(m,n)*fCauchy_stress_tensor_current_n_IP(a,i)
										*fDeformation_Gradient_n(b,K)*fDeformation_Gradient_Inverse(I,i)*dev_Cauchy_stress_tr(a,b)
										*fCauchy_stress_tensor_current_n_IP(m,p)*dev_Cauchy_stress_tr(n,p);
									}
								}
							}
						}
					}
				}
			row++;
			}
		col++;
		}

		row=0;
		col=0;

		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int i=0;i<3;i++)
				{
					for(int n=0;n<3;n++)
					{
						for(int m=0;m<3;m++)
						{
							for(int a=0;a<3;a++)
							{
								for(int b=0;b<3;b++)
								{
									fTemp_matrix_nsd_x_nsd3(row,col)+=fdFYdCauchy_Stress(m,n)*fCauchy_stress_tensor_current_n_IP(a,i)
									*fDeformation_Gradient_n(b,K)*fDeformation_Gradient_Inverse(I,i)*dev_Cauchy_stress_tr(a,b)
									*dev_Cauchy_stress_tr(m,n);
								}
							}
						}
					}
				}
			row++;
			}
		col++;
		}
		fTemp_matrix_nsd_x_nsd3*= fMaterial_Params[kMu];

		row=0;
		col=0;

		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int i=0;i<3;i++)
				{
					for(int n=0;n<3;n++)
					{
						for(int m=0;m<3;m++)
						{
							for(int a=0;a<3;a++)
							{
								for(int b=0;b<3;b++)
								{
									I_temp_11p_11(row,col)+=fdFYdCauchy_Stress(m,n)*fCauchy_stress_tensor_current_n_IP(a,i)
									*fDeformation_Gradient_n(b,K)*fDeformation_Gradient_Inverse(I,i)*dev_Cauchy_stress_tr(a,b)
									*dev_Cauchy_stress_tr(n,m);
								}
							}
						}
					}
				}
			row++;
			}
		col++;
		}
		I_temp_11p_11*= fMaterial_Params[kMu];

		I_temp_11p_11+= fTemp_matrix_nsd_x_nsd;
		I_temp_11p_11+= fTemp_matrix_nsd_x_nsd2;
		I_temp_11p_11+= fTemp_matrix_nsd_x_nsd3;
}

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_12()
{

    int row=0;
    int col=0;
    I_temp_11p_12=0.0;
    fTemp_matrix_nsd_x_nsd = 0.0;
    fTemp_matrix_nsd_x_nsd2 = 0.0;
    fTemp_matrix_nsd_x_nsd3 = 0.0;


    for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int p=0;p<3;p++)
				{
					for(int n=0;n<3;n++)
					{
						for(int m=0;m<3;m++)
						{
							for(int a=0;a<3;a++)
							{
								for(int b=0;b<3;b++)
								{
									fTemp_matrix_nsd_x_nsd(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(a,K)
									*fDeformation_Gradient_Inverse(I,b)*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(m,p)
									*fCauchy_stress_tensor_current_n_IP(p,n);
								}
							}
						}
					}
				}
			row++;
			}
		col++;
		}

		row=0;
		col=0;

		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int p=0;p<3;p++)
				{
					for(int n=0;n<3;n++)
					{
						for(int m=0;m<3;m++)
						{
							for(int a=0;a<3;a++)
							{
								for(int b=0;b<3;b++)
								{
									fTemp_matrix_nsd_x_nsd2(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(a,K)
									*fDeformation_Gradient_Inverse(I,b)*dev_Cauchy_stress_tr(a,b)
									*fCauchy_stress_tensor_current_n_IP(m,p)*dev_Cauchy_stress_tr(n,p);
								}
							}
						}
					}
				}
			row++;
			}
		col++;
		}

		row=0;
		col=0;

		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int n=0;n<3;n++)
				{
					for(int m=0;m<3;m++)
					{
						for(int a=0;a<3;a++)
						{
							for(int b=0;b<3;b++)
							{
								fTemp_matrix_nsd_x_nsd3(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(a,K)
								*fDeformation_Gradient_Inverse(I,b)*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(m,n);
							}
						}
					}
				}
			row++;
			}
		col++;
		}
		fTemp_matrix_nsd_x_nsd3*= fMaterial_Params[kMu];

		row=0;
		col=0;

		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int n=0;n<3;n++)
				{
					for(int m=0;m<3;m++)
					{
						for(int a=0;a<3;a++)
						{
							for(int b=0;b<3;b++)
							{
								I_temp_11p_12(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(a,K)
								*fDeformation_Gradient_Inverse(I,b)*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(n,m);
							}
						}
					}
				}
			row++;
			}
		col++;
		}
		I_temp_11p_12*= fMaterial_Params[kMu];

		I_temp_11p_12+= fTemp_matrix_nsd_x_nsd;
		I_temp_11p_12+= fTemp_matrix_nsd_x_nsd2;
		I_temp_11p_12+= fTemp_matrix_nsd_x_nsd3;
}

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_13()
{

    int row=0;
    int col=0;
    I_temp_11p_13=0.0;
    fTemp_matrix_nsd_x_nsd = 0.0;
    fTemp_matrix_nsd_x_nsd2 = 0.0;
    fTemp_matrix_nsd_x_nsd3 = 0.0;


    for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int p=0;p<3;p++)
				{
					for(int n=0;n<3;n++)
					{
						for(int m=0;m<3;m++)
						{
							for(int a=0;a<3;a++)
							{
								for(int b=0;b<3;b++)
								{
									fTemp_matrix_nsd_x_nsd(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(b,K)
									*fDeformation_Gradient_Inverse(I,a)*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(m,p)
									*fCauchy_stress_tensor_current_n_IP(p,n);
								}
							}
						}
					}
				}
			row++;
			}
		col++;
		}

		row=0;
		col=0;

		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int p=0;p<3;p++)
				{
					for(int n=0;n<3;n++)
					{
						for(int m=0;m<3;m++)
						{
							for(int a=0;a<3;a++)
							{
								for(int b=0;b<3;b++)
								{
									fTemp_matrix_nsd_x_nsd2(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(b,K)
									*fDeformation_Gradient_Inverse(I,a)*dev_Cauchy_stress_tr(a,b)
									*fCauchy_stress_tensor_current_n_IP(m,p)*dev_Cauchy_stress_tr(n,p);
								}
							}
						}
					}
				}
			row++;
			}
		col++;
		}

		row=0;
		col=0;

		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int n=0;n<3;n++)
				{
					for(int m=0;m<3;m++)
					{
						for(int a=0;a<3;a++)
						{
							for(int b=0;b<3;b++)
							{
								fTemp_matrix_nsd_x_nsd3(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(b,K)
								*fDeformation_Gradient_Inverse(I,a)*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(m,n);
							}
						}
					}
				}
			row++;
			}
		col++;
		}
		fTemp_matrix_nsd_x_nsd3*= fMaterial_Params[kMu];
		row=0;
		col=0;

		for(int K=0;K<3;K++)
		{
			row=0;
			for(int I=0;I<3;I++)
			{
				for(int n=0;n<3;n++)
				{
					for(int m=0;m<3;m++)
					{
						for(int a=0;a<3;a++)
						{
							for(int b=0;b<3;b++)
							{
								I_temp_11p_13(row,col)+=fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(b,K)
								*fDeformation_Gradient_Inverse(I,a)*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(n,m);
							}
						}
					}
				}
			row++;
			}
		col++;
		}
		I_temp_11p_13*= fMaterial_Params[kMu];
		I_temp_11p_13+= fTemp_matrix_nsd_x_nsd;
		I_temp_11p_13+= fTemp_matrix_nsd_x_nsd2;
		I_temp_11p_13+= fTemp_matrix_nsd_x_nsd3;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_1_1()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_1_1=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int i=0;i<3;i++)
								{
									for(int p=0;p<3;p++)
									{
										for(int n=0;n<3;n++)
										{
											for(int m=0;m<3;m++)
											{
												I_temp_11p_test_2_1_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
												*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(I,i)
												*fdev_Cauchy_stress_tensor_current_n_IP(m,p)*fCauchy_stress_tensor_current_n_IP(p,n)*fDeformation_Gradient_Inverse(K,l);
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
			col++;
		}
	}
}


void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_1_2()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_1_2=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int i=0;i<3;i++)
								{
									for(int p=0;p<3;p++)
									{
										for(int n=0;n<3;n++)
										{
											for(int m=0;m<3;m++)
											{
												I_temp_11p_test_2_1_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
												*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(i,K)
												*fDeformation_Gradient_Inverse(I,i)*fCauchy_stress_tensor_current_n_IP(m,p)*fdev_Cauchy_stress_tensor_current_n_IP(n,p)*fDeformation_Gradient_Inverse(K,l);
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
			col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_1_3()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_1_3=0.0;

	    for (int I=0;I<3;I++)
		{
			for(int l=0;l<3;l++)
			{
				row=0;
				for(int L=0;L<3;L++)
				{
					for(int k=0;k<3;k++)
					{
						for(int K=0;K<3;K++)
						{
							for(int j=0;j<3;j++)
							{
								for(int i=0;i<3;i++)
								{
									for(int q=0;q<3;q++)
									{
										for(int n=0;n<3;n++)
										{
											for(int m=0;m<3;m++)
											{
													I_temp_11p_test_2_1_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
													*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(i,K)
													*fDeformation_Gradient_Inverse(I,i)*fdev_Cauchy_stress_tensor_current_n_IP(m,n)*fDeformation_Gradient_Inverse(K,l);
											}
										}
									}
								}
							}
						}
					row++;
					}
				}
			col++;
		}
	}
	    I_temp_11p_test_2_1_3*= fMaterial_Params[kMu];
}

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_1_4()
{

	int row=0;
	int col=0;
	I_temp_11p_test_2_1_4=0.0;

	for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							for(int i=0;i<3;i++)
							{
								for(int q=0;q<3;q++)
								{
									for(int n=0;n<3;n++)
									{
										for(int m=0;m<3;m++)
										{
												I_temp_11p_test_2_1_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
												*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(i,K)
												*fDeformation_Gradient_Inverse(I,i)*fdev_Cauchy_stress_tensor_current_n_IP(n,m)*fDeformation_Gradient_Inverse(K,l);
										}
									}
								}
							}
						}
					}
				row++;
				}
			}
		col++;
		}
	}
	I_temp_11p_test_2_1_4*= fMaterial_Params[kMu];
}





void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_2_1()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_2_1=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int i=0;i<3;i++)
								{
									for(int p=0;p<3;p++)
									{
										for(int n=0;n<3;n++)
										{
											for(int m=0;m<3;m++)
											{
												I_temp_11p_test_2_2_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
												*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(m,K)*fDeformation_Gradient_Inverse(I,i)
												*fCauchy_stress_tensor_current_n_IP(i,p)*fCauchy_stress_tensor_current_n_IP(p,n)*fDeformation_Gradient_Inverse(K,l);
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
			col++;
		}
	}
}


void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_2_2()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_2_2=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int i=0;i<3;i++)
								{
									for(int p=0;p<3;p++)
									{
										for(int n=0;n<3;n++)
										{
											for(int m=0;m<3;m++)
											{
												I_temp_11p_test_2_2_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
												*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fCauchy_stress_tensor_current_n_IP(m,p)*fDeformation_Gradient_n(n,K)
												*fDeformation_Gradient_Inverse(I,i)*fCauchy_stress_tensor_current_n_IP(i,p)*fDeformation_Gradient_Inverse(K,l);
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
			col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_2_3()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_2_3=0.0;

	    for (int I=0;I<3;I++)
		{
			for(int l=0;l<3;l++)
			{
				row=0;
				for(int L=0;L<3;L++)
				{
					for(int k=0;k<3;k++)
					{
						for(int K=0;K<3;K++)
						{
							for(int j=0;j<3;j++)
							{
								for(int i=0;i<3;i++)
								{
									for(int q=0;q<3;q++)
									{
										for(int n=0;n<3;n++)
										{
											for(int m=0;m<3;m++)
											{
													I_temp_11p_test_2_2_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
													*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(m,K)
													*fDeformation_Gradient_Inverse(I,i)*fCauchy_stress_tensor_current_n_IP(i,n)*fDeformation_Gradient_Inverse(K,l);
											}
										}
									}
								}
							}
						}
					row++;
					}
				}
			col++;
		}
	}
	    I_temp_11p_test_2_2_3*= fMaterial_Params[kMu];
}

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_2_4()
{

	int row=0;
	int col=0;
	I_temp_11p_test_2_2_4=0.0;

	for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							for(int i=0;i<3;i++)
							{
								for(int q=0;q<3;q++)
								{
									for(int n=0;n<3;n++)
									{
										for(int m=0;m<3;m++)
										{
												I_temp_11p_test_2_2_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
												*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(n,K)
												*fDeformation_Gradient_Inverse(I,i)*fCauchy_stress_tensor_current_n_IP(i,m)*fDeformation_Gradient_Inverse(K,l);
										}
									}
								}
							}
						}
					}
				row++;
				}
			}
		col++;
		}
	}
	I_temp_11p_test_2_2_4*= fMaterial_Params[kMu];
}


void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_3_1()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_3_1=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int i=0;i<3;i++)
								{
									for(int p=0;p<3;p++)
									{
										for(int n=0;n<3;n++)
										{
											for(int m=0;m<3;m++)
											{
												I_temp_11p_test_2_3_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
												*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(I,p)
												*fCauchy_stress_tensor_current_n_IP(p,i)*fCauchy_stress_tensor_current_n_IP(m,n)*fDeformation_Gradient_Inverse(K,l);
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
			col++;
		}
	}
    I_temp_11p_test_2_3_1*= 2.0;
}


void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_3_2()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_3_2=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int i=0;i<3;i++)
								{
									for(int p=0;p<3;p++)
									{
										for(int n=0;n<3;n++)
										{
											for(int m=0;m<3;m++)
											{
												I_temp_11p_test_2_3_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
												*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(I,p)
												*fCauchy_stress_tensor_current_n_IP(p,i)*fIdentity_matrix(m,n)*fDeformation_Gradient_Inverse(K,l);
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
			col++;
		}
	}
    I_temp_11p_test_2_3_2*= 2.0;
    I_temp_11p_test_2_3_2*= fMaterial_Params[kMu];
}



void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_4_1()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_4_1=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int i=0;i<3;i++)
					{
						for(int p=0;p<3;p++)
						{
							for(int n=0;n<3;n++)
							{
								for(int m=0;m<3;m++)
								{
									for(int j=0;j<3;j++)
									{
										for(int q=0;q<3;q++)
										{
											for(int K=0;K<3;K++)
											{
												I_temp_11p_test_2_4_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
												*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fCauchy_stress_tensor_current_n_IP(m,i)*fDeformation_Gradient_n(p,K)
												*fDeformation_Gradient_Inverse(K,l)*fDeformation_Gradient_Inverse(I,i)
												*fCauchy_stress_tensor_current_n_IP(p,n);
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
			col++;
		}
	}
}


void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_4_2()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_4_2=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int i=0;i<3;i++)
					{
						for(int p=0;p<3;p++)
						{
							for(int n=0;n<3;n++)
							{
								for(int m=0;m<3;m++)
								{
									for(int j=0;j<3;j++)
									{
										for(int q=0;q<3;q++)
										{
											for(int K=0;K<3;K++)
											{
												I_temp_11p_test_2_4_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
												*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fCauchy_stress_tensor_current_n_IP(m,p)
												*fCauchy_stress_tensor_current_n_IP(n,i)*fDeformation_Gradient_n(p,K)*fDeformation_Gradient_Inverse(K,l)
												*fDeformation_Gradient_Inverse(I,i);
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
			col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_4_3()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_4_3=0.0;

	    for (int I=0;I<3;I++)
		{
			for(int l=0;l<3;l++)
			{
				row=0;
				for(int L=0;L<3;L++)
				{
					for(int k=0;k<3;k++)
					{
						for(int i=0;i<3;i++)
						{
							for(int K=0;K<3;K++)
							{
								for(int n=0;n<3;n++)
								{
									for(int m=0;m<3;m++)
									{
										for(int j=0;j<3;j++)
										{
											for(int q=0;q<3;q++)
											{
													I_temp_11p_test_2_4_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
													*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fCauchy_stress_tensor_current_n_IP(m,i)
													*fDeformation_Gradient_n(n,K)*fDeformation_Gradient_Inverse(K,l)
													*fDeformation_Gradient_Inverse(I,i);
											}
										}
									}
								}
							}
						}
					row++;
					}
				}
			col++;
			}
		}
	    I_temp_11p_test_2_4_3*= fMaterial_Params[kMu];
}

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_4_4()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_4_4=0.0;

	    for (int I=0;I<3;I++)
		{
			for(int l=0;l<3;l++)
			{
				row=0;
				for(int L=0;L<3;L++)
				{
					for(int k=0;k<3;k++)
					{
						for(int i=0;i<3;i++)
						{
							for(int K=0;K<3;K++)
							{
								for(int n=0;n<3;n++)
								{
									for(int m=0;m<3;m++)
									{
										for(int j=0;j<3;j++)
										{
											for(int q=0;q<3;q++)
											{
													I_temp_11p_test_2_4_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
													*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fCauchy_stress_tensor_current_n_IP(n,i)
													*fDeformation_Gradient_n(m,K)*fDeformation_Gradient_Inverse(K,l)
													*fDeformation_Gradient_Inverse(I,i);
											}
										}
									}
								}
							}
						}
					row++;
					}
				}
			col++;
			}
		}
	    I_temp_11p_test_2_4_4*= fMaterial_Params[kMu];
}



void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_5_1()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_5_1=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int i=0;i<3;i++)
								{
									for(int p=0;p<3;p++)
									{
										for(int n=0;n<3;n++)
										{
											for(int m=0;m<3;m++)
											{
												I_temp_11p_test_2_5_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
												*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fCauchy_stress_tensor_current_n_IP(p,i)
												*fDeformation_Gradient_n(p,K)*fDeformation_Gradient_Inverse(I,i)
												*fCauchy_stress_tensor_current_n_IP(m,n)*fDeformation_Gradient_Inverse(K,l);
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
			col++;
		}
	}
    I_temp_11p_test_2_5_1*= 2.0;
}




void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_5_2()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_5_2=0.0;

	    for (int I=0;I<3;I++)
		{
			for(int l=0;l<3;l++)
			{
				row=0;
				for(int L=0;L<3;L++)
				{
					for(int k=0;k<3;k++)
					{
						for(int K=0;K<3;K++)
						{
							for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									for(int i=0;i<3;i++)
									{
										for(int q=0;q<3;q++)
										{
											for(int n=0;n<3;n++)
											{
												for(int m=0;m<3;m++)
												{
														I_temp_11p_test_2_5_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
														*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fCauchy_stress_tensor_current_n_IP(p,i)*fDeformation_Gradient_n(p,K)
														*fDeformation_Gradient_Inverse(I,i)*fIdentity_matrix(m,n)*fDeformation_Gradient_Inverse(K,l);
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
			col++;
		}
	}
	    I_temp_11p_test_2_5_2*= 2.0;
	    I_temp_11p_test_2_5_2*= fMaterial_Params[kMu];
}



void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_6_1()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_6_1=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int p=0;p<3;p++)
								{
									for(int n=0;n<3;n++)
									{
										for(int m=0;m<3;m++)
										{
											I_temp_11p_test_2_6_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
											*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(m,K)*fDeformation_Gradient_Inverse(I,p)
											*fCauchy_stress_tensor_current_n_IP(p,n)*fDeformation_Gradient_Inverse(K,l);
										}
									}
								}
							}
						}
					}
				row++;
				}
			}
			col++;
		}
	}
}


void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_6_2()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_6_2=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int p=0;p<3;p++)
								{
									for(int n=0;n<3;n++)
									{
										for(int m=0;m<3;m++)
										{
											I_temp_11p_test_2_6_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
											*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fCauchy_stress_tensor_current_n_IP(m,p)*fDeformation_Gradient_n(n,K)
											*fDeformation_Gradient_Inverse(I,p)*fDeformation_Gradient_Inverse(K,l);
										}
									}
								}
							}
						}
					}
				row++;
				}
			}
			col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_6_3()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_6_3=0.0;

	    for (int I=0;I<3;I++)
		{
			for(int l=0;l<3;l++)
			{
				row=0;
				for(int L=0;L<3;L++)
				{
					for(int k=0;k<3;k++)
					{
						for(int K=0;K<3;K++)
						{
							for(int j=0;j<3;j++)
							{
								for(int q=0;q<3;q++)
								{
									for(int n=0;n<3;n++)
									{
										for(int m=0;m<3;m++)
										{
												I_temp_11p_test_2_6_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
												*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(m,K)
												*fDeformation_Gradient_Inverse(I,n)*fDeformation_Gradient_Inverse(K,l);
										}
									}
								}
							}
						}
					row++;
					}
				}
			col++;
		}
	}
	    I_temp_11p_test_2_6_3*= fMaterial_Params[kMu];
}

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_6_4()
{

	int row=0;
	int col=0;
	I_temp_11p_test_2_6_4=0.0;

	for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							for(int q=0;q<3;q++)
							{
								for(int n=0;n<3;n++)
								{
									for(int m=0;m<3;m++)
									{
											I_temp_11p_test_2_6_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
											*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(n,K)
											*fDeformation_Gradient_Inverse(I,m)*fDeformation_Gradient_Inverse(K,l);
									}
								}
							}
						}
					}
				row++;
				}
			}
		col++;
		}
	}
	I_temp_11p_test_2_6_4*= fMaterial_Params[kMu];
}




void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_7_1()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_7_1=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int p=0;p<3;p++)
								{
									for(int n=0;n<3;n++)
									{
										for(int m=0;m<3;m++)
										{
											I_temp_11p_test_2_7_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
											*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)
											*fDeformation_Gradient_n(p,K)*fDeformation_Gradient_Inverse(I,p)
											*fCauchy_stress_tensor_current_n_IP(m,n)*fDeformation_Gradient_Inverse(K,l);
										}
									}
								}
							}
						}
					}
				row++;
				}
			}
			col++;
		}
	}
    I_temp_11p_test_2_7_1*= 2.0;
}




void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_7_2()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_7_2=0.0;

	    for (int I=0;I<3;I++)
		{
			for(int l=0;l<3;l++)
			{
				row=0;
				for(int L=0;L<3;L++)
				{
					for(int k=0;k<3;k++)
					{
						for(int K=0;K<3;K++)
						{
							for(int q=0;q<3;q++)
							{
								for(int j=0;j<3;j++)
								{
									for(int p=0;p<3;p++)
									{
										for(int n=0;n<3;n++)
										{
											for(int m=0;m<3;m++)
											{
												I_temp_11p_test_2_7_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
												*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(p,K)
												*fDeformation_Gradient_Inverse(I,p)*fIdentity_matrix(m,n)*fDeformation_Gradient_Inverse(K,l);
											}
										}
									}
								}
							}
						}
					row++;
					}
				}
			col++;
		}
	}
	    I_temp_11p_test_2_7_2*= 2.0;
	    I_temp_11p_test_2_7_2*= fMaterial_Params[kMu];
}




void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_8_1()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_8_1=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int p=0;p<3;p++)
								{
									for(int n=0;n<3;n++)
									{
										for(int m=0;m<3;m++)
										{
											I_temp_11p_test_2_8_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
											*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(n,K)*fDeformation_Gradient_Inverse(I,p)
											*fCauchy_stress_tensor_current_n_IP(p,m)*fDeformation_Gradient_Inverse(K,l);
										}
									}
								}
							}
						}
					}
				row++;
				}
			}
			col++;
		}
	}
}


void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_8_2()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_8_2=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int p=0;p<3;p++)
								{
									for(int n=0;n<3;n++)
									{
										for(int m=0;m<3;m++)
										{
											I_temp_11p_test_2_8_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
											*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fCauchy_stress_tensor_current_n_IP(n,p)*fDeformation_Gradient_n(m,K)
											*fDeformation_Gradient_Inverse(I,p)*fDeformation_Gradient_Inverse(K,l);
										}
									}
								}
							}
						}
					}
				row++;
				}
			}
			col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_8_3()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_8_3=0.0;

	    for (int I=0;I<3;I++)
		{
			for(int l=0;l<3;l++)
			{
				row=0;
				for(int L=0;L<3;L++)
				{
					for(int k=0;k<3;k++)
					{
						for(int K=0;K<3;K++)
						{
							for(int j=0;j<3;j++)
							{
								for(int q=0;q<3;q++)
								{
									for(int n=0;n<3;n++)
									{
										for(int m=0;m<3;m++)
										{
												I_temp_11p_test_2_8_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
												*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(n,K)
												*fDeformation_Gradient_Inverse(I,m)*fDeformation_Gradient_Inverse(K,l);
										}
									}
								}
							}
						}
					row++;
					}
				}
			col++;
		}
	}
	    I_temp_11p_test_2_8_3*= fMaterial_Params[kMu];
}

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_8_4()
{

	int row=0;
	int col=0;
	I_temp_11p_test_2_8_4=0.0;

	for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							for(int q=0;q<3;q++)
							{
								for(int n=0;n<3;n++)
								{
									for(int m=0;m<3;m++)
									{
											I_temp_11p_test_2_8_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
											*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(m,K)
											*fDeformation_Gradient_Inverse(I,n)*fDeformation_Gradient_Inverse(K,l);
									}
								}
							}
						}
					}
				row++;
				}
			}
		col++;
		}
	}
	I_temp_11p_test_2_8_4*= fMaterial_Params[kMu];
}



void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_9_1()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_9_1=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int i=0;i<3;i++)
								{
									for(int b=0;b<3;b++)
									{
										for(int a=0;a<3;a++)
										{
											for(int p=0;p<3;p++)
											{
												for(int n=0;n<3;n++)
												{
													for(int m=0;m<3;m++)
													{
														I_temp_11p_test_2_9_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
														*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(I,i)
														*fdev_Cauchy_stress_tensor_current_n_IP(a,b)*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(m,p)
														*fCauchy_stress_tensor_current_n_IP(p,n)*fDeformation_Gradient_Inverse(K,l);
													}
												}
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
			col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_9_2()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_9_2=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int i=0;i<3;i++)
								{
									for(int b=0;b<3;b++)
									{
										for(int a=0;a<3;a++)
										{
											for(int p=0;p<3;p++)
											{
												for(int n=0;n<3;n++)
												{
													for(int m=0;m<3;m++)
													{
														I_temp_11p_test_2_9_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
														*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fCauchy_stress_tensor_current_n_IP(m,p)
														*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(I,i)
														*fdev_Cauchy_stress_tensor_current_n_IP(a,b)*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(n,p)
														*fDeformation_Gradient_Inverse(K,l);
													}
												}
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
			col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_9_3()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_9_3=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int i=0;i<3;i++)
								{
									for(int b=0;b<3;b++)
									{
										for(int a=0;a<3;a++)
										{
											for(int n=0;n<3;n++)
											{
												for(int m=0;m<3;m++)
												{
													I_temp_11p_test_2_9_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
													*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)
													*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(I,i)
													*fdev_Cauchy_stress_tensor_current_n_IP(a,b)*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(m,n)
													*fDeformation_Gradient_Inverse(K,l);
												}
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
			col++;
		}
	}
    I_temp_11p_test_2_9_3*= fMaterial_Params[kMu];
}
void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_9_4()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_9_4=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int i=0;i<3;i++)
								{
									for(int b=0;b<3;b++)
									{
										for(int a=0;a<3;a++)
										{
											for(int n=0;n<3;n++)
											{
												for(int m=0;m<3;m++)
												{
													I_temp_11p_test_2_9_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
													*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)
													*fDeformation_Gradient_n(i,K)*fDeformation_Gradient_Inverse(I,i)
													*fdev_Cauchy_stress_tensor_current_n_IP(a,b)*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(n,m)
													*fDeformation_Gradient_Inverse(K,l);
												}
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
			col++;
		}
	}
    I_temp_11p_test_2_9_4*= fMaterial_Params[kMu];
}




void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_10_1()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_10_1=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int i=0;i<3;i++)
								{
									for(int b=0;b<3;b++)
									{
										for(int a=0;a<3;a++)
										{
											for(int p=0;p<3;p++)
											{
												for(int n=0;n<3;n++)
												{
													for(int m=0;m<3;m++)
													{
														I_temp_11p_test_2_10_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
														*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(a,K)*fDeformation_Gradient_Inverse(I,i)
														*fCauchy_stress_tensor_current_n_IP(i,b)*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(m,p)
														*fCauchy_stress_tensor_current_n_IP(p,n)*fDeformation_Gradient_Inverse(K,l);
													}
												}
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
			col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_10_2()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_10_2=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int i=0;i<3;i++)
								{
									for(int b=0;b<3;b++)
									{
										for(int a=0;a<3;a++)
										{
											for(int p=0;p<3;p++)
											{
												for(int n=0;n<3;n++)
												{
													for(int m=0;m<3;m++)
													{
														I_temp_11p_test_2_10_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
														*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fCauchy_stress_tensor_current_n_IP(m,p)
														*fDeformation_Gradient_n(a,K)*fDeformation_Gradient_Inverse(I,i)
														*fCauchy_stress_tensor_current_n_IP(i,b)*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(n,p)
														*fDeformation_Gradient_Inverse(K,l);
													}
												}
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
			col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_10_3()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_10_3=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int i=0;i<3;i++)
								{
									for(int b=0;b<3;b++)
									{
										for(int a=0;a<3;a++)
										{
											for(int n=0;n<3;n++)
											{
												for(int m=0;m<3;m++)
												{
													I_temp_11p_test_2_10_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
													*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)
													*fDeformation_Gradient_n(a,K)*fDeformation_Gradient_Inverse(I,i)
													*fCauchy_stress_tensor_current_n_IP(i,b)*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(m,n)
													*fDeformation_Gradient_Inverse(K,l);
												}
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
			col++;
		}
	}
    I_temp_11p_test_2_10_3*= fMaterial_Params[kMu];
}
void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_10_4()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_10_4=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int i=0;i<3;i++)
								{
									for(int b=0;b<3;b++)
									{
										for(int a=0;a<3;a++)
										{
											for(int n=0;n<3;n++)
											{
												for(int m=0;m<3;m++)
												{
													I_temp_11p_test_2_10_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
													*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)
													*fDeformation_Gradient_n(a,K)*fDeformation_Gradient_Inverse(I,i)
													*fCauchy_stress_tensor_current_n_IP(i,b)*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(n,m)
													*fDeformation_Gradient_Inverse(K,l);
												}
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
			col++;
		}
	}
    I_temp_11p_test_2_10_4*= fMaterial_Params[kMu];
}




void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_11_1()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_11_1=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int i=0;i<3;i++)
								{
									for(int b=0;b<3;b++)
									{
										for(int a=0;a<3;a++)
										{
											for(int p=0;p<3;p++)
											{
												for(int n=0;n<3;n++)
												{
													for(int m=0;m<3;m++)
													{
														I_temp_11p_test_2_11_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
														*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fCauchy_stress_tensor_current_n_IP(a,i)
														*fDeformation_Gradient_n(b,K)*fDeformation_Gradient_Inverse(I,i)
														*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(m,p)
														*fCauchy_stress_tensor_current_n_IP(p,n)*fDeformation_Gradient_Inverse(K,l);
													}
												}
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
			col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_11_2()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_11_2=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int i=0;i<3;i++)
								{
									for(int b=0;b<3;b++)
									{
										for(int a=0;a<3;a++)
										{
											for(int p=0;p<3;p++)
											{
												for(int n=0;n<3;n++)
												{
													for(int m=0;m<3;m++)
													{
														I_temp_11p_test_2_11_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
														*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fCauchy_stress_tensor_current_n_IP(m,p)
														*fCauchy_stress_tensor_current_n_IP(a,i)*fDeformation_Gradient_n(b,K)*fDeformation_Gradient_Inverse(I,i)
														*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(n,p)
														*fDeformation_Gradient_Inverse(K,l);
													}
												}
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
			col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_11_3()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_11_3=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int i=0;i<3;i++)
								{
									for(int b=0;b<3;b++)
									{
										for(int a=0;a<3;a++)
										{
											for(int n=0;n<3;n++)
											{
												for(int m=0;m<3;m++)
												{
													I_temp_11p_test_2_11_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
													*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fCauchy_stress_tensor_current_n_IP(a,i)
													*fDeformation_Gradient_n(b,K)*fDeformation_Gradient_Inverse(I,i)
													*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(m,n)
													*fDeformation_Gradient_Inverse(K,l);
												}
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
			col++;
		}
	}
    I_temp_11p_test_2_11_3*= fMaterial_Params[kMu];
}
void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_11_4()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_11_4=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int i=0;i<3;i++)
								{
									for(int b=0;b<3;b++)
									{
										for(int a=0;a<3;a++)
										{
											for(int n=0;n<3;n++)
											{
												for(int m=0;m<3;m++)
												{
													I_temp_11p_test_2_11_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
													*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)
													*fCauchy_stress_tensor_current_n_IP(a,i)*fDeformation_Gradient_n(b,K)*fDeformation_Gradient_Inverse(I,i)
													*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(n,m)
													*fDeformation_Gradient_Inverse(K,l);
												}
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
			col++;
		}
	}
    I_temp_11p_test_2_11_4*= fMaterial_Params[kMu];
}


void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_12_1()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_12_1=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int b=0;b<3;b++)
								{
									for(int a=0;a<3;a++)
									{
										for(int p=0;p<3;p++)
										{
											for(int n=0;n<3;n++)
											{
												for(int m=0;m<3;m++)
												{
													I_temp_11p_test_2_12_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
													*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(a,K)*fDeformation_Gradient_Inverse(I,b)
													*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(m,p)
													*fCauchy_stress_tensor_current_n_IP(p,n)*fDeformation_Gradient_Inverse(K,l);
												}
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
			col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_12_2()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_12_2=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int b=0;b<3;b++)
								{
									for(int a=0;a<3;a++)
									{
										for(int p=0;p<3;p++)
										{
											for(int n=0;n<3;n++)
											{
												for(int m=0;m<3;m++)
												{
													I_temp_11p_test_2_12_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
													*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fCauchy_stress_tensor_current_n_IP(m,p)
													*fDeformation_Gradient_n(a,K)*fDeformation_Gradient_Inverse(I,b)
													*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(n,p)
													*fDeformation_Gradient_Inverse(K,l);
												}
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
			col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_12_3()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_12_3=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int b=0;b<3;b++)
								{
									for(int a=0;a<3;a++)
									{
										for(int n=0;n<3;n++)
										{
											for(int m=0;m<3;m++)
											{
												I_temp_11p_test_2_12_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
												*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)
												*fDeformation_Gradient_n(a,K)*fDeformation_Gradient_Inverse(I,b)
												*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(m,n)
												*fDeformation_Gradient_Inverse(K,l);
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
			col++;
		}
	}
    I_temp_11p_test_2_12_3*= fMaterial_Params[kMu];
}
void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_12_4()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_12_4=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int b=0;b<3;b++)
								{
									for(int a=0;a<3;a++)
									{
										for(int n=0;n<3;n++)
										{
											for(int m=0;m<3;m++)
											{
												I_temp_11p_test_2_12_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
												*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)
												*fDeformation_Gradient_n(a,K)*fDeformation_Gradient_Inverse(I,b)
												*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(n,m)
												*fDeformation_Gradient_Inverse(K,l);
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
			col++;
		}
	}
    I_temp_11p_test_2_12_4*= fMaterial_Params[kMu];
}


void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_13_1()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_13_1=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int b=0;b<3;b++)
								{
									for(int a=0;a<3;a++)
									{
										for(int p=0;p<3;p++)
										{
											for(int n=0;n<3;n++)
											{
												for(int m=0;m<3;m++)
												{
													I_temp_11p_test_2_13_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
													*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_n(b,K)*fDeformation_Gradient_Inverse(I,a)
													*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(m,p)
													*fCauchy_stress_tensor_current_n_IP(p,n)*fDeformation_Gradient_Inverse(K,l);
												}
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
			col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_13_2()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_13_2=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int b=0;b<3;b++)
								{
									for(int a=0;a<3;a++)
									{
										for(int p=0;p<3;p++)
										{
											for(int n=0;n<3;n++)
											{
												for(int m=0;m<3;m++)
												{
													I_temp_11p_test_2_13_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
													*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)*fCauchy_stress_tensor_current_n_IP(m,p)
													*fDeformation_Gradient_n(b,K)*fDeformation_Gradient_Inverse(I,a)
													*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(n,p)
													*fDeformation_Gradient_Inverse(K,l);
												}
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
			col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_13_3()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_13_3=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int b=0;b<3;b++)
								{
									for(int a=0;a<3;a++)
									{
										for(int n=0;n<3;n++)
										{
											for(int m=0;m<3;m++)
											{
												I_temp_11p_test_2_13_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
												*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)
												*fDeformation_Gradient_n(b,K)*fDeformation_Gradient_Inverse(I,a)
												*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(m,n)
												*fDeformation_Gradient_Inverse(K,l);
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
			col++;
		}
	}
    I_temp_11p_test_2_13_3*= fMaterial_Params[kMu];
}
void FSMicromorphic3DCurrConfigT:: Form_I_temp_11p_test_2_13_4()
{

    int row=0;
    int col=0;
    I_temp_11p_test_2_13_4=0.0;

    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int j=0;j<3;j++)
					{
						for(int q=0;q<3;q++)
						{
							for(int K=0;K<3;K++)
							{
								for(int b=0;b<3;b++)
								{
									for(int a=0;a<3;a++)
									{
										for(int n=0;n<3;n++)
										{
											for(int m=0;m<3;m++)
											{
												I_temp_11p_test_2_13_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,q)
												*fCauchy_stress_tensor_current_n_IP(q,k)*fdFYdCauchy_Stress(m,n)
												*fDeformation_Gradient_n(b,K)*fDeformation_Gradient_Inverse(I,a)
												*dev_Cauchy_stress_tr(a,b)*dev_Cauchy_stress_tr(n,m)
												*fDeformation_Gradient_Inverse(K,l);
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
		col++;
		}
	}
    I_temp_11p_test_2_13_4*= fMaterial_Params[kMu];
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////New Try for the contribution of variation of flow rule on the variation of Del gamma///////////
void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_1_1()
{

    int row=0;
    int col=0;
    II_temp_11p_1_1=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_1_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)
									*fCauchy_stress_tensor_current_n_IP(j,k)*I_temp_11p_1(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_1_2()
{

    int row=0;
    int col=0;
    II_temp_11p_1_2=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_1_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,p)
									*fCauchy_stress_tensor_current_n_IP(p,k)*I_temp_11p_1(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_1_3()
{

    int row=0;
    int col=0;
    II_temp_11p_1_3=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_1_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,p)
									*fdGdCauchy_Stress_tr(k,p)*I_temp_11p_1(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_1_4()
{

    int row=0;
    int col=0;
    II_temp_11p_1_4=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
						{
							for(int j=0;j<3;j++)
							{
								II_temp_11p_1_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)
								*fIdentity_matrix(j,k)*I_temp_11p_1(I,K)*fDeformation_Gradient_Inverse(K,l);
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}




void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_1_5()
{
    int row=0;
    int col=0;
    II_temp_11p_1_5=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							II_temp_11p_1_5(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,k)
							*I_temp_11p_1(I,K)*fDeformation_Gradient_Inverse(K,l);
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_1_6()
{
    int row=0;
    int col=0;
    II_temp_11p_1_6=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							II_temp_11p_1_6(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(k,j)
							*I_temp_11p_1(I,K)*fDeformation_Gradient_Inverse(K,l);
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_2_1()
{

    int row=0;
    int col=0;
    II_temp_11p_2_1=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_2_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)
									*fCauchy_stress_tensor_current_n_IP(j,k)*I_temp_11p_2(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_2_2()
{

    int row=0;
    int col=0;
    II_temp_11p_2_2=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_2_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,p)
									*fCauchy_stress_tensor_current_n_IP(p,k)*I_temp_11p_2(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_2_3()
{

    int row=0;
    int col=0;
    II_temp_11p_2_3=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_2_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,p)
									*fdGdCauchy_Stress_tr(k,p)*I_temp_11p_2(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_2_4()
{

    int row=0;
    int col=0;
    II_temp_11p_2_4=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
						{
							for(int j=0;j<3;j++)
							{
								II_temp_11p_2_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)
								*fIdentity_matrix(j,k)*I_temp_11p_2(I,K)*fDeformation_Gradient_Inverse(K,l);
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}




void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_2_5()
{
    int row=0;
    int col=0;
    II_temp_11p_2_5=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							II_temp_11p_2_5(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,k)
							*I_temp_11p_2(I,K)*fDeformation_Gradient_Inverse(K,l);
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_2_6()
{
    int row=0;
    int col=0;
    II_temp_11p_2_6=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							II_temp_11p_2_6(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(k,j)
							*I_temp_11p_2(I,K)*fDeformation_Gradient_Inverse(K,l);
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_3_1()
{

    int row=0;
    int col=0;
    II_temp_11p_3_1=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_3_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)
									*fCauchy_stress_tensor_current_n_IP(j,k)*I_temp_11p_3(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_3_2()
{

    int row=0;
    int col=0;
    II_temp_11p_3_2=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_3_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,p)
									*fCauchy_stress_tensor_current_n_IP(p,k)*I_temp_11p_3(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_3_3()
{

    int row=0;
    int col=0;
    II_temp_11p_3_3=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_3_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,p)
									*fdGdCauchy_Stress_tr(k,p)*I_temp_11p_3(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_3_4()
{

    int row=0;
    int col=0;
    II_temp_11p_3_4=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
						{
							for(int j=0;j<3;j++)
							{
								II_temp_11p_3_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)
								*fIdentity_matrix(j,k)*I_temp_11p_3(I,K)*fDeformation_Gradient_Inverse(K,l);
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}




void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_3_5()
{
    int row=0;
    int col=0;
    II_temp_11p_3_5=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							II_temp_11p_3_5(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,k)
							*I_temp_11p_3(I,K)*fDeformation_Gradient_Inverse(K,l);
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_3_6()
{
    int row=0;
    int col=0;
    II_temp_11p_3_6=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							II_temp_11p_3_6(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(k,j)
							*I_temp_11p_3(I,K)*fDeformation_Gradient_Inverse(K,l);
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_4_1()
{

    int row=0;
    int col=0;
    II_temp_11p_4_1=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_4_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)
									*fCauchy_stress_tensor_current_n_IP(j,k)*I_temp_11p_4(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_4_2()
{

    int row=0;
    int col=0;
    II_temp_11p_4_2=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_4_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,p)
									*fCauchy_stress_tensor_current_n_IP(p,k)*I_temp_11p_4(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_4_3()
{

    int row=0;
    int col=0;
    II_temp_11p_4_3=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_4_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,p)
									*fdGdCauchy_Stress_tr(k,p)*I_temp_11p_4(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_4_4()
{

    int row=0;
    int col=0;
    II_temp_11p_4_4=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
						{
							for(int j=0;j<3;j++)
							{
								II_temp_11p_4_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)
								*fIdentity_matrix(j,k)*I_temp_11p_4(I,K)*fDeformation_Gradient_Inverse(K,l);
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}




void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_4_5()
{
    int row=0;
    int col=0;
    II_temp_11p_4_5=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							II_temp_11p_4_5(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,k)
							*I_temp_11p_4(I,K)*fDeformation_Gradient_Inverse(K,l);
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_4_6()
{
    int row=0;
    int col=0;
    II_temp_11p_4_6=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							II_temp_11p_4_6(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(k,j)
							*I_temp_11p_4(I,K)*fDeformation_Gradient_Inverse(K,l);
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_5_1()
{

    int row=0;
    int col=0;
    II_temp_11p_5_1=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_5_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)
									*fCauchy_stress_tensor_current_n_IP(j,k)*I_temp_11p_5(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_5_2()
{

    int row=0;
    int col=0;
    II_temp_11p_5_2=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_5_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,p)
									*fCauchy_stress_tensor_current_n_IP(p,k)*I_temp_11p_5(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_5_3()
{

    int row=0;
    int col=0;
    II_temp_11p_5_3=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_5_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,p)
									*fdGdCauchy_Stress_tr(k,p)*I_temp_11p_5(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_5_4()
{

    int row=0;
    int col=0;
    II_temp_11p_5_4=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
						{
							for(int j=0;j<3;j++)
							{
								II_temp_11p_5_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)
								*fIdentity_matrix(j,k)*I_temp_11p_5(I,K)*fDeformation_Gradient_Inverse(K,l);
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}




void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_5_5()
{
    int row=0;
    int col=0;
    II_temp_11p_5_5=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							II_temp_11p_5_5(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,k)
							*I_temp_11p_5(I,K)*fDeformation_Gradient_Inverse(K,l);
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_5_6()
{
    int row=0;
    int col=0;
    II_temp_11p_5_6=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							II_temp_11p_5_6(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(k,j)
							*I_temp_11p_5(I,K)*fDeformation_Gradient_Inverse(K,l);
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_6_1()
{

    int row=0;
    int col=0;
    II_temp_11p_6_1=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_6_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)
									*fCauchy_stress_tensor_current_n_IP(j,k)*I_temp_11p_6(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_6_2()
{

    int row=0;
    int col=0;
    II_temp_11p_6_2=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_6_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,p)
									*fCauchy_stress_tensor_current_n_IP(p,k)*I_temp_11p_6(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_6_3()
{

    int row=0;
    int col=0;
    II_temp_11p_6_3=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_6_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,p)
									*fdGdCauchy_Stress_tr(k,p)*I_temp_11p_6(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_6_4()
{

    int row=0;
    int col=0;
    II_temp_11p_6_4=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
						{
							for(int j=0;j<3;j++)
							{
								II_temp_11p_6_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)
								*fIdentity_matrix(j,k)*I_temp_11p_6(I,K)*fDeformation_Gradient_Inverse(K,l);
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}




void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_6_5()
{
    int row=0;
    int col=0;
    II_temp_11p_6_5=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							II_temp_11p_6_5(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,k)
							*I_temp_11p_6(I,K)*fDeformation_Gradient_Inverse(K,l);
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_6_6()
{
    int row=0;
    int col=0;
    II_temp_11p_6_6=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							II_temp_11p_6_6(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(k,j)
							*I_temp_11p_6(I,K)*fDeformation_Gradient_Inverse(K,l);
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_7_1()
{

    int row=0;
    int col=0;
    II_temp_11p_7_1=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_7_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)
									*fCauchy_stress_tensor_current_n_IP(j,k)*I_temp_11p_7(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_7_2()
{

    int row=0;
    int col=0;
    II_temp_11p_7_2=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_7_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,p)
									*fCauchy_stress_tensor_current_n_IP(p,k)*I_temp_11p_7(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_7_3()
{

    int row=0;
    int col=0;
    II_temp_11p_7_3=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_7_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,p)
									*fdGdCauchy_Stress_tr(k,p)*I_temp_11p_7(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_7_4()
{

    int row=0;
    int col=0;
    II_temp_11p_7_4=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
						{
							for(int j=0;j<3;j++)
							{
								II_temp_11p_7_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)
								*fIdentity_matrix(j,k)*I_temp_11p_7(I,K)*fDeformation_Gradient_Inverse(K,l);
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}




void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_7_5()
{
    int row=0;
    int col=0;
    II_temp_11p_7_5=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							II_temp_11p_7_5(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,k)
							*I_temp_11p_7(I,K)*fDeformation_Gradient_Inverse(K,l);
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_7_6()
{
    int row=0;
    int col=0;
    II_temp_11p_7_6=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							II_temp_11p_7_6(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(k,j)
							*I_temp_11p_7(I,K)*fDeformation_Gradient_Inverse(K,l);
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_8_1()
{

    int row=0;
    int col=0;
    II_temp_11p_8_1=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_8_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)
									*fCauchy_stress_tensor_current_n_IP(j,k)*I_temp_11p_8(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_8_2()
{

    int row=0;
    int col=0;
    II_temp_11p_8_2=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_8_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,p)
									*fCauchy_stress_tensor_current_n_IP(p,k)*I_temp_11p_8(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_8_3()
{

    int row=0;
    int col=0;
    II_temp_11p_8_3=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_8_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,p)
									*fdGdCauchy_Stress_tr(k,p)*I_temp_11p_8(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_8_4()
{

    int row=0;
    int col=0;
    II_temp_11p_8_4=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
						{
							for(int j=0;j<3;j++)
							{
								II_temp_11p_8_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)
								*fIdentity_matrix(j,k)*I_temp_11p_8(I,K)*fDeformation_Gradient_Inverse(K,l);
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}




void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_8_5()
{
    int row=0;
    int col=0;
    II_temp_11p_8_5=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							II_temp_11p_8_5(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,k)
							*I_temp_11p_8(I,K)*fDeformation_Gradient_Inverse(K,l);
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_8_6()
{
    int row=0;
    int col=0;
    II_temp_11p_8_6=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							II_temp_11p_8_6(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(k,j)
							*I_temp_11p_8(I,K)*fDeformation_Gradient_Inverse(K,l);
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_9_1()
{

    int row=0;
    int col=0;
    II_temp_11p_9_1=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_9_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)
									*fCauchy_stress_tensor_current_n_IP(j,k)*I_temp_11p_9(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_9_2()
{

    int row=0;
    int col=0;
    II_temp_11p_9_2=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_9_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,p)
									*fCauchy_stress_tensor_current_n_IP(p,k)*I_temp_11p_9(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_9_3()
{

    int row=0;
    int col=0;
    II_temp_11p_9_3=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_9_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,p)
									*fdGdCauchy_Stress_tr(k,p)*I_temp_11p_9(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_9_4()
{

    int row=0;
    int col=0;
    II_temp_11p_9_4=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
						{
							for(int j=0;j<3;j++)
							{
								II_temp_11p_9_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)
								*fIdentity_matrix(j,k)*I_temp_11p_9(I,K)*fDeformation_Gradient_Inverse(K,l);
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}




void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_9_5()
{
    int row=0;
    int col=0;
    II_temp_11p_9_5=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							II_temp_11p_9_5(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,k)
							*I_temp_11p_9(I,K)*fDeformation_Gradient_Inverse(K,l);
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_9_6()
{
    int row=0;
    int col=0;
    II_temp_11p_9_6=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							II_temp_11p_9_6(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(k,j)
							*I_temp_11p_9(I,K)*fDeformation_Gradient_Inverse(K,l);
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_10_1()
{

    int row=0;
    int col=0;
    II_temp_11p_10_1=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_10_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)
									*fCauchy_stress_tensor_current_n_IP(j,k)*I_temp_11p_10(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_10_2()
{

    int row=0;
    int col=0;
    II_temp_11p_10_2=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_10_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,p)
									*fCauchy_stress_tensor_current_n_IP(p,k)*I_temp_11p_10(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_10_3()
{

    int row=0;
    int col=0;
    II_temp_11p_10_3=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_10_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,p)
									*fdGdCauchy_Stress_tr(k,p)*I_temp_11p_10(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_10_4()
{

    int row=0;
    int col=0;
    II_temp_11p_10_4=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
						{
							for(int j=0;j<3;j++)
							{
								II_temp_11p_10_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)
								*fIdentity_matrix(j,k)*I_temp_11p_10(I,K)*fDeformation_Gradient_Inverse(K,l);
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}




void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_10_5()
{
    int row=0;
    int col=0;
    II_temp_11p_10_5=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							II_temp_11p_10_5(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,k)
							*I_temp_11p_10(I,K)*fDeformation_Gradient_Inverse(K,l);
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_10_6()
{
    int row=0;
    int col=0;
    II_temp_11p_10_6=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							II_temp_11p_10_6(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(k,j)
							*I_temp_11p_10(I,K)*fDeformation_Gradient_Inverse(K,l);
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_11_1()
{

    int row=0;
    int col=0;
    II_temp_11p_11_1=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_11_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)
									*fCauchy_stress_tensor_current_n_IP(j,k)*I_temp_11p_11(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_11_2()
{

    int row=0;
    int col=0;
    II_temp_11p_11_2=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_11_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,p)
									*fCauchy_stress_tensor_current_n_IP(p,k)*I_temp_11p_11(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_11_3()
{

    int row=0;
    int col=0;
    II_temp_11p_11_3=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_11_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,p)
									*fdGdCauchy_Stress_tr(k,p)*I_temp_11p_11(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_11_4()
{

    int row=0;
    int col=0;
    II_temp_11p_11_4=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
						{
							for(int j=0;j<3;j++)
							{
								II_temp_11p_11_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)
								*fIdentity_matrix(j,k)*I_temp_11p_11(I,K)*fDeformation_Gradient_Inverse(K,l);
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}




void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_11_5()
{
    int row=0;
    int col=0;
    II_temp_11p_11_5=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							II_temp_11p_11_5(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,k)
							*I_temp_11p_11(I,K)*fDeformation_Gradient_Inverse(K,l);
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_11_6()
{
    int row=0;
    int col=0;
    II_temp_11p_11_6=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							II_temp_11p_11_6(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(k,j)
							*I_temp_11p_11(I,K)*fDeformation_Gradient_Inverse(K,l);
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_12_1()
{

    int row=0;
    int col=0;
    II_temp_11p_12_1=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_12_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)
									*fCauchy_stress_tensor_current_n_IP(j,k)*I_temp_11p_12(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_12_2()
{

    int row=0;
    int col=0;
    II_temp_11p_12_2=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_12_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,p)
									*fCauchy_stress_tensor_current_n_IP(p,k)*I_temp_11p_12(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_12_3()
{

    int row=0;
    int col=0;
    II_temp_11p_12_3=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_12_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,p)
									*fdGdCauchy_Stress_tr(k,p)*I_temp_11p_12(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_12_4()
{

    int row=0;
    int col=0;
    II_temp_11p_12_4=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
						{
							for(int j=0;j<3;j++)
							{
								II_temp_11p_12_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)
								*fIdentity_matrix(j,k)*I_temp_11p_12(I,K)*fDeformation_Gradient_Inverse(K,l);
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}




void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_12_5()
{
    int row=0;
    int col=0;
    II_temp_11p_12_5=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							II_temp_11p_12_5(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,k)
							*I_temp_11p_12(I,K)*fDeformation_Gradient_Inverse(K,l);
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_12_6()
{
    int row=0;
    int col=0;
    II_temp_11p_12_6=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							II_temp_11p_12_6(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(k,j)
							*I_temp_11p_12(I,K)*fDeformation_Gradient_Inverse(K,l);
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_13_1()
{

    int row=0;
    int col=0;
    II_temp_11p_13_1=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_13_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)
									*fCauchy_stress_tensor_current_n_IP(j,k)*I_temp_11p_13(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_13_2()
{

    int row=0;
    int col=0;
    II_temp_11p_13_2=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_13_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,p)
									*fCauchy_stress_tensor_current_n_IP(p,k)*I_temp_11p_13(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_13_3()
{

    int row=0;
    int col=0;
    II_temp_11p_13_3=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									II_temp_11p_13_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,p)
									*fdGdCauchy_Stress_tr(k,p)*I_temp_11p_13(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_13_4()
{

    int row=0;
    int col=0;
    II_temp_11p_13_4=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
						{
							for(int j=0;j<3;j++)
							{
								II_temp_11p_13_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)
								*fIdentity_matrix(j,k)*I_temp_11p_13(I,K)*fDeformation_Gradient_Inverse(K,l);
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}




void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_13_5()
{
    int row=0;
    int col=0;
    II_temp_11p_13_5=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							II_temp_11p_13_5(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,k)
							*I_temp_11p_13(I,K)*fDeformation_Gradient_Inverse(K,l);
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_II_temp_11p_13_6()
{
    int row=0;
    int col=0;
    II_temp_11p_13_6=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							II_temp_11p_13_6(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(k,j)
							*I_temp_11p_13(I,K)*fDeformation_Gradient_Inverse(K,l);
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FSMicromorphic3DCurrConfigT:: Form_I12p_1()
{

    int row=0;
    int col=0;
    I12p_1=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
							{
								for(int j=0;j<3;j++)
								{
									I12p_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)
									*fCauchy_stress_tensor_current_n_IP(j,k)*I11p_1(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_I12p_2()
{

    int row=0;
    int col=0;
    I12p_2=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									I12p_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,p)
									*fCauchy_stress_tensor_current_n_IP(p,k)*I11p_1(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_I12p_3()
{

    int row=0;
    int col=0;
    I12p_3=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									I12p_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,p)
									*fdGdCauchy_Stress_tr(k,p)*I11p_1(I,K)*fDeformation_Gradient_Inverse(K,l);
								}

							}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_I12p_4()
{

    int row=0;
    int col=0;
    I12p_4=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
						{
							for(int j=0;j<3;j++)
							{
								I12p_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)
								*fIdentity_matrix(j,k)*I11p_1(I,K)*fDeformation_Gradient_Inverse(K,l);
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}




void FSMicromorphic3DCurrConfigT:: Form_I12p_5()
{
    int row=0;
    int col=0;
    I12p_5=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							I12p_5(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,k)
							*I11p_1(I,K)*fDeformation_Gradient_Inverse(K,l);
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_I12p_6()
{
    int row=0;
    int col=0;
    I12p_6=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							I12p_6(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(k,j)
							*I11p_1(I,K)*fDeformation_Gradient_Inverse(K,l);
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

void FSMicromorphic3DCurrConfigT:: Form_I13p_1()
{

    int row=0;
    int col=0;
    I13p_1=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
							{
								for(int j=0;j<3;j++)
								{
									for(int m=0;m<3;m++)
									{
										for(int n=0;n<3;n++)
										{
											for(int i=0;i<3;i++)
											{
												for(int p=0;p<3;p++)
												{
													I13p_1(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)
													*fCauchy_stress_tensor_current_n_IP(j,k)*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_Inverse(I,i)
													*fCauchy_stress_tensor_current_n_IP(m,i)*fDeformation_Gradient_n(p,K)*fCauchy_stress_tensor_current_n_IP(p,n)*fDeformation_Gradient_Inverse(K,l);
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
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_I13p_2()
{

    int row=0;
    int col=0;
    I13p_2=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									for(int m=0;m<3;m++)
									{
										for(int n=0;n<3;n++)
										{
											for(int i=0;i<3;i++)
											{
												for(int q=0;q<3;q++)
												{
														I13p_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,p)
														*fCauchy_stress_tensor_current_n_IP(p,k)*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_Inverse(I,i)
														*fCauchy_stress_tensor_current_n_IP(m,i)*fDeformation_Gradient_n(q,K)*fCauchy_stress_tensor_current_n_IP(q,n)*fDeformation_Gradient_Inverse(K,l);
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
		col++;
		}
	}
}

/*

*/


void FSMicromorphic3DCurrConfigT:: Form_I13p_3()
{

    int row=0;
    int col=0;
    I13p_3=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									for(int m=0;m<3;m++)
									{
										for(int n=0;n<3;n++)
										{
											for(int i=0;i<3;i++)
											{
												for(int q=0;q<3;q++)
												{
													I13p_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fCauchy_stress_tensor_current_n_IP(j,p)
													*fdGdCauchy_Stress_tr(k,p)*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_Inverse(I,i)
													*fCauchy_stress_tensor_current_n_IP(m,i)*fDeformation_Gradient_n(q,K)*fCauchy_stress_tensor_current_n_IP(q,n)*fDeformation_Gradient_Inverse(K,l);
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
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_I13p_4()
{

    int row=0;
    int col=0;
    I13p_4=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int q=0;q<3;q++)
						{
							for(int j=0;j<3;j++)
							{
								for(int m=0;m<3;m++)
								{
									for(int n=0;n<3;n++)
									{
										for(int i=0;i<3;i++)
										{
											for(int p=0;p<3;p++)
											{
												I13p_4(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(q,q)
												*fIdentity_matrix(j,k)*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_Inverse(I,i)
												*fCauchy_stress_tensor_current_n_IP(m,i)*fDeformation_Gradient_n(p,K)*fCauchy_stress_tensor_current_n_IP(p,n)*fDeformation_Gradient_Inverse(K,l);
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
		col++;
		}
	}
}




void FSMicromorphic3DCurrConfigT:: Form_I13p_5()
{
    int row=0;
    int col=0;
    I13p_5=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							for(int m=0;m<3;m++)
							{
								for(int n=0;n<3;n++)
								{
									for(int i=0;i<3;i++)
									{
										for(int p=0;p<3;p++)
										{
											I13p_5(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,k)
											*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_Inverse(I,i)
											*fCauchy_stress_tensor_current_n_IP(m,i)*fDeformation_Gradient_n(p,K)*fCauchy_stress_tensor_current_n_IP(p,n)*fDeformation_Gradient_Inverse(K,l);
										}
									}
								}
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}



void FSMicromorphic3DCurrConfigT:: Form_I13p_6()
{
    int row=0;
    int col=0;
    I13p_6=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int j=0;j<3;j++)
						{
							for(int m=0;m<3;m++)
							{
								for(int n=0;n<3;n++)
								{
									for(int i=0;i<3;i++)
									{
										for(int p=0;p<3;p++)
										{
											I13p_6(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(k,j)
											*fdFYdCauchy_Stress(m,n)*fDeformation_Gradient_Inverse(I,i)
											*fCauchy_stress_tensor_current_n_IP(m,i)*fDeformation_Gradient_n(p,K)*fCauchy_stress_tensor_current_n_IP(p,n)*fDeformation_Gradient_Inverse(K,l);
										}
									}
								}
							}
						}
					}
				 row++;
				}
			}
		col++;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
void FSMicromorphic3DCurrConfigT:: Form_I12p_test_2()
{

    int row=0;
    int col=0;
    I12p_test_2=0.0;
    I12p_test_3=0.0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									I12p_test_2(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,p)
									*fCauchy_stress_tensor_current_n_IP(p,k)*I_temp_11p_4(I,K)*fDeformation_Gradient_Inverse(K,l);
								}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}

    col=0;
    row=0;
    for (int I=0;I<3;I++)
	{
		for(int l=0;l<3;l++)
		{
			row=0;
			for(int L=0;L<3;L++)
			{
				for(int k=0;k<3;k++)
				{
					for(int K=0;K<3;K++)
					{
						for(int p=0;p<3;p++)
							{
								for(int j=0;j<3;j++)
								{
									I12p_test_3(row,col)+=fDeformation_Gradient_Inverse(L,j)*fdGdCauchy_Stress_tr(j,p)
									*fCauchy_stress_tensor_current_n_IP(p,k)*I_temp_11p_4_Transpose(K,I)*fDeformation_Gradient_Inverse(K,l);
								}
							}
					}
				 row++;
				}
			}
		col++;
		}
	}
    I12p_test_3+= I12p_test_2;
}
*/
