/* $Id: FSMicromorphic3DTsaved.h,v 1.1 2010/11/10 17:51:09 isbuga Exp $ */
//DEVELOPMENT
#ifndef _FS_MICROMORPHIC_3D_T_H_
#define _FS_MICROMORPHIC_3D_T_H_

/* base classes */
#include "ElementBaseT.h"
#include "StringT.h"
#include "Traction_CardT.h"
#include "ShapeFunctionT.h"
#include "eIntegratorT.h"
#include "ModelManagerT.h"

#include "ifstreamT.h"
#include "ofstreamT.h"

#include "iAutoArrayT.h"
#include "ScheduleT.h"

/* direct members */
#include "LocalArrayT.h"
#include "GeometryT.h"
#include "dTensor3DT.h"// I added this header file to use Tensor built-in functions!

#include "VariArrayT.h"
#include "nVariArray2DT.h"
#include "VariLocalArrayT.h"


namespace Tahoe {

/* forward declarations */
class ShapeFunctionT;
class Traction_CardT;
class StringT;

/** FSMicromorphic3DT: This class contains a coupled finite deformation
 * micromorphic (displacement and micro-displacement gradient dofs) finite
 * element implementation in 3D.  The model description can be
 * found in Regueiro ASCE JEM 135:178-191 for pressure-sensitive elastoplasticity.
 **/

class FSMicromorphic3DT: public ElementBaseT
{

public:
/* material parameters */
    enum fMaterial_T         {
        kMu,
        kLambda,
        kKappa,
        kNu,
        kSigma_const,
        kTau,
        kEta,
        kRho_0,
        kTau1,
        kTau2,
        kTau3,
        kTau4,
        kTau5,
        kTau6,
        kTau7,
        kTau8,
        kTau9,
        kTau10,
        kTau11,
        kg,
        kg1,
        kg2,
        kg3,
        //add to this list
        kNUM_FMATERIAL_TERMS        };

//  enum fIntegrate_T         {
//      kBeta,
//      kGamma,
//      kNUM_FINTEGRATE_TERMS        };

    /** constructor */
    FSMicromorphic3DT(        const ElementSupportT& support );

    /** destructor */
    ~FSMicromorphic3DT(void);

    /** reference to element shape functions */
    const ShapeFunctionT& ShapeFunctionDispl(void) const;
    const ShapeFunctionT& ShapeFunctionMicro(void) const;

    /** echo input */
    void Echo_Input_Data (void);

    /** return true if the element contributes to the solution of the
     * given group. ElementBaseT::InGroup returns true if group is the
     * same as the group of the FieldT passed in to ElementBaseT::ElementBaseT. */
    virtual bool InGroup(int group) const;

    /* initialize/finalize time increment */
    /*@{*/
    virtual void InitStep(void);
    virtual void CloseStep(void);
    //virtual GlobalT::RelaxCodeT ResetStep(void); // restore last converged state

    /** element level reconfiguration for the current time increment */
    //virtual GlobalT::RelaxCodeT RelaxSystem(void);
    /*@}*/

    /** collecting element group equation numbers. See ElementBaseT::Equations
     * for more information */
    virtual void Equations( AutoArrayT<const iArray2DT*>& eq_d,
                AutoArrayT<const RaggedArray2DT<int>*>& eq_phi);

    /** return a const reference to the run state flag */
    virtual GlobalT::SystemTypeT TangentType(void) const;

    /** accumulate the residual force on the specified node
     * \param node test node
     * \param force array into which to assemble to the residual force */
    virtual void AddNodalForce(const FieldT& field, int node, dArrayT& force);

    /** returns the energy as defined by the derived class types */
    virtual double InternalEnergy(void);

    /** \name writing output */
    /*@{*/
    /** register element for output */
    virtual void RegisterOutput(void);

    /** write element output */
    virtual void WriteOutput(void);

    /** compute specified output parameter and send for smoothing */
    virtual void SendOutput(int kincode);
    /*@}*/

    /** return geometry and number of nodes on each facet */
    void FacetGeometry(ArrayT<GeometryT::CodeT>& facet_geometry, iArrayT& num_facet_nodes) const;

    /** return the geometry code */
    virtual GeometryT::CodeT GeometryCode(void) const;

    /*set active elements*/
    //virtual void SetStatus(const ArrayT<ElementCardT::StatusT>& status);

    /** initial condition/restart functions (per time sequence) */
    virtual void InitialCondition(void);

    /** mass types */
    enum MassTypeT {kNoMass = 0, /**< do not compute mass matrix */
            kConsistentMass = 1, /**< variationally consistent mass matrix */
                kLumpedMass = 2, /**< diagonally lumped mass */
             kAutomaticMass = 3  /**< select the mass type base on the time integration scheme */};
    MassTypeT static int2MassTypeT(int i);

    /** \name restart functions */
    /*@{*/
    /** write restart data to the output stream. Should be paired with
     * the corresponding ElementBaseT::ReadRestart implementation. */
    virtual void WriteRestart(ostream& out) const;

    /** read restart data to the output stream. Should be paired with
     * the corresponding ElementBaseT::WriteRestart implementation. */
    virtual void ReadRestart(istream& in);
    /*@}*/

    /** \name implementation of the ParameterInterfaceT interface */
    /*@{*/
    /** information about subordinate parameter lists */
    virtual void DefineSubs(SubListT& sub_list) const;

    /** return the description of the given inline subordinate parameter list */
    virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order,
                     SubListT& sub_lists) const;

    /** a pointer to the ParameterInterfaceT of the given subordinate */
    virtual ParameterInterfaceT* NewSub(const StringT& name) const;

    /** \name implementation of the ParameterInterfaceT interface */
    /*@{*/
    /** describe the parameters needed by the interface */
    virtual void DefineParameters(ParameterListT& list) const;

    /** accept parameter list */
    virtual void TakeParameterList(const ParameterListT& list);
    /*@}*/

protected:

    /** \name drivers called by ElementBaseT::FormRHS and ElementBaseT::FormLHS */
    /*@{*/
    /** form group contribution to the stiffness matrix */
    virtual void LHSDriver(GlobalT::SystemTypeT);

    /** form group contribution to the residual */
    virtual void RHSDriver(void);
    /*@}*/

    /** compute shape functions and derivatives */
    virtual void SetGlobalShape(void);

    void Select_Equations ( const int &iBalLinMom, const int &iBalFirstMomMom );

private:

    /** \name solution methods.
     * Both of these drivers assemble the LHS as well as the residual.
     */
    /*@{*/
    /** driver for staggered solution */
    void RHSDriver_staggered(void);

    /** driver for monolithic solution */
    void RHSDriver_monolithic(void);
    /*@}*/

protected:

    /* output control */
    iArrayT        fNodalOutputCodes;
    iArrayT        fElementOutputCodes;

private:

    /** Gradients and other matrices */
    dMatrixT fgrad_u, fgrad_u_n;
    dMatrixT fgrad_Phi, fgrad_Phi_n;

    dMatrixT fShapeDispl, fShapeDisplGrad, fShapeDisplGrad_t;
    dMatrixT fShapeDisplGrad_t_Transpose, fShapeDisplGradGrad, fShapeDisplGrad_temp;
    dMatrixT fShapeMicro,fShapeMicroGrad,fShapeMicroGrad_temp;

    dMatrixT fDefGrad, fDefGradInv, fDefGradInvMatrix;

    /** \name  values read from input in the constructor */
    /*@{*/
    /** element geometry */
    GeometryT::CodeT fGeometryCode_displ, fGeometryCodeSurf_displ;
    GeometryT::CodeT fGeometryCode_micro, fGeometryCodeSurf_micro;
    int fGeometryCode_displ_int,fGeometryCode_micro_int, fGeometryCodeSurf_displ_int;

    /** number of integration points */
    int        fNumIP_displ, fNumIPSurf_displ, fNumIP_micro, fNumIPSurf_micro;
    int knum_d_state, knum_i_state, knumstress, knumstrain;//,knumdispl;
    int num_sidesets;

    /*@}*/

    /** \name element displacements in local ordering */
    /*@{*/
    LocalArrayT u;                //solid displacement
    LocalArrayT u_n;         //solid displacement from time t_n
    LocalArrayT u_dot;         //solid velocity
    LocalArrayT u_dot_n;         //solid velocity from time t_n
    LocalArrayT u_dotdot;         //solid acceleration
    LocalArrayT u_dotdot_n;         //solid acceleration from time t_n
    LocalArrayT del_u;        //displacement increment including Newton-R update, i.e. del_u = u_{n+1}^{k+1} - u_n
    LocalArrayT Phi;        //micro-displacement-gradient dof
    LocalArrayT Phi_dot;        //micro-displacement-gradient first time derivative
    LocalArrayT Phi_dot_n;        //micro-displacement-gradient first time derivative from time t_n
    LocalArrayT Phi_dotdot;        //micro-displacement-gradient second derivative
    LocalArrayT Phi_dotdot_n;        //micro-displacement-gradient second derivative from time t_n
    LocalArrayT Phi_n;        //micro-displacement-gradient from time t_n
    LocalArrayT del_Phi;        //micro-displacement-gradient increment
    dArrayT                del_u_vec;          // vector form
    dArrayT                del_Phi_vec;        // vector form
    dArrayT                u_vec;          // solid displacement in vector form
    dArrayT                u_dot_vec;          // solid velocity in vector form
    dArrayT                u_dotdot_vec;          // solid acceleration in vector form
    dArrayT                Phi_vec;        // micro-displacement-gradient in vector form
    dArrayT                Phi_dot_vec;        // first derivative of micro-displacement-gradient in vector form
    dArrayT                Phi_dotdot_vec;        // second derivative of micro-displacement-gradient in vector form

    /*@}*/

    // problem size definitions
    int n_en_displ, n_en_displ_x_n_sd, n_sd_x_n_sd,n_sd_x_n_sd_x_n_sd,n_en_micro_x_n_sd_x_n_sd, n_en_micro_x_n_sd;
    int n_el, n_sd, n_sd_surf, n_en_surf;
    int n_en_micro, ndof_per_nd_micro, n_en_micro_x_ndof_per_nd_micro, ndof_per_nd_micro_x_n_sd;
    int step_number;
    int iConstitutiveModelType;

    //name of output vector
    StringT output;

    dArrayT fForces_at_Node;
    bool bStep_Complete;

    double time, fRho, fRho_0;

    void Get_Fd_ext ( dArrayT &fFd_ext );

    //-- Material Parameters
    dArrayT fMaterial_Params;

    //-- Newmark Time Integration Parameters
    dArrayT fIntegration_Params;

    /** \name shape functions wrt to current coordinates */
    /*@{*/
    /** shape functions and derivatives. The derivatives are wrt to the
      * reference coordinates */
    ShapeFunctionT* fShapes_displ;
    ShapeFunctionT* fShapes_micro;

    /** reference coordinates */
    LocalArrayT fInitCoords_displ, fInitCoords_micro;
    /** current coordinates */
    LocalArrayT fCurrCoords_displ, fCurrCoords_micro;
    /*@}*/

    /* Data Storage */
    ElementMatrixT fKdd, fKdphi;
    ElementMatrixT fKphid, fKphiphi;
    dArrayT         fFd_int;
    dArrayT         fFd_ext;
    dArrayT                fFphi_int;
    dArrayT                fFphi_ext;

    dArrayT                fGrad_disp_vector;
    dMatrixT    fGrad_disp_matrix;
    dArrayT         fDefGradInv_vector;
    dArrayT         fKirchhoff_vector;
    dArrayT         fSecond_Piola_vector;
    dArrayT                fChi_temp_vector;
    dArrayT                fFd_int_N1_vector;
    dArrayT                fTemp_vector_ndof_se;
    dArrayT                fTemp_vector_nen_micro;
    dArrayT                fTemp_vector_9x1;
    dArrayT                fPi_temp_transpose_vector;
    dArrayT                fGrad_1_J_vector;
    dArrayT                fTemp_nsd_vector;
    dArrayT                fFd_int_smallstrain_vector;
    dArrayT                fGravity_vector;
    dArrayT                fFd_int_G4_vector;
    dArrayT         fTemp_nine_values;
////////////////////////////////////
    dArrayT         fE_values;
    dArrayT         fVarepsilon;
////////////////////////////////////
    dArrayT         fTemp_six_values;
    dArrayT         fGradv_vector;
    dArrayT         fgradv_vector;

    dMatrixT        fDeformation_Gradient;
    dMatrixT        fDefGradT_9x9_matrix;
    dMatrixT        fRight_Cauchy_Green_tensor;
    dMatrixT        fRight_Cauchy_Green_tensor_Inverse;
    dMatrixT        fLeft_Cauchy_Green_tensor;
    dMatrixT        fLeft_Cauchy_Green_tensor_Inverse;
    dMatrixT        fDeformation_Gradient_Inverse;
    dMatrixT        fDeformation_Gradient_Transpose;
    dMatrixT        fDeformation_Gradient_Inverse_Transpose;
    dMatrixT        fDefGradInv_Grad_grad;
    dMatrixT        fDefGradInv_Grad_grad_Transpose;
    dMatrixT        fIdentity_matrix;
        dMatrixT    fSecond_Piola_tensor;
        dMatrixT    fTemp_matrix_nsd_x_nsd;
        dMatrixT    fTemp_matrix_nsd_x_1;
        dMatrixT    fTemp_matrix_ndof_se_x_ndof_se;
        dMatrixT    fTemp_matrix1_ndof_se_x_ndof_se;
        dMatrixT    fTemp_matrix_ndof_se_x_nen_micro;
        dMatrixT    fTemp_matrix1_nen_press_x_ndof_se;
        dMatrixT    fTemp_matrix_nsd_x_ndof_se;
        dMatrixT    fTemp_matrix_nsd_x_nen_micro;
        dMatrixT    fKirchhoff_tensor;
        dMatrixT    fIota_temp_matrix;
        dMatrixT    fVarpi_temp_matrix;


        dMatrixT    fIm_temp_matrix;
        dMatrixT    fHbar_temp_matrix;
        dMatrixT    fEll_temp_matrix;
        dMatrixT    fPi_temp_row_matrix;
        dMatrixT    fK_dd_G3_1_matrix;
        dMatrixT    fK_dd_G3_2_matrix;
        dMatrixT    fK_dd_G3_3_matrix;
        dMatrixT    fK_dd_G3_4_matrix;
        dMatrixT    fK_dd_G4_matrix;
        dMatrixT    fI_ij_column_matrix;
        dMatrixT    fShapeMicro_row_matrix;
        dMatrixT    fWp_temp_matrix;
        dMatrixT    fChi_temp_column_matrix;
        dMatrixT    fc_matrix;
        dMatrixT    fC_matrix;
        dMatrixT    fIm_Prim_temp_matrix;
        dMatrixT    fB_matrix;
        dMatrixT    fD_matrix;
        dMatrixT    fK_dd_BTDB_matrix;
    dMatrixT         fDefGradInv_column_matrix;
    dMatrixT         fDefGradInv_column_matrix_Transpose;
    dMatrixT        u_dotdot_column_matrix;
    dMatrixT        fXi_temp_matrix;
    dMatrixT        fVarsigma_temp_matrix;
    dMatrixT        fI_ijkl_matrix;
    dMatrixT        u_dot_column_matrix;
    dMatrixT        u_dot_column_matrix_Transpose;
    dMatrixT        fGravity_column_matrix;
    dMatrixT        fAleph_temp_matrix;
    dMatrixT        micro_dot_column_matrix;
    dMatrixT        fImath_temp_matrix;
    dMatrixT        fPf_0_matrix;
    //////////////////////////////////////////////////////////
    /////DEFINITIONS FOR MICROMORPHIC MATRICES////////////////
    //////////////////////////////////////////////////////////
    double KrDelta[3][3];
    double trdeltad;
    double trdeltaEp;

    double Counter;
    dArray2DT Counter_IPs_el_n;
    dArray2DT Counter_IPs_el;
    dArrayT Counter_IPs;


    //Varitional Matrices coming from the Balance of linear Momentum
    dMatrixT Tsigma_1;
    dMatrixT fG1_1;
    double SigN[3][3];//unsymmetric Cauchy stress tensor found at previous step
    dMatrixT SigN_m;
    dMatrixT Fn_m;
    dMatrixT Finv_m;
    dMatrixT deltaL;
    dMatrixT deltaL_Tr;
    dMatrixT tempSig;
    dMatrixT deltad;
    dMatrixT SigN_ar;

    dMatrixT Mat1;
    dMatrixT Mat2;
    dMatrixT Mat3;
    dMatrixT Mat4;
    dMatrixT Mat5;
    dArrayT RHS;
    dMatrixT Mat5_Inv;
    dMatrixT Sigma1;
    dMatrixT Sigma2;
    dMatrixT Sigma3;
    dMatrixT Sigma4;
    dMatrixT Sigma5;
    dMatrixT Sigma6;
    dMatrixT Sigma6_1;
    dMatrixT Sigma6_2;
    dMatrixT Sigma6_3;
    dMatrixT Sigma6_4;
    dMatrixT Sigma7;
    dMatrixT Sigma8;

    dArray2DT SigN_IPs;
    dArray2DT Sig_IPs;
    dArray2DT SigN_IPs_n;
    dArray2DT SigN_IPs_el;
    dArray2DT Sig_IPs_el;
    dArray2DT SigN_IPs_el_n;
    dArrayT Temp_Identity_array;


    dMatrixT SPiolaN;
    dMatrixT SPiola;

    dArray2DT SPiolaN_IPs;
    dArray2DT SPiola_IPs;
    dArray2DT SPiolaN_IPs_el;
    dArray2DT SPiola_IPs_el;
    dArray2DT SPiolaN_IPs_el_n;



    dMatrixT Sigma; // unsymetric Cauchy stress tensor at current step
    double Fn[3][3];
    dMatrixT Fn_ar;
    dArray2DT Fn_ar_IPs;
    dArray2DT F_ar_IPs;
    dArray2DT Fn_ar_IPs_el;
    dArray2DT F_ar_IPs_el;
    dArray2DT Fn_ar_IPs_el_n;

    double FnInv[3][3];
    dMatrixT FnInv_ar;
    dArray2DT FnInv_ar_IPs;
    dArray2DT FInv_ar_IPs;
//  dArray2DT FInv_ar_IPs_n;
    dArray2DT FnInv_ar_IPs_el;
    dArray2DT FInv_ar_IPs_el;
    dArray2DT FnInv_ar_IPs_el_n;

    double Finv[3][3];

    double Chi[3][3];
    dMatrixT Chi_m;
    double ChiInv[3][3];
    dMatrixT ChiInv_m;
    dMatrixT ChiN_m;
    dMatrixT deltaEp;
    dMatrixT deltaNu;
    double ChiN[3][3];
    dMatrixT ChinN_m;
    dMatrixT ChiN_ar;
    dMatrixT Chi_ar;
    dArray2DT ChiN_ar_IPs;
    dArray2DT Chi_ar_IPs;
    dArray2DT ChiN_ar_IPs_n;
    dArray2DT Chi_ar_IPs_el;
    dArray2DT ChiN_ar_IPs_el;
    dArray2DT ChiN_ar_IPs_el_n;

//  double  Trial[6];
    dMatrixT fIota_w_temp_matrix;
    dMatrixT fTemp_matrix_nudof_x_nchidof;
    dMatrixT fTemp_matrix_nchidof_x_nchidof;
    dMatrixT fTemp_matrix_nchidof_x_nudof;
    dMatrixT fTemp_matrix_nudof_x_nudof;

    dMatrixT fH1_Etagrad;
    dMatrixT TransShapeDisplGrad;
    dMatrixT Var_F;
    dMatrixT GRAD_Nuw;
    dMatrixT GRAD_Nuw_Tr;
    dMatrixT Finv_w; // to create Iota_w which is different than Iota because sequence of the components in wk,l
    dMatrixT Tsigma_2;
    dMatrixT fG1_2;//not being calculated yet
    dMatrixT Tsigma_3;
    dMatrixT fG1_3;// not being calculated yet
    dMatrixT TFn_1;// to be multiplied by (lamda+Tau)
    dMatrixT fG1_4;
    dMatrixT TFn_2;// to be multiplied by (Mu+sigma)
    dMatrixT fG1_5;
    dMatrixT TFn_3;// to be multiplied by (Mu+sigma)
    dMatrixT fG1_6;
    dMatrixT TChi_1;// to be multiplied by eta
    dMatrixT fG1_7;
    dMatrixT TFn_4;// to be multiplied by eta
    dMatrixT fG1_8;
    dMatrixT TChi_2;// to be multiplied by kappa
    dMatrixT fG1_9;
    dMatrixT TFn_5;// to be multiplied by kappa
    dMatrixT fG1_10;
    dMatrixT TChi_3;// to be multiplied by nu
    dMatrixT fG1_11;
    dMatrixT TFn_6;// to be multiplied by nu
    dMatrixT fG1_12;
    dMatrixT SigCurr;
    dMatrixT fG1_13;
    dMatrixT fG1_14;
 // Variational Matrices coming from the Balance of First Moment of Momentum
    double mn[3][3][3];
    dArrayT mn_ar;
    dArray2DT mn_IPs;
    dArray2DT mn_IPs_n;
    dArray2DT mn_IPs_el;
    dArray2DT mn_IPs_el_n;

    double Mnplus1[3][3][3];
    double Gamma[3][3][3];
    double GammaN[3][3][3];
    dArrayT GammaN_ar;

    dArray2DT GammaN_IPs;
    dArray2DT GammaN_IPs_n;
    dArray2DT GammaN_IPs_el;
    dArray2DT GammaN_IPs_el_n;

    double CCof[3][3][3][3][3][3];
    double GRAD_ChiN[3][3][3];// GRADIENT  in reference configuration!
    dArrayT GRAD_Chi_ar;
    dArrayT GRAD_ChiN_ar;
    dArray2DT GRAD_ChiN_ar_IPs;
    dArray2DT GRAD_Chi_ar_IPs;
    dArray2DT GRAD_ChiN_ar_IPs_n;
    dArray2DT GRAD_ChiN_ar_IPs_el;
    dArray2DT GRAD_Chi_ar_IPs_el;
    dArray2DT GRAD_ChiN_ar_IPs_el_n;

    double GRAD_Chi[3][3][3];// GRADIENT in reference configuration!



    dMatrixT fIota_eta_temp_matrix;
    dMatrixT Finv_eta; // to create Iota_eta


    dMatrixT Etagrad;
    dMatrixT Mm_1;
    dMatrixT Mm_2;
    dMatrixT Mm_3;
    dMatrixT Mm_4;
    dMatrixT Mm_5;
    dMatrixT Mm_6;
    dMatrixT Mm_7;
    dMatrixT Mm_71;
    dMatrixT Mm_72;
    dMatrixT Mm_73;
    dMatrixT Mm_74;
    dMatrixT Mm_75;
    dMatrixT Mm_76;
    dMatrixT Mm_77;
    dMatrixT Mm_78;
    dMatrixT Mm_8;
    dMatrixT Mm_9;
    dMatrixT Mm_10;
    dMatrixT Mm_11;
    dMatrixT Mm_12;
    dMatrixT Mm_13;
    dMatrixT Mm_14;
    dMatrixT Ru_1;//u
    dMatrixT Ru_2;
    dMatrixT Ru_3;
    dMatrixT RChi_1;
    dMatrixT Ru_4;
    dMatrixT RChi_2;
    dMatrixT Ru_5;
    dMatrixT Ru_6;
    dMatrixT Ru_7;
    dMatrixT Ru_8;
    dMatrixT Ru_9;
    dMatrixT RChi_3;
    dMatrixT Rs_sigma;
    dMatrixT R_gradu;
    dMatrixT R_Capital_Gamma_Chi;
    dMatrixT CapitalLambda;

    dMatrixT sn_sigman;
    dArray2DT sn_sigman_IPs;
    dArray2DT sn_sigman_IPs_n;
    dArray2DT sn_sigman_IPs_el;
    dArray2DT sn_sigman_IPs_el_n;

    dMatrixT s_sigma;

    dMatrixT fShapeDispl_Tr;


    dMatrixT NCHI;
    dMatrixT NCHI_Tr;
//    dMatrixT NCHI_eta;//  same with the one above no need!
//    dMatrixT GRAD_NCHI_Phi;//no need for this same with GRAD_NCHI
    dMatrixT GRAD_NCHI;


    dMatrixT fH1_1;
    dMatrixT fH1_2;
    dMatrixT fH1_3;
    dMatrixT fH1_4;
    dMatrixT fH1_5;
    dMatrixT fH1_6;
    dMatrixT fH1_7;
    dMatrixT fH1_71;
    dMatrixT fH1_72;
    dMatrixT fH1_73;
    dMatrixT fH1_74;
    dMatrixT fH1_75;
    dMatrixT fH1_76;
    dMatrixT fH1_77;
    dMatrixT fH1_78;
    dMatrixT fH1_8;
    dMatrixT fH1_9;
    dMatrixT fH1_10;
    dMatrixT fH1_11;
    dMatrixT fH1_12;
    dMatrixT fH1_13;
    dMatrixT fH1_14;

    dMatrixT fH2_1;
    dMatrixT fH2_2;
    dMatrixT fH2_3;
    dMatrixT fH2_4;
    dMatrixT fH2_5;
    dMatrixT fH2_6;
    dMatrixT fH2_7;
    dMatrixT fH2_8;
    dMatrixT fH2_9;
    dMatrixT fH2_10;
    dMatrixT fH2_11;
    dMatrixT fH2_12;
    dMatrixT fH2_13;
   // dMatrixT fH2_14;

    dMatrixT fH3_1;

    dArrayT Chi_vec;
    dArrayT GRAD_Chi_vec;
    //internal Force Matrices
    dArrayT G1;
    dArrayT Uint_1;
    dArrayT Uint_1_temp;
    dArrayT Pint_1_temp;
    dArrayT Pint_2_temp;
    dArrayT Pint_3_temp;
    dArrayT Uext_1;
    dArrayT Text;
    dArrayT Gext;

    dArrayT H1;
    dArrayT Pint_1;
    dArrayT H2;
    dArrayT Pint_2;
    dArrayT H3;
    dArrayT Pint_3;
    dArrayT Hext;
    dArrayT Pext;

    dMatrixT Lambda;
    dMatrixT Omega;
    //////////////////////////////////////////////////////////
    //////FINITE STRAIN ELASTICITY MATRICES START HERE/////////
    //////////////////////////////////////////////////////////
    //dMatrixT V_1;
    dMatrixT Jmat;
    dMatrixT KJmat;
    dMatrixT SPK;
    dMatrixT KirchhoffST;// The second Piola-Kirchhoff Matrix
    dMatrixT Temp_SPK;//temporary Matrix used in calculation of SPK
  //  dMatrixT FSF;
  //  dMatrixT LST;//Lagrangian strain tensor used in some functions to get rid of long name
    dMatrixT LagrangianStn;
    dMatrixT MicroStnTensor;//Micro-strain tensor
    dMatrixT PSI;//deformation measure PSI=Transpose(F).chi
    dMatrixT ChiM; //Micro-deformation tensor Chi ( used a different tensor this time )
    dMatrixT I1_1;
    dMatrixT I1_2;
    dMatrixT I1_3;
    dMatrixT I1_4;
    dMatrixT I1_5;
    dMatrixT I1_6;
    dMatrixT I1_7;
    dMatrixT I2_1;
    dMatrixT I1_8;
    dMatrixT I2_2;
    dMatrixT I1_9;
    dMatrixT I2_3;
/*********************/
    dMatrixT fFJ;
    dMatrixT fJF;
    dMatrixT fJ1_1;
    dMatrixT fJ1_2;
    dMatrixT fJ1_3;
    dMatrixT fJ1_4;
    dMatrixT fJ2_1;
    dMatrixT fJ1_5;
    dMatrixT fJ2_2;
    dMatrixT fJ1_6;
    dMatrixT fJ2_3;

    dArrayT Vint_1;
    dArrayT Vint_1_temp;
    dArrayT Vint_2;
    dArrayT Vint_2_temp;
    dArrayT Vint_3_temp;
    dArrayT Vint_3;
    dArrayT fV1;
    dArrayT fV2;
    dArrayT fV3;
    dMatrixT fKu_1;
    dMatrixT fKu_2;
    dMatrixT fKu_3;
    dMatrixT fKu_4;
    dMatrixT fKu_5;
    dMatrixT fKu_6;
    dMatrixT fKu_7;
    dMatrixT fKuphi_1;
    dMatrixT fKu_8;
    dMatrixT fKuphi_2;
    dMatrixT fKu_9;
    dMatrixT fKuphi_3;
    dMatrixT SIGMA_S;
    dMatrixT fKFJu;
    dMatrixT fKJFu;
    dMatrixT fKphiu_1;
    dMatrixT fKphiu_2;
    dMatrixT fKphiu_3;
    dMatrixT fKphiu_4;
    dMatrixT fKphiphi_1;
    dMatrixT fKphiu_5;
    dMatrixT fKphiphi_2;
    dMatrixT fKphiu_6;
    dMatrixT fKphiphi_3;

    dMatrixT fFM;
    dMatrixT fMF;
    dMatrixT fMchi;
    dMatrixT fEtaM;
    dMatrixT fMpu_1;
    dMatrixT fMpp_1;
    dMatrixT fMpu_2;
    dMatrixT fMpp_2;
    dMatrixT fMpu_3;
    dMatrixT fMpp_3;
    dMatrixT fMpu_4;
    dMatrixT fMpp_4;   
    dMatrixT fMpu_6;
    dMatrixT fMpp_6;        
    dMatrixT fKMphiu_3;
    dMatrixT fKMphiphi_3;
    dMatrixT fKMphiu_4;
    dMatrixT fKMphiphi_4;
    dMatrixT fKMphiu_6;
    dMatrixT fKMphiphi_6;

    dTensor3DT fMKLM;
    dTensor3DT GAMMA;
    dTensor3DT GRAD_CHIM;
    dTensor3DT fTemp_tensor_n_sd_x_n_sd_x_nsd;



    dMatrixT fKEtaM;
    dMatrixT fKMFphiu;
    dMatrixT fKMchiphiphi;
    dMatrixT fKMphiu_1;
    dMatrixT fKMphiphi_1;
    dMatrixT fKMphiu_2;
    dMatrixT fKMphiphi_2;
/////stress invariants variables////////
 double Cauchy_inv;
 double Rel_stres_inv;
 double Higher_orderT_inv;
 double temp_inv;
 dMatrixT devsigma;
 dMatrixT devRelsts;
 dTensor3DT  devmklm;
 dMatrixT s_sigma_temp;
 dTensor3DT fmklm;

//////////////////////////////////////////////////



    double trLST;
    double invJ;
////////////////////////////////////////////////
    double lambda_cap;
    double Mu_cap;
    double g1_;
    double g2_;
    double b1_;
    double b2_;
    double b3_;


/////////////////////////////////////////
 //   int element_number;
 //   int el_num;
  //  dMatrixT u_el;
  //  dArrayT u_element;

    double trsigma;
    double trs_sigma;
    double trmklm;
    dArrayT trvecmklm;
    dArrayT ftemp_u_element;
    dArrayT fState_variables;

    //////////////////////////////////////////////////////////
    /////DEFINITIONS FINISH HERE FOR MICROMORPHIC MATRICES////
    //////////////////////////////////////////////////////////

        /* to store fEulerian_strain_tensor_current_IP */
        dMatrixT    fEulerian_strain_tensor_current_IP;
        /* to store fEulerian_strain_IPs for each of the 27 IPs of each element */
        dArray2DT   fEulerian_strain_IPs;
        /* to store fCauchy_stress_tensor_current_IP */
        dMatrixT    fCauchy_stress_tensor_current_IP;
        /* to store fCauchy_stress_IPs for each of the 27 IPs of each element */
        dArray2DT   fCauchy_stress_IPs;
        ////////////////////////////////////
        dArray2DT  fE_values_IPs;
        dArray2DT  fVarepsilon_IPs;
        ////////////////////////////////////
        dArray2DT   fState_variables_IPs;
        dArray2DT   fEulerian_strain_Elements_IPs;
        dArray2DT   fCauchy_stress_Elements_IPs;
        dArray2DT   fState_variables_Elements_IPs;
        /////////////////////////////////////////
        dArray2DT  fE_values_Element_IPs;
        dArray2DT  fVarepsilon_Element_IPs;
        ////////////////////////////////////////

        dArray2DT fDisplacement_Element_IPs;
        dArrayT   fDisplacements_current_IPs;
        /* to store displacements for each of 27  IPs of each element
         ( probably it is better to do it only  for certain nodes but I do not know how to extract the coordinates of the nodes and calculate
         the shape functions a nodes)
         */
       // dArray2DT  fDisplacement_IPs;


    /** the solid displacement field */
    const FieldT* fDispl;

    /** the micro-displacement-gradient field */
    const FieldT* fMicro;

    /** \name state variable storage *
     * State variables are handled ABAQUS-style. For every iteration, the state
     * variables from the previous increment are passed to the element, which
     * updates the values in place. Each row in the array is the state variable
     * storage for all integration points for an element */
    /*@{*/
    dArray2DT fdState_new;
    dArray2DT fdState;

    iArray2DT fiState_new;
    iArray2DT fiState;
    /*@}*/

    /** \name connectivities */
    /*@{*/
    ArrayT<const iArray2DT*> fConnectivities_displ;
    ArrayT<const iArray2DT*> fConnectivities_micro;
    ArrayT<iArray2DT> fConnectivities_reduced;
    /*@}*/

    /** \name equations */
    /*@{*/
    ArrayT<iArray2DT> fEqnos_displ;
    ArrayT<iArray2DT> fEqnos_micro;
    /*@}*/

    /** \name element cards */
    /*@{*/
    AutoArrayT<ElementCardT> fElementCards_displ;
    AutoArrayT<ElementCardT> fElementCards_micro;
    /*@}*/

    /** \name output */
    /*@{*/
    /** output ID */
    int fOutputID;

    /** integration point stresses. Calculated and stored during
     * FSMicromorphic2DT::RHSDriver */
    dArray2DT fIPVariable;
    /*@}*/

    /** \name prescribed plastic gradient side set ID */
    /*@{*/
    ArrayT<StringT> fSideSetID;

    /** prescribed micro-displacement-gradient weight over the side set;
        the direction is defined by {n1,n2,n3} ?? */
    //ArrayT<iArray2DT> fMicroWght;

    /** for each side set, the global nodes on the faces in the set */
    ArrayT<iArray2DT> fMicroFaces;

    /** equation numbers for the nodes on each face */
    ArrayT<iArray2DT> fMicroFaceEqnos;

    /** side set elements */
    ArrayT<iArrayT> fSideSetElements;

    /** side set faces */
    ArrayT<iArrayT> fSideSetFaces;
    /*@}*/

    /** write output for debugging */
    /*@{*/
    /** output file stream */
    ofstreamT fs_micromorph3D_out;
    ofstreamT fs_micromorph3DMn_out;

    /** line output formating variables */
    int outputPrecision, outputFileWidth;
    /*@}*/

    void Form_solid_shape_functions(const double* &shapes_displ_X);
    void Form_Gradient_of_solid_shape_functions(const dMatrixT &fShapeDisplGrad_temp);
    void Form_Gradient_t_of_solid_shape_functions(const dMatrixT &fShapeDisplGrad_temp);
    void Form_micro_shape_functions(const double* &shapes_micro_X);
    void Form_deformation_gradient_tensor(void);
    void Form_Grad_grad_transformation_matrix(void);
    void Form_fDefGradT_9x9_matrix(void);
    void Form_deformation_gradient_inv_vector(void);
    void Form_kirchhoff_stress_vector(void);
    void Form_Varpi_temp_matrix(void);
    void Form_Im_temp_matrix(void);
    void Form_Hbar_temp_matrix(void);
    void Form_Ell_temp_matrix(void);
    void Form_C_matrix(const double& J_Prim);
    void Form_c_matrix(void);
    void Form_Im_Prim_temp_matrix(void);
    void Form_D_matrix(void);
    void Form_B_matrix(void);
    void Extract_six_values_from_symmetric_tensor(const dMatrixT &fTensor,dArrayT& fTemp_six_values);
    void Put_values_In_dArrayT_vector(const dArray2DT &f2DArrayT,const int& e,const int& IP,dArrayT& fArrayT);
    void Put_values_In_Array(const dArray2DT &f2DArrayT,const int& e,const int& IP,dArrayT& fArrayT);
    void Form_gradv_vector(void);
    void Form_Xi_temp_matrix(void);
    void Form_Varsigma_temp_matrix(void);
    void Form_I_ijkl_matrix(void);
    void Compute_norm_of_array(double& norm,const LocalArrayT& B);
    //////////////////////////////////////////////////////////
    /////FUNCTIONS  FOR MICROMORPHIC MATRICES////////////////
    //////////////////////////////////////////////////////////
    //Forming the Matrices coming from the Balance of Linear Momentum
    void Form_Gamma_tensor3D(void);
    void Form_micro_deformation_tensor_Chi(void);
    void Form_Chi_inv_matrix(void);
    void Form_GRAD_Chi_matrix(void);
    void Form_double_Finv_from_Deformation_tensor_inverse(void);
    void Form_GRAD_Nuw_matrix(const dMatrixT &fShapeDisplGrad_temp);
    void Form_Finv_w_matrix(void);

    void Mapping_double_and_Array(const int condition);//
//  void Mapping_double_and_Array(double& dmat, dArrayT& fArrayT,const int& dim,const int& condition);
    void Form_deformation_tensors_arrays(const int condition);//

    void Form_KroneckerDelta_matrix(void);
    void Form_Var_F_tensor(void);
    void Form_Tsigma_1_matrix(void);
    void Form_fG1_1_matrix(void);//not defined yet
    void Form_Tsigma_2_matrix(void);
    void Form_fG1_2_matrix(void);
    void Form_Tsigma_3_matrix(void);
    void Form_fG1_3_matrix(void);
    void Form_TFn_1_matrix(void);
    void Form_fG1_4_matrix(void);
    void Form_TFn_2_matrix(void);
    void Form_fG1_5_matrix(void);
    void Form_TFn_3_matrix(void);
    void Form_fG1_6_matrix(void);
    void Form_TChi_1_matrix(void);
    void Form_fG1_7_matrix(void);
    void Form_TFn_4_matrix(void);
    void Form_fG1_8_matrix(void);
    void Form_TChi_2_matrix(void);
    void Form_fG1_9_matrix(void);
    void Form_TFn_5_matrix(void);
    void Form_fG1_10_matrix(void);
    void Form_TChi_3_matrix(void);
    void Form_fG1_11_matrix(void);
    void Form_TFn_6_matrix(void);
    void Form_fG1_12_matrix(void);
    void Form_SigCurr_matrix(void);
    void Form_fG1_13_matrix(void);
    //Forming the Matrices coming from the Balance of First Moment of Momentum
    void Form_Gradient_of_micro_shape_eta_functions(const dMatrixT &fShapeMicroGrad);
//  void Form_NCHI_eta_matrix(const dMatrixT &fShapeMicro_row_matrix);//same with the one below
    void Form_NCHI_matrix(const dMatrixT &fShapeMicro_row_matrix); //shape function matrice for micro deformatin {ETA}=[NCHI].{alpha}

    void Form_Etagrad_matrix(void);
    void Form_Mm_1_matrix(void);// needs to be multiplied by "-" and J
    void Form_Mm_2_matrix(void);// needs to be multiplied by J
    void Form_Mm_3_matrix(void);// needs to be multiplied by J
    void Form_Mm_4_matrix(void);// needs to be multiplied by J
    void Form_Mm_5_matrix(void);// needs to be multiplied by J
    void Form_Mm_6_matrix(void);// needs to be multiplied by J
    void Form_Mm_7_matrix(void);
    void Form_Mm_71_matrix(void);
    void Form_Mm_72_matrix(void);
    void Form_Mm_73_matrix(void);
    void Form_Mm_74_matrix(void);
    void Form_Mm_75_matrix(void);
    void Form_Mm_76_matrix(void);
    void Form_Mm_77_matrix(void);
    void Form_Mm_78_matrix(void);
    void Form_Mm_8_matrix(void);
    void Form_Mm_9_matrix(void);
    void Form_Mm_10_matrix(void);
    void Form_Mm_11_matrix(void);
    void Form_Mm_12_matrix(void);
    void Form_Mm_13_matrix(void);
    void Form_Mm_14_matrix(void);//something should be changed in the loop due to div(du)!!! changed done ok!
    void Form_Ru_1_matrix(void);
    void Form_Ru_2_matrix(void);
    void Form_Ru_3_matrix(void);
    void Form_RChi_1_matrix(void);// needs to be multiplied by Kappa and J
    void Form_Ru_4_matrix(void);// needs to be multiplied by Kappa and J
    void Form_RChi_2_matrix(void);//needs to be multiplied by Nu and J
    void Form_Ru_5_matrix(void);// needs to be multiplied by Nu and J
    void Form_Ru_6_matrix(void);
    void Form_Ru_7_matrix(void);
    void Form_Ru_8_matrix(void);
    void Form_Ru_9_matrix(void);
    void Form_RChi_3_matrix(void);
    void Form_Rs_sigma_matrix(void);
    void Form_R_Capital_Lambda_Chi_matrix(void);// DO NOT multiply with J !!!
    void Form_Finv_eta_matrix(void);
    void Form_CapitalLambda_matrix(void);
    void Form_CCof_tensor(void);
    void Form_SPiola_matrix(void);

//    void Form_GRAD_NCHI_Phi_matrix(const  dMatrixT &fShapeMicroGrad); no need for this

    void Form_G1_matrix(void);
    void Form_H1_matrix(void);
    void Form_H2_matrix(void);
    void Form_H3_matrix(void);


    //////////////////////////////////////////////////////////
    /////FINITE STRAIN ELASTICITY MATRICES START /////////////
    //////////////////////////////////////////////////////////
    void Form_fV1(void);
    void Form_fV2(void);
    void Form_fV3(void);
    void Form_I1_1_matrix(void);
    void Form_Second_Piola_Kirchhoff_SPK(void);
    void Form_ChiM(void);
    void Form_I1_1(void);
    void Form_I1_2(void);
    void Form_I1_3(void);
    void Form_I1_4(void);
    void Form_I1_5(void);
    void Form_I1_6(void);
    void Form_I1_7(void);
    void Form_I2_1(void);
    void Form_I1_8(void);
    void Form_I2_2(void);
    void Form_I1_9(void);
    void Form_I2_3(void);
    void Form_SIGMA_S(void);
    void Form_fFJ(void);
    void Form_fJF(void);
    void Form_fJ1_1(void);
    void Form_fJ1_2(void);
    void Form_fJ1_3(void);
    void Form_fJ1_4(void);
    void Form_fJ2_1(void);
    void Form_fJ1_5(void);
    void Form_fJ2_2(void);
    void Form_fJ1_6(void);
    void Form_fJ2_3(void);
    void Form_fFM(void);
    void Form_fMF(void);
    void Form_fMKLM(void);
    void Form_GAMMA(void);
    void Form_fEtaM(void);
    void Form_fMchi(void);
    void Form_fMpu_1(void);
    void Form_fMpp_1(void);
    void Form_fMpu_2(void);
    void Form_fMpp_2(void);
    void Form_fMpu_3(void);
    void Form_fMpp_3(void);
    void Form_fMpu_4(void);
    void Form_fMpp_4(void);
    void Form_fMpu_6(void);
    void Form_fMpp_6(void);



    void Form_Jmat(void);


    ////////////////////////////////////////////////////////
    /////////// functions to calculate stress measures ////
    void Calculate_Cauchy_INV(void);
    void Calculate_stress_diff_INV(void);
    void Calculate_higher_order_tensor_INV(void);
    void Calculate_fmklm(void);


    /////////////////////////////////////////////////////////



    //////////////////////////////////////////////////////////
    /////FUNCTIONS FINISH HERE FOR MICROMORPHIC MATRICES////
    //////////////////////////////////////////////////////////

protected:

    /** extract natural boundary condition information */
    void TakeNaturalBC(const ParameterListT& list);

    /** apply traction boundary conditions to displacement equations */
    void ApplyTractionBC(void);

    /** update traction BC data for displacement equations */
    void SetTractionBC(void);

    /* traction data */
    ArrayT<Traction_CardT> fTractionList;
    int fTractionBCSet;

    /** \name arrays with local ordering */
    /*@{*/
    LocalArrayT fLocInitCoords;   /**< initial coords with local ordering */
    LocalArrayT fLocDisp;              /**< solid displacements with local ordering  */
    /*@}*/

    /** \name work space */
    /*@{*/
    dArrayT fNEEvec; /**< work space vector: [element DOF] */
    dArrayT fDOFvec; /**< work space vector: [nodal DOF]   */
    /*@}*/

};


inline const ShapeFunctionT& FSMicromorphic3DT::ShapeFunctionDispl(void) const
{
#if __option(extended_errorcheck)
    if (!fShapes_displ)
        ExceptionT::GeneralFail("FSMicromorphic3DT::ShapeFunctionDispl", "no displ shape functions");
#endif
    return *fShapes_displ;
}

inline const ShapeFunctionT& FSMicromorphic3DT::ShapeFunctionMicro(void) const
{
#if __option(extended_errorcheck)
    if (!fShapes_micro)
        ExceptionT::GeneralFail("FSMicromorphic3DT::ShapeFunctionMicro", "no micro shape functions");
#endif
    return *fShapes_micro;
}

/* return the geometry code */
inline GeometryT::CodeT FSMicromorphic3DT::GeometryCode(void) const
{ return fGeometryCode_displ; }


} // namespace Tahoe
#endif /* _FS_MICROMORPHIC_3D_T_H_ */



