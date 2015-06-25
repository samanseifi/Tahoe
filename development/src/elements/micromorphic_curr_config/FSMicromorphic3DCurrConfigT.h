/* $Id: FSMicromorphic3DCurrConfigT.h,v 1.3 2013/10/01 00:55:05 tahoe.fash5153 Exp $ */
//DEVELOPMENT
#ifndef _FS_MICROMORPHIC_3D_CURR_CONFIG_T_H_
#define _FS_MICROMORPHIC_3D_CURR_CONFIG_T_H_

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

/** FSMicromorphic3DCurrConfigT: This class contains a coupled finite deformation
 * micromorphic (displacement and micro-displacement gradient dofs) finite
 * element implementation in 3D.  The model description can be
 * found in Regueiro ASCE JEM 135:178-191 for pressure-sensitive elastoplasticity.
 **/

class FSMicromorphic3DCurrConfigT: public ElementBaseT
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
        // plasticity parameters
        kHc,
        kc0,
        kHc_chi,
        kc0_chi,
        kHGc_chi,
        kGc0_chi1,
        kGc0_chi2,
        kGc0_chi3,
        //kZ0c,
        kFphi,
        kDpsi,
        kFphi_chi,
        kDpsi_chi,
        kFGphi_chi,
        kDGpsi_chi,
        //
        kg,
        kg1,
        kg2,
        kg3,
        //add to this list
        kNUM_FMATERIAL_TERMS        };

    enum fMaterialState_T         {
    //    kkappa,
        kc,
        kc_chi,
        //    kZkappa,
    //    kZc,
    //    khkappa,
        khc,
        khc_chi,
    //    kIntrinsic_Perm,
    //    kJ,
    //    kJp,
    //    kphi_s,
    //    kphi_f,
    //    kDevSS,
    //    kMeanS,
    //    kEpsVolp,
        kDelgamma,
        kDelgammachi,
        ktrCauchy_Stress,
        kNorm_Dev_Cauchy_Stress,
        ktrRel,
        kRel_inv,
        ktrm,
        km_inv,
        ktrS,
        kinvdevS,
        ktrSIGMA_S,
        kinvdevSIGMA_S,
        kinvtrM,
        kinvdevM,
        kinvPhi,
        kinvGPhi,
        ktreps,
        kdeveps,
	kinvtrgammastn,
	kinvdevgammastn,
    kNUM_FMATERIAL_STATE_TERMS
    };


//  enum fIntegrate_T         {
//      kBeta,
//      kGamma,
//      kNUM_FINTEGRATE_TERMS        };

    /** constructor */
    FSMicromorphic3DCurrConfigT(        const ElementSupportT& support );

    /** destructor */
    ~FSMicromorphic3DCurrConfigT(void);

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
    int step_number,global_iteration;
    int iConstitutiveModelType;
    //double Alpha;

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
    dArrayT         fFphi_int;
    dArrayT         fFphi_ext;

    dArrayT         fGrad_disp_vector;
    dMatrixT        fGrad_disp_matrix;

    dArrayT         fTemp_nine_values;
    dArrayT         fTemp_six_values;

    dMatrixT        fDeformation_Gradient;
    dMatrixT        fRight_Cauchy_Green_tensor;
    dMatrixT        fRight_Cauchy_Green_tensor_Inverse;
    dMatrixT        fRight_Elastic_Cauchy_Green_tensor;
    dMatrixT        fRight_Elastic_Cauchy_Green_tensor_tr;

    dMatrixT        fLeft_Cauchy_Green_tensor;
    dMatrixT        fLeft_Cauchy_Green_tensor_Inverse;
    dMatrixT        fDeformation_Gradient_Inverse;
    dMatrixT        fDeformation_Gradient_Transpose;
    dMatrixT        fDefGradInv_Grad_grad;

    dMatrixT        fIdentity_matrix;
    dMatrixT        fTemp_matrix_nsd_x_nsd;
    dMatrixT        fTemp_matrix_nsd_x_nsd2;
    dMatrixT        fTemp_matrix_nsd_x_nsd3;

    dMatrixT        fShapeMicro_row_matrix;


    dMatrixT        u_dotdot_column_matrix;
    dMatrixT        u_dot_column_matrix;
    dMatrixT        u_dot_column_matrix_Transpose;

    ///////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    ///////////// Plasticity Parameters//////////////////////
    ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////

//////////////////////////// A(delgamma)^2+B(delgamma)+D = 0//////////////////////////
    double			fA_coeff;
    double			fB_coeff;
    double			fD_coeff;
    double			ftemp_D_coeff;
    double			ftemp_D_coeff1;
    double			ftemp_B_coeff;
    double			fDelta_function;
/////////////////////////////////////////////////////////////////////////////////////



    dArray2DT       fdGdCauchy_Stress_Elements_IPs;
    dArray2DT       fdGdCauchy_Stress_n_Elements_IPs;
    dArray2DT		fdGdCauchy_Stress_IPs;
    dArray2DT		fdGdCauchy_Stress_n_IPs;
    dMatrixT		fdGdCauchy_Stress;
    dMatrixT		fdGdCauchy_Stress_n;
    dMatrixT		fdGdCauchy_Stress_n_Transpose;
    dMatrixT		fdGdCauchy_Stress_Transpose;
    dMatrixT		fdGdCauchy_Stress_tr;
    dMatrixT		fdGdCauchy_Stress_tr_Transpose;
    double			fdGdCauchy_Stress_tr_trace;
    double			fdGdCauchy_Stress_n_trace;
    double			fdGdCauchy_Stress_trace;

    dMatrixT		fdFYdCauchy_Stress;
    dArray2DT       fdFYdCauchy_Stress_Elements_IPs;
    dArray2DT       fdFYdCauchy_Stress_n_Elements_IPs;
    dArray2DT		fdFYdCauchy_Stress_IPs;
    dArray2DT		fdFYdCauchy_Stress_n_IPs;
    dMatrixT		fdFYdCauchy_Stress_n;
    dMatrixT		fdFYdCauchy_Stress_n_Transpose;



    dArray2DT   	fCauchy_stress_Elements_IPs;
    dArray2DT   	fCauchy_stress_IPs;
    dMatrixT    	fCauchy_stress_tensor_current_IP;
    dArray2DT   	fCauchy_stress_Elements_n_IPs;
    dArray2DT   	fCauchy_stress_n_IPs;
    dMatrixT    	fCauchy_stress_tensor_current_n_IP;
    dMatrixT    	fdev_Cauchy_stress_tensor_current_IP;
    dMatrixT    	fdev_Cauchy_stress_tensor_current_n_IP;
    dMatrixT		fdCauchy_stressdDelgamma;
    double			fCauchy_stress_tensor_current_IP_trace;


    double			fNorm_dev_Cauchy_stress_tensor_current_IP;
    double			fNorm_dev_Cauchy_stress_tensor_current_IP_tr;



    double			Beta;
    double			Aphi;
    double			Bphi;
    double 			Bpsi;
    double 			Apsi;


    dMatrixT    	fCauchy_stress_tensor_current_IP_tr;
    double       	mean_Cauchy_stress_tr;
    dMatrixT		dev_Cauchy_stress_tr;
    double			Cauchy_Stress_Norm_tr;
    dMatrixT		fdGdS_tr;
    dMatrixT		fdGdS_tr_transpose;
    double			mean_Cauchy_stress_n;
    double			mean_Cauchy_stress;
    double			fdmean_Cauchy_stressdDelgamma;
    dMatrixT		fdDev_Cauchy_stressdDelgamma;

    double			dcdDelgamma;
    double			fNorm_dDev_Cauchy_stressdDelgamma;
    double			dFYdDelgamma;
    dMatrixT		fElastic_Velocity_Gradient_current_IP;
    double			fElastic_Velocity_Gradient_current_IP_trace;
    dMatrixT		fElastic_Velocity_Gradient_current_IP_transpose;

    dMatrixT		dev_Cauchy_stress_n;
    dMatrixT		Predictor_stress_terms;
    double			fVelocity_Gradient_current_IP_trace;
    double			Temp_trace_value;
    dMatrixT		fVelocity_Gradient_current_IP_transpose;
    dMatrixT		Corrector_stress_terms;
    double			Predictor_mean_stress_terms;
    double			Corrector_mean_stress_terms;
    double			fTemp_matrix_one_x_one;
    double			fDelgamma_current_configuration_root_one;
    double			fDelgamma_current_configuration_root_two;
    double			fDelgamma_current_configuration;
    double			fTemp_delgamma_coeff;
    double			fTemp_delgamma_coeff1;
    double			Je;
    double			fYield_function_current_configuration;
    double			cohesion_current_configuration;
    double			fYield_function_tr;
    double			fYield_function;
    double			fF_tr_fact;
    double			fdelDelgamma;
    double			fDelgamma;
    double			iter_count;
    double 			cohesion;
    double			fVariation_Delgamma;
    dMatrixT    	fVelocity_Gradient_current_IP;
    dMatrixT		fSymmetric_part_Velocity_Gradient_current_IP;
    double			fSymmetric_part_Velocity_Gradient_current_IP_trace;
    double			fdFYdc;
    double			Comp11;
    double			Comp22;
    double			Comp33;
    double			Comp44;
    double			J;

    dArray2DT       fDeformation_Gradient_Elements_IPs;
    dArray2DT       fDeformation_Gradient_n_Elements_IPs;
    dArray2DT		fDeformation_Gradient_IPs;
    dArray2DT		fDeformation_Gradient_n_IPs;
    dMatrixT		fDeformation_Gradient_n;
    dMatrixT		fDel_Deformation_Gradient;


    ///Runge Kutta matrices
    dMatrixT		fTemp_Runge_Kutta_K1_nsd_x_nsd;
    dMatrixT		fTemp_Runge_Kutta_K2_nsd_x_nsd;
    dMatrixT		fTemp_Runge_Kutta_K3_nsd_x_nsd;
    dMatrixT		fTemp_Runge_Kutta_K4_nsd_x_nsd;

    dMatrixT		fCauchy_stress_tensor_current_IP_from_piola_stress;
    dMatrixT    	fLeft_Cauchy_Green_tensor_current_IP;
    dMatrixT    	fLeft_Cauchy_Green_tensor_current_IP_transpose;
    dMatrixT		fLeft_Cauchy_Green_tensor_current_IP_Inverse;
    double			fLeft_Cauchy_Green_tensor_current_IP_Trace;

    dMatrixT		fRight_Elastic_Cauchy_Green_tensor_Inverse;

    dMatrixT		fCauchy_stress_tensor_current_IP_Predictor;
    dMatrixT		fCauchy_stress_tensor_current_IP_Corrector;
 ///////////////////// Global Consistent tangent matrices/////////////////////
    //////////////////////////////////////////////////////////////////

    dArrayT Vintp_1;
    dArrayT Vintp_1_temp;
    dArrayT fV1p;

    dArrayT Vinte_1;
    dArrayT Vinte_1_temp;
    dArrayT fV1e;


    dMatrixT		IJe_1;
    dMatrixT		IJe_2;
    dMatrixT		IJe_3;
    dMatrixT		IJe_4;
    dMatrixT		IJe_5;
    dMatrixT		IJe_6;
    dMatrixT		IJe_7;
    dMatrixT		IJe_8;



    dMatrixT		I1p_1;
    dMatrixT		I1p_2;
    dMatrixT		I1p_3;
    dMatrixT		I1p_4;
    dMatrixT		I1p_5;
    dMatrixT		I1p_6;


    dMatrixT		I2p_1;
    dMatrixT		I2p_2;
    dMatrixT		I2p_3;
    dMatrixT		I2p_4;
    dMatrixT		I2p_5;
    dMatrixT		I2p_6;


    dMatrixT		I3p_1;
    dMatrixT		I3p_2;
    dMatrixT		I3p_3;
    dMatrixT		I3p_4;
    dMatrixT		I3p_5;
    dMatrixT		I3p_6;


    dMatrixT		I4p_1;
    dMatrixT		I4p_2;
    dMatrixT		I4p_3;
    dMatrixT		I4p_4;
    dMatrixT		I4p_5;
    dMatrixT		I4p_6;


    dMatrixT		I5p_1;
    dMatrixT		I5p_2;
    dMatrixT		I5p_3;
    dMatrixT		I5p_4;
    dMatrixT		I5p_5;
    dMatrixT		I5p_6;


    dMatrixT		I6p_1;
    dMatrixT		I6p_2;
    dMatrixT		I6p_3;
    dMatrixT		I6p_4;
    dMatrixT		I6p_5;
    dMatrixT		I6p_6;

    dMatrixT		I7p_1;
    dMatrixT		I7p_2;
    dMatrixT		I7p_3;
    dMatrixT		I7p_4;
    dMatrixT		I7p_5;
    dMatrixT		I7p_6;
    dMatrixT		I7p_7;
    dMatrixT		I7p_8;
    dMatrixT		I7p_9;
    dMatrixT		I7p_10;
    dMatrixT		I7p_11;
    dMatrixT		I7p_12;
    dMatrixT		I7p_13;


    dMatrixT		I8p_1;
    dMatrixT		I8p_2;
    dMatrixT		I8p_3;
    dMatrixT		I8p_4;
    dMatrixT		I8p_5;
    dMatrixT		I8p_6;
    dMatrixT		I8p_7;
    dMatrixT		I8p_8;
    dMatrixT		I8p_9;
    dMatrixT		I8p_10;
    dMatrixT		I8p_11;
    dMatrixT		I8p_12;
    dMatrixT		I8p_13;

    dMatrixT		I9p_1;
    dMatrixT		I9p_2;
    dMatrixT		I9p_3;
    dMatrixT		I9p_4;
    dMatrixT		I9p_5;
    dMatrixT		I9p_6;
    dMatrixT		I9p_7;
    dMatrixT		I9p_8;
    dMatrixT		I9p_9;
    dMatrixT		I9p_10;
    dMatrixT		I9p_11;
    dMatrixT		I9p_12;
    dMatrixT		I9p_13;

    dMatrixT		I10p_1;
    dMatrixT		I10p_2;
    dMatrixT		I10p_3;
    dMatrixT		I10p_4;
    dMatrixT		I10p_5;
    dMatrixT		I10p_6;
    dMatrixT		I10p_7;
    dMatrixT		I10p_8;
    dMatrixT		I10p_9;
    dMatrixT		I10p_10;
    dMatrixT		I10p_11;
    dMatrixT		I10p_12;
    dMatrixT		I10p_13;


    dMatrixT		I_temp_11p_1;
    dMatrixT		I_temp_11p_2;
    dMatrixT		I_temp_11p_3;
    dMatrixT		I_temp_11p_4;
    dMatrixT		I_temp_11p_5;
    dMatrixT		I_temp_11p_6;
    dMatrixT		I_temp_11p_7;
    dMatrixT		I_temp_11p_8;
    dMatrixT		I_temp_11p_9;
    dMatrixT		I_temp_11p_10;
    dMatrixT		I_temp_11p_11;
    dMatrixT		I_temp_11p_12;
    dMatrixT		I_temp_11p_13;


    dMatrixT		I11p_1;
    dMatrixT		I11p_2;

    dMatrixT		I12p_1;
    dMatrixT		I12p_2;
    dMatrixT		I12p_3;
    dMatrixT		I12p_4;
    dMatrixT		I12p_5;
    dMatrixT		I12p_6;

    dMatrixT		I13p_1;
    dMatrixT		I13p_2;
    dMatrixT		I13p_3;
    dMatrixT		I13p_4;
    dMatrixT		I13p_5;
    dMatrixT		I13p_6;




    dMatrixT		fKu_IJe_1;
    dMatrixT		fKu_IJe_2;
    dMatrixT		fKu_IJe_3;
    dMatrixT		fKu_IJe_4;
    dMatrixT		fKu_IJe_5;
    dMatrixT		fKu_IJe_6;
    dMatrixT		fKu_IJe_7;
    dMatrixT		fKu_IJe_8;



    dMatrixT		fKu_I1p_1;
    dMatrixT		fKu_I1p_2;
    dMatrixT		fKu_I1p_3;
    dMatrixT		fKu_I1p_4;
    dMatrixT		fKu_I1p_5;
    dMatrixT		fKu_I1p_6;



    dMatrixT		fKu_I2p_1;
    dMatrixT		fKu_I2p_2;
    dMatrixT		fKu_I2p_3;
    dMatrixT		fKu_I2p_4;
    dMatrixT		fKu_I2p_5;
    dMatrixT		fKu_I2p_6;



    dMatrixT		fKu_I3p_1;
    dMatrixT		fKu_I3p_2;
    dMatrixT		fKu_I3p_3;
    dMatrixT		fKu_I3p_4;
    dMatrixT		fKu_I3p_5;
    dMatrixT		fKu_I3p_6;



    dMatrixT		fKu_I4p_1;
    dMatrixT		fKu_I4p_2;
    dMatrixT		fKu_I4p_3;
    dMatrixT		fKu_I4p_4;
    dMatrixT		fKu_I4p_5;
    dMatrixT		fKu_I4p_6;



    dMatrixT		fKu_I5p_1;
    dMatrixT		fKu_I5p_2;
    dMatrixT		fKu_I5p_3;
    dMatrixT		fKu_I5p_4;
    dMatrixT		fKu_I5p_5;
    dMatrixT		fKu_I5p_6;



    dMatrixT		fKu_I6p_1;
    dMatrixT		fKu_I6p_2;
    dMatrixT		fKu_I6p_3;
    dMatrixT		fKu_I6p_4;
    dMatrixT		fKu_I6p_5;
    dMatrixT		fKu_I6p_6;

    dMatrixT		fKu_I7p_1;
    dMatrixT		fKu_I7p_2;
    dMatrixT		fKu_I7p_3;
    dMatrixT		fKu_I7p_4;
    dMatrixT		fKu_I7p_5;
    dMatrixT		fKu_I7p_6;
    dMatrixT		fKu_I7p_7;
    dMatrixT		fKu_I7p_8;
    dMatrixT		fKu_I7p_9;
    dMatrixT		fKu_I7p_10;
    dMatrixT		fKu_I7p_11;
    dMatrixT		fKu_I7p_12;
    dMatrixT		fKu_I7p_13;

    dMatrixT		fKu_I8p_1;
    dMatrixT		fKu_I8p_2;
    dMatrixT		fKu_I8p_3;
    dMatrixT		fKu_I8p_4;
    dMatrixT		fKu_I8p_5;
    dMatrixT		fKu_I8p_6;
    dMatrixT		fKu_I8p_7;
    dMatrixT		fKu_I8p_8;
    dMatrixT		fKu_I8p_9;
    dMatrixT		fKu_I8p_10;
    dMatrixT		fKu_I8p_11;
    dMatrixT		fKu_I8p_12;
    dMatrixT		fKu_I8p_13;

    dMatrixT		fKu_I9p_1;
    dMatrixT		fKu_I9p_2;
    dMatrixT		fKu_I9p_3;
    dMatrixT		fKu_I9p_4;
    dMatrixT		fKu_I9p_5;
    dMatrixT		fKu_I9p_6;
    dMatrixT		fKu_I9p_7;
    dMatrixT		fKu_I9p_8;
    dMatrixT		fKu_I9p_9;
    dMatrixT		fKu_I9p_10;
    dMatrixT		fKu_I9p_11;
    dMatrixT		fKu_I9p_12;
    dMatrixT		fKu_I9p_13;

    dMatrixT		fKu_I10p_1;
    dMatrixT		fKu_I10p_2;
    dMatrixT		fKu_I10p_3;
    dMatrixT		fKu_I10p_4;
    dMatrixT		fKu_I10p_5;
    dMatrixT		fKu_I10p_6;
    dMatrixT		fKu_I10p_7;
    dMatrixT		fKu_I10p_8;
    dMatrixT		fKu_I10p_9;
    dMatrixT		fKu_I10p_10;
    dMatrixT		fKu_I10p_11;
    dMatrixT		fKu_I10p_12;
    dMatrixT		fKu_I10p_13;

    dMatrixT		fKu_I12p_1;
    dMatrixT		fKu_I12p_2;
    dMatrixT		fKu_I12p_3;
    dMatrixT		fKu_I12p_4;
    dMatrixT		fKu_I12p_5;
    dMatrixT		fKu_I12p_6;

    dMatrixT		fKu_I13p_1;
    dMatrixT		fKu_I13p_2;
    dMatrixT		fKu_I13p_3;
    dMatrixT		fKu_I13p_4;
    dMatrixT		fKu_I13p_5;
    dMatrixT		fKu_I13p_6;

    /////////////////////////////Test//////////////////////////////
    dMatrixT		I_temp_11p_test_1;
    dMatrixT		I_temp_11p_test_2;
    dMatrixT		I_temp_11p_test_3;
    dMatrixT		I_temp_11p_test_4;
    dMatrixT		I12p_test_2;
    dMatrixT		I12p_test_3;
    dMatrixT		I12p_test_4;
    dMatrixT		I_temp_11p_20;
    dMatrixT		I_temp_11p_21;
    dMatrixT		I_temp_11p_22;
    dMatrixT		I_temp_11p_23;
    dMatrixT		I_temp_11p_24;

    dMatrixT		I_temp_11p_4_Transpose;

    void Form_I_temp_11p_test_1(void);
    void Form_I_temp_11p_test_2(void);
    void Form_I_temp_11p_test_3(void);
    void Form_I_temp_11p_test_4(void);
    void Form_I12p_test_2(void);


    void Form_I_temp_11p_test_2_1_1(void);
    void Form_I_temp_11p_test_2_1_2(void);
    void Form_I_temp_11p_test_2_1_3(void);
    void Form_I_temp_11p_test_2_1_4(void);

    void Form_I_temp_11p_test_2_2_1(void);
    void Form_I_temp_11p_test_2_2_2(void);
    void Form_I_temp_11p_test_2_2_3(void);
    void Form_I_temp_11p_test_2_2_4(void);

    void Form_I_temp_11p_test_2_3_1(void);
    void Form_I_temp_11p_test_2_3_2(void);


    void Form_I_temp_11p_test_2_4_1(void);
    void Form_I_temp_11p_test_2_4_2(void);
    void Form_I_temp_11p_test_2_4_3(void);
    void Form_I_temp_11p_test_2_4_4(void);

    void Form_I_temp_11p_test_2_5_1(void);
    void Form_I_temp_11p_test_2_5_2(void);


    void Form_I_temp_11p_test_2_6_1(void);
    void Form_I_temp_11p_test_2_6_2(void);
    void Form_I_temp_11p_test_2_6_3(void);
    void Form_I_temp_11p_test_2_6_4(void);

    void Form_I_temp_11p_test_2_7_1(void);
    void Form_I_temp_11p_test_2_7_2(void);


    void Form_I_temp_11p_test_2_8_1(void);
    void Form_I_temp_11p_test_2_8_2(void);
    void Form_I_temp_11p_test_2_8_3(void);
    void Form_I_temp_11p_test_2_8_4(void);

    void Form_I_temp_11p_test_2_9_1(void);
    void Form_I_temp_11p_test_2_9_2(void);
    void Form_I_temp_11p_test_2_9_3(void);
    void Form_I_temp_11p_test_2_9_4(void);

    void Form_I_temp_11p_test_2_10_1(void);
    void Form_I_temp_11p_test_2_10_2(void);
    void Form_I_temp_11p_test_2_10_3(void);
    void Form_I_temp_11p_test_2_10_4(void);

    void Form_I_temp_11p_test_2_11_1(void);
    void Form_I_temp_11p_test_2_11_2(void);
    void Form_I_temp_11p_test_2_11_3(void);
    void Form_I_temp_11p_test_2_11_4(void);

    void Form_I_temp_11p_test_2_12_1(void);
    void Form_I_temp_11p_test_2_12_2(void);
    void Form_I_temp_11p_test_2_12_3(void);
    void Form_I_temp_11p_test_2_12_4(void);

    void Form_I_temp_11p_test_2_13_1(void);
    void Form_I_temp_11p_test_2_13_2(void);
    void Form_I_temp_11p_test_2_13_3(void);
    void Form_I_temp_11p_test_2_13_4(void);


    dMatrixT		I_temp_11p_test_2_1_1;
    dMatrixT		I_temp_11p_test_2_1_2;
    dMatrixT		I_temp_11p_test_2_1_3;
    dMatrixT		I_temp_11p_test_2_1_4;

    dMatrixT		I_temp_11p_test_2_2_1;
    dMatrixT		I_temp_11p_test_2_2_2;
    dMatrixT		I_temp_11p_test_2_2_3;
    dMatrixT		I_temp_11p_test_2_2_4;

    dMatrixT		I_temp_11p_test_2_3_1;
    dMatrixT		I_temp_11p_test_2_3_2;


    dMatrixT		I_temp_11p_test_2_4_1;
    dMatrixT		I_temp_11p_test_2_4_2;
    dMatrixT		I_temp_11p_test_2_4_3;
    dMatrixT		I_temp_11p_test_2_4_4;

    dMatrixT		I_temp_11p_test_2_5_1;
    dMatrixT		I_temp_11p_test_2_5_2;


    dMatrixT		I_temp_11p_test_2_6_1;
    dMatrixT		I_temp_11p_test_2_6_2;
    dMatrixT		I_temp_11p_test_2_6_3;
    dMatrixT		I_temp_11p_test_2_6_4;

    dMatrixT		I_temp_11p_test_2_7_1;
    dMatrixT		I_temp_11p_test_2_7_2;


    dMatrixT		I_temp_11p_test_2_8_1;
    dMatrixT		I_temp_11p_test_2_8_2;
    dMatrixT		I_temp_11p_test_2_8_3;
    dMatrixT		I_temp_11p_test_2_8_4;

    dMatrixT		I_temp_11p_test_2_9_1;
    dMatrixT		I_temp_11p_test_2_9_2;
    dMatrixT		I_temp_11p_test_2_9_3;
    dMatrixT		I_temp_11p_test_2_9_4;

    dMatrixT		I_temp_11p_test_2_10_1;
    dMatrixT		I_temp_11p_test_2_10_2;
    dMatrixT		I_temp_11p_test_2_10_3;
    dMatrixT		I_temp_11p_test_2_10_4;

    dMatrixT		I_temp_11p_test_2_11_1;
    dMatrixT		I_temp_11p_test_2_11_2;
    dMatrixT		I_temp_11p_test_2_11_3;
    dMatrixT		I_temp_11p_test_2_11_4;

    dMatrixT		I_temp_11p_test_2_12_1;
    dMatrixT		I_temp_11p_test_2_12_2;
    dMatrixT		I_temp_11p_test_2_12_3;
    dMatrixT		I_temp_11p_test_2_12_4;

    dMatrixT		I_temp_11p_test_2_13_1;
    dMatrixT		I_temp_11p_test_2_13_2;
    dMatrixT		I_temp_11p_test_2_13_3;
    dMatrixT		I_temp_11p_test_2_13_4;

    dMatrixT		fKu_I_temp_11p_test_2_1_1;
    dMatrixT		fKu_I_temp_11p_test_2_1_2;
    dMatrixT		fKu_I_temp_11p_test_2_1_3;
    dMatrixT		fKu_I_temp_11p_test_2_1_4;

    dMatrixT		fKu_I_temp_11p_test_2_2_1;
    dMatrixT		fKu_I_temp_11p_test_2_2_2;
    dMatrixT		fKu_I_temp_11p_test_2_2_3;
    dMatrixT		fKu_I_temp_11p_test_2_2_4;

    dMatrixT		fKu_I_temp_11p_test_2_3_1;
    dMatrixT		fKu_I_temp_11p_test_2_3_2;


    dMatrixT		fKu_I_temp_11p_test_2_4_1;
    dMatrixT		fKu_I_temp_11p_test_2_4_2;
    dMatrixT		fKu_I_temp_11p_test_2_4_3;
    dMatrixT		fKu_I_temp_11p_test_2_4_4;

    dMatrixT		fKu_I_temp_11p_test_2_5_1;
    dMatrixT		fKu_I_temp_11p_test_2_5_2;


    dMatrixT		fKu_I_temp_11p_test_2_6_1;
    dMatrixT		fKu_I_temp_11p_test_2_6_2;
    dMatrixT		fKu_I_temp_11p_test_2_6_3;
    dMatrixT		fKu_I_temp_11p_test_2_6_4;

    dMatrixT		fKu_I_temp_11p_test_2_7_1;
    dMatrixT		fKu_I_temp_11p_test_2_7_2;


    dMatrixT		fKu_I_temp_11p_test_2_8_1;
    dMatrixT		fKu_I_temp_11p_test_2_8_2;
    dMatrixT		fKu_I_temp_11p_test_2_8_3;
    dMatrixT		fKu_I_temp_11p_test_2_8_4;

    dMatrixT		fKu_I_temp_11p_test_2_9_1;
    dMatrixT		fKu_I_temp_11p_test_2_9_2;
    dMatrixT		fKu_I_temp_11p_test_2_9_3;
    dMatrixT		fKu_I_temp_11p_test_2_9_4;

    dMatrixT		fKu_I_temp_11p_test_2_10_1;
    dMatrixT		fKu_I_temp_11p_test_2_10_2;
    dMatrixT		fKu_I_temp_11p_test_2_10_3;
    dMatrixT		fKu_I_temp_11p_test_2_10_4;

    dMatrixT		fKu_I_temp_11p_test_2_11_1;
    dMatrixT		fKu_I_temp_11p_test_2_11_2;
    dMatrixT		fKu_I_temp_11p_test_2_11_3;
    dMatrixT		fKu_I_temp_11p_test_2_11_4;

    dMatrixT		fKu_I_temp_11p_test_2_12_1;
    dMatrixT		fKu_I_temp_11p_test_2_12_2;
    dMatrixT		fKu_I_temp_11p_test_2_12_3;
    dMatrixT		fKu_I_temp_11p_test_2_12_4;

    dMatrixT		fKu_I_temp_11p_test_2_13_1;
    dMatrixT		fKu_I_temp_11p_test_2_13_2;
    dMatrixT		fKu_I_temp_11p_test_2_13_3;
    dMatrixT		fKu_I_temp_11p_test_2_13_4;

    dMatrixT		fKdd_previous;
    dMatrixT		fKdd_full_implemet;
    dMatrixT		difference;

//////////////////////////////////////////////////////////////////////////


    dMatrixT		II_temp_11p_1_1;
    dMatrixT		II_temp_11p_1_2;
    dMatrixT		II_temp_11p_1_3;
    dMatrixT		II_temp_11p_1_4;
    dMatrixT		II_temp_11p_1_5;
    dMatrixT		II_temp_11p_1_6;

    dMatrixT		II_temp_11p_2_1;
    dMatrixT		II_temp_11p_2_2;
    dMatrixT		II_temp_11p_2_3;
    dMatrixT		II_temp_11p_2_4;
    dMatrixT		II_temp_11p_2_5;
    dMatrixT		II_temp_11p_2_6;

    dMatrixT		II_temp_11p_3_1;
    dMatrixT		II_temp_11p_3_2;
    dMatrixT		II_temp_11p_3_3;
    dMatrixT		II_temp_11p_3_4;
    dMatrixT		II_temp_11p_3_5;
    dMatrixT		II_temp_11p_3_6;

    dMatrixT		II_temp_11p_4_1;
    dMatrixT		II_temp_11p_4_2;
    dMatrixT		II_temp_11p_4_3;
    dMatrixT		II_temp_11p_4_4;
    dMatrixT		II_temp_11p_4_5;
    dMatrixT		II_temp_11p_4_6;

    dMatrixT		II_temp_11p_5_1;
    dMatrixT		II_temp_11p_5_2;
    dMatrixT		II_temp_11p_5_3;
    dMatrixT		II_temp_11p_5_4;
    dMatrixT		II_temp_11p_5_5;
    dMatrixT		II_temp_11p_5_6;

    dMatrixT		II_temp_11p_6_1;
    dMatrixT		II_temp_11p_6_2;
    dMatrixT		II_temp_11p_6_3;
    dMatrixT		II_temp_11p_6_4;
    dMatrixT		II_temp_11p_6_5;
    dMatrixT		II_temp_11p_6_6;

    dMatrixT		II_temp_11p_7_1;
    dMatrixT		II_temp_11p_7_2;
    dMatrixT		II_temp_11p_7_3;
    dMatrixT		II_temp_11p_7_4;
    dMatrixT		II_temp_11p_7_5;
    dMatrixT		II_temp_11p_7_6;

    dMatrixT		II_temp_11p_8_1;
    dMatrixT		II_temp_11p_8_2;
    dMatrixT		II_temp_11p_8_3;
    dMatrixT		II_temp_11p_8_4;
    dMatrixT		II_temp_11p_8_5;
    dMatrixT		II_temp_11p_8_6;

    dMatrixT		II_temp_11p_9_1;
    dMatrixT		II_temp_11p_9_2;
    dMatrixT		II_temp_11p_9_3;
    dMatrixT		II_temp_11p_9_4;
    dMatrixT		II_temp_11p_9_5;
    dMatrixT		II_temp_11p_9_6;

    dMatrixT		II_temp_11p_10_1;
    dMatrixT		II_temp_11p_10_2;
    dMatrixT		II_temp_11p_10_3;
    dMatrixT		II_temp_11p_10_4;
    dMatrixT		II_temp_11p_10_5;
    dMatrixT		II_temp_11p_10_6;

    dMatrixT		II_temp_11p_11_1;
    dMatrixT		II_temp_11p_11_2;
    dMatrixT		II_temp_11p_11_3;
    dMatrixT		II_temp_11p_11_4;
    dMatrixT		II_temp_11p_11_5;
    dMatrixT		II_temp_11p_11_6;

    dMatrixT		II_temp_11p_12_1;
    dMatrixT		II_temp_11p_12_2;
    dMatrixT		II_temp_11p_12_3;
    dMatrixT		II_temp_11p_12_4;
    dMatrixT		II_temp_11p_12_5;
    dMatrixT		II_temp_11p_12_6;

    dMatrixT		II_temp_11p_13_1;
    dMatrixT		II_temp_11p_13_2;
    dMatrixT		II_temp_11p_13_3;
    dMatrixT		II_temp_11p_13_4;
    dMatrixT		II_temp_11p_13_5;
    dMatrixT		II_temp_11p_13_6;


    dMatrixT		fKu_I_temp_11p_1_1;
    dMatrixT		fKu_I_temp_11p_1_2;
    dMatrixT		fKu_I_temp_11p_1_3;
    dMatrixT		fKu_I_temp_11p_1_4;
    dMatrixT		fKu_I_temp_11p_1_5;
    dMatrixT		fKu_I_temp_11p_1_6;

    dMatrixT		fKu_I_temp_11p_2_1;
    dMatrixT		fKu_I_temp_11p_2_2;
    dMatrixT		fKu_I_temp_11p_2_3;
    dMatrixT		fKu_I_temp_11p_2_4;
    dMatrixT		fKu_I_temp_11p_2_5;
    dMatrixT		fKu_I_temp_11p_2_6;

    dMatrixT		fKu_I_temp_11p_3_1;
    dMatrixT		fKu_I_temp_11p_3_2;
    dMatrixT		fKu_I_temp_11p_3_3;
    dMatrixT		fKu_I_temp_11p_3_4;
    dMatrixT		fKu_I_temp_11p_3_5;
    dMatrixT		fKu_I_temp_11p_3_6;

    dMatrixT		fKu_I_temp_11p_4_1;
    dMatrixT		fKu_I_temp_11p_4_2;
    dMatrixT		fKu_I_temp_11p_4_3;
    dMatrixT		fKu_I_temp_11p_4_4;
    dMatrixT		fKu_I_temp_11p_4_5;
    dMatrixT		fKu_I_temp_11p_4_6;

    dMatrixT		fKu_I_temp_11p_5_1;
    dMatrixT		fKu_I_temp_11p_5_2;
    dMatrixT		fKu_I_temp_11p_5_3;
    dMatrixT		fKu_I_temp_11p_5_4;
    dMatrixT		fKu_I_temp_11p_5_5;
    dMatrixT		fKu_I_temp_11p_5_6;

    dMatrixT		fKu_I_temp_11p_6_1;
    dMatrixT		fKu_I_temp_11p_6_2;
    dMatrixT		fKu_I_temp_11p_6_3;
    dMatrixT		fKu_I_temp_11p_6_4;
    dMatrixT		fKu_I_temp_11p_6_5;
    dMatrixT		fKu_I_temp_11p_6_6;

    dMatrixT		fKu_I_temp_11p_7_1;
    dMatrixT		fKu_I_temp_11p_7_2;
    dMatrixT		fKu_I_temp_11p_7_3;
    dMatrixT		fKu_I_temp_11p_7_4;
    dMatrixT		fKu_I_temp_11p_7_5;
    dMatrixT		fKu_I_temp_11p_7_6;

    dMatrixT		fKu_I_temp_11p_8_1;
    dMatrixT		fKu_I_temp_11p_8_2;
    dMatrixT		fKu_I_temp_11p_8_3;
    dMatrixT		fKu_I_temp_11p_8_4;
    dMatrixT		fKu_I_temp_11p_8_5;
    dMatrixT		fKu_I_temp_11p_8_6;

    dMatrixT		fKu_I_temp_11p_9_1;
    dMatrixT		fKu_I_temp_11p_9_2;
    dMatrixT		fKu_I_temp_11p_9_3;
    dMatrixT		fKu_I_temp_11p_9_4;
    dMatrixT		fKu_I_temp_11p_9_5;
    dMatrixT		fKu_I_temp_11p_9_6;

    dMatrixT		fKu_I_temp_11p_10_1;
    dMatrixT		fKu_I_temp_11p_10_2;
    dMatrixT		fKu_I_temp_11p_10_3;
    dMatrixT		fKu_I_temp_11p_10_4;
    dMatrixT		fKu_I_temp_11p_10_5;
    dMatrixT		fKu_I_temp_11p_10_6;

    dMatrixT		fKu_I_temp_11p_11_1;
    dMatrixT		fKu_I_temp_11p_11_2;
    dMatrixT		fKu_I_temp_11p_11_3;
    dMatrixT		fKu_I_temp_11p_11_4;
    dMatrixT		fKu_I_temp_11p_11_5;
    dMatrixT		fKu_I_temp_11p_11_6;

    dMatrixT		fKu_I_temp_11p_12_1;
    dMatrixT		fKu_I_temp_11p_12_2;
    dMatrixT		fKu_I_temp_11p_12_3;
    dMatrixT		fKu_I_temp_11p_12_4;
    dMatrixT		fKu_I_temp_11p_12_5;
    dMatrixT		fKu_I_temp_11p_12_6;

    dMatrixT		fKu_I_temp_11p_13_1;
    dMatrixT		fKu_I_temp_11p_13_2;
    dMatrixT		fKu_I_temp_11p_13_3;
    dMatrixT		fKu_I_temp_11p_13_4;
    dMatrixT		fKu_I_temp_11p_13_5;
    dMatrixT		fKu_I_temp_11p_13_6;

//////////////////////////////////////////////////////////////////////

    void Form_II_temp_11p_1_1(void);
    void Form_II_temp_11p_1_2(void);
    void Form_II_temp_11p_1_3(void);
    void Form_II_temp_11p_1_4(void);
    void Form_II_temp_11p_1_5(void);
    void Form_II_temp_11p_1_6(void);

    void Form_II_temp_11p_2_1(void);
    void Form_II_temp_11p_2_2(void);
    void Form_II_temp_11p_2_3(void);
    void Form_II_temp_11p_2_4(void);
    void Form_II_temp_11p_2_5(void);
    void Form_II_temp_11p_2_6(void);

    void Form_II_temp_11p_3_1(void);
    void Form_II_temp_11p_3_2(void);
    void Form_II_temp_11p_3_3(void);
    void Form_II_temp_11p_3_4(void);
    void Form_II_temp_11p_3_5(void);
    void Form_II_temp_11p_3_6(void);

    void Form_II_temp_11p_4_1(void);
    void Form_II_temp_11p_4_2(void);
    void Form_II_temp_11p_4_3(void);
    void Form_II_temp_11p_4_4(void);
    void Form_II_temp_11p_4_5(void);
    void Form_II_temp_11p_4_6(void);

    void Form_II_temp_11p_5_1(void);
    void Form_II_temp_11p_5_2(void);
    void Form_II_temp_11p_5_3(void);
    void Form_II_temp_11p_5_4(void);
    void Form_II_temp_11p_5_5(void);
    void Form_II_temp_11p_5_6(void);

    void Form_II_temp_11p_6_1(void);
    void Form_II_temp_11p_6_2(void);
    void Form_II_temp_11p_6_3(void);
    void Form_II_temp_11p_6_4(void);
    void Form_II_temp_11p_6_5(void);
    void Form_II_temp_11p_6_6(void);

    void Form_II_temp_11p_7_1(void);
    void Form_II_temp_11p_7_2(void);
    void Form_II_temp_11p_7_3(void);
    void Form_II_temp_11p_7_4(void);
    void Form_II_temp_11p_7_5(void);
    void Form_II_temp_11p_7_6(void);

    void Form_II_temp_11p_8_1(void);
    void Form_II_temp_11p_8_2(void);
    void Form_II_temp_11p_8_3(void);
    void Form_II_temp_11p_8_4(void);
    void Form_II_temp_11p_8_5(void);
    void Form_II_temp_11p_8_6(void);

    void Form_II_temp_11p_9_1(void);
    void Form_II_temp_11p_9_2(void);
    void Form_II_temp_11p_9_3(void);
    void Form_II_temp_11p_9_4(void);
    void Form_II_temp_11p_9_5(void);
    void Form_II_temp_11p_9_6(void);

    void Form_II_temp_11p_10_1(void);
    void Form_II_temp_11p_10_2(void);
    void Form_II_temp_11p_10_3(void);
    void Form_II_temp_11p_10_4(void);
    void Form_II_temp_11p_10_5(void);
    void Form_II_temp_11p_10_6(void);

    void Form_II_temp_11p_11_1(void);
    void Form_II_temp_11p_11_2(void);
    void Form_II_temp_11p_11_3(void);
    void Form_II_temp_11p_11_4(void);
    void Form_II_temp_11p_11_5(void);
    void Form_II_temp_11p_11_6(void);

    void Form_II_temp_11p_12_1(void);
    void Form_II_temp_11p_12_2(void);
    void Form_II_temp_11p_12_3(void);
    void Form_II_temp_11p_12_4(void);
    void Form_II_temp_11p_12_5(void);
    void Form_II_temp_11p_12_6(void);


    void Form_II_temp_11p_13_1(void);
    void Form_II_temp_11p_13_2(void);
    void Form_II_temp_11p_13_3(void);
    void Form_II_temp_11p_13_4(void);
    void Form_II_temp_11p_13_5(void);
    void Form_II_temp_11p_13_6(void);


    /////////////////////////////////////////////////











    void  Form_fV1p(void);
    void  Form_fV1e(void);

    void Form_IJe_1(void);
    void Form_IJe_2(void);
    void Form_IJe_3(void);
    void Form_IJe_4(void);
    void Form_IJe_5(void);
    void Form_IJe_6(void);
    void Form_IJe_7(void);
    void Form_IJe_8(void);


    void Form_I1p_1(void);
    void Form_I1p_2(void);
    void Form_I1p_3(void);
    void Form_I1p_4(void);
    void Form_I1p_5(void);
    void Form_I1p_6(void);


    void Form_I2p_1(void);
    void Form_I2p_2(void);
    void Form_I2p_3(void);
    void Form_I2p_4(void);
    void Form_I2p_5(void);
    void Form_I2p_6(void);


    void Form_I3p_1(void);
    void Form_I3p_2(void);
    void Form_I3p_3(void);
    void Form_I3p_4(void);
    void Form_I3p_5(void);
    void Form_I3p_6(void);


    void Form_I4p_1(void);
    void Form_I4p_2(void);
    void Form_I4p_3(void);
    void Form_I4p_4(void);
    void Form_I4p_5(void);
    void Form_I4p_6(void);


    void Form_I5p_1(void);
    void Form_I5p_2(void);
    void Form_I5p_3(void);
    void Form_I5p_4(void);
    void Form_I5p_5(void);
    void Form_I5p_6(void);


    void Form_I6p_1(void);
    void Form_I6p_2(void);
    void Form_I6p_3(void);
    void Form_I6p_4(void);
    void Form_I6p_5(void);
    void Form_I6p_6(void);


    void Form_I7p_1(void);
    void Form_I7p_2(void);
    void Form_I7p_3(void);
    void Form_I7p_4(void);
    void Form_I7p_5(void);
    void Form_I7p_6(void);
    void Form_I7p_7(void);
    void Form_I7p_8(void);
    void Form_I7p_9(void);
    void Form_I7p_10(void);
    void Form_I7p_11(void);
    void Form_I7p_12(void);
    void Form_I7p_13(void);


    void Form_I8p_1(void);
    void Form_I8p_2(void);
    void Form_I8p_3(void);
    void Form_I8p_4(void);
    void Form_I8p_5(void);
    void Form_I8p_6(void);
    void Form_I8p_7(void);
    void Form_I8p_8(void);
    void Form_I8p_9(void);
    void Form_I8p_10(void);
    void Form_I8p_11(void);
    void Form_I8p_12(void);
    void Form_I8p_13(void);


    void Form_I9p_1(void);
    void Form_I9p_2(void);
    void Form_I9p_3(void);
    void Form_I9p_4(void);
    void Form_I9p_5(void);
    void Form_I9p_6(void);
    void Form_I9p_7(void);
    void Form_I9p_8(void);
    void Form_I9p_9(void);
    void Form_I9p_10(void);
    void Form_I9p_11(void);
    void Form_I9p_12(void);
    void Form_I9p_13(void);


    void Form_I10p_1(void);
    void Form_I10p_2(void);
    void Form_I10p_3(void);
    void Form_I10p_4(void);
    void Form_I10p_5(void);
    void Form_I10p_6(void);
    void Form_I10p_7(void);
    void Form_I10p_8(void);
    void Form_I10p_9(void);
    void Form_I10p_10(void);
    void Form_I10p_11(void);
    void Form_I10p_12(void);
    void Form_I10p_13(void);


    void Form_I_temp_11p_1(void);
    void Form_I_temp_11p_2(void);
    void Form_I_temp_11p_3(void);
    void Form_I_temp_11p_4(void);
    void Form_I_temp_11p_5(void);
    void Form_I_temp_11p_6(void);
    void Form_I_temp_11p_7(void);
    void Form_I_temp_11p_8(void);
    void Form_I_temp_11p_9(void);
    void Form_I_temp_11p_10(void);
    void Form_I_temp_11p_11(void);
    void Form_I_temp_11p_12(void);
    void Form_I_temp_11p_13(void);


    void Form_I12p_1(void);
    void Form_I12p_2(void);
    void Form_I12p_3(void);
    void Form_I12p_4(void);
    void Form_I12p_5(void);
    void Form_I12p_6(void);

    void Form_I13p_1(void);
    void Form_I13p_2(void);
    void Form_I13p_3(void);
    void Form_I13p_4(void);
    void Form_I13p_5(void);
    void Form_I13p_6(void);







    ///////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////



    //////////////////////////////////////////////////////////
    /////DEFINITIONS FOR MICROMORPHIC MATRICES////////////////
    //////////////////////////////////////////////////////////
    double KrDelta[3][3];
    double trdeltad;
    double trdeltaEp;

    dMatrixT PHIMATRIX;
    dMatrixT GPHIMATRIX;
    dTensor3DT GRAD_CHIM;


    double Counter;
    dArray2DT Counter_IPs_el_n;
    dArray2DT Counter_IPs_el;
    dArrayT Counter_IPs;
    dMatrixT fIota_temp_matrix;

    //Varitional Matrices coming from the Balance of linear Momentum
    dMatrixT Tsigma_1;
    dMatrixT fG1_1;
    double SigN[3][3];//unsymmetric Cauchy stress tensor found at previous step
    dMatrixT SigN_m;
    dMatrixT KirchhoffST;
    dMatrixT Temp_SPK;
    dMatrixT SPK;
    dMatrixT Fn_m;
    dMatrixT Finv_m;
    dMatrixT deltaL;
    dMatrixT deltaL_Tr;
    dMatrixT tempSig;
    dMatrixT deltad;
    dMatrixT SigN_ar;


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
    dMatrixT ChiM;
    dMatrixT ChiM_Inverse;
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

    dMatrixT eps;// Micro strain tensor in current config.
    dMatrixT psi;// micro-deformation tensor in current config.

    //plasticity parameters may be needed in the future

   double Apsi_chi;
   double dYieldTrialTol;
   int iPlasticityCheck;

    //////////////////////////////////////////////////////////
    //////FINITE STRAIN ELASTICITY MATRICES START HERE/////////
    //////////////////////////////////////////////////////////





    /***************************************************/
    /*****Micro-scale plasticity matrices **************/
    /***************************************************/



    /***************************************************/
    /************Micro Scale Gradient Plasticity *******/


    /***************************************************/
    /***************************************************/
    /***************************************************/



    dMatrixT fRight_Cauchy_Green_tensor_tr;
    dMatrixT fLagrangian_strain_tensor_tr;
    dMatrixT fMicroRight_Cauchy_Green_tensor;
    dMatrixT fMicroRight_Cauchy_Green_tensor_tr;

    //dMatrixT fMicroRight_Cauchy_Green_tensor;
    //dMatrixT fMicroRight_Cauchy_Green_tensor_tr;


    //dMatrixT fSIGMA_S;





    /* for local Newton-Raphson iteration */
   int iIterationMax;
   double dRelTol, dAbsTol,fConst1;

   int PlasticityCondition;





    dArray2DT   fState_variables_IPs;
    dArray2DT   fState_variables_Elements_IPs;
    dArray2DT   fState_variables_n_IPs;
    dArray2DT   fState_variables_n_Elements_IPs;




/////stress invariants variables////////
 double Cauchy_inv;
 double trS;//trace of SPK
 double invdevS;// Norm of dev. part of SPK
 double Rel_strs_inv;
 double trSIGMA_S;
 double invdevSIGMA_S;// Norm of dev. part of SIGMA_S
 double Higher_orderT_inv;
 double invtrM;//||tr(M)|| trM is a vector
 double invdevMKLM;//||dev(M)||
 double temp_inv;
 double invPhi;
 double invGPhi;
 double deveps; // ||deveps||
 double treps;//trace of epsilon
 double invtrgammastn;
 double invdevgammastn;


 dMatrixT devsigma;
 dMatrixT devRelsts;
 dTensor3DT  devmklm;
 dMatrixT s_sigma_temp;
 dTensor3DT fmklm;
 dTensor3DT gammastn;
 dTensor3DT devgammastn;
 dArrayT trgammastn;


//////////////////////////////////////////////////
    double invJ;

    double trsigma;
    double trs_sigma;
    double trmklm;
    double invtrMKLM;
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
        /* to store fCauchy_stress_IPs for each of the 27 IPs of each element */

        ////////////////////////////////////
        dArray2DT  fE_values_IPs;
        dArray2DT  fVarepsilon_IPs;
        ////////////////////////////////////

        dArray2DT   fEulerian_strain_Elements_IPs;


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
    void Form_micro_shape_functions(const double* &shapes_micro_X);
    void Form_deformation_gradient_tensor(void);
    void Form_Grad_grad_transformation_matrix(void);
    void Form_deformation_gradient_inv_vector(void);
    void Extract_six_values_from_symmetric_tensor(const dMatrixT &fTensor,dArrayT& fTemp_six_values);
    void Put_values_In_dArrayT_vector(const dArray2DT &f2DArrayT,const int& e,const int& IP,dArrayT& fArrayT);
    void Put_values_In_Array(const dArray2DT &f2DArrayT,const int& e,const int& IP,dArrayT& fArrayT);

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
    void Form_ChiM(void);
    void Form_Second_Piola_Kirchhoff_SPK(const dMatrixT& LagStn, const dMatrixT& MicroStn);


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



    ////////////////////////////////////////////////////////
    /////////// functions to calculate stress measures ////
    void Calculate_Cauchy_INV(void);
    void Calculate_stress_diff_INV(void);
    void Calculate_higher_order_tensor_INV(void);
    void Calculate_fmklm(void);
    void Calculate_PHI_GPHI_matrices_INV(void);
    void Caculate_invdevpart_of_Matrix(const dMatrixT &fMatrix,dMatrixT &fdevfMatrix,double devinvariant);
    void Calculate_relative_strain_INV(void);
    void Calculate_HOST_INV(void);

/* Plasticity functions */



// The first and the second terms cancel each other, we start from the thirh term which has 12 matrices
// all the matrices having "e" are the ones which direclty include the term  d(deltau)/dX
// matrices having "p" are the one which include delta(gamma)



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


inline const ShapeFunctionT& FSMicromorphic3DCurrConfigT::ShapeFunctionDispl(void) const
{
#if __option(extended_errorcheck)
    if (!fShapes_displ)
        ExceptionT::GeneralFail("FSMicromorphic3DCurrConfigT::ShapeFunctionDispl", "no displ shape functions");
#endif
    return *fShapes_displ;
}

inline const ShapeFunctionT& FSMicromorphic3DCurrConfigT::ShapeFunctionMicro(void) const
{
#if __option(extended_errorcheck)
    if (!fShapes_micro)
        ExceptionT::GeneralFail("FSMicromorphic3DCurrConfigT::ShapeFunctionMicro", "no micro shape functions");
#endif
    return *fShapes_micro;
}

/* return the geometry code */
inline GeometryT::CodeT FSMicromorphic3DCurrConfigT::GeometryCode(void) const
{ return fGeometryCode_displ; }


} // namespace Tahoe
#endif /* _FS_MICROMORPHIC_3D_CURR_CONFIG_T_H_ */



