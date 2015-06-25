/*$Id: MR3DT.cpp,v 1.13 2011/12/01 20:38:01 beichuan Exp $*/
/* Elastolastic Cohesive Model for Geomaterials*/
#include "MR3DT.h"

#include <iostream>
#include <cmath>

#include "ExceptionT.h"
#include "ifstreamT.h"
#include "ofstreamT.h"
#include "dArrayT.h"
#include "dMatrixT.h"
#include "nMatrixT.h"
#include "ParameterUtils.h"

using namespace Tahoe;

/* class parameters */
const int knumDOF = 3;

/* constructor */
MR3DT::MR3DT(void): SurfacePotentialT(knumDOF)
{
	const char caller[] = "MR3DT::MR3DT";
	
	SetName("elastoplastic_MR_3D");
	
#if 0
	/* Elastic and Fracture Energy parameters */
	in >> fE_n; if (fE_n < 0) throw ExceptionT::kBadInputValue;
	in >> fE_t1; if (fE_t1 < 0) throw ExceptionT::kBadInputValue;
	in >> fE_t2; if (fE_t2 < 0) throw ExceptionT::kBadInputValue;
	in >> fGf_I; if (fGf_I < 0) throw ExceptionT::kBadInputValue;
	in >> fGf_II; if (fGf_II < 0) throw ExceptionT::kBadInputValue;

	/* Inelastic Response initiation parameters */
	in >> fchi_p; if (fchi_p < 0) throw ExceptionT::kBadInputValue;
	in >> fchi_r; if (fchi_r < 0) throw ExceptionT::kBadInputValue;
	in >> fc_p; if (fc_p < 0) throw ExceptionT::kBadInputValue;
	in >> fc_r; if (fc_r < 0) throw ExceptionT::kBadInputValue;
	in >> fphi_p; if (fphi_p < 0) throw ExceptionT::kBadInputValue;
	in >> fphi_r; if (fphi_r < 0) throw ExceptionT::kBadInputValue;
	in >> fpsi_p; if (fpsi_p < 0) throw ExceptionT::kBadInputValue;
	in >> falpha_chi; if (falpha_chi < 0) throw ExceptionT::kBadInputValue;
	in >> falpha_c; if (falpha_c < 0) throw ExceptionT::kBadInputValue;
	in >> falpha_phi; if (falpha_phi < 0) throw ExceptionT::kBadInputValue;
	in >> falpha_psi; if (falpha_psi < 0) throw ExceptionT::kBadInputValue;
	in >> fTol_1; if (fTol_1 < 0) throw ExceptionT::kBadInputValue;
	in >> fTol_2; if (fTol_2 < 0) throw ExceptionT::kBadInputValue;
	
	/*double esp = 0.;
	double enp = 0.;
	fchi = fchi_r + (fchi_p - fchi_r)*exp(-falpha_chi*enp);
	fc   = fc_r + (fc_p - fc_r)*exp(-falpha_c*esp);
	fphi = tan(fphi_r) + (tan(fphi_p) - tan(fphi_r))*exp(-falpha_phi*esp);
	fpsi = tan(fpsi_p)*exp(-falpha_psi*esp);*/
#endif	
}

/* return the number of state variables needed by the model */
int MR3DT::NumStateVariables(void) const { return 23; }  //need to be checked

/* describe the parameters needed by the interface */
void MR3DT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SurfacePotentialT::DefineParameters(list);

	/* model parameters */
	ParameterT E_n(fE_n, "E_n");
	E_n.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(E_n);
	
	ParameterT E_t1(fE_t1, "E_t1");
	E_t1.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(E_t1);

	ParameterT E_t2(fE_t2, "E_t2");
	E_t2.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(E_t2);


	ParameterT Gf_I(fGf_I, "Gf_I");
	Gf_I.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(Gf_I);
	
	ParameterT Gf_II(fGf_II, "Gf_II");
	Gf_II.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(Gf_II);
	
	ParameterT chi_p(fchi_p, "chi_p");
	chi_p.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(chi_p);
	
	ParameterT chi_r(fchi_r, "chi_r");
	chi_r.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(chi_r);
	
	ParameterT c_p(fc_p, "c_p");
	c_p.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(c_p);
	
	ParameterT c_r(fc_r, "c_r");
	c_r.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(c_r);
	
	ParameterT phi_p(fphi_p, "phi_p");
	phi_p.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(phi_p);
	
	ParameterT phi_r(fphi_r, "phi_r");
	phi_r.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(phi_r);
	
	ParameterT psi_p(fpsi_p, "psi_p");
	psi_p.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(psi_p);

	ParameterT alpha_chi(falpha_chi, "alpha_chi");
	alpha_chi.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(alpha_chi);
	
	ParameterT alpha_c(falpha_c, "alpha_c");
	alpha_c.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(alpha_c);
	
	ParameterT alpha_phi(falpha_phi, "alpha_phi");
	alpha_phi.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(alpha_phi);
	
	ParameterT alpha_psi(falpha_psi, "alpha_psi");
	alpha_psi.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(alpha_psi);
	
	ParameterT Tol_1(fTol_1, "Tol_1");
	Tol_1.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(Tol_1);
	
	ParameterT Tol_2(fTol_2, "Tol_2");
	Tol_2.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(Tol_2);
}

/* information about subordinate parameter lists */
void MR3DT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SurfacePotentialT::DefineSubs(sub_list);
}

/* a pointer to the ParameterInterfaceT */
ParameterInterfaceT* MR3DT::NewSub(const StringT& name) const
{
	return SurfacePotentialT::NewSub(name);
}

/* accept parameter list */
void MR3DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SurfacePotentialT::TakeParameterList(list);
	
	/* extract parameters */
	fE_t1 = list.GetParameter("E_t1");
	fE_t2 = list.GetParameter("E_t2");
	fE_n = list.GetParameter("E_n");
	fGf_I = list.GetParameter("Gf_I");
	fGf_II = list.GetParameter("Gf_II");
	fchi_p = list.GetParameter("chi_p");
	fchi_r = list.GetParameter("chi_r");
	fc_p = list.GetParameter("c_p");
	fc_r = list.GetParameter("c_r");
	fphi_p = list.GetParameter("phi_p");
	fphi_r = list.GetParameter("phi_r");
	fpsi_p = list.GetParameter("psi_p");
	falpha_chi = list.GetParameter("alpha_chi");
	falpha_c = list.GetParameter("alpha_c");
	falpha_phi = list.GetParameter("alpha_phi");
	falpha_psi = list.GetParameter("alpha_psi");
	fTol_1 = list.GetParameter("Tol_1");
	fTol_2 = list.GetParameter("Tol_2");

	/* calculate limits on fc_p and fc_r */
	double const1 = tan(fphi_p)*fchi_p;
	if (fc_p < const1) fc_p = const1;
	const1 = tan(fphi_r)*fchi_r;
	if (fc_r < const1) fc_r = const1;

	/* setup output file and format */
	outputPrecision = 10;
	outputFileWidth = outputPrecision + 8;
	mr_ep_3d_out.open("mr_ep_3d.info");
}


/* initialize the state variable array */
void MR3DT::InitStateVariables(ArrayT<double>& state)
{
	int num_state = NumStateVariables();
	if (state.Length() != num_state) {
	#ifndef _SIERRA_TEST_
		cout << "\n SurfacePotentialT::InitStateVariables: expecting state variable array\n"
		     <<   "     length " << num_state << ", found length " << state.Length() << endl;
	#endif
		throw ExceptionT::kSizeMismatch;	
	}

	/* clear */
	if (num_state > 0) state = 0.0;
	
	double enp = state[17];
	double esp = state[18];
	double fchi = fchi_r + (fchi_p - fchi_r)*exp(-falpha_chi*enp);
	double fc   = fc_r + (fc_p - fc_r)*exp(-falpha_c*esp);
	double ftan_phi = tan(fphi_r) + (tan(fphi_p) - tan(fphi_r))*exp(-falpha_phi*esp);
	double ftan_psi = (tan(fpsi_p))*exp(-falpha_psi*esp);	
	state[9] = fchi;
	state[10] = fc ;
	state[11] = ftan_phi;
	state[12] = ftan_psi;
	state[16] = 0.;
}

/* Value of the Yield Function */ 
double MR3DT::YFValue(const ArrayT<double>& state) 
{
	return state[13];
}	

/** dissipated energy. Total amount of energy dissipated reaching
 ** the current state. */
double MR3DT::FractureEnergy(const ArrayT<double>& state)
{
	double FE = state[0]*state[6] + state[1]*state[7] + state[2]*state[8];
	return FE;
}

/** potential energy */
double MR3DT::Potential(const dArrayT& jump_u, const ArrayT<double>& state)
{
	double PT = (jump_u[0]*state[0] + jump_u[1]*state[1] + jump_u[2]*state[2]);
	return PT;
}
/* surface status */
SurfacePotentialT::StatusT MR3DT::Status(const dArrayT& jump_u, 
	const ArrayT<double>& state)
{
	if (state[13]<0.){
		int Status = 0;
	}
	if (state[13]>0) {
		int StatusT = 1;
	}
	if ((jump_u[0]*jump_u[0] + jump_u[1]*jump_u[1] + jump_u[2]*jump_u[2])>100000.) {
		int StatusT = 2;
	}
      
	/*return SurfacePotentialT::Precritical;*/
	return StatusT();
}

void MR3DT::PrintName(ostream& out) const
{
#ifndef _SIERRA_TEST_
	out << "    MR process zone 2D \n";
#endif
}

		
/* traction vector given displacement jump vector */	

const dArrayT& MR3DT::Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma, bool qIntegrate)
{
#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
#endif

	if (! qIntegrate) {
		fTraction[0] = state[0];
		fTraction[1] = state[1];
		fTraction[2] = state[2];
		return fTraction;
	}
	else {
	int i; int j; int kk;

	dMatrixT AA(7,7); dMatrixT KE(3,3); dMatrixT KE_Inv(3,3); dMatrixT I_mat(4,4); dMatrixT I(3,3);
	dMatrixT CMAT(7,7); dMatrixT A_qq(4,4); dMatrixT A_uu(3,3); dMatrixT A_uq(3,4);
	dMatrixT A_qu(4,3); dMatrixT ZMAT(3,4); dMatrixT ZMATP(4,3);
	dMatrixT dQdSig2(3,3); dMatrixT dqbardq(4,4); dMatrixT dQdSigdq(3,4); dMatrixT KE_dQdSig2(3,3); dMatrixT KE_dQdSigdq(3,4);
	dMatrixT dqbardSig(4,3); dMatrixT AA_inv(7,7);

	dArrayT u(3); dArrayT up(3); dArrayT du(3); dArrayT dup(3); dArrayT qn(4);
	dArrayT qo(4); dArrayT Rvec(7); dArrayT Cvec(7); dArrayT upo(3);
	dArrayT R(7); dArrayT Rmod(7); dArrayT Sig(3); dArrayT Sig_I(3);
	dArrayT dQdSig(3); dArrayT dfdq(4); dArrayT qbar(4);
	dArrayT R2(7); dMatrixT X(7,1); dArrayT V_sig(3); dArrayT V_q(4);
	dArrayT dfdSig(3); dArrayT dq(4); dArrayT Y(7); dArrayT KE_dQdSig(3), KE_du(3);

	double ff; double bott; double topp; double dlam; double dlam2; double normr;

		/* Calculate incremental jumps and initialize the neecessary vectors */
		for (i = 0; i<=2; ++i) {
			u[i] = jump_u[i];
			du[i] = u[i] - state[i+3];
			up[i] = state[i+6];
			upo[i] = up[i];
			//Sig_I[i] = 0.;
			Sig_I[i] = state[i];
		}
    
		KE = 0.;
		KE(0,0) = fE_t1;
		KE(1,1) = fE_t2;
		KE(2,2) = fE_n;
		I_mat = 0.;
        	I = 0.;
        	I(0,0) = 1.;
        	I(1,1) = 1.;
		I(2,2) = 1.;
		ZMAT = 0.; ZMATP = 0.;

		KE_Inv.Inverse(KE);
    
		for (i = 0; i<=3; ++i) {
			qn[i] = state[i+9];
			qo[i] = qn[i];
			I_mat(i,i) = 1.;
		}

		dup = 0.;
		//Sig = Sig_I;
		dArrayT ue(3), Sig_e(3);
		ue = u;
		ue -= up;
		KE.MultTx(ue, Sig_e);
		//Sig += Sig_e;
		Sig = Sig_e;


		int iplastic;
		dlam = 0.; dlam2 = 0.; normr = 0.;
    
		/* Check the yield function */
		Yield_f(Sig, qn, ff);

		if (ff < 0.) {
			iplastic = 0;
			state[13] = ff;
			normr = 0.;
			state[16] = normr;
			kk = 0;
		}

		else {
			kk = 0;
			iplastic = 1;
			while (ff > fTol_1 | normr > fTol_2) {
				if (kk > 50000) {
					ExceptionT::GeneralFail("MR3DT::Traction","Too Many Iterations");
				}

				if (kk <= 50) {
#if __option(extended_errorcheck)
				mr_ep_3d_out << setw(outputFileWidth) << "**********" << endl;
				mr_ep_3d_out << setw(outputFileWidth) << "local iteration # = " << kk << endl;
				mr_ep_3d_out << setw(outputFileWidth) << "yield_f = " << ff
					<< setw(outputFileWidth) << "norm = " << normr << endl;
				mr_ep_3d_out << setw(outputFileWidth) << "T_t1 = " << Sig[0] << endl;
				mr_ep_3d_out << setw(outputFileWidth) << "T_t2 = " << Sig[1] << endl;
				mr_ep_3d_out << setw(outputFileWidth) << "T_n = " << Sig[2] << endl;
				mr_ep_3d_out << setw(outputFileWidth) << "u_t1 = " << jump_u[0] << endl;
				mr_ep_3d_out << setw(outputFileWidth) << "u_t2 = " << jump_u[1] << endl;
				mr_ep_3d_out << setw(outputFileWidth) << "u_n = " << jump_u[2] << endl;
				mr_ep_3d_out << setw(outputFileWidth) << "up_t1 = " << up[0] << endl;
				mr_ep_3d_out << setw(outputFileWidth) << "up_t2 = " << up[1] << endl;
				mr_ep_3d_out << setw(outputFileWidth) << "up_n = " << up[2] << endl;
				mr_ep_3d_out << setw(outputFileWidth) << "qn[0] = " << qn[0] << endl;
				mr_ep_3d_out << setw(outputFileWidth) << "qn[0] = " << qn[1] << endl;
				mr_ep_3d_out << setw(outputFileWidth) << "qn[0] = " << qn[2] << endl;
				mr_ep_3d_out << setw(outputFileWidth) << "qn[0] = " << qn[3] << endl;
				mr_ep_3d_out << setw(outputFileWidth) << "**********" << endl;
#endif
				}

				/*
				Sig = Sig_I;
				ue = u;
				ue -= up;
				KE.Multx(ue,Sig_e);
				Sig +=Sig_e;
				Yield_f(Sig, qn, ff);
				*/
				dQdSig_f(Sig, qn, dQdSig);
				qbar_f(Sig, qn, dup, qbar);

                		KE.Multx(du, KE_du);
                		KE.Multx(dQdSig, KE_dQdSig);

				for (i = 0; i<=2; ++i) {
					/*
					R[i]  = upo[i];
					R[i] -= up[i];
					R[i] += dlam*dQdSig[i];
					*/
					R[i]  = Sig[i];
					R[i] -= Sig_I[i];
					R[i] -= KE_du[i];
					R[i] += dlam*KE_dQdSig[i];
				}

				for (i = 0; i<=3; ++i) {
					R[i+3]  = qo[i];
					R[i+3] -= qn[i];
					R[i+3] += dlam*qbar[i];
				}
				/*R[0] = -up[0] + upo[0] + dlam*dQdSig[0];
				  R[1] = -up[1] + upo[1] + dlam*dQdSig[1];
				  R[2] = -up[2] + upo[2] + dlam*dQdSig[2];
				  R[3] = -qn[0] + qo[0] + dlam*qbar[0];
				  R[4] = -qn[1] + qo[1] + dlam*qbar[1];
				  R[5] = -qn[2] + qo[2] + dlam*qbar[2];
				  R[6] = -qn[3] + qo[3] + dlam*qbar[3];*/

				normr = R.Magnitude();
				dQdSig2_f(Sig,qn,dQdSig2);
				dQdSigdq_f(Sig, qn, A_uq);
				dqbardSig_f(Sig, qn, dup, A_qu);
				dqbardq_f(Sig, qn, dup, A_qq);

		                KE_dQdSig2.MultAB(KE, dQdSig2);
                		KE_dQdSigdq.MultAB(KE, A_uq);

				for (i = 0; i<=6; ++i) {
					for (j = 0; j<=6; ++j) {
						/*
						if(i<=2 & j<=2){
							AA_inv(i,j)  = KE_Inv(i,j);
							AA_inv(i,j) += dlam*dQdSig2(i,j);
						}
						if(i<=2 & j >2){
							AA_inv(i,j)  = A_uq(i,j-3);
							AA_inv(i,j) *= dlam;
						}
						if(i >2 & j<=2){
							AA_inv(i,j)  = A_qu(i-3,j);
							AA_inv(i,j) *= dlam;
						}
						if(i >2 & j >2) {
							AA_inv(i,j)  = I_mat(i-3,j-3);
							AA_inv(i,j) *= -1.;
							AA_inv(i,j) += dlam*A_qq(i-3,j-3);
						}
						*/
						if(i<=2 & j<=2){
				                        AA_inv(i,j)  = I(i,j);
                         				AA_inv(i,j) += dlam*KE_dQdSig2(i,j);
						}
						if(i<=2 & j >2){
				                        AA_inv(i,j)  = KE_dQdSigdq(i,j-3);
				                        AA_inv(i,j) *= dlam;
						}
						if(i >2 & j<=2){
							AA_inv(i,j)  = A_qu(i-3,j);
							AA_inv(i,j) *= dlam;
						}
						if(i >2 & j >2) {
							AA_inv(i,j)  = I_mat(i-3,j-3);
							AA_inv(i,j) *= -1.;
							AA_inv(i,j) += dlam*A_qq(i-3,j-3);
						}
					}
				}

				AA.Inverse(AA_inv);
				dfdSig_f(Sig, qn, dfdSig);
				V_sig = dfdSig;
				dfdq_f(Sig, qn, dfdq);
				V_q = dfdq;

				for (i = 0; i<=6; ++i) {
					if (i<=2){
						Rvec[i] = V_sig[i];
						//Cvec[i] = dQdSig[i];
						Cvec[i] = KE_dQdSig[i];
					}
					if (i >2){
						Rvec[i] = V_q[i-3];
						Cvec[i] = qbar[i-3];
					}
				}

				dArrayT tmpVec(7);
				AA.Multx(R,tmpVec);
				topp = ff;
				topp -= dArrayT::Dot(Rvec,tmpVec);
				AA.Multx(Cvec,tmpVec);
				bott = dArrayT::Dot(Rvec,tmpVec);
				dlam2 = topp/bott;

				for (i = 0; i<=6; ++i) {
					for (j = 0; j<=6; ++j) {
						if (i<=2 & j<=2){
							CMAT(i,j) = KE_Inv(i,j);
						}
						if(i<=2 & j >2) {
							CMAT(i,j) = ZMAT(i,j-3);
						}
						if(i >2 & j<=2) {
							CMAT(i,j) = ZMATP(i-3,j);
						}
						if(i >2 & j >2) {
							CMAT(i,j) = -I_mat(i-3,j-3);
						}
					}
				}

				for (i = 0; i<=6; ++i) {
					if (i<=2){
						//Rmod[i] = dQdSig[i];
						Rmod[i] = KE_dQdSig[i];
					}
					if (i >2){
						Rmod[i] = qbar[i-3];
					}
				}

				Rmod *= dlam2;
				R2 = R;
				R2 += Rmod;
				AA.Multx(R2,X);
				CMAT.Multx(X,Y);

				for (i = 0; i<=6; ++i) {
					if (i<=2) {
						dup[i] = Y[i];
					}
					if (i >2) {
						dq[i-3] = Y[i];
					}
				}

				up += dup;
				qn += dq;
				dlam = dlam + dlam2;
				kk = kk + 1;

				// Check the yield function and norm of R
				dup  = up;
				dup -= upo;
				//Sig = Sig_I;
				ue  = u;
				ue -= up;
				KE.Multx(ue, Sig_e);
				//Sig += Sig_e;
				Sig = Sig_e;
				Yield_f(Sig, qn, ff);
				dQdSig_f(Sig, qn, dQdSig);
				qbar_f(Sig, qn, dup, qbar);

				KE.Multx(du, KE_du);
                                KE.Multx(dQdSig, KE_dQdSig);

				for (i = 0; i<=2; ++i) {
					/*
					R[i]  = upo[i];
					R[i] -= up[i];
					R[i] += dlam*dQdSig[i];
					*/
					R[i]  = Sig[i];
                                        R[i] -= Sig_I[i];
                                        R[i] -= KE_du[i];
                                        R[i] += dlam*KE_dQdSig[i];
				}

				for (i = 0; i<=3; ++i) {
					R[i+2]  = qo[i];
					R[i+2] -= qn[i];
					R[i+2] += dlam*qbar[i];
				}

				normr = R.Magnitude();
			}
		}

		state[0] = Sig[0];
		state[1] = Sig[1];
		state[2] = Sig[2];
		fTraction[0] = state[0];
		fTraction[1] = state[1];
		fTraction[2] = state[2];
		state[3] = jump_u[0];
		state[4] = jump_u[1];
		state[5] = jump_u[2];
		state[6] = up[0];
		state[7] = up[1];
		state[8] = up[2];
		state[9] = qn[0];
		state[10] = qn[1];
		state[11] = qn[2];
		state[12] = qn[3];
		state[13] = ff;
		state[14] = dlam;
		state[15] = double(iplastic);
		state[16] = normr;
		dQdSig_f(Sig, qn, dQdSig);
		// need to be checked
		state[17] = Sig[0]*dQdSig[0] + Sig[1]*dQdSig[1];
		state[17] += (Sig[2] + fabs(Sig[2]))*dQdSig[2]/2.;
		state[17] /=fGf_I;
		state[17] *=dlam;
		state[18]  = Sig[0]*dQdSig[0] + Sig[1]*dQdSig[1];
		state[18] *= 1. - fabs((Sig[2] - fabs(Sig[2]))*qn[2]/2.)/sqrt(Sig[0]*Sig[0]+Sig[1]*Sig[1]);
		state[18] /= fGf_II;
		state[18] *=dlam;
		state[19] = double(kk);
		state[20] = dup[0];
		state[21] = dup[1];
		state[22] = dup[2];

#if __option(extended_errorcheck)
		mr_ep_3d_out << setw(outputFileWidth) << "          " << endl;
		mr_ep_3d_out << setw(outputFileWidth) << "----------" << endl;
		mr_ep_3d_out << setw(outputFileWidth) << "yield_f is converged " << "     "
			<< setw(outputFileWidth) << "norm_R is converged " << endl;
		mr_ep_3d_out << setw(outputFileWidth) << "local iteration # = " << kk << endl;
		mr_ep_3d_out << setw(outputFileWidth) << "yield_f = " << ff
			<< setw(outputFileWidth) << "norm = " << normr << endl;
		mr_ep_3d_out << setw(outputFileWidth) << "T_t1 = " << Sig[0] << endl;
		mr_ep_3d_out << setw(outputFileWidth) << "T_t2 = " << Sig[1] << endl;
		mr_ep_3d_out << setw(outputFileWidth) << "T_n = " << Sig[2] << endl;
		mr_ep_3d_out << setw(outputFileWidth) << "u_t1 = " << jump_u[0] << endl;
		mr_ep_3d_out << setw(outputFileWidth) << "u_t2 = " << jump_u[1] << endl;
		mr_ep_3d_out << setw(outputFileWidth) << "u_n = " << jump_u[2] << endl;
		mr_ep_3d_out << setw(outputFileWidth) << "up_t1 = " << up[0] << endl;
		mr_ep_3d_out << setw(outputFileWidth) << "up_t2 = " << up[1] << endl;
		mr_ep_3d_out << setw(outputFileWidth) << "up_n = " << up[2] << endl;
		mr_ep_3d_out << setw(outputFileWidth) << "qn[0] = " << qn[0] << endl;
		mr_ep_3d_out << setw(outputFileWidth) << "qn[0] = " << qn[1] << endl;
		mr_ep_3d_out << setw(outputFileWidth) << "qn[0] = " << qn[2] << endl;
		mr_ep_3d_out << setw(outputFileWidth) << "qn[0] = " << qn[3] << endl;
#endif
		// Check the stiffness
		Stiffness(jump_u, state, sigma);
#if __option(extended_errorcheck)
		mr_ep_3d_out << setw(outputFileWidth) << "KEP(0,0) = " << fStiffness[0] << endl;
		mr_ep_3d_out << setw(outputFileWidth) << "KEP(0,1) = " << fStiffness[1] << endl;
		mr_ep_3d_out << setw(outputFileWidth) << "KEP(0,2) = " << fStiffness[2] << endl;
		mr_ep_3d_out << setw(outputFileWidth) << "KEP(1,0) = " << fStiffness[3] << endl;
		mr_ep_3d_out << setw(outputFileWidth) << "KEP(1,1) = " << fStiffness[4] << endl;
		mr_ep_3d_out << setw(outputFileWidth) << "KEP(1,2) = " << fStiffness[5] << endl;
		mr_ep_3d_out << setw(outputFileWidth) << "KEP(2,0) = " << fStiffness[6] << endl;
		mr_ep_3d_out << setw(outputFileWidth) << "KEP(2,1) = " << fStiffness[7] << endl;
		mr_ep_3d_out << setw(outputFileWidth) << "KEP(2,2) = " << fStiffness[8] << endl;
		mr_ep_3d_out << setw(outputFileWidth) << "----------" << endl;
#endif

		return fTraction;
	}
}

/* calculation of Yield_f */

double& MR3DT::Yield_f(const dArrayT& Sig, const dArrayT& qn, double& ff)
{
	double tmp1, tmp11, tmp22, tmp3, tmp31, tmp32, tmp33, tmp4, tmp5;
  
	tmp1   = qn[1];
	tmp11  = Sig[2];
	tmp11 *= qn[2];
	tmp1  -= tmp11;

	tmp22  = Sig[0];
	tmp22 *= Sig[0];

	tmp33  = Sig[1];
	tmp33 *= Sig[1];

	tmp3   = qn[1];
	tmp31  = qn[0];
	tmp31 *= qn[2];
	tmp3  -= tmp31;
	tmp32  = tmp3;
	tmp32 *= tmp3;

	tmp4  = tmp22 + tmp33;
	tmp4 += tmp32;

	tmp5 = sqrt(tmp4);

	ff  = tmp5;
	ff -= tmp1;

	return ff;
}


/* calculation of qbar_f */

dArrayT& MR3DT::qbar_f(const dArrayT& Sig, const dArrayT& qn, const dArrayT& dup, dArrayT& qbar)
{
/*
	double A1 = -falpha_chi*(qn[0] - fchi_r);
	double B1 = (Sig[2]+fabs(Sig[2]))/2./fGf_I;
	double B2 = Sig[0]/fGf_I;
	double B3 = Sig[1]/fGf_I;
	double DQDN = 2.*qn[3]*(qn[1] - Sig[2]*qn[3]);
	double DQDT1 = 2.*Sig[0];
	double DQDT2 = 2.*Sig[1];
	double A2 = -falpha_c*(qn[1] - fc_r);
	double TNA = (Sig[2] - fabs(Sig[2]))/2.;
	double B4 = (1. - fabs(TNA*qn[2])/sqrt(Sig[0]*Sig[0]+Sig[1]*Sig[1]))*Sig[0]/fGf_II;
	double B5 = (1. - fabs(TNA*qn[2])/sqrt(Sig[0]*Sig[0]+Sig[1]*Sig[1]))*Sig[1]/fGf_II;
	double A3 = -falpha_phi*(qn[2] - tan(fphi_r));
	double A4 = -falpha_psi*qn[3];

	qbar[0] = A1*B1*DQDN + A1*B2*DQDT1 + A1*B3*DQDT2;
	qbar[1] = A2*B4*DQDT1 + A2*B5*DQDT2;
	qbar[2] = A3*B4*DQDT1 + A3*B5*DQDT2;
	qbar[3] = A4*B4*DQDT1 + A4*B5*DQDT2;
*/
	double B4, B5;

	double A1 = -falpha_chi*(qn[0] - fchi_r);
	double A2 = -falpha_c*(qn[1] - fc_r);
	double A3 = -falpha_phi*(qn[2] - tan(fphi_r));
	double A4 = -falpha_psi*qn[3];
	double B1 = (Sig[2] + fabs(Sig[2]))/(2.*fGf_I);
	double B2 = Sig[0]/fGf_I;
	double B3 = Sig[1]/fGf_I;
	double TNA = (Sig[2] - fabs(Sig[2]))/2.;

	if (dup[0] == 0.) {
		B4 = 0.;
	}
	else {
		B4 = (sqrt(Sig[0]*Sig[0] + Sig[1]*Sig[1]) - fabs(TNA*qn[2]))*dup[0]/(sqrt(dup[0]*dup[0] + dup[1]*dup[1])*fGf_II);
	}

	if (dup[1] == 0.) {
		B5 = 0.;
	}
	else {
		B5 = (sqrt(Sig[0]*Sig[0] + Sig[1]*Sig[1]) - fabs(TNA*qn[2]))*dup[1]/(sqrt(dup[0]*dup[0] + dup[1]*dup[1])*fGf_II);
	}

	//double B4 = (sqrt(Sig[0]*Sig[0] + Sig[1]*Sig[1]) - fabs(TNA*qn[2]))*dup[0]/(sqrt(dup[0]*dup[0] + dup[1]*dup[1])*fGf_II);
	//double B5 = (sqrt(Sig[0]*Sig[0] + Sig[1]*Sig[1]) - fabs(TNA*qn[2]))*dup[1]/(sqrt(dup[0]*dup[0] + dup[1]*dup[1])*fGf_II);
	double Shear_Q = Sig[0]*Sig[0] + Sig[1]*Sig[1] + (qn[1] - qn[0]*qn[3])*(qn[1] - qn[0]*qn[3]);
	double DQDT1 = Sig[0]/sqrt(Shear_Q);
	double DQDT2 = Sig[1]/sqrt(Shear_Q);
	double DQDTN = qn[3];

	qbar[0] = A1*B1*DQDTN + A1*B2*DQDT1 + A1*B3*DQDT2;
	qbar[1] = A2*B4*DQDT1 + A2*B5*DQDT2;
	qbar[2] = A3*B4*DQDT1 + A3*B5*DQDT2;
	qbar[3] = A4*B4*DQDT1 + A4*B5*DQDT2;

	return qbar;
 }


/* calculation of dQdSig2_f */

dMatrixT& MR3DT::dQdSig2_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dQdSig2)
{
#pragma unused(Sig)
/*
	dQdSig2(0,0) = 2.;
	dQdSig2(1,1) = 2.;
	dQdSig2(2,2) = -2.*qn[3]*qn[3];
	dQdSig2(0,1) = 0.;
	dQdSig2(0,2) = 0.;
	dQdSig2(1,0) = 0.;
	dQdSig2(1,2) = 0.;
	dQdSig2(2,0) = 0.;
	dQdSig2(2,1) = 0.;
*/
	double Shear_Q = Sig[0]*Sig[0] + Sig[1]*Sig[1] + (qn[1] - qn[0]*qn[3])*(qn[1] - qn[0]*qn[3]);
	dQdSig2(0,0) = (Sig[1]*Sig[1] + (qn[1] - qn[0]*qn[3])*(qn[1] - qn[0]*qn[3]))/(sqrt(Shear_Q)*sqrt(Shear_Q)*sqrt(Shear_Q));
	dQdSig2(0,1) = -Sig[0]*Sig[1]/(sqrt(Shear_Q)*sqrt(Shear_Q)*sqrt(Shear_Q));
	dQdSig2(0,2) = 0.;
	dQdSig2(1,0) = -Sig[0]*Sig[1]/(sqrt(Shear_Q)*sqrt(Shear_Q)*sqrt(Shear_Q));
	dQdSig2(1,1) = (Sig[0]*Sig[0] + (qn[1] - qn[0]*qn[3])*(qn[1] - qn[0]*qn[3]))/(sqrt(Shear_Q)*sqrt(Shear_Q)*sqrt(Shear_Q));
	dQdSig2(1,2) = 0.;
	dQdSig2(2,0) = 0.;
	dQdSig2(2,1) = 0.;
	dQdSig2(2,2) = 0.;

	return dQdSig2;
	return dQdSig2;
}

/* calculation of dfdSig_f */

dArrayT& MR3DT::dfdSig_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& dfdSig)
{
	double Shear = Sig[0]*Sig[0] + Sig[1]*Sig[1] + (qn[1] - qn[0]*qn[2])*(qn[1] - qn[0]*qn[2]);
	dfdSig[0] = Sig[0]/sqrt(Shear);
	dfdSig[1] = Sig[1]/sqrt(Shear);
	dfdSig[2] = qn[2];

	return dfdSig;
}

/* calculation of dQdSig_f */

dArrayT& MR3DT::dQdSig_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& dQdSig)
{
/*
	dQdSig[0] = 2.*Sig[0];
	dQdSig[1] = 2.*Sig[1];
	dQdSig[2] = 2.*qn[3]*(qn[1] - Sig[2]*qn[3]);
*/
	double Shear_Q = Sig[0]*Sig[0] + Sig[1]*Sig[1] + (qn[1] - qn[0]*qn[3])*(qn[1] - qn[0]*qn[3]);
	dQdSig[0] = Sig[0]/sqrt(Shear_Q);
	dQdSig[1] = Sig[1]/sqrt(Shear_Q);
	dQdSig[2] = qn[3];

	return dQdSig;
}


/* calculation of dfdq_f */

dArrayT& MR3DT::dfdq_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& dfdq)
{
	double Shear = Sig[0]*Sig[0] + Sig[1]*Sig[1] + (qn[1] - qn[0]*qn[2])*(qn[1] - qn[0]*qn[2]);
	double zeta = (qn[1] - qn[0]*qn[2])/sqrt(Shear);
	dfdq[0] = -qn[2]*zeta;
	dfdq[1] = zeta - 1.;
	dfdq[2] = -qn[0]*zeta + Sig[2];
	dfdq[3] = 0.;

	return dfdq;
}

/* calculation of dQdSigdq_f */

dMatrixT& MR3DT::dQdSigdq_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dQdSigdq)
{
/*
	dQdSigdq(0,0) = 0.;
	dQdSigdq(0,1) = 0.;
	dQdSigdq(0,2) = 0.;
	dQdSigdq(0,3) = 0.;
	dQdSigdq(1,0) = 0.;
	dQdSigdq(1,1) = 0.;
	dQdSigdq(1,2) = 0.;
	dQdSigdq(1,3) = 0.;
	dQdSigdq(2,0) = 0.;
	dQdSigdq(2,1) = 2.*qn[3];
	dQdSigdq(2,2) = 0.;
	dQdSigdq(2,3) = 2.*qn[1] - 4.*Sig[2]*qn[3];
*/
	double Shear_Q = Sig[0]*Sig[0] + Sig[1]*Sig[1] + (qn[1] - qn[0]*qn[3])*(qn[1] - qn[0]*qn[3]);
	double zeta_q = (qn[1] - qn[0]*qn[3])/(sqrt(Shear_Q)*sqrt(Shear_Q)*sqrt(Shear_Q));
	dQdSigdq(0,0) = Sig[0]*qn[3]*zeta_q;
	dQdSigdq(0,1) = -Sig[0]*zeta_q;
	dQdSigdq(0,2) = 0.;
	dQdSigdq(0,3) = Sig[0]*qn[0]*zeta_q;
	dQdSigdq(1,0) = Sig[1]*qn[3]*zeta_q;
	dQdSigdq(1,1) = -Sig[1]*zeta_q;
	dQdSigdq(1,2) = 0.;
	dQdSigdq(1,3) = Sig[1]*qn[0]*zeta_q;
	dQdSigdq(2,0) = 0.;
	dQdSigdq(2,1) = 0.;
	dQdSigdq(2,2) = 0.;
	dQdSigdq(2,3) = 1.;

	return dQdSigdq;
}

/* calculation of dqbardSig_f */

dMatrixT& MR3DT::dqbardSig_f(const dArrayT& Sig, const dArrayT& qn, const dArrayT& dup, dMatrixT& dqbardSig)
{
/*
	double A1 = -falpha_chi*(qn[0] - fchi_r);
	double B1 = (Sig[2]+fabs(Sig[2]))/2./fGf_I;
	double B2 = Sig[0]/fGf_I;
	double B3 = Sig[1]/fGf_I;
	double DQDN = 2.*qn[3]*(qn[1] - Sig[2]*qn[3]);
	double DQDT1 = 2.*Sig[0];
	double DQDT2 = 2.*Sig[1];
	double A2 = -falpha_c*(qn[1] - fc_r);
	double TNA = (Sig[2] - fabs(Sig[2]))/2.;
	double B4 = (1. - fabs(TNA*qn[2])/sqrt(Sig[0]*Sig[0]+Sig[1]*Sig[1]))*Sig[0]/fGf_II;
	double B5 = (1. - fabs(TNA*qn[2])/sqrt(Sig[0]*Sig[0]+Sig[1]*Sig[1]))*Sig[1]/fGf_II;
	double A3 = -falpha_phi*(qn[2] - tan(fphi_r));
	double A4 = -falpha_psi*qn[3];
	double DB4_DTn = -qn[2]*Sig[0]*signof(TNA)*(1. - signof(Sig[2]))/sqrt(Sig[0]*Sig[0]+Sig[1]*Sig[1])/fGf_II/2.;
	double DB5_DTn = -qn[2]*Sig[1]*signof(TNA)*(1. - signof(Sig[2]))/sqrt(Sig[0]*Sig[0]+Sig[1]*Sig[1])/fGf_II/2.;
	double DB4_DTt1 = (1. - Sig[1]*Sig[1]*fabs(TNA*qn[2])/sqrt(Sig[0]*Sig[0]+Sig[1]*Sig[1])/(Sig[0]*Sig[0]+Sig[1]*Sig[1]))/fGf_II;
	double DB4_DTt2 = (Sig[0]*Sig[1]*fabs(TNA*qn[2])/sqrt(Sig[0]*Sig[0]+Sig[1]*Sig[1])/(Sig[0]*Sig[0]+Sig[1]*Sig[1]))/fGf_II;
	double DB5_DTt1 = (Sig[0]*Sig[1]*fabs(TNA*qn[2])/sqrt(Sig[0]*Sig[0]+Sig[1]*Sig[1])/(Sig[0]*Sig[0]+Sig[1]*Sig[1]))/fGf_II;
	double DB5_DTt2 = (1. - Sig[0]*Sig[0]*fabs(TNA*qn[2])/sqrt(Sig[0]*Sig[0]+Sig[1]*Sig[1])/(Sig[0]*Sig[0]+Sig[1]*Sig[1]))/fGf_II;
	double DQDN2 = -2.*qn[3]*qn[3];
	double DQDT3 = 2.;
	double DQDTN = 0.;
	double DQDNT = 0.;
	double SN = signof(Sig[2]);
	double DB1DN = (SN +fabs(SN))/2./fGf_I;

	dqbardSig(0,0) = A1*B2*DQDT3 + A1*DQDT1/fGf_I;
	dqbardSig(0,1) = A1*B3*DQDT3 + A1*DQDT2/fGf_I;
	dqbardSig(0,2) = A1*B1*DQDN2 + A1*DQDN*DB1DN;
	dqbardSig(1,0) = A2*B4*DQDT3 + A2*DQDT1*DB4_DTt1 + A2*DQDT2*DB5_DTt1;
	dqbardSig(1,1) = A2*B5*DQDT3 + A2*DQDT2*DB5_DTt2 + A2*DQDT1*DB4_DTt2;
	dqbardSig(1,2) = A2*DQDT1*DB4_DTn + A2*DQDT2*DB5_DTn;
	dqbardSig(2,0) = A3*B4*DQDT3 + A3*DQDT1*DB4_DTt1 + A3*DQDT2*DB5_DTt1;
	dqbardSig(2,1) = A3*B5*DQDT3 + A3*DQDT2*DB5_DTt2 + A3*DQDT1*DB4_DTt2;
	dqbardSig(2,2) = A3*DQDT1*DB4_DTn + A3*DQDT2*DB5_DTn;
	dqbardSig(3,0) = A4*B4*DQDT3 + A4*DQDT1*DB4_DTt1 + A4*DQDT2*DB5_DTt1;
	dqbardSig(3,1) = A4*B5*DQDT3 + A4*DQDT2*DB5_DTt2 + A4*DQDT1*DB4_DTt2;
	dqbardSig(3,2) = A4*DQDT1*DB4_DTn + A4*DQDT2*DB5_DTn;
*/
	double B4, DB4DT1, DB4DT2, DB4DTN, B5, DB5DT1, DB5DT2, DB5DTN;

	double A1 = -falpha_chi*(qn[0] - fchi_r);
	double A2 = -falpha_c*(qn[1] - fc_r);
	double A3 = -falpha_phi*(qn[2] - tan(fphi_r));
	double A4 = -falpha_psi*qn[3];
	// double B1 = (Sig[2] + fabs(Sig[2]))/(2.*fGf_I);
	double B2 = Sig[0]/fGf_I;
	double B3 = Sig[1]/fGf_I;
	double TNA = (Sig[2] - fabs(Sig[2]))/2.;
	double DTNADTN = (1. - signof(Sig[2]))/2.;

	if (dup[0] == 0.) {
		B4 = 0.;
		DB4DT1 = 0.;
		DB4DT2 = 0.;
		DB4DTN = 0.;
	}
	else {
		B4 = (sqrt(Sig[0]*Sig[0] + Sig[1]*Sig[1]) - fabs(TNA*qn[2]))*dup[0]/(sqrt(dup[0]*dup[0] + dup[1]*dup[1])*fGf_II);
		DB4DT1 = Sig[0]*dup[0]/(sqrt(Sig[0]*Sig[0] + Sig[1]*Sig[1])*sqrt(dup[0]*dup[0] + dup[1]*dup[1])*fGf_II);
		DB4DT2 = Sig[1]*dup[0]/(sqrt(Sig[0]*Sig[0] + Sig[1]*Sig[1])*sqrt(dup[0]*dup[0] + dup[1]*dup[1])*fGf_II);
		DB4DTN = -signof(TNA*qn[2])*qn[2]*DTNADTN*dup[0]/(sqrt(dup[0]*dup[0] + dup[1]*dup[1])*fGf_II);
	}

	if (dup[1] == 0.) {
		B5 = 0.;
		DB5DT1 = 0.;
		DB5DT2 = 0.;
		DB5DTN = 0.;
	}
	else {
		B5 = (sqrt(Sig[0]*Sig[0] + Sig[1]*Sig[1]) - fabs(TNA*qn[2]))*dup[1]/(sqrt(dup[0]*dup[0] + dup[1]*dup[1])*fGf_II);
		DB5DT1 = Sig[0]*dup[1]/(sqrt(Sig[0]*Sig[0] + Sig[1]*Sig[1])*sqrt(dup[0]*dup[0] + dup[1]*dup[1])*fGf_II);
		DB5DT2 = Sig[1]*dup[1]/(sqrt(Sig[0]*Sig[0] + Sig[1]*Sig[1])*sqrt(dup[0]*dup[0] + dup[1]*dup[1])*fGf_II);
		DB5DTN = -signof(TNA*qn[2])*qn[2]*DTNADTN*dup[1]/(sqrt(dup[0]*dup[0] + dup[1]*dup[1])*fGf_II);
	}

	double Shear_Q = Sig[0]*Sig[0] + Sig[1]*Sig[1] + (qn[1] - qn[0]*qn[3])*(qn[1] - qn[0]*qn[3]);
	double DQDT1 = Sig[0]/sqrt(Shear_Q);
	double DQDT2 = Sig[1]/sqrt(Shear_Q);
	double DQDTN = qn[3];
	double DQDT1_2  = (Sig[1]*Sig[1] + (qn[1] - qn[0]*qn[3])*(qn[1] - qn[0]*qn[3]))/(sqrt(Shear_Q)*sqrt(Shear_Q)*sqrt(Shear_Q));
	double DQDT1DT2 = -Sig[0]*Sig[1]/(sqrt(Shear_Q)*sqrt(Shear_Q)*sqrt(Shear_Q));
	double DQDT2_2  = (Sig[0]*Sig[0] + (qn[1] - qn[0]*qn[3])*(qn[1] - qn[0]*qn[3]))/(sqrt(Shear_Q)*sqrt(Shear_Q)*sqrt(Shear_Q));
	double DB1DTN = (1. + signof(Sig[2]))/(2.*fGf_I);
	double DB2DT1 = 1./fGf_I;
	double DB3DT2 = 1./fGf_I;
	//double DB4DT1 = Sig[0]*dup[0]/(sqrt(Sig[0]*Sig[0] + Sig[1]*Sig[1])*sqrt(dup[0]*dup[0] + dup[1]*dup[1])*fGf_II);
	//double DB4DT2 = Sig[1]*dup[0]/(sqrt(Sig[0]*Sig[0] + Sig[1]*Sig[1])*sqrt(dup[0]*dup[0] + dup[1]*dup[1])*fGf_II);
	//double DB4DTN = -signof(TNA*qn[2])*qn[2]*DTNADTN*dup[0]/(sqrt(dup[0]*dup[0] + dup[1]*dup[1])*fGf_II);
	//double DB5DT1 = Sig[0]*dup[1]/(sqrt(Sig[0]*Sig[0] + Sig[1]*Sig[1])*sqrt(dup[0]*dup[0] + dup[1]*dup[1])*fGf_II);
	//double DB5DT2 = Sig[1]*dup[1]/(sqrt(Sig[0]*Sig[0] + Sig[1]*Sig[1])*sqrt(dup[0]*dup[0] + dup[1]*dup[1])*fGf_II);
	//double DB5DTN = -signof(TNA*qn[2])*qn[2]*DTNADTN*dup[1]/(sqrt(dup[0]*dup[0] + dup[1]*dup[1])*fGf_II);

	dqbardSig(0,0) = A1*DB2DT1*DQDT1 + A1*B2*DQDT1_2 + A1*B3*DQDT1DT2;
	dqbardSig(0,1) = A1*B2*DQDT1DT2 + A1*DB3DT2*DQDT2 + A1*B3*DQDT2_2;
	dqbardSig(0,2) = A1*DB1DTN*DQDTN;
	dqbardSig(1,0) = A2*DB4DT1*DQDT1 + A2*B4*DQDT1_2 + A2*DB5DT1*DQDT2 + A2*B5*DQDT1DT2;
	dqbardSig(1,1) = A2*DB4DT2*DQDT1 + A2*B4*DQDT1DT2 + A2*DB5DT2*DQDT2 + A2*B5*DQDT2_2;
	dqbardSig(1,2) = A2*DB4DTN*DQDT1 + A2*DB5DTN*DQDT2;
	dqbardSig(2,0) = A3*DB4DT1*DQDT1 + A3*B4*DQDT1_2 + A3*DB5DT1*DQDT2 + A3*B5*DQDT1DT2;
	dqbardSig(2,1) = A3*DB4DT2*DQDT1 + A3*B4*DQDT1DT2 + A3*DB5DT2*DQDT2 + A3*B5*DQDT2_2;
	dqbardSig(2,2) = A3*DB4DTN*DQDT1 + A3*DB5DTN*DQDT2;
	dqbardSig(3,0) = A4*DB4DT1*DQDT1 + A4*B4*DQDT1_2 + A4*DB5DT1*DQDT2 + A4*B5*DQDT1DT2;
	dqbardSig(3,1) = A4*DB4DT2*DQDT1 + A4*B4*DQDT1DT2 + A4*DB5DT2*DQDT2 + A4*B5*DQDT2_2;
	dqbardSig(3,2) = A4*DB4DTN*DQDT1 + A4*DB5DTN*DQDT2;

	return dqbardSig;
}
  
/* calculation of dqbardq_f */

dMatrixT& MR3DT::dqbardq_f(const dArrayT& Sig, const dArrayT& qn, const dArrayT& dup, dMatrixT& dqbardq)
{
/*
	double A1 = -falpha_chi*(qn[0] - fchi_r);
	double B1 = (Sig[2]+fabs(Sig[2]))/2./fGf_I;
	double B2 = Sig[0]/fGf_I;
	double B3 = Sig[1]/fGf_I;
	double DQDN = 2.*qn[3]*(qn[1] - Sig[2]*qn[3]);
	double DQDT1 = 2.*Sig[0];
	double DQDT2 = 2.*Sig[1];
	double A2 = -falpha_c*(qn[1] - fc_r);
	double TNA = (Sig[2] - fabs(Sig[2]))/2.;
	double B4 = (1. - fabs(TNA*qn[2])/sqrt(Sig[0]*Sig[0]+Sig[1]*Sig[1]))*Sig[0]/fGf_II;
	double B5 = (1. - fabs(TNA*qn[2])/sqrt(Sig[0]*Sig[0]+Sig[1]*Sig[1]))*Sig[1]/fGf_II;
	double A3 = -falpha_phi*(qn[2] - tan(fphi_r));
	double A4 = -falpha_psi*qn[3];
	double DB4_DTn = -qn[2]*Sig[0]*signof(TNA)*(1. - signof(Sig[2]))/sqrt(Sig[0]*Sig[0]+Sig[1]*Sig[1])/fGf_II/2.;
	double DB5_DTn = -qn[2]*Sig[1]*signof(TNA)*(1. - signof(Sig[2]))/sqrt(Sig[0]*Sig[0]+Sig[1]*Sig[1])/fGf_II/2.;
	double DB4_DTt1 = (1. - Sig[1]*Sig[1]*fabs(TNA*qn[2])/sqrt(Sig[0]*Sig[0]+Sig[1]*Sig[1])/(Sig[0]*Sig[0]+Sig[1]*Sig[1]))/fGf_II;
	double DB4_DTt2 = (Sig[0]*Sig[1]*fabs(TNA*qn[2])/sqrt(Sig[0]*Sig[0]+Sig[1]*Sig[1])/(Sig[0]*Sig[0]+Sig[1]*Sig[1]))/fGf_II;
	double DB5_DTt1 = (Sig[0]*Sig[1]*fabs(TNA*qn[2])/sqrt(Sig[0]*Sig[0]+Sig[1]*Sig[1])/(Sig[0]*Sig[0]+Sig[1]*Sig[1]))/fGf_II;
	double DB5_DTt2 = (1. - Sig[0]*Sig[0]*fabs(TNA*qn[2])/sqrt(Sig[0]*Sig[0]+Sig[1]*Sig[1])/(Sig[0]*Sig[0]+Sig[1]*Sig[1]))/fGf_II;
	double DB4_DTanphi = -fabs(TNA)*Sig[0]/sqrt(Sig[0]*Sig[0]+Sig[1]*Sig[1])/fGf_II;
	double DB5_DTanphi = -fabs(TNA)*Sig[1]/sqrt(Sig[0]*Sig[0]+Sig[1]*Sig[1])/fGf_II;
	double DQDN2 = -2.*qn[3]*qn[3];
	double DQDT3 = 2.;
	double DQDTN = 0.;
	double DQDNT = 0.;
	double SN = signof(Sig[2]);
	double DB1DN = (SN +fabs(SN))/2./fGf_I;

	dqbardq(0,0) = -falpha_chi*(B1*DQDN + B2*DQDT1 + B3*DQDT2);
	dqbardq(0,1) =  A1*B1*(2.*qn[3]);
	dqbardq(0,2) = 0.;
	dqbardq(0,3) =  A1*B1*(2.*qn[1]-4.*Sig[2]*qn[3]);
	dqbardq(1,0) = 0.;
	dqbardq(1,1) = -falpha_c*(B4*DQDT1 + B5*DQDT2);
	dqbardq(1,2) = A2*(DQDT1*DB4_DTanphi + DQDT2*DB5_DTanphi);
	dqbardq(1,3) = 0.;
	dqbardq(2,0) = 0.;
	dqbardq(2,1) = 0.;
	dqbardq(2,2) = -falpha_phi*(B4*DQDT1 + B5*DQDT2) + A3*(DQDT1*DB4_DTanphi + DQDT2*DB5_DTanphi);
	dqbardq(2,3) = 0.;
	dqbardq(3,0) = 0.;
	dqbardq(3,1) = 0.;
	dqbardq(3,2) = A4*(DQDT1*DB4_DTanphi + DQDT2*DB5_DTanphi);
	dqbardq(3,3) = -falpha_psi*(B4*DQDT1 + B5*DQDT2);
*/
	double B4, DB4DTanphi, B5, DB5DTanphi;

	double A1 = -falpha_chi*(qn[0] - fchi_r);
	double A2 = -falpha_c*(qn[1] - fc_r);
	double A3 = -falpha_phi*(qn[2] - tan(fphi_r));
	double A4 = -falpha_psi*qn[3];
	double B1 = (Sig[2] + fabs(Sig[2]))/(2.*fGf_I);
	double B2 = Sig[0]/fGf_I;
	double B3 = Sig[1]/fGf_I;
	double TNA = (Sig[2] - fabs(Sig[2]))/2.;

	if (dup[0] == 0.) {
		B4 = 0.;
		DB4DTanphi = 0.;
	}
	else {
		B4 = (sqrt(Sig[0]*Sig[0] + Sig[1]*Sig[1]) - fabs(TNA*qn[2]))*dup[0]/(sqrt(dup[0]*dup[0] + dup[1]*dup[1])*fGf_II);
		DB4DTanphi = -signof(TNA*qn[2])*TNA*dup[0]/(sqrt(dup[0]*dup[0] + dup[1]*dup[1])*fGf_II);
	}

	if (dup[1] == 0.) {
		B5 = 0.;
		DB5DTanphi = 0.;
	}
	else {
		B5 = (sqrt(Sig[0]*Sig[0] + Sig[1]*Sig[1]) - fabs(TNA*qn[2]))*dup[1]/(sqrt(dup[0]*dup[0] + dup[1]*dup[1])*fGf_II);
		DB5DTanphi = -signof(TNA*qn[2])*TNA*dup[1]/(sqrt(dup[0]*dup[0] + dup[1]*dup[1])*fGf_II);
	}

	double Shear_Q = Sig[0]*Sig[0] + Sig[1]*Sig[1] + (qn[1] - qn[0]*qn[3])*(qn[1] - qn[0]*qn[3]);
	double zeta_q = (qn[1] - qn[0]*qn[3])/(sqrt(Shear_Q)*sqrt(Shear_Q)*sqrt(Shear_Q));
	double DQDT1 = Sig[0]/sqrt(Shear_Q);
	double DQDT2 = Sig[1]/sqrt(Shear_Q);
	double DQDTN = qn[3];
	double DQDT1DChi = Sig[0]*qn[3]*zeta_q;
	double DQDT1DC = -Sig[0]*zeta_q;
	double DQDT1DTanpsi = Sig[0]*qn[0]*zeta_q;
	double DQDT2DChi = Sig[1]*qn[3]*zeta_q;
	double DQDT2DC = -Sig[1]*zeta_q;
	double DQDT2DTanpsi = Sig[1]*qn[0]*zeta_q;
	double DQDTNDTanpsi = 1.;
	//double DB4DTanphi = -signof(TNA*qn[2])*TNA*dup[0]/(sqrt(dup[0]*dup[0] + dup[1]*dup[1])*fGf_II);
	//double DB5DTanphi = -signof(TNA*qn[2])*TNA*dup[1]/(sqrt(dup[0]*dup[0] + dup[1]*dup[1])*fGf_II);

	dqbardq(0,0) = -falpha_chi*(B1*DQDTN + B2*DQDT1 + B3*DQDT2) + A1*B2*DQDT1DChi + A1*B3*DQDT2DChi;
	dqbardq(0,1) = A1*B2*DQDT1DC + A1*B3*DQDT2DC;
	dqbardq(0,2) = 0.;
	dqbardq(0,3) = A1*B1*DQDTNDTanpsi + A1*B2*DQDT1DTanpsi + A1*B3*DQDT2DTanpsi;
	dqbardq(1,0) = A2*B4*DQDT1DChi + A2*B5*DQDT2DChi;
	dqbardq(1,1) = -falpha_c*(B4*DQDT1 + B5*DQDT2) + A2*B4*DQDT1DC + A2*B5*DQDT2DC;
	dqbardq(1,2) = A2*(DB4DTanphi*DQDT1 + DB5DTanphi*DQDT2);
	dqbardq(1,3) = A2*B4*DQDT1DTanpsi + A2*B5*DQDT2DTanpsi;
	dqbardq(2,0) = A3*B4*DQDT1DChi + A3*B5*DQDT2DChi;
	dqbardq(2,1) = A3*B4*DQDT1DC + A3*B5*DQDT2DC;
	dqbardq(2,2) = -falpha_phi*(B4*DQDT1 + B5*DQDT2) + A3*(DB4DTanphi*DQDT1 + DB5DTanphi*DQDT2);
	dqbardq(2,3) = A3*B4*DQDT1DTanpsi + A3*B5*DQDT2DTanpsi;
	dqbardq(3,0) = A4*B4*DQDT1DChi + A4*B5*DQDT2DChi;
	dqbardq(3,1) = A4*B4*DQDT1DC + A4*B5*DQDT2DC;
	dqbardq(3,2) = A4*(DB4DTanphi*DQDT1 + DB5DTanphi*DQDT2);
	dqbardq(3,3) = -falpha_psi*(B4*DQDT1 + B5*DQDT2) + A4*B4*DQDT1DTanpsi + A4*B5*DQDT2DTanpsi;;

	return dqbardq;
}

/* elastoplastic consistent tangent operator*/
const dMatrixT& MR3DT::Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma)
{
#pragma unused(sigma)
#pragma unused(jump_u)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif

	int i, j;

	dMatrixT AA(7,7), I_mat(4,4), CMAT(7,7),AA_inv(7,7),
		A_qq(4,4), A_uu(3,3), A_uq(3,4), A_qu(4,3), ZMAT(3,4),
		ZMATP(4,3), dQdSig2(3,3), dqdbar(4,4), dqbardSig(4,3),
		dQdSigdq(3,4), KP(3,3), KP2(3,3), KEP(3,3), DMAT(7,7), EMAT(7,3), FMAT(7,3);

	dMatrixT I_m(3,3), Rmat(3,3), R_Inv(3,3), KE(3,3), KE_Inv(3,3),
		Ch(4,4), Ch_Inv(4,4), KE1(4,3), KE2(3,3), KE3(3,4), KEE(3,3),
		KEE_Inv(3,3);

	dArrayT  u(3), up(3), du(3), dup(3), qn(4), qo(4), Rvec(7),Cvec(7),
		R(7), Rmod(7), Sig(3), Sig_I(3), dQdSig(3), dfdq(4), qbar(4),
		R2(7), X(7), V_sig(3), V_q(4), dfdSig(3), K1(3), K2(3);

	double bott, dlam, e;

	fStiffness[1] = fStiffness[2] = fStiffness[3] = fStiffness[5] = fStiffness[6] = fStiffness[7] = 0.;
	I_m = 0.;
	I_m(0,0) = 1.; I_m(1,1) = 1.; I_m(2,2) = 1.;
	I_mat = 0.;
	Sig[0] = state[0];
	Sig[1] = state[1];
	Sig[2] = state[2];
	dup[0] = state[20];
	dup[1] = state[21];
	dup[2] = state[22];
	KEE = 0.;
	KEE(0,0) = fE_t1;
	KEE(1,1) = fE_t2;
	KEE(2,2) = fE_n;

	KEE_Inv.Inverse(KEE);

	for (i = 0; i<=3; ++i) {
		qn[i] = state[i+9];
		I_mat(i,i) = 1.;
	}

	if (state[15] == 0.) {
		fStiffness[0] = fE_t1;
		fStiffness[4] = fE_t2;
		fStiffness[8] = fE_n;
	}
	else if (state[15] == 1.) {
		dlam = state[14];
		dQdSig2_f(Sig, qn, dQdSig2);
		dqbardSig_f(Sig, qn, dup, A_qu);
		dqbardq_f(Sig, qn, dup, A_qq);
		dQdSigdq_f(Sig, qn, A_uq);
/*
		Ch  = A_qq;
		Ch *= -dlam;
		Ch += I_mat;
		Ch_Inv.Inverse(Ch);
		KE1.MultAB(Ch_Inv,A_qu);
		KE.MultAB(A_uq,KE1);
		KE *= state[14];
		KE *= state[14];
		//KE = 0.;
		KE2 = dQdSig2;
		KE2 *=state[14];
		KE += KE2;
		KE += KEE_Inv;

		KE_Inv.Inverse(KE);
*/

		for (i = 0; i<=6; ++i) {
			for (j = 0; j<=6; ++j){
				if (i<=2 & j<=2){
					AA_inv(i,j)  = KEE_Inv(i,j);
					AA_inv(i,j) += dlam*dQdSig2(i,j);
				}
				if (i<=2 & j>2){
					AA_inv(i,j)  = A_uq(i,j-3);
					AA_inv(i,j) *= dlam;
				}
				if(i>2 & j<=2){
					AA_inv(i,j)  = A_qu(i-3,j);
					AA_inv(i,j) *= dlam;
				}
				if(i>2 & j >2) {
					AA_inv(i,j)  = I_mat(i-3,j-3);
					AA_inv(i,j) *= -1.;
					AA_inv(i,j) += dlam*A_qq(i-3,j-3);
				}
			}
		}

		AA.Inverse(AA_inv);

		dfdSig_f(Sig, qn, dfdSig);
		dfdq_f(Sig,qn, dfdq);
		dQdSig_f(Sig, qn, dQdSig);
		qbar_f(Sig, qn, dup, qbar);

		for (i = 0; i<=6; ++i) {
			if (i<=2) {
				Rvec[i] = dfdSig[i];
				Cvec[i] = dQdSig[i];
			}
			if (i>2) {
				Rvec[i] = dfdq[i-3];
				Cvec[i] = qbar[i-3];
			}
		}

		dArrayT tmpVec(7), Vvec(3), dVec(3);
		AA.Multx(Cvec,tmpVec);
		//bott = dArrayT::Dot(Rvec,tmpVec);
		e = dArrayT::Dot(Rvec,tmpVec);
/*
		for (i = 0; i<=2; ++i) {
			Vvec[i] = 0.;
			for (j = 0; j<=6; ++j) {
				Vvec[i] += Rvec[j]*AA(j,i);
			}
	        }

		for (i = 0; i<=2; ++i) {
			for (j = 0; j<=2; ++j) {
				KP(i,j) = dQdSig[i]*Vvec[j];
			}
		}

		KE3.MultAB(A_uq, Ch_Inv);
		KE3.Multx(qbar,dVec);
		for (i = 0; i<=2; ++i) {
			for (j = 0; j<=2; ++j) {
				KP2(i,j) = dVec[i]*Vvec[j];
			}
		}

		KP2 *= state[14];
		KP  += KP2;
		KP  /= -bott;
		KP  += I_m;
		KEP.MultAB(KE_Inv, KP);
*/
            	for (i = 0; i<=6; ++i){
			for (j = 0; j<=6; ++j){
                    		CMAT(i,j) = tmpVec[i]*Rvec[j];
                	}
            	}

            	DMAT.MultAB(CMAT, AA);
            	DMAT  /= -e;
            	DMAT  += AA;

            	// Calculate the consistent tangent
            	EMAT = 0.;
            	EMAT(0,0) = 1.;
            	EMAT(1,1) = 1.;
		EMAT(2,2) = 1.;

            	FMAT.MultAB(DMAT, EMAT);

            	for (i = 0; i <= 2; ++i){
                	for (j = 0; j<=2; ++j){
                    		KEP(i,j) = FMAT(i,j);
                	}
            	}



/*
	int i, j;

	dMatrixT KEP(3,3), KEE(3,3);

	dArrayT qn(4), dup(3), Sig(3), dQdSig(3), dfdq(4), qbar(4), dfdSig(3), H1(3), KEP1(3);

	double H2;

	fStiffness[1] = fStiffness[2] = fStiffness[3] = fStiffness[5] = fStiffness[6] = fStiffness[7] = 0.;
	Sig[0] = state[0];
	Sig[1] = state[1];
	Sig[2] = state[2];
	dup[0] = state[20];
	dup[1] = state[21];
	dup[2] = state[22];
	KEE = 0.;
	KEE(0,0) = fE_t1;
	KEE(1,1) = fE_t2;
	KEE(2,2) = fE_n;

	for (i = 0; i<=3; ++i) {
		qn[i] = state[i+9];
	}

	if (state[15] == 0.) {
		fStiffness[0] = fE_t1;
		fStiffness[4] = fE_t2;
		fStiffness[8] = fE_n;
	}
	else if (state[15] == 1.) {

		dfdSig_f(Sig, qn, dfdSig);
		dQdSig_f(Sig, qn, dQdSig);
		qbar_f(Sig, qn, dup, qbar);
		dfdq_f(Sig, qn, dfdq);

		for (i = 0; i <= 2; ++i){
			H1[i] = dfdSig[0]*KEE(0,i) + dfdSig[1]*KEE(1,i) + dfdSig[2]*KEE(2,i);
		}
		// H1[0] = dfdSig[0]*KEE(0,0) + dfdSig[1]*KEE(1,0) + dfdSig[2]*KEE(2,0);
		// H1[1] = dfdSig[0]*KEE(0,1) + dfdSig[1]*KEE(1,1) + dfdSig[2]*KEE(2,1);
		// H1[2] = dfdSig[0]*KEE(0,2) + dfdSig[1]*KEE(1,2) + dfdSig[2]*KEE(2,2);

		H2 = dArrayT::Dot(H1, dQdSig);
		H2 -= dArrayT::Dot(dfdq, qbar);

		for (i = 0; i <= 2; ++i){
			KEP1[i] = KEE(i,0)*dQdSig[0] + KEE(i,1)*dQdSig[1] + KEE(i,2)*dQdSig[2];
		}
		// KEP1[0] = KEE(0,0)*dQdSig[0] + KEE(0,1)*dQdSig[1] + KEE(0,2)*dQdSig[2];
		// KEP1[1] = KEE(1,0)*dQdSig[0] + KEE(1,1)*dQdSig[1] + KEE(1,2)*dQdSig[2];
		// KEP1[1] = KEE(2,0)*dQdSig[0] + KEE(2,1)*dQdSig[1] + KEE(2,2)*dQdSig[2];

		KEP.Outer(KEP1,H1);

		for (i = 0; i <= 2; ++i){
			for (j = 0; j<=2; ++j){
				KEP(i,j) = KEP(i,j)/H2;
				KEP(i,j) = KEE(i,j) - KEP(i,j);
			}
		}
*/
		fStiffness[0] = KEP(0,0);
		fStiffness[1] = KEP(0,1);
		fStiffness[2] = KEP(0,2);
		fStiffness[3] = KEP(1,0);
		fStiffness[4] = KEP(1,1);
		fStiffness[5] = KEP(1,2);
		fStiffness[6] = KEP(2,0);
		fStiffness[7] = KEP(2,1);
		fStiffness[8] = KEP(2,2);
	}

	return fStiffness;
}


/* print parameters to the output stream */
void MR3DT::Print(ostream& out) const
{
#ifndef _SIERRA_TEST_
	out << " Elastic tangential stiffness in 1 direction . . = " << fE_t1 << '\n';
	out << " Elastic tangential stiffness in 2 direction . . = " << fE_t2 << '\n';
	out << " Elastic Normal stiffness. . . . . . . . . . . . = " << fE_n     << '\n';
	out << " Mode_I Fracture Energy. . . . . . . . . . . . . = " << fGf_I     << '\n';
	out << " Mode_II Fracture Energy . . . . . . . . . . . . = " << fGf_II << '\n';
	out << " Peak Cohesion . . . . . . . . . . . . . . . . . = " << fc_p    << '\n';
	out << " Residual Cohesion . . . . . . . . . . . . . . . = " << fc_r    << '\n';
	out << " Peak Tensile Strength . . . . . . . . . . . . . = " << fchi_p << '\n';
	out << " Residual Tensile Strength . . . . . . . . . . . = " << fchi_r << '\n';
	out << " Peak Friction Angle . . . . . . . . . . . . . . = " << fphi_p   << '\n';
	out << " Critical State Friction Angle . . . . . . . . . = " << fphi_r  << '\n';
	out << " Peak Dilation Angle . . . . . . . . . . . . . . = " << fpsi_p   << '\n';
	out << " Coefficient of Tensile Strength Degradation . . = " << falpha_chi   << '\n';
	out << " Coefficient of Cohesion Degradation. .. . . . . = " << falpha_c   << '\n';
	out << " Coefficient for Frictional Angle Degradation. . = " << falpha_phi   << '\n';
	out << " Coefficient for Dilation Angle Degradation  . . = " << falpha_psi   << '\n';
	out << " Error Tolerance for Yield Function. . . . . . . = " << fTol_1 << '\n';
	out << " Error Tolerance for Residual. . . . . . . . . . = " << fTol_2 << '\n';
#endif
}

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default */
int MR3DT::NumOutputVariables(void) const { return 13; }

void MR3DT::OutputLabels(ArrayT<StringT>& labels) const
{
	labels.Dimension(13);
	labels[0] = "u_t1";
	labels[1] = "u_t2";
	labels[2] = "u_n";
	labels[3] = "up_t1";
	labels[4] = "up_t2";
	labels[5] = "up_n";
	labels[6] = "Chi";
	labels[7] = "Cohesion";
	labels[8] = "Friction Angle";
	labels[9] = "Dilation Angle";
	labels[10] = "Yield Function Value";
	labels[11] = "Norm of residuals";
	labels[12] = "No. of Iterations";
}

void MR3DT::ComputeOutput(const dArrayT& jump_u, const ArrayT<double>& state,
	dArrayT& output)
{
#pragma unused(jump_u)
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif	
	output[0] = state[3];
	output[1] = state[4];
	output[2] = state[5];
	output[3] = state[6];
	output[4] = state[7];
	output[5] = state[8];
	output[6] = state[9];
	output[7] = state[10];
	output[8] = atan(state[11]);
	output[9] = atan(state[12]);
	output[10] = state[13];
	output[11] = state[16];
	output[12] = state[19];
}


bool MR3DT::NeedsNodalInfo(void) { return false; }

int MR3DT::NodalQuantityNeeded(void) 
{ 
	return 2;
}

void MR3DT::SetElementGroupsNeeded(iArrayT& iGroups) 
{	
	iGroups[0] = 1;
}

double MR3DT::signof(double r)
{
	if (fabs(r) < kSmall)
		return 0.;
	else
		return fabs(r)/r;
}

