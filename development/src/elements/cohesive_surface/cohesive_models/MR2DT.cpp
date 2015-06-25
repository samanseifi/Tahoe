/*$Id: MR2DT.cpp,v 1.40 2011/12/01 20:38:01 beichuan Exp $*/
/* created by manzari*/
/* Elastolastic Cohesive Model for Geomaterials*/
#include "MR2DT.h"

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
const int knumDOF = 2;

/* constructor */
MR2DT::MR2DT(void): SurfacePotentialT(knumDOF)
{
	const char caller[] = "MR2DT::MR2DT";
	
	SetName("elastoplastic_MR_2D");
	
#if 0
	/* Elastic and Fracture Energy parameters */
	in >> fE_n; if (fE_n < 0) throw ExceptionT::kBadInputValue;
	in >> fE_t; if (fE_t < 0) throw ExceptionT::kBadInputValue;
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
int MR2DT::NumStateVariables(void) const { return 8*knumDOF +1; }

/* describe the parameters needed by the interface */
void MR2DT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SurfacePotentialT::DefineParameters(list);

	/* model parameters */
	ParameterT E_n(fE_n, "E_n");
	E_n.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(E_n);
	
	ParameterT E_t(fE_t, "E_t");
	E_t.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(E_t);
	
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
void MR2DT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SurfacePotentialT::DefineSubs(sub_list);
}

/* a pointer to the ParameterInterfaceT */
ParameterInterfaceT* MR2DT::NewSub(const StringT& name) const
{
	return SurfacePotentialT::NewSub(name);
}

/* accept parameter list */
void MR2DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SurfacePotentialT::TakeParameterList(list);
	
	/* extract parameters */
	fE_t = list.GetParameter("E_t");
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
	mr_ep_2d_out.open("mr_ep_2d.info");

}


/* initialize the state variable array */
void MR2DT::InitStateVariables(ArrayT<double>& state)
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
	
	double enp = state[14];
	double esp = state[15];
	double fchi = fchi_r + (fchi_p - fchi_r)*exp(-falpha_chi*enp);
	double fc   = fc_r + (fc_p - fc_r)*exp(-falpha_c*esp);
	double ftan_phi = tan(fphi_r) + (tan(fphi_p) - tan(fphi_r))*exp(-falpha_phi*esp);
	double ftan_psi = (tan(fpsi_p))*exp(-falpha_psi*esp);	
	state[6] = fchi;
	state[7] = fc ;
	state[8] = ftan_phi;
	state[9] = ftan_psi;
	state[13] = 0.;
}

/* Value of the Yield Function */ 
double MR2DT::YFValue(const ArrayT<double>& state) 
{
	return state[10]; 
}	

/** dissipated energy. Total amount of energy dissipated reaching
	 * the current state. */
double MR2DT::FractureEnergy(const ArrayT<double>& state)
{
	double FE = state[0]*state[4] + state[1]*state[5];
	return FE;
}

/** potential energy */
double MR2DT::Potential(const dArrayT& jump_u, const ArrayT<double>& state)
{
	double PT = (jump_u[0]*state[0] + jump_u[1]*state[1]);
	return PT;
}

/* surface status */
SurfacePotentialT::StatusT MR2DT::Status(const dArrayT& jump_u, 
	const ArrayT<double>& state)
{
	if (state[10]<0.){
		int Status = 0;
	}
	if (state[10]>0.) {
		int StatusT = 1;
	}
	if ((jump_u[0]*jump_u[0] + jump_u[1]*jump_u[1])>100000.) {
		int StatusT = 2;
	}
      
	/*return SurfacePotentialT::Precritical;*/
	return StatusT();
}

void MR2DT::PrintName(ostream& out) const
{
#ifndef _SIERRA_TEST_
	out << "    MR process zone 2D \n";
#endif
}

		
/* traction vector given displacement jump vector */
const dArrayT& MR2DT::Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma, bool qIntegrate)
{
#pragma unused(sigma)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kSizeMismatch;
#endif

	if (! qIntegrate)
	{
		fTraction[0] = state[0];
		fTraction[1] = state[1];
		return fTraction;
	}
	else
	{
		int i; int j; int kk;

		dMatrixT AA(6,6); dMatrixT KE(2,2); dMatrixT KE_Inv(2,2); dMatrixT I_mat(4,4); dMatrixT I(2,2);
		dMatrixT CMAT(6,6); dMatrixT A_qq(4,4); dMatrixT A_uu(2,2); dMatrixT A_uq(2,4);
		dMatrixT A_qu(4,2); dMatrixT ZMAT(2,4); dMatrixT ZMATP(4,2);
		dMatrixT dQdSig2(2,2); dMatrixT dqbardq(4,4); dMatrixT dQdSigdq(2,4); dMatrixT KE_dQdSig2(2,2); dMatrixT KE_dQdSigdq(2,4);
		dMatrixT dqbardSig(4,2); dMatrixT AA_inv(6,6);

		dArrayT u(2); dArrayT up(2); dArrayT du(2); dArrayT dup(2); dArrayT qn(4);
		dArrayT qo(4); dArrayT Rvec(6); dArrayT Cvec(6); dArrayT upo(2);
		dArrayT R(6); dArrayT Rmod(6); dArrayT Sig(2); dArrayT Sig_I(2); dArrayT dSig(2);
		dArrayT dQdSig(2); dArrayT dfdq(4); dArrayT qbar(4);
		dArrayT R2(6); dMatrixT X(6,1); dArrayT V_sig(2); dArrayT V_q(4);
		dArrayT dfdSig(2); dArrayT dq(4); dArrayT Y(6); dArrayT KE_dQdSig(2), KE_du(2);


		double ff; double ff0; double bott; double topp; double dlam; double dlam2; double normr; double normr0;

#if __option(extended_errorcheck)
                mr_ep_2d_out << setw(outputFileWidth) << "************* state_n *************" << endl;
                mr_ep_2d_out << setw(outputFileWidth) << "Tt = state_n[0] " << state[0] << endl;
                mr_ep_2d_out << setw(outputFileWidth) << "Tn = state_n[1] " << state[1] << endl;
                mr_ep_2d_out << setw(outputFileWidth) << "jump_u[0] = state_n[2] " << state[2] << endl;
                mr_ep_2d_out << setw(outputFileWidth) << "jump_u[1] = state_n[3] " << state[3] << endl;
                mr_ep_2d_out << setw(outputFileWidth) << "up_t = state_n[4] = " << state[4] << endl;
                mr_ep_2d_out << setw(outputFileWidth) << "up_n = state_n[5] = " << state[5] << endl;
                mr_ep_2d_out << setw(outputFileWidth) << "chi = state_n[6] = " << state[6] << endl;
                mr_ep_2d_out << setw(outputFileWidth) << "cohesion = state_n[7] = " << state[7] << endl;
                mr_ep_2d_out << setw(outputFileWidth) << "tan(phi) = state_n[8] = " << state[8] << endl;
                mr_ep_2d_out << setw(outputFileWidth) << "tan(psi) = state_n[9] = " << state[9] << endl;
                mr_ep_2d_out << setw(outputFileWidth) << "F = state_n[10] = " << state[10] << endl;
                mr_ep_2d_out << setw(outputFileWidth) << "dlam = state_n[11] = " << state[11] << endl;
                mr_ep_2d_out << setw(outputFileWidth) << "double(iplastic) = state_n[12] = " << state[12] << endl;
                mr_ep_2d_out << setw(outputFileWidth) << "normr = state_n[13] = " << state[13] << endl;
                mr_ep_2d_out << setw(outputFileWidth) << "kk = state_n[16] = " << state[16] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "*************end of state_n *************" << endl;
#endif


		int iplastic;
		dlam = 0.;
		dlam2 = 0.;
		normr = 0.;

		I_mat = 0.;
		I = 0.;

		KE = 0.;
		KE(0,0) = fE_t;
		KE(1,1) = fE_n;
		//ZMAT = 0.; 
		//ZMATP = 0.;
    
		KE_Inv.Inverse(KE);


		/* Calculate incremental jumps and initialize the neecessary vectors */
		for (i = 0; i<=1; ++i) {
			u[i] = jump_u[i];
			du[i] = u[i] - state[i+2];
			up[i] = state[i+4];
			upo[i] = up[i];
			Sig_I[i] = state[i];
			I(i,i) = 1.;
		}

		for (i = 0; i<=3; ++i) {
			qn[i] = state[i+6];
			qo[i] = qn[i];
			I_mat(i,i) = 1.;
		}

		dArrayT ue(2), Sig_tr(2);
		ue  = u;
		ue -= upo;
		KE.MultTx(ue, Sig_tr);		// trial traction
    
		/* Check the yield function */     
		Yield_f(Sig_tr, qn, ff);

		// Check for initial data
#if __option(extended_errorcheck)
		mr_ep_2d_out << setw(outputFileWidth) << "******check for initial data*****" << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "jump_u_t = " << u[0] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "jump_u_n = " << u[1] << endl;	
		mr_ep_2d_out << setw(outputFileWidth) << "T_t_tr = " << Sig_tr[0] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "T_n_tr = " << Sig_tr[1] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "T_t_I = " << Sig_I[0] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "T_n_I = " << Sig_I[1] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "up_t = " << up[0] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "up_n = " << up[1] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "upo_t = " << upo[0] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "upo_n = " << upo[1] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "ue_t = " << ue[0] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "ue_n = " << ue[1] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "qn[0] = " << qn[0] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "qn[1] = " << qn[1] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "qn[2] = " << qn[2] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "qn[3] = " << qn[3] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "qo[0] = " << qo[0] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "qo[1] = " << qo[1] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "qo[2] = " << qo[2] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "qo[3] = " << qo[3] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "F_tr = " << ff << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "******End of check*****" << endl;
#endif

		if (ff < 0.) {
			iplastic = 0;
			Sig = Sig_tr;
			normr = 0.;
			kk = 0;
		}
		else {
			iplastic = 1;
			kk = 0;
			Sig = Sig_I;		// T^0 = T_n

			Yield_f(Sig, qn, ff);
			if (fabs(ff) > 1e-12) {
				ff0 = ff;
			}
			else {
				ff0 = 1.;
			}

			dQdSig_f(Sig, qn, dQdSig);
			qbar_f(Sig, qn, qbar);
			KE.Multx(du, KE_du);
			KE.Multx(dQdSig, KE_dQdSig);

			for (i = 0; i<=1; ++i) {
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
			if (normr > 1e-10) {
				normr0 = normr;
			}
			else {
				normr0 = 1.;
			} 


			//while (ff > fTol_1 || normr > fTol_2) {
			while ((fabs(ff/ff0) > fTol_1 || (normr/normr0) > fTol_2) && (fabs(ff) > fTol_1 || normr > fTol_2)) {

				// Check for data of local iteration
#if __option(extended_errorcheck)
				mr_ep_2d_out << setw(outputFileWidth) << "******check for data of local iteration*****" << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "local iteration # = " << kk << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "jump_u_t = " << u[0] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "jump_u_n = " << u[1] << endl;				
				mr_ep_2d_out << setw(outputFileWidth) << "T_t = " << Sig[0] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "T_n = " << Sig[1] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "up_t = " << up[0] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "up_n = " << up[1] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "ue_t = " << ue[0] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "ue_n = " << ue[1] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "qn[0] = " << qn[0] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "qn[1] = " << qn[1] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "qn[2] = " << qn[2] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "qn[3] = " << qn[3] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dlam = " << dlam << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "F = " << ff << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "norm = " << normr << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "******End of check*****" << endl;
#endif

				if (kk > 50) {
					ExceptionT::GeneralFail("MR2DT::Traction","Too Many Iterations, 50");
				}

				dQdSig_f(Sig, qn, dQdSig);
				qbar_f(Sig, qn, qbar);
				KE.Multx(dQdSig, KE_dQdSig);

#if __option(extended_errorcheck)
				mr_ep_2d_out << setw(outputFileWidth) << "******check for data*****" << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dQdSig[0] = " << dQdSig[0] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dQdSig[1] = " << dQdSig[1] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "qbar[0] = " << qbar[0] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "qbar[1] = " << qbar[1] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "qbar[2] = " << qbar[2] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "qbar[3] = " << qbar[3] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "KE_dQdSig[0] = " << KE_dQdSig[0] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "KE_dQdSig[1] = " << KE_dQdSig[1] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "******End of check*****" << endl;
#endif

				for (i = 0; i<=1; ++i) {
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

				// Check for the residual
#if __option(extended_errorcheck)
				mr_ep_2d_out << setw(outputFileWidth) << "******check for the residual*****" << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "R[0] = " << R[0] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "R[1] = " << R[1] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "R[2] = " << R[2] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "R[3] = " << R[3] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "R[4] = " << R[4] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "R[5] = " << R[5] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "******End of check*****" << endl;
#endif

				dQdSig2_f(Sig,qn,dQdSig2);
				dQdSigdq_f(Sig, qn, dQdSigdq);
				dqbardSig_f(Sig, qn, dqbardSig);
				dqbardq_f(Sig, qn, dqbardq);

				KE_dQdSig2.MultAB(KE, dQdSig2);
				KE_dQdSigdq.MultAB(KE, dQdSigdq);

#if __option(extended_errorcheck)
				mr_ep_2d_out << setw(outputFileWidth) << "******check for data*****" << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dQdSig2(0,0) = " << dQdSig2(0,0) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dQdSig2(0,1) = " << dQdSig2(0,1) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dQdSig2(1,0) = " << dQdSig2(1,0) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dQdSig2(1,1) = " << dQdSig2(1,1) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dQdSigdq(0,0) = " << dQdSigdq(0,0) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dQdSigdq(0,1) = " << dQdSigdq(0,1) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dQdSigdq(0,2) = " << dQdSigdq(0,2) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dQdSigdq(0,3) = " << dQdSigdq(0,3) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dQdSigdq(1,0) = " << dQdSigdq(1,0) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dQdSigdq(1,1) = " << dQdSigdq(1,1) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dQdSigdq(1,2) = " << dQdSigdq(1,2) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dQdSigdq(1,3) = " << dQdSigdq(1,3) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dqbardSig(0,0) = " << dqbardSig(0,0) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dqbardSig(0,1) = " << dqbardSig(0,1) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dqbardSig(1,0) = " << dqbardSig(1,0) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dqbardSig(1,1) = " << dqbardSig(1,1) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dqbardSig(2,0) = " << dqbardSig(2,0) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dqbardSig(2,1) = " << dqbardSig(2,1) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dqbardSig(3,0) = " << dqbardSig(3,0) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dqbardSig(3,1) = " << dqbardSig(3,1) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dqbardq(0,0) = " << dqbardq(0,0) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dqbardq(0,1) = " << dqbardq(0,1) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dqbardq(0,2) = " << dqbardq(0,2) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dqbardq(0,3) = " << dqbardq(0,3) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dqbardq(1,0) = " << dqbardq(1,0) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dqbardq(1,1) = " << dqbardq(1,1) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dqbardq(1,2) = " << dqbardq(1,2) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dqbardq(1,3) = " << dqbardq(1,3) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dqbardq(2,0) = " << dqbardq(2,0) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dqbardq(2,1) = " << dqbardq(2,1) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dqbardq(2,2) = " << dqbardq(2,2) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dqbardq(2,3) = " << dqbardq(2,3) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dqbardq(3,0) = " << dqbardq(3,0) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dqbardq(3,1) = " << dqbardq(3,1) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dqbardq(3,2) = " << dqbardq(3,2) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dqbardq(3,3) = " << dqbardq(3,3) << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "******End of check*****" << endl;
#endif

				for (i = 0; i<=5; ++i) {
					for (j = 0; j<=5; ++j) {
						if (i<=1 & j<=1){
							AA_inv(i,j)  = I(i,j);
							AA_inv(i,j) += dlam*KE_dQdSig2(i,j);
						}
						if (i<=1 & j>1){
							AA_inv(i,j)  = KE_dQdSigdq(i,j-2);
							AA_inv(i,j) *= dlam;
						}
						if(i>1 & j<=1){
							AA_inv(i,j)  = dqbardSig(i-2,j);
							AA_inv(i,j) *= dlam;
						}
						if(i>1 & j >1) {
							AA_inv(i,j)  = I_mat(i-2,j-2);
							AA_inv(i,j) *= -1.;
							AA_inv(i,j) += dlam*dqbardq(i-2,j-2);
						}
					}
				}

				AA.Inverse(AA_inv);
				dfdSig_f(Sig, qn, dfdSig);
				dfdq_f(Sig, qn, dfdq);

#if __option(extended_errorcheck)
				mr_ep_2d_out << setw(outputFileWidth) << "******check for data*****" << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dfdSig[0] = " << dfdSig[0] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dfdSig[1] = " << dfdSig[1] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dfdq[0] = " << dfdq[0] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dfdq[1] = " << dfdq[1] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dfdq[2] = " << dfdq[2] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "dfdq[3] = " << dfdq[3] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "******End of check*****" << endl;
#endif


				for (i = 0; i<=5; ++i) {
					if (i<=1){
						Rvec[i] = dfdSig[i];
						Cvec[i] = KE_dQdSig[i];
					}
					if (i > 1){
						Rvec[i] = dfdq[i-2];
						Cvec[i] = qbar[i-2];
					}
				}

				dArrayT tmpRvec(6), tmpCvec(6);
				AA.Multx(R,tmpRvec);
				topp  = ff;
				topp -= dArrayT::Dot(Rvec,tmpRvec);
				AA.Multx(Cvec,tmpCvec);
				bott  = dArrayT::Dot(Rvec,tmpCvec);
				dlam2 = topp/bott;

				/*
				for (i = 0; i<=5; ++i) {
					for (j = 0; j<=5; ++j) {
						if (i<=1 & j<=1){
							CMAT(i,j) = KE_Inv(i,j);
						}
						if (i<=1 & j>1) {
							CMAT(i,j) = 0.;
						}
						if(i>1 & j<=1) {
							CMAT(i,j) = 0.;
						}
						if(i>1 & j >1) {
							CMAT(i,j) = -I_mat(i-2,j-2);
						}
					}
				}
				*/

				for (i = 0; i<=5; ++i) {
					if (i<=1){
						Rmod[i] = KE_dQdSig[i];
					}
					if (i >1){
						Rmod[i] = qbar[i-2];
					}
				}

				Rmod *= dlam2;
				R2  = R;
				R2 += Rmod;
				AA.Multx(R2,X);
				
				for (i = 0; i<=5; ++i) {
					if (i<=1) {
						//dup[i] = Y[i];
						dSig[i] = -X[i];
					}
					if (i > 1) {
						//dq[i-2] = Y[i];
						dq[i-2] = -X[i];
					}
				}

				// Update the local variables
				Sig += dSig;
				qn += dq;
				dlam = dlam + dlam2;
				kk = kk + 1;

				// Check the yield function and norm of R
				Yield_f(Sig, qn, ff);
				dQdSig_f(Sig, qn, dQdSig);
				qbar_f(Sig, qn, qbar);
                                KE.Multx(dQdSig, KE_dQdSig);

				for (i = 0; i<=1; ++i) {
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

#if __option(extended_errorcheck)
                                mr_ep_2d_out << setw(outputFileWidth) << "******check for data*****" << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "Sig[0] = " << Sig[0] << endl;
				mr_ep_2d_out << setw(outputFileWidth) << "Sig[1] = " << Sig[1] << endl;
                                mr_ep_2d_out << setw(outputFileWidth) << "qn[0] = " << qn[0] << endl;
                                mr_ep_2d_out << setw(outputFileWidth) << "qn[1] = " << qn[1] << endl;
                                mr_ep_2d_out << setw(outputFileWidth) << "qn[2] = " << qn[2] << endl;
                                mr_ep_2d_out << setw(outputFileWidth) << "qn[3] = " << qn[3] << endl;
                                mr_ep_2d_out << setw(outputFileWidth) << "dlam = " << dlam << endl;
                                mr_ep_2d_out << setw(outputFileWidth) << "F = " << ff << endl;
                                mr_ep_2d_out << setw(outputFileWidth) << "dQdSig[0] = " << dQdSig[0] << endl;
                                mr_ep_2d_out << setw(outputFileWidth) << "dQdSig[1] = " << dQdSig[1] << endl;
                                mr_ep_2d_out << setw(outputFileWidth) << "qbar[0] = " << qbar[0] << endl;
                                mr_ep_2d_out << setw(outputFileWidth) << "qbar[1] = " << qbar[1] << endl;
                                mr_ep_2d_out << setw(outputFileWidth) << "qbar[2] = " << qbar[2] << endl;
                                mr_ep_2d_out << setw(outputFileWidth) << "qbar[3] = " << qbar[3] << endl;
                                mr_ep_2d_out << setw(outputFileWidth) << "KE_dQdSig[0] = " << KE_dQdSig[0] << endl;
                                mr_ep_2d_out << setw(outputFileWidth) << "KE_dQdSig[1] = " << KE_dQdSig[1] << endl;
                                mr_ep_2d_out << setw(outputFileWidth) << "******End of check*****" << endl;

                                mr_ep_2d_out << setw(outputFileWidth) << "******check for the residual for next local iteration*****" << endl;
                                mr_ep_2d_out << setw(outputFileWidth) << "R[0] = " << R[0] << endl;
                                mr_ep_2d_out << setw(outputFileWidth) << "R[1] = " << R[1] << endl;
                                mr_ep_2d_out << setw(outputFileWidth) << "R[2] = " << R[2] << endl;
                                mr_ep_2d_out << setw(outputFileWidth) << "R[3] = " << R[3] << endl;
                                mr_ep_2d_out << setw(outputFileWidth) << "R[4] = " << R[4] << endl;
                                mr_ep_2d_out << setw(outputFileWidth) << "R[5] = " << R[5] << endl;
                                mr_ep_2d_out << setw(outputFileWidth) << "******End of check*****" << endl;
#endif

			}

			// Update plastic jump displacement
			for (i = 0; i<=1; ++i) {
				up[i]  = upo[i];
				up[i] += dlam*dQdSig[i];
			}
		}

		state[0] = Sig[0];
		state[1] = Sig[1];
		fTraction[0] = state[0];
		fTraction[1] = state[1];
		state[2] = jump_u[0];
		state[3] = jump_u[1];
		state[4] = up[0];
		state[5] = up[1];
		state[6] = qn[0];
		state[7] = qn[1];
		state[8] = qn[2];
		state[9] = qn[3];
		state[10] = ff;
		state[11] = dlam;
		state[12] = double(iplastic);
		state[13] = normr;
		dQdSig_f(Sig, qn, dQdSig);
		state[14] = Sig[0]*dQdSig[0];
		state[14] += (Sig[1] + fabs(Sig[1]))*dQdSig[1]/2.;
		state[14] /=fGf_I;
		state[14] *=dlam;
		state[15]  = signof(Sig[0]);
		state[15] -= signof(Sig[0])*fabs(Sig[1]*qn[2]);
		state[15] *= dQdSig[0];
		state[15] /= fGf_II;
		state[15] *=dlam;
		state[16] = double(kk);

		// Check for yield and norm_R after convergence is achieved
#if __option(extended_errorcheck)
		mr_ep_2d_out << setw(outputFileWidth) << "*************Results are converged*************" << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "Tt = state[0] " << state[0] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "Tn = state[1] " << state[1] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "jump_u[0] = state[2] " << state[2] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "jump_u[1] = state[3] " << state[3] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "up_t = state[4] = " << state[4] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "up_n = state[5] = " << state[5] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "chi = state[6] = " << state[6] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "cohesion = state[7] = " << state[7] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "tan(phi) = state[8] = " << state[8] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "tan(psi) = state[9] = " << state[9] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "F = state[10] = " << state[10] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "dlam = state[11] = " << state[11] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "double(iplastic) = state[12] = " << state[12] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "normr = state[13] = " << state[13] << endl;
		mr_ep_2d_out << setw(outputFileWidth) << "kk = state[16] = " << state[16] << endl;

		// Check for the stiffness after convergence is achieved
		Stiffness(jump_u, state, sigma);
		mr_ep_2d_out << setw(outputFileWidth) << "KEP(0,0) = " << fStiffness[0]	<< endl;
		mr_ep_2d_out << setw(outputFileWidth) << "KEP(0,1) = " << fStiffness[2]	<< endl;
		mr_ep_2d_out << setw(outputFileWidth) << "KEP(1,0) = " << fStiffness[1]	<< endl;
		mr_ep_2d_out << setw(outputFileWidth) << "KEP(1,1) = " << fStiffness[3]	<< endl;
		mr_ep_2d_out << setw(outputFileWidth) << "*************End of results*************" << endl;
#endif

		return fTraction;
	}
}

/* calculation of Yield_f */
double& MR2DT::Yield_f(const dArrayT& Sig, const dArrayT& qn, double& ff)
{
	double tmp1, tmp11, tmp22, tmp3, tmp31, tmp32, tmp4, tmp5;

	tmp1   = qn[1];
	tmp11  = Sig[1];
	tmp11 *= qn[2];
	tmp1  -= tmp11;

	tmp22  = Sig[0];
	tmp22 *= Sig[0];

	tmp3   = qn[1];
	tmp31  = qn[0];
	tmp31 *= qn[2];
	tmp3  -= tmp31;
	tmp32  = tmp3;
	tmp32 *= tmp3;

	tmp4  = tmp22;
	tmp4 += tmp32;

	tmp5 = sqrt(tmp4);

	ff  = tmp5;
	ff -= tmp1;

	return ff;
}


/* calculation of qbar_f */
dArrayT& MR2DT::qbar_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& qbar)
{
/*
	double A1 = -falpha_chi*(qn[0] - fchi_r);
	double B1 = (Sig[1] + fabs(Sig[1]))/2./fGf_I;
	double B2 = Sig[0]/fGf_I;
	double DQDN = 2.*qn[3]*(qn[1] - Sig[1]*qn[3]);
	double DQDT = 2.*Sig[0];
	double A2 = -falpha_c*(qn[1] - fc_r);
	double TNA = (Sig[1] - fabs(Sig[1]))/2.;
	double B3 = (Sig[0] - fabs(TNA*qn[2])*signof(Sig[0]))/fGf_II;
	double A3 = -falpha_phi*(qn[2] - tan(fphi_r));
	double A4 = -falpha_psi*qn[3];
*/
/*
	double A1 = -falpha_chi*(qn[0] - fchi_r);
	double B1 = (Sig[1] + fabs(Sig[1]))/(2.*fGf_I);
	double B2 = Sig[0]/fGf_I;
	double Shear_Q = Sig[0]*Sig[0] + (qn[1] - qn[0]*qn[3])*(qn[1] - qn[0]*qn[3]);
	double DQDN = qn[3];
	double DQDT = Sig[0]/sqrt(Shear_Q);
	double A2 = -falpha_c*(qn[1] - fc_r);
	double TNA = (Sig[1] - fabs(Sig[1]))/2.;
	double B3 = (Sig[0] - fabs(TNA*qn[2])*signof(Sig[0]))/fGf_II;
	double A3 = -falpha_phi*(qn[2] - tan(fphi_r));
	double A4 = -falpha_psi*qn[3];

	qbar[0] = A1*B1*DQDN + A1*B2*DQDT;
	qbar[1] = A2*B3*DQDT;
	qbar[2] = A3*B3*DQDT;
	qbar[3] = A4*B3*DQDT;
*/

	double A1 = -falpha_chi*(qn[0] - fchi_r);
	double A2 = -falpha_chi*(qn[0] - fchi_r);
	double A3 = -falpha_c*(qn[1] - fc_r);
	double A4 = -falpha_c*(qn[1] - fc_r);
	double A6 = -falpha_phi*(qn[2] - tan(fphi_r));
	double A8 = -falpha_psi*qn[3];

	double B1 = (Sig[1] + fabs(Sig[1]))/(2.*fGf_I);
	double TNA = (Sig[1] - fabs(Sig[1]))/2.;
	//double B3 = (Sig[0] - fabs(TNA*qn[2])*signof(Sig[0]))/fGf_II;
	double Tt_excess = fabs(Sig[0]) - fabs(TNA*qn[2]);
	double B3func = 0.5*(Tt_excess + fabs(Tt_excess));
	double B3 = B3func*signof(Sig[0])/fGf_II;
	
	double Shear_Q = Sig[0]*Sig[0] + (qn[1] - qn[0]*qn[3])*(qn[1] - qn[0]*qn[3]);
	double DQDN = qn[3];
	double DQDT = Sig[0]/sqrt(Shear_Q);

	qbar[0] = A1*B1*DQDN + A2*B3*DQDT;
	qbar[1] = A3*B1*DQDN + A4*B3*DQDT;
	qbar[2] = A6*B3*DQDT;
	qbar[3] = A8*B3*DQDT;

	return qbar;
 }


/* calculation of dQdSig2_f */
dMatrixT& MR2DT::dQdSig2_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dQdSig2)
{
#pragma unused(Sig)
/*
	dQdSig2(0,0) = 2.;
	dQdSig2(1,1) = -2.*qn[3]*qn[3];
	dQdSig2(0,1) = 0.;
	dQdSig2(1,0) = 0.;
*/
	double Shear_Q = Sig[0]*Sig[0] + (qn[1] - qn[0]*qn[3])*(qn[1] - qn[0]*qn[3]);

	dQdSig2(0,0) = (qn[1] - qn[0]*qn[3])*(qn[1] - qn[0]*qn[3])/sqrt(Shear_Q*Shear_Q*Shear_Q);
	dQdSig2(0,1) = 0.;
	dQdSig2(1,0) = 0.;
	dQdSig2(1,1) = 0.;

	return dQdSig2;
}

/* calculation of dfdSig_f */
dArrayT& MR2DT::dfdSig_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& dfdSig)
{
	double Shear = Sig[0]*Sig[0] + (qn[1] - qn[0]*qn[2])*(qn[1] - qn[0]*qn[2]);

	dfdSig[0] = Sig[0]/sqrt(Shear);
	dfdSig[1] = qn[2];

	return dfdSig;
}

/* calculation of dQdSig_f */
dArrayT& MR2DT::dQdSig_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& dQdSig)
{
/*
	dQdSig[0] = 2.*Sig[0];
	dQdSig[1] = 2.*qn[3]*(qn[1] - Sig[1]*qn[3]);
*/
	double Shear_Q = Sig[0]*Sig[0] + (qn[1] - qn[0]*qn[3])*(qn[1] - qn[0]*qn[3]);

	dQdSig[0] = Sig[0]/sqrt(Shear_Q);
	dQdSig[1] = qn[3];

	return dQdSig;
}


/* calculation of dfdq_f */
dArrayT& MR2DT::dfdq_f(const dArrayT& Sig, const dArrayT& qn, dArrayT& dfdq)
{
	double Shear = Sig[0]*Sig[0] + (qn[1] - qn[0]*qn[2])*(qn[1] - qn[0]*qn[2]);
	double zeta = (qn[1] - qn[0]*qn[2])/sqrt(Shear);

	dfdq[0] = -qn[2]*zeta;
	dfdq[1] = zeta - 1.;
	dfdq[2] = -qn[0]*zeta + Sig[1];
	dfdq[3] = 0.;

	return dfdq;
}

/* calculation of dQdSigdq_f */
dMatrixT& MR2DT::dQdSigdq_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dQdSigdq)
{
/*
	dQdSigdq(0,0) = 0.;
	dQdSigdq(0,1) = 0.;
	dQdSigdq(0,2) = 0.;
	dQdSigdq(0,3) = 0.;
	dQdSigdq(1,0) = 0.;
	dQdSigdq(1,1) = 2.*qn[3];
	dQdSigdq(1,2) = 0.;
	dQdSigdq(1,3) = 2.*qn[1] - 4.*Sig[1]*qn[3];
*/
	double Shear_Q = Sig[0]*Sig[0] + (qn[1] - qn[0]*qn[3])*(qn[1] - qn[0]*qn[3]);
	double zeta_Q = (qn[1] - qn[0]*qn[3])/sqrt(Shear_Q*Shear_Q*Shear_Q);

	dQdSigdq(0,0) = Sig[0]*qn[3]*zeta_Q;
	dQdSigdq(0,1) = -Sig[0]*zeta_Q;
	dQdSigdq(0,2) = 0.;
	dQdSigdq(0,3) = Sig[0]*qn[0]*zeta_Q;
	dQdSigdq(1,0) = 0.;
	dQdSigdq(1,1) = 0.;
	dQdSigdq(1,2) = 0.;
	dQdSigdq(1,3) = 1.;

	return dQdSigdq;
}

/* calculation of dqbardSig_f */
dMatrixT& MR2DT::dqbardSig_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dqbardSig)
{
/*
	double A1 = -falpha_chi*(qn[0] - fchi_r);
	double B1 = (Sig[1]+fabs(Sig[1]))/2./fGf_I;
	double B2 = Sig[0]/fGf_I;
	double DQDN = 2.*qn[3]*(qn[1] - Sig[1]*qn[3]);
	double DQDT = 2.*Sig[0];
	double A2 = -falpha_c*(qn[1] - fc_r);
	double TNA = (Sig[1]-fabs(Sig[1]))/2.;
	double B3 = (Sig[0] - fabs(TNA*qn[2])*signof(Sig[0]))/fGf_II;
	double A3 = -falpha_phi*(qn[2] - tan(fphi_r));
	double A4 = -falpha_psi*qn[3];
	double DB3_DTn = -qn[2]*signof(Sig[0])*signof(TNA)*(1. - signof(Sig[1]))/fGf_II/2.;
	double DB3_DTt = 1./fGf_II;
	double DB3_DTanphi = -fabs(TNA)*signof(Sig[0])/fGf_II;
	double DQDN2 = -2.*qn[3]*qn[3];
	double DQDT2 = 2.;
	double DQDTN = 0.;
	double DQDNT = 0.;
	double SN = signof(Sig[1]);
	double DB1DN = (SN +fabs(SN))/2./fGf_I;

	dqbardSig(0,0) = A1*B2*DQDT2 + A1*DQDT/fGf_I;
	dqbardSig(0,1) = A1*B1*DQDN2 + A1*DQDN*DB1DN;
	dqbardSig(1,0) = A2*B3*DQDT2 + A2*DQDT*DB3_DTt;
	dqbardSig(1,1) = A2*DQDT*DB3_DTn;
	dqbardSig(2,0) = A3*B3*DQDT2 + A3*DQDT*DB3_DTt;
	dqbardSig(2,1) = A3*DQDT*DB3_DTn;
	dqbardSig(3,0) = A4*B3*DQDT2 + A4*DQDT*DB3_DTt;
	dqbardSig(3,1) = A4*DQDT*DB3_DTn;
*/
/*
	double A1 = -falpha_chi*(qn[0] - fchi_r);
	double B1 = (Sig[1] + fabs(Sig[1]))/(2.*fGf_I);
	double B2 = Sig[0]/fGf_I;
	double Shear_Q = Sig[0]*Sig[0] + (qn[1] - qn[0]*qn[3])*(qn[1] - qn[0]*qn[3]);
	double DQDN = qn[3];
	double DQDT = Sig[0]/sqrt(Shear_Q);
	double A2 = -falpha_c*(qn[1] - fc_r);
	double TNA = (Sig[1] - fabs(Sig[1]))/2.;
	double B3 = (Sig[0] - fabs(TNA*qn[2])*signof(Sig[0]))/fGf_II;
	double A3 = -falpha_phi*(qn[2] - tan(fphi_r));
	double A4 = -falpha_psi*qn[3];
	double DB3_DTn = -qn[2]*signof(Sig[0])*signof(TNA)*(1. - signof(Sig[1]))/(2.*fGf_II);
	double DB3_DTt = 1./fGf_II;
	double DQDN2 = 0.;
	double DQDT2 = (qn[1] - qn[0]*qn[3])*(qn[1] - qn[0]*qn[3])/sqrt(Shear_Q*Shear_Q*Shear_Q);
	double SN = signof(Sig[1]);
	double DB1DN = (SN +fabs(SN))/(2.*fGf_I);

	dqbardSig(0,0) = A1*B2*DQDT2 + A1*DQDT/fGf_I;
	dqbardSig(0,1) = A1*B1*DQDN2 + A1*DQDN*DB1DN;
	dqbardSig(1,0) = A2*B3*DQDT2 + A2*DQDT*DB3_DTt;
	dqbardSig(1,1) = A2*DQDT*DB3_DTn;
	dqbardSig(2,0) = A3*B3*DQDT2 + A3*DQDT*DB3_DTt;
	dqbardSig(2,1) = A3*DQDT*DB3_DTn;
	dqbardSig(3,0) = A4*B3*DQDT2 + A4*DQDT*DB3_DTt;
	dqbardSig(3,1) = A4*DQDT*DB3_DTn;
*/

	double A1 = -falpha_chi*(qn[0] - fchi_r);
	double A2 = -falpha_chi*(qn[0] - fchi_r);
	double A3 = -falpha_c*(qn[1] - fc_r);
	double A4 = -falpha_c*(qn[1] - fc_r);
	double A6 = -falpha_phi*(qn[2] - tan(fphi_r));
	double A8 = -falpha_psi*qn[3];

	double B1 = (Sig[1] + fabs(Sig[1]))/(2.*fGf_I);
	double TNA = (Sig[1] - fabs(Sig[1]))/2.;
	//double B3 = (Sig[0] - fabs(TNA*qn[2])*signof(Sig[0]))/fGf_II;
	double Tt_excess = fabs(Sig[0]) - fabs(TNA*qn[2]);
	double B3func = 0.5*(Tt_excess + fabs(Tt_excess));
	double B3 = B3func*signof(Sig[0])/fGf_II;
		
	double Shear_Q = Sig[0]*Sig[0] + (qn[1] - qn[0]*qn[3])*(qn[1] - qn[0]*qn[3]);
	double DQDN = qn[3];
	double DQDT = Sig[0]/sqrt(Shear_Q);

	double signfun = (1. + signof(Tt_excess))/2.;
	//double DB3_DTn = -qn[2]*signof(Sig[0])*signof(TNA)*(1. - signof(Sig[1]))/(2.*fGf_II);
	double DB3_DTn = -qn[2]*signof(Sig[0])*signof(TNA)*(1. - signof(Sig[1]))*signfun/(2.*fGf_II);
	//double DB3_DTt = 1./fGf_II;
	double DB3_DTt = 1.*signfun/fGf_II;
	double DQDN2 = 0.;
	double DQDT2 = (qn[1] - qn[0]*qn[3])*(qn[1] - qn[0]*qn[3])/sqrt(Shear_Q*Shear_Q*Shear_Q);
	double SN = signof(Sig[1]);
	double DB1_DTn = (1. + SN)/(2.*fGf_I);

	dqbardSig(0,0) = A2*DB3_DTt*DQDT + A2*B3*DQDT2;
	dqbardSig(0,1) = A1*DB1_DTn*DQDN + A1*B1*DQDN2 + A2*DB3_DTn*DQDT;
	dqbardSig(1,0) = A4*DB3_DTt*DQDT + A4*B3*DQDT2;
	dqbardSig(1,1) = A3*DB1_DTn*DQDN + A3*B1*DQDN2 + A4*DB3_DTn*DQDT;
	dqbardSig(2,0) = A6*DB3_DTt*DQDT + A6*B3*DQDT2;
	dqbardSig(2,1) = A6*DB3_DTn*DQDT;
	dqbardSig(3,0) = A8*DB3_DTt*DQDT + A8*B3*DQDT2;
	dqbardSig(3,1) = A8*DB3_DTn*DQDT;

	return dqbardSig;
}
  
/* calculation of dqbardq_f */
dMatrixT& MR2DT::dqbardq_f(const dArrayT& Sig, const dArrayT& qn, dMatrixT& dqbardq)
{
/*
	double A1 = -falpha_chi*(qn[0] - fchi_r);
	double B1 = (Sig[1]+fabs(Sig[1]))/2./fGf_I;
	double B2 = Sig[0]/fGf_I;
	double DQDN = 2.*qn[3]*(qn[1] - Sig[1]*qn[3]);
	double DQDT = 2.*Sig[0];
	double A2 = -falpha_c*(qn[1] - fc_r);
	double TNA = (Sig[1]-fabs(Sig[1]))/2.;
	double B3 = (Sig[0] - fabs(TNA*qn[2])*signof(Sig[0]))/fGf_II;
	double A3 = -falpha_phi*(qn[2] - tan(fphi_r));
	double A4 = -falpha_psi*qn[3];
	double DB3_DTn = -qn[2]*signof(Sig[0])*signof(TNA)*(1. - signof(Sig[1]))/fGf_II/2.;
	double DB3_DTt = 1./fGf_II;
	double DB3_DTanphi = -fabs(TNA)*signof(Sig[0])/fGf_II;
	double DQDN2 = -2.*qn[3]*qn[3];
	double DQDT2 = 2.;
	double DQDTN = 0.;
	double DQDNT = 0.;
	double SN = signof(Sig[1]);
	double DB1DN = (SN +fabs(SN))/2./fGf_I;

	dqbardq(0,0) = -falpha_chi*(B1*DQDN + B2*DQDT);
	dqbardq(0,1) =  A1*B1*(2.*qn[3]);
	dqbardq(0,2) = 0.;
	dqbardq(0,3) =  A1*B1*(2.*qn[1]-4.*Sig[1]*qn[3]);
	dqbardq(1,0) = 0.;
	dqbardq(1,1) = -falpha_c*B3*DQDT;
	dqbardq(1,2) = A2*DQDT*DB3_DTanphi;
	dqbardq(1,3) = 0.;
	dqbardq(2,0) = 0.;
	dqbardq(2,1) = 0.;
	dqbardq(2,2) = -falpha_phi*B3*DQDT + A3*DQDT*DB3_DTanphi;
	dqbardq(2,3) = 0.;
	dqbardq(3,0) = 0.;
	dqbardq(3,1) = 0.;
	dqbardq(3,2) = A4*DQDT*DB3_DTanphi;
	dqbardq(3,3) = -falpha_psi*B3*DQDT;
*/
/*
	double A1 = -falpha_chi*(qn[0] - fchi_r);
	double B1 = (Sig[1] + fabs(Sig[1]))/(2.*fGf_I);
	double B2 = Sig[0]/fGf_I;
	double Shear_Q = Sig[0]*Sig[0] + (qn[1] - qn[0]*qn[3])*(qn[1] - qn[0]*qn[3]);
	double zeta_Q = (qn[1] - qn[0]*qn[3])/sqrt(Shear_Q*Shear_Q*Shear_Q);
	double DQDN = qn[3];
	double DQDT = Sig[0]/sqrt(Shear_Q);
	double A2 = -falpha_c*(qn[1] - fc_r);
	double TNA = (Sig[1] - fabs(Sig[1]))/2.;
	double B3 = (Sig[0] - fabs(TNA*qn[2])*signof(Sig[0]))/fGf_II;
	double A3 = -falpha_phi*(qn[2] - tan(fphi_r));
	double A4 = -falpha_psi*qn[3];
	double DB3_DTanphi = -fabs(TNA)*signof(Sig[0])/fGf_II;
	double DQDTDChi = Sig[0]*qn[3]*zeta_Q;
	double DQDTDC = -Sig[0]*zeta_Q;
	double DQDTDTanpsi = Sig[0]*qn[0]*zeta_Q;
	double DQDNDTanpsi = 1.;


	dqbardq(0,0) = -falpha_chi*(B1*DQDN + B2*DQDT) + A1*B2*DQDTDChi;
	dqbardq(0,1) = A1*B1*DQDTDC;
	dqbardq(0,2) = 0.;
	dqbardq(0,3) = A1*B1*DQDNDTanpsi + A1*B2*DQDTDTanpsi;
	dqbardq(1,0) = A2*B3*DQDTDChi;
	dqbardq(1,1) = -falpha_c*B3*DQDT + A2*B3*DQDTDC;
	dqbardq(1,2) = A2*DQDT*DB3_DTanphi;
	dqbardq(1,3) = A2*B3*DQDTDTanpsi;
	dqbardq(2,0) = A3*B3*DQDTDChi;
	dqbardq(2,1) = A3*B3*DQDTDC;
	dqbardq(2,2) = -falpha_phi*B3*DQDT + A3*DQDT*DB3_DTanphi;
	dqbardq(2,3) = A3*B3*DQDTDTanpsi;
	dqbardq(3,0) = A4*B3*DQDTDChi;
	dqbardq(3,1) = A4*B3*DQDTDC;
	dqbardq(3,2) = A4*DQDT*DB3_DTanphi;
	dqbardq(3,3) = -falpha_psi*B3*DQDT + A4*B3*DQDTDTanpsi;
*/

	double A1 = -falpha_chi*(qn[0] - fchi_r);
	double A2 = -falpha_chi*(qn[0] - fchi_r);
	double A3 = -falpha_c*(qn[1] - fc_r);
	double A4 = -falpha_c*(qn[1] - fc_r);
	double A6 = -falpha_phi*(qn[2] - tan(fphi_r));
	double A8 = -falpha_psi*qn[3];

	double DA1_DChi = -falpha_chi;
	double DA2_DChi = -falpha_chi;
	double DA3_DC = -falpha_c;
	double DA4_DC = -falpha_c;
	double DA6_DTanphi = -falpha_phi;
	double DA8_DTanpsi = -falpha_psi;

	double B1 = (Sig[1] + fabs(Sig[1]))/(2.*fGf_I);
	double TNA = (Sig[1] - fabs(Sig[1]))/2.;
	//double B3 = (Sig[0] - fabs(TNA*qn[2])*signof(Sig[0]))/fGf_II;
	double Tt_excess = fabs(Sig[0]) - fabs(TNA*qn[2]);
	double B3func = 0.5*(Tt_excess + fabs(Tt_excess));
	double B3 = B3func*signof(Sig[0])/fGf_II;
		
	double Shear_Q = Sig[0]*Sig[0] + (qn[1] - qn[0]*qn[3])*(qn[1] - qn[0]*qn[3]);
	double DQDN = qn[3];
	double DQDT = Sig[0]/sqrt(Shear_Q);

	double signfun = (1. + signof(Tt_excess))/2.;
	double zeta_Q = (qn[1] - qn[0]*qn[3])/sqrt(Shear_Q*Shear_Q*Shear_Q);
	//double DB3_DTanphi = -fabs(TNA)*signof(Sig[0])/fGf_II;
	double DB3_DTanphi = -fabs(TNA)*signof(Sig[0])*signfun/fGf_II;
	double DQDTDChi = Sig[0]*qn[3]*zeta_Q;
	double DQDTDC = -Sig[0]*zeta_Q;
	double DQDTDTanpsi = Sig[0]*qn[0]*zeta_Q;
	double DQDNDTanpsi = 1.;

	//dqbardq(0,0) = -falpha_chi*(B1*DQDN + B3*DQDT) + A2*B3*DQDTDChi;
	dqbardq(0,0) = DA1_DChi*B1*DQDN + DA2_DChi*B3*DQDT + A2*B3*DQDTDChi;
	dqbardq(0,1) = A2*B3*DQDTDC;
	dqbardq(0,2) = A2*DB3_DTanphi*DQDT;
	dqbardq(0,3) = A1*B1*DQDNDTanpsi + A2*B3*DQDTDTanpsi;
	dqbardq(1,0) = A4*B3*DQDTDChi;
	//dqbardq(1,1) = -falpha_c*(B1*DQDN + B3*DQDT) + A4*B3*DQDTDC;
	dqbardq(1,1) = DA3_DC*B1*DQDN + DA4_DC*B3*DQDT + A4*B3*DQDTDC;
	dqbardq(1,2) = A4*DB3_DTanphi*DQDT;
	dqbardq(1,3) = A3*B1*DQDNDTanpsi + A4*B3*DQDTDTanpsi;
	dqbardq(2,0) = A6*B3*DQDTDChi;
	dqbardq(2,1) = A6*B3*DQDTDC;
	//dqbardq(2,2) = -falpha_phi*B3*DQDT + A6*DB3_DTanphi*DQDT;
	dqbardq(2,2) = DA6_DTanphi*B3*DQDT + A6*DB3_DTanphi*DQDT;
	dqbardq(2,3) = A6*B3*DQDTDTanpsi;
	dqbardq(3,0) = A8*B3*DQDTDChi;
	dqbardq(3,1) = A8*B3*DQDTDC;
	dqbardq(3,2) = A8*DB3_DTanphi*DQDT;
	//dqbardq(3,3) = -falpha_psi*B3*DQDT + A8*B3*DQDTDTanpsi;
	dqbardq(3,3) = DA8_DTanpsi*B3*DQDT + A8*B3*DQDTDTanpsi;

	return dqbardq;
}

/* elastoplastic consistent tangent operator*/
const dMatrixT& MR2DT::Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma)
{
#pragma unused(sigma)
#pragma unused(jump_u)
#if __option(extended_errorcheck)
	if (jump_u.Length() != knumDOF) throw ExceptionT::kSizeMismatch;
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif

	int i, j;

	dMatrixT AA(6,6), I_mat(4,4), CMAT(6,6), AA_inv(6,6),
		A_qq(4,4), A_uu(2,2), A_uq(2,4), A_qu(4,2), ZMAT(2,4),
		ZMATP(4,2), dQdSig2(2,2), dqdbar(4,4), dqbardSig(4,2), dqbardq(4,4),
		dQdSigdq(2,4), KP(2,2), KP2(2,2), KEP(2,2), DMAT(6,6), EMAT(6,2), FMAT(6,2);

	dMatrixT I_m(2,2), Rmat(2,2), R_Inv(2,2), KE(2,2), KE_Inv(2,2),
		Ch(4,4), Ch_Inv(4,4), KE1(4,2), KE2(2,2), KE3(2,4), KEE(2,2),
		KEE_Inv(2,2);

	dArrayT u(2), up(2), du(2), dup(2), qn(4), qo(4), Rvec(6),Cvec(6),
		R(6), Rmod(6), Sig(2), Sig_I(2), dQdSig(2), dfdq(4), qbar(4),
		R2(6), X(6), V_sig(2), V_q(4), dfdSig(2), K1(2), K2(2);

	double bott, dlam, e;

	//I_m(0,0) = 1.;
	//I_m(0,1) = 0.;
	//I_m(1,0) = 0.;
	//I_m(1,1) = 1.;
	I_mat = 0.;

	if (state[12] == 0.){
		fStiffness[0] = fE_t;
		fStiffness[1] = 0.;
		fStiffness[2] = 0.;
		fStiffness[3] = fE_n;
	}
	else
		if (state[12] == 1.){
			for (i = 0; i<=3; ++i){
				qn[i] = state[i+6];
				I_mat(i,i) = 1.;
			}

			Sig[0] = state[0];
			Sig[1] = state[1];

			KEE = 0.;
			KEE(0,0) = fE_t;
			KEE(1,1) = fE_n;
			KEE_Inv.Inverse(KEE);

			dlam = state[11];
			dQdSig2_f(Sig, qn, dQdSig2);
			dQdSigdq_f(Sig, qn, dQdSigdq);
			dqbardSig_f(Sig, qn, dqbardSig);
			dqbardq_f(Sig, qn, dqbardq);

			for (i = 0; i<=5; ++i){
				for (j = 0; j<=5; ++j) {
					if (i<=1 & j<=1){
						AA_inv(i,j)  = KEE_Inv(i,j);
						AA_inv(i,j) += dlam*dQdSig2(i,j);
					}
					if (i<=1 & j>1){
						AA_inv(i,j)  = dQdSigdq(i,j-2);
						AA_inv(i,j) *= dlam;
					}
					if(i>1 & j<=1){
						AA_inv(i,j)  = dqbardSig(i-2,j);
						AA_inv(i,j) *= dlam;
					}
					if(i>1 & j >1){
						AA_inv(i,j)  = I_mat(i-2,j-2);
						AA_inv(i,j) *= -1.;
						AA_inv(i,j) += dlam*dqbardq(i-2,j-2);
					}
				}
			}

			AA.Inverse(AA_inv);
			
			dfdSig_f(Sig, qn, dfdSig);
			dfdq_f(Sig,qn, dfdq);
			dQdSig_f(Sig, qn, dQdSig);
			qbar_f(Sig, qn, qbar);

			for (i = 0; i<=5; ++i){
				if (i<=1){
					Rvec[i] = dfdSig[i];
					Cvec[i] = dQdSig[i];
				}
				if (i>1){
					Rvec[i] = dfdq[i-2];
					Cvec[i] = qbar[i-2];
				}
			}

			dArrayT tmpVec(6), Vvec(2), dVec(2);
			AA.Multx(Cvec,tmpVec);
			e = dArrayT::Dot(Rvec,tmpVec);

			for (i = 0; i<=5; ++i){
				for (j = 0; j<=5; ++j){
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

			FMAT.MultAB(DMAT, EMAT);
			for (i = 0; i <= 1; ++i){
				for (j = 0; j<=1; ++j){
					KEP(i,j) = FMAT(i,j);
				}
			}

	/*
	// Check for the consistent tangent
#if __option(extended_errorcheck)
	mr_ep_2d_out << setw(outputFileWidth) << "**********Check for consistent tangent**********" << endl;
	mr_ep_2d_out << setw(outputFileWidth) << "Tt = Sig[0] = state[0] = " << state[0] << endl;
	mr_ep_2d_out << setw(outputFileWidth) << "Tn = Sig[1] = state[1] = " << state[1] << endl;
	mr_ep_2d_out << setw(outputFileWidth) << "jump_u_t = " << state[2] << endl;
	mr_ep_2d_out << setw(outputFileWidth) << "jump_u_n = " << state[3] << endl;
	mr_ep_2d_out << setw(outputFileWidth) << "up_t = " << state[4] << endl;
	mr_ep_2d_out << setw(outputFileWidth) << "up_n = " << state[5] << endl;
	mr_ep_2d_out << setw(outputFileWidth) << "chi = " << state[6] << endl;
	mr_ep_2d_out << setw(outputFileWidth) << "cohesion = " << state[7] << endl;
	mr_ep_2d_out << setw(outputFileWidth) << "tan(phi) = " << state[8] << endl;
	mr_ep_2d_out << setw(outputFileWidth) << "tan(psi) = " << state[9] << endl;
	mr_ep_2d_out << setw(outputFileWidth) << "F = " << state[10] << endl;
	mr_ep_2d_out << setw(outputFileWidth) << "dlam = state[11] = " << state[11] << endl;
	mr_ep_2d_out << setw(outputFileWidth) << "double(iplastic) = " << state[12] << endl;
	mr_ep_2d_out << setw(outputFileWidth) << "normr = " << state[13] << endl;
	mr_ep_2d_out << setw(outputFileWidth) << "KEE(0,0) = fE_t = " << fE_t << endl;
	mr_ep_2d_out << setw(outputFileWidth) << "KEE(1,1) = fE_n = " << fE_n << endl;
	mr_ep_2d_out << setw(outputFileWidth) << "KEP(0,0) = " << KEP(0,0) << endl;
	mr_ep_2d_out << setw(outputFileWidth) << "KEP(0,1) = " << KEP(0,1) << endl;
	mr_ep_2d_out << setw(outputFileWidth) << "KEP(1,0) = " << KEP(1,0) << endl;
	mr_ep_2d_out << setw(outputFileWidth) << "KEP(1,1) = " << KEP(1,1) << endl;
	mr_ep_2d_out << setw(outputFileWidth) << "**********End of check**********" << endl;
#endif
	*/

		fStiffness[0] = KEP(0,0);
		fStiffness[2] = KEP(0,1);
		fStiffness[1] = KEP(1,0);
		fStiffness[3] = KEP(1,1);
	}

	return fStiffness;
}


/* print parameters to the output stream */
void MR2DT::Print(ostream& out) const
{
#ifndef _SIERRA_TEST_
	out << " Elastic tangential stiffness. . . . . . . . . . = " << fE_t << '\n';
	out << " Elastic Normal stiffness . . .. . . . . . . . . = " << fE_n     << '\n';
	out << " Mode_I Fracture Energy            . . . . . . . = " << fGf_I     << '\n';
	out << " Mode_II Fracture Energy            . . .  . . . = " << fGf_II << '\n';
	out << " Peak Cohesion                 . . . . . . . . . = " << fc_p    << '\n';
	out << " Residual Cohesion             . . . . . . . . . = " << fc_r    << '\n';
	out << " Peak Tensile Strength   . . . . . . . . . . . . = " << fchi_p << '\n';
	out << " Residual Tensile Strength . . . . . . . . . . . = " << fchi_r << '\n';
	out << " Peak Friction Angle         . . . . . . . . . . = " << fphi_p   << '\n';
	out << " Critical State Friction Angle       . . . . . . = " << fphi_r  << '\n';
	out << " Peak Dilation Angle.. . . . . . . . . . . . . . = " << fpsi_p   << '\n';
	out << " Coefficient of Tensile Strength Degradation ..  = " << falpha_chi   << '\n';
	out << " Coefficient of Cohesion Degradation. .. . . . . = " << falpha_c   << '\n';
	out << " Coefficient for Frictional Angle Degradation .  = " << falpha_phi   << '\n';
	out << " Coefficient for Dilation Angle Degradation  . . = " << falpha_psi   << '\n';
	out << " Error Tolerance for Yield Function. . . . . . . = " << fTol_1 << '\n';
	out << " Error Tolerance for Residual    . . . . . . . . = " << fTol_2 << '\n';
#endif
}

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default */
int MR2DT::NumOutputVariables(void) const { return 11; }

void MR2DT::OutputLabels(ArrayT<StringT>& labels) const
{
	labels.Dimension(11);
	labels[0] = "u_t";
	labels[1] = "u_n";
	labels[2] = "up_t";
	labels[3] = "up_n";
	labels[4] = "Chi";
	labels[5] = "Cohesion";
	labels[6] = "Friction Angle";
	labels[7] = "Dilation Angle";
	labels[8] = "Yield Function Value";
	labels[9] = "Norm of residuals";
	labels[10] = "No. of Iterations";
}

void MR2DT::ComputeOutput(const dArrayT& jump_u, const ArrayT<double>& state,
	dArrayT& output)
{
#pragma unused(jump_u)
#if __option(extended_errorcheck)
	if (state.Length() != NumStateVariables()) throw ExceptionT::kGeneralFail;
#endif
	output[0] = state[2];
	output[1] = state[3];
	output[2] = state[4];
	output[3] = state[5];
	output[4] = state[6];
	output[5] = state[7];
	// output[6] = state[8];
	output[6] = atan(state[8]); // state[8] = ftan_phi
	output[7] = atan(state[9]);
	output[8] = state[10];
	output[9] = state[13];
	output[10] = state[16];
	
}


bool MR2DT::NeedsNodalInfo(void) { return false; }

int MR2DT::NodalQuantityNeeded(void) 
{
	return 2;
}

void MR2DT::SetElementGroupsNeeded(iArrayT& iGroups) 
{	
	iGroups[0] = 1;
}

double MR2DT::signof(double r)
{
	if (fabs(r) < kSmall)
		return 0.;
	else
		return fabs(r)/r;
}

