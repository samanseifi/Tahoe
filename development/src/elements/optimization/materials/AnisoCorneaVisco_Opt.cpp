 /* $Id: AnisoCorneaVisco_Opt.cpp,v 1.5 2011/12/01 20:38:06 beichuan Exp $ */
#include "AnisoCorneaVisco_Opt.h"
#include "ParameterContainerT.h"
#include "ExceptionT.h"

#include <cmath>
#include "bessel.h"
#include "toolboxConstants.h"
#include "C1FunctionT.h"
#include "ParameterContainerT.h"
#include "ModelManagerT.h"
#include "ofstreamT.h"

#ifdef VIB_MATERIAL

/* point generator */
#include "EvenSpacePtsT.h"

#include "FungType.h"
#include "VWType.h"
#include "LanirFiber.h"

/*viscosity functions*/
#include "ScaledCsch.h"

const double Pi = acos(-1.0);
static const int perm[3][3] = {0,1,2,1,2,0,2,0,1};

using namespace Tahoe;

AnisoCorneaVisco_Opt::AnisoCorneaVisco_Opt(void):
	ParameterInterfaceT("aniso_cornea_visco_opto")
{
}

void AnisoCorneaVisco_Opt::DefineParameters(ParameterListT& list) const
{
	/*inherited*/
	AnisoCorneaVisco::DefineParameters(list);
	list.SetDescription("parameter order: [gamma, mu, beta, alpha_eq, alpha_neq_k, eta0_k, tau0_k, n, (phi)]");
}

/* information about subordinate parameter lists */
void AnisoCorneaVisco_Opt::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FSSolidMatT::DefineSubs(sub_list);

	/*material parameters for matrix*/
	sub_list.AddSub("matrix_opto_params", ParameterListT::Once);

	/* choice of energy potential for fibrils */
	sub_list.AddSub("eq_fibril_opto_params", ParameterListT::Once);
	
	sub_list.AddSub("neq_fibril_opto_params", ParameterListT::Any);

	/* choice of energy potential for fibrils */
	sub_list.AddSub("viscosity_opto_params", ParameterListT::Any);

	/* choice of fibril distribution funcion */
	sub_list.AddSub("fibril_distr_opto_params", ParameterListT::Once);	
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* AnisoCorneaVisco_Opt::NewSub(const StringT& name) const
{
	/* inherited */
	ParameterInterfaceT* sub = FSSolidMatT::NewSub(name);
	if (sub) 
	{
		return sub;
	}
	else if (name == "matrix_opto_params")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);

		/* bound */
		LimitT lower(0.0, LimitT::Lower);
		
		/* exponential functions*/
		ParameterContainerT matrix("Neo-Hookean");
		ParameterT mu(ParameterT::Double, "shear_modulus");
		mu.AddLimit(lower);
		ParameterT gamma(ParameterT::Double, "bulk_modulus");
		gamma.AddLimit(lower);
		matrix.AddParameter(mu);
		matrix.AddParameter(gamma);
		matrix.SetDescription("param order, [kappa, mu]");
		choice->AddSub(matrix);
		return choice;
	}
	else if (name == "eq_fibril_opto_params" )
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);

		ParameterContainerT fung("eq_fung_type");		
		{
			LimitT lower(0.0, LimitT::Lower);

			ParameterT alpha(ParameterT::Double, "alpha");
			ParameterT beta(ParameterT::Double, "beta");

			fung.AddParameter(alpha);
			fung.AddParameter(beta);
			alpha.AddLimit(lower);
			beta.AddLimit(lower);

			/* set the description */
			fung.SetDescription("param order, [alpha_eq, beta]");	
			choice->AddSub(fung);
		}

		ParameterContainerT fung0("eq_vw_type");		
		{
			LimitT lower(0.0, LimitT::Lower);
			
			ParameterT alpha(ParameterT::Double, "alpha");
			ParameterT beta(ParameterT::Double, "beta");
			
			fung0.AddParameter(alpha);
			fung0.AddParameter(beta);
			alpha.AddLimit(lower);
			beta.AddLimit(lower);
			
			/* set the description */
			fung.SetDescription("param order, [alpha_eq, beta]");	
			choice->AddSub(fung0);
		}
		return(choice);
	}
	else if (name == "neq_fibril_opto_params" )
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);

		ParameterContainerT fung("neq_fung_type");		
		{
			LimitT lower(0.0, LimitT::Lower);

			ParameterT alpha(ParameterT::Double, "alpha");

			fung.AddParameter(alpha);
			alpha.AddLimit(lower);

			/* set the description */
			fung.SetDescription("param order, [alpha_neq]");	
			choice->AddSub(fung);
		}
		ParameterContainerT fung0("neq_vw_type");		
		{
			LimitT lower(0.0, LimitT::Lower);
			
			ParameterT alpha(ParameterT::Double, "alpha");
			
			fung0.AddParameter(alpha);
			alpha.AddLimit(lower);
			
			/* set the description */
			fung0.SetDescription("param order, [alpha_neq]");	
			choice->AddSub(fung0);
		}
		
		return(choice);
	}
	else if (name == "viscosity_opto_params")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);

		ParameterContainerT sinh("scaled_csch");		
		{
			LimitT lower(0.0, LimitT::Lower);

			ParameterT eta0(ParameterT::Double, "eta0");
			ParameterT tau0(ParameterT::Double, "tau0");

			sinh.AddParameter(eta0);
			sinh.AddParameter(tau0);
			eta0.AddLimit(lower);
//			tau0.SetDefault(100000);
//			tau0.AddLimit(100000, LimitT::Only);
			tau0.AddLimit(0.0, LimitT::Lower);

			/* set the description */
			sinh.SetDescription("param oder, [eta0, tau0]");	
			choice->AddSub(sinh);
		}
		return(choice);
	}
	else if (name == "fibril_distr_opto_params")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);
		LimitT lower(0.0, LimitT::LowerInclusive);
		LimitT upper(1.0, LimitT::UpperInclusive);

		/* sine & cosine raised to a power */
		ParameterContainerT circ("circumferential");		
		{
			ParameterT m(ParameterT::Integer, "number_of_peaks");
			ParameterT k(ParameterT::Double, "concentration_k");
			ParameterT phi(ParameterT::Double, "location_phi");

			m.AddLimit(2, LimitT::LowerInclusive);
			k.AddLimit(-1.0e-13, LimitT::Lower);
			phi.AddLimit(-3.14159, LimitT::LowerInclusive);
			phi.AddLimit(3.14159, LimitT::UpperInclusive);
			
			circ.AddParameter(m);
			circ.AddParameter(k);
			circ.AddParameter(phi);
		}
		/* set the description */	
		circ.SetDescription("D(theta) = 1/(2 Pi I0(k)) exp(k*cos*(m*(theta-phi))) ");	
		choice->AddSub(circ);

		/* sine & cosine raised to a power */
		ParameterContainerT orthogonal("meridional");		
		{
			ParameterT m(ParameterT::Integer, "number_of_peaks");
			ParameterT k(ParameterT::Double, "concentration_k");
			ParameterT phi(ParameterT::Double, "location_phi");

			m.AddLimit(2, LimitT::LowerInclusive);
			k.AddLimit(-1.0e-13, LimitT::Lower);
			phi.AddLimit(-3.14159, LimitT::LowerInclusive);
			phi.AddLimit(3.14159, LimitT::UpperInclusive);
			
			orthogonal.AddParameter(m);
			orthogonal.AddParameter(k);
			orthogonal.AddParameter(phi);
		}
		orthogonal.SetDescription("D(theta) = 1/(2 Pi I0(k)) exp(k*cos*(m*(theta-phi))) ");	
		choice->AddSub(orthogonal);
		
		return(choice);
	}
}

/* accept parameter list */
void AnisoCorneaVisco_Opt::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "AnisoCorneaVisco_Opt::TakeParameterList";

	/* inherited */
	FSFiberMatT::TakeParameterList(list);

	/*initializing some parameters for cornea_mod fibril distribution*/
	int num_neq_pot = list.NumLists("neq_fibril_opto_params");
	int num_visc = list.NumLists("viscosity_opto_params");

	if (num_visc != num_neq_pot)
		ExceptionT::GeneralFail("AnisoCorneaVisco::TakeParameterList", 
			"number of viscosity functions does not match number of nonequilibrium potentials");

	int b = list.GetParameter("pressure_model");
	fVolType = (b == kBlatz) ? kBlatz : kOgden;

	fNumFibProcess = num_neq_pot;
	fPotential.Dimension(fNumFibProcess+1);
	if(fNumFibProcess > 0)
		fViscosity.Dimension(fNumFibProcess);
		
	/*count parameters*/
	fnum_params = 2; //matrix params
	fnum_params += 2; // fiber eq params: alpha_eq, beta
		
	for (int i = 0; i < fNumFibProcess; i++)
	{
		fnum_params ++; //fiber neq params: alpha_neq;
		fnum_params += 2; //viscosity func: eta0, tau0;
	}

	fnum_params += 2; //distribution function: k, phi
	
	fParamGrads.Dimension(dSymMatrixT::NumValues(NumSD()), fnum_params);
	fconstraint_grad.Dimension(fnum_params);
	fconstraint_grad = 0.0;
	flabels.Dimension(fnum_params);
	fparams.Dimension(fnum_params);

	flabels[0] = "gamma";
	flabels[1] = "mu";
	flabels[3] = "alpha_eq";
	flabels[2] = "beta";
	int dex = 4;
	for (int i = 0; i < fNumFibProcess; i++)
	{
		StringT alpha, eta, tau;
		alpha = "alpha_neq_";
		eta = "eta0_";
		tau = "tau0_";
		flabels[dex++] = alpha.Append(i,3);
		flabels[dex++] = eta.Append(i,3);
		flabels[dex++] = tau.Append(i,3);
	}

	flabels[dex++] = "k";
	flabels[dex++] = "phi";

	const ParameterListT& matrix = list.GetListChoice(*this, "matrix_opto_params");
	if (matrix.Name() == "Neo-Hookean")
	{
		fMu = matrix.GetParameter("shear_modulus");
		fGamma = matrix.GetParameter("bulk_modulus");
	}
	fparams[0] = fGamma;
	fparams[1] = fMu;

	const ParameterListT& potential = list.GetListChoice(*this, "eq_fibril_opto_params");
	if (potential.Name() == "eq_fung_type")
	{
		double alpha_eq = potential.GetParameter("alpha");
		double beta = potential.GetParameter("beta");
		fPotential[0] = new FungType(alpha_eq, beta);
		if (!fPotential[0]) throw ExceptionT::kOutOfMemory;
		
		fparams[2] = beta;
		fparams[3] = alpha_eq;
		ffiber_type = 1;
	}
	else if (potential.Name() == "eq_vw_type")
	{
		double alpha_eq = potential.GetParameter("alpha");
		double beta = potential.GetParameter("beta");
		fPotential[0] = new VWType(alpha_eq, beta);
		if (!fPotential[0]) throw ExceptionT::kOutOfMemory;
		
		fparams[2] = beta;
		fparams[3] = alpha_eq;
		ffiber_type = 0;
	}
	else 
		ExceptionT::GeneralFail(caller, "no such potential");
	dex = 4;
	
	for (int i = 0; i < fNumFibProcess; i++)
	{
		const ParameterListT& neq_potential = list.GetListChoice(*this, "neq_fibril_opto_params", i);
		if (neq_potential.Name() == "neq_fung_type")
		{
			double alpha_neq = neq_potential.GetParameter("alpha");

			fPotential[i+1] = new FungType(alpha_neq,fparams[2]);
			if (!fPotential[i+1]) throw ExceptionT::kOutOfMemory;
			
			fparams[dex++] = alpha_neq;
		}		
		else if (neq_potential.Name() == "neq_vw_type")
		{
			double alpha_neq = neq_potential.GetParameter("alpha");
			
			fPotential[i+1] = new VWType(alpha_neq,fparams[2]);
			if (!fPotential[i+1]) throw ExceptionT::kOutOfMemory;
			
			fparams[dex++] = alpha_neq;
		}		
		const ParameterListT& visc = list.GetListChoice(*this, "viscosity_opto_params", i);
		if (visc.Name() == "scaled_csch")
		{
			double a = visc.GetParameter("eta0");
			double b = visc.GetParameter("tau0");
			fViscosity[i] = new ScaledCsch(a,b);
			if (!fViscosity[i]) throw ExceptionT::kOutOfMemory;
		
			fparams[dex++] = a;
			fparams[dex++] = b;
		}
	}
	
	const ParameterListT& distr = list.GetListChoice(*this, "fibril_distr_opto_params");
	if (distr.Name() == "circumferential")
	{
		fm = distr.GetParameter("number_of_peaks");
		fk = distr.GetParameter("concentration_k");
		fphi = distr.GetParameter("location_phi");
		 
		fparams[dex++] = fk;
		fparams[dex++] = fphi;

		fDType = kCircumferential;
	}
	else if (distr.Name() == "meridional")
	{
		fm = distr.GetParameter("number_of_peaks");
		fk = distr.GetParameter("concentration_k");
		fphi = distr.GetParameter("location_phi");
		 
		fparams[dex++] = fk;
		fparams[dex++] = fphi;

		fDType = kMeridional;
	}

		
	/* allocate memory */
	/*2D fiber stress and modulus*/
	fNumFibStress = dSymMatrixT::NumValues(fNumSD-1);
	fNumFibModuli = dSymMatrixT::NumValues(fNumFibStress);
	fFiberStretch.Dimension(fNumSD-1);
	fFiberStress.Dimension(fNumSD-1);
	fFiberMod.Dimension(fNumFibStress);

	/*viscous fiber stretch at time step n, viscous stretch at time step n and vn*/
	fFiberStretch_v.Dimension(fNumSD-1);
	fFiberStretch_vn.Dimension(fNumSD-1);

	/*Dimension work spaces*/
	fInverse.Dimension(fNumSD);

	fCalg.Dimension(fNumFibStress);	
	/* allocate memory */
	/*dimension invserse viscosity matrix*/
	fiVisc.Dimension(fNumFibStress);

	/*Dimension work spaces*/
	fFlowStress.Dimension(fNumSD-1);
	fResidual.Dimension(fNumSD-1);
	fiK.Dimension(fNumFibStress);
	fG.Dimension(fNumFibStress);
	fMod1.Dimension(fNumFibStress);
	fMod2.Dimension(fNumFibStress);
	fMod3.Dimension(fNumFibStress);
	fVec.Dimension(fNumFibStress);

	/*dimension state variable storage arrays*/
	fNumMatProcess = 0;
	int numprocess = fNumFibProcess;
	SetStateVariables(numprocess);

	/* point generator */
	int points = list.GetParameter("n_points");
	fCircle = new EvenSpacePtsT(points);
	Construct();
	int numbonds = fI4.Length();
	fI4e.Dimension(numbonds);
	
	/*workspace*/
	fb.Dimension(NumSD());
	fdsf_dalpha.Dimension(fNumSD-1);
	fdsf_dbeta.Dimension(fNumSD-1);
	fdsf_deta0.Dimension(fNumSD-1);
	fdsf_dtau0.Dimension(fNumSD-1);
	fdsf_dk.Dimension(fNumSD-1);
	fdsf_dphi.Dimension(fNumSD-1);

	fdS.Dimension(fNumSD);
	fdS = 0.0;
}



const dArray2DT& AnisoCorneaVisco_Opt::ds_ij_dlambda_q(void)
{
	fParamGrads = 0.0;
	
	const dMatrixT& F = F_mechanical();
		
	fb.MultAAT(F);
	fC.MultATA(F);
	double J = F.Det();

	/*matrix params*/
	if (fVolType == kOgden)
	{
		double third = 1.0/3.0;
		double I3 = fC.Det();
		double I3rthird = pow(I3,-third);
		double I1 = fC[0]+fC[1]+fC[2];
	
		// dsig_dgamma
		double coeff = (1.0/J)*0.5*fGamma*(I3-1.0);
		fParamGrads(0,0) = coeff;
		fParamGrads(1,0) = coeff;
		fParamGrads(2,0) = coeff;
		fParamGrads(3,0) = 0.0;
		fParamGrads(4,0) = 0.0;
		fParamGrads(5,0) = 0.0;	
	
		// dsig_dmu
		fParamGrads(0,1) = (1.0/J)*I3rthird *(fb[0]  - third*I1);
		fParamGrads(1,1) = (1.0/J)*I3rthird *(fb[1]  - third*I1);
		fParamGrads(2,1) = (1.0/J)*I3rthird *(fb[2]  - third*I1);
		fParamGrads(3,1) = (1.0/J)*I3rthird * fb[3];
		fParamGrads(4,1) = (1.0/J)*I3rthird * fb[4];
		fParamGrads(5,1) = (1.0/J)*I3rthird * fb[5];
	}
	else if(fVolType == kBlatz)
	{
		double coeff = pow(J, -2.0*fGamma);
		double coeff2 = coeff* 2.0*log(J);

		//	dsig_dgamma
		fParamGrads(0,0) = fMu/J *(coeff2);
		fParamGrads(1,0) = fMu/J *(coeff2);
		fParamGrads(2,0) = fMu/J *(coeff2);
		fParamGrads(3,0) = 0.0;
		fParamGrads(4,0) = 0.0;
		fParamGrads(5,0) = 0.0;	
	
		//	dsig_dmu
		fParamGrads(0,1) = 1.0/J *(fb[0]  - coeff);
		fParamGrads(1,1) = 1.0/J *(fb[1]  - coeff);
		fParamGrads(2,1) = 1.0/J *(fb[2]  - coeff);
		fParamGrads(3,1) = 1.0/J*fb[3];
		fParamGrads(4,1) = 1.0/J*fb[4];
		fParamGrads(5,1) = 1.0/J*fb[5];
	}
	/////////////////////////////////////////////////////////
	/*eq fiber*/
	
	/*****initialize********/
	/* initialize kernel pointers */
	double elnum = fFSFiberMatSupport->CurrElementNumber();
	const dArrayT& pj = fjacobians[elnum];
//	const dArrayT& pj_drho = fjac_drho[elnum];
	const dArrayT& pj_dk = fjac_dk[elnum];
	const dArrayT& pj_dphi = fjac_dphi[elnum];
	
	double *p1  = fStressTable(0);
	double *p2  = fStressTable(1);
	double *p3  = fStressTable(2);
	
	fdsf_dalpha = 0.0;
	double& ss1 = fdsf_dalpha[0]; /*dsf11_dalpha*/
	double& ss2 = fdsf_dalpha[1]; /*dsf22_dalpha*/
	double& ss3 = fdsf_dalpha[2]; /*dsf12_dalpha*/

	fdsf_dbeta = 0.0;
	double& t1 = fdsf_dbeta[0];
	double& t2 = fdsf_dbeta[1];
	double& t3 = fdsf_dbeta[2];
	
	fdsf_dk = 0.0;
	double& k1 = fdsf_dk[0];
	double& k2 = fdsf_dk[1];
	double& k3 = fdsf_dk[2];

	fdsf_dphi = 0.0;
	double& z1 = fdsf_dphi[0];
	double& z2 = fdsf_dphi[1];
	double& z3 = fdsf_dphi[2];
	
	/*integrate w.r.t in-plane orientation theta*/
	ComputeFiberStretch(fC, fFiberStretch);
	CompI4(fFiberStretch);
	fPotential[0]->MapDFunction(fI4, fdU);

	double alpha_eq = fparams[3];
	double beta = fparams[2];

	for (int i = 0; i < fI4.Length(); i++)
	{
		/*dseq_dalpha*/
		double delta = fI4[i] - 1.0;
		
		double factors;
		if(ffiber_type == 0)
//OLD			factors=2.0*pj[i]*(exp(beta*delta) -1.0/(fI4[i]*fI4[i]));
			factors=2.0*pj[i]*(exp(beta*delta)-1.0);
		else if (ffiber_type ==1)
			factors=2.0*delta*(exp(beta*delta*delta))*pj[i];

		ss1 += factors * p1[i];
		ss2 += factors * p2[i];
		ss3 += factors * p3[i];
	

		/*dseq_dbeta*/
		double factort;
		if(ffiber_type == 0)
			factort = 2.0*pj[i]*(alpha_eq*delta*exp(beta*delta));
		else if (ffiber_type ==1)
			factort = 2.0*(alpha_eq*delta*delta*delta*exp(beta*delta*delta))*pj[i];

		t1 += factort * p1[i];
		t2 += factort * p2[i];
		t3 += factort * p3[i];
		
		double seq = 2.0*fdU[i];
		/*dseq_dk*/
		double factork = seq*pj_dk[i];
		k1 += factork * p1[i];
		k2 += factork * p2[i];
		k3 += factork * p3[i];
	
		/*dseq_dk*/
		double factorz = seq*pj_dphi[i];
		z1 += factorz * p1[i];
		z2 += factorz * p2[i];
		z3 += factorz * p3[i];

	}
	/* rotate and assemble dS_dalpha to lab coordinates */
	AssembleFiberStress(fdsf_dalpha, fdS, dSymMatrixT::kOverwrite);

	/*push forward to current configuration*/
	fdS.SetToScaled(1.0/J, PushForward(F, fdS));
	
	fParamGrads(0,3) = fdS[0];
	fParamGrads(1,3) = fdS[1];
	fParamGrads(2,3) = fdS[2];
	fParamGrads(3,3) = fdS[3];
	fParamGrads(4,3) = fdS[4];
	fParamGrads(5,3) = fdS[5];
	
	////////////////// NOT YET IMPLEMENTED CORRECTLY  ////////////////////
	/*neq fiber*/
	/*get time step*/
	const double dt = fFSFiberMatSupport->TimeStep();
	/*retrieve history*/
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());

	/*initialize*/	
	double& v1 = fdsf_deta0[0];
	double& v2 = fdsf_deta0[1];
	double& v3 = fdsf_deta0[2];

	double& u1 = fdsf_dtau0[0];
	double& u2 = fdsf_dtau0[1];
	double& u3 = fdsf_dtau0[2];

	int dex = 4;
	for (int j = 0;  j< fNumFibProcess & fNumFibProcess > 0; j++)
	{
		fdsf_dalpha = 0.0;
		fdsf_deta0 = 0.0;
		fdsf_dtau0 = 0.0;

		ComputeFiberStretch(fC_v[j], fFiberStretch_v);
		CompI4(fFiberStretch, fFiberStretch_v);
		fI4e = fI4;
		fI4e /= fI4v; /*fI4 = I4e*/

		/* derivatives of the potential */
		fPotential[j+1]->MapDFunction(fI4e, fdU);
		fPotential[j+1]->MapDDFunction(fI4e, fddU);

		/* initialize kernel pointers */
		double alpha_neq = fparams[dex];
		double eta0 = fparams[dex+1];
		double tau0 = fparams[dex+2];
		/*integrate w.r.t in-plane orientation theta*/
		for (int i = 0; i < fI4e.Length(); i++)
		{

			/* c: stiffness of neq stress response*/
			double cneq_b = (4.0*fddU[i]*fI4e[i] + 4.0*fdU[i])*fI4e[i];
			double tauneq = (2.0*fdU[i]*fI4e[i]);
			double sneq = tauneq/fI4[i];
			double eta = fViscosity[j]->Function(tauneq);
			double coeff =  0.00;
			coeff = cosh(tauneq/tau0);
			
			/* k: stiffness of residual for lambda_v_k*/
			//assuming constant viscosity eta
			double k = 1.0 + dt/eta* (cneq_b *eta/eta0 *coeff - tauneq);
			
			/*dsneq_dalphaneq*/			
			double factors = sneq/alpha_neq * pj[i];
			factors *= (1.0 - dt* cneq_b/(eta0*k) * coeff);
			ss1 += factors * p1[i];
			ss2 += factors * p2[i];
			ss3 += factors * p3[i];
			
			
			/*dsneq_dbeta*/
			double r = fI4e[i] - 1.0;
			double factort =2.0*pj[i]/fI4v[i] * (alpha_neq*r*exp(beta*r));
			factort *= (1.0 - dt* cneq_b/(eta0*k) * coeff);
			t1 += factort * p1[i];
			t2 += factort * p2[i];
			t3 += factort * p3[i];
			
			/*dsneq_deta0*/
			double factorv = sneq/eta * pj[i];
			factorv *= dt*cneq_b/(eta0*k);
			v1 += factorv* p1[i];
			v2 += factorv* p2[i];
			v3 += factorv* p3[i];

			/*dsneq_tau0*/
			double factoru = -sneq/eta * pj[i];
			factoru *= dt*cneq_b/(tau0*k)*(1.0-eta/eta0*coeff);
			u1 += factoru* p1[i];
			u2 += factoru* p2[i];
			u3 += factoru* p3[i];

			/*dseq_dk*/
			double factork = sneq * pj_dk[i];
			k1 += factork * p1[i];
			k2 += factork * p2[i];
			k3 += factork * p3[i];
	
			/*dseq_dk*/
			double factorz = sneq * pj_dphi[i];
			z1 += factorz * p1[i];
			z2 += factorz * p2[i];
			z3 += factorz * p3[i];
		}

		AssembleFiberStress(fdsf_dalpha, fdS, dSymMatrixT::kOverwrite);
		/*push forward to current configuration*/
		fdS.SetToScaled(1.0/F.Det(), PushForward(F, fdS));
	
		fParamGrads(0,dex) = fdS[0];
		fParamGrads(1,dex) = fdS[1];
		fParamGrads(2,dex) = fdS[2];
		fParamGrads(3,dex) = fdS[3];
		fParamGrads(4,dex) = fdS[4];
		fParamGrads(5,dex) = fdS[5];
		dex++;

		
		/* rotate and assemble dS_deta to lab coordinates */
		AssembleFiberStress(fdsf_deta0, fdS, dSymMatrixT::kOverwrite);

		/*push forward to current configuration*/
		fdS.SetToScaled(1.0/F.Det(), PushForward(F, fdS));
	
		fParamGrads(0,dex) = fdS[0];
		fParamGrads(1,dex) = fdS[1];
		fParamGrads(2,dex) = fdS[2];
		fParamGrads(3,dex) = fdS[3];
		fParamGrads(4,dex) = fdS[4];
		fParamGrads(5,dex) = fdS[5];
		dex++;
		
		/* rotate and assemble dS_deta to lab coordinates */
		AssembleFiberStress(fdsf_dtau0, fdS, dSymMatrixT::kOverwrite);

		/*push forward to current configuration*/
		fdS.SetToScaled(1.0/F.Det(), PushForward(F, fdS));
	
		fParamGrads(0,dex) = fdS[0];
		fParamGrads(1,dex) = fdS[1];
		fParamGrads(2,dex) = fdS[2];
		fParamGrads(3,dex) = fdS[3];
		fParamGrads(4,dex) = fdS[4];
		fParamGrads(5,dex) = fdS[5];
		dex++;
	} //for j < fNumFibProcess;

	///////////////////////////////////////////////////

	/* rotate and assemble stress_gradient to lab coordinates */
	AssembleFiberStress(fdsf_dbeta, fdS, dSymMatrixT::kOverwrite);

	/*push forward to current configuration*/
	fdS.SetToScaled(1.0/J, PushForward(F, fdS));
	
	fParamGrads(0,2) = fdS[0];
	fParamGrads(1,2) = fdS[1];
	fParamGrads(2,2) = fdS[2];
	fParamGrads(3,2) = fdS[3];
	fParamGrads(4,2) = fdS[4];
	fParamGrads(5,2) = fdS[5];
	
	/* rotate and assemble stress_gradient to lab coordinates */
/*	AssembleFiberStress(fdsf_drho, fdS, dSymMatrixT::kOverwrite);

	fdS.SetToScaled(1.0/J, PushForward(F, fdS));
	
	fParamGrads(0,dex) = fdS[0];
	fParamGrads(1,dex) = fdS[1];
	fParamGrads(2,dex) = fdS[2];
	fParamGrads(3,dex) = fdS[3];
	fParamGrads(4,dex) = fdS[4];
	fParamGrads(5,dex) = fdS[5];
	dex++;
*/	
	/* rotate and assemble stress_gradient to lab coordinates */
	AssembleFiberStress(fdsf_dk, fdS, dSymMatrixT::kOverwrite);

	/*push forward to current configuration*/
	fdS.SetToScaled(1.0/J, PushForward(F, fdS));
	
	fParamGrads(0,dex) = fdS[0];
	fParamGrads(1,dex) = fdS[1];
	fParamGrads(2,dex) = fdS[2];
	fParamGrads(3,dex) = fdS[3];
	fParamGrads(4,dex) = fdS[4];
	fParamGrads(5,dex) = fdS[5];
	dex++;

	/* rotate and assemble stress_gradient to lab coordinates */
	AssembleFiberStress(fdsf_dphi, fdS, dSymMatrixT::kOverwrite);

	/*push forward to current configuration*/
	fdS.SetToScaled(1.0/F.Det(), PushForward(F, fdS));
	
	fParamGrads(0,dex) = fdS[0];
	fParamGrads(1,dex) = fdS[1];
	fParamGrads(2,dex) = fdS[2];
	fParamGrads(3,dex) = fdS[3];
	fParamGrads(4,dex) = fdS[4];
	fParamGrads(5,dex) = fdS[5];
	
	return(fParamGrads);
}

/* Initialize angle tables */
void AnisoCorneaVisco_Opt::Construct(void)
{
	/* fetch points */
	const dArray2DT& points = fCircle->CirclePoints(0.0);
	int numbonds = points.MajorDim();
		
	int nel = NumElements();     
	fjacobians.Dimension(nel);
	fjac_dk.Dimension(nel);
	fjac_dphi.Dimension(nel);

	/* length table */
	fI4.Dimension(numbonds);
	fI4v.Dimension(numbonds);

	/* potential tables */
	fU.Dimension(numbonds);
	fdU.Dimension(numbonds);
	fddU.Dimension(numbonds);

	/* STRESS angle tables - by associated stress component */
	fStressTable.Dimension(fNumFibStress, numbonds);
	  	
	/* MODULI angle tables - using Cauchy symmetry */ 	
	fModuliTable.Dimension(fNumFibModuli, numbonds);	

	/* set pointers */
	double *s1 = fStressTable(0);
	double *s2 = fStressTable(1);
	double *s3 = fStressTable(2);

	double *c11 = fModuliTable(0);
	double *c22 = fModuliTable(1);
	double *c33 = fModuliTable(2);
	double *c23 = fModuliTable(3);
	double *c13 = fModuliTable(4);
	double *c12 = fModuliTable(5);

	StringT name("distributions.dat");
	ofstreamT dist_out(name);

	/*modified von-mises distribution function [0,Pi/2]. integral[D(theta),{theta,-Pi,Pi}] = 2 Pi */
	if (fDType == kMeridional)
	{			 		
		/*integrate dtemp - min to normalize distribution function*/
		dArrayT& jac = fjacobians[0];
		dArrayT& jac_dk = fjac_dk[0];
		dArrayT& jac_dphi = fjac_dphi[0];
		
		jac.Dimension(numbonds);
		jac_dk.Dimension(numbonds);
		jac_dphi.Dimension(numbonds);
		
		const dArrayT& angles = fCircle->CircleAngles(0.0);
		const dArray2DT& points = fCircle->CirclePoints(0.0);  /*set jacobians*/
		const dArrayT& weights = fCircle->Jacobians();
		
		double I0 = besseli0(fk);
		double I1 = besseli1(fk);

		for (int i = 0; i < numbonds; i++)
		{
			/*function*/
			double theta = angles[i];
			double y = fk*cos(fm*(theta-fphi));

			double f = exp(y)/(I0*2.0*Pi);
			jac[i] = f * weights[i];
			
			/*dfunction_dk*/
			double df_dk = f*(I0*cos(fm*(theta-fphi))-I1)/(I0);
			jac_dk[i] = df_dk*weights[i];
			
			/*dfunction_dhi*/
			double df_dphi = fm*fk*f*sin(fm*(theta-fphi));
			jac_dphi[i] = df_dphi*weights[i];
		}
		
		for (int el = 1; el < nel; el++)
		{
			fjacobians[el].Dimension(numbonds);
			fjac_dk[el].Dimension(numbonds);
			fjac_dphi[el].Dimension(numbonds);
			
			fjacobians[el] = jac;
			fjac_dk[el] = jac_dk;
			fjac_dphi[el] = jac_dphi;
		}
	}
	/*modified von-mises distribution function [0,Pi]. integral[D(theta),{theta,-Pi,Pi}] = 2 Pi */
	else if(fDType == kCircumferential)
	{
		dArrayT xc(3);
		const dArray2DT& coordinates = fFSFiberMatSupport->InitialCoordinates();

		fFSFiberMatSupport->TopElement();
		int ielm = fFSFiberMatSupport->CurrElementNumber();
	
		const dArrayT& angles = fCircle->CircleAngles(0.0);
		const dArray2DT& points = fCircle->CirclePoints(0.0);  /*set jacobians*/
		const dArrayT& weights = fCircle->Jacobians();
		
		double I0 = besseli0(fk);
		double I1 = besseli1(fk);
		
		while (fFSFiberMatSupport->NextElement()) 
		{
			int ielm = fFSFiberMatSupport->CurrElementNumber();
		
			/* calculate element centroid */
			iArrayT nodes = fFSFiberMatSupport->CurrentElement()->NodesX();
			int nen = NumElementNodes();       
			xc = 0.0;
			for (int i = 0; i < nen; i++) 
			{
				for (int j = 0; j < coordinates.MinorDim(); j++)
				xc[j] += coordinates(nodes[i],j);
			}
			xc /= nen;
		
			int inormal = 2; // HACK
			double x1 = xc[perm[inormal][1]];
			double x2 = xc[perm[inormal][2]];
		
			// projected polar coordinates
			double r = sqrt(x1*x1 + x2*x2);
			double xi = atan2(x2,x1); //angle measured from radius axis
			xi -= 0.5*Pi; //rotate to hoop axis
				
			dArrayT& jac = fjacobians[ielm];
			dArrayT& jac_dk = fjac_dk[ielm];
			dArrayT& jac_dphi = fjac_dphi[ielm];
			
			jac.Dimension(numbonds);
			jac_dk.Dimension(numbonds);
			jac_dphi.Dimension(numbonds);
			
			const dArrayT& angles = fCircle->CircleAngles(0.0);
			const dArray2DT& points = fCircle->CirclePoints(0.0);  /*set jacobians*/
			const dArrayT& weights = fCircle->Jacobians();
		
			double I0 = besseli0(fk);
			double I1 = besseli1(fk);

			for (int i = 0; i < numbonds; i++)
			{
				/*function*/
				double theta = angles[i];
				double y = fk*cos(fm*(theta-xi-fphi));

				double f = exp(y)/(I0*2.0*Pi);
				jac[i] = f * weights[i];
			
				/*dfunction_dk*/
				double df_dk = f*(I0*cos(fm*(theta-fphi-xi))-I1)/(I0);
				jac_dk[i] = df_dk*weights[i];
			
				/*dfunction_dphi*/
				double df_dphi = fm*fk*f*sin(fm*(theta-fphi-xi));
				jac_dphi[i] = df_dphi*weights[i];
			}
		}
	}	

	for (int i = 0; i < numbonds; i++)
	{
		/* direction cosines */
		const double *xsi = points(i);
		double cosi = xsi[0];
		double sini = xsi[1];
		
		/* stress angle tables */
		s1[i] = cosi*cosi;      
		s2[i] = sini*sini;
		s3[i] = sini*cosi;
	
		/* moduli angle tables */
		c11[i] = s1[i]*s1[i];
		c22[i] = s2[i]*s2[i];
		c33[i] = s3[i]*s3[i];
		c23[i] = s2[i]*s3[i];
		c13[i] = s1[i]*s3[i];
		c12[i] = s2[i]*s1[i];
	}
}


	/*modified von-mises distribution function [0,Pi/2]. integral[D(theta),{theta,-Pi,Pi}] = 2 Pi */
/*	if (fDType == kOrthogonal)
	{			 		
		dArrayT& jac = fjacobians[0];
		dArrayT& jac_drho = fjac_drho[0];
		dArrayT& jac_dn = fjac_dn[0];
		dArrayT& jac_dxi = fjac_dxi[0];
		
		jac.Dimension(numbonds);
		jac_drho.Dimension(numbonds);			
		jac_dn.Dimension(numbonds);
		jac_dxi.Dimension(numbonds);
		
		const dArrayT& angles = fCircle->CircleAngles(0.0);
		const dArray2DT& points = fCircle->CirclePoints(0.0); 
		const dArrayT& weights = fCircle->Jacobians();
		
		double I0 = besseli0(n1);
		double I1 = besseli1(n1);
		double denom = exp(n1)*I0 - 1.0;

		for (int i = 0; i < numbonds; i++)
		{
			double theta = angles[i];
			double y = 2.0*n1*cos(2.0*theta)*cos(2.0*theta);
			double f = (exp(y) -1.0)/denom;
			jac[i] = (c1 + (1.0-c1)*f) * weights[i];
			
			jac_drho[i] =(1.0 - f) *weights[i];
			
			double df_dn = 2.0*exp(y)*cos(2.0*theta)*cos(2.0*theta);
			df_dn += (1.0 - exp(y))*exp(n1)*(I0 + I1)/denom;
			df_dn /= denom;
			jac_dn[i] = (1.0-c1)*df_dn*weights[i];
			
			jac_dxi[i] = 0.0;
		}
		
		for (int el = 1; el < nel; el++)
		{
			fjacobians[el].Dimension(numbonds);
			fjac_drho[el].Dimension(numbonds);			
			fjac_dn[el].Dimension(numbonds);
			fjac_dxi[el].Dimension(numbonds);
			
			fjacobians[el] = jac;
			fjac_drho[el] = jac_drho;
			fjac_dn[el] = jac_dn;
			fjac_dxi[el] = jac_dxi;
		}
	}
	else if(fDType == kCircumferential)
	{
		dArrayT xc(3);
		const dArray2DT& coordinates = fFSFiberMatSupport->InitialCoordinates();

		fFSFiberMatSupport->TopElement();
		int ielm = fFSFiberMatSupport->CurrElementNumber();
	
		const dArrayT& angles = fCircle->CircleAngles(0.0);
		const dArray2DT& points = fCircle->CirclePoints(0.0); 
		const dArrayT& weights = fCircle->Jacobians();
		
		double I0 = besseli0(n1);
		double I1 = besseli1(n1);
		double denom = exp(n1)*I0 - 1.0;

		while (fFSFiberMatSupport->NextElement()) 
		{
			int ielm = fFSFiberMatSupport->CurrElementNumber();
		
			iArrayT nodes = fFSFiberMatSupport->CurrentElement()->NodesX();
			int nen = NumElementNodes();       
			xc = 0.0;
			for (int i = 0; i < nen; i++) 
			{
				for (int j = 0; j < coordinates.MinorDim(); j++)
				xc[j] += coordinates(nodes[i],j);
			}
			xc /= nen;
		
			int inormal = 2; // HACK
			double x1 = xc[perm[inormal][1]];
			double x2 = xc[perm[inormal][2]];
		
			// projected polar coordinates
			double r = sqrt(x1*x1 + x2*x2);
			double phi = atan2(x2,x1);
			
			//shift not yet implemented//
			phi -= 0.5*Pi;
//			phi -= xi;
				
			dArrayT& jac = fjacobians[ielm];
			dArrayT& jac_drho = fjac_drho[ielm];
			dArrayT& jac_dn = fjac_dn[ielm];
			dArrayT& jac_dxi = fjac_dxi[ielm];
			
			jac.Dimension(numbonds);
			jac_drho.Dimension(numbonds);			
			jac_dn.Dimension(numbonds);
			jac_dxi.Dimension(numbonds);
			
			for (int i = 0; i < numbonds; i++)
			{
				double theta = angles[i];
				double y = 2.0*n1*cos(theta-phi)*cos(theta-phi);
				double f = (exp(y) -1.0)/denom;
				jac[i] = (c1 + (1.0-c1)*f) * weights[i];
			
				jac_drho[i] =(1.0- f) *weights[i];
			
				double df_dn = 2.0*exp(y)*cos(theta-phi)*cos(theta-phi);
				df_dn += (1.0 - exp(y))*exp(n1)*(I0 + I1)/denom;
				df_dn /= denom;
				jac_dn[i] = (1.0-c1)*df_dn*weights[i];
			
//				double df_dxi = -4.0*exp(y)*n1*cos(theta-phi)*sin(theta-phi);
//				df_dxi /= denom;
//				jac_dxi[i] = (1.0-c1)*df_dxi*weights[i];
				jac_dxi[i] = 0.0;

			}
		}
	}
	else if (fDType == kOneFiber)
	{
		dArrayT& jac = fjacobians[0];
		dArrayT& jac_drho = fjac_drho[0];
		dArrayT& jac_dn = fjac_dn[0];
		dArrayT& jac_dxi = fjac_dxi[0];
		
		jac.Dimension(numbonds);
		jac_drho.Dimension(numbonds);			
		jac_dn.Dimension(numbonds);
		jac_dxi.Dimension(numbonds);
		
		const dArrayT& angles = fCircle->CircleAngles(0.0);
		const dArray2DT& points = fCircle->CirclePoints(0.0);  
		const dArrayT& weights = fCircle->Jacobians();
		
		double I0 = besseli0(n1);
		double I1 = besseli1(n1);
		double denom = exp(n1)*I0 - 1.0;

		for (int i = 0; i < numbonds; i++)
		{
			double theta = angles[i];
			double y = 2.0*n1*cos(theta-xi)*cos(theta-xi);
			double f = (exp(y) -1.0)/denom;
			jac[i] = (c1 + (1.0-c1)*f) * weights[i];
			
			jac_drho[i] =(1.0- f) *weights[i];
			
			double df_dn = 2.0*exp(y)*cos(theta-xi)*cos(theta-xi);
			df_dn += (1.0 - exp(y))*exp(n1)*(I0 + I1)/denom;
			df_dn /= denom;
			jac_dn[i] = (1.0-c1)*df_dn*weights[i];
			
			double df_dxi = 4.0*exp(y)*n1*cos(theta-xi)*sin(theta-xi);
			df_dxi /= denom;
			jac_dxi[i] = (1.0-c1)*df_dxi*weights[i];
		}
		
		for (int el = 1; el < nel; el++)
		{
			fjacobians[el].Dimension(numbonds);
			fjac_drho[el].Dimension(numbonds);			
			fjac_dn[el].Dimension(numbonds);
			fjac_dxi[el].Dimension(numbonds);
			
			fjacobians[el] = jac;
			fjac_drho[el] = jac_drho;
			fjac_dn[el] = jac_dn;
			fjac_dxi[el] = jac_dxi;
		}		
	}
*/

#endif /*VIB_MATERIAL*/
