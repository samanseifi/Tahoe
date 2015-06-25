/* $Id: TersoffPairT.cpp,v 1.4 2011/12/01 21:11:39 bcyansfn Exp $ */
#include "TersoffPairT.h"
#include <iostream>
#include <cmath>
#include "dMatrixT.h"
#include "dArrayT.h"

using namespace Tahoe;

// initialize static parameters
double TersoffPairT::s_A = 0.0;
double TersoffPairT::s_B = 0.0;
double TersoffPairT::s_lambda = 0.0;
double TersoffPairT::s_mu = 0.0;
double TersoffPairT::s_beta = 0.0;
double TersoffPairT::s_n = 0.0;
double TersoffPairT::s_c = 0.0;
double TersoffPairT::s_d = 0.0;
double TersoffPairT::s_h = 0.0;
double TersoffPairT::s_chi = 0.0;
double TersoffPairT::s_R = 0.0;
double TersoffPairT::s_S = 0.0;

const double Pi = acos(-1.0);
				 
				 
// constructor
// Uses Tersoff potential as defined in Tersoff PRB 1989
TersoffPairT::TersoffPairT(void):
	f_A(0.0),
	f_B(0.0),
	f_lambda(0.0),
	f_mu(0.0),
	f_beta(0.0),
	f_n(0.0),
	f_c(0.0),
	f_d(0.0),
	f_h(0.0),
	f_chi(0.0),
	f_R(0.0),
	f_S(0.0)
{
	SetName("Tersoff");
}

// return a pointer to the energy function
TersoffPropertyT::EnergyFunction TersoffPairT::getEnergyFunction(void)
{
	// copy my data to static
	s_A = f_A;
	s_lambda = f_lambda;
	s_R = f_R;
    s_S = f_S;
    s_B = f_B;
	s_mu = f_mu;
	s_beta = f_beta;
	s_n = f_n;
	s_c = f_c;
	s_d = f_d;
	s_h = f_h;
	s_chi = f_chi;

    
	// return function pointer
	return TersoffPairT::Energy;
}

// return a pointer to the force function
TersoffPropertyT::ForceFunction TersoffPairT::getForceFunction_ij(void)
{
	// copy my data to static
	s_A = f_A;
	s_lambda = f_lambda;
	s_R = f_R;
    s_S = f_S;
    s_B = f_B;
	s_mu = f_mu;
	s_beta = f_beta;
	s_n = f_n;
	s_c = f_c;
	s_d = f_d;
	s_h = f_h;
	s_chi = f_chi;

	// return function pointer
	return TersoffPairT::Force_ij;
}

// return a pointer to the force function
TersoffPropertyT::ForceFunction TersoffPairT::getForceFunction_ik(void)
{
	// copy my data to static
	s_A = f_A;
	s_lambda = f_lambda;
	s_R = f_R;
    s_S = f_S;
    s_B = f_B;
	s_mu = f_mu;
	s_beta = f_beta;
	s_n = f_n;
	s_c = f_c;
	s_d = f_d;
	s_h = f_h;
	s_chi = f_chi;

	// return function pointer
	return TersoffPairT::Force_ik;
}

// return a pointer to the force function
TersoffPropertyT::ForceFunction TersoffPairT::getForceFunction_jk(void)
{
	// copy my data to static
	s_A = f_A;
	s_lambda = f_lambda;
	s_R = f_R;
    s_S = f_S;
    s_B = f_B;
	s_mu = f_mu;
	s_beta = f_beta;
	s_n = f_n;
	s_c = f_c;
	s_d = f_d;
	s_h = f_h;
	s_chi = f_chi;

	// return function pointer
	return TersoffPairT::Force_jk;
}


// return a pointer to the stiffness function
TersoffPropertyT::StiffnessFunction TersoffPairT::getStiffnessFunction(void)
{
	// copy my data to static
	s_A = f_A;
	s_lambda = f_lambda;
	s_R = f_R;
    s_S = f_S;
    s_B = f_B;
	s_mu = f_mu;
	s_beta = f_beta;
	s_n = f_n;
	s_c = f_c;
	s_d = f_d;
	s_h = f_h;
	s_chi = f_chi;

    
	// return function pointer
	return TersoffPairT::Stiffness;
}

// describe the parameters needed by the interface
void TersoffPairT::DefineParameters(ParameterListT& list) const
{
	// inherited
	ParameterInterfaceT::DefineParameters(list);

	ParameterT mass(fMass, "mass");
	mass.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(mass, ParameterListT::ZeroOrOnce);

	ParameterT A(f_A, "rep_energy_scale_Aij");
	A.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(A, ParameterListT::ZeroOrOnce);
	
	ParameterT B(f_B, "attr_energy_scale_Bij");
	B.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(B, ParameterListT::ZeroOrOnce);
	
	ParameterT lambda(f_lambda, "rep_energy_exponent_lambdaij");
	lambda.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(lambda, ParameterListT::ZeroOrOnce);
	
	ParameterT mu(f_mu, "attr_energy_exponent_muij");
	mu.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(mu, ParameterListT::ZeroOrOnce);
	
	ParameterT beta(f_beta, "bond_order_coeff1_betai");
	beta.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(beta, ParameterListT::ZeroOrOnce);
	
	ParameterT n(f_n, "bond_order_exponent_ni");
	n.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(n, ParameterListT::ZeroOrOnce);
	
	ParameterT c(f_c, "bond_order_coeff2_ci");
	c.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(c, ParameterListT::ZeroOrOnce);
	
	ParameterT d(f_d, "bond_order_coeff3_di");
	d.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(d, ParameterListT::ZeroOrOnce);
	
	ParameterT h(f_h, "bond_order_coeff4_hi");
	list.AddParameter(h, ParameterListT::ZeroOrOnce);
	
	ParameterT chi(f_chi, "bond_order_scaling_chi_ij");
	chi.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(chi, ParameterListT::ZeroOrOnce);
	
	ParameterT R(f_R, "cutoff_func_length_1_Rij");
	R.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(R, ParameterListT::ZeroOrOnce);
	
	ParameterT S(f_S, "cutoff_func_length_2_Sij");
	S.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(S, ParameterListT::ZeroOrOnce);

}

// accept parameter list
void TersoffPairT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	/* All parameters default to Si unless specified
	 * Si Parameters taken from Devanathan et al, Journal of Nuclear Materials, 1998
	 */
	const ParameterT* opt_param = list.Parameter("mass");
	if (opt_param)
		fMass = *opt_param;
	else
		fMass = 28.0855; // mass in amu (from webelements)
	
	opt_param = list.Parameter("rep_energy_scale_Aij");
	if (opt_param)
		f_A = *opt_param;
	else
		f_A = 1830.8; // scaling term in eV
	
	opt_param = list.Parameter("attr_energy_scale_Bij");
	if (opt_param)
		f_B = *opt_param;
	else
		f_B = 471.18; // scaling term in eV
		
	opt_param = list.Parameter("rep_energy_exponent_lambdaij");
	if (opt_param)
		f_lambda = *opt_param;
	else
		f_lambda = 2.4799; // exponent in angstrom^-1
		
	opt_param = list.Parameter("attr_energy_exponent_muij");
	if (opt_param)
		f_mu = *opt_param;
	else
		f_mu = 1.7322; // exponent in angstrom^-1
		
	opt_param = list.Parameter("bond_order_coeff1_betai");
	if (opt_param)
		f_beta = *opt_param;
	else
		f_beta = 1.1e-6; // bond order coeff 1 - pure number
	
	opt_param = list.Parameter("bond_order_exponent_ni");
	if (opt_param)
		f_n = *opt_param;
	else
		f_n = .78734; // bond order exponent - pure number
	
	opt_param = list.Parameter("bond_order_coeff2_ci");
	if (opt_param)
		f_c = *opt_param;
	else
		f_c = 100390.0; // bond order coeff 2 - pure number
	
	opt_param = list.Parameter("bond_order_coeff3_di");
	if (opt_param)
		f_d = *opt_param;
	else
		f_d = 16.217; // bond order coeff 3 - pure number
	
	opt_param = list.Parameter("bond_order_coeff4_hi");
	if (opt_param)
		f_h = *opt_param;
	else
		f_h = -0.59825; // bond order coeff 4 - pure number
		
	opt_param = list.Parameter("bond_order_scaling_chi_ij");
	if (opt_param)
		f_chi = *opt_param;
	else
		f_chi = 1.0; // bond order scaling coeff (for identical atoms = 1)- pure number
		
	opt_param = list.Parameter("cutoff_func_length_1_Rij");
	if (opt_param)
		f_R = *opt_param;
	else
		f_R = 2.7; // Cutoff length parameter 1 - angstrom
	
	opt_param = list.Parameter("cutoff_func_length_2_Sij");
	if (opt_param)
		f_S = *opt_param;
	else
		f_S = 3.0; // Cutoff length parameter 2 - angstrom	


	SetRange(1.1*f_S);					// range, with a tolerance
	SetNearestNeighbor(1.1*f_R);	// nearest neighbor, with an tolerance

}

/***********************************************************************
 * Private
 ***********************************************************************/

double TersoffPairT::Energy(double rij, iArrayT neighbors, const int j, const AutoArrayT<int> type, ArrayT<TersoffPropertyT*> tersoff_properties, nMatrixT<int>& properties_map, const dArray2DT& coords)
{
	// Determine the energy
	
	// Determine value of cutoff function
	double FCij = 0.0;
	
	// figure out some needed values
	if (rij < s_R)
		FCij = 1.0;	
	else if (rij >= s_R && rij < s_S)
		FCij = .5 + .5*cos(Pi * (rij - s_R)/(s_S - s_R));
	else
		FCij = 0.0;
	
	// Determine the repulsive part
	double FR = s_A * exp(-s_lambda * rij);
	
	// Determine the attractive part
	double FA = -s_B * exp(-s_mu * rij);
	
	// Determine the bond order parameter
	double ksi_ij = 0.0;
	for (int k = 1; k < neighbors.Length(); k++)
	{
		// skip if atom_k == atom_j or atom_i
		if (k == j || neighbors[k] == neighbors[0])
			continue;
		
		// determine some needed values
		const double* x_i = coords(neighbors[0]);
		const double* x_j = coords(neighbors[j]);
		const double* x_k = coords(neighbors[k]);
		
		double r_ij[3], r_ik[3];
		r_ij[0] = x_j[0] - x_i[0];
		r_ij[1] = x_j[1] - x_i[1];
		r_ij[2] = x_j[2] - x_i[2];
		
		r_ik[0] = x_k[0] - x_i[0];
		r_ik[1] = x_k[1] - x_i[1];
		r_ik[2] = x_k[2] - x_i[2];
		double rik = sqrt(r_ik[0]*r_ik[0] + r_ik[1]*r_ik[1] + r_ik[2]*r_ik[2]);
		
		double costheta = (r_ij[0]*r_ik[0] + r_ij[1]*r_ik[1] + r_ij[2]*r_ik[2])/(rij*rik);
		
		// now get the needed interaction info for i-k bonds (the cutoff distance parameters)
		double Rik = tersoff_properties[properties_map(type[neighbors[0]], type[neighbors[k]])]->GetR();
		double Sik = tersoff_properties[properties_map(type[neighbors[0]], type[neighbors[k]])]->GetS();
		
		// Determine the parts of the bond order parameter
		double c2 = s_c*s_c;
		double d2 = s_d*s_d;
		
		double g = 1 + c2/d2 - c2/(d2 + pow((s_h - costheta),2.0));
		
		// Determine value of cutoff function
		double FCik = 0.0;
		if (rik < Rik)
		{
			ksi_ij += 1.0 * g;
		}
		else if (rik >= Rik && rik < Sik)
		{
			FCik = .5 + .5*cos(Pi * (rik - Rik)/(Sik - Rik));
			ksi_ij += FCik * g;
		}
		else
			ksi_ij += 0.0;		
	}
	
	// Assemble bond order term
	double bij = s_chi * pow((1 + (pow(s_beta,s_n)*pow(ksi_ij,s_n))), -1/(2*s_n));
	if (neighbors.Length() == 3 || ksi_ij == 0.0)
		bij = 1.0;	// if no ikth bond, assume bond order = 1, results in morselike potential
	
	// return the energy
	double Vij = FCij*FR + bij*FCij*FA;
	return (.5*Vij);
}

double TersoffPairT::Force_ij(double rij, iArrayT neighbors, const int j, const int kk, const AutoArrayT<int> type, ArrayT<TersoffPropertyT*> tersoff_properties, nMatrixT<int>& properties_map, const dArray2DT& coords)
{
	// Determine the i,j force contribution
	
	// Determine value of cutoff function and its derivative w/respect to rij
	double FCij = 0.0;
	double dFCijdr = 0.0;
	
	// figure out some needed values
	
	if (rij < s_R)
	{
		FCij = 1.0;
		dFCijdr = 0.0;
	}	
	else if (rij >= s_R && rij < s_S)
	{
		FCij = .5 + .5*cos(Pi * (rij - s_R)/(s_S - s_R));
		dFCijdr = -Pi/(2*(s_S - s_R)) * sin(Pi * (rij - s_R)/(s_S - s_R));
	}
	else
	{
		FCij = 0.0;
		dFCijdr = 0.0;
	}
		
	// Determine the repulsive part and its derivative
	double FR = s_A * exp(-s_lambda * rij);
	double dFRdr = - s_lambda * FR;
	
	// Determine the attractive  part and its derivative
	double FA = -s_B * exp(-s_mu * rij);
	double dFAdr = -s_mu * FA;
	
	// Determine the bond order parameter and its derivative w/respect to rij
	// derivative information obtained from Miejie Tang's Thesis, MIT, 1995
	double ksi_ij = 0.0;
	double dksi_ijdrij = 0.0;
	for (int k = 1; k < neighbors.Length(); k++)
	{
		// skip if atom_k == atom_j or atom_i
		if (k == j || neighbors[k] == neighbors[0])
			continue;
		
		// determine some needed values
		const double* x_i = coords(neighbors[0]);
		const double* x_j = coords(neighbors[j]);
		const double* x_k = coords(neighbors[k]);
		
		double r_ij[3], r_ik[3];
		r_ij[0] = x_j[0] - x_i[0];
		r_ij[1] = x_j[1] - x_i[1];
		r_ij[2] = x_j[2] - x_i[2];
		
		r_ik[0] = x_k[0] - x_i[0];
		r_ik[1] = x_k[1] - x_i[1];
		r_ik[2] = x_k[2] - x_i[2];
		double rik = sqrt(r_ik[0]*r_ik[0] + r_ik[1]*r_ik[1] + r_ik[2]*r_ik[2]);
		
		double costheta = (r_ij[0]*r_ik[0] + r_ij[1]*r_ik[1] + r_ij[2]*r_ik[2])/(rij*rik);
		
		// now get the needed interaction info for i-k bonds (the cutoff distance parameters)
		double Rik = tersoff_properties[properties_map(type[neighbors[0]], type[neighbors[k]])]->GetR();
		double Sik = tersoff_properties[properties_map(type[neighbors[0]], type[neighbors[k]])]->GetS();
		
		// Determine the parts of the bond order parameter
		double c2 = s_c*s_c;
		double d2 = s_d*s_d;
		
		double g = 1 + c2/d2 - c2/(d2 + pow((s_h - costheta),2.0));
		double dgdrij = (-2.0*c2*(s_h - costheta))/pow((d2 + pow((s_h - costheta),2)),2.0) * ((1/rik) - (costheta/rij)); 
		
		// Determine value of cutoff function
		double FCik = 0.0;
		if (rik < Rik)
		{
			ksi_ij += 1.0 * g;
			dksi_ijdrij += 1.0 * dgdrij;
		}
		else if (rik >= Rik && rik < Sik)
		{
			FCik = .5 + .5*cos(Pi * (rik - Rik)/(Sik - Rik));
			ksi_ij += FCik * g;
			dksi_ijdrij += FCik * dgdrij;
		}
		else
		{
			ksi_ij += 0.0;
			dksi_ijdrij += 0.0;
		}		
	}
	
	// Assemble bond order term
	double bij = s_chi * pow((1 + (pow(s_beta,s_n)*pow(ksi_ij,s_n))), -1/(2*s_n));
	
	double dbijdZij = (-.5*s_chi* pow((1 + (pow(s_beta,s_n)*pow(ksi_ij,s_n))), -(1+(2*s_n))/(2*s_n))) * (pow(s_beta,s_n-1.0)*pow(ksi_ij,s_n-1.0));
	double dbijdrij = s_beta * dbijdZij * dksi_ijdrij;
	
	if (neighbors.Length() == 3 || ksi_ij == 0.0)
	{
		bij = 1.0;	// if no ikth bond, assume bond order function = 1, results in morselike potential, no third body effects
		dbijdrij = 0.0;
	}
		
	// return force
	double Fij = (dFCijdr*FR) + (FCij*dFRdr) + (dFCijdr*bij*FA) + (FCij*dbijdrij*FA) + (FCij*bij*dFAdr);
	return (.5*Fij);
}

double TersoffPairT::Force_ik(double rij, iArrayT neighbors, const int j, const int kk, const AutoArrayT<int> type, ArrayT<TersoffPropertyT*> tersoff_properties, nMatrixT<int>& properties_map, const dArray2DT& coords)
{
	// Determine the i,k force contribution
	
	// Determine value of cutoff function (derivative w/respect to rik == 0)
	double FCij = 0.0;
	
	// figure out some needed values
	
	if (rij < s_R)
	{
		FCij = 1.0;
	}	
	else if (rij >= s_R && rij < s_S)
	{
		FCij = .5 + .5*cos(Pi * (rij - s_R)/(s_S - s_R));
	}
	else
	{
		FCij = 0.0;
	}
	
	// Determine the attractive part (derivatives and other parts are zero)
	double FA = -s_B * exp(-s_mu * rij);
	
	// Determine the bond order parameter and its derivative w/respect to rij
	// derivative information obtained from Miejie Tang's Thesis, MIT, 1995
	double ksi_ij = 0.0;
	double dksi_ijdrik = 0.0;
	for (int k = 1; k < neighbors.Length(); k++)
	{
		// skip if atom_k == atom_j or atom_i
		if (k == j || neighbors[k] == neighbors[0])
			continue;
		
		// determine some needed values
		const double* x_i = coords(neighbors[0]);
		const double* x_j = coords(neighbors[j]);
		const double* x_k = coords(neighbors[k]);
		
		double r_ij[3], r_ik[3];
		r_ij[0] = x_j[0] - x_i[0];
		r_ij[1] = x_j[1] - x_i[1];
		r_ij[2] = x_j[2] - x_i[2];
		
		r_ik[0] = x_k[0] - x_i[0];
		r_ik[1] = x_k[1] - x_i[1];
		r_ik[2] = x_k[2] - x_i[2];
		double rik = sqrt(r_ik[0]*r_ik[0] + r_ik[1]*r_ik[1] + r_ik[2]*r_ik[2]);
		
		double costheta = (r_ij[0]*r_ik[0] + r_ij[1]*r_ik[1] + r_ij[2]*r_ik[2])/(rij*rik);
		
		// now get the needed interaction info for i-k bonds (the cutoff distance parameters)
		double Rik = tersoff_properties[properties_map(type[neighbors[0]], type[neighbors[k]])]->GetR();
		double Sik = tersoff_properties[properties_map(type[neighbors[0]], type[neighbors[k]])]->GetS();
		
		// Determine the parts of the bond order parameter
		double c2 = s_c*s_c;
		double d2 = s_d*s_d;
		
		double g = 1 + c2/d2 - c2/(d2 + pow((s_h - costheta),2.0));
		
		// Determine value of cutoff function
		double FCik = 0.0;
		if (rik < Rik)
			ksi_ij += 1.0 * g;
		else if (rik >= Rik && rik < Sik)
		{
			FCik = .5 + .5*cos(Pi * (rik - Rik)/(Sik - Rik));
			ksi_ij += FCik * g;
		}
		else
			ksi_ij += 0.0;
		
		// Figure out the i,k derivative of ksi
		if (k == kk)
		{
			double dgdrik = (-2.0*c2*(s_h - costheta))/pow((d2 + pow((s_h - costheta),2)),2.0) * ((1/rij) - (costheta/rik));
			
			if (rik < Rik)
			{
				dksi_ijdrik = 1.0 * dgdrik;
			}
			else if (rik >= Rik && rik < Sik)
			{
				FCik = .5 + .5*cos(Pi * (rik - Rik)/(Sik - Rik));
				dksi_ijdrik = FCik * dgdrik - (Pi/(2*(Sik-Rik)))*sin(Pi * (rik - Rik)/(Sik - Rik))*g;
			}
			else
			{
				dksi_ijdrik = 0.0;
			}
		}
		
	}
	
	// Assemble bond order term
	double bij = s_chi * pow((1 + (pow(s_beta,s_n)*pow(ksi_ij,s_n))), -1/(2*s_n));
	
	double dbijdZij = (-.5*s_chi* pow((1 + (pow(s_beta,s_n)*pow(ksi_ij,s_n))), -(1+(2*s_n))/(2*s_n))) * (pow(s_beta,s_n-1.0)*pow(ksi_ij,s_n-1.0));
	double dbijdrik = s_beta * dbijdZij * dksi_ijdrik;
	
	if (neighbors.Length() == 3 || ksi_ij == 0.0)
	{
		bij = 1.0;	// if no ikth bond, assume bond order function = 1, results in morselike potential, no third body effects
		dbijdrik = 0.0;
	}
	
	// return force
	double Fik = (FCij*dbijdrik*FA);
	return (.5*Fik);
}

double TersoffPairT::Force_jk(double rij, iArrayT neighbors, const int j, const int kk, const AutoArrayT<int> type, ArrayT<TersoffPropertyT*> tersoff_properties, nMatrixT<int>& properties_map, const dArray2DT& coords)
{
	// Determine the j,k force contribution
	
	// Determine value of cutoff function (derivative w/respect to rik == 0)
	double FCij = 0.0;
	
	// figure out some needed values
	
	if (rij < s_R)
	{
		FCij = 1.0;
	}	
	else if (rij >= s_R && rij < s_S)
	{
		FCij = .5 + .5*cos(Pi * (rij - s_R)/(s_S - s_R));
	}
	else
	{
		FCij = 0.0;
	}
	
	// Determine the attractive part (derivatives and other parts are zero)
	double FA = -s_B * exp(-s_mu * rij);
	
	// Determine the bond order parameter and its derivative w/respect to rij
	// derivative information obtained from Miejie Tang's Thesis, MIT, 1995
	double ksi_ij = 0.0;
	double dksi_ijdrjk = 0.0;
	for (int k = 1; k < neighbors.Length(); k++)
	{
		// skip if atom_k == atom_j or atom_i
		if (k == j || neighbors[k] == neighbors[0])
			continue;
		
		// determine some needed values
		const double* x_i = coords(neighbors[0]);
		const double* x_j = coords(neighbors[j]);
		const double* x_k = coords(neighbors[k]);
		
		double r_ij[3], r_ik[3], r_jk[3];
		r_ij[0] = x_j[0] - x_i[0];
		r_ij[1] = x_j[1] - x_i[1];
		r_ij[2] = x_j[2] - x_i[2];
		
		r_ik[0] = x_k[0] - x_i[0];
		r_ik[1] = x_k[1] - x_i[1];
		r_ik[2] = x_k[2] - x_i[2];
		
		r_jk[0] = x_k[0] - x_j[0];
		r_jk[1] = x_k[1] - x_j[1];
		r_jk[2] = x_k[2] - x_j[2];
		
		double rik = sqrt(r_ik[0]*r_ik[0] + r_ik[1]*r_ik[1] + r_ik[2]*r_ik[2]);
		double rjk = sqrt(r_jk[0]*r_jk[0] + r_jk[1]*r_jk[1] + r_jk[2]*r_jk[2]);
		
		double costheta = (r_ij[0]*r_ik[0] + r_ij[1]*r_ik[1] + r_ij[2]*r_ik[2])/(rij*rik);
		
		// now get the needed interaction info for i-k bonds (the cutoff distance parameters)
		double Rik = tersoff_properties[properties_map(type[neighbors[0]], type[neighbors[k]])]->GetR();
		double Sik = tersoff_properties[properties_map(type[neighbors[0]], type[neighbors[k]])]->GetS();
		
		// Determine the parts of the bond order parameter
		double c2 = s_c*s_c;
		double d2 = s_d*s_d;
		
		double g = 1 + c2/d2 - c2/(d2 + pow((s_h - costheta),2.0));
		
		// Determine value of cutoff function
		double FCik = 0.0;
		if (rik < Rik)
			ksi_ij += 1.0 * g;
		else if (rik >= Rik && rik < Sik)
		{
			FCik = .5 + .5*cos(Pi * (rik - Rik)/(Sik - Rik));
			ksi_ij += FCik * g;
		}
		else
			ksi_ij += 0.0;
		
		// Figure out the j,k derivative of ksi
		if (k == kk)
		{
			double dgdrjk = (-2.0*c2*(s_h - costheta))/pow((d2 + pow((s_h - costheta),2)),2.0) * (-rjk/(rij*rik));
			
			if (rik < Rik)
			{
				dksi_ijdrjk = 1.0 * dgdrjk;
			}
			else if (rik >= Rik && rik < Sik)
			{
				FCik = .5 + .5*cos(Pi * (rik - Rik)/(Sik - Rik));
				dksi_ijdrjk = FCik * dgdrjk;
			}
			else
			{
				dksi_ijdrjk = 0.0;
			}
		}
		
	}
	
	// Assemble bond order term
	double bij = s_chi * pow((1 + (pow(s_beta,s_n)*pow(ksi_ij,s_n))), -1/(2*s_n));
	
	double dbijdZij = (-.5*s_chi* pow((1 + (pow(s_beta,s_n)*pow(ksi_ij,s_n))), -(1+(2*s_n))/(2*s_n))) * (pow(s_beta,s_n-1.0)*pow(ksi_ij,s_n-1.0));
	double dbijdrjk = s_beta * dbijdZij * dksi_ijdrjk;
	
	if (neighbors.Length() == 3 || ksi_ij == 0.0)
	{
		bij = 1.0;	// if no ikth bond, assume bond order function = 1, results in morselike potential, no third body effects
		dbijdrjk = 0.0;
	}
	
	// return force
	double Fjk = (FCij*dbijdrjk*FA);
	return (.5*Fjk);
}

double TersoffPairT::Stiffness(double rij, iArrayT neighbors, const int j, const AutoArrayT<int> type, ArrayT<TersoffPropertyT*> tersoff_properties, nMatrixT<int>& properties_map, const dArray2DT& coords)
{
	// Determine the stiffness
	
	// Determine value of cutoff function and its derivatives w/respect to rij
	double FCij = 0.0;
	double dFCijdr = 0.0;
	double d2FCijdr2 = 0.0;
	
	if (rij < s_R)
	{
		FCij = 1.0;
		dFCijdr = 0.0;
		d2FCijdr2 = 0.0;
	}	
	else if (rij >= s_R && rij < s_S)
	{
		FCij = .5 + .5*cos(Pi * (rij - s_R)/(s_S - s_R));
		dFCijdr = -Pi/(2*(s_S - s_R)) * sin(Pi * (rij - s_R)/(s_S - s_R));
		d2FCijdr2 = -.5 * (Pi/((s_S - s_R))) * (Pi/((s_S - s_R))) * cos(Pi * (rij - s_R)/(s_S - s_R));
	}
	else
	{
		FCij = 0.0;
		dFCijdr = 0.0;
		d2FCijdr2 = 0.0;
	}
	
	// Determine the repulsive part and its derivatives
	double FR = s_A * exp(-s_lambda * rij);
	double dFRdr = - s_lambda * FR;
	double d2FRdr2 = s_lambda * s_lambda * FR;
	
	// Determine the attractive  part and its derivatives
	double FA = -s_B * exp(-s_mu * rij);
	double dFAdr = -s_mu * FA;
	double d2FAdr2 = s_mu * s_mu * FA;
	
	// Determine the bond order parameter and its derivative w/respect to rij
	// derivative information obtained from Miejie Tang's Thesis, MIT, 1995
	
	double ksi_ij = 0.0;
	double dksi_ijdrij = 0.0;
	double d2ksi_ijdrij2 = 0.0;
	for (int k = 1; k < neighbors.Length(); k++)
	{
		// skip if atom_k == atom_j or atom_i
		if (k == j || neighbors[k] == neighbors[0])
			continue;
		
		// determine some needed values
		const double* x_i = coords(neighbors[0]);
		const double* x_j = coords(neighbors[j]);
		const double* x_k = coords(neighbors[k]);
		
		double r_ij[3], r_ik[3];
		r_ij[0] = x_j[0] - x_i[0];
		r_ij[1] = x_j[1] - x_i[1];
		r_ij[2] = x_j[2] - x_i[2];
		
		r_ik[0] = x_k[0] - x_i[0];
		r_ik[1] = x_k[1] - x_i[1];
		r_ik[2] = x_k[2] - x_i[2];
		double rik = sqrt(r_ik[0]*r_ik[0] + r_ik[1]*r_ik[1] + r_ik[2]*r_ik[2]);
		
		double costheta = (r_ij[0]*r_ik[0] + r_ij[1]*r_ik[1] + r_ij[2]*r_ik[2])/(rij*rik);
		
		// now get the needed interaction info for i-k bonds (the cutoff distance parameters)
		double Rik = tersoff_properties[properties_map(type[neighbors[0]], type[neighbors[k]])]->GetR();
		double Sik = tersoff_properties[properties_map(type[neighbors[0]], type[neighbors[k]])]->GetS();
		
		// Determine the parts of the bond order parameter
		double c2 = s_c*s_c;
		double d2 = s_d*s_d;
		
		double g = 1 + c2/d2 - c2/(d2 + pow((s_h - costheta),2.0));
		double dgdrij = (-2.0*c2*(s_h - costheta))/pow((d2 + pow((s_h - costheta),2)),2.0) * ((1/rik) - (costheta/rij)); 
		double d2gdrij2_pt1 = (((2.0*c2)/pow((d2 + pow((s_h - costheta),2.0)),2.0)) - ((8.0*c2*pow(s_h-costheta,2))/pow((d2 + pow((s_h - costheta),2.0)),3.0))) * pow(((1/rik) - (costheta/rij)),2.0);
		double d2gdrij2_pt2 = (-2.0*c2*(s_h - costheta))/pow((d2 + pow((s_h - costheta),2)),2.0) * ((2*costheta/(rij*rij)) - 1/(rij*rik));
		double d2gdrij2 = d2gdrij2_pt1 + d2gdrij2_pt2;
		
		// Determine value of cutoff function
		double FCik = 0.0;
		if (rik < Rik)
		{
			ksi_ij += 1.0 * g;
			dksi_ijdrij += 1.0 * dgdrij;
			d2ksi_ijdrij2 += 1.0 * d2gdrij2;
		}
		else if (rik >= Rik && rik < Sik)
		{
			FCik = .5 + .5*cos(Pi * (rik - Rik)/(Sik - Rik));
			ksi_ij += FCik * g;
			dksi_ijdrij += FCik * dgdrij;
			d2ksi_ijdrij2 += FCik * d2gdrij2;
		}
		else
		{
			ksi_ij += 0.0;
			dksi_ijdrij += 0.0;
			d2ksi_ijdrij2 += 0.0;
		}
		
	}
	
	// Assemble bond order term
	double bij = s_chi * pow((1 + (pow(s_beta,s_n)*pow(ksi_ij,s_n))), -1/(2*s_n));
	
	double dbijdZij = (-.5*s_chi* pow((1 + (pow(s_beta,s_n)*pow(ksi_ij,s_n))), -(1+(2*s_n))/(2*s_n))) * (pow(s_beta,s_n-1.0)*pow(ksi_ij,s_n-1.0));
	double dbijdrij = s_beta * dbijdZij * dksi_ijdrij;
	 
	double d2bijdZij2_pt1a = s_chi * ((2*s_n + 1)/(4*s_n) *pow((1 + (pow(s_beta,s_n)*pow(ksi_ij,s_n))), -(1+(4*s_n))/(2*s_n))*(pow(s_beta,2*s_n-2)*pow(ksi_ij,2*s_n-2)));
	double d2bijdZij2_pt1b = s_chi * (.5*(s_n - 1)*pow((1 + (pow(s_beta,s_n)*pow(ksi_ij,s_n))), -(1+(2*s_n))/(2*s_n))*(pow(s_beta,s_n-2)*pow(ksi_ij,s_n-2)));
	double d2bijdrij2_pt1 = pow(s_beta,2) * (d2bijdZij2_pt1a+d2bijdZij2_pt1b) * pow(dksi_ijdrij,2);
	double d2bijdrij2_pt2 = s_beta * dbijdZij * d2ksi_ijdrij2;
	double d2bijdrij2 = d2bijdrij2_pt1 + d2bijdrij2_pt2;
	
	if (neighbors.Length() == 3 || ksi_ij == 0.0)
	{
		bij = 1.0;	// if no ikth bond, assume bond order function = 1, results in morselike potential, no third body effects
		dbijdrij = 0.0;
		d2bijdrij2 = 0.0;
	}
	
	// return stiffness
	double Kij = (d2FCijdr2*FR) + (dFCijdr*dFRdr) + (FCij*d2FRdr2) + (dFCijdr*dFRdr) + (d2FCijdr2*bij*FA) + (dFCijdr*dbijdrij*FA) + (dFCijdr*bij*dFAdr) + (FCij*d2bijdrij2*FA) + (FCij*dbijdrij*dFAdr) + (dFCijdr*dbijdrij*dFAdr) + (FCij*dbijdrij*dFAdr) + (FCij*bij*d2FAdr2);
	return (.5*Kij);
}
