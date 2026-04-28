/* $Id: PenaltyContact3DT.cpp,v 1.17 2011/12/01 21:11:36 bcyansfn Exp $ */
/* created: paklein (02/09/2000) */
#include "PenaltyContact3DT.h"

#include <cmath>
#include <iostream>
#include <iomanip>

#include "ifstreamT.h"
#include "eIntegratorT.h"

using namespace Tahoe;

/* vector functions */
inline static void CrossProduct(const double* A, const double* B, double* AxB)
{   AxB[0] = A[1]*B[2] - A[2]*B[1];
	AxB[1] = A[2]*B[0] - A[0]*B[2];
	AxB[2] = A[0]*B[1] - A[1]*B[0];
};

inline static double Dot(const double* A, const double* B)
{ return A[0]*B[0] + A[1]*B[1] + A[2]*B[2]; };

inline static void Vector(const double* start, const double* end, double* v)
{
	v[0] = end[0] - start[0];
	v[1] = end[1] - start[1];
	v[2] = end[2] - start[2];
};

/* constructor */
PenaltyContact3DT::PenaltyContact3DT(const ElementSupportT& support):
	Contact3DT(support),
	fK(0.0),
	fMu(0.0),
	fFrictionEps(1.0e-6)
{
	SetName("contact_3D_penalty");
}

/* describe the parameters needed by the interface */
void PenaltyContact3DT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	Contact3DT::DefineParameters(list);

	/* penalty stiffness */
	ParameterT stiffness(ParameterT::Double, "penalty_stiffness");
	stiffness.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(stiffness);

	/* Coulomb friction coefficient (0 = frictionless, default).
	 * When >0, regularised kinetic Coulomb is applied:
	 *   f_t = -mu * |f_n| * v_t / sqrt(|v_t|^2 + eps^2)  */
	ParameterT mu(ParameterT::Double, "friction_coefficient");
	mu.SetDefault(0.0);
	mu.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(mu, ParameterListT::ZeroOrOnce);

	ParameterT eps(ParameterT::Double, "friction_epsilon_velocity");
	eps.SetDefault(1.0e-6);
	eps.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(eps, ParameterListT::ZeroOrOnce);
}

/* accept parameter list */
void PenaltyContact3DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	Contact3DT::TakeParameterList(list);

	/* contact stiffness */
	fK = list.GetParameter("penalty_stiffness");

	/* Coulomb friction (defaults: mu=0 frictionless, eps=1e-6 regularisation) */
	const ParameterT* p_mu = list.Parameter("friction_coefficient");
	if (p_mu) fMu = *p_mu;
	const ParameterT* p_eps = list.Parameter("friction_epsilon_velocity");
	if (p_eps) fFrictionEps = *p_eps;
	if (fMu > 0.0)
		std::cout << "PenaltyContact3DT: Coulomb friction enabled, mu="
		          << fMu << ", eps=" << fFrictionEps << std::endl;

	/* dimension workspace */
	fElCoord.Dimension(fNumFacetNodes + 1, NumSD());
	fElRefCoord.Dimension(fNumFacetNodes + 1, NumSD());
	fElDisp.Dimension(fNumFacetNodes + 1, NumDOF());

	fdc_du.Dimension(NumSD(), fElDisp.Length());
	fdn_du.Dimension(NumSD(), fElDisp.Length());
	fM1.Dimension(NumSD());
	fM2.Dimension(NumSD(), fElDisp.Length());
	fV1.Dimension(fElDisp.Length());

	double third = 1.0/3.0;
	double* p = fdc_du.Pointer();
	*p++ =-third;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ =-third;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ =-third;
	*p++ =-third;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ =-third;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ =-third;
	*p++ =-third;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ =-third;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ =-third;
	*p++ = 1;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ = 1;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p   = 1;
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* called by FormRHS and FormLHS */
void PenaltyContact3DT::LHSDriver(GlobalT::SystemTypeT)
{
	/* time integrator parameters */
	double constK = 0.0;
	int formK = fIntegrator->FormK(constK);
	if (!formK) return;

	/* references to global nodal data */
	const dArray2DT& init_coords = ElementSupport().InitialCoordinates();
	const dArray2DT& disp = Field()[0]; /* displacements */

	/* work space */
	dArrayT c(3), n(3);
	double a[3], b[3];
	iArrayT eqnos;
	
	/* loop over active elements */
	const iArray2DT& connects = *(fConnectivities[0]);
	for (int i = 0; i < connects.MajorDim(); i++)
	{
		const int* pelem = connects(i);

		/* collect element configuration */
		fElRefCoord.RowCollect(pelem, init_coords);
		fElDisp.RowCollect(pelem, disp);

		/* current configuration using effective displacement */
		fElCoord.SumOf(fElRefCoord, fElDisp);
	
		/* get facet and striker coords */
		fElCoord.RowAlias(0, fx1);
		fElCoord.RowAlias(1, fx2);
		fElCoord.RowAlias(2, fx3);
		fElCoord.RowAlias(3, fStriker);

		/* facet normal (direction) = a x b */
		Vector(fx1.Pointer(), fx2.Pointer(), a);
		Vector(fx1.Pointer(), fx3.Pointer(), b);
		CrossProduct(a, b, n.Pointer());
		double mag = sqrt(Dot(n.Pointer(), n.Pointer()));
		n[0] /= mag;
		n[1] /= mag;
		n[2] /= mag;

		/* (store) distance to facet */
		c[0] = fStriker[0] - (fx1[0] + fx2[0] + fx3[0])/3.0;
		c[1] = fStriker[1] - (fx1[1] + fx2[1] + fx3[1])/3.0;
		c[2] = fStriker[2] - (fx1[2] + fx2[2] + fx3[2])/3.0;
		double h = Dot(n.Pointer(), c.Pointer());

		/* contact */
		if (h < 0.0)
		{
			/* initialize */
			fRHS = 0.0;
			fLHS = 0.0;
	
			/* second variation of gap */
			DDg_tri_facet(
				fElRefCoord(0), fElRefCoord(1), fElRefCoord(2), fElRefCoord(3),
				fElDisp(0), fElDisp(1), fElDisp(2), fElDisp(3),
				fLHS);
			int striker_index = fStrikerTags_map.Map(pelem[3]);
			fLHS *= fK*h*constK*fStrikerArea[striker_index];

			/* d_c */
			fdc_du.MultTx(n, fV1);
			fRHS.SetToScaled(1.0, fV1);

			/* d_normal */
			Set_dn_du(fElCoord, fdn_du);					
			fM1.Outer(n, n);
			fM1.PlusIdentity(-1.0);
			fM2.MultATB(fM1, fdn_du);
			fM2.MultTx(c, fV1);
			fRHS.AddScaled(-1.0/mag, fV1);

			/* add term g^T g */
			fLHS.Outer(fRHS, fRHS, fK*constK*fStrikerArea[striker_index], dMatrixT::kAccumulate);

			/* get equation numbers */
			fEqnos[0].RowAlias(i, eqnos);
			
			/* assemble */
			ElementSupport().AssembleLHS(Group(), fLHS, eqnos);
		}
	}
}

void PenaltyContact3DT::RHSDriver(void)
{
	/* time integration parameters */
	double constKd = 0.0;
	int     formKd = fIntegrator->FormKd(constKd);
	if (!formKd) return;

	/* references to global nodal data */
	const dArray2DT& init_coords = ElementSupport().InitialCoordinates();
	const dArray2DT& disp = Field()[0]; /* displacements */

	/* work space */
	dArrayT c(3), n(3);
	iArrayT eqnos;
	double a[3], b[3];

	/* reset tracking data */
	int num_contact = 0;
	double h_max = 0.0;

	/* clear force */
	fStrikerForce2D = 0.0;

	const int* pelem = fConnectivities[0]->Pointer();
	int rowlength = fConnectivities[0]->MinorDim();
	for (int i = 0; i < fConnectivities[0]->MajorDim(); i++, pelem += rowlength)
	{
		/* collect element configuration */
		fElCoord.RowCollect(pelem, init_coords);
		fElDisp.RowCollect(pelem, disp);

		/* current configuration using effective displacement */
		fElCoord.AddScaled(constKd, fElDisp);
	
		/* get facet and striker coords */
		fElCoord.RowAlias(0, fx1);
		fElCoord.RowAlias(1, fx2);
		fElCoord.RowAlias(2, fx3);
		fElCoord.RowAlias(3, fStriker);

		/* facet normal (direction) = a x b */
		Vector(fx1.Pointer(), fx2.Pointer(), a);
		Vector(fx1.Pointer(), fx3.Pointer(), b);
		CrossProduct(a, b, n.Pointer());
		double mag = sqrt(Dot(n.Pointer(), n.Pointer()));
		n[0] /= mag;
		n[1] /= mag;
		n[2] /= mag;

		/* (store) distance to facet */
		c[0] = fStriker[0] - (fx1[0] + fx2[0] + fx3[0])/3.0;
		c[1] = fStriker[1] - (fx1[1] + fx2[1] + fx3[1])/3.0;
		c[2] = fStriker[2] - (fx1[2] + fx2[2] + fx3[2])/3.0;
		double h = Dot(n.Pointer(), c.Pointer());

		/* contact */
		if (h < 0.0)
		{
			/* tracking data */
			num_contact++;
			h_max = (h < h_max) ? h : h_max;
		
			/* penetration force */
			int striker_index = fStrikerTags_map.Map(pelem[3]);			
			double dphi =-fK*h*fStrikerArea[striker_index];

			/* d_c */
			fdc_du.MultTx(n, fV1);
			fRHS.SetToScaled(dphi, fV1);

			/* d_normal */
			Set_dn_du(fElCoord, fdn_du);
			fM1.Outer(n, n);
			fM1.PlusIdentity(-1.0);
			fM2.MultATB(fM1, fdn_du);
			fM2.MultTx(c, fV1);
			fRHS.AddScaled(-dphi/mag, fV1);

			/* ----- Coulomb friction (regularised kinetic) -----
			 * f_t = -mu * |f_n| * v_t / sqrt(|v_t|^2 + eps^2)
			 *
			 * Distribution to nodal RHS:
			 *   striker:  -f_t       (resists striker tangential motion)
			 *   facet 1..3: +f_t / 3 (Newton's 3rd law, equal split)
			 * The 1/3 split is a simplification of the proper shape-function
			 * weighting; preserves linear momentum conservation, sufficient
			 * for kinetic friction in penalty contact.
			 */
			if (fMu > 0.0) {
				const dArray2DT& vel = Field()[1];
				double v_str[3] = {vel(pelem[3], 0), vel(pelem[3], 1), vel(pelem[3], 2)};
				double v_fc[3]  = {
					(vel(pelem[0], 0) + vel(pelem[1], 0) + vel(pelem[2], 0)) / 3.0,
					(vel(pelem[0], 1) + vel(pelem[1], 1) + vel(pelem[2], 1)) / 3.0,
					(vel(pelem[0], 2) + vel(pelem[1], 2) + vel(pelem[2], 2)) / 3.0
				};
				double v_rel[3] = {v_str[0]-v_fc[0], v_str[1]-v_fc[1], v_str[2]-v_fc[2]};
				double vrn = v_rel[0]*n[0] + v_rel[1]*n[1] + v_rel[2]*n[2];
				double vt[3] = {v_rel[0]-vrn*n[0], v_rel[1]-vrn*n[1], v_rel[2]-vrn*n[2]};
				double vt_mag2 = vt[0]*vt[0] + vt[1]*vt[1] + vt[2]*vt[2];
				double inv_smooth = 1.0 / sqrt(vt_mag2 + fFrictionEps*fFrictionEps);
				double f_n_mag = fabs(dphi);
				double s = -fMu * f_n_mag * inv_smooth;  /* opposes v_t */
				double ft[3] = {s*vt[0], s*vt[1], s*vt[2]};

				/* fRHS layout: 4 nodes * 3 dofs = 12 doubles, node order
				 * [facet1, facet2, facet3, striker] */
				double third = 1.0 / 3.0;
				/* facet 1 — opposite reaction */
				fRHS[0]  += -ft[0] * third;
				fRHS[1]  += -ft[1] * third;
				fRHS[2]  += -ft[2] * third;
				/* facet 2 */
				fRHS[3]  += -ft[0] * third;
				fRHS[4]  += -ft[1] * third;
				fRHS[5]  += -ft[2] * third;
				/* facet 3 */
				fRHS[6]  += -ft[0] * third;
				fRHS[7]  += -ft[1] * third;
				fRHS[8]  += -ft[2] * third;
				/* striker — friction acts ON it */
				fRHS[9]  += ft[0];
				fRHS[10] += ft[1];
				fRHS[11] += ft[2];
			}

			/* get equation numbers */
			fEqnos[0].RowAlias(i, eqnos);

			/* assemble */
			ElementSupport().AssembleRHS(Group(), fRHS, eqnos);

			/* store force vector for output */
			fStrikerForce2D(striker_index,0) = dphi*n[0];
			fStrikerForce2D(striker_index,1) = dphi*n[1];
			fStrikerForce2D(striker_index,2) = dphi*n[2];
		}
	}

	/* set tracking */
	SetTrackingData(num_contact, h_max);
}
