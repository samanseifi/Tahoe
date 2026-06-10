/* PlaneStressJ2.h — plane-stress J2 (von Mises) radial-return stress update for the meshfree
 * Kirchhoff-Love shell (epic #59, Fig 18 elasto-plastic track).
 *
 * Each through-thickness Gauss point of the shell is a PLANE-STRESS material point (sigma33 = 0).
 * A 3D J2 model cannot be used directly: the sigma33 = 0 constraint couples into the return, so
 * the deviatoric flow is non-radial in the in-plane stress space. This is the standard
 * Simo & Hughes plane-stress return, written as a stress update only (no consistent tangent --
 * the explicit central-difference solver needs only the updated stress).
 *
 * In-plane Voigt convention: s = [s11, s22, s12], strain = [e11, e22, gamma12] (ENGINEERING shear).
 * von Mises in plane stress: seq^2 = s11^2 - s11*s22 + s22^2 + 3*s12^2 = s^T P s, with
 *   P = [[1, -1/2, 0], [-1/2, 1, 0], [0, 0, 3]].
 * Isotropic linear hardening: Y(ep) = Y0 + H*ep.
 */
#ifndef _PLANE_STRESS_J2_H_
#define _PLANE_STRESS_J2_H_

#include <cmath>
#include "SecantMethodT.h"

namespace Tahoe {
namespace KLShell {

/* seq^2 = s^T P s */
inline double J2pq(const double s[3])
{
	return s[0]*s[0] - s[0]*s[1] + s[1]*s[1] + 3.0*s[2]*s[2];
}

/* flow stress Y(ep) = Y0 + H*ep + (Ysat-Y0)*(1-exp(-delta*ep)).  Linear hardening: Ysat<=Y0 or
 * delta<=0 (the saturation term vanishes). Necking (paper Fig 12): Y0=343,H=300,Ysat=680,delta=16.93. */
inline double J2yield(double ep, double Y0, double H, double Ysat, double delta)
{
	double Y = Y0 + H*ep;
	if (Ysat > Y0 && delta > 0.0) Y += (Ysat - Y0)*(1.0 - std::exp(-delta*ep));
	return Y;
}

/* solve [I + b*C*P] x = rhs for the in-plane stress (3x3, P with engineering-shear scaling) */
inline void J2solve(double b, double c, double nu, const double rhs[3], double x[3])
{
	/* C (plane stress, engineering shear) = c*[[1,nu,0],[nu,1,0],[0,0,(1-nu)/2]], c=E/(1-nu^2).
	 * M = I + b*C*P, with P=[[1,-1/2,0],[-1/2,1,0],[0,0,3]]. */
	double cs = c*(1.0-nu)/2.0;
	double CP[3][3] = {
		{ c*(1.0 - 0.5*nu),  c*(nu - 0.5),     0.0 },
		{ c*(nu - 0.5),      c*(1.0 - 0.5*nu), 0.0 },
		{ 0.0,               0.0,              3.0*cs } };
	double M[3][3];
	for (int i=0;i<3;i++) for (int j=0;j<3;j++) M[i][j] = (i==j?1.0:0.0) + b*CP[i][j];
	/* the 2x2 (11,22) block is decoupled from shear (12); solve directly */
	double det = M[0][0]*M[1][1] - M[0][1]*M[1][0];
	x[0] = ( M[1][1]*rhs[0] - M[0][1]*rhs[1]) / det;
	x[1] = (-M[1][0]*rhs[0] + M[0][0]*rhs[1]) / det;
	x[2] = rhs[2] / M[2][2];
}

/* Plane-stress J2 radial return. In/out: sig (in-plane stress), ep (equivalent plastic strain).
 * In: deps (in-plane engineering strain increment), E, nu, Y0 (initial yield), H (hardening mod). */
inline void PlaneStressJ2Return(double sig[3], const double deps[3], double& ep,
                                double E, double nu, double Y0, double H,
                                double Ysat=0.0, double delta=0.0)
{
	double c = E/(1.0 - nu*nu);
	/* elastic trial stress */
	double st[3];
	st[0] = sig[0] + c*(deps[0] + nu*deps[1]);
	st[1] = sig[1] + c*(nu*deps[0] + deps[1]);
	st[2] = sig[2] + c*(1.0-nu)/2.0*deps[2];

	double seq_tr = std::sqrt(J2pq(st));
	double Yn = J2yield(ep, Y0, H, Ysat, delta);
	if (seq_tr <= Yn || seq_tr < 1.0e-30) {       /* elastic step */
		sig[0]=st[0]; sig[1]=st[1]; sig[2]=st[2];
		return;
	}

	/* plastic: find dgamma so that seq(sig(dgamma)) = Y(ep + dgamma). At a candidate dgamma the
	 * consistency target Yc = Y0 + H*(ep+dgamma) sets b = dgamma/Yc and sig = [I + b C P]^-1 st;
	 * the residual r = seq(sig) - Yc is monotone -> bracket + bisection (robust, derivative-free). */
	double lo = 0.0, hi = (seq_tr - Yn)/ (1.5*E/(1.0+nu) + (H>0.0?H:0.0)); /* elastic-predictor guess */
	if (hi <= 0.0) hi = seq_tr/E + 1.0e-8;
	/* expand hi until r(hi) < 0 */
	double sg[3];
	for (int it=0; it<60; it++) {
		double Yc = J2yield(ep+hi, Y0, H, Ysat, delta);
		J2solve(hi/Yc, c, nu, st, sg);
		if (std::sqrt(J2pq(sg)) - Yc < 0.0) break;
		hi *= 2.0;
	}
	double dgamma = 0.5*(lo+hi);
	for (int it=0; it<80; it++) {
		dgamma = 0.5*(lo+hi);
		double Yc = J2yield(ep+dgamma, Y0, H, Ysat, delta);
		J2solve(dgamma/Yc, c, nu, st, sg);
		double r = std::sqrt(J2pq(sg)) - Yc;
		if (r > 0.0) lo = dgamma; else hi = dgamma;
		if (hi-lo < 1.0e-14*(1.0+hi)) break;
	}
	sig[0]=sg[0]; sig[1]=sg[1]; sig[2]=sg[2];
	ep += dgamma;
}

/* ---------------------------------------------------------------------------------------------- *
 * Paper Algorithm 3: enforce sigma33=0 with a FULL 3D J2 law + Newton/secant on the through-
 * thickness rate of deformation D33. Unlike the 2D plane-stress condensation above, this returns
 * the actual D33 (through-thickness strain increment) -> the physical thickness change that drives
 * necking localization.
 * ---------------------------------------------------------------------------------------------- */

/* 3D J2 radial return. sig[6] Voigt (11,22,33,23,13,12) tensor shear; deps[6] engineering shear. */
inline void J2_3D_return(double sig[6], const double deps[6], double& ep,
                         double E, double nu, double Y0, double H, double Ysat, double delta)
{
	double lam = E*nu/((1.0+nu)*(1.0-2.0*nu)), G = E/(2.0*(1.0+nu));
	double tr = deps[0]+deps[1]+deps[2];
	double st[6];
	st[0]=sig[0]+lam*tr+2.0*G*deps[0]; st[1]=sig[1]+lam*tr+2.0*G*deps[1]; st[2]=sig[2]+lam*tr+2.0*G*deps[2];
	st[3]=sig[3]+G*deps[3]; st[4]=sig[4]+G*deps[4]; st[5]=sig[5]+G*deps[5];   /* engineering shear -> G, not 2G */
	double p=(st[0]+st[1]+st[2])/3.0;
	double s[6]={st[0]-p,st[1]-p,st[2]-p,st[3],st[4],st[5]};
	double q=std::sqrt(1.5*(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]+2.0*(s[3]*s[3]+s[4]*s[4]+s[5]*s[5])));
	double Y=J2yield(ep,Y0,H,Ysat,delta);
	if (q<=Y || q<1.0e-30) { for(int i=0;i<6;i++) sig[i]=st[i]; return; }
	/* dgamma: q - 3G dgamma = Y(ep+dgamma), monotone -> bisection */
	double lo=0.0, hi=(q-Y)/(3.0*G)+1.0e-12, dg=0.0;
	for (int it=0;it<60;it++){ dg=0.5*(lo+hi); double Yc=J2yield(ep+dg,Y0,H,Ysat,delta);
		if (q-3.0*G*dg-Yc>0.0) lo=dg; else hi=dg; if (hi-lo<1.0e-14*(1.0+hi)) break; }
	double fac=1.0-3.0*G*dg/q;
	sig[0]=fac*s[0]+p; sig[1]=fac*s[1]+p; sig[2]=fac*s[2]+p;
	sig[3]=fac*s[3]; sig[4]=fac*s[4]; sig[5]=fac*s[5];
	ep+=dg;
}

/* sigma33=0 via secant on D33. sig[3]=in-plane (s11,s22,s12); dip[3]=in-plane strain incr (de11,de22,dg12).
 * Returns de33 (through-thickness strain increment); updates sig (in-plane) + ep. */
inline double PlaneStressJ2_D33(double sig[3], const double dip[3], double& ep,
                                double E, double nu, double Y0, double H, double Ysat, double delta)
{
	double s_in[6]={sig[0],sig[1],0.0,0.0,0.0,sig[2]};   /* stored 3D stress (sigma33=0 from last step) */
	/* residual r(de33) = sigma33 after the 3D return; bracket two guesses (elastic plane-stress + perturb) */
	double g0 = -nu/(1.0-nu)*(dip[0]+dip[1]);            /* elastic plane-stress de33 */
	double g1 = g0 - 0.3*(std::fabs(dip[0])+std::fabs(dip[1])) - 1.0e-9;
	double sa[6]; for(int i=0;i<6;i++) sa[i]=s_in[i]; double epa=ep; double da[6]={dip[0],dip[1],g0,0,0,dip[2]};
	J2_3D_return(sa,da,epa,E,nu,Y0,H,Ysat,delta); double ea=sa[2];
	double sb[6]; for(int i=0;i<6;i++) sb[i]=s_in[i]; double epb=ep; double db[6]={dip[0],dip[1],g1,0,0,dip[2]};
	J2_3D_return(sb,db,epb,E,nu,Y0,H,Ysat,delta); double eb=sb[2];
	SecantMethodT secant(30, 1.0e-10);
	secant.Reset(g0,ea,g1,eb);
	double de33=g0, epf=epa; double sf[6]; for(int i=0;i<6;i++) sf[i]=sa[i];
	for (int it=0; it<30; it++){
		de33=secant.NextGuess();
		double sc[6]; for(int i=0;i<6;i++) sc[i]=s_in[i]; double epc=ep; double dc[6]={dip[0],dip[1],de33,0,0,dip[2]};
		J2_3D_return(sc,dc,epc,E,nu,Y0,H,Ysat,delta);
		SecantMethodT::StatusT st=secant.NextPoint(de33,sc[2]);
		epf=epc; for(int i=0;i<6;i++) sf[i]=sc[i];
		if (st==SecantMethodT::kConverged || st==SecantMethodT::kFail) break;
	}
	sig[0]=sf[0]; sig[1]=sf[1]; sig[2]=sf[5]; ep=epf;
	return de33;
}

} /* namespace KLShell */
} /* namespace Tahoe */

#endif /* _PLANE_STRESS_J2_H_ */
