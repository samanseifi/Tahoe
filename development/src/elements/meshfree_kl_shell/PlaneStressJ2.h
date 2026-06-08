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

namespace Tahoe {
namespace KLShell {

/* seq^2 = s^T P s */
inline double J2pq(const double s[3])
{
	return s[0]*s[0] - s[0]*s[1] + s[1]*s[1] + 3.0*s[2]*s[2];
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
                                double E, double nu, double Y0, double H)
{
	double c = E/(1.0 - nu*nu);
	/* elastic trial stress */
	double st[3];
	st[0] = sig[0] + c*(deps[0] + nu*deps[1]);
	st[1] = sig[1] + c*(nu*deps[0] + deps[1]);
	st[2] = sig[2] + c*(1.0-nu)/2.0*deps[2];

	double seq_tr = std::sqrt(J2pq(st));
	double Yn = Y0 + H*ep;
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
		double Yc = Y0 + H*(ep+hi);
		J2solve(hi/Yc, c, nu, st, sg);
		if (std::sqrt(J2pq(sg)) - Yc < 0.0) break;
		hi *= 2.0;
	}
	double dgamma = 0.5*(lo+hi);
	for (int it=0; it<80; it++) {
		dgamma = 0.5*(lo+hi);
		double Yc = Y0 + H*(ep+dgamma);
		J2solve(dgamma/Yc, c, nu, st, sg);
		double r = std::sqrt(J2pq(sg)) - Yc;
		if (r > 0.0) lo = dgamma; else hi = dgamma;
		if (hi-lo < 1.0e-14*(1.0+hi)) break;
	}
	sig[0]=sg[0]; sig[1]=sg[1]; sig[2]=sg[2];
	ep += dgamma;
}

} /* namespace KLShell */
} /* namespace Tahoe */

#endif /* _PLANE_STRESS_J2_H_ */
