/* FlanaganTaylor.h — Flanagan & Taylor (1987) incremental rotation/stretch integration
 * for the Green-Naghdi objective co-rotational stress update (Wang & Bazilevs 2025, Algorithm 2).
 *
 * Advances the proper orthogonal rotation R and the left-stretch V over a step given the
 * velocity-gradient increment dL = L*dt (L = grad v in the current configuration). The material
 * response is then integrated in the co-rotational frame: D~ = R^T D R, sigma~ += dt*F(D~),
 * sigma = R sigma~ R^T -- this is what makes the stress objective under large rotation.
 *
 * Conventions:
 *   axial vector of a skew matrix S:  S x = s x x  ->  s = (S32, S13, S21)
 *   skew(s)_ij = eps_ikj s_k  ->  the matrix with that axial vector.
 */
#ifndef _FLANAGAN_TAYLOR_H_
#define _FLANAGAN_TAYLOR_H_

#include <cmath>

namespace Tahoe {
namespace KLShell {

/* C = A*B (3x3) */
inline void M3_mul(const double A[3][3], const double B[3][3], double C[3][3])
{
	for (int i=0;i<3;i++) for (int j=0;j<3;j++){
		double s=0.0; for(int k=0;k<3;k++) s+=A[i][k]*B[k][j]; C[i][j]=s; }
}

/* C = A*B^T (3x3) */
inline void M3_mulT(const double A[3][3], const double B[3][3], double C[3][3])
{
	for (int i=0;i<3;i++) for (int j=0;j<3;j++){
		double s=0.0; for(int k=0;k<3;k++) s+=A[i][k]*B[j][k]; C[i][j]=s; }
}

/* C = A^T*B (3x3) */
inline void M3_TmulT(const double A[3][3], const double B[3][3], double C[3][3])
{
	for (int i=0;i<3;i++) for (int j=0;j<3;j++){
		double s=0.0; for(int k=0;k<3;k++) s+=A[k][i]*B[k][j]; C[i][j]=s; }
}

/* 3x3 inverse via cofactors; returns det (0 -> singular, Ainv left unset) */
inline double M3_inv(const double A[3][3], double Ai[3][3])
{
	double c00=A[1][1]*A[2][2]-A[1][2]*A[2][1];
	double c01=A[1][2]*A[2][0]-A[1][0]*A[2][2];
	double c02=A[1][0]*A[2][1]-A[1][1]*A[2][0];
	double det=A[0][0]*c00+A[0][1]*c01+A[0][2]*c02;
	if (std::fabs(det)<1.0e-300) return 0.0;
	double id=1.0/det;
	Ai[0][0]=c00*id;
	Ai[0][1]=(A[0][2]*A[2][1]-A[0][1]*A[2][2])*id;
	Ai[0][2]=(A[0][1]*A[1][2]-A[0][2]*A[1][1])*id;
	Ai[1][0]=c01*id;
	Ai[1][1]=(A[0][0]*A[2][2]-A[0][2]*A[2][0])*id;
	Ai[1][2]=(A[0][2]*A[1][0]-A[0][0]*A[1][2])*id;
	Ai[2][0]=c02*id;
	Ai[2][1]=(A[0][1]*A[2][0]-A[0][0]*A[2][1])*id;
	Ai[2][2]=(A[0][0]*A[1][1]-A[0][1]*A[1][0])*id;
	return det;
}

/* Flanagan-Taylor step. dL = velocity-gradient increment (L*dt). R, V updated in place.
 * R must be a proper rotation on entry (init = reference frame), V the left stretch (init = I). */
inline void FlanaganTaylorStep(const double dL[3][3], double R[3][3], double V[3][3])
{
	/* 1) D = sym(dL), W = skew(dL) */
	double D[3][3], W[3][3];
	for (int i=0;i<3;i++) for (int j=0;j<3;j++){
		D[i][j]=0.5*(dL[i][j]+dL[j][i]); W[i][j]=0.5*(dL[i][j]-dL[j][i]); }

	/* axial vector of W:  w = (W32, W13, W21) */
	double w[3] = { W[2][1], W[0][2], W[1][0] };

	/* 2) z_i = eps_ikj D_jl V_lk = (DV)_cyc - (DV)_anticyc;  z = 2*axial(skew(DV)) */
	double DV[3][3]; M3_mul(D,V,DV);
	double z[3] = { DV[2][1]-DV[1][2], DV[0][2]-DV[2][0], DV[1][0]-DV[0][1] };

	/* [I tr(V) - V]^{-1} z */
	double trV=V[0][0]+V[1][1]+V[2][2];
	double M[3][3];
	for (int i=0;i<3;i++) for (int j=0;j<3;j++) M[i][j]=(i==j?trV:0.0)-V[i][j];
	double Mi[3][3]; double det=M3_inv(M,Mi);
	double q[3];
	if (det!=0.0) for(int i=0;i<3;i++){ double s=0.0; for(int k=0;k<3;k++) s+=Mi[i][k]*z[k]; q[i]=w[i]+s; }
	else          for(int i=0;i<3;i++) q[i]=w[i];

	/* Omega = skew(q):  Omega_ij = eps_ikj q_k */
	double Om[3][3] = {{0.0,-q[2],q[1]},{q[2],0.0,-q[0]},{-q[1],q[0],0.0}};

	/* 3) Q = (I - Om/2)^{-1}(I + Om/2);  R <- Q R */
	double Ap[3][3], Am[3][3];
	for (int i=0;i<3;i++) for (int j=0;j<3;j++){
		double e=(i==j?1.0:0.0); Ap[i][j]=e+0.5*Om[i][j]; Am[i][j]=e-0.5*Om[i][j]; }
	double Ami[3][3]; M3_inv(Am,Ami);
	double Q[3][3]; M3_mul(Ami,Ap,Q);
	double Rn[3][3]; M3_mul(Q,R,Rn);
	for (int i=0;i<3;i++) for (int j=0;j<3;j++) R[i][j]=Rn[i][j];

	/* 4-5) Vdot = (D+W)V - V*Omega ;  V <- V + Vdot */
	double DW[3][3]; for(int i=0;i<3;i++)for(int j=0;j<3;j++) DW[i][j]=D[i][j]+W[i][j];
	double DWV[3][3]; M3_mul(DW,V,DWV);
	double VOm[3][3]; M3_mul(V,Om,VOm);
	for (int i=0;i<3;i++) for (int j=0;j<3;j++) V[i][j]+=DWV[i][j]-VOm[i][j];
}

} /* namespace KLShell */
} /* namespace Tahoe */

#endif /* _FLANAGAN_TAYLOR_H_ */
