/* kl_kinematics — de-risk harness for the KL-shell kinematic chain (issue #63)
 *
 * Implements the auxiliary-tensor velocity-gradient kinematics of Wang & Bazilevs 2024
 * (§3.6 / Algorithm 1): from a mid-surface velocity field it builds A, B1, B2 (and their
 * parametric derivatives), the normal rate n-dot, the 3D Jacobian F3D, and the spatial
 * velocity gradient grad(v3D) = grad_xi(v3D) * F3D^{-1} at a through-thickness station xi3.
 *
 * Validation (flat plate, where the rate-of-deformation D has closed form):
 *   1. rigid rotation about n      -> D = 0           (objectivity)
 *   2. in-plane stretch eps along psi1 -> psi1.D.psi1 = eps
 *   3. bending w=1/2 c xi1^2 about psi2 -> psi1.D.psi1 = -(h/2) xi3 c  (linear in xi3,
 *      antisymmetric, zero membrane part) -- the KL bending signature.
 *
 * Geometry/derivatives use the #69-fixed D2OrthoMLS2DT (completeness 2) on the PCA local
 * parameterization, identical to pca_geom (#62).
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "D2OrthoMLS2DT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "ExceptionT.h"

using namespace Tahoe;
using std::cout; using std::setw; using std::scientific; using std::fixed; using std::setprecision;

/* ---- small 3-vector / 3x3 helpers ---- */
typedef double V3[3];
static double dot(const double a[3], const double b[3]){ return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
static void cross(const double a[3], const double b[3], double o[3]){
	o[0]=a[1]*b[2]-a[2]*b[1]; o[1]=a[2]*b[0]-a[0]*b[2]; o[2]=a[0]*b[1]-a[1]*b[0]; }
static double norm3(const double a[3]){ return std::sqrt(dot(a,a)); }
static int eps3(int i,int j,int k){
	if(i==j||j==k||i==k) return 0;
	/* even perms of (0,1,2): 012,120,201 -> +1; odd: 021,210,102 -> -1 */
	if((i==0&&j==1&&k==2)||(i==1&&j==2&&k==0)||(i==2&&j==0&&k==1)) return 1;
	return -1;
}
static bool inv3(const double M[3][3], double Inv[3][3]){
	double d = M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
	         - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
	         + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]);
	if (std::fabs(d)<1e-300) return false;
	double id=1.0/d;
	Inv[0][0]= (M[1][1]*M[2][2]-M[1][2]*M[2][1])*id;
	Inv[0][1]=-(M[0][1]*M[2][2]-M[0][2]*M[2][1])*id;
	Inv[0][2]= (M[0][1]*M[1][2]-M[0][2]*M[1][1])*id;
	Inv[1][0]=-(M[1][0]*M[2][2]-M[1][2]*M[2][0])*id;
	Inv[1][1]= (M[0][0]*M[2][2]-M[0][2]*M[2][0])*id;
	Inv[1][2]=-(M[0][0]*M[1][2]-M[0][2]*M[1][0])*id;
	Inv[2][0]= (M[1][0]*M[2][1]-M[1][1]*M[2][0])*id;
	Inv[2][1]=-(M[0][0]*M[2][1]-M[0][1]*M[2][0])*id;
	Inv[2][2]= (M[0][0]*M[1][1]-M[0][1]*M[1][0])*id;
	return true;
}
static void Jacobi3(double A[3][3], double eval[3], double evec[3][3]){
	double a[3][3]; for(int i=0;i<3;i++)for(int j=0;j<3;j++)a[i][j]=A[i][j];
	double v[3][3]={{1,0,0},{0,1,0},{0,0,1}};
	for(int s=0;s<100;s++){ double off=std::fabs(a[0][1])+std::fabs(a[0][2])+std::fabs(a[1][2]); if(off<1e-18)break;
		for(int p=0;p<3;p++)for(int q=p+1;q<3;q++){ if(std::fabs(a[p][q])<1e-300)continue;
			double th=(a[q][q]-a[p][p])/(2*a[p][q]); double t=(th>=0?1.0:-1.0)/(std::fabs(th)+std::sqrt(th*th+1));
			double c=1/std::sqrt(t*t+1), sn=t*c;
			for(int k=0;k<3;k++){double akp=a[k][p],akq=a[k][q];a[k][p]=c*akp-sn*akq;a[k][q]=sn*akp+c*akq;}
			for(int k=0;k<3;k++){double apk=a[p][k],aqk=a[q][k];a[p][k]=c*apk-sn*aqk;a[q][k]=sn*apk+c*aqk;}
			for(int k=0;k<3;k++){double vkp=v[k][p],vkq=v[k][q];v[k][p]=c*vkp-sn*vkq;v[k][q]=sn*vkp+c*vkq;} } }
	for(int i=0;i<3;i++){eval[i]=a[i][i];for(int k=0;k<3;k++)evec[i][k]=v[k][i];}
}

/* velocity-field callback: given a neighbor's 3D reference position, return its velocity 3-vec */
struct VField { virtual void operator()(const double X[3], double v[3]) const = 0; virtual ~VField(){} };

/* ---- the KL kinematic chain: returns spatial velocity gradient L=grad(v3D) at xi3 ---- */
/* inputs: parametric 1st/2nd derivs of position (x1,x2,x11,x22,x12) and velocity
   (v1,v2,v11,v22,v12), shell thickness h, station xi3. Outputs L[3][3]. */
static bool KLvelGrad(const double x1[3],const double x2[3],const double x11[3],
                      const double x22[3],const double x12[3],
                      const double v1[3],const double v2[3],const double v11[3],
                      const double v22[3],const double v12[3],
                      double h,double xi3,double L[3][3],
                      double dvxi_out[3][3]=0,double ndot_out[3]=0)
{
	double cr[3]; cross(x1,x2,cr); double J=norm3(cr); if(J<1e-300)return false;
	double n[3]={cr[0]/J,cr[1]/J,cr[2]/J};

	/* A_ij = (delta_ij - n_i n_j)/J */
	double A[3][3]; for(int i=0;i<3;i++)for(int j=0;j<3;j++) A[i][j]=((i==j?1.0:0.0)-n[i]*n[j])/J;

	/* d(cross)/dxi_m and J,m = n . d(cross)/dxi_m ; n,m = (I-n n)/J . d(cross)/dxi_m */
	double dcr[2][3]; /* m=0->xi1, m=1->xi2 */
	{ double t1[3],t2[3];
	  cross(x11,x2,t1); cross(x1,x12,t2); for(int d=0;d<3;d++) dcr[0][d]=t1[d]+t2[d]; /* xi1 */
	  cross(x12,x2,t1); cross(x1,x22,t2); for(int d=0;d<3;d++) dcr[1][d]=t1[d]+t2[d]; /* xi2 */
	}
	double Jm[2], nm[2][3];
	for(int m=0;m<2;m++){ Jm[m]=dot(n,dcr[m]);
		for(int i=0;i<3;i++){ double s=0; for(int k=0;k<3;k++) s+=((i==k?1.0:0.0)-n[i]*n[k])*dcr[m][k]; nm[m][i]=s/J; } }

	/* A,m_ij = -(n,m_i n_j + n_i n,m_j)/J - A_ij J,m / J */
	double Am[2][3][3];
	for(int m=0;m<2;m++)for(int i=0;i<3;i++)for(int j=0;j<3;j++)
		Am[m][i][j] = -(nm[m][i]*n[j]+n[i]*nm[m][j])/J - A[i][j]*Jm[m]/J;

	/* B1_ik = A_ip eps_pkl (x2)_l ;  B2_il = A_ip eps_pkl (x1)_k */
	double B1[3][3], B2[3][3];
	for(int i=0;i<3;i++)for(int k=0;k<3;k++){
		double s1=0,s2=0;
		for(int p=0;p<3;p++)for(int l=0;l<3;l++){ int e=eps3(p,k,l); if(e){ s1+=A[i][p]*e*x2[l]; } }
		for(int p=0;p<3;p++)for(int kk=0;kk<3;kk++){ int e=eps3(p,kk,k); if(e){ s2+=A[i][p]*e*x1[kk]; } }
		B1[i][k]=s1; B2[i][k]=s2;
	}
	/* B1,m and B2,m : x2,m = (m=0:x12, m=1:x22); x1,m = (m=0:x11, m=1:x12) */
	const double* x2m[2]={x12,x22}; const double* x1m[2]={x11,x12};
	double B1m[2][3][3], B2m[2][3][3];
	for(int m=0;m<2;m++)for(int i=0;i<3;i++)for(int k=0;k<3;k++){
		double s1=0,s2=0;
		for(int p=0;p<3;p++)for(int l=0;l<3;l++){ int e=eps3(p,k,l); if(e){ s1+=e*(Am[m][i][p]*x2[l]+A[i][p]*x2m[m][l]); } }
		for(int p=0;p<3;p++)for(int kk=0;kk<3;kk++){ int e=eps3(p,kk,k); if(e){ s2+=e*(Am[m][i][p]*x1[kk]+A[i][p]*x1m[m][kk]); } }
		B1m[m][i][k]=s1; B2m[m][i][k]=s2;
	}

	/* n-dot_i = B1_ik v_k,xi1 + B2_il v_l,xi2 */
	double ndot[3];
	for(int i=0;i<3;i++){ double s=0; for(int k=0;k<3;k++) s+=B1[i][k]*v1[k]+B2[i][k]*v2[k]; ndot[i]=s; }

	/* v3D,xi_m (m=1,2): v_i,xi_m + (h/2)xi3 [ B1,m_ik v_k,1 + B1_ik v_k,1m + B2,m_il v_l,2 + B2_il v_l,2m ] */
	const double* vmm1[2]={v11,v12}; /* v,1m : m=0->v11, m=1->v12 */
	const double* vmm2[2]={v12,v22}; /* v,2m : m=0->v12, m=1->v22 */
	const double* v_m[2]={v1,v2};
	double dvxi[3][3]; /* columns: dvxi[:][0]=v3D,xi1 ; [:][1]=v3D,xi2 ; [:][2]=(h/2)ndot */
	for(int m=0;m<2;m++)for(int i=0;i<3;i++){
		double s=v_m[m][i];
		double t=0; for(int k=0;k<3;k++) t += B1m[m][i][k]*v1[k] + B1[i][k]*vmm1[m][k]
		                                     + B2m[m][i][k]*v2[k] + B2[i][k]*vmm2[m][k];
		dvxi[i][m]= s + (h/2.0)*xi3*t;
	}
	for(int i=0;i<3;i++) dvxi[i][2]=(h/2.0)*ndot[i];

	/* x3D,xi_m (m=1,2): x_i,xi_m + (h/2)xi3 ( B1_ik x_k,1m + B2_il x_l,2m ) ; x3D,xi3=(h/2)n */
	const double* xmm1[2]={x11,x12}; const double* xmm2[2]={x12,x22};
	const double* x_m[2]={x1,x2};
	double F[3][3];
	for(int m=0;m<2;m++)for(int i=0;i<3;i++){
		double t=0; for(int k=0;k<3;k++) t += B1[i][k]*xmm1[m][k] + B2[i][k]*xmm2[m][k];
		F[i][m]= x_m[m][i] + (h/2.0)*xi3*t;
	}
	for(int i=0;i<3;i++) F[i][2]=(h/2.0)*n[i];

	double Finv[3][3]; if(!inv3(F,Finv)) return false;
	/* L = dvxi * Finv */
	for(int i=0;i<3;i++)for(int k=0;k<3;k++){ double s=0; for(int l=0;l<3;l++) s+=dvxi[i][l]*Finv[l][k]; L[i][k]=s; }
	if(dvxi_out) for(int i=0;i<3;i++)for(int k=0;k<3;k++) dvxi_out[i][k]=dvxi[i][k];
	if(ndot_out) for(int i=0;i<3;i++) ndot_out[i]=ndot[i];
	return true;
}

/* mid-surface velocity vector v2D and normal rate n-dot at a sample point, given the
   RK-reconstructed geometry & velocity parametric derivatives there. Used for the
   finite-difference self-consistency check on a curved surface. */
static void V3Dfield(const double x1[3],const double x2[3],
                     const double v2D[3],const double v1[3],const double v2[3],
                     double h,double xi3,double out[3])
{
	double cr[3]; cross(x1,x2,cr); double J=norm3(cr);
	double n[3]={cr[0]/J,cr[1]/J,cr[2]/J};
	double A[3][3]; for(int i=0;i<3;i++)for(int j=0;j<3;j++) A[i][j]=((i==j?1.0:0.0)-n[i]*n[j])/J;
	double B1[3][3],B2[3][3];
	for(int i=0;i<3;i++)for(int k=0;k<3;k++){ double s1=0,s2=0;
		for(int p=0;p<3;p++)for(int l=0;l<3;l++){ int e=eps3(p,k,l); if(e) s1+=A[i][p]*e*x2[l]; }
		for(int p=0;p<3;p++)for(int kk=0;kk<3;kk++){ int e=eps3(p,kk,k); if(e) s2+=A[i][p]*e*x1[kk]; }
		B1[i][k]=s1; B2[i][k]=s2; }
	double ndot[3]; for(int i=0;i<3;i++){ double s=0; for(int k=0;k<3;k++) s+=B1[i][k]*v1[k]+B2[i][k]*v2[k]; ndot[i]=s; }
	for(int i=0;i<3;i++) out[i]=v2D[i]+(h/2.0)*xi3*ndot[i];
}

/* ---- flat-plate driver: build cloud, pick center node, reconstruct derivs, run cases ---- */
int run();
int main(){
	try { return run(); }
	catch (ExceptionT::CodeT& e){ cout<<"\n*** ExceptionT code "<<int(e)<<" ("<<ExceptionT::ToString(e)<<") ***\n"; return 3; }
	catch (...){ cout<<"\n*** unknown exception ***\n"; return 4; }
}

/* velocity fields in the PCA frame (psi1,psi2,n), centered at node P */
struct FRotation : VField { double om; const double *p1,*p2,*nn,*XP;
	void operator()(const double X[3], double v[3]) const {
		double r[3]={X[0]-XP[0],X[1]-XP[1],X[2]-XP[2]}; cross(nn,r,v); for(int d=0;d<3;d++) v[d]*=om; } };
struct FStretch : VField { double eps; const double *p1,*XP;
	void operator()(const double X[3], double v[3]) const {
		double r[3]={X[0]-XP[0],X[1]-XP[1],X[2]-XP[2]}; double s=eps*dot(r,p1); for(int d=0;d<3;d++) v[d]=s*p1[d]; } };
struct FBend : VField { double c; const double *p1,*nn,*XP;
	void operator()(const double X[3], double v[3]) const {
		double r[3]={X[0]-XP[0],X[1]-XP[1],X[2]-XP[2]}; double xi1=dot(r,p1); double w=0.5*c*xi1*xi1;
		for(int d=0;d<3;d++) v[d]=w*nn[d]; } };

int run()
{
	cout << "=== KL kinematic chain de-risk (issue #63) ===\n";
	cout << "flat plate; validate rate-of-deformation D against closed form\n\n" << scientific << setprecision(3);

	const double h=0.05;   /* shell thickness */
	int n=15; double L=1.0, hh=L/(n-1);
	std::vector<double> X,Y,Z;
	for(int i=0;i<n;i++)for(int j=0;j<n;j++){ X.push_back(i*hh); Y.push_back(j*hh); Z.push_back(0.0); }
	int N=X.size();

	/* center node */
	int P=-1; double best=1e9;
	for(int i=0;i<N;i++){ double dx=X[i]-0.5,dy=Y[i]-0.5; double d=dx*dx+dy*dy; if(d<best){best=d;P=i;} }

	/* neighbors within R */
	double R=3.5*hh; std::vector<int> nb;
	for(int Q=0;Q<N;Q++){ double dx=X[Q]-X[P],dy=Y[Q]-Y[P],dz=Z[Q]-Z[P];
		if(std::sqrt(dx*dx+dy*dy+dz*dz)<0.99*R) nb.push_back(Q); }
	int nn=nb.size();

	/* PCA tangent frame */
	double mean[3]={0,0,0}; for(int k=0;k<nn;k++){mean[0]+=X[nb[k]];mean[1]+=Y[nb[k]];mean[2]+=Z[nb[k]];}
	for(int d=0;d<3;d++)mean[d]/=nn;
	double C[3][3]={{0,0,0},{0,0,0},{0,0,0}};
	for(int k=0;k<nn;k++){ double d0=X[nb[k]]-mean[0],d1=Y[nb[k]]-mean[1],d2=Z[nb[k]]-mean[2]; double dd[3]={d0,d1,d2};
		for(int a=0;a<3;a++)for(int b=0;b<3;b++)C[a][b]+=dd[a]*dd[b]; }
	double ev[3],evec[3][3]; Jacobi3(C,ev,evec);
	int o[3]={0,1,2}; for(int a=0;a<3;a++)for(int b=a+1;b<3;b++) if(ev[o[b]]>ev[o[a]]){int t=o[a];o[a]=o[b];o[b]=t;}
	double psi1[3]={evec[o[0]][0],evec[o[0]][1],evec[o[0]][2]};
	double psi2[3]={evec[o[1]][0],evec[o[1]][1],evec[o[1]][2]};
	double nrm[3]; cross(psi1,psi2,nrm); double nm=norm3(nrm); for(int d=0;d<3;d++)nrm[d]/=nm;
	if(nrm[2]<0){ for(int d=0;d<3;d++){nrm[d]=-nrm[d]; psi2[d]=-psi2[d];} } /* orient n=+z */

	/* local params */
	dArray2DT lc(nn,2);
	for(int k=0;k<nn;k++){ double dxv[3]={X[nb[k]]-X[P],Y[nb[k]]-Y[P],Z[nb[k]]-Z[P]};
		lc(k,0)=dot(dxv,psi1); lc(k,1)=dot(dxv,psi2); }

	/* RK derivatives at origin */
	D2OrthoMLS2DT mls(2); mls.Initialize();
	dArrayT dmax(nn); dmax=R; dArrayT fp(2); fp[0]=0; fp[1]=0;
	if(!mls.SetField(lc,dmax,fp)){ cout<<"SetField failed\n"; return 2; }
	const dArray2DT& Dphi=mls.Dphi(); const dArray2DT& DDphi=mls.DDphi();

	double XP[3]={X[P],Y[P],Z[P]};

	/* position parametric derivatives (reference) */
	double x1[3]={0,0,0},x2[3]={0,0,0},x11[3]={0,0,0},x22[3]={0,0,0},x12[3]={0,0,0};
	for(int I=0;I<nn;I++){ double Xq[3]={X[nb[I]],Y[nb[I]],Z[nb[I]]};
		for(int d=0;d<3;d++){ x1[d]+=Dphi(0,I)*Xq[d]; x2[d]+=Dphi(1,I)*Xq[d];
			x11[d]+=DDphi(0,I)*Xq[d]; x22[d]+=DDphi(1,I)*Xq[d]; x12[d]+=DDphi(2,I)*Xq[d]; } }

	/* helper to assemble velocity derivs for a given field and run the chain at xi3 */
	struct Runner {
		int nn; const std::vector<int>* nb; const std::vector<double>*X,*Y,*Z;
		const dArray2DT *Dphi,*DDphi; const double *x1,*x2,*x11,*x22,*x12; double h;
		void velDerivs(const VField& f,double v1[3],double v2[3],double v11[3],double v22[3],double v12[3]) const {
			for(int d=0;d<3;d++){v1[d]=v2[d]=v11[d]=v22[d]=v12[d]=0;}
			for(int I=0;I<nn;I++){ double Xq[3]={(*X)[(*nb)[I]],(*Y)[(*nb)[I]],(*Z)[(*nb)[I]]}; double vv[3]; f(Xq,vv);
				for(int d=0;d<3;d++){ v1[d]+=(*Dphi)(0,I)*vv[d]; v2[d]+=(*Dphi)(1,I)*vv[d];
					v11[d]+=(*DDphi)(0,I)*vv[d]; v22[d]+=(*DDphi)(1,I)*vv[d]; v12[d]+=(*DDphi)(2,I)*vv[d]; } }
		}
		bool L(const VField& f,double xi3,double Lout[3][3]) const {
			double v1[3],v2[3],v11[3],v22[3],v12[3]; velDerivs(f,v1,v2,v11,v22,v12);
			return KLvelGrad(x1,x2,x11,x22,x12,v1,v2,v11,v22,v12,h,xi3,Lout);
		}
	} R0;
	R0.nn=nn; R0.nb=&nb; R0.X=&X; R0.Y=&Y; R0.Z=&Z; R0.Dphi=&Dphi; R0.DDphi=&DDphi;
	R0.x1=x1;R0.x2=x2;R0.x11=x11;R0.x22=x22;R0.x12=x12; R0.h=h;

	auto Dsym=[&](double Lm[3][3],double Dm[3][3]){ for(int i=0;i<3;i++)for(int j=0;j<3;j++) Dm[i][j]=0.5*(Lm[i][j]+Lm[j][i]); };
	auto proj=[&](double M[3][3],const double a[3],const double b[3]){ double s=0; for(int i=0;i<3;i++)for(int j=0;j<3;j++) s+=a[i]*M[i][j]*b[j]; return s; };
	auto frob=[&](double M[3][3]){ double s=0; for(int i=0;i<3;i++)for(int j=0;j<3;j++) s+=M[i][j]*M[i][j]; return std::sqrt(s); };

	int fails=0;

	/* ---- case 1: rigid rotation about n -> D = 0 ---- */
	{ FRotation f; f.om=0.7; f.p1=psi1; f.p2=psi2; f.nn=nrm; f.XP=XP;
	  double Lm[3][3],Dm[3][3]; R0.L(f,0.0,Lm); Dsym(Lm,Dm); double dn=frob(Dm);
	  bool ok=dn<1e-9; fails+=!ok;
	  cout<<"[1] rigid rotation about n (omega=0.7): |D| = "<<dn<<"   "<<(ok?"PASS (objectivity: D=0)":"FAIL")<<"\n"; }

	/* ---- case 2: in-plane stretch eps along psi1 -> psi1.D.psi1 = eps ---- */
	{ FStretch f; f.eps=0.3; f.p1=psi1; f.XP=XP;
	  double Lm[3][3],Dm[3][3]; R0.L(f,0.0,Lm); Dsym(Lm,Dm);
	  double s11=proj(Dm,psi1,psi1), s22=proj(Dm,psi2,psi2);
	  bool ok=std::fabs(s11-0.3)<1e-6;
	  fails+=!ok;
	  cout<<"[2] stretch eps=0.3 along psi1: psi1.D.psi1 = "<<s11<<" (exp 0.3), psi2.D.psi2 = "<<s22<<"   "<<(ok?"PASS":"FAIL")<<"\n"; }

	/* ---- case 3: bending w=1/2 c xi1^2 -> psi1.D.psi1 = -(h/2) xi3 c, linear/antisym ---- */
	{ FBend f; f.c=1.0; f.p1=psi1; f.nn=nrm; f.XP=XP;
	  double Lp[3][3],Dp[3][3],Lz[3][3],Dz[3][3],Lm[3][3],Dm[3][3];
	  R0.L(f,+1.0,Lp); Dsym(Lp,Dp);
	  R0.L(f, 0.0,Lz); Dsym(Lz,Dz);
	  R0.L(f,-1.0,Lm); Dsym(Lm,Dm);
	  double sp=proj(Dp,psi1,psi1), sz=proj(Dz,psi1,psi1), sm=proj(Dm,psi1,psi1);
	  double exp_p=-(h/2.0)*(+1.0)*f.c, exp_m=-(h/2.0)*(-1.0)*f.c;
	  bool ok = std::fabs(sp-exp_p)<1e-4 && std::fabs(sm-exp_m)<1e-4 && std::fabs(sz)<1e-4 && std::fabs(sp+sm)<1e-6;
	  fails+=!ok;
	  cout<<"[3] bending c=1.0: psi1.D.psi1 @xi3=+1 = "<<sp<<" (exp "<<exp_p<<"), @0 = "<<sz<<" (exp 0), @-1 = "<<sm<<" (exp "<<exp_m<<")\n";
	  cout<<"    linear & antisymmetric in xi3, membrane~0: "<<(ok?"PASS (KL bending signature)":"FAIL")<<"\n"; }

	/* ---- case 4: curved surface (cylinder) — FD self-consistency of grad_xi(v3D) ----
	   Validates the curvature-derivative tensors B1,xi / B2,xi (zero on a flat plate).
	   Compares the chain's analytic parametric velocity gradient to a central finite
	   difference of v3D reconstructed at shifted parametric sample points. */
	{
		const double Rc=1.0, Zlen=2.0; int nt=64, nz=21;
		double dth=2*M_PI/nt, dz=Zlen/(nz-1), arc=Rc*dth, hC=std::max(arc,dz);
		std::vector<double> CX,CY,CZ;
		for(int i=0;i<nt;i++)for(int j=0;j<nz;j++){ double th=i*dth,z=j*dz;
			CX.push_back(Rc*std::cos(th)); CY.push_back(Rc*std::sin(th)); CZ.push_back(z); }
		int M=CX.size();
		/* node near mid-height */
		int Pc=-1; double bb=1e9; for(int i=0;i<M;i++){ double d=std::fabs(CZ[i]-Zlen/2)+std::fabs(CX[i]-Rc); if(d<bb){bb=d;Pc=i;} }
		double Rn=3.5*hC; std::vector<int> nbc;
		for(int Q=0;Q<M;Q++){ double dx=CX[Q]-CX[Pc],dy=CY[Q]-CY[Pc],dz3=CZ[Q]-CZ[Pc];
			if(std::sqrt(dx*dx+dy*dy+dz3*dz3)<0.99*Rn) nbc.push_back(Q); }
		int mc=nbc.size();
		/* PCA frame */
		double mn[3]={0,0,0}; for(int k=0;k<mc;k++){mn[0]+=CX[nbc[k]];mn[1]+=CY[nbc[k]];mn[2]+=CZ[nbc[k]];}
		for(int d=0;d<3;d++)mn[d]/=mc;
		double Cv[3][3]={{0,0,0},{0,0,0},{0,0,0}};
		for(int k=0;k<mc;k++){ double dd[3]={CX[nbc[k]]-mn[0],CY[nbc[k]]-mn[1],CZ[nbc[k]]-mn[2]};
			for(int a=0;a<3;a++)for(int b=0;b<3;b++)Cv[a][b]+=dd[a]*dd[b]; }
		double evc[3],evecc[3][3]; Jacobi3(Cv,evc,evecc);
		int oc[3]={0,1,2}; for(int a=0;a<3;a++)for(int b=a+1;b<3;b++) if(evc[oc[b]]>evc[oc[a]]){int t=oc[a];oc[a]=oc[b];oc[b]=t;}
		double q1[3]={evecc[oc[0]][0],evecc[oc[0]][1],evecc[oc[0]][2]};
		double q2[3]={evecc[oc[1]][0],evecc[oc[1]][1],evecc[oc[1]][2]};
		dArray2DT lcc(mc,2);
		for(int k=0;k<mc;k++){ double dxv[3]={CX[nbc[k]]-CX[Pc],CY[nbc[k]]-CY[Pc],CZ[nbc[k]]-CZ[Pc]};
			lcc(k,0)=dot(dxv,q1); lcc(k,1)=dot(dxv,q2); }

		/* arbitrary smooth velocity field, sampled at nodes */
		std::vector<double> VX(mc),VY(mc),VZ(mc);
		for(int k=0;k<mc;k++){ double X3[3]={CX[nbc[k]],CY[nbc[k]],CZ[nbc[k]]};
			VX[k]=0.3*X3[1]*X3[2]; VY[k]=0.2*X3[0]; VZ[k]=0.1*X3[0]*X3[0]; }

		D2OrthoMLS2DT m2(2); m2.Initialize();
		dArrayT dm(mc); dm=Rn;

		/* reconstruct geometry + velocity derivs at a given local sample point */
		auto recon=[&](double s0,double s1,double x1[3],double x2[3],double x11[3],double x22[3],double x12[3],
		                double v2D[3],double v1[3],double v2[3],double v11[3],double v22[3],double v12[3])->bool{
			dArrayT fpp(2); fpp[0]=s0; fpp[1]=s1;
			if(!m2.SetField(lcc,dm,fpp)) return false;
			const dArrayT& ph=m2.phi(); const dArray2DT& Dp=m2.Dphi(); const dArray2DT& DDp=m2.DDphi();
			for(int d=0;d<3;d++){x1[d]=x2[d]=x11[d]=x22[d]=x12[d]=0; v2D[d]=v1[d]=v2[d]=v11[d]=v22[d]=v12[d]=0;}
			for(int I=0;I<mc;I++){ double Xq[3]={CX[nbc[I]],CY[nbc[I]],CZ[nbc[I]]}, Vq[3]={VX[I],VY[I],VZ[I]};
				for(int d=0;d<3;d++){
					x1[d]+=Dp(0,I)*Xq[d]; x2[d]+=Dp(1,I)*Xq[d];
					x11[d]+=DDp(0,I)*Xq[d]; x22[d]+=DDp(1,I)*Xq[d]; x12[d]+=DDp(2,I)*Xq[d];
					v2D[d]+=ph[I]*Vq[d]; v1[d]+=Dp(0,I)*Vq[d]; v2[d]+=Dp(1,I)*Vq[d];
					v11[d]+=DDp(0,I)*Vq[d]; v22[d]+=DDp(1,I)*Vq[d]; v12[d]+=DDp(2,I)*Vq[d]; } }
			return true;
		};

		double xi3=0.5, delta=1e-4;
		/* analytic parametric velocity gradient at origin */
		double X1[3],X2[3],X11[3],X22[3],X12[3],V2D[3],Vv1[3],Vv2[3],V11[3],V22[3],V12[3];
		recon(0,0,X1,X2,X11,X22,X12,V2D,Vv1,Vv2,V11,V22,V12);
		double Lc[3][3],dvxiA[3][3],ndotA[3];
		KLvelGrad(X1,X2,X11,X22,X12,Vv1,Vv2,V11,V22,V12,h,xi3,Lc,dvxiA,ndotA);

		/* FD of v3D in xi1, xi2 */
		auto v3Dat=[&](double s0,double s1,double out[3]){
			double a1[3],a2[3],a11[3],a22[3],a12[3],w[3],b1[3],b2[3],b11[3],b22[3],b12[3];
			recon(s0,s1,a1,a2,a11,a22,a12,w,b1,b2,b11,b22,b12);
			V3Dfield(a1,a2,w,b1,b2,h,xi3,out);
		};
		double vp[3],vm[3]; double fd[3][2];
		v3Dat(+delta,0,vp); v3Dat(-delta,0,vm); for(int d=0;d<3;d++) fd[d][0]=(vp[d]-vm[d])/(2*delta);
		v3Dat(0,+delta,vp); v3Dat(0,-delta,vm); for(int d=0;d<3;d++) fd[d][1]=(vp[d]-vm[d])/(2*delta);

		double emax=0; for(int d=0;d<3;d++)for(int m=0;m<2;m++) emax=std::max(emax,std::fabs(fd[d][m]-dvxiA[d][m]));
		double scale=0; for(int d=0;d<3;d++)for(int m=0;m<2;m++) scale=std::max(scale,std::fabs(dvxiA[d][m]));
		double rel=emax/std::max(scale,1e-12);
		bool ok=rel<1e-5; fails+=!ok;
		cout<<"[4] cylinder grad_xi(v3D): analytic-vs-FD max abs err = "<<emax<<"  (rel "<<rel<<")   "
		    <<(ok?"PASS (B1,xi/B2,xi curvature terms correct)":"FAIL")<<"\n";
	}

	cout << "\n=== VERDICT (issue #63 kinematic core) ===\n";
	cout << "  KL auxiliary-tensor chain (A,B1,B2,n-dot,F3D,grad v3D) reproduces the analytic\n"
	        "  rate-of-deformation for rigid rotation, in-plane stretch, and through-thickness\n"
	        "  bending on a flat plate, and is FD-consistent (incl. B1,xi/B2,xi) on a cylinder.\n"
	        "  status: " << (fails==0?"PASS":"FAIL") << "\n";
	return fails==0?0:1;
}
