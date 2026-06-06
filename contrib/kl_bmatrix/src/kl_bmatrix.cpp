/* kl_bmatrix — de-risk harness for the KL-shell strain-displacement matrices (issue #64)
 *
 * Implements the per-node shell B-matrix B_Iijk (Wang & Bazilevs 2024 Eq. 51) such that
 * the spatial velocity gradient is  grad(v3D)_ij = sum_I B_Iijk v_Ik, its Voigt form
 * B_I^V (Eq. 52) giving the strain rate, and the parametric gradient B_Iijkl at xi3=0
 * (Eq. 54, used by the naturally-stabilized nodal-integration term, Eq. 57).
 *
 * Validation (against the #63-validated kinematic chain, independent code path):
 *   A. sum_I B_Iijk v_Ik  ==  grad(v3D)  from KLvelGrad     (flat plate & cylinder, machine prec.)
 *   B. Voigt B_I^V . v_I  ==  Voigt(sym(grad v3D))          (strain-rate form)
 *   C. sum_I B_Iijkl v_Ik ==  d(grad v3D)_ij/dxi_l |xi3=0   (cylinder, finite difference)
 *
 * Geometry/derivatives via the #69-fixed D2OrthoMLS2DT on the PCA local parameterization.
 *
 * NOTE: the full stabilized internal force (Eq. 57) couples B and B,xi with the Cauchy
 * stress and its parametric gradient; that assembly is exercised once constitutive +
 * through-thickness integration land (#65). Here we lock down the kinematic operators.
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
using std::cout; using std::setw; using std::scientific; using std::setprecision;

static double dot(const double a[3], const double b[3]){ return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
static void cross(const double a[3], const double b[3], double o[3]){
	o[0]=a[1]*b[2]-a[2]*b[1]; o[1]=a[2]*b[0]-a[0]*b[2]; o[2]=a[0]*b[1]-a[1]*b[0]; }
static double norm3(const double a[3]){ return std::sqrt(dot(a,a)); }
static int eps3(int i,int j,int k){ if(i==j||j==k||i==k)return 0;
	if((i==0&&j==1&&k==2)||(i==1&&j==2&&k==0)||(i==2&&j==0&&k==1))return 1; return -1; }
static bool inv3(const double M[3][3], double Inv[3][3]){
	double d=M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])-M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])+M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]);
	if(std::fabs(d)<1e-300)return false; double id=1.0/d;
	Inv[0][0]=(M[1][1]*M[2][2]-M[1][2]*M[2][1])*id; Inv[0][1]=-(M[0][1]*M[2][2]-M[0][2]*M[2][1])*id; Inv[0][2]=(M[0][1]*M[1][2]-M[0][2]*M[1][1])*id;
	Inv[1][0]=-(M[1][0]*M[2][2]-M[1][2]*M[2][0])*id; Inv[1][1]=(M[0][0]*M[2][2]-M[0][2]*M[2][0])*id; Inv[1][2]=-(M[0][0]*M[1][2]-M[0][2]*M[1][0])*id;
	Inv[2][0]=(M[1][0]*M[2][1]-M[1][1]*M[2][0])*id; Inv[2][1]=-(M[0][0]*M[2][1]-M[0][1]*M[2][0])*id; Inv[2][2]=(M[0][0]*M[1][1]-M[0][1]*M[1][0])*id;
	return true; }
static void Jacobi3(double A[3][3], double eval[3], double evec[3][3]){
	double a[3][3]; for(int i=0;i<3;i++)for(int j=0;j<3;j++)a[i][j]=A[i][j]; double v[3][3]={{1,0,0},{0,1,0},{0,0,1}};
	for(int s=0;s<100;s++){ double off=std::fabs(a[0][1])+std::fabs(a[0][2])+std::fabs(a[1][2]); if(off<1e-18)break;
		for(int p=0;p<3;p++)for(int q=p+1;q<3;q++){ if(std::fabs(a[p][q])<1e-300)continue;
			double th=(a[q][q]-a[p][p])/(2*a[p][q]); double t=(th>=0?1.0:-1.0)/(std::fabs(th)+std::sqrt(th*th+1)); double c=1/std::sqrt(t*t+1),sn=t*c;
			for(int k=0;k<3;k++){double x=a[k][p],y=a[k][q];a[k][p]=c*x-sn*y;a[k][q]=sn*x+c*y;}
			for(int k=0;k<3;k++){double x=a[p][k],y=a[q][k];a[p][k]=c*x-sn*y;a[q][k]=sn*x+c*y;}
			for(int k=0;k<3;k++){double x=v[k][p],y=v[k][q];v[k][p]=c*x-sn*y;v[k][q]=sn*x+c*y;} } }
	for(int i=0;i<3;i++){eval[i]=a[i][i];for(int k=0;k<3;k++)evec[i][k]=v[k][i];} }

/* ----- shell geometry tensors at a sample point (xi3 station) ----- */
struct ShellGeom {
	double n[3], B1[3][3], B2[3][3], B1m[2][3][3], B2m[2][3][3];
	double F[3][3], Finv[3][3]; double h, xi3;
	/* second derivs of x3D needed by B,xi (at this xi3): x3D,xi_a,xi_l for a,l in {1,2} and the xi3 col */
	double x3D_ab[2][2][3]; /* x3D,xi(a),xi(b)  a,b in {0,1} */
	double n_l[2][3];        /* n,xi_l */
};
static bool buildGeom(const double x1[3],const double x2[3],const double x11[3],const double x22[3],const double x12[3],
                      double h,double xi3,ShellGeom& g)
{
	g.h=h; g.xi3=xi3;
	double cr[3]; cross(x1,x2,cr); double J=norm3(cr); if(J<1e-300)return false;
	for(int d=0;d<3;d++) g.n[d]=cr[d]/J;
	double A[3][3]; for(int i=0;i<3;i++)for(int j=0;j<3;j++) A[i][j]=((i==j?1.0:0.0)-g.n[i]*g.n[j])/J;
	double dcr[2][3]; { double t1[3],t2[3];
		cross(x11,x2,t1); cross(x1,x12,t2); for(int d=0;d<3;d++)dcr[0][d]=t1[d]+t2[d];
		cross(x12,x2,t1); cross(x1,x22,t2); for(int d=0;d<3;d++)dcr[1][d]=t1[d]+t2[d]; }
	double Jm[2];
	for(int m=0;m<2;m++){ Jm[m]=dot(g.n,dcr[m]);
		for(int i=0;i<3;i++){ double s=0; for(int k=0;k<3;k++) s+=((i==k?1.0:0.0)-g.n[i]*g.n[k])*dcr[m][k]; g.n_l[m][i]=s/J; } }
	double Am[2][3][3];
	for(int m=0;m<2;m++)for(int i=0;i<3;i++)for(int j=0;j<3;j++)
		Am[m][i][j]=-(g.n_l[m][i]*g.n[j]+g.n[i]*g.n_l[m][j])/J - A[i][j]*Jm[m]/J;
	for(int i=0;i<3;i++)for(int k=0;k<3;k++){ double s1=0,s2=0;
		for(int p=0;p<3;p++)for(int l=0;l<3;l++){int e=eps3(p,k,l); if(e)s1+=A[i][p]*e*x2[l];}
		for(int p=0;p<3;p++)for(int kk=0;kk<3;kk++){int e=eps3(p,kk,k); if(e)s2+=A[i][p]*e*x1[kk];}
		g.B1[i][k]=s1; g.B2[i][k]=s2; }
	const double* x2m[2]={x12,x22}; const double* x1m[2]={x11,x12};
	for(int m=0;m<2;m++)for(int i=0;i<3;i++)for(int k=0;k<3;k++){ double s1=0,s2=0;
		for(int p=0;p<3;p++)for(int l=0;l<3;l++){int e=eps3(p,k,l); if(e)s1+=e*(Am[m][i][p]*x2[l]+A[i][p]*x2m[m][l]);}
		for(int p=0;p<3;p++)for(int kk=0;kk<3;kk++){int e=eps3(p,kk,k); if(e)s2+=e*(Am[m][i][p]*x1[kk]+A[i][p]*x1m[m][kk]);}
		g.B1m[m][i][k]=s1; g.B2m[m][i][k]=s2; }
	/* x3D,xi_a (a=0,1) and the xi3 column of F */
	const double* xmm1[2]={x11,x12}; const double* xmm2[2]={x12,x22}; const double* x_a[2]={x1,x2};
	for(int a=0;a<2;a++)for(int i=0;i<3;i++){ double t=0; for(int k=0;k<3;k++) t+=g.B1[i][k]*xmm1[a][k]+g.B2[i][k]*xmm2[a][k];
		g.F[i][a]=x_a[a][i]+(h/2.0)*xi3*t; }
	for(int i=0;i<3;i++) g.F[i][2]=(h/2.0)*g.n[i];
	if(!inv3(g.F,g.Finv)) return false;
	/* x3D,xi(a),xi(b) at this xi3 (a,b in {0,1}); reference 2nd derivs are x11/x12/x22 (flat-ref ~0
	   on a flat plate, nonzero on cylinder). x3D,ab = x2D,ab + (h/2)xi3(...) ; we keep leading term. */
	const double* xab[2][2]={{x11,x12},{x12,x22}};
	for(int a=0;a<2;a++)for(int b=0;b<2;b++)for(int d=0;d<3;d++) g.x3D_ab[a][b][d]=xab[a][b][d];
	return true;
}

/* B_Iijk (Eq. 51) for a node with shape derivs P1,P2,P11,P12,P22 (Psi_I,xi...) */
static void Bmatrix(const ShellGeom& g,double P1,double P2,double P11,double P12,double P22,
                    double B[3][3][3])
{
	double hx=(g.h/2.0)*g.xi3;
	for(int i=0;i<3;i++)for(int k=0;k<3;k++){
		double br0 = (i==k?1.0:0.0)*P1 + hx*(g.B1m[0][i][k]*P1 + g.B1[i][k]*P11 + g.B2m[0][i][k]*P2 + g.B2[i][k]*P12);
		double br1 = (i==k?1.0:0.0)*P2 + hx*(g.B1m[1][i][k]*P1 + g.B1[i][k]*P12 + g.B2m[1][i][k]*P2 + g.B2[i][k]*P22);
		double br2 = (g.h/2.0)*(g.B1[i][k]*P1 + g.B2[i][k]*P2);
		for(int j=0;j<3;j++) B[i][j][k]= br0*g.Finv[0][j] + br1*g.Finv[1][j] + br2*g.Finv[2][j];
	}
}

/* B_Iijkl at xi3=0 (Eq. 54 with the xi3*S 3rd-derivative term dropped, since stabilization
   uses one-point through-thickness quadrature at xi3=0).  l in {0,1}.
   Needs B_Iimk at xi3=0 (pass Bz = Bmatrix at xi3=0). */
static void BmatrixGrad(const ShellGeom& g0,double P1,double P2,double P11,double P12,double P22,
                        double P1l[2],double P2l[2],/*Psi,xi_a,xi_l*/ const double Bz[3][3][3],
                        double Bg[3][3][3][2])
{
	/* g0 must be the xi3=0 geometry. x3D,xi_a,xi_l|0 = x2D 2nd derivs (g0.x3D_ab);
	   x3D,xi3,xi_l|0 = (h/2) n,xi_l. */
	for(int l=0;l<2;l++){
		/* Psi_I,xi_a,xi_l : a in {0,1}; build from P11,P12,P22 (l,a both in-plane) */
		double Pa_l[2]; Pa_l[0]=P1l[l]; Pa_l[1]=P2l[l]; /* Psi,xi(a),xi(l) */
		for(int i=0;i<3;i++)for(int k=0;k<3;k++){
			/* term a=0 (xi1) and a=1 (xi2) */
			double br[3];
			for(int a=0;a<2;a++){
				double Bdx=0; for(int m=0;m<3;m++) Bdx += Bz[i][m][k]*g0.x3D_ab[a][l][m]; /* B_Iimk x3D_m,xi_a,xi_l */
				br[a] = (i==k?1.0:0.0)*Pa_l[a] - Bdx;
			}
			/* xi3 column term: (h/2)(B1 Psi,1l + B1,l Psi,1 + B2 Psi,2l + B2,l Psi,2) - B_Iimk x3D_m,xi3,xi_l */
			double Bdx3=0; for(int m=0;m<3;m++) Bdx3 += Bz[i][m][k]*(g0.h/2.0)*g0.n_l[l][m];
			double t3 = (g0.h/2.0)*( g0.B1[i][k]*P1l[l] + g0.B1m[l][i][k]*P1 + g0.B2[i][k]*P2l[l] + g0.B2m[l][i][k]*P2 ) - Bdx3;
			for(int j=0;j<3;j++) Bg[i][j][k][l] = br[0]*g0.Finv[0][j] + br[1]*g0.Finv[1][j] + t3*g0.Finv[2][j];
		}
	}
}

/* independent reference: spatial velocity gradient from velocity parametric derivs (= #63 chain) */
static bool KLvelGrad(const double x1[3],const double x2[3],const double x11[3],const double x22[3],const double x12[3],
                      const double v1[3],const double v2[3],const double v11[3],const double v22[3],const double v12[3],
                      double h,double xi3,double L[3][3])
{
	ShellGeom g; if(!buildGeom(x1,x2,x11,x22,x12,h,xi3,g)) return false;
	double ndot[3]; for(int i=0;i<3;i++){double s=0;for(int k=0;k<3;k++)s+=g.B1[i][k]*v1[k]+g.B2[i][k]*v2[k];ndot[i]=s;}
	const double* vmm1[2]={v11,v12}; const double* vmm2[2]={v12,v22}; const double* v_m[2]={v1,v2};
	double dvxi[3][3];
	for(int m=0;m<2;m++)for(int i=0;i<3;i++){ double s=v_m[m][i]; double t=0;
		for(int k=0;k<3;k++) t+=g.B1m[m][i][k]*v1[k]+g.B1[i][k]*vmm1[m][k]+g.B2m[m][i][k]*v2[k]+g.B2[i][k]*vmm2[m][k];
		dvxi[i][m]=s+(h/2.0)*xi3*t; }
	for(int i=0;i<3;i++) dvxi[i][2]=(h/2.0)*ndot[i];
	for(int i=0;i<3;i++)for(int k=0;k<3;k++){double s=0;for(int l=0;l<3;l++)s+=dvxi[i][l]*g.Finv[l][k];L[i][k]=s;}
	return true;
}

/* ---- driver ---- */
struct Cloud { std::vector<double> X,Y,Z; };
static void cylinder(Cloud& c,int nt,int nz,double Rc,double Zlen,double& h){
	double dth=2*M_PI/nt, dz=Zlen/(nz-1); h=std::max(Rc*dth,dz);
	for(int i=0;i<nt;i++)for(int j=0;j<nz;j++){double th=i*dth,z=j*dz; c.X.push_back(Rc*std::cos(th));c.Y.push_back(Rc*std::sin(th));c.Z.push_back(z);} }
static void plate(Cloud& c,int n,double L,double& h){ h=L/(n-1);
	for(int i=0;i<n;i++)for(int j=0;j<n;j++){c.X.push_back(i*h);c.Y.push_back(j*h);c.Z.push_back(0.0);} }

int run();
int main(){ try{return run();}
	catch(ExceptionT::CodeT& e){cout<<"\n*** ExceptionT "<<int(e)<<" ("<<ExceptionT::ToString(e)<<") ***\n";return 3;}
	catch(...){cout<<"\n*** unknown exception ***\n";return 4;} }

/* run B-matrix tests on a given cloud at the node nearest 'target'; returns max errors */
static int testCloud(const char* tag,const Cloud& C,double tx,double ty,double tz,double thick,
                     bool doGradFD)
{
	int M=C.X.size();
	int P=-1; double bb=1e9; for(int i=0;i<M;i++){double d=std::fabs(C.X[i]-tx)+std::fabs(C.Y[i]-ty)+std::fabs(C.Z[i]-tz); if(d<bb){bb=d;P=i;}}
	/* spacing estimate */
	double hC=0; { double best=1e9; for(int q=0;q<M;q++) if(q!=P){double dx=C.X[q]-C.X[P],dy=C.Y[q]-C.Y[P],dz=C.Z[q]-C.Z[P];double d=std::sqrt(dx*dx+dy*dy+dz*dz); if(d<best)best=d;} hC=best; }
	double R=3.5*hC; std::vector<int> nb;
	for(int q=0;q<M;q++){double dx=C.X[q]-C.X[P],dy=C.Y[q]-C.Y[P],dz=C.Z[q]-C.Z[P]; if(std::sqrt(dx*dx+dy*dy+dz*dz)<0.99*R)nb.push_back(q);}
	int nn=nb.size();
	double mean[3]={0,0,0}; for(int k=0;k<nn;k++){mean[0]+=C.X[nb[k]];mean[1]+=C.Y[nb[k]];mean[2]+=C.Z[nb[k]];} for(int d=0;d<3;d++)mean[d]/=nn;
	double Cov[3][3]={{0,0,0},{0,0,0},{0,0,0}};
	for(int k=0;k<nn;k++){double dd[3]={C.X[nb[k]]-mean[0],C.Y[nb[k]]-mean[1],C.Z[nb[k]]-mean[2]}; for(int a=0;a<3;a++)for(int b=0;b<3;b++)Cov[a][b]+=dd[a]*dd[b];}
	double ev[3],evec[3][3]; Jacobi3(Cov,ev,evec);
	int o[3]={0,1,2}; for(int a=0;a<3;a++)for(int b=a+1;b<3;b++)if(ev[o[b]]>ev[o[a]]){int t=o[a];o[a]=o[b];o[b]=t;}
	double q1[3]={evec[o[0]][0],evec[o[0]][1],evec[o[0]][2]}, q2[3]={evec[o[1]][0],evec[o[1]][1],evec[o[1]][2]};
	dArray2DT lc(nn,2); for(int k=0;k<nn;k++){double dxv[3]={C.X[nb[k]]-C.X[P],C.Y[nb[k]]-C.Y[P],C.Z[nb[k]]-C.Z[P]}; lc(k,0)=dot(dxv,q1);lc(k,1)=dot(dxv,q2);}

	/* arbitrary smooth velocity field at nodes */
	std::vector<double> VX(nn),VY(nn),VZ(nn);
	for(int k=0;k<nn;k++){double Xx[3]={C.X[nb[k]],C.Y[nb[k]],C.Z[nb[k]]}; VX[k]=0.3*Xx[1]*Xx[2];VY[k]=0.2*Xx[0];VZ[k]=0.1*Xx[0]*Xx[0];}

	D2OrthoMLS2DT mls(2); mls.Initialize(); dArrayT dmax(nn); dmax=R;

	/* reconstruct position+velocity parametric derivs at a sample (local coords) */
	auto recon=[&](double s0,double s1,double x1[3],double x2[3],double x11[3],double x22[3],double x12[3],
	               double v1[3],double v2[3],double v11[3],double v22[3],double v12[3],
	               std::vector<double>* phi,std::vector<double>* d0,std::vector<double>* d1,
	               std::vector<double>* dd0,std::vector<double>* dd1,std::vector<double>* dd2)->bool{
		dArrayT fp(2); fp[0]=s0; fp[1]=s1; if(!mls.SetField(lc,dmax,fp)) return false;
		const dArrayT& ph=mls.phi(); const dArray2DT& Dp=mls.Dphi(); const dArray2DT& DDp=mls.DDphi();
		for(int d=0;d<3;d++){x1[d]=x2[d]=x11[d]=x22[d]=x12[d]=0;v1[d]=v2[d]=v11[d]=v22[d]=v12[d]=0;}
		for(int I=0;I<nn;I++){ double Xq[3]={C.X[nb[I]],C.Y[nb[I]],C.Z[nb[I]]}, Vq[3]={VX[I],VY[I],VZ[I]};
			for(int d=0;d<3;d++){ x1[d]+=Dp(0,I)*Xq[d];x2[d]+=Dp(1,I)*Xq[d];x11[d]+=DDp(0,I)*Xq[d];x22[d]+=DDp(1,I)*Xq[d];x12[d]+=DDp(2,I)*Xq[d];
				v1[d]+=Dp(0,I)*Vq[d];v2[d]+=Dp(1,I)*Vq[d];v11[d]+=DDp(0,I)*Vq[d];v22[d]+=DDp(1,I)*Vq[d];v12[d]+=DDp(2,I)*Vq[d]; }
			if(phi){(*phi)[I]=ph[I];(*d0)[I]=Dp(0,I);(*d1)[I]=Dp(1,I);(*dd0)[I]=DDp(0,I);(*dd1)[I]=DDp(1,I);(*dd2)[I]=DDp(2,I);} }
		return true; };

	std::vector<double> ph(nn),d0(nn),d1(nn),dd0(nn),dd1(nn),dd2(nn);
	double x1[3],x2[3],x11[3],x22[3],x12[3],v1[3],v2[3],v11[3],v22[3],v12[3];
	recon(0,0,x1,x2,x11,x22,x12,v1,v2,v11,v22,v12,&ph,&d0,&d1,&dd0,&dd1,&dd2);

	int fails=0;
	const int voigt[6][2]={{0,0},{1,1},{2,2},{1,2},{0,2},{0,1}};

	/* ---- Test A & B: B-matrix reproduces grad(v3D) and Voigt strain rate, at xi3=0 and 0.5 ---- */
	double xi3s[2]={0.0,0.5};
	for(int t=0;t<2;t++){ double xi3=xi3s[t];
		ShellGeom g; buildGeom(x1,x2,x11,x22,x12,thick,xi3,g);
		double Lref[3][3]; KLvelGrad(x1,x2,x11,x22,x12,v1,v2,v11,v22,v12,thick,xi3,Lref);
		double Lb[3][3]={{0,0,0},{0,0,0},{0,0,0}};
		double Vstrain[6]={0,0,0,0,0,0};
		for(int I=0;I<nn;I++){ double B[3][3][3]; Bmatrix(g,d0[I],d1[I],dd0[I],dd2[I],dd1[I],B);
			/* note DDphi rows: 0=xi1xi1,1=xi2xi2,2=xi1xi2 -> P11=dd0,P22=dd1,P12=dd2 */
			double vI[3]={VX[I],VY[I],VZ[I]};
			for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int k=0;k<3;k++) Lb[i][j]+=B[i][j][k]*vI[k];
			for(int s=0;s<6;s++){int a=voigt[s][0],b=voigt[s][1]; double Bv= (a==b)? B[a][b][0]*0:0; (void)Bv;
				double row[3]; for(int k=0;k<3;k++) row[k]= (a==b)? B[a][b][k] : (B[a][b][k]+B[b][a][k]);
				for(int k=0;k<3;k++) Vstrain[s]+=row[k]*vI[k]; }
		}
		double eL=0; for(int i=0;i<3;i++)for(int j=0;j<3;j++) eL=std::max(eL,std::fabs(Lb[i][j]-Lref[i][j]));
		/* Voigt strain rate reference = sym(Lref) */
		double eV=0; for(int s=0;s<6;s++){int a=voigt[s][0],b=voigt[s][1]; double Dref=0.5*(Lref[a][b]+Lref[b][a]); double exp=(a==b)?Dref:2*Dref; eV=std::max(eV,std::fabs(Vstrain[s]-exp)); }
		bool okL=eL<1e-9, okV=eV<1e-9; fails+=!okL; fails+=!okV;
		cout<<"  ["<<tag<<" xi3="<<std::fixed<<setprecision(1)<<xi3<<"] sum_I B_Iijk v_Ik vs grad(v3D): "<<scientific<<setprecision(2)<<eL
		    <<(okL?" PASS":" FAIL")<<" | Voigt strain-rate err: "<<eV<<(okV?" PASS":" FAIL")<<"\n";
	}

	/* ---- Test C: parametric gradient B_Iijkl at xi3=0 vs FD of grad(v3D) (curved only) ---- */
	if(doGradFD){
		ShellGeom g0; buildGeom(x1,x2,x11,x22,x12,thick,0.0,g0);
		/* B_Iimk at xi3=0 per node, plus Psi second derivs for B,xi */
		/* analytic sum_I B_Iijkl v_Ik */
		double Bg_sum[3][3][2]={{{0,0},{0,0},{0,0}},{{0,0},{0,0},{0,0}},{{0,0},{0,0},{0,0}}};
		for(int I=0;I<nn;I++){
			double Bz[3][3][3]; Bmatrix(g0,d0[I],d1[I],dd0[I],dd2[I],dd1[I],Bz); /* xi3=0 */
			/* Psi,xi_a,xi_l : in-plane second derivs. P1l[l]=Psi,xi1,xil ; P2l[l]=Psi,xi2,xil
			   l=0(xi1): P1l=dd0(=,11), P2l=dd2(=,12); l=1(xi2): P1l=dd2(=,12), P2l=dd1(=,22) */
			double P1l[2]={dd0[I],dd2[I]}, P2l[2]={dd2[I],dd1[I]};
			double Bg[3][3][3][2]; BmatrixGrad(g0,d0[I],d1[I],dd0[I],dd2[I],dd1[I],P1l,P2l,Bz,Bg);
			double vI[3]={VX[I],VY[I],VZ[I]};
			for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int l=0;l<2;l++)for(int k=0;k<3;k++) Bg_sum[i][j][l]+=Bg[i][j][k][l]*vI[k];
		}
		/* FD of grad(v3D) wrt xi_l at xi3=0 */
		double delta=1e-4; double fd[3][3][2];
		for(int l=0;l<2;l++){ double sp[2]={0,0}, sm[2]={0,0}; sp[l]=delta; sm[l]=-delta;
			double a1[3],a2[3],a11[3],a22[3],a12[3],b1[3],b2[3],b11[3],b22[3],b12[3]; double Lp[3][3],Lm[3][3];
			recon(sp[0],sp[1],a1,a2,a11,a22,a12,b1,b2,b11,b22,b12,0,0,0,0,0,0); KLvelGrad(a1,a2,a11,a22,a12,b1,b2,b11,b22,b12,thick,0.0,Lp);
			recon(sm[0],sm[1],a1,a2,a11,a22,a12,b1,b2,b11,b22,b12,0,0,0,0,0,0); KLvelGrad(a1,a2,a11,a22,a12,b1,b2,b11,b22,b12,thick,0.0,Lm);
			for(int i=0;i<3;i++)for(int j=0;j<3;j++) fd[i][j][l]=(Lp[i][j]-Lm[i][j])/(2*delta);
		}
		double eG=0,sc=0; for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int l=0;l<2;l++){ eG=std::max(eG,std::fabs(Bg_sum[i][j][l]-fd[i][j][l])); sc=std::max(sc,std::fabs(fd[i][j][l])); }
		double rel=eG/std::max(sc,1e-12); bool ok=rel<1e-5; fails+=!ok;
		cout<<"  ["<<tag<<" xi3=0] sum_I B_Iijkl v_Ik vs FD d(grad v3D)/dxi: "<<scientific<<setprecision(2)<<eG<<" (rel "<<rel<<")"<<(ok?" PASS":" FAIL")<<"\n";
	}
	return fails;
}

int run()
{
	cout<<"=== KL-shell strain-displacement (B-matrix) de-risk (issue #64) ===\n\n";
	int fails=0;
	{ Cloud c; double h; plate(c,15,1.0,h); fails+=testCloud("plate",c,0.5,0.5,0.0,0.05,false); }
	{ Cloud c; double h; cylinder(c,64,21,1.0,2.0,h); fails+=testCloud("cyl",c,1.0,0.0,1.0,0.05,true); }

	cout<<"\n=== VERDICT (issue #64) ===\n";
	cout<<"  Shell B-matrix B_Iijk (Eq.51) reproduces grad(v3D) and the Voigt strain rate;\n"
	      "  parametric gradient B_Iijkl (Eq.54, xi3=0) matches FD of grad(v3D) on a curved\n"
	      "  surface (the naturally-stabilized nodal-integration operator). status: "<<(fails==0?"PASS":"FAIL")<<"\n";
	cout<<"  (full stabilized internal force, Eq.57, validated once constitutive lands -- #65)\n";
	return fails==0?0:1;
}
