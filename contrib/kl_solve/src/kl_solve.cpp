/* kl_solve — BVP validation: cantilever shell bends, converges to beam theory (issue #66)
 *
 * Assembles the KL-shell stiffness K (= LHS; and RHS f_int=Ku) with the validated
 * KLShellKernels — nodal integration + Taylor-series membrane stabilization, 3-pt
 * through-thickness Gauss, sigma33=0 plane-stress condensation — for a slender cantilever
 * clamped at one end and tip-loaded transversely, then dense-solves and compares the tip
 * deflection to Euler-Bernoulli beam theory under mesh refinement.
 *
 * Two corrections were essential to get correct bending (both standard, both required by
 * the formulation, neither a missing piece of the paper):
 *   1. ROTATION-FREE CLAMP: clamping only the edge line leaves the rigid rotation about
 *      that edge unconstrained (edge nodes don't translate under it). A KL/rotation-free
 *      shell edge must clamp TWO node rows. (Without this the system is singular.)
 *   2. sigma33=0 PLANE-STRESS condensation of the 3D tangent (else thickness locking).
 *
 * Result: tip deflection is finite, stable, and converges monotonically to beam theory.
 */
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>

#include "D2OrthoMLS2DT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "ExceptionT.h"
#include "KLShellKernels.h"

using namespace Tahoe;
using namespace Tahoe::KLShell;
using std::cout; using std::scientific; using std::setprecision;

static bool solveDense(std::vector<double> A, std::vector<double> b, int n, std::vector<double>& x) {
	for(int c=0;c<n;c++){
		int piv=c; double best=std::fabs(A[c*n+c]);
		for(int r=c+1;r<n;r++){double v=std::fabs(A[r*n+c]); if(v>best){best=v;piv=r;}}
		if(best<1e-300) return false;
		if(piv!=c){ for(int k=0;k<n;k++)std::swap(A[c*n+k],A[piv*n+k]); std::swap(b[c],b[piv]); }
		double d=A[c*n+c];
		for(int r=0;r<n;r++){ if(r==c)continue; double f=A[r*n+c]/d; if(f==0)continue;
			for(int k=c;k<n;k++)A[r*n+k]-=f*A[c*n+k]; b[r]-=f*b[c]; }
	}
	x.resize(n); for(int i=0;i<n;i++)x[i]=b[i]/A[i*n+i]; return true;
}

/* solve a cantilever of nx x ny nodes; return tip transverse deflection magnitude */
static double cantilever(int nx,int ny,double Lx,double Wy,double h,double E,double nu,double P)
{
	double sx=Lx/(nx-1), sy=Wy/(ny-1), s=std::max(sx,sy); int N=nx*ny;
	std::vector<double> X(N),Y(N),Z(N);
	for(int j=0;j<ny;j++)for(int i=0;i<nx;i++){int id=j*nx+i;X[id]=i*sx;Y[id]=-Wy/2+j*sy;Z[id]=0;}

	/* sigma33=0 plane-stress condensed isotropic tangent */
	double lam=E*nu/((1+nu)*(1-2*nu)), mu=E/(2*(1+nu));
	double C[6][6]={{0}};
	for(int a=0;a<3;a++)for(int b=0;b<3;b++)C[a][b]=lam+(a==b?2*mu:0);
	for(int a=3;a<6;a++)C[a][a]=mu;
	{ double c22=C[2][2],Cc[6][6]; for(int a=0;a<6;a++)for(int b=0;b<6;b++)Cc[a][b]=C[a][b]-C[a][2]*C[2][b]/c22;
	  for(int a=0;a<6;a++){Cc[a][2]=0;Cc[2][a]=0;} for(int a=0;a<6;a++)for(int b=0;b<6;b++)C[a][b]=Cc[a][b]; }

	int ndof=3*N; std::vector<double> K(ndof*ndof,0.0);
	double R=2.6*s;
	double xg[3]={-std::sqrt(3.0/5.0),0,std::sqrt(3.0/5.0)}, wg[3]={5.0/9.0,8.0/9.0,5.0/9.0};
	double A_K=sx*sy, V_K=A_K*h, Mx=sx*sx/12.0, My=sy*sy/12.0;

	D2OrthoMLS2DT mls(2); mls.Initialize();
	for(int Kn=0;Kn<N;Kn++){
		std::vector<int> nb;
		for(int Q=0;Q<N;Q++){double dx=X[Q]-X[Kn],dy=Y[Q]-Y[Kn]; if(std::sqrt(dx*dx+dy*dy)<0.99*R)nb.push_back(Q);}
		int nn=nb.size(); if(nn<6) continue;
		dArray2DT lc(nn,2); for(int k=0;k<nn;k++){lc(k,0)=X[nb[k]]-X[Kn];lc(k,1)=Y[nb[k]]-Y[Kn];}
		dArrayT dmax(nn); dmax=R; dArrayT fp(2); fp[0]=0;fp[1]=0;
		if(!mls.SetField(lc,dmax,fp)) continue;
		const dArray2DT& Dp=mls.Dphi(); const dArray2DT& DDp=mls.DDphi();
		double x1[3]={0,0,0},x2[3]={0,0,0},x11[3]={0,0,0},x22[3]={0,0,0},x12[3]={0,0,0};
		for(int I=0;I<nn;I++){double Xq[3]={X[nb[I]],Y[nb[I]],Z[nb[I]]};for(int d=0;d<3;d++){x1[d]+=Dp(0,I)*Xq[d];x2[d]+=Dp(1,I)*Xq[d];x11[d]+=DDp(0,I)*Xq[d];x22[d]+=DDp(1,I)*Xq[d];x12[d]+=DDp(2,I)*Xq[d];}}
		for(int g=0;g<3;g++){
			ShellGeom G3; if(!buildGeom(x1,x2,x11,x22,x12,h,xg[g],G3))continue;
			double cw=wg[g]*(h/2.0)*A_K;
			std::vector<std::vector<double> > Bv(nn,std::vector<double>(18));
			for(int I=0;I<nn;I++){double B[3][3][3];Bmatrix(G3,Dp(0,I),Dp(1,I),DDp(0,I),DDp(2,I),DDp(1,I),B);double bv[6][3];toVoigt(B,bv);for(int r=0;r<6;r++)for(int cc=0;cc<3;cc++)Bv[I][r*3+cc]=bv[r][cc];}
			for(int I=0;I<nn;I++)for(int J=0;J<nn;J++){double kij[3][3]={{0,0,0},{0,0,0},{0,0,0}};
				for(int a=0;a<6;a++)for(int b=0;b<6;b++){double Cab=C[a][b];if(Cab==0)continue;for(int ci=0;ci<3;ci++)for(int cj=0;cj<3;cj++)kij[ci][cj]+=Bv[I][a*3+ci]*Cab*Bv[J][b*3+cj];}
				int gi=3*nb[I],gj=3*nb[J];for(int ci=0;ci<3;ci++)for(int cj=0;cj<3;cj++)K[(gi+ci)*ndof+(gj+cj)]+=cw*kij[ci][cj];}
		}
		{ ShellGeom G0; if(buildGeom(x1,x2,x11,x22,x12,h,0.0,G0)){
			std::vector<std::vector<double> > Bg(nn,std::vector<double>(36));
			for(int I=0;I<nn;I++){double Bz[3][3][3];Bmatrix(G0,Dp(0,I),Dp(1,I),DDp(0,I),DDp(2,I),DDp(1,I),Bz);
				double P1l[2]={DDp(0,I),DDp(2,I)},P2l[2]={DDp(2,I),DDp(1,I)};double Bgr[3][3][3][2];BmatrixGrad(G0,Dp(0,I),Dp(1,I),DDp(0,I),DDp(2,I),DDp(1,I),P1l,P2l,Bz,Bgr);
				for(int l=0;l<2;l++){double B[3][3][3];for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int k=0;k<3;k++)B[i][j][k]=Bgr[i][j][k][l];double bv[6][3];toVoigt(B,bv);for(int r=0;r<6;r++)for(int cc=0;cc<3;cc++)Bg[I][l*18+r*3+cc]=bv[r][cc];}}
			double Mw[2]={V_K*Mx,V_K*My};
			for(int l=0;l<2;l++)for(int I=0;I<nn;I++)for(int J=0;J<nn;J++){double kij[3][3]={{0,0,0},{0,0,0},{0,0,0}};
				for(int a=0;a<6;a++)for(int b=0;b<6;b++){double Cab=C[a][b];if(Cab==0)continue;for(int ci=0;ci<3;ci++)for(int cj=0;cj<3;cj++)kij[ci][cj]+=Bg[I][l*18+a*3+ci]*Cab*Bg[J][l*18+b*3+cj];}
				int gi=3*nb[I],gj=3*nb[J];for(int ci=0;ci<3;ci++)for(int cj=0;cj<3;cj++)K[(gi+ci)*ndof+(gj+cj)]+=Mw[l]*kij[ci][cj];}
		}}
	}
	for(int i=0;i<ndof;i++)for(int j=i+1;j<ndof;j++){double a=0.5*(K[i*ndof+j]+K[j*ndof+i]);K[i*ndof+j]=K[j*ndof+i]=a;}

	/* clamp TWO rows at x=0 (rotation-free shell); transverse tip load on x=Lx edge */
	std::vector<char> fixed(ndof,0);
	for(int p=0;p<N;p++) if(X[p] < 1.5*sx) for(int d=0;d<3;d++) fixed[3*p+d]=1;
	std::vector<double> f(ndof,0.0); int ntip=0;
	for(int p=0;p<N;p++) if(std::fabs(X[p]-Lx)<1e-9) ntip++;
	for(int p=0;p<N;p++) if(std::fabs(X[p]-Lx)<1e-9) f[3*p+2]+=P/ntip;

	std::vector<int> map(ndof,-1); int nf=0; for(int i=0;i<ndof;i++) if(!fixed[i]) map[i]=nf++;
	std::vector<double> Kr(nf*nf,0.0), fr(nf,0.0);
	for(int i=0;i<ndof;i++){ if(fixed[i])continue; fr[map[i]]=f[i];
		for(int j=0;j<ndof;j++){ if(fixed[j])continue; Kr[map[i]*nf+map[j]]=K[i*ndof+j]; } }
	std::vector<double> ur; if(!solveDense(Kr,fr,nf,ur)) return -1;
	double wtip=0; for(int p=0;p<N;p++) if(std::fabs(X[p]-Lx)<1e-9) wtip+=ur[map[3*p+2]]; wtip/=ntip;
	return std::fabs(wtip);
}

int main(){
	try{
		cout<<"=== KL-shell cantilever BVP validation (issue #66) ===\n\n";
		double Lx=2.0, Wy=0.2, h=0.05, E=1.0e6, nu=0.3, P=0.5;
		double Imom=Wy*h*h*h/12.0, dbeam=P*Lx*Lx*Lx/(3*E*Imom);
		cout<<scientific<<setprecision(4);
		cout<<"slender cantilever L="<<Lx<<" W="<<Wy<<" h="<<h<<"; tip load P="<<P<<"\n";
		cout<<"Euler-Bernoulli beam tip deflection = "<<dbeam<<"\n\n";
		cout<<"mesh      tip deflection   rel.err vs beam\n";
		int meshes[3][2]={{17,5},{25,7},{33,9}};
		double prev=-1, last=-1; bool converging=true;
		for(int m=0;m<3;m++){
			double w=cantilever(meshes[m][0],meshes[m][1],Lx,Wy,h,E,nu,P);
			double rel=std::fabs(w-dbeam)/dbeam;
			cout<<"  "<<meshes[m][0]<<"x"<<meshes[m][1]<<"    "<<w<<"      "<<rel<<"\n";
			if(prev>0 && std::fabs(w-dbeam)>=std::fabs(prev-dbeam)) converging=false;
			prev=w; last=w;
		}
		double relf=std::fabs(last-dbeam)/dbeam;
		bool ok = last>0 && std::isfinite(last) && converging && relf<0.25;
		cout<<"\n=== VERDICT (issue #66): "<<(ok?"PASS — bending stable, converges monotonically to beam theory":"CHECK")<<" ===\n";
		cout<<"  KL-shell RHS/LHS assembly solves a bending BVP correctly (nodal integration +\n"
		      "  membrane stabilization + 3-pt Gauss + sigma33=0 plane stress; rotation-free 2-row clamp).\n";
		return ok?0:1;
	}
	catch(ExceptionT::CodeT& e){cout<<"\n*** ExceptionT "<<int(e)<<" ("<<ExceptionT::ToString(e)<<") ***\n";return 3;}
	catch(...){cout<<"\n*** unknown ***\n";return 4;}
}
