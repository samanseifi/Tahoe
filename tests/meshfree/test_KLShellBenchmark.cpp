/* test_KLShellBenchmark.cpp — meshfree KL-shell obstacle-course benchmark.
 *
 * Reproduces the Scordelis-Lo roof from the Belytschko shell obstacle course
 * (Wang & Bazilevs 2025, section 4.2.1): a short cylindrical roof under gravity,
 * supported by rigid end diaphragms, with free longitudinal edges. The reference
 * vertical deflection at the mid-point of a free edge is 0.3006.
 *
 * Driven through the validated KL-shell kernels (KLShellKernels.h) with a general
 * curved-surface assembler: per-node PCA tangent frame, nodal integration + membrane
 * stabilization + 3-pt through-thickness Gauss, sigma33=0 plane-stress condensation,
 * dense linear solve. The benchmark checks the deflection converges toward 0.3006
 * under mesh refinement.
 */

#include "gtest/gtest.h"

#include <cmath>
#include <vector>
#include <algorithm>

#include "MLSSolverT.h"
#include "MeshFreeT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "KLShellKernels.h"

using namespace Tahoe;
using namespace Tahoe::KLShell;

namespace {

/* sigma33=0 plane-stress condensation of the isotropic 3D Voigt tangent */
void PlaneStressTangent(double E, double nu, double C[6][6])
{
	double lam = E*nu/((1+nu)*(1-2*nu)), mu = E/(2*(1+nu));
	for (int a=0;a<6;a++) for (int b=0;b<6;b++) C[a][b] = 0.0;
	for (int a=0;a<3;a++) for (int b=0;b<3;b++) C[a][b] = lam + (a==b ? 2*mu : 0.0);
	for (int a=3;a<6;a++) C[a][a] = mu;
	double c22 = C[2][2];
	double Cc[6][6];
	for (int a=0;a<6;a++) for (int b=0;b<6;b++) Cc[a][b] = C[a][b] - C[a][2]*C[2][b]/c22;
	for (int a=0;a<6;a++) { Cc[a][2] = 0.0; Cc[2][a] = 0.0; }
	for (int a=0;a<6;a++) for (int b=0;b<6;b++) C[a][b] = Cc[a][b];
}

bool SolveDense(std::vector<double> A, std::vector<double> b, int n, std::vector<double>& x)
{
	for (int c=0;c<n;c++) {
		int piv=c; double best=std::fabs(A[c*n+c]);
		for (int r=c+1;r<n;r++) { double v=std::fabs(A[r*n+c]); if (v>best) { best=v; piv=r; } }
		if (best < 1e-300) return false;
		if (piv!=c) { for (int k=0;k<n;k++) std::swap(A[c*n+k],A[piv*n+k]); std::swap(b[c],b[piv]); }
		double d=A[c*n+c];
		for (int r=0;r<n;r++) { if (r==c) continue; double f=A[r*n+c]/d; if (f==0) continue;
			for (int k=c;k<n;k++) A[r*n+k]-=f*A[c*n+k]; b[r]-=f*b[c]; }
	}
	x.resize(n);
	for (int i=0;i<n;i++) x[i] = b[i]/A[i*n+i];
	return true;
}

/* general curved-surface shell stiffness assembler.
 * X is a 3D point cloud (N x 3, row-major); nodalArea is the per-node integration weight;
 * h is the shell thickness; C is the (plane-stress) Voigt tangent. K is dense [3N x 3N]. */
void AssembleShell(const std::vector<double>& X, const std::vector<double>& nodalArea,
	double h, const double C[6][6], double support, std::vector<double>& K, double alpha = 1.0)
{
	int N = int(X.size())/3;
	int ndof = 3*N;
	K.assign((size_t)ndof*ndof, 0.0);

	double xg[3] = {-std::sqrt(3.0/5.0), 0.0, std::sqrt(3.0/5.0)};
	double wg[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0};

	/* RKPM (reproducing kernel), quadratic completeness, Gaussian window.
	 * Gaussian (not cubic-spline) because its 2nd-derivative path is robust at one-sided
	 * boundary neighborhoods — the cubic-spline window produces NaN DDphi there. */
	dArrayT gwin(3);
	gwin[0] = 1.5;  /* support scaling */
	gwin[1] = 0.4;  /* sharpening */
	gwin[2] = 3.0;  /* cutoff */
	MLSSolverT rkpm(2, 2, false, MeshFreeT::kGaussian, gwin);
	rkpm.Initialize();

	int nskip = 0;
	for (int P=0;P<N;P++) {

		/* neighbors within the support radius (3D distance on the surface) */
		std::vector<int> nb;
		for (int Q=0;Q<N;Q++) {
			double d0=X[3*Q]-X[3*P], d1=X[3*Q+1]-X[3*P+1], d2=X[3*Q+2]-X[3*P+2];
			if (std::sqrt(d0*d0+d1*d1+d2*d2) < 0.99*support) nb.push_back(Q);
		}
		int nn = int(nb.size());
		if (nn < 6) continue;

		/* PCA tangent frame at P */
		std::vector<double> nbX(3*nn);
		for (int k=0;k<nn;k++) { nbX[3*k]=X[3*nb[k]]; nbX[3*k+1]=X[3*nb[k]+1]; nbX[3*k+2]=X[3*nb[k]+2]; }
		double psi1[3], psi2[3], n0[3];
		PCAFrame(&nbX[0], nn, psi1, psi2, n0);

		/* local parametric coordinates of the neighbors */
		dArray2DT lc(nn, 2);
		for (int k=0;k<nn;k++) {
			double dxv[3] = {X[3*nb[k]]-X[3*P], X[3*nb[k]+1]-X[3*P+1], X[3*nb[k]+2]-X[3*P+2]};
			lc(k,0) = Dot(dxv, psi1);
			lc(k,1) = Dot(dxv, psi2);
		}
		dArray2DT np(nn, 1); np = support;          /* RKPM support size per neighbor */
		dArrayT vol(nn);
		for (int k=0;k<nn;k++) vol[k] = nodalArea[nb[k]]; /* nodal volumes for the moment matrix */
		dArrayT sample(2); sample[0]=0.0; sample[1]=0.0;
		if (!rkpm.SetField(lc, np, vol, sample, 3)) { nskip++; continue; } /* order 3 for the bending stabilization */
		const dArray2DT& Dp = rkpm.Dphi();
		const dArray2DT& DDp = rkpm.DDphi();
		const dArray2DT& DDDp = rkpm.DDDphi(); /* 2D comps: 0:111 1:122 2:112 3:222 */

		/* reference position parametric derivatives at P */
		double x1[3]={0,0,0},x2[3]={0,0,0},x11[3]={0,0,0},x22[3]={0,0,0},x12[3]={0,0,0};
		for (int I=0;I<nn;I++) {
			double Xq[3] = {X[3*nb[I]], X[3*nb[I]+1], X[3*nb[I]+2]};
			for (int d=0;d<3;d++) {
				x1[d]+=Dp(0,I)*Xq[d]; x2[d]+=Dp(1,I)*Xq[d];
				x11[d]+=DDp(0,I)*Xq[d]; x22[d]+=DDp(1,I)*Xq[d]; x12[d]+=DDp(2,I)*Xq[d];
			}
		}

		double A_K = nodalArea[P];
		double V_K = A_K*h;
		/* nodal cell length scale for the stabilization moment (~ sqrt(area)) */
		double cell = std::sqrt(A_K);
		double Mmom = cell*cell/12.0;

		/* nodal-integration term (3-pt Gauss through thickness) */
		for (int g=0;g<3;g++) {
			ShellGeom G;
			if (!BuildGeom(x1,x2,x11,x22,x12,h,xg[g],G)) continue;
			double e1[3],e2[3]; OrthoTangents(G.n, e1, e2); /* local frame, plane stress along normal */
			double cw = wg[g]*(h/2.0)*A_K;
			std::vector<std::vector<double> > Bv(nn, std::vector<double>(18));
			for (int I=0;I<nn;I++) {
				double B[3][3][3];
				BMatrix(G, Dp(0,I),Dp(1,I),DDp(0,I),DDp(2,I),DDp(1,I), B);
				double bv[6][3]; ToVoigtLocal(B, e1, e2, G.n, bv);
				for (int r=0;r<6;r++) for (int c=0;c<3;c++) Bv[I][r*3+c]=bv[r][c];
			}
			for (int I=0;I<nn;I++) for (int J=0;J<nn;J++) {
				double kij[3][3]={{0,0,0},{0,0,0},{0,0,0}};
				for (int a=0;a<6;a++) for (int b=0;b<6;b++) {
					double Cab=C[a][b]; if (Cab==0.0) continue;
					for (int ci=0;ci<3;ci++) for (int cj=0;cj<3;cj++) kij[ci][cj]+=Bv[I][a*3+ci]*Cab*Bv[J][b*3+cj];
				}
				int gi=3*nb[I], gj=3*nb[J];
				for (int ci=0;ci<3;ci++) for (int cj=0;cj<3;cj++) K[(size_t)(gi+ci)*ndof+(gj+cj)] += cw*kij[ci][cj];
			}
		}

		/* membrane stabilization (xi3=0) */
		ShellGeom G0;
		if (BuildGeom(x1,x2,x11,x22,x12,h,0.0,G0)) {
			double e1[3],e2[3]; OrthoTangents(G0.n, e1, e2);
			std::vector<std::vector<double> > Bg(nn, std::vector<double>(36));
			for (int I=0;I<nn;I++) {
				double Bz[3][3][3];
				BMatrix(G0, Dp(0,I),Dp(1,I),DDp(0,I),DDp(2,I),DDp(1,I), Bz);
				double P1l[2]={DDp(0,I),DDp(2,I)}, P2l[2]={DDp(2,I),DDp(1,I)};
				double Bgr[3][3][3][2];
				BMatrixGradient(G0, Dp(0,I),Dp(1,I), P1l,P2l,Bz,Bgr);
				for (int l=0;l<2;l++) {
					double B[3][3][3];
					for (int i=0;i<3;i++) for (int j=0;j<3;j++) for (int k=0;k<3;k++) B[i][j][k]=Bgr[i][j][k][l];
					double bv[6][3]; ToVoigtLocal(B, e1, e2, G0.n, bv);
					for (int r=0;r<6;r++) for (int c=0;c<3;c++) Bg[I][l*18+r*3+c]=bv[r][c];
				}
			}
			double Mw = alpha*V_K*Mmom; /* alpha (Eq.107) scales membrane stabilization to relieve locking */
			for (int l=0;l<2;l++) for (int I=0;I<nn;I++) for (int J=0;J<nn;J++) {
				double kij[3][3]={{0,0,0},{0,0,0},{0,0,0}};
				for (int a=0;a<6;a++) for (int b=0;b<6;b++) {
					double Cab=C[a][b]; if (Cab==0.0) continue;
					for (int ci=0;ci<3;ci++) for (int cj=0;cj<3;cj++) kij[ci][cj]+=Bg[I][l*18+a*3+ci]*Cab*Bg[J][l*18+b*3+cj];
				}
				int gi=3*nb[I], gj=3*nb[J];
				for (int ci=0;ci<3;ci++) for (int cj=0;cj<3;cj++) K[(size_t)(gi+ci)*ndof+(gj+cj)] += Mw*kij[ci][cj];
			}

			/* bending (curvature-gradient) stabilization: penalize d(kappa)/d(xi_l), the part
			 * of the strain gradient that lives in the xi3-linear (bending) term and is invisible
			 * to the membrane (xi3=0) stabilization. Coefficient M_aa * (h^3/12) is the bending
			 * counterpart of the membrane M_aa * h; uses 3rd derivatives (DDDp). */
			std::vector<std::vector<double> > Bkv(nn, std::vector<double>(36));
			for (int I=0;I<nn;I++) {
				double Bk[3][3][3][2];
				/* DDDp comps (2D): 0:Psi,111  1:Psi,122  2:Psi,112  3:Psi,222 */
				BMatrixCurvatureGradient(G0, DDp(0,I),DDp(2,I),DDp(1,I),
					DDDp(0,I),DDDp(2,I),DDDp(1,I),DDDp(3,I), Bk);
				for (int l=0;l<2;l++) {
					double B[3][3][3];
					for (int i=0;i<3;i++) for (int j=0;j<3;j++) for (int k=0;k<3;k++) B[i][j][k]=Bk[i][j][k][l];
					double bv[6][3]; ToVoigtLocal(B, e1, e2, G0.n, bv);
					for (int r=0;r<6;r++) for (int c=0;c<3;c++) Bkv[I][l*18+r*3+c]=bv[r][c];
				}
			}
			double Mw_bend = (A_K*Mmom)*(h*h*h/12.0);
			for (int l=0;l<2;l++) for (int I=0;I<nn;I++) for (int J=0;J<nn;J++) {
				double kij[3][3]={{0,0,0},{0,0,0},{0,0,0}};
				for (int a=0;a<6;a++) for (int b=0;b<6;b++) {
					double Cab=C[a][b]; if (Cab==0.0) continue;
					for (int ci=0;ci<3;ci++) for (int cj=0;cj<3;cj++) kij[ci][cj]+=Bkv[I][l*18+a*3+ci]*Cab*Bkv[J][l*18+b*3+cj];
				}
				int gi=3*nb[I], gj=3*nb[J];
				for (int ci=0;ci<3;ci++) for (int cj=0;cj<3;cj++) K[(size_t)(gi+ci)*ndof+(gj+cj)] += Mw_bend*kij[ci][cj];
			}
		}
	}
	if (nskip > 0) printf("  [AssembleShell: %d/%d nodes skipped (SetField failed)]\n", nskip, N);
	for (long i=0;i<ndof;i++) for (long j=i+1;j<ndof;j++) {
		double a = 0.5*(K[i*ndof+j]+K[j*ndof+i]);
		K[i*ndof+j] = K[j*ndof+i] = a;
	}
}

/* cyclic Jacobi eigenvalues (ascending) of a symmetric n x n matrix */
void JacobiEig(std::vector<double> A, int n, std::vector<double>& ev)
{
	for (int sweep=0;sweep<80;sweep++) {
		double off=0.0;
		for (int p=0;p<n;p++) for (int q=p+1;q<n;q++) off += A[(size_t)p*n+q]*A[(size_t)p*n+q];
		if (off < 1e-26) break;
		for (int p=0;p<n;p++) for (int q=p+1;q<n;q++) {
			double apq=A[(size_t)p*n+q];
			if (std::fabs(apq) < 1e-300) continue;
			double th=(A[(size_t)q*n+q]-A[(size_t)p*n+p])/(2*apq);
			double t=(th>=0?1.0:-1.0)/(std::fabs(th)+std::sqrt(th*th+1));
			double c=1/std::sqrt(t*t+1), sn=t*c;
			for (int k=0;k<n;k++){double a1=A[(size_t)k*n+p],a2=A[(size_t)k*n+q];A[(size_t)k*n+p]=c*a1-sn*a2;A[(size_t)k*n+q]=sn*a1+c*a2;}
			for (int k=0;k<n;k++){double a1=A[(size_t)p*n+k],a2=A[(size_t)q*n+k];A[(size_t)p*n+k]=c*a1-sn*a2;A[(size_t)q*n+k]=sn*a1+c*a2;}
		}
	}
	ev.resize(n);
	for (int i=0;i<n;i++) ev[i]=A[(size_t)i*n+i];
	std::sort(ev.begin(), ev.end());
}

/* number of near-zero eigenvalues of a free cylinder patch (expect 6 rigid-body modes) */
int CurvedPatchZeroModes(int nt, int nz, std::vector<double>& smallest9, double supportFac = 3.0,
	double thick = 0.05)
{
	const double Rc = 1.0;
	double arc = 1.0/(nt-1);          /* small arc span ~1 rad total */
	double dz = 1.0/(nz-1);
	int N = nt*nz;
	std::vector<double> X(3*N), area(N);
	for (int i=0;i<nt;i++) for (int j=0;j<nz;j++) {
		int id=j*nt+i;
		double th = i*arc;
		X[3*id]=Rc*std::cos(th); X[3*id+1]=Rc*std::sin(th); X[3*id+2]=j*dz;
		area[id]=(Rc*arc)*dz;
	}
	double C[6][6]; PlaneStressTangent(1.0e4, 0.0, C);
	double support = supportFac*std::max(Rc*arc, dz);
	std::vector<double> K;
	AssembleShell(X, area, thick, C, support, K);
	int ndof=3*N;
	std::vector<double> ev;
	JacobiEig(K, ndof, ev);
	double lmax=ev.back();
	int nz0=0;
	for (int i=0;i<ndof;i++) if (ev[i] < 1e-8*lmax) nz0++;
	smallest9.assign(ev.begin(), ev.begin()+std::min(9,ndof));
	for (double& e : smallest9) e/=lmax;
	return nz0;
}

/* Diagnostic: does the assembler's RKPM (MLSSolverT+Gaussian) recover a cylinder's
 * principal curvatures {1/R, 0} from nodal coordinates? Reconstructs x1,x2,x11,x22,x12 at
 * an interior node and computes curvatures via the fundamental forms. */
void CylinderCurvatureRKPM(double R, int nt, int nz, double supportFac,
	double& kmax, double& kmin)
{
	const double Lz = 4.0*R, phiMax = 40.0*M_PI/180.0;
	int N = nt*nz;
	std::vector<double> X(3*N);
	double dphi=(2*phiMax)/(nt-1), dz=Lz/(nz-1), arc=R*dphi;
	for (int i=0;i<nt;i++) for (int j=0;j<nz;j++) {
		int id=j*nt+i; double phi=-phiMax+i*dphi;
		X[3*id]=R*std::sin(phi); X[3*id+1]=R*std::cos(phi); X[3*id+2]=j*dz;
	}
	double support = supportFac*std::max(arc,dz);
	dArrayT gwin(3); gwin[0]=1.5; gwin[1]=0.4; gwin[2]=3.0;
	MLSSolverT rkpm(2,2,false,MeshFreeT::kGaussian,gwin); rkpm.Initialize();

	int P = (nz/2)*nt + nt/2; /* interior node */
	std::vector<int> nb;
	for (int Q=0;Q<N;Q++){double d0=X[3*Q]-X[3*P],d1=X[3*Q+1]-X[3*P+1],d2=X[3*Q+2]-X[3*P+2];
		if (std::sqrt(d0*d0+d1*d1+d2*d2)<0.99*support) nb.push_back(Q);}
	int nn=int(nb.size());
	std::vector<double> nbX(3*nn);
	for (int k=0;k<nn;k++){nbX[3*k]=X[3*nb[k]];nbX[3*k+1]=X[3*nb[k]+1];nbX[3*k+2]=X[3*nb[k]+2];}
	double psi1[3],psi2[3],n0[3]; PCAFrame(&nbX[0],nn,psi1,psi2,n0);
	dArray2DT lc(nn,2);
	for (int k=0;k<nn;k++){double dv[3]={X[3*nb[k]]-X[3*P],X[3*nb[k]+1]-X[3*P+1],X[3*nb[k]+2]-X[3*P+2]};
		lc(k,0)=Dot(dv,psi1); lc(k,1)=Dot(dv,psi2);}
	dArray2DT np(nn,1); np=support; dArrayT vol(nn); vol=arc*dz; dArrayT s(2); s[0]=s[1]=0.0;
	rkpm.SetField(lc,np,vol,s,2);
	const dArray2DT& Dp=rkpm.Dphi(); const dArray2DT& DDp=rkpm.DDphi();
	double x1[3]={0,0,0},x2[3]={0,0,0},x11[3]={0,0,0},x22[3]={0,0,0},x12[3]={0,0,0};
	for (int I=0;I<nn;I++){double Xq[3]={X[3*nb[I]],X[3*nb[I]+1],X[3*nb[I]+2]};
		for(int d=0;d<3;d++){x1[d]+=Dp(0,I)*Xq[d];x2[d]+=Dp(1,I)*Xq[d];
			x11[d]+=DDp(0,I)*Xq[d];x22[d]+=DDp(1,I)*Xq[d];x12[d]+=DDp(2,I)*Xq[d];}}
	double nv[3]; Cross(x1,x2,nv); double nm=Norm(nv); for(int d=0;d<3;d++) nv[d]/=nm;
	double E=Dot(x1,x1),F=Dot(x1,x2),G=Dot(x2,x2);
	double Lf=Dot(x11,nv),Mf=Dot(x12,nv),Nf=Dot(x22,nv);
	double a2=E*G-F*F, a1=-(E*Nf-2*F*Mf+G*Lf), a0=Lf*Nf-Mf*Mf;
	double disc=a1*a1-4*a2*a0; if(disc<0)disc=0;
	double k1=std::fabs((-a1+std::sqrt(disc))/(2*a2)), k2=std::fabs((-a1-std::sqrt(disc))/(2*a2));
	kmax=std::max(k1,k2); kmin=std::min(k1,k2);
}

/* Curved-shell stiffness by BACKGROUND-CELL GAUSS integration over a structured nt x nz
 * node grid (classic meshfree Galerkin). Each grid cell gets 2x2 in-plane Gauss points;
 * the RKPM shape functions are evaluated AT each Gauss point (not at nodes), so the
 * integration is consistent and stable -- no nodal-integration hourglass, no stabilization.
 * Through-thickness: 3-pt Gauss with sigma33=0 plane stress (via the condensed C). */
void AssembleShellGauss(const std::vector<double>& X, int nt, int nz, double h,
	const double C[6][6], double support, std::vector<double>& K)
{
	int N = nt*nz, ndof = 3*N;
	K.assign((size_t)ndof*ndof, 0.0);

	double xg3[3] = {-std::sqrt(3.0/5.0), 0.0, std::sqrt(3.0/5.0)};
	double wg3[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0};
	double g2[2] = {-1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0)};

	dArrayT gwin(3); gwin[0]=1.5; gwin[1]=0.4; gwin[2]=3.0;
	MLSSolverT rkpm(2, 2, false, MeshFreeT::kGaussian, gwin);
	rkpm.Initialize();

	for (int ci=0; ci<nt-1; ci++)
		for (int cj=0; cj<nz-1; cj++) {

			int corner[4] = { cj*nt+ci, cj*nt+ci+1, (cj+1)*nt+ci+1, (cj+1)*nt+ci };

			for (int gi=0; gi<2; gi++)
				for (int gj=0; gj<2; gj++) {

					double xi=g2[gi], eta=g2[gj];
					double Nc[4]   = {(1-xi)*(1-eta)/4,(1+xi)*(1-eta)/4,(1+xi)*(1+eta)/4,(1-xi)*(1+eta)/4};
					double dNx[4]  = {-(1-eta)/4,(1-eta)/4,(1+eta)/4,-(1+eta)/4};
					double dNe[4]  = {-(1-xi)/4,-(1+xi)/4,(1+xi)/4,(1-xi)/4};
					double xgp[3]={0,0,0}, dxi[3]={0,0,0}, deta[3]={0,0,0};
					for (int k=0;k<4;k++) for (int d=0;d<3;d++) {
						xgp[d]  += Nc[k]*X[3*corner[k]+d];
						dxi[d]  += dNx[k]*X[3*corner[k]+d];
						deta[d] += dNe[k]*X[3*corner[k]+d];
					}
					double cr[3]; Cross(dxi, deta, cr);
					double wIP = Norm(cr); /* area element; 2pt Gauss weights are 1 */

					/* neighbors of the Gauss point within the support */
					std::vector<int> nb;
					for (int Q=0;Q<N;Q++) {
						double d0=X[3*Q]-xgp[0], d1=X[3*Q+1]-xgp[1], d2=X[3*Q+2]-xgp[2];
						if (std::sqrt(d0*d0+d1*d1+d2*d2) < 0.99*support) nb.push_back(Q);
					}
					int nn=int(nb.size());
					if (nn < 6) continue;

					std::vector<double> nbX(3*nn);
					for (int k=0;k<nn;k++){nbX[3*k]=X[3*nb[k]];nbX[3*k+1]=X[3*nb[k]+1];nbX[3*k+2]=X[3*nb[k]+2];}
					double psi1[3],psi2[3],n0[3];
					PCAFrame(&nbX[0], nn, psi1, psi2, n0);

					dArray2DT lc(nn,2);
					for (int k=0;k<nn;k++) {
						double dv[3]={X[3*nb[k]]-xgp[0],X[3*nb[k]+1]-xgp[1],X[3*nb[k]+2]-xgp[2]};
						lc(k,0)=Dot(dv,psi1); lc(k,1)=Dot(dv,psi2);
					}
					dArray2DT np(nn,1); np=support;
					dArrayT vol(nn); vol=wIP; /* nominal nodal volume for the moment matrix */
					dArrayT sample(2); sample[0]=0.0; sample[1]=0.0;
					if (!rkpm.SetField(lc, np, vol, sample, 2)) continue;
					const dArray2DT& Dp = rkpm.Dphi();
					const dArray2DT& DDp = rkpm.DDphi();

					double x1[3]={0,0,0},x2[3]={0,0,0},x11[3]={0,0,0},x22[3]={0,0,0},x12[3]={0,0,0};
					for (int I=0;I<nn;I++) {
						double Xq[3]={X[3*nb[I]],X[3*nb[I]+1],X[3*nb[I]+2]};
						for (int d=0;d<3;d++){
							x1[d]+=Dp(0,I)*Xq[d]; x2[d]+=Dp(1,I)*Xq[d];
							x11[d]+=DDp(0,I)*Xq[d]; x22[d]+=DDp(1,I)*Xq[d]; x12[d]+=DDp(2,I)*Xq[d];
						}
					}

					for (int g=0; g<3; g++) {
						ShellGeom G;
						if (!BuildGeom(x1,x2,x11,x22,x12,h,xg3[g],G)) continue;
						double cw = wg3[g]*(h/2.0)*wIP;
						std::vector<std::vector<double> > Bv(nn, std::vector<double>(18));
						for (int I=0;I<nn;I++) {
							double B[3][3][3];
							BMatrix(G, Dp(0,I),Dp(1,I),DDp(0,I),DDp(2,I),DDp(1,I), B);
							double bv[6][3]; ToVoigt(B, bv);
							for (int r=0;r<6;r++) for (int c=0;c<3;c++) Bv[I][r*3+c]=bv[r][c];
						}
						for (int I=0;I<nn;I++) for (int J=0;J<nn;J++) {
							double kij[3][3]={{0,0,0},{0,0,0},{0,0,0}};
							for (int a=0;a<6;a++) for (int b=0;b<6;b++) {
								double Cab=C[a][b]; if (Cab==0.0) continue;
								for (int ci2=0;ci2<3;ci2++) for (int cj2=0;cj2<3;cj2++) kij[ci2][cj2]+=Bv[I][a*3+ci2]*Cab*Bv[J][b*3+cj2];
							}
							int gI=3*nb[I], gJ=3*nb[J];
							for (int ci2=0;ci2<3;ci2++) for (int cj2=0;cj2<3;cj2++) K[(size_t)(gI+ci2)*ndof+(gJ+cj2)] += cw*kij[ci2][cj2];
						}
					}
				}
		}
	for (long i=0;i<ndof;i++) for (long j=i+1;j<ndof;j++) {
		double a=0.5*(K[i*ndof+j]+K[j*ndof+i]); K[i*ndof+j]=K[j*ndof+i]=a;
	}
}

/* Curved-shell stiffness using a SINGLE GLOBAL parametric chart (the structured (i,j) grid)
 * instead of a per-evaluation-point PCA frame. Because every Gauss point uses the same
 * consistent parameterization, the RKPM shape functions form a consistent global Galerkin
 * basis -> no spurious rank deficiency. The 3D geometry/curvature is reconstructed from the
 * nodal 3D positions via the RKPM derivatives in (u,v); the area element is |x,u x x,v|.
 * 2x2 in-plane Gauss per parametric cell, 3-pt through-thickness Gauss. */
void AssembleShellParam(const std::vector<double>& X, int nt, int nz, double h,
	const double C[6][6], double supportFac, std::vector<double>& K)
{
	int N=nt*nz, ndof=3*N; K.assign((size_t)ndof*ndof,0.0);
	double xg3[3]={-std::sqrt(3.0/5.0),0.0,std::sqrt(3.0/5.0)}, wg3[3]={5.0/9.0,8.0/9.0,5.0/9.0};
	double g2[2]={-1.0/std::sqrt(3.0),1.0/std::sqrt(3.0)};
	double support = supportFac; /* in parametric (index) units */

	dArrayT gwin(3); gwin[0]=1.5; gwin[1]=0.4; gwin[2]=3.0;
	MLSSolverT rkpm(2,2,false,MeshFreeT::kGaussian,gwin); rkpm.Initialize();

	/* node parametric coords (u,v) = (i,j) */
	for (int cu=0; cu<nt-1; cu++)
		for (int cv=0; cv<nz-1; cv++)
			for (int gu=0; gu<2; gu++)
				for (int gv=0; gv<2; gv++) {

					double ug = cu + 0.5*(1.0+g2[gu]);  /* parametric position of the Gauss point */
					double vg = cv + 0.5*(1.0+g2[gv]);

					std::vector<int> nb;
					for (int j=0;j<nz;j++) for (int i=0;i<nt;i++) {
						double du=i-ug, dv=j-vg;
						if (std::sqrt(du*du+dv*dv) < 0.99*support) nb.push_back(j*nt+i);
					}
					int nn=int(nb.size());
					if (nn<6) continue;

					dArray2DT lc(nn,2);
					for (int k=0;k<nn;k++) { int i=nb[k]%nt, j=nb[k]/nt; lc(k,0)=i-ug; lc(k,1)=j-vg; }
					dArray2DT np(nn,1); np=support; dArrayT vol(nn); vol=1.0; dArrayT s(2); s[0]=s[1]=0.0;
					if (!rkpm.SetField(lc,np,vol,s,2)) continue;
					const dArray2DT& Dp=rkpm.Dphi(); const dArray2DT& DDp=rkpm.DDphi();

					double x1[3]={0,0,0},x2[3]={0,0,0},x11[3]={0,0,0},x22[3]={0,0,0},x12[3]={0,0,0};
					for (int I=0;I<nn;I++){double Xq[3]={X[3*nb[I]],X[3*nb[I]+1],X[3*nb[I]+2]};
						for(int d=0;d<3;d++){x1[d]+=Dp(0,I)*Xq[d];x2[d]+=Dp(1,I)*Xq[d];
							x11[d]+=DDp(0,I)*Xq[d];x22[d]+=DDp(1,I)*Xq[d];x12[d]+=DDp(2,I)*Xq[d];}}
					double cr[3]; Cross(x1,x2,cr); double jac=Norm(cr); /* physical area per unit (u,v) */
					double wIP = jac*0.25; /* 2x2 Gauss on a unit parametric cell: weight 1*1, du dv = 1/4 */

					for (int g=0;g<3;g++) {
						ShellGeom G;
						if (!BuildGeom(x1,x2,x11,x22,x12,h,xg3[g],G)) continue;
						double e1[3],e2[3]; OrthoTangents(G.n, e1, e2); /* local frame, e3 = normal */
						double cw = wg3[g]*(h/2.0)*wIP;
						std::vector<std::vector<double> > Bv(nn, std::vector<double>(18));
						for (int I=0;I<nn;I++){double B[3][3][3];
							BMatrix(G,Dp(0,I),Dp(1,I),DDp(0,I),DDp(2,I),DDp(1,I),B);
							double bv[6][3]; ToVoigtLocal(B, e1, e2, G.n, bv); /* plane stress along normal */
							for(int r=0;r<6;r++)for(int c=0;c<3;c++) Bv[I][r*3+c]=bv[r][c];}
						for (int I=0;I<nn;I++) for (int J=0;J<nn;J++){
							double kij[3][3]={{0,0,0},{0,0,0},{0,0,0}};
							for(int a=0;a<6;a++)for(int b=0;b<6;b++){double Cab=C[a][b];if(Cab==0.0)continue;
								for(int ci=0;ci<3;ci++)for(int cj=0;cj<3;cj++) kij[ci][cj]+=Bv[I][a*3+ci]*Cab*Bv[J][b*3+cj];}
							int gI=3*nb[I], gJ=3*nb[J];
							for(int ci=0;ci<3;ci++)for(int cj=0;cj<3;cj++) K[(size_t)(gI+ci)*ndof+(gJ+cj)]+=cw*kij[ci][cj];
						}
					}
				}
	for (long i=0;i<ndof;i++) for (long j=i+1;j<ndof;j++){double a=0.5*(K[i*ndof+j]+K[j*ndof+i]);K[i*ndof+j]=K[j*ndof+i]=a;}
}

/* Quarter-cylinder under internal pressure (pure MEMBRANE, plane strain u_z=0).
 * Known: radial expansion u_r = p R^2 / (E h) (nu=0). Returns the radial displacement at an
 * interior node. Tests curved-membrane stiffness in isolation. */
double QuarterCylinderPressure(int nt, int nz)
{
	const double R = 10.0, Lz = 20.0, h = 0.5, E = 1.0e4, nu = 0.0, p = 1.0;
	int N = nt*nz;
	std::vector<double> X(3*N), area(N);
	double dth = (0.5*M_PI)/(nt-1), dz = Lz/(nz-1), arc = R*dth;
	for (int i=0;i<nt;i++) for (int j=0;j<nz;j++) {
		int id=j*nt+i; double th=i*dth;
		X[3*id]=R*std::cos(th); X[3*id+1]=R*std::sin(th); X[3*id+2]=j*dz;
		double wa=(i==0||i==nt-1)?0.5:1.0, wb=(j==0||j==nz-1)?0.5:1.0;
		area[id]=arc*dz*wa*wb;
	}
	double C[6][6]; PlaneStressTangent(E, nu, C);
	(void)arc;
	std::vector<double> K; AssembleShellParam(X, nt, nz, h, C, 3.0, K);
	int ndof=3*N;

	std::vector<char> fixed(ndof,0);
	for (int i=0;i<nt;i++) for (int j=0;j<nz;j++){ int id=j*nt+i; fixed[3*id+2]=1; } /* u_z=0 plane strain */
	for (int j=0;j<nz;j++){ fixed[3*(j*nt+0)+1]=1; fixed[3*(j*nt+nt-1)+0]=1; } /* symmetry edges */

	std::vector<double> f(ndof,0.0);
	for (int i=0;i<nt;i++) for (int j=0;j<nz;j++){ int id=j*nt+i; double th=i*dth;
		f[3*id+0]+=p*area[id]*std::cos(th); f[3*id+1]+=p*area[id]*std::sin(th); }

	std::vector<int> map(ndof,-1); int nf=0;
	for (int i=0;i<ndof;i++) if(!fixed[i]) map[i]=nf++;
	std::vector<double> Kr((size_t)nf*nf,0.0), fr(nf,0.0);
	for (int i=0;i<ndof;i++){ if(fixed[i])continue; fr[map[i]]=f[i];
		for(int j=0;j<ndof;j++) if(!fixed[j]) Kr[(size_t)map[i]*nf+map[j]]=K[(size_t)i*ndof+j]; }
	std::vector<double> ur;
	if(!SolveDense(Kr,fr,nf,ur)) return -1.0;
	int id=(nz/2)*nt + nt/2; double th=(nt/2)*dth;
	double ux=ur[map[3*id+0]], uy=ur[map[3*id+1]];
	return ux*std::cos(th)+uy*std::sin(th); /* radial component */
}

/* FULL pinched cylinder (no symmetry planes -> avoids the rotation-free symmetry-BC issue;
 * the 3D-distance neighbor search handles the circumferential seam automatically). Diaphragms
 * at both ends, two opposite inward point loads at mid-span. R=300,L=600,h=3,E=3e6,nu=0.3,
 * P=1; reference radial deflection under the load = 1.8248e-5. nt must be even (load at
 * theta=0 and pi), nz odd (load at z=L/2). */
double FullPinchedCylinder(int nt, int nz, double supportFac = 3.0)
{
	const double R=300.0, L=600.0, h=3.0, E=3.0e6, nu=0.3, P=1.0;
	int N=nt*nz;
	std::vector<double> X(3*N), area(N);
	double dth=(2.0*M_PI)/nt, dz=L/(nz-1), arc=R*dth;
	for (int i=0;i<nt;i++) for (int j=0;j<nz;j++){
		int id=j*nt+i; double th=i*dth;
		X[3*id]=R*std::cos(th); X[3*id+1]=R*std::sin(th); X[3*id+2]=j*dz;
		double wb=(j==0||j==nz-1)?0.5:1.0;
		area[id]=arc*dz*wb; /* full circumference: no circumferential edge halving */
	}
	double C[6][6]; PlaneStressTangent(E, nu, C);
	double support=supportFac*std::max(arc,dz);
	std::vector<double> K; AssembleShell(X, area, h, C, support, K);
	int ndof=3*N;

	std::vector<char> fixed(ndof,0);
	for (int i=0;i<nt;i++){ int id0=0*nt+i, idL=(nz-1)*nt+i;  /* diaphragms at z=0 and z=L */
		fixed[3*id0+0]=1; fixed[3*id0+1]=1; fixed[3*idL+0]=1; fixed[3*idL+1]=1; }
	fixed[3*(0*nt+0)+2]=1; /* one u_z restraint -> remove axial rigid translation */

	std::vector<double> f(ndof,0.0);
	int jMid=(nz-1)/2;
	int idTop=jMid*nt + 0;       /* theta=0  -> (R,0,zmid),  load inward = -x */
	int idBot=jMid*nt + nt/2;    /* theta=pi -> (-R,0,zmid), load inward = +x */
	f[3*idTop+0] += -P;
	f[3*idBot+0] += +P;

	std::vector<int> map(ndof,-1); int nf=0;
	for (int i=0;i<ndof;i++) if(!fixed[i]) map[i]=nf++;
	std::vector<double> Kr((size_t)nf*nf,0.0), fr(nf,0.0);
	for (int i=0;i<ndof;i++){ if(fixed[i])continue; fr[map[i]]=f[i];
		for(int j=0;j<ndof;j++) if(!fixed[j]) Kr[(size_t)map[i]*nf+map[j]]=K[(size_t)i*ndof+j]; }
	std::vector<double> ur;
	if(!SolveDense(Kr,fr,nf,ur)) return -1.0;
	return -ur[map[3*idTop+0]]; /* inward (-x) deflection magnitude at the load */
}

/* Pinched cylinder with end diaphragms (Belytschko obstacle course). Two opposite radial
 * point loads at mid-span; ends are rigid diaphragms. Modeled as 1/8 by symmetry. Standard
 * data: R=300, L=600, h=3, E=3e6, nu=0.3, P=1; reference radial deflection under the load
 * = 1.8248e-5. Returns the deflection magnitude under the load. */
double PinchedCylinder(int nt, int nz, double supportFac = 3.0)
{
	const double R=300.0, L=600.0, h=3.0, E=3.0e6, nu=0.3, P=1.0;
	int N=nt*nz;
	std::vector<double> X(3*N), area(N);
	double dphi=(0.5*M_PI)/(nt-1), dz=(0.5*L)/(nz-1), arc=R*dphi;
	for (int i=0;i<nt;i++) for (int j=0;j<nz;j++){
		int id=j*nt+i; double phi=i*dphi;            /* phi=0 at top (load), phi=90 at side */
		X[3*id]=R*std::sin(phi); X[3*id+1]=R*std::cos(phi); X[3*id+2]=j*dz; /* z=0 mid-span */
		double wa=(i==0||i==nt-1)?0.5:1.0, wb=(j==0||j==nz-1)?0.5:1.0;
		area[id]=arc*dz*wa*wb;
	}
	double C[6][6]; PlaneStressTangent(E, nu, C);
	double support=supportFac*std::max(arc,dz);
	std::vector<double> K; AssembleShell(X, area, h, C, support, K);
	int ndof=3*N;

	std::vector<char> fixed(ndof,0);
	for (int i=0;i<nt;i++) { int id=(nz-1)*nt+i; fixed[3*id+0]=1; fixed[3*id+1]=1; } /* z=L/2 diaphragm */
	for (int i=0;i<nt;i++) { int id=0*nt+i; fixed[3*id+2]=1; }                       /* z=0 symmetry: u_z=0 */
	for (int j=0;j<nz;j++) { int id=j*nt+0;      fixed[3*id+0]=1; }                  /* phi=0 symmetry: u_x=0 */
	for (int j=0;j<nz;j++) { int id=j*nt+(nt-1); fixed[3*id+1]=1; }                  /* phi=90 symmetry: u_y=0 */

	std::vector<double> f(ndof,0.0);
	f[3*(0*nt+0)+1] = -0.25*P; /* load corner (top, mid-span): P/4 radial inward (-y) for the 1/8 model */

	std::vector<int> map(ndof,-1); int nf=0;
	for (int i=0;i<ndof;i++) if(!fixed[i]) map[i]=nf++;
	std::vector<double> Kr((size_t)nf*nf,0.0), fr(nf,0.0);
	for (int i=0;i<ndof;i++){ if(fixed[i])continue; fr[map[i]]=f[i];
		for(int j=0;j<ndof;j++) if(!fixed[j]) Kr[(size_t)map[i]*nf+map[j]]=K[(size_t)i*ndof+j]; }
	std::vector<double> ur;
	if(!SolveDense(Kr,fr,nf,ur)) return -1.0;
	return -ur[map[3*(0*nt+0)+1]]; /* inward radial (-y) deflection magnitude at the load */
}

/* FLAT cantilever strip via the SAME curved assembler (AssembleShell): clamped at x=0,
 * transverse tip load at x=L. Returns tip deflection. Beam theory: w = P L^3/(3 E I),
 * I = W h^3/12. Isolates whether AssembleShell itself is consistent (flat) before blaming
 * curvature. */
double FlatCantilever(int nx, int nz, double h, double supportFac = 5.0)
{
	const double L = 8.0, W = 8.0, E = 1.0e4, nu = 0.0, P = 1.0; /* square -> isotropic spacing */
	double dx = L/(nx-1), dz = W/(nz-1);
	int N = nx*nz;
	std::vector<double> X(3*N), area(N);
	for (int i=0;i<nx;i++) for (int j=0;j<nz;j++) {
		int id=j*nx+i;
		X[3*id]=i*dx; X[3*id+1]=0.0; X[3*id+2]=j*dz; /* flat plate in x-z plane, normal = y */
		double wa=(i==0||i==nx-1)?0.5:1.0, wb=(j==0||j==nz-1)?0.5:1.0;
		area[id]=dx*dz*wa*wb;
	}
	double C[6][6]; PlaneStressTangent(E, nu, C);
	(void)dx; (void)dz;
	std::vector<double> K; AssembleShellParam(X, nx, nz, h, C, supportFac, K);
	int ndof=3*N;

	std::vector<char> fixed(ndof,0);
	for (int j=0;j<nz;j++) for (int ir=0; ir<2; ir++) { /* clamp 2 rows at x=0 */
		int id=j*nx+ir; fixed[3*id]=fixed[3*id+1]=fixed[3*id+2]=1;
	}
	std::vector<double> f(ndof,0.0);
	for (int j=0;j<nz;j++) { int id=j*nx+(nx-1); double wb=(j==0||j==nz-1)?0.5:1.0; f[3*id+1]=-P*wb/(nz-1); }

	std::vector<int> map(ndof,-1); int nf=0;
	for (int i=0;i<ndof;i++) if(!fixed[i]) map[i]=nf++;
	std::vector<double> Kr((size_t)nf*nf,0.0), fr(nf,0.0);
	for (int i=0;i<ndof;i++){ if(fixed[i])continue; fr[map[i]]=f[i];
		for(int j=0;j<ndof;j++) if(!fixed[j]) Kr[(size_t)map[i]*nf+map[j]]=K[(size_t)i*ndof+j]; }
	std::vector<double> ur;
	if (!SolveDense(Kr,fr,nf,ur)) return 0.0;
	int idTip = ((nz-1)/2)*nx + (nx-1);
	return -ur[map[3*idTip+1]];
}

/* Scordelis-Lo roof: returns the downward deflection at the free-edge midpoint.
 * nt = nodes across the 80-degree arc, nz = nodes along the length. */
double ScordelisLo(int nt, int nz, double supportFac = 3.0, double alpha = 1.0)
{
	const double R = 25.0, Lz = 50.0, h = 0.25;
	const double E = 4.32e8, nu = 0.0, grav = 90.0; /* load per unit area, downward */
	const double phiMax = 40.0*M_PI/180.0;          /* +/-40 deg from the crown */

	int N = nt*nz;
	std::vector<double> X(3*N);
	std::vector<double> area(N, 0.0);
	double dphi = (2*phiMax)/(nt-1);
	double dz = Lz/(nz-1);
	double arc = R*dphi;                             /* circumferential spacing */
	for (int i=0;i<nt;i++)
		for (int j=0;j<nz;j++) {
			int id = j*nt+i;
			double phi = -phiMax + i*dphi;
			X[3*id]   = R*std::sin(phi);             /* x */
			X[3*id+1] = R*std::cos(phi);             /* y (up); crown at y=R */
			X[3*id+2] = j*dz;                        /* z (length) */
			/* nodal area = (circumferential spacing) x (axial spacing), halved at edges */
			double wa = (i==0||i==nt-1) ? 0.5 : 1.0;
			double wb = (j==0||j==nz-1) ? 0.5 : 1.0;
			area[id] = (arc*dz)*wa*wb;
		}

	double C[6][6];
	PlaneStressTangent(E, nu, C);
	double support = supportFac*std::max(arc, dz);
	std::vector<double> K;
	AssembleShell(X, area, h, C, support, K, alpha); /* nodal integration (reduced) + stabilization */
	int ndof = 3*N;

	/* boundary conditions (standard Scordelis-Lo):
	 *   rigid end diaphragms at z=0 and z=Lz: u_x = u_y = 0, u_z FREE (a simple support, not a
	 *   clamp -- a diaphragm does not restrain rotation, so one node ring per end is correct);
	 *   remove the remaining axial rigid translation by fixing u_z at a single node. */
	std::vector<char> fixed(ndof, 0);
	for (int i=0;i<nt;i++) {
		int id0 = 0*nt+i;          /* z=0 diaphragm */
		int idL = (nz-1)*nt+i;     /* z=Lz diaphragm */
		fixed[3*id0+0]=1; fixed[3*id0+1]=1;
		fixed[3*idL+0]=1; fixed[3*idL+1]=1;
	}
	fixed[3*(0*nt+0)+2]=1;         /* single u_z restraint: remove axial rigid translation */

	/* gravity: nodal force = -grav * area in the y (vertical) direction */
	std::vector<double> f(ndof, 0.0);
	for (int p=0;p<N;p++) f[3*p+1] = -grav*area[p];

	std::vector<int> map(ndof, -1);
	int nf = 0;
	for (int i=0;i<ndof;i++) if (!fixed[i]) map[i]=nf++;
	std::vector<double> Kr((size_t)nf*nf, 0.0), fr(nf, 0.0);
	for (int i=0;i<ndof;i++) {
		if (fixed[i]) continue;
		fr[map[i]] = f[i];
		for (int j=0;j<ndof;j++) if (!fixed[j]) Kr[(size_t)map[i]*nf+map[j]] = K[(size_t)i*ndof+j];
	}
	std::vector<double> ur;
	if (!SolveDense(Kr, fr, nf, ur)) return 0.0;

	/* vertical deflection at the free-edge midpoint (phi = +phiMax, z = Lz/2) */
	int iEdge = nt-1;
	int jMid = (nz-1)/2;
	int idEdge = jMid*nt + iEdge;
	double uy = ur[map[3*idEdge+1]];
	return -uy; /* downward magnitude */
}

} /* anonymous namespace */

/* A FREE curved (cylinder) patch must have exactly 6 zero-energy (rigid-body) modes.
 *
 * History: with EFG and membrane-only stabilization this gave 9 modes at 7x7 and 16 at 9x9
 * (a curved-shell BENDING hourglass that grows with refinement). The cure is the bending
 * (curvature-gradient) stabilization (BMatrixCurvatureGradient, using DDDphi) PLUS a
 * sufficiently large support (~5x spacing) so the quadratic RKPM strain operator is rich
 * enough — small supports (e.g. 3x) leave a resonance with spurious modes. With both, the
 * patch is rank-clean (exactly 6) across refinement. */
/* DISABLED diagnostic. The pinched cylinder (R/h=100, two opposite point loads) is
 * BENDING / inextensional dominated -- unlike Scordelis-Lo (membrane dominated) which works.
 * Both 1/8-symmetry and full models come out ~100-700x too soft and grow with refinement:
 * a spurious inextensional-bending hourglass.
 *
 * Reviewer-guided audit (2026-06-07) -- two levers pulled:
 *   - SUPPORT SWEEP (lever 1): the softness is strongly support-driven (11x at supportFac 1.8,
 *     ~840x at 3.4). Larger support flatlines the kernel -> the higher-order stabilization
 *     derivatives vanish -> hourglass uncontrolled. (Paper section 5.1 / Approach 1 limit.)
 *   - THROUGH-THICKNESS B,xi COUPLING (lever 2 / Check 1): added the dropped auxiliary-tensor
 *     gradient terms (dB1/dxi_l = B1m[l]) to BMatrixCurvatureGradient. Negligible effect
 *     (identical to 5 sig figs) -- and a 50x boost of the whole bending-stab term only helps
 *     ~4x. So the first-gradient bending stabilization is simply too weak (~h^3) to control
 *     this inextensional hourglass, regardless of its exact form.
 *
 * Conclusion: needs the paper's section 5.2 (Approach 2) stabilization or a higher-completeness
 * (cubic) basis -- not a tweak to the current scheme. Membrane curved shells work
 * (Scordelis-Lo); bending-dominated curved shells remain blocked. Tracked in #66/#68. */
TEST(KLShellBenchmark, DISABLED_PinchedCylinderAudit)
{
	double ref = 1.8248e-5;
	printf("FULL pinched cylinder (ref %.4e) -- 32x17, support sweep (added B1m,l terms):\n", ref);
	for (double sf = 1.8; sf <= 3.21; sf += 0.4) {
		double w = FullPinchedCylinder(32, 17, sf);
		printf("  supportFac=%.1f : w=%.4e  (%.0f%% of ref)\n", sf, w, 100.0*w/ref);
	}
}

TEST(KLShellBenchmark, DISABLED_CurvedPatchBendingStabilized)
{
	std::vector<double> s9;
	for (double sf = 3.0; sf <= 5.01; sf += 1.0) {
		int z7  = CurvedPatchZeroModes(7,  7,  s9, sf);
		int z9  = CurvedPatchZeroModes(9,  9,  s9, sf);
		printf("support=%.1f: 7x7=%d 9x9=%d  (7th/lmax=%.2e)\n", sf, z7, z9, s9[6]);
	}
}

/* Diagnostic audits used to root-cause the curved-shell BVP (now RESOLVED -- see
 * ScordelisLoRoof). The decisive findings were: the assembler is correct on flat problems
 * (FlatCantileverAudit), curvature and strain are reconstructed correctly (CylinderCurvature
 * Audit, RadialHoopStrain), but the curved BVP was garbage (CurvedMembraneAudit) -- which
 * traced to the plane-stress condensation being applied in the GLOBAL frame instead of the
 * shell-normal frame (fixed by ToVoigtLocal), plus an over-clamped Scordelis-Lo end (fixed
 * to a proper diaphragm). These remain as DISABLED regression diagnostics.
 *
 * Run with:  ./test_KLShellBenchmark --gtest_also_run_disabled_tests */
TEST(KLShellBenchmark, DISABLED_GaussCurvedSpectrum)
{
	/* free cylinder patch assembled by GAUSS integration -> count zero modes (expect 6) */
	for (int nn = 7; nn <= 11; nn += 2) {
		const double Rc=1.0; double arc=1.0/(nn-1), dz=1.0/(nn-1);
		int N=nn*nn; std::vector<double> X(3*N);
		for (int i=0;i<nn;i++) for (int j=0;j<nn;j++){int id=j*nn+i;double th=i*arc;
			X[3*id]=Rc*std::cos(th);X[3*id+1]=Rc*std::sin(th);X[3*id+2]=j*dz;}
		double C[6][6]; PlaneStressTangent(1.0e4,0.0,C);
		double support=3.0*std::max(Rc*arc,dz);
		std::vector<double> K; AssembleShellGauss(X,nn,nn,0.05,C,support,K);
		int ndof=3*N; std::vector<double> ev; JacobiEig(K,ndof,ev);
		double lmax=ev.back(); int z=0; for(int i=0;i<ndof;i++) if(ev[i]<1e-8*lmax) z++;
		printf("Gauss free cylinder %dx%d: zero modes=%d (expect 6)  ev[6]/lmax=%.2e\n", nn,nn,z, ev[6]/lmax);
	}
}

TEST(KLShellBenchmark, DISABLED_RadialHoopStrain)
{
	const double R = 10.0, w = 0.1;
	int nt=15, nz=15; double Lz=20.0;
	int N=nt*nz;
	std::vector<double> X(3*N);
	double dth=(0.5*M_PI)/(nt-1), dz=Lz/(nz-1), arc=R*dth;
	for (int i=0;i<nt;i++) for (int j=0;j<nz;j++){int id=j*nt+i;double th=i*dth;
		X[3*id]=R*std::cos(th);X[3*id+1]=R*std::sin(th);X[3*id+2]=j*dz;}
	double support=3.0*std::max(arc,dz);
	dArrayT gwin(3);gwin[0]=1.5;gwin[1]=0.4;gwin[2]=3.0;
	MLSSolverT rkpm(2,2,false,MeshFreeT::kGaussian,gwin); rkpm.Initialize();

	int P=(nz/2)*nt+nt/2;
	std::vector<int> nb;
	for(int Q=0;Q<N;Q++){double d0=X[3*Q]-X[3*P],d1=X[3*Q+1]-X[3*P+1],d2=X[3*Q+2]-X[3*P+2];
		if(std::sqrt(d0*d0+d1*d1+d2*d2)<0.99*support) nb.push_back(Q);}
	int nn=int(nb.size());
	std::vector<double> nbX(3*nn);
	for(int k=0;k<nn;k++){nbX[3*k]=X[3*nb[k]];nbX[3*k+1]=X[3*nb[k]+1];nbX[3*k+2]=X[3*nb[k]+2];}
	double psi1[3],psi2[3],n0[3]; PCAFrame(&nbX[0],nn,psi1,psi2,n0);
	dArray2DT lc(nn,2);
	for(int k=0;k<nn;k++){double dv[3]={X[3*nb[k]]-X[3*P],X[3*nb[k]+1]-X[3*P+1],X[3*nb[k]+2]-X[3*P+2]};
		lc(k,0)=Dot(dv,psi1);lc(k,1)=Dot(dv,psi2);}
	dArray2DT np(nn,1);np=support; dArrayT vol(nn);vol=arc*dz; dArrayT s(2);s[0]=s[1]=0.0;
	rkpm.SetField(lc,np,vol,s,2);
	const dArray2DT& Dp=rkpm.Dphi(); const dArray2DT& DDp=rkpm.DDphi();
	double x1[3]={0,0,0},x2[3]={0,0,0},x11[3]={0,0,0},x22[3]={0,0,0},x12[3]={0,0,0};
	for(int I=0;I<nn;I++){double Xq[3]={X[3*nb[I]],X[3*nb[I]+1],X[3*nb[I]+2]};
		for(int d=0;d<3;d++){x1[d]+=Dp(0,I)*Xq[d];x2[d]+=Dp(1,I)*Xq[d];
			x11[d]+=DDp(0,I)*Xq[d];x22[d]+=DDp(1,I)*Xq[d];x12[d]+=DDp(2,I)*Xq[d];}}
	ShellGeom G; BuildGeom(x1,x2,x11,x22,x12,0.5,0.0,G);

	/* impose v_I = w * (radial direction at node I) and accumulate grad(v3D) = sum B_I v_I */
	double L[3][3]={{0,0,0},{0,0,0},{0,0,0}};
	for(int I=0;I<nn;I++){
		double B[3][3][3]; BMatrix(G,Dp(0,I),Dp(1,I),DDp(0,I),DDp(2,I),DDp(1,I),B);
		double th=std::atan2(X[3*nb[I]+1],X[3*nb[I]]);
		double vI[3]={w*std::cos(th), w*std::sin(th), 0.0}; /* radial */
		for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int k=0;k<3;k++) L[i][j]+=B[i][j][k]*vI[k];
	}
	/* hoop strain = e_theta . sym(L) . e_theta, e_theta = circumferential tangent (psi that is hoop) */
	double D[3][3]; for(int i=0;i<3;i++)for(int j=0;j<3;j++) D[i][j]=0.5*(L[i][j]+L[j][i]);
	double thp=(nt/2)*dth; double et[3]={-std::sin(thp),std::cos(thp),0.0};
	double hoop=0; for(int i=0;i<3;i++)for(int j=0;j<3;j++) hoop+=et[i]*D[i][j]*et[j];
	double tr = D[0][0]+D[1][1]+D[2][2];
	printf("radial w=%.2f R=%.1f : hoop strain=%.5f (exact w/R=%.5f)  trace=%.5f\n",
		w, R, hoop, w/R, tr);
	EXPECT_NEAR(hoop, w/R, 0.1*w/R);
}

TEST(KLShellBenchmark, DISABLED_CurvedMembraneAudit)
{
	double exact = 1.0*10.0*10.0/(1.0e4*0.5); /* p R^2/(E h) = 0.02 */
	printf("Quarter-cylinder pressure (exact u_r=%.4f):\n", exact);
	for (int n = 9; n <= 29; n += 4)
		printf("  %dx%d : u_r=%.5f\n", n, n, QuarterCylinderPressure(n, n));
}

TEST(KLShellBenchmark, DISABLED_CylinderCurvatureAudit)
{
	for (double sf = 2.5; sf <= 4.01; sf += 0.5) {
		double kmax, kmin;
		CylinderCurvatureRKPM(25.0, 15, 15, sf, kmax, kmin);
		printf("R=25 supportFac=%.1f : kmax=%.5f (exact 0.04) kmin=%.5f (exact 0)\n", sf, kmax, kmin);
	}
}

TEST(KLShellBenchmark, DISABLED_FlatCantileverAudit)
{
	for (double h = 0.5; h >= 0.079; h *= 0.5) {
		double I = 8.0*h*h*h/12.0, wbeam = 512.0/(3.0*1.0e4*I);
		printf("FLAT h=%.3f (L/h=%.0f) beam=%.4f :", h, 8.0/h, wbeam);
		for (int n = 9; n <= 17; n += 4)
			printf("  %dx%d=%.4f", n, n, FlatCantilever(n, n, h, 3.0));
		printf("\n");
	}
}

/* Scordelis-Lo roof (Belytschko shell obstacle course; paper section 4.2.1). Reference
 * vertical deflection at the free-edge midpoint = 0.3006. REPRODUCED once two bugs were
 * fixed: (1) the plane-stress sigma33=0 condensation must act along the shell NORMAL, not the
 * global z-axis (ToVoigtLocal); (2) the diaphragm is a SIMPLE support (u_x=u_y=0, u_z free)
 * at both ends -- not a clamp. With both, the deflection converges to ~0.30 (within a few
 * percent) under refinement. */
TEST(KLShellBenchmark, ScordelisLoRoof)
{
	const double reference = 0.3006;

	double wCoarse = ScordelisLo(15, 15, 3.0);
	double wFine   = ScordelisLo(23, 23, 3.0);

	ASSERT_GT(wCoarse, 0.0);
	ASSERT_GT(wFine, 0.0);
	/* converges toward the reference under refinement, and the refined mesh is within
	 * engineering tolerance of 0.3006 */
	EXPECT_LT(std::fabs(wFine - reference), std::fabs(wCoarse - reference));
	EXPECT_LT(std::fabs(wFine - reference)/reference, 0.08);
}
