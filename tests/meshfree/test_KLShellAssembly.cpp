/* test_KLShellAssembly.cpp — meshfree KL-shell assembly & boundary-value-problem tests.
 *
 * Assembles the shell tangent stiffness K = sum_K int B^T C B (nodal integration +
 * Taylor-series membrane stabilization + 3-pt through-thickness Gauss + sigma33=0
 * plane-stress condensation) using KLShellKernels, and checks:
 *
 *   - rigid-body consistency: K . r = 0 for the 6 rigid-body modes;
 *   - rank sufficiency: the free patch has exactly 6 zero-energy modes (no spurious modes);
 *   - cantilever bending: with the rotation-free two-row clamp, the tip deflection converges
 *     monotonically to Euler-Bernoulli beam theory under mesh refinement.
 */

#include "gtest/gtest.h"

#include <cmath>
#include <vector>
#include <algorithm>

#include "D2OrthoMLS2DT.h"
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

/* Voigt 6x3 of B[3][3][3] -> flat 18 */
void VoigtFlat(const double B[3][3][3], double out[18])
{
	double bv[6][3];
	toVoigt(B, bv);
	for (int r=0;r<6;r++) for (int c=0;c<3;c++) out[r*3+c] = bv[r][c];
}

/* assemble the dense stiffness for a flat nx x ny cloud on [0,Lx] x [-Wy/2,Wy/2] */
void Assemble(int nx, int ny, double Lx, double Wy, double h, const double C[6][6],
	std::vector<double>& X, std::vector<double>& Y, std::vector<double>& K)
{
	double sx = Lx/(nx-1), sy = Wy/(ny-1), s = std::max(sx, sy);
	int N = nx*ny;
	X.assign(N, 0.0); Y.assign(N, 0.0);
	for (int j=0;j<ny;j++) for (int i=0;i<nx;i++) { int id=j*nx+i; X[id]=i*sx; Y[id]=-Wy/2+j*sy; }

	int ndof = 3*N;
	K.assign(ndof*ndof, 0.0);
	double R = 2.6*s;
	double xg[3] = {-std::sqrt(3.0/5.0), 0.0, std::sqrt(3.0/5.0)};
	double wg[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0};
	double A_K = sx*sy, V_K = A_K*h, Mx = sx*sx/12.0, My = sy*sy/12.0;

	D2OrthoMLS2DT efg(2);
	efg.Initialize();

	for (int Kn=0;Kn<N;Kn++) {

		std::vector<int> nb;
		for (int Q=0;Q<N;Q++) {
			double dx=X[Q]-X[Kn], dy=Y[Q]-Y[Kn];
			if (std::sqrt(dx*dx+dy*dy) < 0.99*R) nb.push_back(Q);
		}
		int nn = int(nb.size());
		if (nn < 6) continue;

		dArray2DT lc(nn, 2);
		for (int k=0;k<nn;k++) { lc(k,0)=X[nb[k]]-X[Kn]; lc(k,1)=Y[nb[k]]-Y[Kn]; }
		dArrayT dmax(nn); dmax = R;
		dArrayT sample(2); sample[0]=0.0; sample[1]=0.0;
		if (!efg.SetField(lc, dmax, sample)) continue;
		const dArray2DT& Dp = efg.Dphi();
		const dArray2DT& DDp = efg.DDphi();

		double x1[3]={0,0,0},x2[3]={0,0,0},x11[3]={0,0,0},x22[3]={0,0,0},x12[3]={0,0,0};
		for (int I=0;I<nn;I++) {
			double Xq[3]={X[nb[I]],Y[nb[I]],0.0};
			for (int d=0;d<3;d++) {
				x1[d]+=Dp(0,I)*Xq[d]; x2[d]+=Dp(1,I)*Xq[d];
				x11[d]+=DDp(0,I)*Xq[d]; x22[d]+=DDp(1,I)*Xq[d]; x12[d]+=DDp(2,I)*Xq[d];
			}
		}

		/* nodal-integration term: 3-pt Gauss through thickness */
		for (int g=0;g<3;g++) {
			ShellGeom G;
			if (!buildGeom(x1,x2,x11,x22,x12,h,xg[g],G)) continue;
			double cw = wg[g]*(h/2.0)*A_K;
			std::vector<std::vector<double> > Bv(nn, std::vector<double>(18));
			for (int I=0;I<nn;I++) {
				double B[3][3][3];
				Bmatrix(G, Dp(0,I),Dp(1,I),DDp(0,I),DDp(2,I),DDp(1,I), B);
				VoigtFlat(B, &Bv[I][0]);
			}
			for (int I=0;I<nn;I++) for (int J=0;J<nn;J++) {
				double kij[3][3]={{0,0,0},{0,0,0},{0,0,0}};
				for (int a=0;a<6;a++) for (int b=0;b<6;b++) {
					double Cab = C[a][b];
					if (Cab == 0.0) continue;
					for (int ci=0;ci<3;ci++) for (int cj=0;cj<3;cj++)
						kij[ci][cj] += Bv[I][a*3+ci]*Cab*Bv[J][b*3+cj];
				}
				int gi=3*nb[I], gj=3*nb[J];
				for (int ci=0;ci<3;ci++) for (int cj=0;cj<3;cj++)
					K[(gi+ci)*ndof+(gj+cj)] += cw*kij[ci][cj];
			}
		}

		/* membrane stabilization (xi3=0) */
		ShellGeom G0;
		if (buildGeom(x1,x2,x11,x22,x12,h,0.0,G0)) {
			std::vector<std::vector<double> > Bg(nn, std::vector<double>(36));
			for (int I=0;I<nn;I++) {
				double Bz[3][3][3];
				Bmatrix(G0, Dp(0,I),Dp(1,I),DDp(0,I),DDp(2,I),DDp(1,I), Bz);
				double P1l[2]={DDp(0,I),DDp(2,I)}, P2l[2]={DDp(2,I),DDp(1,I)};
				double Bgr[3][3][3][2];
				BmatrixGrad(G0, Dp(0,I),Dp(1,I),DDp(0,I),DDp(2,I),DDp(1,I), P1l,P2l,Bz,Bgr);
				for (int l=0;l<2;l++) {
					double B[3][3][3];
					for (int i=0;i<3;i++) for (int j=0;j<3;j++) for (int k=0;k<3;k++) B[i][j][k]=Bgr[i][j][k][l];
					VoigtFlat(B, &Bg[I][l*18]);
				}
			}
			double Mw[2] = {V_K*Mx, V_K*My};
			for (int l=0;l<2;l++) for (int I=0;I<nn;I++) for (int J=0;J<nn;J++) {
				double kij[3][3]={{0,0,0},{0,0,0},{0,0,0}};
				for (int a=0;a<6;a++) for (int b=0;b<6;b++) {
					double Cab = C[a][b];
					if (Cab == 0.0) continue;
					for (int ci=0;ci<3;ci++) for (int cj=0;cj<3;cj++)
						kij[ci][cj] += Bg[I][l*18+a*3+ci]*Cab*Bg[J][l*18+b*3+cj];
				}
				int gi=3*nb[I], gj=3*nb[J];
				for (int ci=0;ci<3;ci++) for (int cj=0;cj<3;cj++)
					K[(gi+ci)*ndof+(gj+cj)] += Mw[l]*kij[ci][cj];
			}
		}
	}
	for (int i=0;i<ndof;i++) for (int j=i+1;j<ndof;j++) {
		double a = 0.5*(K[i*ndof+j]+K[j*ndof+i]);
		K[i*ndof+j] = K[j*ndof+i] = a;
	}
}

/* cyclic Jacobi eigenvalues (ascending) of a symmetric n x n matrix */
void JacobiEig(std::vector<double> A, int n, std::vector<double>& ev)
{
	for (int sweep=0;sweep<80;sweep++) {
		double off = 0.0;
		for (int p=0;p<n;p++) for (int q=p+1;q<n;q++) off += A[p*n+q]*A[p*n+q];
		if (off < 1e-26) break;
		for (int p=0;p<n;p++) for (int q=p+1;q<n;q++) {
			double apq = A[p*n+q];
			if (std::fabs(apq) < 1e-300) continue;
			double th = (A[q*n+q]-A[p*n+p])/(2*apq);
			double t = (th>=0 ? 1.0 : -1.0)/(std::fabs(th)+std::sqrt(th*th+1));
			double c = 1/std::sqrt(t*t+1), sn = t*c;
			for (int k=0;k<n;k++) { double a1=A[k*n+p],a2=A[k*n+q]; A[k*n+p]=c*a1-sn*a2; A[k*n+q]=sn*a1+c*a2; }
			for (int k=0;k<n;k++) { double a1=A[p*n+k],a2=A[q*n+k]; A[p*n+k]=c*a1-sn*a2; A[q*n+k]=sn*a1+c*a2; }
		}
	}
	ev.resize(n);
	for (int i=0;i<n;i++) ev[i] = A[i*n+i];
	std::sort(ev.begin(), ev.end());
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

double CantileverTipDeflection(int nx, int ny, double Lx, double Wy, double h, double E, double nu, double P)
{
	double C[6][6];
	PlaneStressTangent(E, nu, C);
	std::vector<double> X, Y, K;
	Assemble(nx, ny, Lx, Wy, h, C, X, Y, K);
	int N = nx*ny, ndof = 3*N;
	double sx = Lx/(nx-1);

	/* clamp TWO rows at x=0 (rotation-free shell); transverse tip load at x=Lx */
	std::vector<char> fixed(ndof, 0);
	for (int p=0;p<N;p++) if (X[p] < 1.5*sx) for (int d=0;d<3;d++) fixed[3*p+d]=1;
	int ntip = 0;
	for (int p=0;p<N;p++) if (std::fabs(X[p]-Lx) < 1e-9) ntip++;
	std::vector<double> f(ndof, 0.0);
	for (int p=0;p<N;p++) if (std::fabs(X[p]-Lx) < 1e-9) f[3*p+2] += P/ntip;

	std::vector<int> map(ndof, -1);
	int nf = 0;
	for (int i=0;i<ndof;i++) if (!fixed[i]) map[i]=nf++;
	std::vector<double> Kr(nf*nf, 0.0), fr(nf, 0.0);
	for (int i=0;i<ndof;i++) {
		if (fixed[i]) continue;
		fr[map[i]] = f[i];
		for (int j=0;j<ndof;j++) if (!fixed[j]) Kr[map[i]*nf+map[j]] = K[i*ndof+j];
	}
	std::vector<double> ur;
	if (!SolveDense(Kr, fr, nf, ur)) return -1.0;

	double wtip = 0.0;
	for (int p=0;p<N;p++) if (std::fabs(X[p]-Lx) < 1e-9) wtip += ur[map[3*p+2]];
	return std::fabs(wtip / ntip);
}

} /* anonymous namespace */

TEST(KLShellAssembly, RigidBodyModesAndRank)
{
	const int n = 7;
	const double L = 1.0, h = 0.05, E = 100.0, nu = 0.3;
	double C[6][6];
	PlaneStressTangent(E, nu, C);
	std::vector<double> X, Y, K;
	Assemble(n, n, L, L, h, C, X, Y, K);
	int N = n*n, ndof = 3*N;

	/* eigenvalue spectrum: exactly 6 zero (rigid) modes, no spurious zero modes */
	std::vector<double> ev;
	JacobiEig(K, ndof, ev);
	double lmax = ev.back();
	int nzero = 0;
	for (int i=0;i<ndof;i++) if (ev[i] < 1e-8*lmax) nzero++;
	EXPECT_EQ(nzero, 6);

	/* rigid-body consistency: K . r = 0 for translations and rotations */
	double maxResid = 0.0;
	for (int mode=0; mode<6; mode++) {
		std::vector<double> u(ndof, 0.0);
		for (int p=0;p<N;p++) {
			double Xp[3] = {X[p], Y[p], 0.0};
			double uu[3] = {0,0,0};
			if (mode < 3) uu[mode] = 1.0;
			else {
				double w[3] = {0,0,0};
				w[mode-3] = 1.0;
				double r[3] = {Xp[0]-0.5, Xp[1]-0.5, Xp[2]};
				cross3(w, r, uu);
			}
			for (int d=0;d<3;d++) u[3*p+d] = uu[d];
		}
		double un = 0.0;
		for (int i=0;i<ndof;i++) un = std::max(un, std::fabs(u[i]));
		for (int i=0;i<ndof;i++) {
			double s = 0.0;
			for (int j=0;j<ndof;j++) s += K[i*ndof+j]*u[j];
			maxResid = std::max(maxResid, std::fabs(s)/un);
		}
	}
	EXPECT_LT(maxResid, 1e-7*lmax);
}

TEST(KLShellAssembly, CantileverConvergesToBeamTheory)
{
	const double Lx = 2.0, Wy = 0.2, h = 0.05, E = 1.0e6, nu = 0.3, P = 0.5;
	double Imom = Wy*h*h*h/12.0;
	double dbeam = P*Lx*Lx*Lx/(3.0*E*Imom);

	double wCoarse = CantileverTipDeflection(17, 5, Lx, Wy, h, E, nu, P);
	double wFine   = CantileverTipDeflection(25, 7, Lx, Wy, h, E, nu, P);

	ASSERT_GT(wCoarse, 0.0);
	ASSERT_GT(wFine, 0.0);

	/* finite, sensible, and converging toward beam theory under refinement */
	EXPECT_LT(std::fabs(wFine - dbeam), std::fabs(wCoarse - dbeam));
	EXPECT_LT(std::fabs(wFine - dbeam)/dbeam, 0.30);
}
