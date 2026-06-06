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
	double h, const double C[6][6], double support, std::vector<double>& K)
{
	int N = int(X.size())/3;
	int ndof = 3*N;
	K.assign((size_t)ndof*ndof, 0.0);

	double xg[3] = {-std::sqrt(3.0/5.0), 0.0, std::sqrt(3.0/5.0)};
	double wg[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0};

	D2OrthoMLS2DT efg(2);
	efg.Initialize();

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
		dArrayT dmax(nn); dmax = support;
		dArrayT sample(2); sample[0]=0.0; sample[1]=0.0;
		if (!efg.SetField(lc, dmax, sample)) continue;
		const dArray2DT& Dp = efg.Dphi();
		const dArray2DT& DDp = efg.DDphi();

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
			double cw = wg[g]*(h/2.0)*A_K;
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
					for (int ci=0;ci<3;ci++) for (int cj=0;cj<3;cj++) kij[ci][cj]+=Bv[I][a*3+ci]*Cab*Bv[J][b*3+cj];
				}
				int gi=3*nb[I], gj=3*nb[J];
				for (int ci=0;ci<3;ci++) for (int cj=0;cj<3;cj++) K[(size_t)(gi+ci)*ndof+(gj+cj)] += cw*kij[ci][cj];
			}
		}

		/* membrane stabilization (xi3=0) */
		ShellGeom G0;
		if (BuildGeom(x1,x2,x11,x22,x12,h,0.0,G0)) {
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
					double bv[6][3]; ToVoigt(B, bv);
					for (int r=0;r<6;r++) for (int c=0;c<3;c++) Bg[I][l*18+r*3+c]=bv[r][c];
				}
			}
			double Mw = V_K*Mmom;
			for (int l=0;l<2;l++) for (int I=0;I<nn;I++) for (int J=0;J<nn;J++) {
				double kij[3][3]={{0,0,0},{0,0,0},{0,0,0}};
				for (int a=0;a<6;a++) for (int b=0;b<6;b++) {
					double Cab=C[a][b]; if (Cab==0.0) continue;
					for (int ci=0;ci<3;ci++) for (int cj=0;cj<3;cj++) kij[ci][cj]+=Bg[I][l*18+a*3+ci]*Cab*Bg[J][l*18+b*3+cj];
				}
				int gi=3*nb[I], gj=3*nb[J];
				for (int ci=0;ci<3;ci++) for (int cj=0;cj<3;cj++) K[(size_t)(gi+ci)*ndof+(gj+cj)] += Mw*kij[ci][cj];
			}
		}
	}
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
int CurvedPatchZeroModes(int nt, int nz, std::vector<double>& smallest9)
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
	double support = 3.0*std::max(Rc*arc, dz);
	std::vector<double> K;
	AssembleShell(X, area, 0.05, C, support, K);
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

/* Scordelis-Lo roof: returns the downward deflection at the free-edge midpoint.
 * nt = nodes across the 80-degree arc, nz = nodes along the length. */
double ScordelisLo(int nt, int nz)
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
	double support = 3.0*std::max(arc, dz);
	std::vector<double> K;
	AssembleShell(X, area, h, C, support, K);
	int ndof = 3*N;

	/* boundary conditions:
	 *   rigid end diaphragms at z=0 and z=Lz: u_x = u_y = 0 (in the diaphragm plane);
	 *   remove the axial rigid mode: u_z = 0 at the z=0 diaphragm. */
	std::vector<char> fixed(ndof, 0);
	for (int i=0;i<nt;i++)
		for (int jr=0; jr<2; jr++) {          /* two rings near each end (rotation-free support) */
			int id0 = jr*nt+i;                /* z=0 end */
			int idL = (nz-1-jr)*nt+i;         /* z=Lz end */
			fixed[3*id0+0]=1; fixed[3*id0+1]=1; fixed[3*id0+2]=1;  /* diaphragm + axial restraint */
			fixed[3*idL+0]=1; fixed[3*idL+1]=1;                    /* diaphragm only */
		}

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

/* Diagnostic: does a FREE curved (cylinder) patch have exactly 6 zero-energy modes?
 * If > 6, the curved-surface nodal integration has a spurious interior hourglass
 * (basis/stabilization problem); if == 6, the Scordelis-Lo blow-up is a boundary/BC issue. */
TEST(KLShellBenchmark, DISABLED_CurvedPatchSpectrum)
{
	std::vector<double> s9a, s9b;
	int z1 = CurvedPatchZeroModes(7, 7, s9a);
	int z2 = CurvedPatchZeroModes(9, 9, s9b);
	printf("curved 7x7 zero modes: %d  (smallest9:", z1);
	for (double e : s9a) printf(" %.2e", e);
	printf(")\ncurved 9x9 zero modes: %d  (smallest9:", z2);
	for (double e : s9b) printf(" %.2e", e);
	printf(")\n");
	EXPECT_EQ(z1, 6);
	EXPECT_EQ(z2, 6);
}

/* DISABLED — work in progress. The geometry, diaphragm boundary conditions, gravity load
 * and curved-surface assembly are in place, but the solve currently exhibits a spurious
 * nodal-integration mode that grows under mesh refinement on the CURVED surface (the flat
 * cantilever in test_KLShellAssembly converges correctly, so curvature is the trigger).
 * The membrane-only Taylor-series stabilization, sufficient for the flat case, does not
 * control this curved-shell hourglass with the EFG orthogonal-MLS basis used here (Tahoe's
 * RKPM PolyBasis2DT is capped at linear completeness, issue #61). Resolving it likely needs
 * either an RKPM quadratic basis or an additional bending/curvature stabilization (#66).
 *
 * Run with:  ./test_KLShellBenchmark --gtest_also_run_disabled_tests */
TEST(KLShellBenchmark, DISABLED_ScordelisLoRoof)
{
	const double reference = 0.3006;

	double wCoarse = ScordelisLo(15, 15);
	double wFine   = ScordelisLo(21, 21);

	ASSERT_GT(wCoarse, 0.0);
	ASSERT_GT(wFine, 0.0);
	/* target behaviour (not yet met): monotonic convergence to the reference 0.3006 */
	EXPECT_LT(std::fabs(wFine - reference), std::fabs(wCoarse - reference));
	EXPECT_LT(std::fabs(wFine - reference)/reference, 0.20);
}
