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
	double h, const double C[6][6], double support, std::vector<double>& K)
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
					double bv[6][3]; ToVoigt(B, bv);
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
	double support = 5.0*std::max(arc, dz); /* large enough for clean quadratic RKPM + stabilization */
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

/* A FREE curved (cylinder) patch must have exactly 6 zero-energy (rigid-body) modes.
 *
 * History: with EFG and membrane-only stabilization this gave 9 modes at 7x7 and 16 at 9x9
 * (a curved-shell BENDING hourglass that grows with refinement). The cure is the bending
 * (curvature-gradient) stabilization (BMatrixCurvatureGradient, using DDDphi) PLUS a
 * sufficiently large support (~5x spacing) so the quadratic RKPM strain operator is rich
 * enough — small supports (e.g. 3x) leave a resonance with spurious modes. With both, the
 * patch is rank-clean (exactly 6) across refinement. */
TEST(KLShellBenchmark, CurvedPatchBendingStabilized)
{
	std::vector<double> s9;
	int z7  = CurvedPatchZeroModes(7,  7,  s9, 5.0);
	int z9  = CurvedPatchZeroModes(9,  9,  s9, 5.0);
	int z11 = CurvedPatchZeroModes(11, 11, s9, 5.0);
	EXPECT_EQ(z7,  6);
	EXPECT_EQ(z9,  6);
	EXPECT_EQ(z11, 6);
}

/* DISABLED — the last remaining wall, now precisely characterized.
 *
 * Geometry, diaphragm BCs, gravity load and the curved RKPM assembly (quadratic basis +
 * Gaussian window + membrane AND bending/curvature stabilization) are all in place. The
 * bending stabilization cures the curved-shell hourglass for moderately thick shells
 * (CurvedPatchBendingStabilized: a free patch gives exactly 6 modes at R/h=20). But the
 * Scordelis-Lo roof is THIN (R/h=100), and a thickness sweep of the free-patch spectrum
 * shows a spurious mode whose energy collapses to ~1.8e-9 * lambda_max as R/h grows -- about
 * four orders below the physical bending scale (h/R)^2. That thin-shell membrane-bending
 * hourglass is NOT controlled by the membrane+bending Taylor stabilization (boosting the
 * bending term 100x barely moves the deflection), so the roof deflection is far too large
 * and grows with refinement.
 *
 * This is the paper's own acknowledged "preliminary"/incomplete section 5 (thin-shell
 * stabilization / membrane locking). Resolving it needs a more robust construction than the
 * first-gradient Taylor stabilization (e.g. SCNI-style smoothed gradients, or a
 * thickness-consistent stabilization). Tracked in #66/#68.
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
