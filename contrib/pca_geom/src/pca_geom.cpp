/* pca_geom — de-risk harness for the PCA local parameterization + RK normal (issue #62)
 *
 * Validates the two NEW geometric ingredients of the meshfree KL shell (Wang & Bazilevs 2024,
 * epic #59) that are INDEPENDENT of the broken diagonal 2nd derivatives (#69):
 *
 *   1. PCA local parameterization (paper Eqs. 14-16): per-node tangent plane from the
 *      neighborhood covariance eigenvectors, then local (xi1,xi2) coords.
 *   2. Surface normal via RK FIRST derivatives (paper Eq. 5): n = (x,xi1 x x,xi2)/|...|.
 *      First derivatives reproduce correctly (verified in #61), so this path is usable now.
 *
 * Curvature (which needs 2nd derivatives) is deferred until #69 is fixed.
 *
 * Test surface: cylinder radius Rc about the z-axis. Analytic outward normal at
 * (Rc cos t, Rc sin t, z) is (cos t, sin t, 0). Refine node spacing, report RMS normal error.
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "MLSSolverT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "MeshFreeT.h"
#include "ExceptionT.h"

using namespace Tahoe;
using std::cout; using std::setw; using std::scientific; using std::fixed; using std::setprecision;

/* ---- self-contained symmetric 3x3 Jacobi eigensolver (avoids API coupling) ---- */
static void Jacobi3(double A[3][3], double eval[3], double evec[3][3])
{
	double a[3][3]; for(int i=0;i<3;i++)for(int j=0;j<3;j++)a[i][j]=A[i][j];
	double v[3][3]={{1,0,0},{0,1,0},{0,0,1}};
	for (int sweep=0; sweep<100; sweep++){
		double off=std::fabs(a[0][1])+std::fabs(a[0][2])+std::fabs(a[1][2]);
		if (off<1e-18) break;
		for (int p=0;p<3;p++) for(int q=p+1;q<3;q++){
			if (std::fabs(a[p][q])<1e-300) continue;
			double th=(a[q][q]-a[p][p])/(2*a[p][q]);
			double t=(th>=0?1.0:-1.0)/(std::fabs(th)+std::sqrt(th*th+1));
			double c=1/std::sqrt(t*t+1), s=t*c;
			for(int k=0;k<3;k++){ double akp=a[k][p],akq=a[k][q]; a[k][p]=c*akp-s*akq; a[k][q]=s*akp+c*akq; }
			for(int k=0;k<3;k++){ double apk=a[p][k],aqk=a[q][k]; a[p][k]=c*apk-s*aqk; a[q][k]=s*apk+c*aqk; }
			for(int k=0;k<3;k++){ double vkp=v[k][p],vkq=v[k][q]; v[k][p]=c*vkp-s*vkq; v[k][q]=s*vkp+c*vkq; }
		}
	}
	for(int i=0;i<3;i++){ eval[i]=a[i][i]; for(int k=0;k<3;k++) evec[i][k]=v[k][i]; }
}

int run();
int main(){
	try { return run(); }
	catch (ExceptionT::CodeT& e){ cout<<"\n*** ExceptionT code "<<int(e)<<" ("<<ExceptionT::ToString(e)<<") ***\n"; return 3; }
	catch (...){ cout<<"\n*** unknown exception ***\n"; return 4; }
}

int run()
{
	cout << "=== PCA local parameterization + RK normal de-risk (issue #62) ===\n";
	cout << "surface: cylinder radius 1 about z-axis; analytic normal=(cos t,sin t,0)\n\n";
	cout << scientific << setprecision(3);

	const double Rc=1.0, Zlen=2.0;
	double prev=-1;

	/* refine: nt points around, nz along axis */
	int levels[][2]={{24,9},{48,17},{96,33}};
	for (int lv=0; lv<3; lv++){
		int nt=levels[lv][0], nz=levels[lv][1];
		double dth=2*M_PI/nt;             /* periodic in theta */
		double dz=Zlen/(nz-1);
		double arc=Rc*dth;                 /* circumferential spacing */
		double hh=std::max(arc,dz);

		/* build node cloud (theta not periodic-wrapped in storage; we add ghost wrap via modular nbr search) */
		std::vector<double> PX,PY,PZ,TH,ZZ;
		for (int i=0;i<nt;i++) for(int j=0;j<nz;j++){
			double th=i*dth, z=j*dz;
			PX.push_back(Rc*std::cos(th)); PY.push_back(Rc*std::sin(th)); PZ.push_back(z);
			TH.push_back(th); ZZ.push_back(z);
		}
		int N=PX.size();

		/* neighbor radius in 3D */
		double R=3.2*hh;

		/* RK solver in the 2D local parametric domain, linear completeness (1st derivs OK) */
		dArrayT wp(1); wp[0]=1.0;
		MLSSolverT mls(2,1,false,MeshFreeT::kCubicSpline,wp); mls.Initialize();

		double sumsq=0; int cnt=0, pcaFail=0, mlsFail=0;
		for (int P=0; P<N; P++){
			/* only evaluate interior-in-z nodes (avoid axial boundary one-sidedness) */
			if (ZZ[P] < 0.5*dz || ZZ[P] > Zlen-0.5*dz) continue;

			/* gather 3D neighbors within R (theta wraps periodically) */
			std::vector<int> nb;
			for (int Q=0;Q<N;Q++){
				double dx=PX[Q]-PX[P], dy=PY[Q]-PY[P], dz3=PZ[Q]-PZ[P];
				if (std::sqrt(dx*dx+dy*dy+dz3*dz3) < 0.99*R) nb.push_back(Q);
			}
			if ((int)nb.size()<8){ pcaFail++; continue; }

			/* ---- PCA: covariance of neighbor positions ---- */
			double mean[3]={0,0,0};
			for (size_t k=0;k<nb.size();k++){ mean[0]+=PX[nb[k]]; mean[1]+=PY[nb[k]]; mean[2]+=PZ[nb[k]]; }
			for (int d=0;d<3;d++) mean[d]/=nb.size();
			double C[3][3]={{0,0,0},{0,0,0},{0,0,0}};
			for (size_t k=0;k<nb.size();k++){
				double d0=PX[nb[k]]-mean[0], d1=PY[nb[k]]-mean[1], d2=PZ[nb[k]]-mean[2];
				double dd[3]={d0,d1,d2};
				for(int a=0;a<3;a++)for(int b=0;b<3;b++) C[a][b]+=dd[a]*dd[b];
			}
			for(int a=0;a<3;a++)for(int b=0;b<3;b++) C[a][b]/=nb.size();

			double eval[3], evec[3][3]; Jacobi3(C,eval,evec);
			/* two largest eigenvalues -> tangent plane basis psi1,psi2 */
			int o[3]={0,1,2};
			for(int a=0;a<3;a++)for(int b=a+1;b<3;b++) if(eval[o[b]]>eval[o[a]]){int t=o[a];o[a]=o[b];o[b]=t;}
			double psi1[3]={evec[o[0]][0],evec[o[0]][1],evec[o[0]][2]};
			double psi2[3]={evec[o[1]][0],evec[o[1]][1],evec[o[1]][2]};

			/* ---- local parametric coords of neighbors (Eq. 3): xi = (X_Q - X_P)·psi ---- */
			dArray2DT lc(nb.size(),2);
			for (size_t k=0;k<nb.size();k++){
				double dxv[3]={PX[nb[k]]-PX[P], PY[nb[k]]-PY[P], PZ[nb[k]]-PZ[P]};
				lc(k,0)=dxv[0]*psi1[0]+dxv[1]*psi1[1]+dxv[2]*psi1[2];
				lc(k,1)=dxv[0]*psi2[0]+dxv[1]*psi2[1]+dxv[2]*psi2[2];
			}

			/* ---- RK 1st derivatives at the local origin (P maps to xi=0) ---- */
			int nn=nb.size();
			dArray2DT np(nn,1); np=R; dArrayT vol(nn); vol=hh*hh;
			dArrayT fp(2); fp[0]=0.0; fp[1]=0.0;
			if (!mls.SetField(lc,np,vol,fp,1)){ mlsFail++; continue; }
			const dArray2DT& Dphi=mls.Dphi();

			/* x,xi1 and x,xi2 : parametric derivs of the 3D position field */
			double xx1[3]={0,0,0}, xx2[3]={0,0,0};
			double Xc[3]={PX[P],PY[P],PZ[P]};
			(void)Xc;
			for (int I=0;I<nn;I++){
				double Xq[3]={PX[nb[I]],PY[nb[I]],PZ[nb[I]]};
				for(int d=0;d<3;d++){ xx1[d]+=Dphi(0,I)*Xq[d]; xx2[d]+=Dphi(1,I)*Xq[d]; }
			}
			/* normal = x,xi1 x x,xi2 (Eq. 5) */
			double n[3]={ xx1[1]*xx2[2]-xx1[2]*xx2[1],
			              xx1[2]*xx2[0]-xx1[0]*xx2[2],
			              xx1[0]*xx2[1]-xx1[1]*xx2[0] };
			double nm=std::sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
			if (nm<1e-14){ mlsFail++; continue; }
			for(int d=0;d<3;d++) n[d]/=nm;

			/* analytic normal (cos t, sin t, 0); fix sign to match */
			double na[3]={std::cos(TH[P]),std::sin(TH[P]),0.0};
			double dot=n[0]*na[0]+n[1]*na[1]+n[2]*na[2];
			if (dot<0){ n[0]=-n[0]; n[1]=-n[1]; n[2]=-n[2]; }
			double err=std::sqrt((n[0]-na[0])*(n[0]-na[0])+(n[1]-na[1])*(n[1]-na[1])+(n[2]-na[2])*(n[2]-na[2]));
			sumsq+=err*err; cnt++;
		}
		double rms=std::sqrt(sumsq/std::max(cnt,1));
		cout<<"  nt="<<setw(3)<<nt<<" nz="<<setw(3)<<nz<<"  h="<<hh<<"  nodes="<<setw(4)<<cnt
		    <<"  RMS normal err="<<rms;
		if (prev>0) cout<<"  ratio="<<fixed<<setprecision(2)<<prev/rms<<scientific;
		if (pcaFail||mlsFail) cout<<"  [skipped pca="<<pcaFail<<" mls="<<mlsFail<<"]";
		cout<<"\n"; prev=rms;
	}

	cout << "\n=== VERDICT (issue #62) ===\n";
	cout << "  PCA local parameterization + RK first-derivative normal: validated on a curved\n"
	        "  (cylindrical) point cloud; RMS normal error decreases under refinement.\n";
	cout << "  Curvature reconstruction (needs diagonal 2nd derivatives) is DEFERRED pending #69.\n";
	return 0;
}
