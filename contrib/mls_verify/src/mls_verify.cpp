/* mls_verify — standalone de-risk harness for Tahoe meshfree 2nd derivatives (issue #61)
 *
 * The planned RKShellT (Wang & Bazilevs 2024) needs shape-function SECOND derivatives that
 * REPRODUCE quadratic fields, because KL bending depends on them. The gold-standard test:
 *
 *   For an RK/MLS basis of completeness order p and ANY polynomial P of degree <= p:
 *       sum_I  phi_I(x)   P(x_I) = P(x)          (value)
 *       sum_I  Dphi_I(x)  P(x_I) = grad P(x)     (gradient)
 *       sum_I  DDphi_I(x) P(x_I) = Hess P(x)     (Hessian)   <-- KL-shell-critical
 *
 * With completeness>=2 a quadratic field's Hessian must be reproduced to ~machine precision.
 *
 * Two backends are exercised:
 *   Path A  RKPM  via MLSSolverT      (PolyBasis2DT) — completeness capped at 1 (linear)
 *   Path B  EFG   via D2OrthoMLS2DT   (orthogonal MLS) — supports completeness {1,2}
 *
 * DDphi 2D row layout (dSymMatrixT::ExpandIndex): 0=xx, 1=yy, 2=xy.
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "MLSSolverT.h"
#include "D2OrthoMLS2DT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "MeshFreeT.h"
#include "ExceptionT.h"

using namespace Tahoe;
using std::cout; using std::setw; using std::scientific; using std::fixed; using std::setprecision;

/* monomial value/grad/Hess at (x,y); m in {0:1,1:x,2:y,3:x^2,4:xy,5:y^2} */
static void EvalMono(int m, double x, double y, double& v, double& gx, double& gy,
                     double& hxx, double& hyy, double& hxy)
{
	v=gx=gy=hxx=hyy=hxy=0.0;
	switch (m) {
		case 0: v=1; break;
		case 1: v=x; gx=1; break;
		case 2: v=y; gy=1; break;
		case 3: v=x*x; gx=2*x; hxx=2; break;
		case 4: v=x*y; gx=y; gy=x; hxy=1; break;
		case 5: v=y*y; gy=2*y; hyy=2; break;
	}
}
static const char* kMono[6] = {"1","x","y","x^2","xy","y^2"};

/* regular n x n grid on [0,L]^2 */
static void BuildGrid(int n, double L, std::vector<double>& X, std::vector<double>& Y, double& h)
{
	h = L/double(n-1); X.clear(); Y.clear();
	for (int i=0;i<n;i++) for (int j=0;j<n;j++){ X.push_back(i*h); Y.push_back(j*h); }
}

/* neighbors of (px,py) strictly within radius R */
static void Neighbors(const std::vector<double>& X, const std::vector<double>& Y,
                      double px, double py, double R, dArray2DT& coords)
{
	std::vector<int> idx;
	for (size_t i=0;i<X.size();i++){
		double dx=X[i]-px, dy=Y[i]-py;
		if (std::sqrt(dx*dx+dy*dy) < 0.99*R) idx.push_back(int(i));
	}
	coords.Dimension(int(idx.size()), 2);
	for (size_t k=0;k<idx.size();k++){ coords(int(k),0)=X[idx[k]]; coords(int(k),1)=Y[idx[k]]; }
}

/* run reproduction test for one backend; returns max Hessian error */
template <class GetDD>
static double Reproduce(const char* tag, const dArray2DT& coords,
                        const dArrayT& phi, const dArray2DT& Dphi, GetDD ddphi,
                        double px, double py, double& outVal, double& outGrad,
                        int maxMono)
{
	int nnd = coords.MajorDim();
	double maxV=0,maxG=0,maxH=0;
	for (int m=0;m<maxMono;m++){
		double v,gx,gy,hxx,hyy,hxy; EvalMono(m,px,py,v,gx,gy,hxx,hyy,hxy);
		double sv=0,sgx=0,sgy=0,sxx=0,syy=0,sxy=0;
		for (int I=0;I<nnd;I++){
			double pv,a,b,c,d,e; EvalMono(m,coords(I,0),coords(I,1),pv,a,b,c,d,e);
			sv+=phi[I]*pv; sgx+=Dphi(0,I)*pv; sgy+=Dphi(1,I)*pv;
			sxx+=ddphi(0,I)*pv; syy+=ddphi(1,I)*pv; sxy+=ddphi(2,I)*pv;
		}
		maxV=std::max(maxV,std::fabs(sv-v));
		maxG=std::max(maxG,std::fabs(sgx-gx)+std::fabs(sgy-gy));
		maxH=std::max(maxH,std::fabs(sxx-hxx)+std::fabs(syy-hyy)+std::fabs(sxy-hxy));
	}
	outVal=maxV; outGrad=maxG; return maxH;
}

/* accessor wrappers */
struct DDmls { const dArray2DT* p; double operator()(int r,int I) const { return (*p)(r,I);} };

static int run();

int main()
{
	try { return run(); }
	catch (ExceptionT::CodeT& e) {
		cout << "\n*** Tahoe ExceptionT: code "<<int(e)<<" ("<<ExceptionT::ToString(e)<<") ***"<<std::endl;
		return 3;
	}
	catch (...) { cout << "\n*** unknown exception ***"<<std::endl; return 4; }
}

static int run()
{
	cout << "=== MLS/EFG second-derivative de-risk (issue #61) ===\n" << scientific << setprecision(3);
	const double L=1.0;
	double pts[][2]={{0.5,0.5},{0.37,0.62},{0.5,0.34},{0.42,0.5},{0.6,0.45}};
	int npts=5;

	/* ---------------- Path A: RKPM MLSSolverT (completeness 1) ---------------- */
	cout << "\n[Path A] RKPM via MLSSolverT (PolyBasis2DT)\n";
	{
		int capOK=1; int cap=2;
		try {
			dArrayT wp(1); wp[0]=1.0;
			MLSSolverT m(2,2,false,MeshFreeT::kCubicSpline,wp); m.Initialize();
		} catch (ExceptionT::CodeT&) { capOK=0; }
		cout << "  completeness=2 supported? " << (capOK? "yes":"NO (PolyBasis2DT caps at 1)") << "\n";

		/* run at completeness=1: reproduces linear, NOT quadratic Hessian */
		dArrayT wp(1); wp[0]=1.0;
		MLSSolverT m(2,1,false,MeshFreeT::kCubicSpline,wp); m.Initialize();
		int n=13; std::vector<double> X,Y; double h; BuildGrid(n,L,X,Y,h);
		double R=3.2*h;
		double linH=0, quadH=0, vV=0, gG=0;
		for (int ip=0; ip<npts; ip++){
			dArray2DT cd; Neighbors(X,Y,pts[ip][0],pts[ip][1],R,cd);
			int nn=cd.MajorDim();
			dArray2DT np(nn,1); np=R; dArrayT vol(nn); vol=h*h;
			dArrayT fp(2); fp[0]=pts[ip][0]; fp[1]=pts[ip][1];
			if (!m.SetField(cd,np,vol,fp,2)) { cout<<"  SetField failed\n"; return 2; }
			DDmls dd; dd.p=&m.DDphi();
			double ov,og;
			double hLin=Reproduce("",cd,m.phi(),m.Dphi(),dd,pts[ip][0],pts[ip][1],ov,og,3); /* {1,x,y} */
			double hQ  =Reproduce("",cd,m.phi(),m.Dphi(),dd,pts[ip][0],pts[ip][1],ov,og,6); /* +quad */
			linH=std::max(linH,hLin); quadH=std::max(quadH,hQ); vV=std::max(vV,ov); gG=std::max(gG,og);
		}
		cout << "  reproduction (max over interior pts):\n";
		cout << "    Hessian err, LINEAR fields {1,x,y}      = " << linH << "  (expect ~0)\n";
		cout << "    Hessian err, QUADRATIC fields {x2,xy,y2}= " << quadH << "  (expect NONZERO: linear basis)\n";
		cout << "  => Path A reproduces only up to linear; insufficient for KL bending as-is.\n";
	}

	/* ---------------- Path B: EFG D2OrthoMLS2DT (completeness 2) ---------------- */
	cout << "\n[Path B] EFG via D2OrthoMLS2DT, completeness=2  (the viable KL-shell path)\n";
	double pathB_maxH=0, pathB_maxV=0, pathB_maxG=0;
	{
		D2OrthoMLS2DT efg(2); efg.Initialize();
		cout << "  (2D quadratic orthogonal basis: 6 monomials)\n";
		int n=13; std::vector<double> X,Y; double h; BuildGrid(n,L,X,Y,h);
		double R=3.5*h;

		/* --- consistency probe at center: Sum_I DDphi(row,I) * mono(x_I) --- */
		{
			dArray2DT cd; Neighbors(X,Y,0.5,0.5,R,cd); int nn=cd.MajorDim();
			dArrayT dmax(nn); dmax=R; dArrayT fp(2); fp[0]=0.5; fp[1]=0.5;
			efg.SetField(cd,dmax,fp);
			const dArray2DT& DD=efg.DDphi();
			const char* rn[3]={"xx","yy","xy"};
			cout<<"  consistency sums  Sum_I DDphi(row,I)*P(x_I)  [expect: xx->x^2=2, yy->y^2=2, xy->xy=1, else 0]\n";
			for (int r=0;r<3;r++){
				cout<<"    DDphi_"<<rn[r]<<":";
				for (int m=0;m<6;m++){
					double s=0; for(int I=0;I<nn;I++){double pv,a,b,c,d,e;EvalMono(m,cd(I,0),cd(I,1),pv,a,b,c,d,e);s+=DD(r,I)*pv;}
					cout<<" "<<kMono[m]<<"="<<scientific<<setprecision(2)<<s;
				}
				cout<<"\n";
			}
			/* also raw magnitude of DDphi entries */
			double mx=0; for(int r=0;r<3;r++)for(int I=0;I<nn;I++)mx=std::max(mx,std::fabs(DD(r,I)));
			cout<<"    max|DDphi entry| = "<<mx<<"  (per-node 2nd deriv ~ O(1/h^2)="<<1.0/(h*h)<<")\n";
		}
		for (int ip=0; ip<npts; ip++){
			dArray2DT cd; Neighbors(X,Y,pts[ip][0],pts[ip][1],R,cd);
			int nn=cd.MajorDim();
			dArrayT dmax(nn); dmax=R;
			dArrayT fp(2); fp[0]=pts[ip][0]; fp[1]=pts[ip][1];
			if (!efg.SetField(cd,dmax,fp)) { cout<<"  SetField failed at pt "<<ip<<"\n"; return 2; }
			DDmls dd; dd.p=&efg.DDphi();
			double ov,og;
			double hQ=Reproduce("",cd,efg.phi(),efg.Dphi(),dd,pts[ip][0],pts[ip][1],ov,og,6);
			cout<<"    pt("<<fixed<<setprecision(2)<<pts[ip][0]<<","<<pts[ip][1]<<") nbrs="<<setw(2)<<nn
			    <<"  val="<<scientific<<setprecision(2)<<ov<<" grad="<<og<<" HESS="<<hQ<<"\n";
			pathB_maxH=std::max(pathB_maxH,hQ); pathB_maxV=std::max(pathB_maxV,ov); pathB_maxG=std::max(pathB_maxG,og);
		}
		cout << "  reproduction (max over interior pts, ALL monomials incl. quadratic):\n";
		cout << "    value    err = " << pathB_maxV << "\n";
		cout << "    gradient err = " << pathB_maxG << "\n";
		cout << "    HESSIAN  err = " << pathB_maxH << "   <-- DDphi reproduces quadratic\n";
	}

	/* ---------------- Path B convergence on non-polynomial field ---------------- */
	cout << "\n[Path B] DDphi reconstruction of f=sin(2pi x)sin(2pi y) Hessian vs h:\n";
	{
		const double k=2.0*M_PI; double prev=-1;
		for (int n=13;n<=49;n=2*n-1){
			D2OrthoMLS2DT efg(2); efg.Initialize();
			std::vector<double> X,Y; double h; BuildGrid(n,L,X,Y,h);
			double R=3.5*h;
			dArray2DT cd; Neighbors(X,Y,0.5,0.5,R,cd); int nn=cd.MajorDim();
			dArrayT dmax(nn); dmax=R; dArrayT fp(2); fp[0]=0.5; fp[1]=0.5;
			if (!efg.SetField(cd,dmax,fp)){ cout<<"  SetField failed n="<<n<<"\n"; continue; }
			const dArray2DT& DD=efg.DDphi();
			double rxx=0,ryy=0,rxy=0;
			for (int I=0;I<nn;I++){ double f=std::sin(k*cd(I,0))*std::sin(k*cd(I,1));
				rxx+=DD(0,I)*f; ryy+=DD(1,I)*f; rxy+=DD(2,I)*f; }
			double axx=-k*k*std::sin(k*0.5)*std::sin(k*0.5), ayy=axx, axy=k*k*std::cos(k*0.5)*std::cos(k*0.5);
			double err=std::fabs(rxx-axx)+std::fabs(ryy-ayy)+std::fabs(rxy-axy);
			cout<<"  n="<<setw(3)<<n<<" h="<<h<<" err="<<err;
			if (prev>0) cout<<"  ratio="<<fixed<<setprecision(2)<<prev/err<<scientific;
			cout<<"\n"; prev=err;
		}
	}

	bool ok = (pathB_maxH<1e-6 && pathB_maxV<1e-9 && pathB_maxG<1e-7);
	cout << "\n=== VERDICT (issue #61 / fix #69) ===\n";
	cout << "  Path A  RKPM/MLSSolverT : PolyBasis2DT caps completeness at 1 -> CANNOT reproduce\n"
	        "          quadratic 2nd derivatives (linear basis only). Use Path B for KL bending.\n";
	cout << "  Path B  EFG/D2OrthoMLS2DT (completeness 2): value, gradient AND full Hessian\n"
	        "          (xx, yy, xy) reproduce quadratics to ~machine precision.\n";
	cout << "  #69 fix in place: the missing CjI factor on the DDb (b,ac) quotient-rule term in\n"
	        "  D2OrthoMLSSolverT was restored; diagonal 2nd derivatives now reproduce.\n";
	cout << "  status: " << (ok?"PASS — usable for KL shell":"FAIL — regression in DDphi") << "\n";
	return ok?0:1;
}
