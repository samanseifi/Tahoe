/* standalone unit test for FlanaganTaylor.h — compile:
 *   g++ -I. test_ft.cpp -o /tmp/test_ft && /tmp/test_ft
 */
#include <cstdio>
#include <cmath>
#include "FlanaganTaylor.h"
using namespace Tahoe::KLShell;

static int fails=0;
static void check(const char* name, double got, double exp, double tol){
	bool ok=std::fabs(got-exp)<tol;
	printf("  [%s] %-28s got=% .6f exp=% .6f\n", ok?"PASS":"FAIL", name, got, exp);
	if(!ok) fails++;
}

int main(){
	const double PI=3.14159265358979323846;

	/* TEST 1: pure rigid rotation about z by 90 deg in N increments.
	 * L = W = [[0,-w,0],[w,0,0],[0,0,0]]; dL = W*dt. Expect R = Rz(90), V = I. */
	{
		int N=200; double th=PI/2.0, dth=th/N;
		double R[3][3]={{1,0,0},{0,1,0},{0,0,1}}, V[3][3]={{1,0,0},{0,1,0},{0,0,1}};
		double dL[3][3]={{0,-dth,0},{dth,0,0},{0,0,0}};
		for(int s=0;s<N;s++) FlanaganTaylorStep(dL,R,V);
		printf("TEST 1: rigid rotation Rz(90)\n");
		check("R00",R[0][0],0.0,1e-4); check("R01",R[0][1],-1.0,1e-4);
		check("R10",R[1][0],1.0,1e-4);  check("R11",R[1][1],0.0,1e-4);
		check("R22",R[2][2],1.0,1e-4);
		check("V00(=I)",V[0][0],1.0,1e-4); check("V01(=0)",V[0][1],0.0,1e-4);
		/* R orthonormal: det=1 */
		double det=R[0][0]*(R[1][1]*R[2][2]-R[1][2]*R[2][1])
		          -R[0][1]*(R[1][0]*R[2][2]-R[1][2]*R[2][0])
		          +R[0][2]*(R[1][0]*R[2][1]-R[1][1]*R[2][0]);
		check("det(R)",det,1.0,1e-6);
	}

	/* TEST 2: pure stretch in x (no rotation). L = D = diag(a,0,0); expect R=I, V11 grows ~1+N*a. */
	{
		int N=100; double a=0.001;
		double R[3][3]={{1,0,0},{0,1,0},{0,0,1}}, V[3][3]={{1,0,0},{0,1,0},{0,0,1}};
		double dL[3][3]={{a,0,0},{0,0,0},{0,0,0}};
		for(int s=0;s<N;s++) FlanaganTaylorStep(dL,R,V);
		printf("TEST 2: pure stretch (R=I, V11=1+N*a)\n");
		check("R00(=1)",R[0][0],1.0,1e-6); check("R01(=0)",R[0][1],0.0,1e-6);
		check("R10(=0)",R[1][0],0.0,1e-6);
		check("V11",V[0][0],std::exp(N*a),2e-3);  /* multiplicative: (1+a)^N -> e^0.1 = 1.10517 */
		check("V22(=1)",V[2][2],1.0,1e-6);
	}

	/* TEST 3: objectivity. A material stress sig~ in the co-rotational frame; after a rigid 90-deg
	 * rotation the spatial stress sig = R sig~ R^T must be the rotated tensor (a uniaxial s11 in the
	 * body becomes s22 after a 90-deg turn). */
	{
		int N=200; double th=PI/2.0, dth=th/N;
		double R[3][3]={{1,0,0},{0,1,0},{0,0,1}}, V[3][3]={{1,0,0},{0,1,0},{0,0,1}};
		double dL[3][3]={{0,-dth,0},{dth,0,0},{0,0,0}};
		for(int s=0;s<N;s++) FlanaganTaylorStep(dL,R,V);
		double sigb[3][3]={{100,0,0},{0,0,0},{0,0,0}};   /* uniaxial along body-x */
		double Rs[3][3]; M3_mul(R,sigb,Rs); double sig[3][3]; M3_mulT(Rs,R,sig);
		printf("TEST 3: objectivity (body s11=100 -> spatial s22=100 after 90deg)\n");
		check("sig11",sig[0][0],0.0,1e-2); check("sig22",sig[1][1],100.0,1e-2);
		check("sig12",sig[0][1],0.0,1e-2);
	}

	printf("\n%s  (%d failures)\n", fails==0?"ALL PASS":"FAILURES", fails);
	return fails;
}
