/* BonetTetT.cpp */
#include "BonetTetT.h"

#include "ANPHelperT.h"
#include "ElementSupportT.h"
#include "ModelManagerT.h"
#include "FiniteStrainT.h"
#include "ElementCardT.h"
#include "iArrayT.h"
#include "dArray2DT.h"
#include "dMatrixT.h"

#include <iostream>
#include <cmath>

using namespace Tahoe;

BonetTetT::BonetTetT(const ElementSupportT& support)
	: UpdatedLagrangianT(support),
	  fANP(NULL),
	  fFlatConn(NULL),
	  fVrefE(NULL),
	  fJe(NULL),
	  fJbarE(NULL),
	  fNelem(0)
{
	SetName("bonet_tet");
}

BonetTetT::~BonetTetT(void)
{
	delete fANP;
	delete[] fFlatConn;
	delete[] fVrefE;
	delete[] fJe;
	delete[] fJbarE;
}

void BonetTetT::TakeParameterList(const ParameterListT& list)
{
	/* base class — sets up connectivity, shape funcs, materials, etc. */
	UpdatedLagrangianT::TakeParameterList(list);

	/* this element only makes sense for Tet4 */
	if (NumElementNodes() != 4 || NumSD() != 3) {
		std::cout << "BonetTetT: WARNING — only valid for 3D Tet4 (nen=4, nsd=3); "
		          << "got nen=" << NumElementNodes() << " nsd=" << NumSD()
		          << ".  ANP F-bar disabled." << std::endl;
		return;
	}

	BuildANPData();
	std::cout << "BonetTetT: ANP-Tet4 (LS-DYNA ELFORM=13) active for "
	          << fNelem << " tets" << std::endl;
}

/* Build flat connectivity, reference volumes, and instantiate ANP helper. */
void BonetTetT::BuildANPData(void)
{
	/* count total elements */
	fNelem = 0;
	for (int b = 0; b < fBlockData.Length(); b++)
		fNelem += fBlockData[b].Dimension();

	const int nen = 4;
	delete[] fFlatConn; delete[] fVrefE; delete[] fJe; delete[] fJbarE;
	fFlatConn = new int[fNelem * nen];
	fVrefE    = new double[fNelem];
	fJe       = new double[fNelem]();
	fJbarE    = new double[fNelem]();

	const dArray2DT& ref = ElementSupport().InitialCoordinates();
	for (int e = 0; e < fNelem; e++) {
		const ElementCardT& card = ElementCard(e);
		const iArrayT& nodes = card.NodesX();
		for (int n = 0; n < nen; n++)
			fFlatConn[e * nen + n] = nodes[n];

		double x[4], y[4], z[4];
		for (int n = 0; n < nen; n++) {
			int gn = fFlatConn[e * nen + n];
			x[n] = ref(gn, 0); y[n] = ref(gn, 1); z[n] = ref(gn, 2);
		}
		double a1=x[0]-x[2], a2=y[0]-y[2], a3=z[0]-z[2];
		double b1=x[1]-x[2], b2=y[1]-y[2], b3=z[1]-z[2];
		double c1=x[3]-x[2], c2=y[3]-y[2], c3=z[3]-z[2];
		fVrefE[e] = std::fabs(a1*(b2*c3-b3*c2)
		                    - a2*(b1*c3-b3*c1)
		                    + a3*(b1*c2-b2*c1)) / 6.0;
	}

	int numnod = ref.MajorDim();
	delete fANP;
	fANP = new ANPHelperT();
	fANP->Init(fNelem, numnod, nen, fFlatConn, fVrefE);
}

/* Compute J_e = V_curr / V_ref for every element. */
void BonetTetT::ComputeAllJe(void)
{
	const int nen = 4;
	const dArray2DT& cur = ElementSupport().CurrentCoordinates();
	#pragma omp parallel for if(fNelem > 1024)
	for (int e = 0; e < fNelem; e++) {
		const int* ec = fFlatConn + e * nen;
		double x[4], y[4], z[4];
		for (int n = 0; n < nen; n++) {
			x[n] = cur(ec[n], 0); y[n] = cur(ec[n], 1); z[n] = cur(ec[n], 2);
		}
		double a1=x[0]-x[2], a2=y[0]-y[2], a3=z[0]-z[2];
		double b1=x[1]-x[2], b2=y[1]-y[2], b3=z[1]-z[2];
		double c1=x[3]-x[2], c2=y[3]-y[2], c3=z[3]-z[2];
		double V = std::fabs(a1*(b2*c3-b3*c2)
		                   - a2*(b1*c3-b3*c1)
		                   + a3*(b1*c2-b2*c1)) / 6.0;
		fJe[e] = (fVrefE[e] > 0.0) ? V / fVrefE[e] : 1.0;
	}
}

/* Override SetGlobalShape — call parent first, then apply F-bar in-place. */
void BonetTetT::SetGlobalShape(void)
{
	UpdatedLagrangianT::SetGlobalShape();
	if (!fANP || NumElementNodes() != 4) return;

	int e = CurrElementNumber();
	if (e < 0 || e >= fNelem) return;

	double Jb = fJbarE[e];
	int nip = fF_List.Length();
	for (int ip = 0; ip < nip; ip++) {
		dMatrixT& F = fF_List[ip];
		double J = F.Det();
		if (J <= 0.0) continue;
		double scale = std::cbrt(Jb / J);
		double* p = F.Pointer();
		for (int k = 0; k < 9; k++) p[k] *= scale;
	}
}

/* RHS pre-pass: compute J_e for every element, run the ANP helper, then
 * dispatch to the inherited element loop (which calls our SetGlobalShape). */
void BonetTetT::RHSDriver(void)
{
	if (fANP) {
		ComputeAllJe();
		fANP->ComputeJBar(fJe, fJbarE);
	}
	UpdatedLagrangianT::RHSDriver();
}

/* Same pre-pass for the stiffness assembly so the numerical / analytical
 * tangent uses the same F_bar that the residual was based on. */
void BonetTetT::LHSDriver(GlobalT::SystemTypeT sys_type)
{
	if (fANP) {
		ComputeAllJe();
		fANP->ComputeJBar(fJe, fJbarE);
	}
	UpdatedLagrangianT::LHSDriver(sys_type);
}
