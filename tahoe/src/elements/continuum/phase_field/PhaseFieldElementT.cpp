/* Phase-field fracture element (AT2 model) */
#include "PhaseFieldElementT.h"

#include <iostream>
#include <iomanip>
#include <cmath>

#include "ifstreamT.h"
#include "ElementCardT.h"
#include "ShapeFunctionT.h"
#include "eIntegratorT.h"
#include "iAutoArrayT.h"
#include "ParameterContainerT.h"

/* materials */
#include "PhaseFieldMaterialT.h"
#include "PhaseFieldMatSupportT.h"
#include "PhaseFieldMatListT.h"

using namespace Tahoe;

/* initialize static data */
const int PhaseFieldElementT::NumNodalOutputCodes = 3;
static const char* NodalOutputNames[] = {
	"coordinates",
	"phase_field",
	"material_output"};

/* phase-field has 1 DOF per node */
const int kPhaseFieldNDOF = 1;

/* constructor */
PhaseFieldElementT::PhaseFieldElementT(const ElementSupportT& support):
	ContinuumElementT(support),
	fLocDisplacement(NULL),
	fPhaseFieldMatSupport(NULL),
	fMechanicalCoupling(false),
	fTotalNumIP(0),
	fMu(0.0),
	fLambda(0.0)
{
	SetName("phase_field");
}

/* destructor */
PhaseFieldElementT::~PhaseFieldElementT(void)
{
	delete fPhaseFieldMatSupport;
	delete fLocDisplacement;
}

/* compute nodal force */
void PhaseFieldElementT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
#pragma unused(field)
#pragma unused(node)
#pragma unused(force)
}

/* returns the stored energy (crack surface energy) */
double PhaseFieldElementT::InternalEnergy(void)
{
	double energy = 0.0;

	Top();
	while (NextElement())
	{
		SetGlobalShape();
		SetLocalU(fLocDisp);

		const double* Det    = fShapes->IPDets();
		const double* Weight = fShapes->IPWeights();

		double Gc  = fCurrMaterial->FractureToughness();
		double ell = fCurrMaterial->LengthScale();

		fShapes->TopIP();
		while (fShapes->NextIP())
		{
			/* interpolate d at this IP */
			fShapes->InterpolateU(fLocDisp, fDOFvec);
			double d = fDOFvec[0];

			/* gradient of d */
			dMatrixT grad(1, NumSD());
			IP_ComputeGradient(fLocDisp, grad);

			/* |grad d|^2 */
			double grad_d_sq = 0.0;
			for (int i = 0; i < NumSD(); i++)
				grad_d_sq += grad(0,i)*grad(0,i);

			/* AT2 crack surface energy: Gc/(2*ell)*d^2 + Gc*ell/2*|grad d|^2 */
			energy += (Gc/(2.0*ell)*d*d + Gc*ell/2.0*grad_d_sq)*(*Det++)*(*Weight++);
		}
	}
	return energy;
}

void PhaseFieldElementT::SendOutput(int kincode)
{
	iArrayT flags(fNodalOutputCodes.Length());
	flags = IOBaseT::kAtNever;
	switch (kincode)
	{
		case iNodalDisp:
		    flags[iNodalDisp] = NumDOF();
			break;
		default:
			cout << "\n PhaseFieldElementT::SendOutput: invalid output code: ";
			cout << kincode << endl;
	}

	iArrayT n_counts;
	SetNodalOutputCodes(IOBaseT::kAtInc, flags, n_counts);
	ElementSupport().ResetAverage(n_counts.Sum());

	iArrayT e_counts(fElementOutputCodes.Length());
	e_counts = 0;

	dArray2DT e_values, n_values;
	ComputeOutput(n_counts, n_values, e_counts, e_values);
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* initialize local arrays */
void PhaseFieldElementT::SetLocalArrays(void)
{
	/* inherited */
	ContinuumElementT::SetLocalArrays();

	int nen = NumElementNodes();

	/* look for displacement field for mechanical coupling */
	const FieldT* fDisplacementVectorField = ElementSupport().Field("displacement");
	if (fDisplacementVectorField) {
		fLocDisplacement = new LocalArrayT(LocalArrayT::kDisp, nen, fDisplacementVectorField->NumDOF());
		fDisplacementVectorField->RegisterLocal(*fLocDisplacement);
	}
}

/* set the correct shape functions */
void PhaseFieldElementT::SetShape(void)
{
	fShapes = new ShapeFunctionT(GeometryCode(), NumIP(), fLocInitCoords);
	if (!fShapes) throw ExceptionT::kOutOfMemory;
	fShapes->Initialize();
}

/* compute deformation gradient from displacement field */
void PhaseFieldElementT::SetGlobalShape(void)
{
	/* inherited */
	ContinuumElementT::SetGlobalShape();

	/* mechanical coupling: compute deformation gradient at IPs */
	if (fMechanicalCoupling && fLocDisplacement)
	{
		SetLocalU(*fLocDisplacement);

		for (int i = 0; i < NumIP(); i++)
		{
			dMatrixT& mat = fF_List[i];
			fShapes->GradU(*fLocDisplacement, mat, i);
			mat.PlusIdentity();
		}
	}
}

/* construct output labels array */
void PhaseFieldElementT::SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
	iArrayT& counts) const
{
	counts.Dimension(flags.Length());
	counts = 0;

	if (flags[iNodalCoord] == mode)
		counts[iNodalCoord] = NumSD();
	if (flags[iNodalDisp] == mode)
		counts[iNodalDisp] = NumDOF();
	if (flags[iMaterialData] == mode)
		counts[iMaterialData] = (*fMaterialList)[0]->NumOutputVariables();
}

void PhaseFieldElementT::SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
	iArrayT& counts) const
{
#pragma unused(mode)
#pragma unused(flags)
	if (counts.Sum() != 0)
		ExceptionT::BadInputValue("PhaseFieldElementT::SetElementOutputCodes", "not implemented");
}

/* construct the effective mass matrix */
void PhaseFieldElementT::LHSDriver(GlobalT::SystemTypeT sys_type)
{
	/* inherited */
	ContinuumElementT::LHSDriver(sys_type);

	double constK = 0.0;
	int formK = fIntegrator->FormK(constK);

	Top();
	while (NextElement())
	{
		fLHS = 0.0;
		SetGlobalShape();

		if (formK) FormStiffness(constK);

		AssembleLHS();
	}
}

void PhaseFieldElementT::RHSDriver(void)
{
	/* inherited */
	ContinuumElementT::RHSDriver();

	double constKd = 0.0;
	int formKd = fIntegrator->FormKd(constKd);

	Top();
	while (NextElement())
	{
		fRHS = 0.0;
		SetGlobalShape();

		if (formKd)
		{
			SetLocalU(fLocDisp);
			FormKd(-constKd);
		}

		AssembleRHS();
	}
}

/* set the B matrix at the specified integration point */
static void PhaseField_B(const ShapeFunctionT* shapes, int ip, dMatrixT& B_matrix)
{
	const dArray2DT& DNa = shapes->Derivatives_U(ip);
	int nnd = DNa.MinorDim();
	double* pB = B_matrix.Pointer();

	if (DNa.MajorDim() == 2)
	{
		const double* pNax = DNa(0);
		const double* pNay = DNa(1);
		for (int i = 0; i < nnd; i++)
		{
			*pB++ = *pNax++;
			*pB++ = *pNay++;
		}
	}
	else
	{
		const double* pNax = DNa(0);
		const double* pNay = DNa(1);
		const double* pNaz = DNa(2);
		for (int i = 0; i < nnd; i++)
		{
			*pB++ = *pNax++;
			*pB++ = *pNay++;
			*pB++ = *pNaz++;
		}
	}
}

/* form the element stiffness matrix
 *
 * K_dd = integral[ Gc*ell*B^T*B + (Gc/ell + 2*H)*N^T*N ] dV
 */
void PhaseFieldElementT::FormStiffness(double constK)
{
	dMatrixT::SymmetryFlagT format =
		(fLHS.Format() == ElementMatrixT::kNonSymmetric) ?
		dMatrixT::kWhole :
		dMatrixT::kUpperOnly;

	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	double Gc  = fCurrMaterial->FractureToughness();
	double ell = fCurrMaterial->LengthScale();

	/* global element index for history variable lookup */
	int elem = CurrElementNumber();
	int nip  = NumIP();

	fShapes->TopIP();
	while (fShapes->NextIP())
	{
		int ip = fShapes->CurrIP();
		double scale = constK*(*Det++)*(*Weight++);

		/* gradient B matrix */
		PhaseField_B(fShapes, ip, fB);

		/* diffusion-like term: scale*Gc*ell*B^T*B */
		fD.Identity(scale * Gc * ell);
		fLHS.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);

		/* reaction term: (Gc/ell + 2*H)*N^T*N */
		int H_index = elem * nip + ip;
		double H = fH_current[H_index];
		double reaction_coeff = scale * (Gc/ell + 2.0*H);

		/* use FormMass-like assembly for N^T*N with scalar coefficient */
		const double* Na = fShapes->IPShapeU(ip);
		int nen = NumElementNodes();
		for (int a = 0; a < nen; a++)
		{
			for (int b = a; b < nen; b++)
			{
				double val = reaction_coeff * Na[a] * Na[b];
				fLHS(a, b) += val;
				if (a != b && format == dMatrixT::kWhole)
					fLHS(b, a) += val;
			}
		}
	}
}

/* calculate the internal force contribution
 *
 * f_d = integral[ Gc*ell*B^T*grad(d) + (Gc/ell + 2*H)*N^T*d - 2*H*N^T ] dV
 */
void PhaseFieldElementT::FormKd(double constK)
{
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	int nsd = NumSD();
	int nen = NumElementNodes();

	double Gc  = fCurrMaterial->FractureToughness();
	double ell = fCurrMaterial->LengthScale();

	int elem = CurrElementNumber();
	int nip  = NumIP();

	dMatrixT grad;
	dArrayT grad_d(nsd);

	fShapes->TopIP();
	while (fShapes->NextIP())
	{
		int ip = fShapes->CurrIP();
		double jw = (*Det++) * (*Weight++);

		/* compute gradient of d at this IP */
		grad.Set(1, nsd, fGradient_list[ip].Pointer());
		IP_ComputeGradient(fLocDisp, grad);
		for (int i = 0; i < nsd; i++)
			grad_d[i] = grad(0, i);

		/* interpolate d at this IP */
		fShapes->InterpolateU(fLocDisp, fDOFvec);
		double d = fDOFvec[0];

		/* get H at this integration point */
		int H_index = elem * nip + ip;
		double H = fH_current[H_index];

		/* update history variable: H = max(H_last, psi(F)) */
		if (fMechanicalCoupling)
		{
			double psi = ComputeStrainEnergyDensity(ip);
			if (psi > H) {
				H = psi;
				fH_current[H_index] = H;
			}
		}

		/* B matrix */
		PhaseField_B(fShapes, ip, fB);

		/* diffusion term: Gc*ell*B^T*grad(d) */
		fB.MultTx(grad_d, fNEEvec);
		fRHS.AddScaled(constK * jw * Gc * ell, fNEEvec);

		/* reaction + source terms via shape functions */
		const double* Na = fShapes->IPShapeU(ip);
		double reaction = (Gc/ell + 2.0*H) * d - 2.0*H;
		for (int a = 0; a < nen; a++)
			fRHS[a] += constK * jw * reaction * Na[a];
	}
}

/* compute strain energy density from deformation gradient at given IP.
 * Uses compressible Neo-Hookean model:
 *   psi = mu/2*(I1 - nsd) - mu*ln(J) + lambda/2*(ln(J))^2
 */
double PhaseFieldElementT::ComputeStrainEnergyDensity(int ip) const
{
	if (!fMechanicalCoupling) return 0.0;

	int nsd = NumSD();
	const dMatrixT& F = fF_List[ip];

	/* right Cauchy-Green tensor C = F^T*F */
	dMatrixT C(nsd);
	C.MultATB(F, F);

	/* I1 = tr(C) */
	double I1 = 0.0;
	for (int i = 0; i < nsd; i++)
		I1 += C(i,i);

	/* J = det(F) */
	double J = F.Det();
	double lnJ = log(J);

	/* compressible Neo-Hookean */
	return fMu/2.0*(I1 - nsd) - fMu*lnJ + fLambda/2.0*lnJ*lnJ;
}

/* construct a new material support and return a pointer */
MaterialSupportT* PhaseFieldElementT::NewMaterialSupport(MaterialSupportT* p) const
{
	if (!p) p = new PhaseFieldMatSupportT(NumDOF(), NumIP());

	/* inherited initializations */
	ContinuumElementT::NewMaterialSupport(p);

	PhaseFieldMatSupportT* ps = TB_DYNAMIC_CAST(PhaseFieldMatSupportT*, p);
	if (ps) {
		ps->SetContinuumElement(this);
		ps->SetGradient(&fGradient_list);
	}

	return p;
}

/* return a pointer to a new material list */
MaterialListT* PhaseFieldElementT::NewMaterialList(const StringT& name, int size)
{
	if (name != "phase_field_material")
		return NULL;

	if (size > 0)
	{
		if (!fPhaseFieldMatSupport) {
			fPhaseFieldMatSupport = TB_DYNAMIC_CAST(PhaseFieldMatSupportT*, NewMaterialSupport());
			if (!fPhaseFieldMatSupport)
				ExceptionT::GeneralFail("PhaseFieldElementT::NewMaterialList");
		}
		return new PhaseFieldMatListT(size, *fPhaseFieldMatSupport);
	}
	else
		return new PhaseFieldMatListT;
}

/* current element operations */
bool PhaseFieldElementT::NextElement(void)
{
	bool result = ContinuumElementT::NextElement();

	if (result)
	{
		ContinuumMaterialT* pcont_mat = (*fMaterialList)[CurrentElement().MaterialNumber()];
		fCurrMaterial = (PhaseFieldMaterialT*) pcont_mat;
	}

	return result;
}

/* driver for calculating output values */
void PhaseFieldElementT::ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
	const iArrayT& e_codes, dArray2DT& e_values)
{
	int n_out = n_codes.Sum();
	int e_out = e_codes.Sum();

	if (n_out == 0 && e_out == 0) return;

#pragma unused(e_values)
	if (e_out > 0)
		ExceptionT::GeneralFail("PhaseFieldElementT::ComputeOutput", "element output not supported");

	int nen = NumElementNodes();
	int nsd = NumSD();

	ElementSupport().ResetAverage(n_out);

	dArray2DT nodal_space(nen, n_out);
	dArray2DT nodal_all(nen, n_out);
	dArray2DT coords, disp;
	dArray2DT matdat;

	dArrayT ipmat(n_codes[iMaterialData]);

	double* pall = nodal_space.Pointer();
	coords.Set(nen, n_codes[iNodalCoord], pall);
	pall += coords.Length();
	disp.Set(nen, n_codes[iNodalDisp], pall);
	pall += disp.Length();
	matdat.Set(nen, n_codes[iMaterialData], pall);

	Top();
	while (NextElement())
	{
		nodal_space = 0.0;
		SetGlobalShape();
		SetLocalU(fLocDisp);

		if (n_codes[iNodalCoord]) fLocInitCoords.ReturnTranspose(coords);
		if (n_codes[iNodalDisp])  fLocDisp.ReturnTranspose(disp);

		fShapes->TopIP();
		while (fShapes->NextIP())
		{
			if (n_codes[iMaterialData])
			{
				fCurrMaterial->ComputeOutput(ipmat);
				fShapes->Extrapolate(ipmat, matdat);
			}
		}

		int colcount = 0;
		nodal_all.BlockColumnCopyAt(disp,   colcount); colcount += disp.MinorDim();
		nodal_all.BlockColumnCopyAt(coords, colcount); colcount += coords.MinorDim();
		nodal_all.BlockColumnCopyAt(matdat, colcount);

		ElementSupport().AssembleAverage(CurrentElement().NodesX(), nodal_all);
	}

	ElementSupport().OutputUsedAverage(n_values);
}

/* information about subordinate parameter lists */
void PhaseFieldElementT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ContinuumElementT::DefineSubs(sub_list);

	/* nodal output codes (optional) */
	sub_list.AddSub("phase_field_element_nodal_output", ParameterListT::ZeroOrOnce);

	/* element block/material specification */
	sub_list.AddSub("phase_field_element_block", ParameterListT::OnePlus);
}

/* return the description of the given inline subordinate parameter list */
ParameterInterfaceT* PhaseFieldElementT::NewSub(const StringT& name) const
{
	if (name == "phase_field_element_nodal_output")
	{
		ParameterContainerT* node_output = new ParameterContainerT(name);

		for (int i = 0; i < NumNodalOutputCodes; i++) {
			ParameterT output(ParameterT::Integer, NodalOutputNames[i]);
			output.SetDefault(1);
			node_output->AddParameter(output, ParameterListT::ZeroOrOnce);
		}

		return node_output;
	}
	else if (name == "phase_field_element_block")
	{
		ParameterContainerT* block = new ParameterContainerT(name);

		/* list of element block ID's */
		block->AddSub("block_ID_list", ParameterListT::Once);

		/* choice of materials lists (inline) */
		block->AddSub("phase_field_material", ParameterListT::Once);

		/* set this as source of subs */
		block->SetSubSource(this);

		return block;
	}
	else /* inherited */
		return ContinuumElementT::NewSub(name);
}

/* accept parameter list */
void PhaseFieldElementT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ContinuumElementT::TakeParameterList(list);

	/* dimensions */
	int nsd = NumSD();
	int nen = NumElementNodes();
	int nip = NumIP();
	fD.Dimension(nsd);
	fB.Dimension(nsd, nen);
	fGradient_list.Dimension(nip);
	for (int i = 0; i < fGradient_list.Length(); i++)
		fGradient_list[i].Dimension(nsd);

	/* nodal output codes */
	fNodalOutputCodes.Dimension(NumNodalOutputCodes);
	fNodalOutputCodes = 0;
	const ParameterListT* node_output = list.List("phase_field_element_nodal_output");
	if (node_output)
		for (int i = 0; i < NumNodalOutputCodes; i++)
		{
			const ParameterT* nodal_value = node_output->Parameter(NodalOutputNames[i]);
			if (nodal_value) {
				int do_write = *nodal_value;
				if (do_write)
					fNodalOutputCodes[i] = 1;
			}
		}

	/* look for displacement field (mechanical coupling) */
	const FieldT* fDisplacementVectorField = ElementSupport().Field("displacement");
	fMechanicalCoupling = false;
	if (fDisplacementVectorField) {
		fMechanicalCoupling = true;
		std::cout << "PhaseFieldElementT: mechanical coupling detected." << std::endl;
	} else {
		std::cout << "PhaseFieldElementT: standalone mode (no mechanical coupling)." << std::endl;
	}

	/* deformation gradient storage */
	fF_all.Dimension(nip*nsd*nsd);
	fF_List.Dimension(nip);
	fF_all = 0.0;
	for (int i = 0; i < nip; i++)
		fF_List[i].Alias(nsd, nsd, fF_all.Pointer(i*nsd*nsd));

	/* allocate history variable H for all elements and IPs.
	 * Total IPs = number_of_elements * nip.
	 * We count elements across all element blocks. */
	fTotalNumIP = NumElements() * nip;
	fH_current.Dimension(fTotalNumIP);
	fH_last.Dimension(fTotalNumIP);
	fH_current = 0.0;
	fH_last = 0.0;

	/* read mechanical parameters if coupling is active */
	const ParameterT* mu_param = list.Parameter("shear_modulus");
	const ParameterT* lambda_param = list.Parameter("lame_lambda");
	if (mu_param) fMu = *mu_param;
	if (lambda_param) fLambda = *lambda_param;
}

/* extract the list of material parameters */
void PhaseFieldElementT::CollectMaterialInfo(const ParameterListT& all_params,
	ParameterListT& mat_params) const
{
	const char caller[] = "PhaseFieldElementT::CollectMaterialInfo";

	mat_params.Clear();
	mat_params.SetName("phase_field_material");

	int num_blocks = all_params.NumLists("phase_field_element_block");
	for (int i = 0; i < num_blocks; i++) {
		const ParameterListT& block = all_params.GetList("phase_field_element_block", i);
		const ParameterListT& mat_list = block.GetList(mat_params.Name());
		const ArrayT<ParameterListT>& mat = mat_list.Lists();
		mat_params.AddList(mat[0]);
	}
}

/***********************************************************************
 * Private
 ***********************************************************************/

/* construct output labels array */
void PhaseFieldElementT::GenerateOutputLabels(const iArrayT& n_codes,
	ArrayT<StringT>& n_labels, const iArrayT& e_codes,
	ArrayT<StringT>& e_labels) const
{
#pragma unused(e_labels)

	n_labels.Dimension(n_codes.Sum());
	int count = 0;

	if (n_codes[iNodalDisp])
	{
		const ArrayT<StringT>& labels = Field().Labels();
		for (int i = 0; i < labels.Length(); i++)
			n_labels[count++] = labels[i];
	}

	if (n_codes[iNodalCoord])
	{
		const char* xlabels[] = {"x1", "x2", "x3"};
		for (int i = 0; i < NumSD(); i++)
			n_labels[count++] = xlabels[i];
	}

	if (n_codes[iMaterialData])
	{
		ArrayT<StringT> matlabels;
		(*fMaterialList)[0]->OutputLabels(matlabels);

		for (int i = 0; i < n_codes[iMaterialData]; i++)
			n_labels[count++] = matlabels[i];
	}

	if (e_codes.Sum() != 0)
		ExceptionT::GeneralFail("PhaseFieldElementT::GenerateOutputLabels",
			"not expecting any element output codes");
}
