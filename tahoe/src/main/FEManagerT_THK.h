/* $Id: FEManagerT_THK.h,v 1.3 2007/08/14 16:22:44 d-farrell2 Exp $ */

#ifndef _FE_MANAGER_THK_H_
#define _FE_MANAGER_THK_H_

/* element configuration header */
#include "ElementsConfig.h"
#if defined(BRIDGING_ELEMENT)

/* base class */
#include "FEManagerT_bridging.h"
#include "dMatrixT.h"

/* File related to calculation using the Wagner-Karpov-Liu Bridging Scale Method
 * If you make use of this code, please cite the following publications (they are also handy references for the method and implementation):
 * 
 *		1)	Wagner GJ, WK Liu. "Coupling of Atomistic and Continuum Simulations using a Bridging Scale Decomposition", Journal of Computational Physics, 190:249-274 (2003)
 *		2)	Wagner GJ, EG Karpov, WK Liu. "Molecular Dynamics Boundary Conditions for Regular Crystal Lattices", CMAME, 193(17-20):1579-1601 (2004)
 *		3)	Park HS, WK Liu. "An Introduction and Tutorial on Multiple Scale Analysis in Solids", CMAME, 193(17-20):1733-1772 (2004)
 *		4)	Park HS, EG Karpov, PA Klein, WK Liu. "Three-Dimensional Bridging Scale Analysis of Dynamic Fracture", Journal of Computational Physics, 207(2):588-609 (2005)
 *		5)	Park HS, EG Karpov, WK Liu, PA Klein. "The Bridging Scale for Two-Dimensional Atomistic/Continuum Coupling", Philosophical Magazine, 85(1):79-113 (2005)
 *		6)	Farrell DE, HS Park, WK Liu. "Implementation Aspects of the Bridging Scale Method and Application to Intersonic Crack Propagation", IJNME 71:583-605 (2007)
 *		7)	Farrell DE, EG Karpov, WK Liu. "Algorithms for Bridging Scale Method Parameters", Computational Mechanics, in print DOI: 10.1007/s00466-007-0156-z (2007)
 */


namespace Tahoe {

/* forward declarations */
class ParticlePairT;

/** FEManagerT to support the time history kernel (THK) formulation */
class FEManagerT_THK: public FEManagerT_bridging
{
public:

	/** constructor */
	FEManagerT_THK(const StringT& input, ofstreamT& output, CommunicatorT& comm,
		const ArrayT<StringT>& argv, TaskT task);

	/** return array containing atom numbers of boundary and ghost atoms **/
	const iArrayT& InterpolationNodes(void);
	
	/** predictor routine for FEM solution interpolated to MD boundary atoms.
		predictor and corrector combined because of constant acceleration assumption.  **/
	void BAPredictAndCorrect(double timestep, dArray2DT& badisp, dArray2DT& bavel, dArray2DT& baacc);
	
	/** calculate THK displacement for ghost atoms for 2/3D disp formulation **/
	const dArray2DT& THKDisp(const StringT& bridging_field, const dArray2DT& badisp);
	
	/** calculate THK displacement for ghost atoms for 2/3D disp formulation, use Beta form **/
	const dArray2DT& BetaTHKDisp(const StringT& bridging_field, const dArray2DT& badisp, const dArray2DT& bavel);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;
	
	/** set up new subordinate parameter list */
	ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	
	/** 2D/3D MD/THK and BSM THK Initialization */
	void InitializeTHK(bool ignore_continuum);
	
	/** accessor for the ghostoffmapping */
	const nMatrixT<int> GetGhostMap(void) { return fghostoffmap;};
	
	/** accessor for the THK type */
	const StringT GetTHKType(void) { return fTHK_type;};

	/*@}*/

	

private:
	
	/** perform neighbor search for THK boundary atoms, 2D */
	void DoNeighSearch2D(void);
	
	/** perform neighbor search for THK boundary atoms, 3D */
	void DoNeighSearch3D(void);
	
	/** find the ghost atom properties map */
	void DoGhostMap(void);
	
	/** compute theta tables for 2D/3D disp/disp or force/disp formulation (doesn't matter, its all the same) */
	void ComputeThetaTables(void);
	
	/** compute beta tables for 2D/3D disp/vel or force/vel formulation (doesn't matter, its all the same) */
	void ComputeBetaTables(void);

private:

	/** \name input parameters */
	/*@{*/
	int fNcrit;
	double fTcut, fLatticeParameter, fSearchParameter;
	StringT fThetaFile, fGhostMapFile;
	ArrayT<StringT> fTHKNodes, fTHKGhostNodes;
	/*@}*/

	int fN_times, fNumstep_crit, fNeighbors;
	dArray2DT fTHKforce, fTHKdisp;
	iArrayT fShift;

	/** interpolation points */
	iArrayT fInterpolationNodes;
	
// DEF added these:
	
	// ghostoffmap matrix
	nMatrixT<int> fghostoffmap;
	iArrayT fthk_bound_lengths;
	iArray2DT fthk_boundary_atoms;
	int fmax_thk_bound_length;
	int fnumsets;
	
	int ftotal_b_atoms, ftotal_g_atoms; // total number of boundary atoms, ghost atoms
	
	// array of THK BC plane normals
	ArrayT<dArrayT> fTHK_normals; 
	
	/** atoms in each node set (ghost atoms): [boundary_n] x [n_boundary_atoms]  (no repeats) */
	ArrayT<iArrayT> fghost_set_atoms;
	
	/** atoms in each node set (real atoms): [boundary_n] x [n_boundary_atoms]  (has repeats) */
	ArrayT<iArrayT> fbound_set_atoms;
	
	/** displacement history: [boundary_n] x [[n_boundary_atoms] x [time x ndof]] (has repeats) */
	ArrayT< ArrayT<dArray2DT> > fHistoryTable;
	
	/** THK values: [boundary_n] x [[neighbor] x [time x ndof*ndof]] (has repeats)*/
	ArrayT< ArrayT<dArray2DT> > fThetaTable_array;
	
	/** boundary neighbors: [boundary_n] x [[n_boundary_atoms] x [n_neighbors]] (has repeats)*/
	ArrayT<iArray2DT> fbound_neighbor_atoms;
	
	/** file containing the fourier coefficients for the sine series used to calculate Theta matrices */
	// size: [n_sets] - 1 entry per set
	ArrayT<StringT> fThetaFile_array;
	
	// parameter which allows the scaling of Tcrit and the Bn's for the sine series (defaults to 1)
	// This allows fourier sine coeffs calculated for k,m = 1 to be used for any k, m it is sqrt(k/m) for the desired values
	double fOmega_sys;
	
	// parameter which allows choice of theta type (displacement-displacement) or beta type (displacement-velocity) THK
	// choices are beta, theta. default is theta
	StringT fTHK_type;
};

} /* namespace Tahoe */

#endif /* BRIDGING_ELEMENT */
#endif /* _FE_MANAGER_THK_H_ */
