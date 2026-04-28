/* ANPHelperT.h — average nodal pressure (Bonet-Burton 1998) F-bar helper.
 *
 * Equivalent to LS-DYNA ELFORM=13 for tet4: averages the dilatation
 * J = det(F) at nodes (volume-weighted) to eliminate volumetric locking.
 * Designed to be reusable between integrator types — the explicit batched
 * driver and the classic-Tahoe implicit element call this same helper.
 *
 * Algorithm (Bonet & Burton 1998):
 *   Per element e:  J_e = det(F_e),  V_e^ref = reference volume (constant)
 *   Per node n:     J_n = (sum_{e ∋ n} V_e^ref * J_e) / (sum_{e ∋ n} V_e^ref)
 *   Per element e:  J_bar_e = (1/nen) * sum_{n in e} J_n
 *   F_bar_e = (J_bar_e / J_e)^(1/nsd) * F_e
 *
 * Pure data interface (no inheritance, no virtual dispatch).  All memory
 * is owned externally except the per-node accumulators allocated by Init.
 */

#ifndef _ANP_HELPER_T_H_
#define _ANP_HELPER_T_H_

namespace Tahoe {

class ANPHelperT
{
public:
	ANPHelperT(void);
	~ANPHelperT(void);

	/** Initialize the helper.  Called once at simulation start.
	 *  \param nelem   total number of elements
	 *  \param nnod    total number of nodes (for sizing nodal accumulators)
	 *  \param nen     nodes per element (4 for tet, 8 for hex, ...)
	 *  \param conn    flat connectivity, [nelem * nen], 0-based node IDs
	 *  \param V_ref_e reference volume per element, [nelem]
	 *
	 * The helper retains pointers to conn and V_ref_e but does not own
	 * them; the caller must keep them valid for the lifetime of the helper.
	 */
	void Init(int nelem, int nnod, int nen,
	          const int* conn, const double* V_ref_e);

	/** Compute J_bar per element from raw J per element.
	 *  \param J_e     input: J = det(F) per element, [nelem]
	 *  \param J_bar_e output: averaged J per element, [nelem]
	 *
	 *  Two-pass internal:
	 *    pass 1 — scatter V_e^ref * J_e into nodal numerator;
	 *             nodal denominator (sum V_e^ref) is precomputed at Init.
	 *    pass 2 — divide; gather J_n back to J_bar_e per element.
	 */
	void ComputeJBar(const double* J_e, double* J_bar_e);

	/** Apply F-bar correction: F_bar = (J_bar/J)^(1/3) * F  (3D).
	 *
	 *  Layout matches the SoA batch convention: F is 9 components, each
	 *  stored as a contiguous row of length `stride` (typically MVSIZ).
	 *  Component order: F11, F12, F13, F21, F22, F23, F31, F32, F33.
	 *
	 *  \param nelem    number of elements to process (≤ stride)
	 *  \param stride   minor-dimension stride of the F arrays
	 *  \param F_e_ptr  pointer to F_e[0][0] — total length 9*stride
	 *  \param J_e      [nelem] det(F) per element
	 *  \param J_bar_e  [nelem] nodal-averaged J per element
	 *  \param F_bar_e_ptr  output, may alias F_e_ptr in-place
	 */
	static void ApplyFBar3D(int nelem, int stride,
	                        const double* F_e_ptr,
	                        const double* J_e,
	                        const double* J_bar_e,
	                        double* F_bar_e_ptr);

	/** Number of elements registered. */
	int NumElements(void) const { return fNelem; }

private:
	int fNelem;
	int fNnod;
	int fNen;
	const int*    fConn;       /* not owned */
	const double* fVrefE;      /* not owned */
	double*       fSumVrefN;   /* owned: per-node sum_{e ∋ n} V_e^ref */
	double*       fJN;         /* owned: per-node J after averaging */
};

} /* namespace Tahoe */
#endif /* _ANP_HELPER_T_H_ */
