# Closing the Wang & Bazilevs (2025) gaps — implementation notes

This documents the formulation gaps closed in `RKShellT` (epic #59) to make the element
paper-faithful, and how each was tested. Source: `development/src/elements/meshfree_kl_shell/`.

The five gaps were identified by comparing the paper section-by-section against the code
(see the review in the session log). They are addressed below in the order of leverage for
the necking benchmark (§4.3).

---

## Gap 1 + 3 — natural stabilization via the accumulated stress gradient σ,ξ (paper §3.10, Eqs 67–72)

**Was:** the membrane stabilizer used `C_ps^alg · (B,ξl · u_total)` — a consistent-tangent times
the *total* strain gradient — or a penalty. Both keep full perpendicular/volumetric stiffness in the
yielding neck, so localization is over-smoothed (diffuse neck, ε_p ~0.5 vs the paper's 2.0).

**Now (`stab_siggrad=1`, default):** the stabilization force uses the **accumulated Cauchy stress
gradient** σ,ξl as a history variable, exactly as the paper does:

- per node, `fSigGrad[i]` holds `[σ,ξ1 ; σ,ξ2]` (two in-plane symmetric tensors, local frame);
- each step it is advanced by `σ,ξl += C̃^P_mid · (B,ξl · du)` where `C̃^P_mid` is the **mid-surface
  (ξ3=0) plane-stress ALGORITHMIC tangent** (`PlaneStressJ2Calg`), captured at the g==1 Gauss station;
- it is co-rotated with the Flanagan–Taylor frame each step (Eq 72);
- the stab force is `Σ_l (B,ξl)^T σ,ξl · V_K M_Kξl`.

**Why it sharpens the neck:** once the mid-surface yields, `C̃^P_mid` collapses along the plastic-flow
normal, so the *increment* to σ,ξl along that direction → 0 and σ,ξl **saturates**. The stabilizer
therefore stops adding stiffness in the localizing direction (no elastic clamp) while still controlling
hourglass in the orthogonal/volumetric directions. This is the single high-leverage item the prior
`DEVIATIONS_FROM_PAPER.md` / `RESULTS_SUMMARY.md` root-caused as blocking ε_p → 2.0.

**Approximation noted:** the stored `B,ξl` operators are reference-config (built once); the per-step
gradient strain increment uses them with the current `du`, and objectivity of the *accumulated* σ,ξl is
restored by the co-rotation. A fully current-config `B,ξl` rebuild each step is a possible refinement.

## Gap 2 — objective stress update is now the default (paper §3.8, Algorithm 2)

`corotational` now **defaults to 1**. The Green–Naghdi / Flanagan–Taylor co-rotational frame
(`FlanaganTaylor.h`) was already implemented and unit-verified but was opt-in; the legacy path
re-derived an arbitrary in-plane frame each step (not objective for in-plane shear under finite
rotation). Linear / small-strain runs are unaffected (the frame only acts in `InternalForceFS`).

## Gap 4 — rotational-continuity penalty coupling (paper §3.12, Eqs 80–82)

New, **off by default** (`penalty_coupling=0`). Preserves the kink angle between adjacent shell patches
at a C0 interface by penalizing the relative **normal-rotation jump** of paired interface nodes. The
linearized normal change `ṅ(u) = Σ_I (B1 Ψ,ξ1_I + B2 Ψ,ξ2_I) u_I` uses each side's own one-sided PCA
chart (the auxiliary tensors already in `KLShellKernels.h`). Interface pairs come from two ordered
node sets: `couple_node_ID_a`, `couple_node_ID_b`. Force `(h³/12) C E/h_pl` is scattered to both sides'
stencils (`AddCouplingForce`), distributed per node in `RHSDriver`.

Unit self-test (`KLSHELL_FEATURETEST=1`): `ṅ(rigid translation) = 0` to 1e-16 — a translated patch does
not rotate its normal, confirming the operator. **Full-problem validation needs the multi-patch
tube-crush case (§4.5)** — not exercised by the smooth necking cylinder.

## Gap 5 — pinball / volumetric-potential self-contact (paper §3.13, Eqs 83–86)

New, **off by default** (`contact_stiffness=0`). A node repels every NON-neighbor node within the
contact band `[r_in, r_out]` using current positions: `f = Σ_K ψ(r_iK) (x_i−x_K)/‖·‖ V_K V_i`, with the
pinball force-density `ψ` (Eq 85) and the `c1/c2` C1-smoothing of Eq 86 (p=2). Parameters:
`contact_stiffness` (kc), `contact_r_in`, `contact_r_out` (auto-set to 1.44 dx / 2.74 dx if omitted).

Validated two ways:
1. unit self-test: `ψ` positive, monotone-decreasing, `ψ(r_out)=0`.
2. end-to-end `twosheets_contact.xml` (gen `gen_twosheets_geom.py`): a top sheet falls under load onto a
   fixed bottom sheet. **Contact ON → arrested at gap ~1.8 (max|u|≈1.2, no interpenetration); OFF → falls
   straight through (max|u|≈40 and growing).**

---

## How to reproduce

- Feature self-tests: `KLSHELL_FEATURETEST=1 tahoe -f <any deck>` (runs at setup).
- Self-contact demo: `tahoe -f twosheets_contact.xml` vs `twosheets_nocontact.xml`.
- Necking, coarse A/B (sharp vs diffuse): `_neckAB_new.xml` (`stab_siggrad=1`) vs `_neckAB_old.xml` (`=0`).
- Necking, paper resolution (overnight): `_neck_paper_siggrad.xml` (8262 nodes, real density).

## Parameter summary (new / changed defaults)

| param | default | meaning |
|---|---|---|
| `stab_siggrad` | 1 | accumulated σ,ξ stabilization (paper §3.10) |
| `corotational` | 1 | Flanagan–Taylor objective stress update (was 0) |
| `contact_stiffness` | 0 | self-contact kc (>0 enables, §3.13) |
| `contact_r_in/out` | 0 | contact band radii (auto from spacing if 0) |
| `penalty_coupling` | 0 | kink-angle penalty C (>0 enables, §3.12) |
| `couple_node_ID_a/b` | none | ordered interface node sets for coupling |
