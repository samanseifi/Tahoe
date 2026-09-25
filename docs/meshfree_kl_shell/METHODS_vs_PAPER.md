# Method-by-method: our RKShellT vs Wang & Bazilevs (2025) — and a ratification plan

> **Superseded (2026-09-25).** The current, verified account of the formulation, deviations and results is
> `technical_report.pdf` in this directory (source `technical_report.tex`). This note is kept for history;
> several of its statements are out of date (Scordelis–Lo, hemisphere, stabilization type, necking, Fig 18).


Goal: enumerate every methodological component, mark whether it matches the paper, rate its impact on
the quantitative result, and give a concrete test/fix to **ratify** each one (i.e. prove it matches, or
make it match). "Ratified" = verified equal to the paper either by a unit test or by closing the gap.

## A. Comparison table

| # | Component | Paper (Wang & Bazilevs 2025) | Our RKShellT | Match? | Impact | Ratification action |
|---|-----------|------------------------------|--------------|--------|--------|---------------------|
| 1 | KL kinematics (aux tensors A, B1, B2; Eqs 39–55) | analytic auxiliary-tensor form | same (`KLShellKernels.h`) | ✅ | — | **DONE** — FD-validated to machine precision (objectivity, stretch, bending, cylinder κ) |
| 2 | Local parameterization | PCA tangent plane (neighbor cloud) | PCA (`PCAFrame`) | ✅ | low | **DONE** — vase normal/curvature convergence replicated (§4.1) |
| 3 | RK completeness | quadratic (cubic optional) | quadratic (`completeness=2`) | ✅ | med | regression: confirm completeness/reproducing conditions to 1e-12 |
| 4 | **RK kernel** | **cubic B-spline, C²** | cubic B-spline by default (`kernel=1`); Gaussian compatibility mode | 🟡 | **high** | on-node singularity fixed, but dilation 1.0 gives Scordelis −0.194 vs −0.292; recalibrate without breaking pinch |
| 5 | Base mid-surface quadrature | nodal (value at node) × nodal area | same (direct nodal point-sample) | ✅ | med | confirm A_K (nodal area) definition matches (Voronoi vs spacing²) |
| 6 | Through-thickness quadrature | 3-pt Gauss | 3-pt Gauss | ✅ | low | **DONE** |
| 7 | **Membrane stabilization** | natural Taylor σ,ξ (Eq 33/57): BV,ξ : σ,ξ · M_K | σ,ξ stress-gradient (`stab_siggrad`) | 🟡 | **high** | see #8 — structure matches, fidelity differs |
| 8 | **Stress-gradient update** | full co-rotational σ̃,ξ = σ̃,ξ + Δt·C̃ᴾ·D̃,ξ (Eqs 67–72), **current-config** D,ξ | current-config B,ξ + **mid-surface** tangent only | 🟡 | **high** | use per-station tangent; FD-check σ,ξ vs Eq 68 |
| 9 | Extra stabilizers (SCNI penalty, bending penalty) | **none** | present in code, **set to 0** in necking decks | 🟡 | — | confirm OFF in every paper-faithful deck; gate them out of the default path |
| 10 | Co-rotational stress update | Green–Naghdi + Flanagan–Taylor (Alg 2) | same (`FlanaganTaylor.h`, `corotational=1`) | ✅ | med | **DONE** — unit test (rotation/stretch) + objectivity check |
| 11 | Plane stress σ33=0 | Newton on D̃33 + 3D return (Alg 3) | secant on D33 + 3D J2 (`PlaneStressJ2_D33`) | ✅ | med | **DONE** — uniaxial J2 reproduces 1.0× at finite strain |
| 12 | Thickness update | Padé (1+ΔtD/2)/(1−ΔtD/2), Eq 79 | same (`ThicknessStretchPade`) | ✅ | negl | **DONE** — exact paper update, independently tested |
| 13 | Material (plane-stress J2 + sat. hardening) | standard return map | bisection consistency return | 🟡 | low | benchmark a single point vs a reference return (Simo–Hughes) to 1e-8 |
| 14 | Time integration | explicit central difference (velocity Verlet) | central_difference (Tahoe) | ✅ | low | confirm the lumped-mass/half-step convention matches |
| 15 | Loading | velocity control, 30 m/s ramp over 5e-4 s | `D_u` velocity BC, same ramp | ✅ | med | **DONE** (this turn) — velocity BC verified |
| 16 | Lumped mass | translational only (rotational h² dropped, Eq 28) | same | ✅ | low | **DONE** |
| 17 | Mesh | 8284 nodes | 8262 nodes (102×81) | 🟡 | med | regenerate at exactly the paper's node count/layout |
| 18 | Essential BC enforcement | penalty (static) / nodal collocation (dynamic) | Tahoe KBC on (non-interpolatory) coefficients | ❌ | med | for dynamics the paper uses collocation; verify our KBC ≈ collocation on the driven ring |
| 19 | Seed/imperfection | **none** (thickness-update triggers necking) | none (removed this turn) | ✅ | high (location) | **DONE** — exact-paper deck has no seed |

Legend: ✅ matches & verified · 🟡 same intent, detail differs · ❌ genuinely different.

## B. Where the quantitative gap most likely lives
Ranked by expected effect on the necking ultimate/peak (currently ~9% low) and the softening shape:
1. **#8 σ,ξ fidelity** (mid-surface-tangent approximation) — highest leverage, fully in our control.
2. **#18 BC enforcement** + **#17 mesh** — shift localization onset/band.
3. Dynamic ringing from 30 m/s loading — biases the measured reaction (not a formulation error).

Note: the **post-peak ε_p (the "2.0") and softening rate are regularization-dependent** for local
plasticity and are NOT expected to bit-match any independent code; judge ratification on the **rise +
ultimate**, which are physically determined.

## C. Ratification plan (phased, each with a pass criterion)

**Phase 0 — lock the verified matches (1–2 h).**
- Re-run the existing unit tests for #1,#2,#6,#10,#11,#16 and record pass. Add a reproducing-condition
  test for #3 (RK reproduces 1, ξ, ξ² to 1e-12). *Pass:* all green, logged.

**Phase 1 — close #8 (the σ,ξ stress-gradient), highest leverage.**
1. Keep rebuilding B,ξl in the **current configuration** each step (now implemented like the base Bv).
2. Use the **per-through-thickness-station** algorithmic tangent (or the paper's plane-stress C̃ᴾ at the
   mid-surface but applied to D̃,ξ properly), and verify the increment equals Eq 68 by finite difference.
- *Pass:* a patch-test where a known σ,ξ field is reproduced to ~1e-6; necking peak moves toward 605.

**Phase 2 — #4 kernel parity. Reopened after fresh rerun.**
- The **cubic B-spline** is now the default in `MLSSolverT`; Gaussian remains available as `kernel=0`.
- The on-node spline derivative limit is finite, but the current default dilation 1.0 gives Scordelis
  −0.194 on the checked-in mesh. Dilation 1.4 improves it to −0.269, but destabilizes the quadratic
  pinched-cylinder configuration. A single paper-faithful normalization has not yet been ratified.

**Phase 3 — #18 BC + #17 mesh (isolate discretization).**
- Build the exact 8284-node layout; for dynamics, enforce the driven ring by **collocation** and
  compare to the KBC path. *Pass:* peak/location insensitive to which (≤1–2%), or adopt the paper's.

**Phase 4 — separate dynamics from formulation.**
- Run the necking **quasi-statically** (slow ramp + light damping) to remove ringing, AND at the paper's
  30 m/s. *Pass:* the quasi-static ultimate is rate-independent and matches the IGA refs' ultimate;
  the 30 m/s curve carries the paper's ringing. This tells us how much of the 9% is dynamic vs formulation.

**Phase 5 — clean the default path (#9 and #13; #12 complete).**
- Make the natural σ,ξ the sole default stabilizer (SCNI/bending penalties off unless explicitly asked);
  cross-check the J2 return against a reference. *Pass:* paper-faithful
  deck runs with a single stabilizer, regression suite still green.

**Acceptance for "ratified against the paper":**
- Rise + **ultimate within ~3%** of the IGA reference (Alaydin/Ambati Fig 8) after Phases 1–4, AND
- qualitative localization (central neck, deep thinning) reproduced, AND
- every row above either ✅ (unit-verified) or its residual difference quantified and shown immaterial.
</content>
