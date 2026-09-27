# Brinell indentation — rigid ball on elastoplastic block (#47)

Quasi-static implicit indentation of a near-rigid hemispherical indenter
into a J2-plastic block, with light Coulomb friction.  Same quarter-
symmetry topology as the Hertz benchmark (`level.5/hertz/`), but the
base block is elastoplastic and the indenter is ~10× stiffer, so plastic
flow concentrates in the block.

The benchmark exercises three capabilities together for the first time:

| Capability | First landed in | Used here as |
|------------|-----------------|--------------|
| `cubic_spline` isotropic hardening for J2 | original Tahoe | takes 6 (ε_p, σ_y) points → C¹ spline through them |
| Implicit Coulomb friction in `PenaltyContact3DT` | PR #45 (#40) | μ = 0.1 sliding friction at the contact pair |
| Newton + line search (`<nonlinear_solver_LS>`) | original Tahoe | required for robust convergence under changing contact set + plasticity |

## Files

| File | Notes |
|------|-------|
| [`generate_brinell_mesh.py`](generate_brinell_mesh.py) | Writes both mesh variants (smoke, fine). |
| [`brinell_smoke.xml`](brinell_smoke.xml) | **Smoke variant** — 2 304 Hex8, 10 steps to δ = 0.15 mm, converges serial in ~55 s.  This is what CI / autonomous loops should run. |
| [`brinell.xml`](brinell.xml) | Fine variant — 8 800 Hex8, 20 steps to δ = 0.30 mm.  Slower (~1 h serial); deeper into the fully-plastic regime, suitable for Tabor-relation validation. |
| `brinell_smoke.geom`, `brinell.geom` | Generated meshes (text Tahoe geom format). |
| [`compare_to_tabor.py`](compare_to_tabor.py) | Post-processing — extracts δ, P, contact radius `a`, mean pressure `p_m` from the Exodus output; writes a CSV plus two PNGs (P-δ vs Hertz reference, and `p_m/σ_y0` vs `δ/R` with the Tabor 2.8 line). |

## Setup

Both variants share the same geometry shape and BCs:

```
Indenter (block 1): hemisphere R = 5 mm, quarter-symmetry,
                    near-rigid elastic (E = 2 000 000 MPa, ν = 0.3)
Block    (block 2): R_box × R_box × H_base = 2.5 × 2.5 × 3 mm,
                    J2 plastic, E = 200 000 MPa, ν = 0.3,
                    σ_y0 = 250 MPa
                    cubic-spline hardening through:
                        (ε_p, σ_y) = (0.000, 250), (0.005, 320),
                                     (0.020, 420), (0.050, 500),
                                     (0.150, 580), (0.500, 650)  [MPa]
```

Boundary conditions:
- indenter top (NS1): prescribed `u_z = -δ_max · schedule(t)`
- indenter and block symmetry on x=0 and y=0 faces
- block bottom (NS4): fully clamped
- contact pair: `<contact_3D_penalty>` with `μ = 0.1`, regularised slip
  scale `friction_epsilon_velocity = 1e-4`, `penalty_stiffness = 1e6` (see *Contact penalty* below)

Solver: `<nonlinear_solver_LS>` + SPOOLES, `max_iterations=40`,
`search_iterations=5`.

## Running

```bash
cd benchmark_XML/level.5/brinell
python3 generate_brinell_mesh.py        # writes brinell_smoke.geom + brinell.geom
../../../build/bin/tahoe -f brinell_smoke.xml   # ~3.5 min on serial SPOOLES
# or for the fine variant:
../../../build/bin/tahoe -f brinell.xml         # ~1.6 h on one core
```

## Convergence (brinell_smoke.xml)

All 10 steps converge in 6–7 Newton iterations.  Step 10 (deepest load, δ = 0.15 mm, fully plastic):

```
init: 0 LS: 3 ... | 0: Rel error = 1.37e-02
                  | 1: Rel error = 1.55e-03
                  | 2: Rel error = 5.19e-04
                  | 3: Rel error = 2.02e-04
                  | 4: Rel error = 1.98e-05
                  | 5: Rel error = 2.21e-07    ← converged
```

Contact patch grows monotonically as plastic deformation accumulates:
48 active strikers at step 1 → 80 at step 9 (plastic impression spreading).

## Smoke results (from `compare_to_tabor.py brinell_smoke`)

10 frames, δ from 0.015 → 0.150 mm (`penalty_stiffness = 1e6`). `p_m` is the quarter load over the
quarter contact disc, `Y_r` the flow stress at Tabor's representative strain `ε_r = 0.2 a/R`:

| Frame | δ [mm] | P_quarter [N] | a [mm] | p_m [MPa] | p_m / σ_y0 | p_m / Y_r |
|------:|-------:|--------------:|-------:|----------:|-----------:|----------:|
| 0 | 0.015 |   93.3 | 0.308 | 1253.7 | 5.02 | 3.26 |
| 4 | 0.075 |  681.8 | 0.791 | 1386.6 | 5.55 | 3.03 |
| 9 | 0.150 | 1317.9 | 1.116 | 1346.6 | 5.39 | 2.75 |

Over δ ≥ 0.06 mm the coarse mesh gives `p_m / Y_r` = 2.95 on average (2.75–3.16), 5 % above Tabor.

## Contact penalty (2026-09-27, issue #47)

The penalty force is `k · g · A` per striker, so `k` is a pressure per unit penetration and the
penetration is `g ≈ p/k`. At Tabor pressure (~700 MPa) that is 7e-5 mm for `k = 1e7` and 7e-4 mm
for `k = 1e6`.

The decks used `1e7` (50× the block modulus) until 2026-09-27. Every change of the active contact
set then cost a burst of Newton iterations. `brinell_smoke.xml` needed 7, 19, 10, 6, 21, 11, 7, 7, 7
and 8 iterations per step at `1e7`, and 6, 6, 6, 6, 6, 6, 7, 7, 6 and 6 at `1e6`, 35 % less wall
time. The load and the Tabor ratio barely move:

| δ [mm] | P_quarter, k = 1e6 [N] | P_quarter, k = 1e7 [N] | ΔP | p_m/Y_r (1e6 / 1e7) |
|------:|------:|------:|------:|------:|
| 0.015 |   93.3 |   97.2 | −4.1 % | 3.26 / 3.40 |
| 0.060 |  513.8 |  520.5 | −1.3 % | 3.16 / 3.20 |
| 0.090 |  813.4 |  821.5 | −1.0 % | 3.14 / 3.17 |
| 0.150 | 1317.9 | 1325.1 | −0.6 % | 2.75 / 2.77 |

The difference is a fixed penetration offset, so it matters only at the shallowest frames and is
below 1 % in the δ ≳ 0.1 mm range where the Tabor ratio is read. Both decks now use `1e6`, and
`brinell.xml` allows five load-step cuts.

The iteration counts are not caused by the J2 tangent. A finite-difference global stiffness gives
the same Newton history as the analytic `Simo_J2` tangent (#78, which fixed a small term found by
that check).

## Fine-run results (2026-09-27, issue #47)

`brinell.xml` (8 800 Hex8, 20 steps to δ = 0.30 mm, δ/R = 6 %) now runs to the end.

| Run | Wall time (1 core) | Newton iterations per step | Step cuts |
|-----|------:|------|------:|
| `penalty_stiffness = 1e6`, `max_step_cuts = 5` (committed deck) | 97 min | 6–10 | 0 |
| `penalty_stiffness = 1e7`, `max_step_cuts = 5` | 217 min | 6–42 | 1 (at δ = 0.09 mm) |
| `penalty_stiffness = 1e7`, no cuts (2026-09-26) | stalled at δ = 0.09 mm | 10–42 | – |

Tabor check for the committed deck (`compare_to_tabor.py brinell`):

| δ [mm] | a/R | P_quarter [N] | p_m [MPa] | ε_r = 0.2 a/R | Y_r [MPa] | p_m / Y_r |
|------:|----:|------:|------:|------:|------:|------:|
| 0.030 | 0.100 |  218.9 | 1107.9 | 0.020 | 420.2 | 2.64 |
| 0.060 | 0.143 |  507.9 | 1267.5 | 0.029 | 448.2 | 2.83 |
| 0.090 | 0.178 |  790.6 | 1276.6 | 0.036 | 467.6 | 2.73 |
| 0.120 | 0.204 | 1041.2 | 1273.7 | 0.041 | 480.5 | 2.65 |
| 0.150 | 0.220 | 1277.4 | 1338.8 | 0.044 | 487.9 | 2.74 |
| 0.180 | 0.244 | 1503.4 | 1289.5 | 0.049 | 497.5 | 2.59 |
| 0.210 | 0.256 | 1706.5 | 1324.5 | 0.051 | 502.4 | 2.64 |
| 0.240 | 0.275 | 1911.9 | 1283.1 | 0.055 | 509.5 | 2.52 |
| 0.270 | 0.288 | 2099.9 | 1293.7 | 0.058 | 513.6 | 2.52 |
| 0.300 | 0.302 | 2269.1 | 1262.9 | 0.061 | 518.6 | 2.44 |

Over the fully plastic frames (δ ≥ 0.06 mm) `p_m / Y_r` averages **2.63 (2.44–2.83), 6 % below
Tabor's 2.8**; the `1e7` run gives 2.67 (2.54–2.87). The coarse smoke mesh gives 2.95, so the two
meshes bracket 2.8. The frame-to-frame scatter comes from the contact radius `a`, which moves in
steps of one node ring: at δ = 0.30 mm the two penalties agree on the load to 0.4 % but the `1e6` run
picks up one more ring, lowering its `p_m` by 5 %. The slow drift below 2.8 at depth is consistent
with pile-up widening the measured contact ring.

**Tabor is validated to within 10 % on average**, with the deepest frame at −13 %.

The first report of this run (2026-09-26, `p_m / σ_y0 ≈ 1.3`) was wrong for two reasons, both fixed
in `compare_to_tabor.py`: the quarter-model load was divided by the full contact disc (a factor of
4), and a hardening material has to be compared at Tabor's representative flow stress, not at the
initial yield stress. With `σ_y0` the ratio is about 5.

## Expected physics

1. **Elastic phase** (δ ≲ 0.02 mm) — `P-δ` matches Hertz `P = (4/3) E* √R δ^{3/2}`.
2. **Yield onset** — maximum von Mises beneath the indenter reaches σ_y0 = 250 MPa
   at indentation depth `δ_y ≈ 0.012 mm` (Hertz analytical `p₀,y = 1.6 σ_y0`).
3. **Plastic phase** — contact patch grows faster than Hertz; mean pressure
   `p_m = P / πa²` flattens.
4. **Tabor relation** — in the fully plastic regime `p_m ≈ 2.8 Y_r`, with `Y_r` the flow stress
   at the representative strain `ε_r ≈ 0.2 a/R` (Tabor 1951). Checked above.

## See also
- Issue #47 — Brinell benchmark tracking
- `level.5/hertz/` — same topology, elastic-only Hertz validation
- `level.5/implicit_friction/sliding_cubes.xml` — bare Coulomb friction test (#40)
- `level.5/tet_classic/tet4_hyperelastic_anp.xml` — implicit Newton + ANP-Tet4 (#29)

## Limitations

- The "rigid" indenter is actually a stiff elastic body (E = 2 × 10⁶ MPa).
  Its compliance contributes ~0.5 % of total displacement.  A true rigid-body
  element would be cleaner — tracked as #30.
- The fine variant (`brinell.xml`) takes about 1.6 h serial with line-search
  Newton.  Each step is about 5 min, dominated by SPOOLES factorisation of
  the 29 022-DOF tangent.  For repeated runs, MUMPS or
  SuperLU would speed this up; not changed here so the benchmark is
  reproducible on a default build with no optional flags.
- Unloading / residual depth is not in this benchmark; would require an
  additional load schedule and a small code-side change to record
  per-step results without rewriting from scratch.

## Plastic strain + von Mises field on the block top (last frame)

Material output from `<Simo_J2>` (enabled by `material_output="1"` on the
plastic-block element group) gives per-node `alpha` (equivalent plastic
strain ε_p), `VM_Kirch` (von Mises Kirchhoff), `press`, and `norm_beta`
(back-stress norm — zero here since isotropic hardening only).
`compare_to_tabor.py` slices the block top surface (z = 0 in the
reference config) at the last frame and writes `*_field.png`.

| Quantity | Range (smoke run, δ = 0.15 mm) |
|----------|--------------------------------|
| `α` (ε_p) | 0.000 → 0.106 — significant plastic flow at the impression centre |
| `σ_VM`    | 50 → 566 MPa  — saturates at ≈ σ_y(α = 0.1) per the cubic-spline curve |

The 566 MPa peak matching the spline value σ_y(0.10) ≈ 560 MPa is direct
evidence that the radial-return is converging to the yield surface — a
useful sanity check on the J2 update under contact loading.

> **Output-channel layout (after splitting the two materials into separate
> element groups so the J2 outputs can be emitted without clashing with
> the elastic indenter).**
> ```
> *.io0.exo   — indenter (block 1) — D_X/Y/Z, Cauchy stress
> *.io1.exo   — plastic block (block 2) — same + alpha, norm_beta, VM_Kirch, press
> *.io2.exo   — contact group — D_X/Y/Z and F_X/Y/Z on strikers
> ```
