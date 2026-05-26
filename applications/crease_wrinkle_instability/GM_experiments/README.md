# Gurtin-Murdoch surface elasticity — experiments (2D)

Strain-dependent surface stress
`σ_s(ε) = γ₀ + E_s · ε_eng,  ε_eng = (L − L₀)/L₀`
applied to selected crease/wrinkle/buckling setups from
[../2D_young_laplace_experiments/](../2D_young_laplace_experiments/) and [../prestretch_workflow/](../prestretch_workflow/).
The underlying element implementation (`SimoQ1P0_Surface`) was added
under issue [#54](https://github.com/samanseifi/Tahoe/issues/54);
research questions for the paper live under issue
[#55](https://github.com/samanseifi/Tahoe/issues/55).

YL (Young-Laplace, `E_s = 0`) is the default and is preserved
bit-identically in the existing 2D_young_laplace_experiments XMLs — these GM variants
just add `E_s` on the `<surface_tension>` element.

## XMLs

| File                                                              | Seeded from                                         | E_s |
| ----------------------------------------------------------------- | --------------------------------------------------- | --: |
| [staggered_explicit_2D_GM.xml](staggered_explicit_2D_GM.xml)       | [../2D_young_laplace_experiments/staggered_explicit_2D.xml](../2D_young_laplace_experiments/staggered_explicit_2D.xml) | 10  |
| [stage2_explicit_eps20_GM.xml](stage2_explicit_eps20_GM.xml)       | [../prestretch_workflow/stage2_explicit_eps20.xml](../prestretch_workflow/stage2_explicit_eps20.xml) | 10  |

### Why these two first

* `staggered_explicit_2D_GM.xml` is the no-prestretch wrinkle case
  (γ̄ = 5).  Top-fiber engineering strain only develops *during* the
  voltage ramp, so the GM contribution `E_s · ε_eng` grows in step
  with the instability — this is the cleanest GM-vs-YL comparison.

* `stage2_explicit_eps20_GM.xml` carries a pre-stretch of ε_pre = 20%
  *into* the surface energy.  At t = 0 of Stage 2 the top fiber is
  already at ε_eng ≈ 0.20 → `E_s · ε_eng = 2.0` ≈ γ₀, so the GM term
  is comparable to Young-Laplace from the start.  Expect a noticeable
  V_crit shift relative to the YL baseline.

## What to expect

| Quantity                          | YL (E_s = 0)              | GM (E_s > 0)                                |
| --------------------------------- | ------------------------- | ------------------------------------------- |
| σ_s before voltage                | γ₀ (constant)             | γ₀ + E_s · ε_eng (grows with stretch)       |
| In-plane stiffness of top fiber   | 0 to lowest order         | E_s (added to bulk)                         |
| Wavelength selection past V_crit  | set by γ₀/μ               | shifted by `E_s` relative to γ₀             |
| V_crit (crease/wrinkle)           | Wang-Zhao 2013            | open — measure and compare                   |

## First result — no-prestretch wrinkle (γ̄ = 5, E_s = 10)

Ran [`staggered_explicit_2D_GM.xml`](staggered_explicit_2D_GM.xml) and
its YL counterpart
[`../2D_young_laplace_experiments/staggered_explicit_2D.xml`](../2D_young_laplace_experiments/staggered_explicit_2D.xml)
side-by-side, both with γ = 20, V(t) = 0.05 · t, no pre-stretch.
Bifurcation thresholds extracted from the top-edge peak-to-peak D_Y
amplitude history:

| amp crosses | YL       | GM       | shift     |
| :---------: | :------: | :------: | :-------: |
| 1·10⁻⁷       | V ≈ 14.88 | V ≈ 15.42 | **+0.54 V** |
| 1·10⁻⁴       | V ≈ 15.42 | V ≈ 15.90 | +0.48 V     |

Both runs subsequently lose elements to inversion past V ≈ 16, but
YL saturates the wrinkle ~2.5× harder than GM at the same voltage
(peak-to-peak |D_Y| = 0.99 vs 6.2·10⁻⁴ at V = 16).  The strain-dependent
surface stress (`E_s · ε_eng`) acts as an extra restoring force on
the top fiber and delays nucleation by ~3.6 % in V_crit relative to
Wang-Zhao Young-Laplace.

Outputs in this folder:

* `gm_vs_yl_no_prestretch.png` — top-edge profile at V ∈ {8, 12, 14, 15, 16}.
* `gm_vs_yl_no_prestretch_amplitude.png` — semilog amp(V) bifurcation diagram.
* `gm_vs_yl_no_prestretch_amplitude.csv` — raw amp(t, V) for both runs.

Reproduce with:

```bash
(cd 2D_young_laplace_experiments && ../../../build/bin/tahoe -f staggered_explicit_2D.xml)
(cd GM_experiments               && ../../../build/bin/tahoe -f staggered_explicit_2D_GM.xml)
(cd GM_experiments && python3 ../scripts/plot_top_edge_evolution.py \
    --gm staggered_explicit_2D_GM.io1.exo \
    --yl ../2D_young_laplace_experiments/staggered_explicit_2D.io1.exo \
    --voltages 8 12 14 15 16 --out gm_vs_yl_no_prestretch.png)
```

## Second result — droplet hypothesis on the fineY mesh

[`droplet_hypothesis_GM.xml`](droplet_hypothesis_GM.xml) and
[`droplet_hypothesis_YL_control.xml`](droplet_hypothesis_YL_control.xml)
use a paper-aspect strip (Lx=80, H=4) refined 4× through-thickness
(Nx=80, Ny=16, dy=0.25, mesh
[`../meshes/bar_2D_fineY.geom`](../meshes/bar_2D_fineY.geom))
with dt=0.005 to test whether the strain-dependent surface law
drives wrinkles toward droplet-like pinching that pure surface tension
cannot reach.

What the runs show:

| V | YL amp | GM amp | ratio | comment |
| ---: | ---: | ---: | ---: | --- |
| 14.0 | 6·10⁻¹¹ | 9·10⁻¹³ | 65×  | both still in linear noise |
| 15.0 | 2·10⁻⁶ | 4·10⁻⁹ | 430× | YL bifurcation visible |
| 15.5 | 1·10⁻³ | 2·10⁻⁶ | 660× | YL fully past V_crit |
| 15.8 | 0.10 | 1·10⁻⁴ | **800×** | clean sinusoidal YL wrinkles |
| 16.0 | 0.49 | 3·10⁻³ | 180× | YL died at element 880 (side, mid-thickness) |

Outcome: a clean **+0.5 V** shift in the wrinkle threshold and an
**~800× amplitude suppression** at matched V (V=15.8) for GM vs YL.
But **no droplet morphology emerged in either case** before the
Lagrangian mesh inverted.  Top-edge profile
([`droplet_fineY_topedge.png`](droplet_fineY_topedge.png)) shows
smooth ~10-period sinusoidal wrinkles in YL up to V=15.9 and the same
smooth pattern emerging in GM starting at V≈16.2.  No necking, no
pinching, no localized severance.

This isn't a falsification of the droplet hypothesis — Lagrangian
quad4 elements simply cannot represent a topological pinch-off.  Both
runs cap out at V ≈ 16 from element inversion at a side / top corner
where the side roller BC pins D_X while the top wants to fold.  4×
refinement in y compared to the original 80×4 paper mesh did *not*
push the V_max ceiling outward — it just relocated the failing
element.  A genuine droplet test would need adaptive remeshing
(unavailable in Tahoe) or a mesh-free method (MPM, SPH, ALE,
phase-field).

Artifacts:

* `droplet_fineY_amplitude.png` — bifurcation diagram (semilog),
  GM curve sits ~0.5 V to the right of YL.
* `droplet_fineY_topedge.png` — top-edge profiles at V ∈ {15.5, 15.8,
  15.9, 16.0, 16.2}, makes the GM suppression visually obvious.
* `droplet_fineY.png` — same in plot_top_edge_evolution.py's per-V
  panel style.
* `droplet_fineY_amplitude.csv` — raw amp(t, V) for both runs.

## Out of scope here (blocked on #54 Phase 5)

* 3D GM (`SimoQ1P0_3D_Surface` needs `E_s`).
* Monolithic GM (`dielectric_elastomer_Q1P0Elastocapillary` needs
  the strain-dependent surface stress integrated into its coupled
  tangent).

## Running

```bash
cd GM_experiments
../../../build/bin/tahoe -f staggered_explicit_2D_GM.xml
```

Compare the resulting `*.io1.exo` against the corresponding YL run in
`../2D_young_laplace_experiments/` — same mesh, same BCs, same time grid, only the
surface law differs.

## Reference

Element-level validation of the analytical residual and tangent lives
under [benchmark_XML/level.4/surface_tension/verify_GM_surface.py](../../../benchmark_XML/level.4/surface_tension/verify_GM_surface.py)
and the integrated benchmark in
[benchmark_XML/level.4/surface_tension/gm_strip/](../../../benchmark_XML/level.4/surface_tension/gm_strip/).
