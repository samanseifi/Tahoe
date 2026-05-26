# Crease & wrinkle instability of a dielectric elastomer (DE) film

Voltage-driven surface instabilities (crease, wrinkle, and buckling) of
a constrained DE film, with **Young-Laplace** and **Gurtin-Murdoch**
surface elasticity.  Tracked under issue
[#55](https://github.com/samanseifi/Tahoe/issues/55); the underlying
strain-dependent surface stress lives in issue
[#54](https://github.com/samanseifi/Tahoe/issues/54).

## What's here

```
crease_wrinkle_instability/
├── meshes/                  bar_2D.geom, bar_3D.geom, plate_3D.geom
├── scripts/                 mesh generators, sweep driver
├── 2D_paper/                Seifi-Park 2016 IJSS replication (YL only)
├── 3D_paper/                3D bar and plate variants
├── prestretch_workflow/     two-stage pre-stretch + voltage XMLs
├── GM_experiments/          NEW — selected runs with E_s ≠ 0
└── references/              papers (PDFs)
```

Every XML in `2D_paper/`, `3D_paper/`, `prestretch_workflow/` and
`GM_experiments/` references the shared meshes via
`geometry_file="../meshes/*.geom"`.  Run each XML from inside its own
subdirectory so the relative path resolves.

## Physics

Replication and extension of:

* **Seifi & Park 2016**, IJSS — "Surface creasing in constrained 2D
  strip" — paper-baseline geometry and BCs.  [PDF](references/seifiIJSS2016.pdf)
* **Wang & Zhao 2013** — analytic V_crease and V_wrinkle thresholds.
* **Gurtin & Murdoch 1975** — strain-dependent surface stress used in
  `GM_experiments/` (`σ_s = γ₀ + E_s · ε_eng`); validated by the
  benchmark suite at [benchmark_XML/level.4/surface_tension/](../../benchmark_XML/level.4/surface_tension).

Material is nearly-incompressible neo-Hookean dielectric with
μ = ε = ρ = 1, κ ≈ 1000.67 (ν ≈ 0.4995), all in the paper's
non-dimensional units.

Elastocapillary number `γ̄ = γ/(μH)` selects the instability mode:

| γ̄ | regime | analytic V_crit (H = 4) |
| ---: | --- | ---: |
| 0    | crease       | V ≈ 4.1  |
| 0.5  | crease       | V ≈ 9.4  |
| 1    | transition   | V ≈ 11.3 |
| 2    | wrinkle      | V ≈ 12.6 |
| 5    | wrinkle      | V ≈ 13.9 |

## Two solver paths

* **Monolithic** (implicit) — single coupled element
  `dielectric_elastomer_Q1P0Elastocapillary`, 3 DOF/node, surface
  tension carried as a material parameter, nonlinear-HHT in time.
  Used by [2D_paper/monolithic_2D.xml](2D_paper/monolithic_2D.xml).
* **Staggered** (explicit or static) — two passes: a `<diffusion>`
  Poisson solve for the potential Ψ (MUMPS, one factorization cached),
  followed by `<updated_lagrangian_Q1P0_surface>` mechanics with the
  Maxwell stress added by the SimoQ1P0 base and Young-Laplace /
  Gurtin-Murdoch force applied on the top side set.
  Used by [2D_paper/staggered_explicit_2D.xml](2D_paper/staggered_explicit_2D.xml)
  and everything in [GM_experiments/](GM_experiments/).

Explicit staggered is the recommended path for chasing the crease
bifurcation — no Newton iterations to stall at threshold.

## Quick start

```bash
# 1. (Re)generate meshes
python3 scripts/generate_bar_2D.py           # 80x4 (paper default)
python3 scripts/generate_bar_3D.py           # 80x4x4 bar
python3 scripts/generate_bar_3D.py --Lx 40 --Lz 40 --Nx 40 --Nz 40 \
    --out ../meshes/plate_3D.geom            # plate

# 2. Run a paper-replication case from its own subdir
(cd 2D_paper && ../../../build/bin/tahoe -f staggered_explicit_2D.xml)

# 3. Or run a GM experiment
(cd GM_experiments && ../../../build/bin/tahoe -f staggered_explicit_2D_GM.xml)
```

Outputs (per run):

* `*.io0.exo` — Ψ (electrostatic) history.
* `*.io1.exo` — `D_X, D_Y, s11, s22, s12` (mechanical) history.
* `*.out`, `*.echo.xml`, `*.valid.xml`, `*.log` — convergence + parsed XML.

Open `.io1.exo` in ParaView, warp by displacement, look at the top
edge for the surface morphology.

## GM vs YL — what to compare

Each pair of XMLs in [GM_experiments/](GM_experiments/) was seeded
from the corresponding YL XML in `2D_paper/` or `prestretch_workflow/`
and differs only by adding `E_s="..."` to `<surface_tension>`.  Useful
diagnostics:

* **V_crit shift** — does `E_s > 0` raise or lower the crease/wrinkle
  threshold?  Sweep `E_s ∈ {0, 5, 10, 20}` at fixed γ.
* **Wavelength selection** — Fourier-transform the top edge profile
  just past nucleation and look for a shift relative to the YL case.
* **Hysteresis under pre-stretch** — the prestretch workflow already
  imposes ε_eng ≈ ε_pre on the top fiber; the GM term then biases
  σ_s by E_s · ε_pre *before* the voltage even turns on.

## Status (open work tracked under #55)

* Monolithic GM extension and 3D GM extension are blocked on issue #54
  Phase 5.  Until that lands, `GM_experiments/` is 2D staggered only.
* The sweep driver `scripts/sweep_prestretch_phase_diagram.sh` was
  reworked to use absolute mesh paths; its prior `../bar_2D.geom`
  resolution was broken from inside `sweep_runs/<tag>/`.
