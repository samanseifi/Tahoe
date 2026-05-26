# Pre-stretch workflow — stretch first, then drive the voltage

Two-stage runs that pre-stretch the film mechanically, then apply the
voltage on the pre-stretched state.  Used for the phase-diagram axis
of paper #1 (V_crit as a function of ε_pre and γ̄) and for voltage-
driven buckling experiments.

Workflow:

* **Stage 1** — static (or quasi-static) mechanical pre-stretch via
  side-edge BCs, no voltage.  Restart file (`*.rsN_of_N`) saved at
  end of run.
* **Stage 2** — restart from Stage 1, declare both Ψ and displacement
  fields, ramp `γ` over `t_0`, then ramp `V`.  Explicit central-
  difference for the mechanics, MUMPS-cached diffusion for Ψ.

## Stage 1 — pre-stretch only

| File                                                              | Setup |
| ----------------------------------------------------------------- | ----- |
| [prestretch_static_2D.xml](prestretch_static_2D.xml)               | minimal static stretch, ε_pre = 5% |
| [prestretch_staggered_static.xml](prestretch_staggered_static.xml) | both fields, static, ε_pre = 5% (paired with `voltage_after_stretch_staggered.xml`) |
| [prestretch_staggered_static_eps10.xml](prestretch_staggered_static_eps10.xml) | static stretch ε_pre = 10% |
| [prestretch_staggered_static_eps20.xml](prestretch_staggered_static_eps20.xml) | static stretch ε_pre = 20% |
| [prestretch_eps10_stage1.xml](prestretch_eps10_stage1.xml)         | explicit Stage 1 at ε_pre = 10% |
| [prestretch_eps10_stage1_dt001.xml](prestretch_eps10_stage1_dt001.xml) | same at smaller dt (Newton robustness) |
| [prestretch_eps20_stage1.xml](prestretch_eps20_stage1.xml)         | explicit Stage 1 at ε_pre = 20% |
| [stage1_plainQ1P0_eps20.xml](stage1_plainQ1P0_eps20.xml)           | plain Q1P0 (no surface tension) Stage 1, ε_pre = 20% |

## Stage 2 — voltage on pre-stretched state

| File                                                              | Restarts from                                | Notes |
| ----------------------------------------------------------------- | -------------------------------------------- | ----- |
| [voltage_after_stretch_staggered.xml](voltage_after_stretch_staggered.xml) | `prestretch_staggered_static.rs2000of2000`     | HHT mech, γ̄ = 0.5 |
| [stage2_explicit_eps10.xml](stage2_explicit_eps10.xml)             | `prestretch_eps10_stage1_dt001.rs4000of4000`   | explicit, ε_pre = 10%, γ = 2 |
| [stage2_explicit_eps20.xml](stage2_explicit_eps20.xml)             | `prestretch_eps20_stage1.rs2000of2000`         | explicit, ε_pre = 20%, γ = 2 |
| [voltage_eps20_gamma0_stage2.xml](voltage_eps20_gamma0_stage2.xml) | `prestretch_eps20_stage1.rs2000of2000`         | ε_pre = 20%, γ = 0 control |

## Monolithic alternative

[prestretched_monolithic_2D.xml](prestretched_monolithic_2D.xml) bundles
the stretch and voltage into one schedule on the coupled element — no
restart, but slower and more fragile near V_crit (see issue #53 for the
monolithic Newton stall behavior).

## Sweep driver

[`../scripts/sweep_prestretch_phase_diagram.sh`](../scripts/sweep_prestretch_phase_diagram.sh)
loops over (ε_pre, γ) pairs, writing per-point Stage 1 / Stage 2 XMLs
into `sweep_runs/<tag>/` and parsing the final `D_Y` spread back into
a CSV summary.

## Running

```bash
cd prestretch_workflow
../../../build/bin/tahoe -f prestretch_eps20_stage1.xml      # Stage 1
../../../build/bin/tahoe -f stage2_explicit_eps20.xml        # Stage 2 (uses Stage 1 restart)
```

The Stage 1 restart files (`*.rs*`) land in the prestretch_workflow/
directory, which is where Stage 2 looks for them.
