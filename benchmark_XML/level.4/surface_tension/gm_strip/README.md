# #54 Gurtin-Murdoch surface elasticity — 2D plane-strain strip

A small (8x2 quad4) plane-strain strip with surface elasticity on the top
edge.  Verifies the Gurtin-Murdoch (1975) **strain-dependent** surface
stress

```
sigma_s(eps_eng) = gamma_0 + E_s * eps_eng
eps_eng = (L - L_0) / L_0
```

against the Young-Laplace baseline (`E_s = 0`).  Reduces to the classical
Young-Laplace case by setting `E_s = 0`, which keeps existing surface_tension
XMLs bit-identical (verified by `diag_Q1P0_surface_stretch.xml` smoke test).

## Files

| File                          | Purpose                                          |
| ----------------------------- | ------------------------------------------------ |
| `generate_strip.py`           | builds the `strip.geom` mesh                     |
| `gm_strip_implicit.xml`       | nonlinear_HHT, gamma=2, E_s=10                   |
| `gm_strip_explicit.xml`       | central_difference, gamma=2, E_s=10              |
| `yl_baseline_implicit.xml`    | nonlinear_HHT, gamma=2, E_s=0  (Young-Laplace)   |
| `yl_baseline_explicit.xml`    | central_difference, gamma=2, E_s=0               |
| `verify_gm_strip.py`          | runs all four XMLs and cross-checks them         |
| `strip.geom`                  | generated mesh (27 nodes, 16 elements)           |

## Geometry and BCs

```
        D_X = -0.2           top: side set 1, surface_tension here
        D_Y free                    │ │ │ │ │ │ │ │
              ┌─────────────────────┴─┴─┴─┴─┴─┴─┴─┴─────────────┐
              │ . . . . . . . . . . . . . . . . . . . . . . . . │  D_X = +0.2
   x = 0 ───► │ . . . . . . . . . . . . . . . . . . . . . . . . │ ◄── x = 8
              │_________________________________________________│
                              D_Y = 0  (substrate-bonded bottom)
```

* Plane strain, neo-Hookean Q1P0 mixed (RG_split_general, kappa=1000.67, mu=1).
* Side stretch +/- 0.2 ramped over t in [0, 20], held to t_f = 40 (5% strain).
* Surface stress lives only on the top edge (side set 1).

## How to run

```bash
python3 generate_strip.py            # build mesh
python3 verify_gm_strip.py           # runs all four XMLs, prints checks
```

Set `TAHOE=/path/to/tahoe` to override the binary path; the default looks
in `../../../../build/bin/tahoe`.

## What the verifier checks

1. **Implicit / explicit consistency**: the top-edge length agrees between
   the two integrators to better than 1e-6 in either E_s configuration.
2. **GM vs YL kinematic signature**: with E_s > 0 the surface is stiffer,
   so under lateral stretch the top sags less; mean top-y differs by ~6e-4
   between the two runs for E_s = 10 vs 0.
3. **Analytical GM law**: the measured top-edge engineering strain
   reproduces `sigma_s = gamma_0 + E_s * eps_eng` to machine precision.
4. **Backward compatibility**: YL baseline (E_s=0) yields
   `sigma_s = gamma_0` exactly.

The JAX residual/stiffness verifier `../verify_GM_surface.py` covers the
element-level math separately.

## References

* Gurtin, M.E. and Murdoch, A.I. (1975), *A continuum theory of elastic
  material surfaces*. Arch. Rational Mech. Anal. 57, 291-323.
* Henann, D.L. and Bertoldi, K. (2014), *Modeling of elasto-capillary
  phenomena*. Soft Matter 10, 709-717.
