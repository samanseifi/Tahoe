# Fresh pinched-cylinder comparison (2026-09-16)

This comparison was regenerated with the current `RKShellT`, rather than reusing the historical
curves in this directory. The raw data and reproducible plot are in:

- `pinched_current_comparison.csv` — linear pinched-cylinder convergence;
- `pinched_plastic_current.csv` — two fresh nonlinear runs through approximately 150 mm physical crush;
- `plot_pinched_paper_comparison.py` — standard-library SVG generator;
- `pinched_paper_comparison.svg` — the resulting vector two-panel comparison;
- `pinched_paper_comparison.png` — dependency-free raster rendering generated locally by the script.
  It is intentionally Git-ignored because the web PR system rejects binary patches; the tracked SVG
  is the reviewable source artifact.

## Linear pinched cylinder

The paper's displacement reference is `1.8248e-5`. The checked-in deck had drifted to **cubic RK
completeness**, although the paper comparison specifies quadratic completeness. That configuration
does not converge: the displacement/reference ratio grows from 3.02 at 468 nodes to 6.87 at 2160
nodes. This is the clearest problem found by rerunning the examples.

With quadratic completeness, the paper-default cubic B-spline kernel, `support_factor=2.5`, natural
stabilization on, and bending penalty off, the current implementation gives:

| nodes | displacement | computed / paper |
|---:|---:|---:|
| 816 | 1.231711e-5 | 0.675 |
| 1,344 | 1.644981e-5 | 0.901 |
| 2,160 | 1.646691e-5 | 0.902 |
| 3,168 | 1.617772e-5 | 0.887 |
| 4,680 | 1.745230e-5 | 0.956 |
| 8,480 | 1.736834e-5 | **0.952** |

The corrected result is within 4.8% at the finest freshly run mesh, but is not monotonically
converged. The benchmark deck now explicitly selects that quadratic configuration.

## Elasto-plastic pinch

> **2026-09-21:** resolved, see `../fig18_validation/README.md` (factor-4 load convention). The
> remainder of this section is the state as of 2026-09-16.


The existing input decks used the obsolete `stab_membrane` penalty and coefficient-space displacement
constraints. A fresh 40-by-15-node run exposed why the older force plots are not directly comparable
to the paper:

1. With ordinary KBCs, a prescribed coefficient displacement of 150 mm produces only about 52 mm of
   **physical** displacement because RK shape functions are non-interpolatory.
2. With collocation, the physical displacement is correct, but `[RKShell-react]` still reports the
   coefficient/generalized reaction. That reaction must be transformed through the transpose of the
   collocation operator before it is a physical load that can be overlaid on Fig. 18.

Therefore the nonlinear panel intentionally labels both curves as **raw generalized reactions**. The
paper information presently available in the repository only states that Fig. 18 is below about 1000
at 150 mm; this is plotted as an upper-bound marker, not fabricated into a digitized paper curve.

## Conclusion

- The **linear pinch is substantially repaired by correcting the deck to quadratic completeness** and
  is within 4.8% of the paper at 8,480 nodes.
- The old nonlinear force plot is not a valid paper comparison because its horizontal axis was a meshfree
  coefficient rather than physical displacement.
- The next code fix should expose a **physical collocation reaction**. Until that transformation is
  implemented, a force-displacement overlay against the paper would be quantitatively misleading.

## Scordelis-Lo cross-check

As a control, the checked-in 225-node Scordelis-Lo example was also rerun. With the new default spline
dilation of 1.0 it gives `min(u_y)=-0.1940`, versus the paper value `-0.292` (33.6% too stiff). The same
mesh gives `-0.2695` with spline dilation 1.4 and `-0.2698` with the legacy Gaussian kernel. Thus the
old documented `-0.291` result is **not reproduced by the current default deck**, and the claim that
the kernel-default change retained that result must be withdrawn pending a controlled mesh/dilation
convergence study.
