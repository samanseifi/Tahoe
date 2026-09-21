# Paper reference curves for the necking & pinch decks

Digitized comparison data for the two elasto-plastic benchmarks in

> *A general-purpose meshfree Kirchhoff–Love shell formulation*,
> Engineering with Computers (2025) **41**:1379–1410
> (`ref/s00366-024-01989-x.pdf`)

so that `paper_neck_v2.xml` (necking) and `pinch_s4p4.xml` (pinch) can be
compared against **the same quantities the paper plots**.

## Files

| file | what |
|------|------|
| `digitize_paper_figs.py` | regenerates the two CSVs by colour-segmenting a 10× render of Fig 15 / Fig 18 |
| `fig15_necking_effstress_vs_Unorm.csv` | Fig 15 — effective stress vs `U_norm` |
| `fig18_pinch_force_vs_disp.csv` | Fig 18 — reaction force vs pinch displacement |
| `digitized_overlay.png` | all digitized curves, for visual QC against the PDF |
| `compare_run.py` | parses a run log and overlays the run on the paper curves |

## The two output quantities (must match for a fair comparison)

### Fig 15 — necking of a cylindrical shell (paper §4.3)

- **x = `U_norm`** [mm], a nodal elongation measure
  `U_norm = sqrt( Σ_I u2D_I·u2D_I / NP )` over **all** shell nodes.
- **y = effective stress** [MPa] = `|R_axial| / A0`, where `R_axial` is the
  total reaction on the driven edge and `A0 = 2·π·R·h` is the **undeformed**
  cross-section (`R = 10 mm` mid-surface radius, `h = 1 mm` → `A0 = 62.832 mm²`).

Curves in the CSV (columns):
- `rkpm_membrane_3pt` — paper's recommended result (3 pts thickness + membrane
  stabilization). **This is the target our run should reproduce.**
- `nostab_3pt` — 3 pts, no membrane stab (matches until `U_norm≈8`, then softens).
- `ambati_2018`, `alaydin_2021` — IGA literature references ([33],[41]); plotted
  only to `U_norm≈7` in the paper.
- `1pt_nostab`, `1pt_bending` — under-integrated, *unstable* illustrative cases;
  intentionally noisy, digitized only approximately.

### Fig 18 — pinched elasto-plastic cylinder (paper §4.4)

- **x = pinch displacement** [mm], the prescribed radial displacement (0 → 300).
- **y = reaction force** at the pinch (paper units; with `E=3000 MPa` and mm,
  forces are in N).

Curves in the CSV (columns):
- `rkpm_stab_3pt` — paper's stabilized result. **Target for our run.**
- `nostab_3pt` — no membrane stab; softens badly past ~200 mm (plateaus ≈2150).
- `areias_2010`, `ambati_2018`, `alaydin_2021` — literature references ([79],[33],[41]).

> **Heads-up — the comment in `pinch_s4p4.xml` is wrong about Fig 18.** It says
> "reaction vs displacement 0..150 mm, reaching ~7000–7500 N at 150 mm." In the
> actual figure the x-axis runs to **300 mm**; at 150 mm the force is only the
> small local dip (~600–700), and the curves reach ~5200–6500 only at 300 mm.
> Compare against `fig18_pinch_force_vs_disp.csv`, not that comment.

## How the run emits these quantities

`RKShellT::WriteOutput` prints one parseable line per output step to stdout:

```
necking (Fig 15):  [RKShell-fig15] <Unorm> <Rx> <Ry> <Rz> <|R|>
pinch   (Fig 18):  [RKShell-react] <node> <ux> <uy> <uz> <Rx> <Ry> <Rz>
```

- Necking drives the **axial (Z)** dof → `R_axial = Rz`; effective stress = `|Rz|/A0`.
- Pinch drives the **X** dof → displacement = `|ux|`, reaction = `|Rx|`.

## Workflow

```bash
# from the meshfree_kl_shell/ directory
../../build/bin/tahoe -f paper_neck_v2.xml | tee neck.log
../../build/bin/tahoe -f pinch_s4p4.xml    | tee pinch.log

python3 reference_data/compare_run.py necking neck.log    # -> compare_necking.png + run_necking.csv
python3 reference_data/compare_run.py pinch   pinch.log    # -> compare_pinch.png   + run_pinch.csv
```

`compare_run.py` overlays the black "this run" curve on the dashed paper
references using exactly the axis definitions above.

## Digitization accuracy

Curves recovered by colour segmentation of a 10×-DPI render, calibrated on the
printed axis ticks (necking: x 2→12 mm, y 400→700 MPa; pinch: x 0→300 mm,
y 0→8000). Smooth curves are accurate to roughly the plotted line thickness
(±~10 MPa / ±~50 force units); the noisy 1-pt necking curves are indicative
only. Re-run `digitize_paper_figs.py` (needs `pymupdf numpy pillow`) to rebuild.
