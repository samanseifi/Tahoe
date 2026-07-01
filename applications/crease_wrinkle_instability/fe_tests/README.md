# Finite-element tests — organized by case

Surface instabilities of dielectric-elastomer films with surface tension.
Tahoe decks (`stage1_*` static pre-strain, `stage2_*` explicit voltage ramp),
ExodusII results (`*.io1.exo`), and onsets extracted by `common/extract_onsets.py`.
Meshes are shared in `../meshes/`; decks reference them as `../../meshes/...`.

| dir | case | field | pre-strain | layers |
|-----|------|-------|-----------|--------|
| `1_mechanical`            | mechanical wrinkling          | none | compression | single |
| `2_em_no_prestrain`       | electromechanical             | yes  | none (λ=1)  | single |
| `3_em_precompressed`      | electromechanical             | yes  | λ<1         | single |
| `4_em_prestretched`       | electromechanical             | yes  | λ>1         | single |
| `5_mechanical_multilayer` | mechanical wrinkling          | none | compression | multi  |
| `6_em_multilayer`         | electromechanical             | yes  | varies      | multi  |

`common/`  shared runners (`run_sweep2.py`, `drive_sweep2.sh`, `dyn_test.py`),
the extractor (`extract_onsets.py`), `lam_onset.json` (pre-strain sweep spanning
cases 2-4), and `wnl_landau.json` (weakly-nonlinear data, parked for later).
`_archive/`  superseded coarse (80×4) runs and earlier exploratory studies.

Refined production mesh: `bar_2D_fineY.geom` (80×16). Explicit stage-2 needs
DT≤0.005 (CFL: κ=1000 ⇒ dilatational c≈32).
