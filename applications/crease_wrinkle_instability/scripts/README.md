# Scripts

Mesh generators and the phase-diagram sweep driver.

## Mesh generators

| Script                                       | Default output                | Geometry           |
| -------------------------------------------- | ----------------------------- | ------------------ |
| [`generate_bar_2D.py`](generate_bar_2D.py)   | `../meshes/bar_2D.geom`       | 80 × 4 quad4 strip |
| [`generate_bar_3D.py`](generate_bar_3D.py)   | `../meshes/bar_3D.geom`       | 80 × 4 × 4 hex8 bar |

Both accept `--Lx`, `--H`, `--Lz` (3D only), and `--Nx`, `--Ny`,
`--Nz` for refinement, plus `--out` to override the output path.

To regenerate the plate mesh used by
[../3D_young_laplace_experiments/staggered_explicit_plate_3D.xml](../3D_young_laplace_experiments/staggered_explicit_plate_3D.xml):

```bash
python3 generate_bar_3D.py --Lx 40 --Lz 40 --Nx 40 --Nz 40 \
    --out ../meshes/plate_3D.geom
```

## Sweep driver

[`sweep_prestretch_phase_diagram.sh`](sweep_prestretch_phase_diagram.sh)
sweeps `(ε_pre, γ)` to map the crease/wrinkle V_crit phase diagram for
paper #1.  It generates Stage 1 / Stage 2 XMLs per point in
`sweep_runs/<eps>_<gam>/` and writes a summary CSV.

Usage:

```bash
./sweep_prestretch_phase_diagram.sh                # full grid
./sweep_prestretch_phase_diagram.sh 0.10 2.0       # one point only
TAHOE=/path/to/tahoe ./sweep_prestretch_phase_diagram.sh
```

The script now uses absolute paths for the binary and the mesh, so it
runs correctly from any CWD.

Per-point wall time is ~6-10 min; the full grid is ~3 hr on a single
workstation core.  Output lives under `sweep_runs/` (gitignored).
