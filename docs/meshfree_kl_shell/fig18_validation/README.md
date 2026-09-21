# Fig 18 pinched elasto-plastic cylinder — resolution (2026-09-20/21)

**Result:** the meshfree KL shell reproduces Wang & Bazilevs Fig 18 (and Ambati 2018, Areias 2010,
Alaydin 2021) once the load convention is accounted for. Every literature curve comes from a
quarter/octant symmetry model (Ambati: "only one quarter of the shell is modeled"; Alaydin: "only one
octant of the full cylinder"), where the load point lies on two symmetry planes and the plotted
"Reaction Force" is the applied nodal force, i.e. **one quarter of the physical pinch load**. Our runs
report the physical point load `lambda`, so compare against the digitized curves times 4:

| u [mm] | paper x4 | Ambati x4 | Areias x4 | shell 160x53 (paper mesh) | Tahoe hex20 solid |
|---:|---:|---:|---:|---:|---:|
| 50  |  874 | 1156 | 1069 |  814 |  852 |
| 100 | 1995 | 2502 | 2356 | 2182 | 2131 |
| 200 | 5036 | 5870 | 5425 | 5580 |    – |
| 300 | 21733 | 23595 | – | 25692 |   – |

`fig18_x4_overlay.png` is the overlay; `fig18_pin_compare.png` the same runs against the unscaled
curves (which is what the earlier "4-9x over-stiff" hunts were comparing against).

## Loading and measurement
The pinch is applied with `penalty_displacement_meshfree` (MFPenaltyDisplacementT): a penalty
enforcement of the PHYSICAL displacement `sum_J Phi_J(x_A) d_J = ubar`, with the force
`Phi_J lambda` on every support coefficient and `lambda` = the physical point reaction (critically
damped dashpot to suppress the penalty-spring ringing in explicit dynamics). A coefficient KBC or the
collocation KBC gives the same load within ~10% (verified on 40x15), but only `lambda` is the
work-conjugate quantity. `plot_pin_compare.py label:run.log[:ctrl|pin|react]` overlays runs.

## Runs (decks kept; .exo/.log/.out not committed)
- `pin40_slow10.xml`, `pin80_slow10.xml`, `pin160_slow10.xml`: 40x15 / 80x27 / 160x53 (paper mesh,
  8480 nodes) full crush to 300 mm, dt=3, 50000-step ramp (quasi-static past ~60 mm, KE/W < 0.05).
- `pin160_nostab.xml` (Taylor stabilization off, to 100 mm: identical to stabilized),
  `pin160_elastic.xml` (yield off, to 36 mm: identical to elastoplastic — plasticity only bites later),
  `pin320_dt15.xml` (320x105, 33,600 nodes, dt=1.5: within 5% of 160x53 -> mesh-converged).
- `static_plastic160.xml`: static load control on the paper mesh, linear to 195 N where yield starts.
- `plate/`: simply supported plate under pressure, w_c and bending stress checks (98-100%).
- `fe_solid/`: independent Tahoe reference, hex20 solid 1/8 model with Simo J2
  (`gen_cyl_hex20.py NT NZ NR [grade]`, `fe_coarse.xml`, `fe_linear.xml` = linear obstacle course
  check, 1.8339e-5 vs 1.8248e-5). Note Tahoe's `linear_function` is `a*x + b`: a = H, b = sigma_y.
- Meshes: `python3 generate_pinched_geom.py NT NZ > pinchedNT.geom` (load nodes 281/301, 1041/1061,
  4161/4241, 16641/16721 for 40/80/160/320).

## What was ruled out along the way (all at the paper mesh)
constraint type and reported force; mesh (converged at 160x53); Taylor stabilization; the through-
thickness plasticity path (elastic == plastic to 35 mm in BOTH the shell and the solid model);
resolution of the point-load stress concentration (graded solid meshes); end boundary conditions
(rigid diaphragm u_x=u_y=0 with u_z free is correct; Ambati Fig 16).
