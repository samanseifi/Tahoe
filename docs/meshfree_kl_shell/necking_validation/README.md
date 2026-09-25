# Necking of a cylindrical shell — re-verification (Wang & Bazilevs §4.3, Fig 15), 2026-09-22

Run with the merged develop binary (meshfree shell + the feature/55 explicit-integrator changes) at
the paper's settings: R = 10, L = 50, h = 1, E = 189e3, ν = 0.29, Y = 343 + 300 ε_p + 337 (1 − e^(−16.93 ε_p)),
normalized support 2.4, ρ = 7.8e-9 (7800 kg/m³ in mm–MPa units), Δt = 4e-8 s, axial velocity ramped to
30 m/s over 0.5 ms, thickness update on, natural + stress-gradient stabilization, 3 through-thickness
points. Cloud: 8,262 nodes (paper: 8,284), `gen_necking_geom.py`. Deck: `necking_full.xml`
(`monitor_count` = the driven ring; `[RKShell-fig15]` prints U_norm and the total reaction).

Effective stress = |R_z| / A0, A0 = 2πRh = 62.83 mm², against the digitized Fig 15
(`../fig18_validation/reference_data/fig15_necking_effstress_vs_Unorm.csv`). `plot_fig15.py run.log`
regenerates `fig15_compare.png`.

| U_norm [mm] | paper RKPM, 3PT + membrane stab | Ambati 2018 | Tahoe |
|---:|---:|---:|---:|
| 0.5 | 436 | 420 | 419 |
| 1.0 | 489 | 477 | 492 |
| 2.0 | 569 | 555 | 564 |
| 4.0 | 609 | 602 | 613 |
| 6.0 | 598 | 588 | 593 |
| 8.0 | 543 | – | 542 |
| 10.0 | 444 | – | 447 |
| 12.0 | 379 | – | – |

Peak 614 MPa vs the paper's 609. Up to U_norm = 10.7 mm the run follows the paper's stabilized curve
within 1–4 % (and clearly not its unstabilized branch, which is at ~300 MPa by 10.7 mm). The earlier
note that the post-peak branch "tracks the unstabilized curve" (2026-06) no longer holds with the
current code.

**Open item:** at U_norm ≈ 10.9 mm (step ~21,000 of 30,000) the explicit run produces NaN. At the last
valid output the neck at z ≈ 26 mm has max ε_p = 1.28 (growing ~0.11 per 500 steps), t/t0 = 0.45 and a
band width of 3.75 mm; the paper continues to U_norm = 12.2 mm with ε_p up to ~2 (its Fig 14). The
instability sits in the deep-localization regime (thickness < 0.5 h, strain rate in the neck at its
maximum). Candidates to test: the explicit step at reduced thickness / large stretch (rerun with
Δt = 2e-8), the Padé thickness update near its pole guard, and the stabilization scaling with the
current thickness. Not a regression of the merged changes: the pre-merge runs never reached this
point because they localized earlier.
