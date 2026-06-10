# Pinched-cylinder investigation — run index

**How to read the results:**
- `*.io0.exo` — open in ParaView (displacement field `D_X/D_Y/D_Z`, plus `EQ_PLASTIC_STRAIN` when plastic).
- `*_log.txt` — lines `[RKShell-react] node ux uy uz Rx Ry Rz` → **crush = |ux| (col 3)**, **force = |Rx| (col 6)**;
  and `[RKShell-energy] KE max|v|` for runs with the energy diagnostic.
- `fig18_*.txt` — two columns `ux  Rx` (crush, force; both negative).
- All crush runs: nt=160×nz=52 cylinder (R=300, L=600, t=3), E=3e3, ν=0.3, Y0=24.3, H=300, support 2.4, penalty stab.

---

## THE KEY RUNS (this investigation)

| run (deck) | exo / log | what it is | result |
|---|---|---|---|
| **`_dyn.xml`** | `_dyn.io0.exo` / `dyn_log.txt` | **BASELINE penalty-config crush** (H=300), + KE diagnostic | clean 4-lobe butterfly, propagates (localize 292→18), **quasi-static (max\|v\|~1× load rate)**, force 381→2570 (12→84mm) |
| **`_h0.xml`** | `_h0.io0.exo` / `h0_log.txt` | crush with **H=0 (perfect plasticity)** | force plateaus 348→756 → **~paper scale**. Isolates: 4× is the hardening |
| **`_elastic.xml`** | `_elastic.io0.exo` / `el_log.txt` | crush **elastic (yield=0)** | force 801→3818 (plasticity-uncapped). Kinematics is not the direct culprit |
| **`_tu.xml`** | `_tu.io0.exo` / `tu_log.txt` | crush **with thickness update** | only ~10% lower (0.90 ratio @84mm) — real but NOT the 4× fix (folds bending-dominated) |

## VALIDATION BENCHMARKS

| run (deck) | exo / log | what it is | result |
|---|---|---|---|
| `uniaxial.xml` / `_uni2.xml` | `_uni2.io0.exo` / `uni_log.txt` | **uniaxial J2 strip** (base material test) | **= 1.0× analytical** → base membrane/J2/kinematics correct (use _uni2: dt=0.5) |
| `pinched_cylinder.xml` | `pinched_cylinder.io0.exo` | **linear pinched cylinder** (ref 1.8248e-5) | converges with quadratic+stab_bending=20: 2.40→1.35× (nt 24→80) — **element validated** |
| `scordelis_lo.xml` | `scordelis_lo.io0.exo` | Scordelis-Lo roof (ref -0.292) | **PASS -0.2910** |
| `_lpc.xml` | `_lpc.io0.exo` | linear pinched cyl sweep scratch (last = nt80) | (convergence study working file) |

## STABILIZATION SWEEPS (showed 4× is insensitive → not the stabilization)

| run | log | result |
|---|---|---|
| `_sb5.0.xml` `_sb10.0.xml` | `sb5.0_log.txt` `sb10.0_log.txt` | stab_bending sweep (vs 20) → force ~unchanged |
| `_sm5.0.xml` `_sm10.0.xml` | `sm5.0_log.txt` `sm10.0_log.txt` | stab_membrane sweep (vs 20) → force ~unchanged |

## EARLIER DETOURS (didn't pan out — kept for the record)
These were Eq.33 / SCNI / bending-Taylor experiments that I later reverted (they gave the dimple/wrinkle/sawtooth
that you flagged). The element source has since been reverted to the penalty config.

| run | exo | what it was |
|---|---|---|
| `_test1.xml` `_test2_line.xml` `_test3_pinch.xml` | `_test1.io0.exo` … | penalty-free Eq.33 + support 2.4 (line hourglassed 51 lobes; pinch dimpled) |
| `_eq33line.xml` `_eq33pinch.xml` | `_eq33*.io0.exo` | Eq.33 + damping (line wrinkled, pinch localized) |
| `_scni_*.xml` | `_scni_*.io0.exo` | SCNI base experiments (axial wrinkle) |
| `_natline _natpatch _nomemb` | `_nat*.io0.exo` | natural-Taylor + load-footprint experiments |

## Self-test / unit harness (not crush runs)
- `patch_selftest.xml` + `KLSHELL_SELFTEST=1` → flat-patch stabilization mode-energy unit test.
- `gen_strip_geom.py`, `gen_patch_geom.py`, `generate_pinched_geom.py` → geometry generators.
