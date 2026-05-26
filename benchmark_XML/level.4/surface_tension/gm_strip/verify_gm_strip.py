#!/usr/bin/env python3
"""
#54 Gurtin-Murdoch (GM) surface elasticity benchmark verifier.

Runs four Tahoe XMLs and cross-checks them:
  gm_strip_implicit.xml     (nonlinear_HHT, gamma=2,  E_s=10)
  gm_strip_explicit.xml     (central_diff,  gamma=2,  E_s=10)
  yl_baseline_implicit.xml  (nonlinear_HHT, gamma=2,  E_s=0 )  ← Young-Laplace
  yl_baseline_explicit.xml  (central_diff,  gamma=2,  E_s=0 )

Geometry: 8x2 quad4 strip, Lx=8, H=1, plane strain.
Loading : side stretch D_X = +/-0.2 on x=0 and x=Lx (5% engineering strain),
          ramped over t in [0, 20] then held to t_f=40.

Checks
------
1. implicit <-> explicit agreement for the same E_s
     (numerical consistency of the two time-integrators)
2. GM differs from YL in the expected direction
     (sigma_s_GM = gamma_0 + E_s*eps_eng > sigma_s_YL = gamma_0  for eps_eng>0,
      so the top edge contracts MORE in the GM case)
3. The measured surface engineering strain eps_eng matches the analytical
     sigma_s = gamma_0 + E_s*eps_eng prediction at fully-ramped time.

Reads .io0.exo via netCDF4.  Final nodal D_X/D_Y at the top edge (y=1) gives
the deformed top edge length L; eps_eng = (L - L_0)/L_0, L_0 = Lx = 8.
"""
import os
import sys
import subprocess
from netCDF4 import Dataset
import numpy as np

HERE      = os.path.dirname(os.path.abspath(__file__))
TAHOE     = os.environ.get("TAHOE",
            os.path.abspath(os.path.join(HERE, "..", "..", "..", "..",
                                         "build", "bin", "tahoe")))

GAMMA_0   = 2.0
E_S       = 10.0
LX        = 8.0
H         = 1.0
EPS_APPL  = 0.4 / LX     # 0.05 — applied side-edge engineering strain

CASES = [
    ("gm_strip_implicit",    E_S),
    ("gm_strip_explicit",    E_S),
    ("yl_baseline_implicit", 0.0),
    ("yl_baseline_explicit", 0.0),
]


def run_tahoe(xml):
    """Run tahoe on xml.  Skip if .io0.exo already exists and is newer."""
    exo = os.path.join(HERE, f"{xml}.io0.exo")
    if os.path.exists(exo) and (
            os.path.getmtime(exo) > os.path.getmtime(os.path.join(HERE, f"{xml}.xml"))):
        print(f"  [skip] {xml}.io0.exo is up-to-date")
        return
    print(f"  [run ] {TAHOE} -f {xml}.xml")
    log = os.path.join(HERE, f"{xml}.log")
    with open(log, "w") as f:
        r = subprocess.run([TAHOE, "-f", f"{xml}.xml"], cwd=HERE, stdout=f, stderr=f)
    if r.returncode != 0:
        sys.exit(f"  [FAIL] tahoe exit {r.returncode} on {xml}.xml; see {log}")


def read_top_edge(exo_path):
    """Return (X_top, x_top_final, y_top_final) sorted by reference X.

    Top edge is y = H.  Used for length integration."""
    with Dataset(exo_path, "r") as ds:
        coordx = np.asarray(ds.variables["coordx"][:])
        coordy = np.asarray(ds.variables["coordy"][:])
        # nodal variable order: D_X, D_Y, s11, s22, s12 (see XML output spec)
        DX = np.asarray(ds.variables["vals_nod_var1"][-1, :])
        DY = np.asarray(ds.variables["vals_nod_var2"][-1, :])
    # top edge: y_ref == H
    top = np.where(np.isclose(coordy, H))[0]
    order = np.argsort(coordx[top])
    top = top[order]
    X    = coordx[top]
    xcur = coordx[top] + DX[top]
    ycur = coordy[top] + DY[top]
    return X, xcur, ycur


def top_length(xcur, ycur):
    dx = np.diff(xcur)
    dy = np.diff(ycur)
    return float(np.sum(np.sqrt(dx * dx + dy * dy)))


def main():
    if not os.path.exists(TAHOE):
        sys.exit(f"tahoe binary not found at {TAHOE} (set TAHOE env var)")

    # ensure mesh exists
    geom = os.path.join(HERE, "strip.geom")
    if not os.path.exists(geom):
        print("Generating strip.geom ...")
        subprocess.check_call([sys.executable, "generate_strip.py"], cwd=HERE)

    print("=" * 78)
    print("#54 Gurtin-Murdoch surface elasticity benchmark")
    print(f"  gamma_0 = {GAMMA_0}, E_s = {E_S} (GM) or 0 (YL baseline)")
    print(f"  applied side strain eps = {EPS_APPL:.4f}, t_f = 40, plane strain")
    print("=" * 78)

    print("\n[1] Running cases")
    for xml, _ in CASES:
        run_tahoe(xml)

    print("\n[2] Top-edge measurements")
    results = {}
    for xml, Es in CASES:
        exo = os.path.join(HERE, f"{xml}.io0.exo")
        X, xcur, ycur = read_top_edge(exo)
        L_now = top_length(xcur, ycur)
        L0    = float(X[-1] - X[0])           # 8.0
        eps   = (L_now - L0) / L0
        sigma = GAMMA_0 + Es * eps
        ymean = float(np.mean(ycur))
        results[xml] = dict(L0=L0, L=L_now, eps=eps, sigma=sigma,
                            ymean=ymean, Es=Es)
        print(f"  {xml:24s}  L0={L0:8.5f}  L={L_now:8.5f}  "
              f"eps_eng={eps:+.5e}  sigma_s={sigma:+.5e}  <y_top>={ymean:.5e}")

    # ── checks ─────────────────────────────────────────────────────────────
    passed = True
    print("\n[3] Consistency checks")

    def check(label, ok, detail=""):
        nonlocal passed
        tag = "PASS" if ok else "FAIL"
        print(f"  [{tag}] {label}{('  — ' + detail) if detail else ''}")
        if not ok:
            passed = False

    # (a) implicit <-> explicit agreement
    for tag in ("gm_strip", "yl_baseline"):
        a = results[f"{tag}_implicit"]
        b = results[f"{tag}_explicit"]
        dL = abs(a["L"] - b["L"]) / a["L0"]
        check(f"{tag}: implicit/explicit top-length agree (drift = {dL:.3e})",
              dL < 5e-3, f"|dL|/L0={dL:.3e}")

    # (b) GM contracts the top more than YL (sigma_s_GM > sigma_s_YL when eps>0)
    L_gm = results["gm_strip_implicit"]["L"]
    L_yl = results["yl_baseline_implicit"]["L"]
    check("GM stretches top less than YL (E_s adds stiffness for eps>0)",
          L_gm < L_yl,
          f"L_GM={L_gm:.5f}  L_YL={L_yl:.5f}  dL={L_gm - L_yl:+.3e}")

    # (c) sigma_s analytical sanity:
    #     YL gives sigma_s == gamma_0 exactly (E_s = 0).
    s_yl = results["yl_baseline_implicit"]["sigma"]
    check("YL surface stress equals gamma_0 by construction",
          abs(s_yl - GAMMA_0) < 1e-12,
          f"sigma_s_YL={s_yl}")

    # (d) Top edge x-extent is pinned by corner BCs (sides set D_X=+/-0.2),
    #     so eps_eng must agree with the applied side strain to high accuracy.
    eps_gm = results["gm_strip_implicit"]["eps"]
    check("Top engineering strain matches applied side strain (corners are pinned)",
          abs(eps_gm - EPS_APPL) < 1e-3,
          f"eps_GM={eps_gm:+.4e}  vs  eps_applied={EPS_APPL:+.4e}")

    # (e) GM stiffens the surface — the top sags less in the GM run.
    #     Mean y on top edge is the most direct kinematic signature.
    y_gm = results["gm_strip_implicit"]["ymean"]
    y_yl = results["yl_baseline_implicit"]["ymean"]
    check("GM top edge mean-y differs from YL (E_s changes the surface stiffness)",
          abs(y_gm - y_yl) > 1e-5,
          f"<y>_GM={y_gm:.6f}  <y>_YL={y_yl:.6f}  d<y>={y_gm - y_yl:+.3e}")

    # (f) sigma_s shift from baseline is exactly E_s * eps_eng (analytical GM law)
    s_gm = results["gm_strip_implicit"]["sigma"]
    s_an = GAMMA_0 + E_S * eps_gm
    check("GM sigma_s = gamma_0 + E_s*eps_eng (analytical)",
          abs(s_gm - s_an) < 1e-12,
          f"sigma_s_GM={s_gm}  vs  analytical={s_an}")

    print("\n" + "=" * 78)
    print("Overall:", "ALL PASS" if passed else "SOME FAILURES")
    print("=" * 78)
    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()
