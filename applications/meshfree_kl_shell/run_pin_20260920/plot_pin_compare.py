#!/usr/bin/env python3
"""Overlay pinch runs on digitized Wang & Bazilevs Fig 18.

Series spec  label:run.log[:mode]   mode = react (x=[RKShell-uphys], y=[RKShell-react] coefficient
reaction, the OLD measure), pin (x,y from the legacy element-internal [RKShell-pin] 0 line) or ctrl
(default: x,y = u_mean, lambda_mean of the FIRST penalty_displacement_meshfree controller, read from the
run's .out file next to the .log). KE/W panel uses the same (x,y) pair.
"""
import csv, re, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REF = "/home/samanseifi/codes/tahoe/applications/meshfree_kl_shell/reference_data/fig18_pinch_force_vs_disp.csv"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#8e44ad"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def parse(path, mode):
    lines = open(path).read().split("\n")
    K = [float(re.search(r"KE=(\S+)", l).group(1)) for l in lines if l.startswith("[RKShell-energy]")]
    if mode == "ctrl":
        out = open(path[:-4] + ".out").read().split("\n") if path.endswith(".log") else lines
        P = [l.split() for l in out if l.strip().startswith("penalty_displacement_meshfree: node")]
        first = P[0][2] if P else None
        P = [p for p in P if p[2] == first]
        xy = [(abs(float(p[11])), abs(float(p[14]))) for p in P]
    elif mode == "pin":
        P = [l.split() for l in lines if l.startswith("[RKShell-pin] 0 ")]
        # prefer the interval means ($6,$7) when present, else the instantaneous values
        xy = [(abs(float(p[5])), abs(float(p[6]))) if len(p) >= 7 else (abs(float(p[3])), abs(float(p[4]))) for p in P]
    else:
        R = [l.split() for l in lines if l.startswith("[RKShell-react]")]
        U = [l.split() for l in lines if l.startswith("[RKShell-uphys]")]
        xy = [(abs(float(uu[2])), abs(float(r[5]))) for r, uu in zip(R, U)]
    u, f, ratio = [], [], []
    W, up, fp = 0.0, 0.0, 0.0
    for i, (x, y) in enumerate(xy):
        W += 2 * 0.5 * (y + fp) * (x - up)   # two pinch points
        up, fp = x, y
        u.append(x); f.append(y)
        ratio.append(K[i] / W if W > 0 and i < len(K) else float("nan"))
    return u, f, ratio


def main():
    runs = []
    for a in sys.argv[1:]:
        parts = a.split(":")
        runs.append((parts[0], parts[1], parts[2] if len(parts) > 2 else "ctrl"))
    ref = list(csv.DictReader(open(REF)))
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={"width_ratios": [1.6, 1]})
    for key, style, lab in [("rkpm_stab_3pt", "-", "paper RKPM (target)"),
                            ("ambati_2018", "--", "Ambati 2018"), ("areias_2010", ":", "Areias 2010")]:
        pts = [(float(r["disp_mm"]), float(r[key])) for r in ref if r[key]]
        ax.plot(*zip(*pts), style, color=INK2, lw=1.5, label=lab)
    for c, (lab, path, mode) in zip(SERIES, runs):
        u, f, ratio = parse(path, mode)
        ax.plot(u, f, color=c, lw=2, label=lab)
        bx.plot(u, ratio, color=c, lw=2, label=lab)
    ax.set_yscale("log")
    ax.set_xlabel("physical pinch displacement [mm]", color=INK2)
    ax.set_ylabel("pinch load [N]  (log)", color=INK2)
    ax.set_title("Fig 18 comparison", color=INK, fontsize=11, loc="left")
    bx.axhline(0.05, color=INK2, lw=1, ls="--")
    bx.text(5, 0.06, "quasi-static ~ < 0.05", color=INK2, fontsize=8)
    bx.set_xlabel("physical pinch displacement [mm]", color=INK2)
    bx.set_ylabel("kinetic energy / external work", color=INK2)
    bx.set_title("Inertia share", color=INK, fontsize=11, loc="left")
    for a in (ax, bx):
        a.set_xlim(0, 320); a.grid(True, color=GRID, lw=0.6)
        for s in ("top", "right"): a.spines[s].set_visible(False)
        for s in ("left", "bottom"): a.spines[s].set_color(GRID)
        a.tick_params(colors=INK2)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.tight_layout(); fig.savefig("fig18_pin_compare.png", dpi=150); print("wrote fig18_pin_compare.png")
    # table at reference displacements
    refd = {float(r["disp_mm"]): float(r["rkpm_stab_3pt"]) for r in ref if r["rkpm_stab_3pt"]}
    print("%-22s" % "u[mm]", *["%10.0f" % d for d in (30, 60, 90, 150, 240, 300)])
    print("%-22s" % "paper", *["%10.0f" % refd.get(d, float("nan")) for d in (30, 60, 90, 150, 240, 300)])
    for lab, path, mode in runs:
        u, f, _ = parse(path, mode)
        row = []
        for d in (30, 60, 90, 150, 240, 300):
            i = min(range(len(u)), key=lambda j: abs(u[j] - d)) if u else None
            row.append("%10.0f" % f[i] if i is not None and abs(u[i] - d) < 5 else "%10s" % "-")
        print("%-22s" % lab[:22], *row)


if __name__ == "__main__":
    main()
