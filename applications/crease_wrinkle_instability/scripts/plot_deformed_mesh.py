#!/usr/bin/env python3
"""
Visualize the deformed 2D mesh from a Tahoe `*.io1.exo` run at several
voltage snapshots, side-by-side for two surface laws (GM and YL).
Useful for spotting non-wrinkle morphology — sharp pinches, beading,
or droplet-like detachment in the post-bifurcation regime.

Reads coords + connectivity + nodal displacements; plots the deformed
quad4 mesh as a wireframe, optionally tinted by the local top-fiber
engineering strain.

Usage
-----
    python3 plot_deformed_mesh.py
        --gm droplet_hypothesis_GM.io1.exo
        --yl droplet_hypothesis_YL_control.io1.exo
        --voltages 12 14 15 16 18
        --vdot 0.05
        --out droplet_morphology.png
"""
import argparse
import os
import numpy as np
from netCDF4 import Dataset

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
except ImportError:
    plt = None


def read_mesh_and_history(exo_path):
    """Return (X, Y, conn(0-indexed Mx4), DX[T,N], DY[T,N], times)."""
    with Dataset(exo_path, "r") as ds:
        X    = np.asarray(ds.variables["coordx"][:])
        Y    = np.asarray(ds.variables["coordy"][:])
        conn = np.asarray(ds.variables["connect1"][:]) - 1  # 1- → 0-indexed
        DX   = np.asarray(ds.variables["vals_nod_var1"][:])
        DY   = np.asarray(ds.variables["vals_nod_var2"][:])
        t    = np.asarray(ds.variables["time_whole"][:])
    return X, Y, conn, DX, DY, t


def deformed_edges(X, Y, conn, DX_f, DY_f):
    """Return Mx2x2 array of edges (start[xy], end[xy]) for plotting."""
    xc = X + DX_f
    yc = Y + DY_f
    edges = []
    for e in conn:                                       # CCW order
        for i in range(4):
            a, b = e[i], e[(i + 1) % 4]
            edges.append(((xc[a], yc[a]), (xc[b], yc[b])))
    return np.asarray(edges)


def nearest_frames(times, vdot, voltages):
    return [int(np.argmin(np.abs(times - V / vdot))) for V in voltages]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gm", required=True)
    ap.add_argument("--yl", required=True)
    ap.add_argument("--vdot", type=float, default=0.05)
    ap.add_argument("--voltages", type=float, nargs="+",
                    default=[12, 14, 15, 16, 18])
    ap.add_argument("--out", default="deformed_morphology.png")
    ap.add_argument("--linewidth", type=float, default=0.15)
    ap.add_argument("--y_aspect", type=float, default=1.0,
                    help="multiplier for the y axis vs equal aspect (>1 stretches y)")
    ap.add_argument("--xwin", type=float, nargs=2, default=None,
                    help="x window to zoom in on (e.g. --xwin 20 60)")
    args = ap.parse_args()

    runs = {"YL": args.yl, "GM": args.gm}

    data = {}
    for tag, path in runs.items():
        if not os.path.exists(path):
            raise SystemExit(f"missing exodus file: {path}")
        X, Y, conn, DX, DY, t = read_mesh_and_history(path)
        data[tag] = dict(X=X, Y=Y, conn=conn, DX=DX, DY=DY, t=t,
                         frames=nearest_frames(t, args.vdot, args.voltages))

    nV = len(args.voltages)
    fig, axes = plt.subplots(nV, 2, figsize=(11.5, 2.0 * nV),
                             sharex=True, sharey=True)
    if nV == 1:
        axes = np.atleast_2d(axes)

    H = float(data["YL"]["Y"].max())
    if args.xwin:
        xmin, xmax = args.xwin
    else:
        xmin = -0.5
        xmax = float(data["YL"]["X"].max()) + 0.5

    for i, V in enumerate(args.voltages):
        for col, tag in enumerate(("YL", "GM")):
            d = data[tag]
            f = d["frames"][i]
            edges = deformed_edges(d["X"], d["Y"], d["conn"], d["DX"][f], d["DY"][f])
            ax = axes[i, col]
            lc = LineCollection(edges, linewidths=args.linewidth,
                                colors=("#1f77b4" if tag == "YL" else "#d62728"))
            ax.add_collection(lc)
            ax.axhline(H, ls=":", lw=0.5, color="gray")
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(-0.5, H + 1.5)
            # equal aspect, optionally exaggerated in y
            ax.set_aspect(1.0 / args.y_aspect)
            ax.set_title(f"{tag}  V = {V:.1f}  (t = {d['t'][f]:.1f})", fontsize=9)
            ax.grid(True, alpha=0.25)
        axes[i, 0].set_ylabel("y")
    axes[-1, 0].set_xlabel("x")
    axes[-1, 1].set_xlabel("x")

    fig.suptitle("Deformed mesh: YL (left) vs GM E_s=10 (right)", fontsize=11)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"wrote {args.out}")

    # ── side-by-side top-edge profile zoom (one panel per voltage) ────────
    fig2, axs2 = plt.subplots(nV, 1, figsize=(10, 1.6 * nV), sharex=True)
    if nV == 1:
        axs2 = [axs2]
    colors = {"YL": "#1f77b4", "GM": "#d62728"}
    for ax2, V, i in zip(axs2, args.voltages, range(nV)):
        for tag in ("YL", "GM"):
            d = data[tag]
            f = d["frames"][i]
            top = np.where(np.isclose(d["Y"], H))[0]
            order = np.argsort(d["X"][top])
            top = top[order]
            xt = d["X"][top] + d["DX"][f, top]
            yt = d["Y"][top] + d["DY"][f, top]
            ax2.plot(xt, yt, color=colors[tag], lw=1.0,
                     label=f"{tag} (t={d['t'][f]:.0f})")
        ax2.axhline(H, ls=":", lw=0.6, color="gray")
        ax2.set_ylabel(f"V={V:.1f}\ny_top")
        ax2.legend(loc="upper right", fontsize=8, frameon=False)
        ax2.grid(True, alpha=0.3)
        if args.xwin:
            ax2.set_xlim(args.xwin)
    axs2[-1].set_xlabel("x (deformed)")
    fig2.suptitle("Top-edge profiles: YL vs GM", fontsize=11)
    fig2.tight_layout()
    top_out = os.path.splitext(args.out)[0] + "_topedge.png"
    fig2.savefig(top_out, dpi=150)
    print(f"wrote {top_out}")

    # Also print a per-frame top-edge curvature / pinch indicator.
    print("\nTop-edge geometry indicators (per voltage):")
    print(f"  {'V':>6}  {'tag':>3}  {'amp':>10}  {'max(-y′′)':>10}  {'L_top':>9}")
    for i, V in enumerate(args.voltages):
        for tag in ("YL", "GM"):
            d = data[tag]
            f = d["frames"][i]
            top = np.where(np.isclose(d["Y"], H))[0]
            order = np.argsort(d["X"][top])
            top = top[order]
            xt = d["X"][top] + d["DX"][f, top]
            yt = d["Y"][top] + d["DY"][f, top]
            amp = float(yt.max() - yt.min())
            # discrete second derivative (proxy for curvature/pinch)
            if len(yt) >= 3:
                d2 = np.gradient(np.gradient(yt, xt), xt)
                max_neg_curv = float(np.max(-d2))   # large positive = sharp upward pinch
            else:
                max_neg_curv = 0.0
            L_top = float(np.sum(np.hypot(np.diff(xt), np.diff(yt))))
            print(f"  {V:>6.1f}  {tag:>3}  {amp:>10.3e}  "
                  f"{max_neg_curv:>10.3e}  {L_top:>9.4f}")


if __name__ == "__main__":
    main()
