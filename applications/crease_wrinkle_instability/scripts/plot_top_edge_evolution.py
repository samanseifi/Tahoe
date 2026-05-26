#!/usr/bin/env python3
"""
Plot top-edge evolution for a 2D staggered DE wrinkle/crease run.

Reads `<basename>.io1.exo` (mechanical channel: D_X, D_Y, s11, s22, s12),
finds the y=H top edge, and plots the deformed profile

    (X + D_X, H + D_Y)

at several voltage snapshots.  Optionally overlays two runs (e.g. GM
vs YL) on the same axes.

Usage
-----
    python3 plot_top_edge_evolution.py
        --gm  ../GM_experiments/staggered_explicit_2D_GM.io1.exo
        --yl  ../2D_young_laplace_experiments/staggered_explicit_2D.io1.exo
        --voltages 2 6 10 14 18
        --vdot 0.05
        --out gm_vs_yl_top_edge.png
"""
import argparse
import os
import numpy as np
from netCDF4 import Dataset

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    plt = None


# ── helpers ──────────────────────────────────────────────────────────────────
def read_top_edge_history(exo_path, H=4.0):
    """Return (times, X_top, DX_top[T,N], DY_top[T,N])."""
    with Dataset(exo_path, "r") as ds:
        coordx = np.asarray(ds.variables["coordx"][:])
        coordy = np.asarray(ds.variables["coordy"][:])
        DX_all = np.asarray(ds.variables["vals_nod_var1"][:])  # T x N
        DY_all = np.asarray(ds.variables["vals_nod_var2"][:])
        times  = np.asarray(ds.variables["time_whole"][:])

    top = np.where(np.isclose(coordy, H))[0]
    order = np.argsort(coordx[top])
    top = top[order]
    return times, coordx[top], DX_all[:, top], DY_all[:, top]


def nearest_frames(times, vdot, voltages):
    """For each target V, return the frame index whose t = V/vdot is closest."""
    return [int(np.argmin(np.abs(times - V / vdot))) for V in voltages]


def amplitude(DY):
    """Peak-to-peak top-edge D_Y magnitude (proxy for instability amplitude)."""
    return float(DY.max() - DY.min())


# ── main ─────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gm", required=True, help="path to GM run .io1.exo")
    ap.add_argument("--yl", default=None,  help="path to YL run .io1.exo (optional)")
    ap.add_argument("--H",  type=float, default=4.0)
    ap.add_argument("--vdot", type=float, default=0.05,
                    help="V/t ramp rate; V(t) = vdot * t")
    ap.add_argument("--voltages", type=float, nargs="+",
                    default=[2.0, 6.0, 10.0, 14.0, 18.0])
    ap.add_argument("--out", default="top_edge_evolution.png")
    args = ap.parse_args()

    runs = {"GM": args.gm}
    if args.yl:
        runs["YL"] = args.yl

    # collect histories
    data = {}
    for tag, path in runs.items():
        if not os.path.exists(path):
            raise SystemExit(f"missing exodus file: {path}")
        t, X, DX, DY = read_top_edge_history(path, H=args.H)
        data[tag] = dict(t=t, X=X, DX=DX, DY=DY)

    # pick frames
    for tag, d in data.items():
        d["frames"] = nearest_frames(d["t"], args.vdot, args.voltages)

    # ── numerical summary ───────────────────────────────────────────────────
    print(f"\nTop-edge amplitudes (peak-to-peak D_Y):\n")
    header = "  V   " + "    ".join(f"{tag:>9s}" for tag in data)
    print(header)
    print("  " + "-" * (len(header) - 2))
    for i, V in enumerate(args.voltages):
        row = f"  {V:5.1f}"
        for tag, d in data.items():
            f = d["frames"][i]
            row += f"  {amplitude(d['DY'][f]):9.4e}"
        print(row)

    if data.get("YL") and data.get("GM"):
        print("\nFinal-frame top-edge stats:")
        for tag in ("YL", "GM"):
            d = data[tag]
            f = d["frames"][-1]
            xcur = d["X"] + d["DX"][f]
            ycur = args.H + d["DY"][f]
            print(f"  {tag}: t={d['t'][f]:.2f}  V={args.voltages[-1]:.1f}"
                  f"  L_top={np.sum(np.hypot(np.diff(xcur), np.diff(ycur))):.5f}"
                  f"  max|D_Y|={np.max(np.abs(d['DY'][f])):.4e}")

    # ── plot ────────────────────────────────────────────────────────────────
    if plt is None:
        print("\n(matplotlib unavailable — skipped plot)")
        return

    colors = {"GM": "#d62728", "YL": "#1f77b4"}

    # ── (1) profile snapshots ─────────────────────────────────────────────
    fig, axes = plt.subplots(len(args.voltages), 1,
                             figsize=(8, 2.2 * len(args.voltages)),
                             sharex=True)
    if len(args.voltages) == 1:
        axes = [axes]
    for ax, V, i in zip(axes, args.voltages, range(len(args.voltages))):
        for tag, d in data.items():
            f = d["frames"][i]
            xcur = d["X"] + d["DX"][f]
            ycur = args.H + d["DY"][f]
            ax.plot(xcur, ycur, color=colors.get(tag, "k"),
                    lw=1.3, label=f"{tag}  t={d['t'][f]:.0f}")
        ax.axhline(args.H, ls=":", lw=0.6, color="gray")
        ax.set_ylabel(f"V={V:.1f}\ny_top")
        ax.legend(loc="upper right", fontsize=8, frameon=False)
        ax.grid(True, alpha=0.3)
    axes[-1].set_xlabel("x (deformed)")
    fig.suptitle("Top-edge evolution under voltage ramp", fontsize=11)
    fig.tight_layout()
    fig.savefig(args.out, dpi=140)
    print(f"\nwrote {args.out}")

    # ── (2) amplitude vs V (bifurcation diagram) ──────────────────────────
    fig2, ax2 = plt.subplots(figsize=(6.5, 4.0))
    for tag, d in data.items():
        V_t = args.vdot * d["t"]              # voltage history
        amp = np.array([amplitude(frame) for frame in d["DY"]])
        ax2.semilogy(V_t, np.maximum(amp, 1e-15),
                     color=colors.get(tag, "k"), lw=1.5, label=tag)
    ax2.set_xlabel("V")
    ax2.set_ylabel(r"top-edge $D_Y$ peak-to-peak")
    ax2.set_title("Wrinkle bifurcation: YL vs GM (E_s = 10)")
    ax2.grid(True, which="both", alpha=0.3)
    ax2.legend()
    amp_path = os.path.splitext(args.out)[0] + "_amplitude.png"
    fig2.tight_layout()
    fig2.savefig(amp_path, dpi=140)
    print(f"wrote {amp_path}")

    # ── (3) CSV dump of amplitude(t) ──────────────────────────────────────
    csv_path = os.path.splitext(args.out)[0] + "_amplitude.csv"
    with open(csv_path, "w") as fcsv:
        fcsv.write("t,V," + ",".join(f"amp_{tag}" for tag in data) + "\n")
        # use the shorter history (in case runs differ in length)
        T = min(len(d["t"]) for d in data.values())
        for k in range(T):
            t = next(iter(data.values()))["t"][k]
            row = [f"{t:.4f}", f"{args.vdot * t:.4f}"]
            for tag, d in data.items():
                row.append(f"{amplitude(d['DY'][k]):.6e}")
            fcsv.write(",".join(row) + "\n")
    print(f"wrote {csv_path}")


if __name__ == "__main__":
    main()
