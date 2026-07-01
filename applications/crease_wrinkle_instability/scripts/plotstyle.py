#!/usr/bin/env python3
"""Shared figure style for the crease/wrinkle DE paper.

Aesthetic modelled on Yang, Zhao & Sharma (JAM 2017, 84:031008): white
background, full box, NO gridlines, serif/LaTeX math fonts, thin curves,
distinct-colour-per-curve palette (black -> blue -> magenta -> red -> teal),
solid = primary branch / dashed = secondary, and classic open markers for the
finite-element data points.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# distinct-colour-per-curve palette (MATLAB-classic flavour, not matplotlib default)
PALETTE = ["#101010", "#1b3fb0", "#b5179e", "#c01622", "#1e9e84", "#8a6d00"]
# pair for two-curve comparison plots (compression vs stretch etc.)
PAIR = ["#101010", "#c01622"]
# classic open markers for FE data (no fill, dark edge)
FE_MARKERS = ["o", "s", "^", "D", "v"]


def apply_style():
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "text.latex.preamble": r"\usepackage{amsmath}",
        "font.size": 12,
        "axes.labelsize": 13,
        "axes.titlesize": 13,
        "legend.fontsize": 10,
        "lines.linewidth": 1.1,          # thin curves
        "axes.linewidth": 0.7,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,             # full box with inward ticks
        "axes.grid": False,              # NO gridlines
        "savefig.bbox": "tight",
        "savefig.dpi": 300,
        "figure.dpi": 150,
        "axes.prop_cycle": plt.cycler(color=PALETTE),
    })


def fe_marker(ax, x, y, color="#101010", marker="o", label=None, ms=5):
    """Plot FE data as classic open markers (white fill, coloured edge, no line).
    Drawn on top (high zorder) and unclipped so markers at the axis edges stay
    fully visible instead of hiding behind the spines."""
    return ax.plot(x, y, linestyle="none", marker=marker, ms=ms,
                   mfc="white", mec=color, mew=1.1, label=label,
                   zorder=6, clip_on=False)
