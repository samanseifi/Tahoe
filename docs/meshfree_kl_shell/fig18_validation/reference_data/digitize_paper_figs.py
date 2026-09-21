#!/usr/bin/env python3
"""
Digitize Fig. 15 and Fig. 18 of

  "A general-purpose meshfree Kirchhoff-Love shell formulation",
  Engineering with Computers (2025) 41:1379-1410,
  ref/s00366-024-01989-x.pdf

so the necking (Fig. 15) and pinched-cylinder (Fig. 18) decks in this
directory can be compared against the same quantities the paper reports.

Fig. 15 -- necking of a cylindrical shell (Sec. 4.3)
    x: U_norm = sqrt( sum_I u2D_I . u2D_I / NP )   [mm]   (nodal elongation measure)
    y: effective stress = (reaction force on driven edge) / (undeformed
       cross-section area = 2*pi*R*h)               [MPa]

Fig. 18 -- pinched elasto-plastic cylinder (Sec. 4.4)
    x: prescribed radial pinching displacement      [mm]   (ramped 0 -> 300)
    y: reaction force at the pinch                   [N-equivalent units of the paper]

Curves recovered by colour segmentation of a high-DPI render of each plot.
The smooth "consensus" curves (RKPM 3-pt membrane/stabilized + the Ambati,
Alaydin, Areias literature references) are the comparison targets. The
1-pt (under-integrated) curves in Fig. 15 are noisy, illustrative
instability cases and are digitized only approximately.

Requires: PyMuPDF (fitz), numpy, Pillow.  Run from anywhere; paths are
resolved relative to the repo root inferred from this file's location.
"""
import os, csv, numpy as np
from PIL import Image
import fitz

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
PDF  = os.path.join(REPO, "ref", "s00366-024-01989-x.pdf")
W, H = 595.276, 790.866          # page size in points
DPI  = fitz.Matrix(10, 10)


def render(page, clip, out):
    fitz.open(PDF)[page].get_pixmap(matrix=DPI, clip=clip).save(out)


def col_value(mask, c, r0, r1, gap=18):
    """Largest contiguous cluster of masked rows in column c -> its centroid row.
    Rejects thin anti-aliased contamination from neighbouring curves."""
    rows = np.where(mask[r0:r1, c])[0]
    if len(rows) == 0:
        return None
    rows = rows + r0
    groups, cur = [], [rows[0]]
    for v in rows[1:]:
        if v - cur[-1] <= gap:
            cur.append(v)
        else:
            groups.append(cur); cur = [v]
    groups.append(cur)
    best = max(groups, key=len)
    return float(np.median(best))


def extract(img, specs, xcal, ycal, cols, legend=None):
    (cx0, vx0, cx1, vx1) = xcal
    (ry0, vy0, ry1, vy1) = ycal
    X = lambda c: vx0 + (c - cx0) * (vx1 - vx0) / (cx1 - cx0)
    Y = lambda r: vy0 + (r - ry0) * (vy1 - vy0) / (ry1 - ry0)
    R, G, B = img[:, :, 0], img[:, :, 1], img[:, :, 2]
    out = {}
    for name, (cond, r0, r1) in specs.items():
        m = cond(R, G, B)
        if legend is not None:
            m = m & ~legend
        xs, ys = [], []
        for c in cols:
            rr = col_value(m, c, r0, r1)
            if rr is not None:
                xs.append(X(c)); ys.append(Y(rr))
        xs, ys = np.array(xs), np.array(ys)
        if not name.startswith("1pt"):          # keep the noisy 1-pt curves raw
            xs, ys = despike(xs, ys)
        out[name] = (xs, ys)
    return out


def despike(xs, ys, win=9, tol=120.0):
    """Drop points that deviate from a local rolling median by > tol.
    Used only for the smooth (stable) curves; dashed-line gaps occasionally
    yield a stray pixel that the cluster picker latches onto."""
    if len(ys) < win:
        return xs, ys
    keep = np.ones(len(ys), bool)
    for i in range(len(ys)):
        lo, hi = max(0, i - win), min(len(ys), i + win)
        med = np.median(ys[lo:hi])
        if abs(ys[i] - med) > tol:
            keep[i] = False
    return xs[keep], ys[keep]


def resample(out, grid):
    cols = list(out.keys())
    table = []
    for x in grid:
        row = [round(float(x), 3)]
        for n in cols:
            cx, cy = out[n]
            if len(cx) > 3 and cx[0] <= x <= cx[-1]:
                row.append(round(float(np.interp(x, cx, cy)), 1))
            else:
                row.append("")
        table.append(row)
    return ["x"] + cols, table


def save(path, header, table):
    with open(path, "w", newline="") as f:
        w = csv.writer(f); w.writerow(header); w.writerows(table)
    print("wrote", path)


# ----------------------------------------------------------------------- Fig 15
def fig15():
    f = os.path.join(HERE, "_f15.png")
    render(19, fitz.Rect(0.345*W, 0.335*H, 0.83*W, 0.665*H), f)
    img = np.array(Image.open(f).convert("RGB")).astype(int)
    legend = np.zeros(img.shape[:2], bool); legend[1360:1835, 1850:] = True
    specs = {
        "rkpm_membrane_3pt": (lambda R,G,B:(R>205)&(G<70)&(B<70),                                 150,1955),
        "nostab_3pt":        (lambda R,G,B:(R<60)&(G<60)&(B>175),                                  150,1955),
        "ambati_2018":       (lambda R,G,B:(G>160)&(R>100)&(R<195)&(B<95)&((G-R)>25),              150,1955),
        "alaydin_2021":      (lambda R,G,B:(R>210)&(G>165)&(B<130),                                150,1955),
        "1pt_nostab":        (lambda R,G,B:(R>110)&(R<195)&(G<65)&(B>20)&(B<85),                  1550,1955),
        "1pt_bending":       (lambda R,G,B:(G>55)&(G<135)&(R<75)&(B>35)&(B<95)&((G-R)>15),         150,1955),
    }
    out = extract(img, specs,
                  xcal=(549, 2.0, 2664, 12.0),     # x ticks: col 549->2mm, col 2664->12mm
                  ycal=(229, 700.0, 951, 400.0),   # y ticks: row 229->700,  row 951->400 MPa
                  cols=range(120, 2700), legend=legend)
    hdr, tab = resample(out, np.arange(0, 12.01, 0.25))
    hdr[0] = "U_norm_mm"
    save(os.path.join(HERE, "fig15_necking_effstress_vs_Unorm.csv"), hdr, tab)
    os.remove(f)


# ----------------------------------------------------------------------- Fig 18
def fig18():
    f = os.path.join(HERE, "_f18.png")
    render(22, fitz.Rect(0.46*W, 0.02*H, 0.99*W, 0.34*H), f)
    img = np.array(Image.open(f).convert("RGB")).astype(int)
    legend = np.zeros(img.shape[:2], bool); legend[460:960, 183:820] = True
    specs = {
        "areias_2010":   (lambda R,G,B:(R<110)&(G>90)&(G<185)&(B>150),               460,2437),
        "ambati_2018":   (lambda R,G,B:(G>150)&(B<120)&(R<175)&((G-R)>25),            460,2437),
        "alaydin_2021":  (lambda R,G,B:(R>200)&(G>150)&(B<125),                       460,2437),
        "rkpm_stab_3pt": (lambda R,G,B:(R>185)&(G<100)&(B<100),                       460,2437),
        "nostab_3pt":    (lambda R,G,B:(R<70)&(G<70)&(B>150),                         460,2437),
    }
    out = extract(img, specs,
                  xcal=(182, 0.0, 2530, 300.0),    # x ticks: col 182->0mm, col 2530->300mm
                  ycal=(2436, 0.0, 460, 8000.0),   # y ticks: row 2436->0, row 460->8000
                  cols=range(184, 2531), legend=legend)
    hdr, tab = resample(out, np.arange(0, 300.5, 5.0))
    hdr[0] = "disp_mm"
    save(os.path.join(HERE, "fig18_pinch_force_vs_disp.csv"), hdr, tab)
    os.remove(f)


if __name__ == "__main__":
    fig15()
    fig18()
