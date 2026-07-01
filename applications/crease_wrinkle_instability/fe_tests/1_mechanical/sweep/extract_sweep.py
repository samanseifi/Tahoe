#!/usr/bin/env python3
"""Extract critical strain eps_c and wrinkle wavelength l/H vs gbar from the
symmetric-compression elastocapillary sweep (sweep_g*.io0.exo).

eps_c -- the surface amplitude A (RMS of the cubic-detrended top profile, /H) is
  tracked vs strain; eps_c is obtained by extrapolating A->0 (fit eps vs A^2 over
  the small-amplitude window, intercept at A^2=0).  This is threshold-independent,
  so gamma=0 recovers Biot rather than overshooting at a finite amplitude.
strain -- eps = (DX_left - DX_right)/L  (symmetric push: both edges move inward).
wavelength -- TWO independent estimates, cross-checked:
  (1) zero-padded FFT peak of the detrended top surface,
  (2) first peak of its autocorrelation.
  Computed a few frames past onset (small amplitude).  Skipped for gamma=0, where
  the Biot problem is scale-free and no wavelength is selected."""
import glob, json, os, re
import numpy as np
from netCDF4 import Dataset

HERE = os.path.dirname(os.path.abspath(__file__))
H = 4.0


def series(f):
    with Dataset(f) as ds:
        cx = np.array(ds["coordx"][:]); cy = np.array(ds["coordy"][:])
        DX = np.array(ds["vals_nod_var1"][:]); DY = np.array(ds["vals_nod_var2"][:])
        t = np.array(ds["time_whole"][:])
    top = np.where(np.isclose(cy, cy.max(), atol=1e-6))[0]
    o = np.argsort(cx[top]); top = top[o]; xr = cx[top]
    L = xr.max() - xr.min()
    rr = np.where(np.isclose(cx, cx.max(), atol=1e-6))[0]
    ll = np.where(np.isclose(cx, cx.min(), atol=1e-6))[0]
    eps = (DX[:, ll].mean(1) - DX[:, rr].mean(1)) / L      # symmetric-safe true strain
    I = (xr > 0.12 * L) & (xr < 0.88 * L)                  # interior (skip edge layers)
    A = np.zeros(len(t)); prof = []
    for k in range(len(t)):
        y = cy[top] + DY[k, top]
        yd = y - np.polyval(np.polyfit(xr, y, 3), xr)      # cubic detrend (remove bulk thickening)
        A[k] = np.sqrt(np.mean(yd[I] ** 2)) / H
        prof.append(yd)
    return eps, A, prof, xr, I, L


def eps_c_extrap(eps, A, lo=0.006, hi=0.05):
    """eps at A->0 by fitting eps vs A^2 over the first small-amplitude window."""
    m = (A > lo) & (A < hi) & (eps > 0.05)
    if m.sum() >= 3:
        p = np.polyfit(A[m] ** 2, eps[m], 1)               # eps = p0*A^2 + p1
        return float(p[1]), "extrap", int(m.sum())
    on = np.where((A > 0.03) & (eps > 0.05))[0]            # fallback: fixed threshold
    return (float(eps[on[0]]) if len(on) else float("nan")), "thr", int(len(on))


def wavelengths(prof, xr, I, k0):
    """FFT-peak and autocorrelation wavelength, median over a few onset frames."""
    xi = xr[I]; dx = xi[1] - xi[0]
    lam_fft, lam_ac = [], []
    for k in range(k0, min(k0 + 5, len(prof))):
        yd = prof[k][I]; yd = yd - yd.mean()
        if np.sqrt((yd ** 2).mean()) < 1e-9:
            continue
        F = np.abs(np.fft.rfft(yd, n=8 * len(yd))); fr = np.fft.rfftfreq(8 * len(yd), d=dx)
        j = 1 + int(np.argmax(F[1:]))
        if fr[j] > 0:
            lam_fft.append(1.0 / fr[j])
        ac = np.correlate(yd, yd, "full")[len(yd) - 1:]; ac = ac / ac[0]
        d = np.diff(ac); up = np.where(d > 0)[0]
        if len(up):
            s = up[0]
            pk = s + int(np.argmax(ac[s:min(s + len(yd), len(ac))]))
            if pk > 0:
                lam_ac.append(pk * dx)
    f = float(np.median(lam_fft)) / H if lam_fft else np.nan
    a = float(np.median(lam_ac)) / H if lam_ac else np.nan
    return f, a


def analyze(f, gbar):
    eps, A, prof, xr, I, L = series(f)
    eps_c, how, npts = eps_c_extrap(eps, A)
    on = np.where((A > 0.02) & (eps > 0.05))[0]
    k0 = on[0] if len(on) else len(eps) - 1
    found = bool(len(on))
    if gbar == 0.0 or not found:
        lf = la = np.nan                                   # no wavelength for scale-free gamma=0
    else:
        lf, la = wavelengths(prof, xr, I, k0)
    return eps_c, lf, la, float(A.max()), how, npts, found


def main():
    rows = []
    for f in sorted(glob.glob(os.path.join(HERE, "sweep_g*.io0.exo"))):
        gbar = int(re.search(r"sweep_g(\d+)", f).group(1)) / 100.0
        eps_c, lf, la, amax, how, npts, found = analyze(f, gbar)
        nn = lambda v: round(float(v), 3) if not np.isnan(v) else None
        rows.append({"gbar": gbar, "eps_c": nn(eps_c), "lH_fft": nn(lf), "lH_ac": nn(la),
                     "amp_max": round(amax, 4), "onset_method": how, "fit_pts": npts,
                     "wrinkled": found})
        wl = "n/a (scale-free)" if gbar == 0.0 else f"FFT={lf:.2f} ac={la:.2f}"
        print(f"  gbar={gbar:>5}: eps_c={eps_c:.3f} [{how},{npts}pt]  l/H: {wl}  (amax={amax:.2f})")
    rows.sort(key=lambda r: r["gbar"])
    json.dump(rows, open(os.path.join(HERE, "sweep_results.json"), "w"), indent=2)
    print(f"\nwrote sweep_results.json ({len(rows)} runs)")


if __name__ == "__main__":
    main()
