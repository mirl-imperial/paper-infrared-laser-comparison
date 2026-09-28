#!/usr/bin/env python3
"""Chloride-adduct isotope envelope of m/z 862.65 (Supplementary Note S1, Fig. S7).

Tests the assignment of m/z 862.65 (the green channel of Fig. 2) as
[HexCer 42:1;O3 + 35Cl]- by comparing its measured precursor isotope envelope
with the theoretical one. Chlorine's two isotopes add a large M+2 peak that no
non-chlorinated composition of the same nominal mass can reproduce, so
M+2/M+0 is the discriminating ratio.

    python chloride_isotope_envelope.py laser_comp_data

Method
------
1. Theoretical envelope of the anion C48H93NO9·35Cl- (neutral HexCer 42:1;O3
   plus one Cl, electron mass added): the exact-mass isotope distribution is
   binned to nominal shifts M+0 to M+3, giving each peak's abundance relative
   to M+0, its abundance-weighted m/z and the share of it carried by 37Cl.
2. Theoretical M+2/M+0 of a non-chlorinated composition of the same nominal
   mass, C51H92NO10, as the control.
3. Measured envelope, from each of the 18 homogenate profile (summed)
   spectra (6 per laser). For each peak, a ±0.05 Da window around the
   theoretical m/z is baseline-subtracted (baseline = median of the two
   outermost points at each end) and integrated by the trapezoidal rule. The
   observed m/z is the intensity-weighted mean of the points at or above half
   maximum. Ratios are areas relative to M+0.
4. Summary over replicates: mean and sample s.d. (ddof = 1) of each ratio,
   and the mass error of M+0 in ppm, (observed - theoretical) / theoretical.

The profile spectra are used, not the peak lists, because the peak-list
intensities are averages that distort isotope ratios.

The calculation functions take NumPy arrays (or dicts of formulae) and return
values. File reading is confined to ``read_sum_spectrum``.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

NEUTRAL = {"C": 48, "H": 93, "N": 1, "O": 9}     # HexCer 42:1;O3
ADDUCT_EL = "Cl"
ISOBAR = {"C": 51, "H": 92, "N": 1, "O": 10}     # non-chlorinated control
PEAK_HALFWIDTH_DA = 0.05
N_PEAKS = 4                                      # M+0 .. M+3
E_MASS = 0.00054858                              # electron mass, Da

# (exact mass, natural abundance) per isotope.
ISO = {
    "C":  [(12.0, 0.9893), (13.003355, 0.0107)],
    "H":  [(1.0078250319, 0.999885), (2.0141018, 0.000115)],
    "N":  [(14.0030740052, 0.99632), (15.0001089, 0.00368)],
    "O":  [(15.9949146221, 0.99757), (16.9991315, 0.00038), (17.9991604, 0.00205)],
    "Cl": [(34.96885271, 0.7576), (36.9659026, 0.2424)],
}

# Laser label -> (data folder, file stem). Longest to shortest pulse duration.
LASERS = {
    "2.94 µm · 3 ns":   ("2.94um_3ns_Montfort",   "Montfort"),
    "2.72 µm · 2 ns":   ("2.72um_2ns_Ivy",        "Ivy"),
    "2.94 µm · 300 ps": ("2.94um_300ps_Innolas",  "Innolas"),
}
N_REPLICATES = 6
RAW_DIR = Path("03_brain_homogenate_MS") / "raw_sum_spectra"


# ================================================================ reader ===

def read_sum_spectrum(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Two-column tab-separated profile spectrum: m/z, intensity."""
    arr = np.loadtxt(path)
    return arr[:, 0], arr[:, 1]


# ======================================================= theoretical envelope ===

def element_isotope_dist(el: str, n: int) -> dict:
    """Exact-mass distribution {mass: probability} of n atoms of one element."""
    dist = {0.0: 1.0}
    for _ in range(n):
        nxt: dict = {}
        for m, p in dist.items():
            for dm, pi in ISO[el]:
                nxt[m + dm] = nxt.get(m + dm, 0.0) + p * pi
        dist = {m: p for m, p in nxt.items() if p > 1e-13}
    return dist


def combine_dists(a: dict, b: dict) -> dict:
    """Convolve two exact-mass distributions."""
    out: dict = {}
    for m1, p1 in a.items():
        for m2, p2 in b.items():
            out[m1 + m2] = out.get(m1 + m2, 0.0) + p1 * p2
    return {m: p for m, p in out.items() if p > 1e-13}


def formula_dist(formula: dict) -> dict:
    """Exact-mass distribution of an elemental formula {element: count}."""
    dist = {0.0: 1.0}
    for el, n in formula.items():
        dist = combine_dists(dist, element_isotope_dist(el, n))
    return dist


def nominal_bin(dist: dict, n_peaks: int) -> tuple[dict, dict]:
    """Bin an exact-mass distribution to nominal shifts 0..n_peaks-1.

    Returns ({k: abundance / M+0 abundance}, {k: abundance-weighted mass}).
    """
    base = min(dist)
    ab, wm = {}, {}
    for m, p in dist.items():
        k = int(round(m - base))
        if k < n_peaks:
            ab[k] = ab.get(k, 0.0) + p
            wm[k] = wm.get(k, 0.0) + p * m
    return ({k: ab[k] / ab[0] for k in sorted(ab)},
            {k: wm[k] / ab[k] for k in sorted(ab)})


def theoretical_chloride_envelope(neutral: dict = NEUTRAL, n_peaks: int = N_PEAKS
                                  ) -> tuple[dict, dict, dict]:
    """Theoretical envelope of the [neutral + Cl]- anion.

    Returns
    -------
    rel : {k: abundance relative to M+0}
    cen : {k: m/z of peak M+k (electron mass added)}
    cl37_pct : {k: % of peak M+k carried by 37Cl}
    """
    rest = formula_dist(neutral)
    cl = element_isotope_dist(ADDUCT_EL, 1)
    rel, cen = nominal_bin(combine_dists(rest, cl), n_peaks)
    cen = {k: v + E_MASS for k, v in cen.items()}
    rest_rel, _ = nominal_bin(rest, n_peaks)
    cl37_over_cl35 = cl[max(cl)] / cl[min(cl)]
    cl37_pct = {k: 100.0 * cl37_over_cl35 * rest_rel.get(k - 2, 0.0) / rel[k]
                for k in rel}
    return rel, cen, cl37_pct


def non_chlorinated_isobar_m2_ratio(isobar: dict = ISOBAR) -> float:
    """Theoretical M+2/M+0 of the non-chlorinated control composition."""
    rel, _ = nominal_bin(formula_dist(isobar), 3)
    return float(rel[2])


# ======================================================= measured envelope ===

def integrate_peak(mz: np.ndarray, intensity: np.ndarray, centre: float,
                   halfwidth: float = PEAK_HALFWIDTH_DA) -> tuple[float, float]:
    """Baseline-subtracted trapezoidal area and half-maximum centroid of one peak.

    Returns (area, centroid m/z); NaN where the window holds < 3 points or no signal.
    """
    w = (mz > centre - halfwidth) & (mz < centre + halfwidth)
    if w.sum() < 3:
        return float("nan"), float("nan")
    x, y = mz[w], intensity[w].astype(float)
    base = np.median(np.concatenate([y[:2], y[-2:]]))
    y = np.clip(y - base, 0.0, None)
    trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    area = float(trapz(y, x))
    if y.max() <= 0:
        return area, float("nan")
    top = y >= 0.5 * y.max()
    return area, float((x[top] * y[top]).sum() / y[top].sum())


def measure_envelope(mz: np.ndarray, intensity: np.ndarray,
                     centres: list[float]) -> tuple[np.ndarray, np.ndarray]:
    """Area ratios to M+0 and observed m/z of each peak of one spectrum.

    ``centres`` are the theoretical m/z of M+0, M+1, ... Returns
    (ratios, observed m/z), each of length len(centres).
    """
    areas, cents = zip(*(integrate_peak(mz, intensity, c) for c in centres))
    areas = np.asarray(areas)
    return areas / areas[0], np.asarray(cents)


def summarise_replicates(ratios: np.ndarray, mz_obs_m0: np.ndarray,
                         mz_theo_m0: float) -> dict:
    """Mean and sample s.d. over replicates of the ratios and the M+0 mass error.

    Parameters
    ----------
    ratios : ndarray (n_replicates, n_peaks), M+k / M+0 per replicate
    mz_obs_m0 : ndarray (n_replicates,), observed m/z of M+0
    mz_theo_m0 : theoretical m/z of M+0

    Returns
    -------
    dict with ``ratio_mean``, ``ratio_sd`` (arrays over peaks),
    ``mz_obs_mean``, ``ppm_mean``, ``ppm_sd`` and ``n``.
    """
    ratios = np.asarray(ratios, dtype=float)
    mz_obs_m0 = np.asarray(mz_obs_m0, dtype=float)
    ppm = (mz_obs_m0 - mz_theo_m0) / mz_theo_m0 * 1e6
    return dict(ratio_mean=ratios.mean(axis=0), ratio_sd=ratios.std(axis=0, ddof=1),
                mz_obs_mean=float(mz_obs_m0.mean()), ppm_mean=float(ppm.mean()),
                ppm_sd=float(ppm.std(ddof=1)), n=int(ratios.shape[0]))


# =================================================================== demo ===

def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("data_dir", type=Path, help="unzipped data record, e.g. laser_comp_data")
    args = p.parse_args()

    rel, cen, cl37 = theoretical_chloride_envelope()
    iso_m2 = non_chlorinated_isobar_m2_ratio()
    centres = [cen[k] for k in range(N_PEAKS)]

    print("Per-replicate M+k / M+0 (area ratios) and observed m/z of M+0:")
    ratios, mz_m0 = [], []
    for label, (folder, stem) in LASERS.items():
        for i in range(1, N_REPLICATES + 1):
            mz, inten = read_sum_spectrum(args.data_dir / RAW_DIR / folder / f"{stem}_{i}.txt")
            r, obs = measure_envelope(mz, inten, centres)
            ratios.append(r)
            mz_m0.append(obs[0])
            print(f"  {label:<18} {stem}_{i:<3}"
                  + "".join(f"  M+{k} {r[k]:.4f}" for k in range(1, N_PEAKS))
                  + f"   m/z {obs[0]:.4f}")

    s = summarise_replicates(np.array(ratios), np.array(mz_m0), cen[0])
    print(f"\nPeak   m/z theo   % M+0 obs (n = {s['n']})   % M+0 theo   obs - theo   37Cl share")
    for k in range(N_PEAKS):
        obs, sd, theo = 100 * s["ratio_mean"][k], 100 * s["ratio_sd"][k], 100 * rel[k]
        print(f"M+{k}   {cen[k]:9.4f}   {obs:6.1f} ± {sd:3.1f}{'':12}{theo:6.1f}"
              f"{obs - theo:+13.1f}{cl37[k]:11.0f}%")
    m2, m2_sd = s["ratio_mean"][2], s["ratio_sd"][2]
    print(f"\nM+2/M+0: measured {m2:.2f} ± {m2_sd:.3f}; one chlorine {rel[2]:.2f}; "
          f"non-chlorinated C51H92NO10 {iso_m2:.2f}")
    print(f"Separation of the measured M+2/M+0 from the non-chlorinated value: "
          f"{(m2 - iso_m2) / m2_sd:.0f} x the between-replicate s.d.")
    print(f"Largest |measured - theoretical| over M+0 to M+3: "
          f"{100 * max(abs(s['ratio_mean'][k] - rel[k]) for k in range(N_PEAKS)):.1f} "
          f"percentage points")
    print(f"M+0: theoretical {cen[0]:.4f}, observed {s['mz_obs_mean']:.4f} "
          f"(mean of {s['n']}), error {s['ppm_mean']:+.1f} ± {s['ppm_sd']:.1f} ppm")
    return 0


if __name__ == "__main__":
    sys.exit(main())
