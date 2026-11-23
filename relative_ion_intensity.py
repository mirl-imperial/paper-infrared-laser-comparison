#!/usr/bin/env python3
"""Relative ion intensities of the consensus ions in brain homogenate (Fig. 4).

Measures the intensity of every consensus ion in each of the 18 homogenate
summed spectra (6 per laser), normalises them, and compares the lasers pair by
pair, giving the fold changes and q-values of Fig. 4 and of the fold-change
columns of Supplementary Table S1.

    python relative_ion_intensity.py laser_comp_data

Method
------
0. Resolving power. The 18 spectra are each scaled to unit total intensity,
   interpolated onto the m/z axis of the first spectrum and averaged. Peaks
   of this mean spectrum with a signal-to-noise ratio of at least 10 are
   found (noise = 1.4826 x the running median, over 301 points, of the
   absolute deviation from a 301-point running median), and R is the median
   over those peaks of m/z / FWHM, the FWHM measured at half height. This
   gives about 17,840, used as R = 17,800 below.
1. Ion intensity. Each summed profile spectrum is used on its own m/z axis
   (no recalibration or resampling). For each consensus m/z, the intensities
   of all data points within one FWHM centred on it, [m - m/2R, m + m/2R]
   with resolving power R = 17,800, are summed. Every ion therefore has a
   measured value in every spectrum, whether or not that laser detected it in
   the consensus list.
2. Normalisation. Each spectrum's ion intensities are scaled to sum to 1,
   then probabilistic quotient normalisation (PQN) divides each spectrum by
   the median of its ratios to the median spectrum. Values are log2
   transformed after adding the 1st percentile of the positive values.
3. Comparison. For each ion, the three pairwise differences of the mean log2
   intensity are tested with Student's t-test using the residual variance
   pooled over all three lasers (18 spectra, 15 degrees of freedom). The
   difference is the log2 fold change.
4. Multiple testing. Benjamini-Hochberg q-values are computed for each pair
   over all consensus ions. A difference is reported when it is at least
   2-fold with q < 0.05.

The calculation functions take NumPy arrays and return values. File reading
is confined to ``read_sum_spectrum`` and ``read_consensus_mz``.

Inputs (deposited data)
-----------------------
    03_brain_homogenate_MS/raw_sum_spectra/<laser folder>/<Montfort|Ivy|Innolas>_1..6.txt
    05_consensus_ion_list/master_ion_list_laser_comp.csv
    05_consensus_ion_list/id_subset_laser_comp.csv (demo only, for the names)

Demo output (written to --out-dir, default the working directory)
-----------------------------------------------------------------
    relative_ion_intensity.csv   one row per consensus ion: m/z, annotation,
                                 mean normalised intensity per laser, and per
                                 pair the log2 fold change, p and q
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.ndimage import median_filter
from scipy.signal import find_peaks, peak_widths

RESOLVING_POWER = 17800        # measured by resolving_power(), rounded
MIN_SNR = 10
NOISE_WINDOW = 301             # points, running-median noise estimate
FOLD_THRESHOLD = 2.0
Q_THRESHOLD = 0.05
N_REPLICATES = 6

# Laser label -> (data folder, file stem). Longest to shortest pulse duration.
LASERS = {
    "2.94 µm · 3 ns":   ("2.94um_3ns_Montfort",   "Montfort"),
    "2.72 µm · 2 ns":   ("2.72um_2ns_Ivy",        "Ivy"),
    "2.94 µm · 300 ps": ("2.94um_300ps_Innolas",  "Innolas"),
}
RAW_DIR = Path("03_brain_homogenate_MS") / "raw_sum_spectra"
ION_DIR = Path("05_consensus_ion_list")


# ================================================================ readers ===

def read_sum_spectrum(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Two-column tab-separated profile spectrum: m/z, intensity."""
    arr = np.loadtxt(path)
    return arr[:, 0], arr[:, 1]


def read_consensus_mz(csv_path: Path) -> np.ndarray:
    """Consensus m/z values (column "m/z") of the master ion list, sorted."""
    return np.sort(pd.read_csv(csv_path)["m/z"].to_numpy(dtype=float))


# ============================================================ calculations ===

def mean_spectrum(spectra: list[tuple[np.ndarray, np.ndarray]]) -> tuple[np.ndarray, np.ndarray]:
    """Mean of spectra scaled to unit total intensity, on the first spectrum's m/z axis.

    ``spectra`` is [(mz, intensity), ...]; the others are linearly
    interpolated onto the first m/z axis. Returns (mz, mean intensity).
    """
    mz = np.asarray(spectra[0][0], dtype=float)
    Y = np.array([np.interp(mz, m, y) for m, y in spectra], dtype=float)
    return mz, (Y / Y.sum(axis=1, keepdims=True)).mean(axis=0)


def resolving_power(mz: np.ndarray, intensity: np.ndarray, min_snr: float = MIN_SNR,
                    noise_window: int = NOISE_WINDOW) -> tuple[float, int]:
    """Median m/z / FWHM over the peaks of a profile spectrum with S/N >= min_snr.

    Noise is 1.4826 x the running median (``noise_window`` points) of the
    absolute deviation from a running median of the same width. FWHM is
    measured at half height (scipy ``peak_widths``) and converted to m/z with
    the local point spacing. Returns (R, number of peaks used).
    """
    mz = np.asarray(mz, dtype=float)
    y = np.asarray(intensity, dtype=float)
    noise = median_filter(np.abs(y - median_filter(y, noise_window)), noise_window) * 1.4826
    peaks, _ = find_peaks(y, distance=3)
    peaks = peaks[y[peaks] >= min_snr * (noise[peaks] + 1e-30)]
    width_pts = peak_widths(y, peaks, rel_height=0.5)[0]
    r = mz[peaks] / (width_pts * np.gradient(mz)[peaks])
    return float(np.median(r)), int(peaks.size)


def window_intensity(mz: np.ndarray, intensity: np.ndarray, centres: np.ndarray,
                     resolving_power: float = RESOLVING_POWER) -> np.ndarray:
    """Summed intensity within one FWHM, [m - m/2R, m + m/2R], of each centre.

    ``mz`` must be sorted ascending. Returns an array of len(centres).
    """
    centres = np.asarray(centres, dtype=float)
    half = centres / (2.0 * resolving_power)
    lo = np.searchsorted(mz, centres - half)
    hi = np.searchsorted(mz, centres + half)
    cum = np.concatenate([[0.0], np.cumsum(intensity, dtype=float)])
    return cum[hi] - cum[lo]


def pqn_log2(X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Total-signal scaling, probabilistic quotient normalisation and log2.

    Parameters
    ----------
    X : ndarray (n_spectra, n_ions) of raw ion intensities

    Returns
    -------
    (normalised intensities, their log2 after adding the 1st percentile of
    the positive normalised values)
    """
    X = np.asarray(X, dtype=float)
    T = X / X.sum(axis=1, keepdims=True)
    ref = np.median(T, axis=0)
    P = T / np.median(T / ref, axis=1, keepdims=True)
    return P, np.log2(P + np.percentile(P[P > 0], 1))


def pairwise_t_tests(log_values: np.ndarray, groups: np.ndarray,
                     pairs: list[tuple[str, str]]) -> dict:
    """Pooled-variance t-test of each ion for each pair of groups.

    Parameters
    ----------
    log_values : ndarray (n_spectra, n_ions), log2 normalised intensities
    groups : ndarray (n_spectra,), group label of each spectrum
    pairs : [(a, b), ...]; the fold change is a over b

    Returns
    -------
    {(a, b): (log2 fold change, p)}, each an array over ions. The residual
    variance is pooled over all groups present in ``groups``.
    """
    log_values = np.asarray(log_values, dtype=float)
    groups = np.asarray(groups)
    labels = list(dict.fromkeys(groups))
    means = {g: log_values[groups == g].mean(axis=0) for g in labels}
    counts = {g: int((groups == g).sum()) for g in labels}
    rss = sum(((log_values[groups == g] - means[g]) ** 2).sum(axis=0) for g in labels)
    dof = len(groups) - len(labels)
    s2 = rss / dof
    out = {}
    for a, b in pairs:
        diff = means[a] - means[b]
        se = np.sqrt(s2 * (1.0 / counts[a] + 1.0 / counts[b]))
        out[(a, b)] = (diff, 2.0 * stats.t.sf(np.abs(diff / se), dof))
    return out


def bh_q_values(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values."""
    return stats.false_discovery_control(np.asarray(p, dtype=float))


def is_reported(log2_fc: np.ndarray, q: np.ndarray,
                fold: float = FOLD_THRESHOLD, q_max: float = Q_THRESHOLD) -> np.ndarray:
    """True where the difference is at least ``fold``-fold with q < ``q_max``."""
    return (np.abs(log2_fc) >= np.log2(fold)) & (np.asarray(q) < q_max)


# =================================================================== demo ===

def _fold_range(fc: np.ndarray) -> str:
    f = 2.0 ** np.abs(fc)
    return f"{f.min():.1f}-{f.max():.1f}-fold"


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("data_dir", type=Path, help="unzipped data record, e.g. laser_comp_data")
    p.add_argument("--out-dir", type=Path, default=Path("."),
                   help="folder for the output CSV (default: working directory)")
    args = p.parse_args()

    centres = read_consensus_mz(args.data_dir / ION_DIR / "master_ion_list_laser_comp.csv")
    X, groups, spectra = [], [], []
    for label, (folder, stem) in LASERS.items():
        for i in range(1, N_REPLICATES + 1):
            mz, inten = read_sum_spectrum(args.data_dir / RAW_DIR / folder / f"{stem}_{i}.txt")
            spectra.append((mz, inten))
            X.append(window_intensity(mz, inten, centres))
            groups.append(label)
    X, groups = np.array(X), np.array(groups)

    r_measured, n_peaks = resolving_power(*mean_spectrum(spectra))
    print(f"Resolving power of the mean spectrum: R = {r_measured:.0f} "
          f"(median over {n_peaks} peaks at S/N >= {MIN_SNR:g}); "
          f"windows use R = {RESOLVING_POWER} (+/- {1e6 / (2 * RESOLVING_POWER):.0f} ppm)")

    P, L = pqn_log2(X)
    labels = list(LASERS)
    pairs = [(labels[1], labels[0]), (labels[2], labels[0]), (labels[2], labels[1])]
    res = pairwise_t_tests(L, groups, pairs)

    ids = pd.read_csv(args.data_dir / ION_DIR / "id_subset_laser_comp.csv")
    names = pd.Series(ids["Molecule"].to_numpy(), index=ids["m/z"].round(6))
    table = pd.DataFrame({"m/z": centres})
    table["Molecule"] = table["m/z"].round(6).map(names)
    for label in labels:
        table[f"{label} mean normalised intensity"] = P[groups == label].mean(axis=0)

    print(f"{len(centres)} consensus ions x {len(X)} summed spectra "
          f"(R = {RESOLVING_POWER}, {N_REPLICATES} per laser)")
    print(f"Reported differences (>= {FOLD_THRESHOLD:g}-fold, q < {Q_THRESHOLD}):")
    for a, b in pairs:
        fc, pv = res[(a, b)]
        q = bh_q_values(pv)
        table[f"log2 FC {a} / {b}"], table[f"p {a} / {b}"], table[f"q {a} / {b}"] = fc, pv, q
        rep = is_reported(fc, q)
        print(f"  {a} vs {b}: {int(rep.sum())} ions "
              f"({int((rep & (fc > 0)).sum())} higher with {a}, "
              f"{int((rep & (fc < 0)).sum())} higher with {b}); "
              f"largest difference {2.0 ** np.abs(fc).max():.2f}-fold")

    print("Annotated classes (fold change, all q < 0.05 unless stated):")
    mol = table["Molecule"].fillna("")
    groups_of_interest = {
        "PA": mol.str.match(r"^PA \d"),
        "cLPA": mol.str.match(r"^cLPA"),
        "fumarate, FA 16:1, NAA, ascorbate, malate": mol.str.contains(
            r"Fumarate|^FA 16:1|acetylaspart|Ascorb|Malate", case=False),
    }
    for name, sel in groups_of_interest.items():
        sel = sel.to_numpy()
        if not sel.any():
            continue
        for a, b in pairs[1:]:
            fc = table.loc[sel, f"log2 FC {a} / {b}"].to_numpy()
            q = table.loc[sel, f"q {a} / {b}"].to_numpy()
            note = "" if (q < Q_THRESHOLD).all() else f" ({int((q >= Q_THRESHOLD).sum())} with q >= 0.05)"
            print(f"  {name:<42} {a} vs {b}: n = {int(sel.sum())}, {_fold_range(fc)}{note}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out = args.out_dir / "relative_ion_intensity.csv"
    table.to_csv(out, index=False)
    print(f"\nWrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
