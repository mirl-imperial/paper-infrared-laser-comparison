#!/usr/bin/env python3
"""Histogram of the homogenate summed spectra (Supplementary Fig. S9).

Bins each of the 18 brain homogenate profile (summed) spectra, 6 per laser,
into 5-Da m/z bins over m/z 50-1000 and expresses each bin as a percentage
of that replicate's total signal over the same range, so each replicate's
bins sum to 100%. The histogram for each laser is the mean and sample s.d.
(ddof = 1) of each bin over its 6 replicates. The raw summed intensities are
used directly: no baseline subtraction, noise-floor subtraction or peak
picking.

    python homogenate_histogram.py laser_comp_data

The demo writes homogenate_histogram.csv to the working directory, one row
per laser and bin (laser, bin start and end m/z, mean % and s.d. %, number
of replicates), and prints a short summary.

The calculation functions take NumPy arrays and return arrays. File reading
is confined to ``read_sum_spectrum``.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

MZ_MIN, MZ_MAX = 50.0, 1000.0
BIN_WIDTH = 5.0
N_REPLICATES = 6
RAW_DIR = Path("03_brain_homogenate_MS") / "raw_sum_spectra"
OUT_CSV = "homogenate_histogram.csv"

# Laser label -> (data folder, file stem). Longest to shortest pulse duration.
LASERS = {
    "2.94 µm · 3 ns":   ("2.94um_3ns_Montfort",   "Montfort"),
    "2.72 µm · 2 ns":   ("2.72um_2ns_Ivy",        "Ivy"),
    "2.94 µm · 300 ps": ("2.94um_300ps_Innolas",  "Innolas"),
}


def read_sum_spectrum(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Two-column tab-separated profile spectrum: m/z, intensity."""
    df = pd.read_csv(path, sep="\t", header=None, names=["mz", "intensity"])
    return df["mz"].to_numpy(float), df["intensity"].to_numpy(float)


def bin_edges(mz_min: float = MZ_MIN, mz_max: float = MZ_MAX,
              width: float = BIN_WIDTH) -> np.ndarray:
    """Edges of uniform bins of ``width`` from ``mz_min`` to ``mz_max``."""
    n_bins = int(round((mz_max - mz_min) / width))
    return mz_min + width * np.arange(n_bins + 1)


def percent_signal_per_bin(mz: np.ndarray, intensity: np.ndarray,
                           edges: np.ndarray) -> np.ndarray:
    """% of the total signal in [edges[0], edges[-1]) falling in each bin."""
    in_range = (mz >= edges[0]) & (mz < edges[-1])
    total = intensity[in_range].sum()
    per_bin, _ = np.histogram(mz[in_range], bins=edges, weights=intensity[in_range])
    return 100.0 * per_bin / total


def mean_and_sd(percent_by_replicate: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Per-bin mean and sample s.d. (ddof = 1) over replicates (rows)."""
    a = np.asarray(percent_by_replicate, dtype=float)
    return a.mean(axis=0), a.std(axis=0, ddof=1)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("data_dir", type=Path, help="unzipped data record, e.g. laser_comp_data")
    args = p.parse_args()

    edges = bin_edges()
    rows = []
    print(f"{len(edges) - 1} bins of {BIN_WIDTH:g} Da over m/z {MZ_MIN:g}-{MZ_MAX:g}")
    for label, (folder, stem) in LASERS.items():
        pct = np.array([percent_signal_per_bin(
            *read_sum_spectrum(args.data_dir / RAW_DIR / folder / f"{stem}_{i}.txt"), edges)
            for i in range(1, N_REPLICATES + 1)])
        mean, sd = mean_and_sd(pct)
        for lo, hi, m, s in zip(edges[:-1], edges[1:], mean, sd):
            rows.append(dict(laser=label, mz_bin_start=lo, mz_bin_end=hi,
                             mean_pct_of_total_signal=m, sd_pct_of_total_signal=s,
                             n_replicates=len(pct)))
        top = np.argsort(mean)[::-1][:3]
        print(f"  {label:<18} n = {len(pct)}, replicate totals "
              f"{pct.sum(axis=1).min():.4f}-{pct.sum(axis=1).max():.4f}%; largest bins: "
              + ", ".join(f"m/z {edges[i]:g}-{edges[i + 1]:g} ({mean[i]:.2f} ± {sd[i]:.2f}%)"
                          for i in top))

    pd.DataFrame(rows).to_csv(OUT_CSV, index=False)
    print(f"Wrote {OUT_CSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
