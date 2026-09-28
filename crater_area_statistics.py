#!/usr/bin/env python3
"""Brain crater area statistics per laser (Fig. 3; Supplementary Fig. S2).

Summarises the SEM crater areas measured on the imaged mouse brain sections
from the per-crater ImageJ measurement CSVs in the data record:

    05_crater_micrographs/brain_SEM_measurements/<laser stem>_brain.csv

Each CSV has one row per crater and an "Area" column in µm², from an ellipse
fitted to each crater outline in ImageJ (the fitting was done in ImageJ, not
here). For each laser the areas are summarised as n, mean, sample s.d.
(ddof = 1) and the equivalent-circle diameter of the mean area,
2 * sqrt(mean / pi).

    python crater_area_statistics.py laser_comp_data

``summarise_areas`` takes a NumPy array and returns values; file reading is
confined to ``read_crater_areas``.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Laser label -> file stem. Longest to shortest pulse duration.
LASERS = {
    "2.94 µm · 3 ns":   "2.94um_3ns_Montfort",
    "2.72 µm · 2 ns":   "2.72um_2ns_Ivy",
    "2.94 µm · 300 ps": "2.94um_300ps_Innolas",
}
MEASUREMENTS = Path("05_crater_micrographs") / "brain_SEM_measurements"


def read_crater_areas(csv_path: Path) -> np.ndarray:
    """Read the "Area" column (µm²) of one ImageJ measurement CSV."""
    df = pd.read_csv(csv_path)
    df.columns = [c.strip() for c in df.columns]
    return df["Area"].to_numpy(float)


def brain_measurement_csv(data_dir: Path, stem: str) -> Path:
    """Path of the brain SEM measurement CSV for one laser."""
    return Path(data_dir) / MEASUREMENTS / f"{stem}_brain.csv"


def summarise_areas(areas_um2: np.ndarray) -> dict:
    """n, mean, sample s.d. (ddof = 1) and equivalent diameter of crater areas.

    Parameters
    ----------
    areas_um2 : ndarray, shape (n_craters,)
        Per-crater areas in µm².

    Returns
    -------
    dict with ``n``, ``mean_area_um2``, ``sd_area_um2`` and
    ``equiv_diameter_um`` (2 * sqrt(mean / pi), from the mean area).
    """
    a = np.asarray(areas_um2, dtype=float)
    if a.size < 2:
        raise ValueError("need at least two crater areas")
    mean_a = float(a.mean())
    return dict(n=int(a.size), mean_area_um2=mean_a,
                sd_area_um2=float(a.std(ddof=1)),
                equiv_diameter_um=float(2 * np.sqrt(mean_a / np.pi)))


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("data_dir", type=Path, help="unzipped data record, e.g. laser_comp_data")
    args = p.parse_args()

    print(f"{'laser':<18}{'n':>4}{'mean area (µm²)':>18}{'s.d. (µm²)':>12}"
          f"{'equiv. d (µm)':>15}   rounded")
    for label, stem in LASERS.items():
        s = summarise_areas(read_crater_areas(brain_measurement_csv(args.data_dir, stem)))
        print(f"{label:<18}{s['n']:>4}{s['mean_area_um2']:>18.3f}{s['sd_area_um2']:>12.3f}"
              f"{s['equiv_diameter_um']:>15.2f}   "
              f"{s['mean_area_um2']:.0f} ± {s['sd_area_um2']:.0f} µm²")
    return 0


if __name__ == "__main__":
    sys.exit(main())
