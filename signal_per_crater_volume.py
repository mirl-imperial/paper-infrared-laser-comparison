#!/usr/bin/env python3
"""Signal per crater volume for the three mouse brain images (Results, Mouse brain imaging).

Reproduces the numbers in the "Mouse brain imaging" subsection of Results
(and the abstract): the mean signal per pixel of each image, the crater
areas and volumes, the signal recovered per µm³ of tissue ablated, and the
scatter of that value between lasers.

    python signal_per_crater_volume.py laser_comp_data

Method
------
1. Common ion set. The ions present in all three images (their 150-peak m/z
   axes, matched within 10 ppm) and also on the consensus ion list (within
   10 ppm). Each m/z of the first image's axis is tested in turn; ppm is
   relative to that m/z. One set shared by all three lasers keeps the sum
   from rewarding whichever laser detects more ions.
2. Mean signal per pixel. For each image, the intensities of the common ions
   (the nearest m/z on that image's axis) are summed over every pixel of the
   image and divided by the number of pixels. The whole image is used, with
   no tissue mask or region selection.
3. Crater area. Mean and sample s.d. (ddof = 1) of the brain SEM crater
   areas measured on the same sections (crater_area_statistics.py).
4. Crater volume = mean area x depth, with ± = area s.d. x depth. The depth is
   the full section depth reached by the craters (2 µm by default).
5. Signal per volume = mean signal per pixel / crater volume, with ± = value x
   (area s.d. / area mean). The image itself contributes no error term: there
   is one image per laser.
6. Scatter = sample coefficient of variation (ddof = 1) of the three
   signal-per-volume values, in %. Fold ranges are max / min over the lasers.

The calculation functions take NumPy arrays and return values. File reading
is done by ``imzml_to_intensity_cube`` (imzml_to_intensity_cube.py),
``read_crater_areas`` (crater_area_statistics.py) and
``read_consensus_mz`` below.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from crater_area_statistics import (LASERS, brain_measurement_csv,
                                    read_crater_areas, summarise_areas)
from imzml_to_intensity_cube import imzml_to_intensity_cube

PPM = 10.0
DEPTH_UM = 2.0
CONSENSUS_CSV = Path("06_consensus_ion_list") / "master_ion_list_laser_comp.csv"
IMZML_DIR = Path("02_MSI_IMZML_data")


# ================================================================ readers ===

def read_consensus_mz(csv_path: Path) -> np.ndarray:
    """The "m/z" column of the consensus ion list."""
    return pd.read_csv(csv_path)["m/z"].to_numpy(float)


# ==================================================== calculation functions ===

def within_ppm(target: float, pool: np.ndarray, ppm: float = PPM) -> bool:
    """True if any value in ``pool`` is within ``ppm`` of ``target`` (relative to target)."""
    return bool(np.min(np.abs(np.asarray(pool) - target)) / target * 1e6 <= ppm)


def common_ion_set(mz_axes: list[np.ndarray], reference_mz: np.ndarray,
                   ppm: float = PPM) -> np.ndarray:
    """m/z values of the first axis found in every other axis and on the reference list."""
    first, *rest = mz_axes
    return np.sort(np.array([m for m in first
                             if all(within_ppm(m, axis, ppm) for axis in rest)
                             and within_ppm(m, reference_mz, ppm)]))


def mean_signal_per_pixel(intensities: np.ndarray, mz: np.ndarray,
                          ion_set: np.ndarray) -> float:
    """Summed intensity of ``ion_set`` over all pixels, divided by the pixel count.

    Parameters
    ----------
    intensities : ndarray (n_pixels, n_peaks)
    mz : ndarray (n_peaks,), m/z axis of ``intensities``
    ion_set : ndarray, target m/z values; each takes the nearest column
    """
    cols = [int(np.argmin(np.abs(mz - t))) for t in ion_set]
    return float(intensities[:, cols].sum() / intensities.shape[0])


def fold_range(values: np.ndarray) -> float:
    """max / min."""
    values = np.asarray(values, dtype=float)
    return float(values.max() / values.min())


def crater_volume(area_mean_um2: float, area_sd_um2: float,
                  depth_um: float = DEPTH_UM) -> tuple[float, float]:
    """(volume, ±) in µm³: area x depth, with the area s.d. carried through."""
    return area_mean_um2 * depth_um, area_sd_um2 * depth_um


def signal_per_volume(signal_per_pixel: float, area_mean_um2: float,
                      area_sd_um2: float, depth_um: float = DEPTH_UM
                      ) -> tuple[float, float]:
    """(counts per µm³, ±): signal per pixel / crater volume, ± from the area s.d."""
    value = signal_per_pixel / (area_mean_um2 * depth_um)
    return value, value * area_sd_um2 / area_mean_um2


def sample_cv_percent(values: np.ndarray) -> float:
    """Sample coefficient of variation (ddof = 1), in %."""
    values = np.asarray(values, dtype=float)
    return float(100 * values.std(ddof=1) / values.mean())


# =================================================================== demo ===

def _sig3(x: float) -> str:
    """Three significant figures in the manuscript's style (7.08e3, 1.24e4)."""
    mantissa, exponent = f"{x:.2e}".split("e")
    return f"{mantissa}e{int(exponent)}"


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("data_dir", type=Path, help="unzipped data record, e.g. laser_comp_data")
    p.add_argument("--depth", type=float, default=DEPTH_UM,
                   help=f"crater depth in µm (default {DEPTH_UM:g})")
    p.add_argument("--ppm", type=float, default=PPM,
                   help=f"m/z matching tolerance in ppm (default {PPM:g})")
    args = p.parse_args()

    images = {}
    for label, stem in LASERS.items():
        (imzml,) = sorted((args.data_dir / IMZML_DIR / stem).glob("*.imzML"))
        intensities, mz, _ = imzml_to_intensity_cube(imzml)
        images[label] = (intensities, mz)
        print(f"{label:<18} {imzml.name}: {intensities.shape[0]} pixels, "
              f"{intensities.shape[1]} peaks")

    reference = read_consensus_mz(args.data_dir / CONSENSUS_CSV)
    ions = common_ion_set([mz for _, mz in images.values()], reference, args.ppm)
    print(f"\nCommon ion set: {len(ions)} ions in all three images and on the "
          f"consensus list (within {args.ppm:g} ppm)")

    rows = []
    for label, stem in LASERS.items():
        intensities, mz = images[label]
        area = summarise_areas(read_crater_areas(brain_measurement_csv(args.data_dir, stem)))
        a, a_sd = area["mean_area_um2"], area["sd_area_um2"]
        sig = mean_signal_per_pixel(intensities, mz, ions)
        vol, vol_sd = crater_volume(a, a_sd, args.depth)
        spv, spv_sd = signal_per_volume(sig, a, a_sd, args.depth)
        rows.append(dict(label=label, sig=sig, n=area["n"], a=a, a_sd=a_sd,
                         vol=vol, vol_sd=vol_sd, spv=spv, spv_sd=spv_sd))

    sig = np.array([r["sig"] for r in rows])
    areas = np.array([r["a"] for r in rows])
    spv = np.array([r["spv"] for r in rows])

    print(f"\nFull precision (crater depth {args.depth:g} µm):")
    for r in rows:
        print(f"  {r['label']:<18} signal/pixel {r['sig']:.4f}  "
              f"area {r['a']:.4f} ± {r['a_sd']:.4f} µm² (n = {r['n']})  "
              f"volume {r['vol']:.4f} ± {r['vol_sd']:.4f} µm³  "
              f"signal/volume {r['spv']:.4f} ± {r['spv_sd']:.4f} counts/µm³")
    print(f"  fold range, signal per pixel:    {fold_range(sig):.4f}")
    print(f"  fold range, mean crater area:    {fold_range(areas):.4f}")
    print(f"  scatter (CV) of signal/volume:   {sample_cv_percent(spv):.4f} %")

    print("\nRounded as in the manuscript:")
    for r in rows:
        print(f"  {r['label']:<18} signal/pixel {_sig3(r['sig'])} counts  "
              f"area {r['a']:.0f} ± {r['a_sd']:.0f} µm²  "
              f"volume {r['vol']:.0f} ± {r['vol_sd']:.0f} µm³  "
              f"signal/volume {r['spv']:.1f} ± {r['spv_sd']:.1f} counts/µm³")
    print(f"  signal per pixel: {fold_range(sig):.1f}-fold range")
    print(f"  crater size:      {fold_range(areas):.1f}-fold range")
    print(f"  signal per volume scatters by {sample_cv_percent(spv):.1f}%")
    return 0


if __name__ == "__main__":
    sys.exit(main())
