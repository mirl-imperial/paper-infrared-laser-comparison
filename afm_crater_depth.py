#!/usr/bin/env python3
"""Crater depth from atomic force microscopy (Supplementary Fig. S3; Table 1).

Levels the AFM height map of a single 2.72 µm · 2 ns crater in intact mouse
brain to the undisturbed tissue surface and reports the crater floor, rim and
peak-to-valley height. The ~2 µm floor depth is the crater depth used for the
crater volumes and the signal per crater volume of Table 1
(``signal_per_crater_volume.py``, ``--depth``).

    python afm_crater_depth.py laser_comp_data

Method
------
1. Height map. Layer 0 (HeightRetrace, m) of the Igor binary wave; the scan
   size (m) is read from the wave note, and the pixel size is the scan size
   divided by the number of rows.
2. Levelling. A plane z = a x + b y + c is fitted by least squares to the
   pixels of the undisturbed surface. These are found iteratively: after each
   fit, pixels whose residual lies more than 2.5 x the robust s.d.
   (1.4826 x MAD) from the median residual are excluded, so the crater and
   its rim do not tilt the plane. The loop stops when the pixel set no longer
   changes (at most 6 iterations). Heights are then referenced to the median
   of the undisturbed surface, so 0 is the tissue surface.
3. Summary. Deepest point, mean height of the crater floor (pixels deeper
   than -1.5 µm, where the floor forms a plateau before the sloped wall),
   highest rim point, and peak-to-valley height of the levelled map. The
   line profile of panel c is the middle row of the map after rotating it
   90° anticlockwise to the orientation of the figure.

The calculation functions take NumPy arrays and return values. File reading
is confined to ``read_afm_height``, which needs the ``igor2`` package.

Inputs (deposited data)
-----------------------
    06_AFM_topography/2.72um_2ns_Ivy_brain_25um_scan.ibw   single crater
    06_AFM_topography/2.72um_2ns_Ivy_brain_90um_scan.ibw   crater array
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

MAD_THRESHOLD = 2.5
MAX_ITER = 6
FLOOR_THRESHOLD_UM = -1.5
AFM_DIR = Path("06_AFM_topography")
SINGLE_SCAN = "2.72um_2ns_Ivy_brain_25um_scan.ibw"
ARRAY_SCAN = "2.72um_2ns_Ivy_brain_90um_scan.ibw"


# ================================================================ reader ===

def read_afm_height(path: Path) -> tuple[np.ndarray, float]:
    """HeightRetrace layer (m) of an Asylum .ibw file and its pixel size (m)."""
    from igor2 import binarywave

    wave = binarywave.load(str(path))["wave"]
    height = wave["wData"][:, :, 0].astype(np.float64)
    note = wave["note"].decode("latin1", "ignore")
    scan = float(next(l for l in note.splitlines() if l.startswith("ScanSize:")).split(":")[1])
    return height, scan / height.shape[0]


# ============================================================ calculations ===

def level_to_undisturbed_surface(height: np.ndarray, n_iter: int = MAX_ITER,
                                 mad_threshold: float = MAD_THRESHOLD
                                 ) -> tuple[np.ndarray, np.ndarray]:
    """Plane-level a height map to its undisturbed surface (step 2).

    Returns (levelled height, same units as ``height``, with 0 at the median
    of the undisturbed surface; boolean mask of the undisturbed pixels).
    """
    height = np.asarray(height, dtype=float)
    ny, nx = height.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    mask = np.ones_like(height, dtype=bool)
    plane = np.zeros_like(height)
    for _ in range(n_iter):
        A = np.column_stack([xx[mask], yy[mask], np.ones(mask.sum())])
        coef, *_ = np.linalg.lstsq(A, height[mask], rcond=None)
        plane = coef[0] * xx + coef[1] * yy + coef[2]
        resid = height - plane
        med = np.median(resid[mask])
        mad = np.median(np.abs(resid[mask] - med)) * 1.4826
        new_mask = np.abs(resid - med) < mad_threshold * mad
        converged = new_mask.sum() == mask.sum()
        mask = new_mask
        if converged:
            break
    levelled = height - plane
    return levelled - np.median(levelled[mask]), mask


def crater_depth_summary(levelled_um: np.ndarray,
                         floor_threshold_um: float = FLOOR_THRESHOLD_UM) -> dict:
    """Floor, rim and peak-to-valley of a levelled height map in µm (step 3)."""
    z = np.asarray(levelled_um, dtype=float)
    floor = z[z < floor_threshold_um]
    return dict(deepest_um=float(z.min()),
                floor_mean_um=float(floor.mean()) if floor.size else float("nan"),
                rim_um=float(z.max()),
                peak_to_valley_um=float(z.max() - z.min()))


def middle_row_profile(levelled_um: np.ndarray) -> np.ndarray:
    """Height along the middle row after a 90° anticlockwise rotation (panel c)."""
    disp = np.rot90(np.asarray(levelled_um), k=1)
    return disp[disp.shape[0] // 2, :]


# =================================================================== demo ===

def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("data_dir", type=Path, help="unzipped data record, e.g. laser_comp_data")
    args = p.parse_args()

    for name in (SINGLE_SCAN, ARRAY_SCAN):
        height, px = read_afm_height(args.data_dir / AFM_DIR / name)
        levelled, mask = level_to_undisturbed_surface(height)
        z = levelled * 1e6
        print(f"{name}: {z.shape[0]} x {z.shape[1]} px, {px * 1e9:.1f} nm/px, "
              f"undisturbed-surface s.d. {z[mask].std() * 1000:.0f} nm")
        if name == SINGLE_SCAN:
            s = crater_depth_summary(z)
            prof = middle_row_profile(z)
            print(f"  whole map:    floor mean {s['floor_mean_um']:.2f} µm "
                  f"(pixels below {FLOOR_THRESHOLD_UM} µm), deepest {s['deepest_um']:.2f} µm, "
                  f"rim {s['rim_um']:+.2f} µm, peak-to-valley {s['peak_to_valley_um']:.2f} µm")
            print(f"  line profile: deepest {prof.min():.2f} µm, highest {prof.max():+.2f} µm")
        else:
            print(f"  levelled range {z.min():.2f} to {z.max():+.2f} µm")
    return 0


if __name__ == "__main__":
    sys.exit(main())
