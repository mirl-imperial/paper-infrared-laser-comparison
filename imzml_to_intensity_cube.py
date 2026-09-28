#!/usr/bin/env python3
"""Read a continuous-mode centroid imzML/ibd pair into a pixel-by-peak array.

A data-conversion tool whose job is file I/O. It reads one imzML + .ibd pair
from 02_MSI_IMZML_data/ with pyimzML and returns (or, from the command line,
saves to a compressed .npz):

    intensities   float32, shape (n_pixels, n_peaks)
    mz            float64, shape (n_peaks,), the m/z axis shared by every pixel
    coords        int32,   shape (n_pixels, 3), imzML (x, y, z) of each pixel

The deposited images are continuous-mode files, in which every pixel shares
one m/z axis, so the image is a plain two-dimensional table indexed by pixel.
The m/z axis is kept in the order stored in the file (not sorted). The pixel
layout can be recovered from ``coords``.

    python imzml_to_intensity_cube.py \\
        laser_comp_data/02_MSI_IMZML_data/2.94um_3ns_Montfort/2026_04_29_MB1_M2_33_Montfort_25um_5_centroid.imzML \\
        --out 2.94um_3ns_Montfort.npz

``signal_per_crater_volume.py`` imports ``imzml_to_intensity_cube`` to read
the three images.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np


def imzml_to_intensity_cube(imzml_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read a continuous-mode imzML/ibd pair.

    Parameters
    ----------
    imzml_path : path to the .imzML file; the .ibd must sit next to it.

    Returns
    -------
    intensities : ndarray (n_pixels, n_peaks), float32
    mz : ndarray (n_peaks,), float64, shared m/z axis in file order
    coords : ndarray (n_pixels, 3), int32, (x, y, z) per pixel

    Raises ValueError for a processed-mode file (per-pixel m/z axes of
    different lengths).
    """
    from pyimzml.ImzMLParser import ImzMLParser

    parser = ImzMLParser(str(imzml_path))
    coords = np.asarray(parser.coordinates, dtype=np.int32)
    n_pixels = len(coords)
    if n_pixels == 0:
        raise ValueError(f"{Path(imzml_path).name}: no pixels")

    mz0, i0 = parser.getspectrum(0)
    mz0 = np.asarray(mz0, dtype=np.float64)
    intensities = np.empty((n_pixels, mz0.size), dtype=np.float32)
    intensities[0] = i0
    for k in range(1, n_pixels):
        mz_k, i_k = parser.getspectrum(k)
        if len(mz_k) != mz0.size:
            raise ValueError(f"{Path(imzml_path).name}: pixel {k} has {len(mz_k)} peaks, "
                             f"pixel 0 has {mz0.size} (not a continuous-mode file)")
        intensities[k] = i_k
    return intensities, mz0, coords


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("imzml", type=Path, help=".imzML file (.ibd next to it)")
    p.add_argument("--out", type=Path, required=True, help="output .npz file")
    args = p.parse_args()

    intensities, mz, coords = imzml_to_intensity_cube(args.imzml)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.out, intensities=intensities, mz=mz, coords=coords)

    x, y = coords[:, 0], coords[:, 1]
    print(f"Read {args.imzml.name}")
    print(f"  pixels: {intensities.shape[0]}  (x {x.min()}-{x.max()}, y {y.min()}-{y.max()})")
    print(f"  peaks:  {intensities.shape[1]}  (m/z {mz.min():.4f}-{mz.max():.4f})")
    print(f"Wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
