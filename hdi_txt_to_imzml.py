"""
Convert a Waters HDImaging (Maldichrom) peak-picked pixel export (.txt) into a
centroid imzML + .ibd pair, readable by any imzML-consuming tool.

Background
----------
Waters HDImaging's own "export to imzML" (HDImaging v1.4 / ImzmlConverter
v0.5) re-reads the raw file and writes full profile spectra. It does not use
the peak-picked data, and the resulting files are very large (tens of GB and
more).

HDImaging's Process tab (Maldichrom) separately produces a tab-delimited .txt
export: one shared m/z peak list for the whole image (the N most intense
peaks, chosen during processing), with each pixel's intensity at each of those
peaks. That is centroided data and is much smaller. This script converts the
.txt export directly into a centroid imzML.

It performs no centroiding, peak picking or binning of its own and applies no
smoothing, filtering, resampling or intensity transform. The peak list and the
per-pixel intensities are written unchanged. The m/z values are written in
the order of the export, which is HDImaging's peak rank order, not increasing
m/z, so the file does not declare an increasing-m/z scan.

.txt format (Waters HDImaging v1.4 Maldichrom output)
------------------------------------------------------
Line 1: blank, or a free-text title (e.g. "Default file"); unused.
Line 2: summed/reference-spectrum row, index "0"; ignored.
Line 3: peak index header (1..N); used only to count the N peaks.
Line 4: the N peak m/z values: the single m/z axis shared by every pixel.
Line 5+: one row per pixel:
    col 0       : pixel index (1-based, sequential in acquisition order)
    col 1       : stage x position (mm). It resets to zero at the start of
                  every scan line, so the gap between consecutive resets is
                  the raster width. The script uses this to cross-check
                  --width and warns on a mismatch (see
                  _count_and_detect_width). --width is still required,
                  because this reset pattern is not guaranteed for every
                  export.
    col 2       : stage y position (mm), constant along a scan line; ignored.
    col 3..3+N-1: the N intensities, aligned to the line-4 m/z list.
    last 2 cols : MassLynx function / scan-number bookkeeping. These are not
                  pixel coordinates, although they look like them: the
                  function number is constant and the scan number repeats
                  col 0. Do not use them for pixel position.

Pixel (x, y) coordinates are reconstructed from the acquisition order (col 0)
and the raster geometry, using the scan-order conventions of imzML's
scan-settings terms (--scan-direction, --line-scan-direction, --scan-pattern,
--scan-type). The defaults match the standard HDImaging raster (top down, one
way, horizontal lines, left to right). --width is required and is
cross-checked against the col-1 reset pattern when available. --height is
inferred from pixel count / width if not given, and checked for an exact
match.

After pyimzML has written the files, the imzML header is edited in three
places only: the pixel size (IMS:1000046/47) is added to the scan settings,
the "increasing m/z scan" parameter (MS:1000093) that pyimzML writes by
default is removed, and the <run> id is set to the output file's base name
(pyimzML writes the full output path there).

Usage
-----
    python hdi_txt_to_imzml.py INPUT.txt OUTPUT_STEM \\
        --width 296 --pixel-size 25 --polarity negative

Produces OUTPUT_STEM.imzML + OUTPUT_STEM.ibd (continuous mode, centroid
spectra). Requires numpy and pyimzML.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import time

import numpy as np
from pyimzml.ImzMLWriter import ImzMLWriter


def pixel_index_to_xy(i0, width, height, scan_direction, line_scan_direction,
                      scan_pattern, scan_type):
    """i0: 0-based pixel index in acquisition order. Returns 1-based (x, y)."""
    if scan_type == "horizontal_line":
        line_len, n_lines = width, height
    elif scan_type == "vertical_line":
        line_len, n_lines = height, width
    else:
        raise ValueError(f"unknown scan_type: {scan_type}")

    line_no = i0 // line_len          # 0-based line number, in acquisition order
    pos_in_line = i0 % line_len       # 0-based position along the line

    base_reversed = line_scan_direction in ("line_right_left", "line_bottom_up")
    if scan_pattern == "meandering" and line_no % 2 == 1:
        reversed_this_line = not base_reversed
    else:
        reversed_this_line = base_reversed
    pos = (line_len - 1 - pos_in_line) if reversed_this_line else pos_in_line

    line_reversed = scan_direction == "bottom_up"
    line = (n_lines - 1 - line_no) if line_reversed else line_no

    if scan_type == "horizontal_line":
        x, y = pos, line
    else:
        x, y = line, pos
    return x + 1, y + 1


def _count_and_detect_width(f):
    """Count data rows and cross-check the raster width from column 1.

    Column 1 (stage x position, mm) resets to zero at the start of every scan
    line, so the gap between consecutive resets is the line width. Returns
    (n_pixels, detected_width_or_None); None if the resets are not perfectly
    uniform, in which case the caller relies on --width alone.
    """
    n_pixels = 0
    prev_t = None
    resets = []
    for line in f:
        if not line.strip():
            continue
        t = float(line.split("\t", 3)[1])
        if prev_t is not None and t < prev_t - 1e-9:
            resets.append(n_pixels)
        prev_t = t
        n_pixels += 1
    if len(resets) < 2:
        return n_pixels, None
    gaps = {resets[i + 1] - resets[i] for i in range(len(resets) - 1)}
    if len(gaps) != 1:
        return n_pixels, None
    return n_pixels, resets[0]


def convert(txt_path, out_path, width, pixel_size, height=None,
            pixel_size_y=None, polarity="negative", scan_direction="top_down",
            line_scan_direction="line_left_right", scan_pattern="one_way",
            scan_type="horizontal_line", limit=None):
    t0 = time.time()

    def open_data_rows():
        """Open the export and read its header; returns (file, n_peaks, mzs).

        The file is reopened for each pass (count, then write) so it is
        never held in memory.
        """
        f = open(txt_path, "r", encoding="utf-8", errors="strict")
        f.readline()  # line 1: blank, or a title; unused
        f.readline()  # line 2: summed / reference row; unused
        row_idx = f.readline().rstrip("\n").split("\t")
        row_mz = f.readline().rstrip("\n").split("\t")
        n_peaks = len([c for c in row_idx[3:] if c != ""])
        mzs = np.array([float(v) for v in row_mz[3:3 + n_peaks]], dtype=np.float64)
        assert len(mzs) == n_peaks
        return f, n_peaks, mzs

    f, n_peaks, mzs = open_data_rows()
    print(f"n_peaks={n_peaks}  mz range={mzs.min():.4f}-{mzs.max():.4f}")

    if limit is not None:
        n_pixels = limit
        for _ in range(limit):
            f.readline()
    else:
        n_pixels, detected_width = _count_and_detect_width(f)
        if detected_width is not None and detected_width != width:
            print(
                f"WARNING: --width {width} does not match the raster width "
                f"({detected_width}) detected from column 1 (stage x "
                f"position) resetting to zero at the start of each scan line. "
                f"Check --width; {detected_width} is more likely to be "
                f"correct."
            )
    f.close()

    if height is None:
        assert n_pixels % width == 0, (
            f"{n_pixels} pixels is not an exact multiple of --width {width}; "
            "pass the correct raster width (or --height explicitly)."
        )
        height = n_pixels // width
    else:
        assert width * height == n_pixels, (
            f"width*height ({width}*{height}={width*height}) != pixel count "
            f"({n_pixels}); check --width/--height."
        )
    print(f"n_pixels={n_pixels}  width={width}  height={height}")

    pixel_size_y = pixel_size if pixel_size_y is None else pixel_size_y

    max_x = max_y = 0
    f, _, _ = open_data_rows()
    with f, ImzMLWriter(
        out_path,
        mode="continuous",
        spec_type="centroid",
        polarity=polarity,
        scan_direction=scan_direction,
        line_scan_direction=line_scan_direction,
        scan_pattern=scan_pattern,
        scan_type=scan_type,
    ) as writer:
        i0 = 0
        for line in f:
            if not line.strip():
                continue
            if i0 >= n_pixels:
                break
            cols = line.rstrip("\n").split("\t")
            intens = np.array(cols[3:3 + n_peaks], dtype=np.float64)
            x, y = pixel_index_to_xy(
                i0, width, height, scan_direction, line_scan_direction,
                scan_pattern, scan_type,
            )
            max_x, max_y = max(max_x, x), max(max_y, y)
            writer.addSpectrum(mzs, intens, (x, y, 1))
            i0 += 1
            if i0 % 200000 == 0:
                print(f"  {i0} pixels written ({time.time() - t0:.0f}s)")

    print(f"Done: {n_pixels} pixels, max_x={max_x}, max_y={max_y}, "
          f"{time.time() - t0:.0f}s")
    _finalise_header(out_path, pixel_size, pixel_size_y)


def _finalise_header(out_path, pixel_size_x, pixel_size_y):
    """Edit the imzML header written by pyimzML.

    1. Add the pixel size (IMS:1000046/47), which pyimzML does not write, to
       the scanSettings block.
    2. Remove the "increasing m/z scan" parameter (MS:1000093) that pyimzML
       writes by default: the m/z array is in the export's peak rank order.
    3. Set the <run> id to the output file's base name instead of the full
       output path that pyimzML writes.
    """
    stem = os.path.splitext(out_path)[0]   # the same stem pyimzML uses
    imzml_path = stem + ".imzML"
    with open(imzml_path, "r", encoding="ISO-8859-1") as f:
        xml = f.read()

    marker = 'name="max count of pixels y"'
    idx = xml.index(marker)
    line_end = xml.index("\n", idx)
    insertion = (
        f'\n      <cvParam cvRef="IMS" accession="IMS:1000046" name="pixel size (x)" '
        f'value="{pixel_size_x}" unitCvRef="UO" unitAccession="UO:0000015" unitName="micrometer"/>'
        f'\n      <cvParam cvRef="IMS" accession="IMS:1000047" name="pixel size y" '
        f'value="{pixel_size_y}" unitCvRef="UO" unitAccession="UO:0000015" unitName="micrometer"/>'
    )
    xml = xml[:line_end] + insertion + xml[line_end:]

    xml, n_removed = re.subn(
        r'\n[ \t]*<cvParam cvRef="MS" accession="MS:1000093" name="increasing m/z scan"/>',
        "", xml)
    assert n_removed == 1, f"expected one increasing-m/z parameter, found {n_removed}"

    run_id = os.path.basename(stem)
    xml, n_run = re.subn(r'(<run defaultInstrumentConfigurationRef="IC1" id=")[^"]*(")',
                         lambda m: m.group(1) + run_id + m.group(2), xml, count=1)
    assert n_run == 1, "could not find the <run> element"

    with open(imzml_path, "w", encoding="ISO-8859-1") as f:
        f.write(xml)
    print(f"Header: pixel size {pixel_size_x} x {pixel_size_y} um, run id "
          f'"{run_id}", increasing-m/z parameter removed ({imzml_path})')


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input_txt", help="HDImaging Maldichrom .txt pixel export")
    p.add_argument("output_stem", help="output path without extension "
                                       "(writes OUTPUT_STEM.imzML + .ibd)")
    p.add_argument("--width", type=int, required=True,
                   help="raster width in pixels (required; cannot be "
                        "recovered reliably from the .txt file)")
    p.add_argument("--height", type=int, default=None,
                   help="raster height in pixels (default: inferred from "
                        "pixel count / width, and checked)")
    p.add_argument("--pixel-size", type=float, required=True,
                   help="pixel size in micrometres (x)")
    p.add_argument("--pixel-size-y", type=float, default=None,
                   help="pixel size in micrometres (y), default: same as --pixel-size")
    p.add_argument("--polarity", choices=["positive", "negative"], required=True)
    p.add_argument("--scan-direction", choices=["top_down", "bottom_up"],
                   default="top_down")
    p.add_argument("--line-scan-direction",
                   choices=["line_left_right", "line_right_left",
                            "line_top_down", "line_bottom_up"],
                   default="line_left_right")
    p.add_argument("--scan-pattern", choices=["one_way", "meandering"],
                   default="one_way")
    p.add_argument("--scan-type", choices=["horizontal_line", "vertical_line"],
                   default="horizontal_line")
    p.add_argument("--limit", type=int, default=None,
                   help="only convert the first N pixels (for testing)")
    args = p.parse_args()

    convert(args.input_txt, args.output_stem, args.width, args.pixel_size,
            height=args.height, pixel_size_y=args.pixel_size_y,
            polarity=args.polarity, scan_direction=args.scan_direction,
            line_scan_direction=args.line_scan_direction,
            scan_pattern=args.scan_pattern, scan_type=args.scan_type,
            limit=args.limit)


if __name__ == "__main__":
    sys.exit(main())
