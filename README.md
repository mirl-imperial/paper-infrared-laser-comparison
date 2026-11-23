# Code for: Multimodal characterisation of three infrared lasers for ambient mass spectrometry imaging

William J. C. Francis, Lucy R. Noyes, Milena Micic, Tugba Temel, Ronan A. Battle,
Daniel Simon, Zoltan Takats, Robert T. Murray

This repository contains the custom code behind specific numbers and tables in the manuscript.
The data is in Zenodo at https://doi.org/10.5281/zenodo.23025495.

**Scope.** The repository holds only the calculations whose results cannot be read directly
from the deposited data: the consensus ion list and its sensitivity to the *m/z* tolerance,
the relative ion intensities of the consensus ions and the resolving power behind them, the
signal per crater volume, the crater area statistics, the AFM crater depth and the
chloride-adduct isotope envelope, plus the conversion that produced
the deposited imzML files. It is not a one-command
figure-reproduction pipeline. Figure assembly, layout, labelling and styling code is not
included. Two kinds of file are provided:

- **Calculation functions** (`build_consensus_ion_list.py`, `relative_ion_intensity.py`,
  `signal_per_crater_volume.py`, `crater_area_statistics.py`, `afm_crater_depth.py`,
  `chloride_isotope_envelope.py`).
  NumPy arrays (or a pandas DataFrame of results) in, values out, with no file I/O or
  plotting inside the functions (the readers `read_peak_list` and `read_consensus_inputs` in
  `build_consensus_ion_list.py`, `read_crater_areas` in `crater_area_statistics.py`,
  `read_consensus_mz` in `signal_per_crater_volume.py` and `relative_ion_intensity.py`, and
  `read_sum_spectrum` in `relative_ion_intensity.py` and `chloride_isotope_envelope.py`, and
  `read_afm_height` in `afm_crater_depth.py` excepted). A demo at the
  bottom of each file, under `if __name__ == "__main__":`, runs it on the deposited data.
- **Data-conversion utilities** (`hdi_txt_to_imzml.py`, `imzml_to_intensity_cube.py`).
  Command-line tools whose job is file I/O, reading one format and writing another.
  `signal_per_crater_volume.py` uses `imzml_to_intensity_cube.py` to read the imzML files.

The three lasers are labelled by wavelength and pulse duration, from longest to shortest pulse
duration. The manuscript also gives each a symbol:

| Symbol in the manuscript | Label in the code output | Folder and file names |
|---|---|---|
| Θ (2.94 µm, 3 ns) | 2.94 µm · 3 ns | `2.94um_3ns_Montfort`, `Montfort` |
| Ψ (2.72 µm, 2 ns) | 2.72 µm · 2 ns | `2.72um_2ns_Ivy`, `Ivy` |
| Ω (2.94 µm, 300 ps) | 2.94 µm · 300 ps | `2.94um_300ps_Innolas`, `Innolas` |

Figure numbers below refer to the submitted manuscript.

## Getting the data

Download the Zenodo record (https://doi.org/10.5281/zenodo.23025495) and unzip each section into one
folder, e.g. `laser_comp_data/`, giving `laser_comp_data/01_MSI_HDI_processed_data/`,
`laser_comp_data/02_MSI_IMZML_data/`, `laser_comp_data/03_brain_homogenate_MS/`,
`laser_comp_data/04_crater_micrographs/`, `laser_comp_data/05_consensus_ion_list/`, and so
on. The paths below refer to that layout, and each demo takes the `laser_comp_data` folder as
its argument.

## Installation

```
pip install -r requirements.txt
```

Tested with Python 3.12.1, NumPy 2.4.6, SciPy 1.16.3, pandas 3.0.3, matchms 0.33.1 and
pyimzML 1.5.5 (Windows 10).

## Files

### `build_consensus_ion_list.py` (Supplementary Table S1; Fig. 3; Supplementary Fig. S2)
```
read_peak_list(path) -> (mz, avg_intensity)
read_consensus_inputs(data_dir) -> (replicates, backgrounds, intact)
per_file_noise_floor(avg_intensity) -> floor
global_noise_floor(background_intensities) -> floor
cluster_mz(mz_values, ppm_tolerance=10) -> labels
is_background_ion(mz, exclusion_mz, ppm_tolerance=10) -> bool
ipa_background_ions(backgrounds, noise_floor) -> exclusion_mz
consensus_ion_list(replicates, noise_floor, exclusion_mz, ppm_tolerance=10) -> DataFrame
tolerance_check(replicates, noise_floor, exclusion_mz, ppm_tolerance) -> dict
add_intact_overlay(master, intact, exclusion_mz) -> (DataFrame, stats)
detection_overlap(master, labels) -> {lasers: n_ions}
```
Builds the consensus ion list of brain homogenate from the Waters HDImaging peak lists (500
peaks each) of the six replicates per laser and six acquisitions of neat 2-propanol (IPA)
with no tissue. The noise floor of each IPA acquisition is the mean + 3 s.d. of its peaks
ranked 201–500, and the highest of the six is the global noise floor, reported as 54.58
counts. IPA peaks are clustered by *m/z* at 10 ppm, and a cluster above the floor in at least
four of the six IPA acquisitions is a solvent-background ion (128 ions). All 18 homogenate
peak lists are then pooled and clustered at 10 ppm with matchms (connected components of the
pairwise match graph; the consensus *m/z* is the cluster mean). A laser detects an ion if it
is above the floor in at least four of its six replicates, and an ion is kept if at least one
laser detects it and it is not within 10 ppm of a solvent-background ion. This gives the 280
ions of Supplementary Table S1, of which 130, 122 and 231 are detected by 2.94 µm · 3 ns,
2.72 µm · 2 ns and 2.94 µm · 300 ps. `detection_overlap` counts the ions in each region of
the Fig. 3 diagram: 75 detected by all three lasers, 38 by 2.94 µm · 3 ns and 2.72 µm · 2 ns
only (so the two share 113 ions), 9 by 2.94 µm · 3 ns and 2.94 µm · 300 ps only, 6 by
2.72 µm · 2 ns and 2.94 µm · 300 ps only, and 8, 3 and 141 by 2.94 µm · 3 ns, 2.72 µm · 2 ns
and 2.94 µm · 300 ps alone. Each intact-tissue image's 150-peak list is overlaid as an extra
column by matching its peaks (solvent background removed, no noise floor, since the floor was
measured for the 500-peak homogenate replicates) to the nearest consensus *m/z* within 10 ppm.

`tolerance_check` rebuilds the list with the homogenate peaks clustered at any other
tolerance, keeping the noise floor, the detection rule and the 10 ppm solvent-background list
unchanged, and returns one point of Supplementary Fig. S2: the number of consensus ions and of
ions per laser, the number of kept ions whose cluster joins two or more peaks of one
acquisition, the widest kept cluster, and the spread of each above-floor peak's *m/z* about its
consensus *m/z*. At 10 ppm it gives 280 ions, none merged, with a deviation s.d. of 3.6 ppm and
95% of peaks within ±7.4 ppm; at 3 ppm, 230 ions; at 16 ppm, 285 ions with 2 merged. The
annotated-ion count of Fig. S2 is not produced here (see the annotations below).

`05_consensus_ion_list/master_ion_list_laser_comp.csv` and `ipa_background_ions.csv` in the
data record were produced by this script. The annotations in
`05_consensus_ion_list/id_subset_laser_comp.csv` were made by accurate-mass matching of the
consensus *m/z* values against LIPID MAPS and HMDB (10 ppm, [M−H]⁻ and [M+Cl]⁻ adducts; sum
composition, MSI Level 3) and are not produced by code in this repository.

The demo writes `master_ion_list_laser_comp.csv` and `ipa_background_ions.csv` to the
working directory (about 20 s). `--ppm-tolerance` also prints the `tolerance_check` values at
that tolerance:
```
python build_consensus_ion_list.py laser_comp_data
python build_consensus_ion_list.py laser_comp_data --ppm-tolerance 3
```

### `relative_ion_intensity.py` (Fig. 4; Supplementary Table S1)
```
read_sum_spectrum(path) -> (mz, intensity)
read_consensus_mz(csv_path) -> mz
mean_spectrum(spectra) -> (mz, mean_intensity)
resolving_power(mz, intensity, min_snr=10, noise_window=301) -> (R, n_peaks)
window_intensity(mz, intensity, centres, resolving_power=17800) -> intensities
pqn_log2(X) -> (normalised, log2_normalised)
pairwise_t_tests(log_values, groups, pairs) -> {(a, b): (log2_fold_change, p)}
bh_q_values(p) -> q
is_reported(log2_fc, q, fold=2, q_max=0.05) -> bool
```
Compares the relative intensities of the 280 consensus ions between lasers, using the 18
homogenate summed profile spectra in `03_brain_homogenate_MS/raw_sum_spectra/`. The resolving
power *R* is measured on the mean of the 18 spectra (each scaled to unit total intensity and
interpolated onto the first spectrum's *m/z* axis): over the peaks with a signal-to-noise
ratio of at least 10 (noise = 1.4826 × the 301-point running median of the absolute deviation
from a 301-point running median), *R* is the median of *m/z* / FWHM, which gives 17,844,
used as 17,800. For the ion intensities each spectrum is used on its own *m/z* axis, with no
recalibration or resampling. The intensity of each ion
is the sum of the data points within one FWHM of its consensus *m/z*,
[*m* − *m*/2*R*, *m* + *m*/2*R*] with *R* = 17,800, so every ion has a measured value in every
spectrum. The 280 values of each spectrum are scaled to sum to 1, normalised by probabilistic
quotient normalisation and log2 transformed after adding the 1st percentile of the positive
values. For each ion, each pair of lasers is compared by Student's t-test on the difference
of the mean log2 intensities, with the residual variance pooled over all three lasers (15
degrees of freedom); the difference is the log2 fold change. Benjamini–Hochberg q-values are
computed per pair over the 280 ions, and a difference is reported when it is at least 2-fold
with q < 0.05.

This gives 101 reported differences between 2.94 µm · 300 ps and 2.94 µm · 3 ns, 80 between
2.94 µm · 300 ps and 2.72 µm · 2 ns, and none between the two nanosecond lasers, whose largest
difference is 1.85-fold. The demo also prints the fold changes of the annotated PA and cLPA
species and of fumarate, FA 16:1, *N*-acetylaspartate, ascorbate and malate (names from
`05_consensus_ion_list/id_subset_laser_comp.csv`), and writes `relative_ion_intensity.csv`
(per ion: mean normalised intensity per laser and, per pair, log2 fold change, p and q) to
the working directory:
```
python relative_ion_intensity.py laser_comp_data
```

### `signal_per_crater_volume.py` (Table 1; Results, Mouse brain imaging)
```
read_consensus_mz(csv_path) -> mz
within_ppm(target, pool, ppm=10) -> bool
common_ion_set(mz_axes, reference_mz, ppm=10) -> mz
mean_signal_per_pixel(intensities, mz, ion_set) -> signal
fold_range(values) -> max / min
crater_volume(area_mean_um2, area_sd_um2, depth_um=2) -> (volume, sd)
signal_per_volume(signal_per_pixel, area_mean_um2, area_sd_um2, depth_um=2) -> (value, sd)
sample_cv_percent(values) -> cv
```
Compares the signal of the three mouse brain images of Fig. 2 against the tissue volume each
laser ablates. The common ion set is the consensus ions with a peak within 10 ppm of the
consensus *m/z* in all three images: 35 ions, the same ions as the rows of Table S1 with all
three Intact values. For each image the intensities of those
ions are summed over every pixel of the image (no tissue mask or region selection) and
divided by the number of pixels, giving the mean signal per pixel, reported as 7.35e3, 3.40e3
and 1.33e4 counts for 2.94 µm · 3 ns, 2.72 µm · 2 ns and 2.94 µm · 300 ps, a 3.9-fold range.
The crater areas are the brain SEM measurements of `crater_area_statistics.py` (159 ± 16,
87 ± 7 and 309 ± 60 µm², a 3.5-fold range). The crater volume is the mean area times the
crater depth (2 µm by default, `--depth`, from `afm_crater_depth.py`), with the area s.d. times the depth as its ±:
317 ± 33, 175 ± 15 and 617 ± 121 µm³. The signal per volume is the mean signal per pixel
divided by the crater volume, with ± equal to the value times the relative s.d. of the area:
23.2 ± 2.4, 19.5 ± 1.6 and 21.6 ± 4.2 counts µm⁻³. Their scatter, the sample coefficient of
variation (ddof = 1) of the three values, is 8.7%. The image carries no error term of its
own, because there is one image per laser.

The demo reads the three imzML files (about 20 s) and prints the values at full precision
and rounded as in the manuscript:
```
python signal_per_crater_volume.py laser_comp_data
```

### `crater_area_statistics.py` (Table 1; Fig. 2d–f; Supplementary Fig. S6)
```
read_crater_areas(csv_path) -> areas_um2
brain_measurement_csv(data_dir, stem) -> path
summarise_areas(areas_um2) -> {"n", "mean_area_um2", "sd_area_um2", "equiv_diameter_um"}
```
Summarises the brain crater areas measured by SEM on the imaged sections, from the per-crater
ImageJ ellipse-fit measurements in `04_crater_micrographs/brain_SEM_measurements/` ("Area"
column, µm²). For each laser it returns n, the mean, the sample s.d. (ddof = 1) and the
diameter of a circle of the mean area. The areas are reported as 159 ± 16, 87 ± 7 and
309 ± 60 µm² (mean ± s.d.; n = 12, 12 and 11) for 2.94 µm · 3 ns, 2.72 µm · 2 ns and
2.94 µm · 300 ps.
```
python crater_area_statistics.py laser_comp_data
```

### `afm_crater_depth.py` (Supplementary Fig. S3; Table 1)
```
read_afm_height(path) -> (height_m, pixel_size_m)
level_to_undisturbed_surface(height, n_iter=6, mad_threshold=2.5) -> (levelled, surface_mask)
crater_depth_summary(levelled_um, floor_threshold_um=-1.5) -> {"deepest_um", "floor_mean_um", "rim_um", "peak_to_valley_um"}
middle_row_profile(levelled_um) -> heights
```
Measures the crater depth used for the crater volumes of Table 1, from the AFM height map of a
single 2.72 µm · 2 ns crater in intact mouse brain
(`06_AFM_topography/2.72um_2ns_Ivy_brain_25um_scan.ibw`, HeightRetrace layer). A plane is
fitted to the undisturbed tissue surface, iteratively excluding pixels more than 2.5 robust
s.d. (1.4826 × MAD) from the median residual, and heights are referenced to the median of that
surface. The crater floor (pixels deeper than −1.5 µm) averages −1.91 µm, the deepest point is
−2.22 µm, the highest rim point is +2.03 µm and the peak-to-valley height is 4.25 µm, reported
as a ~2 µm deep crater with a ~2 µm rim (peak-to-valley ~4.3 µm). The line profile of
Supplementary Fig. S3c (`middle_row_profile`) reaches −2.18 µm. The demo also levels the 90 µm
scan of the crater array:
```
python afm_crater_depth.py laser_comp_data
```

### `chloride_isotope_envelope.py` (Supplementary Note S1; Supplementary Fig. S5)
```
read_sum_spectrum(path) -> (mz, intensity)
theoretical_chloride_envelope(neutral=NEUTRAL, n_peaks=4) -> (rel, centroid_mz, cl37_pct)
non_chlorinated_isobar_m2_ratio(isobar=ISOBAR) -> ratio
integrate_peak(mz, intensity, centre, halfwidth=0.05) -> (area, centroid_mz)
measure_envelope(mz, intensity, centres) -> (ratios, observed_mz)
summarise_replicates(ratios, mz_obs_m0, mz_theo_m0) -> {"ratio_mean", "ratio_sd", ...}
```
Tests the assignment of *m/z* 862.65, the green channel of Fig. 2, as
[HexCer 42:1;O3 + ³⁵Cl]⁻. The theoretical envelope of C48H93NO9·³⁵Cl⁻ is binned to nominal
peaks M+0 to M+3 (*m/z* 862.6544, 863.6578, 864.6546 and 865.6564; 100, 53.7, 48.0 and
20.6% of M+0, with ³⁷Cl carrying 67% of M+2 and 83% of M+3). The control is the M+2/M+0
ratio of a non-chlorinated composition of the same nominal mass, C51H92NO10 (0.18). The
measured envelope comes from the 18 homogenate profile spectra in
`03_brain_homogenate_MS/raw_sum_spectra/`: each peak is integrated over ±0.05 Da around its
theoretical *m/z* after subtracting a baseline (median of the two outermost points at each
end), and its *m/z* is the intensity-weighted mean of the points at or above half maximum.
Over the 18 replicates M+1, M+2 and M+3 are 53.1 ± 1.7, 48.0 ± 1.3 and 19.4 ± 0.3% of
M+0, within 1.3 percentage points of the prediction at every peak, and the measured M+2/M+0 of
0.48 lies 22 times its between-replicate s.d. from 0.18. The M+0 *m/z* is observed at
862.6521, −2.7 ± 2.3 ppm from theory.
```
python chloride_isotope_envelope.py laser_comp_data
```

### `imzml_to_intensity_cube.py` (`02_MSI_IMZML_data/`)
```
imzml_to_intensity_cube(imzml_path) -> (intensities, mz, coords)
python imzml_to_intensity_cube.py IMZML --out OUTPUT.npz
```
Reads a continuous-mode centroid imzML + `.ibd` pair into a pixel-by-peak intensity array
(float32), the *m/z* axis shared by every pixel (float64, in file order) and the (x, y, z)
position of each pixel, and from the command line saves them to a compressed `.npz`.
```
python imzml_to_intensity_cube.py laser_comp_data/02_MSI_IMZML_data/2.94um_3ns_Montfort/2026_04_29_MB1_M2_33_Montfort_25um_5_centroid.imzML --out 2.94um_3ns_Montfort.npz
```

### `hdi_txt_to_imzml.py` (`02_MSI_IMZML_data/`)
```
python hdi_txt_to_imzml.py INPUT.txt OUTPUT_STEM --width W --pixel-size P --polarity {positive,negative}
```
Converts a Waters HDImaging (Maldichrom) peak-picked pixel `.txt` export into a centroid,
continuous-mode imzML + `.ibd` pair. The imzML files in `02_MSI_IMZML_data/` were made this
way from the exports in `01_MSI_HDI_processed_data/`. HDImaging's own imzML export re-reads
the raw file and writes full, uncentroided profile spectra, so this script works from the
much smaller peak-picked export instead.

**The script performs no centroiding, peak picking or binning of its own, and applies no
smoothing, filtering, resampling or intensity transform.** Peak picking was done upstream in
HDImaging (Experimental: MS resolution 20 000, 0.02 Da *m/z* window, 150 most intense
peaks per image), and the export already holds one shared peak list with each pixel's
intensity at each of those *m/z* values. The script writes that peak list and the per-pixel
intensities to imzML/ibd unchanged, in the export's order: HDImaging lists the peaks by
intensity rank, so the *m/z* array is not in increasing order and the files do not declare
an increasing-*m/z* scan. Its only computation is each pixel's (x, y) raster position,
reconstructed from the acquisition order. After writing, it adds the pixel size to the imzML
header and sets the `<run>` id to the output file's base name.

`--width` is always required. It is cross-checked against the raster width detected from the
stage x column, which resets at the start of each scan line, and a mismatch prints a warning.
Optional flags cover other raster geometries (scan direction, line-scan direction, scan
pattern, scan type). The module docstring describes the `.txt` format, including the two
trailing columns that look like pixel coordinates but are MassLynx bookkeeping. The deposited
files used:

| Dataset | `--width` | `--pixel-size` |
|---|---|---|
| 2.94 µm · 3 ns (`2.94um_3ns_Montfort`) | 296 | 25 |
| 2.72 µm · 2 ns (`2.72um_2ns_Ivy`) | 294 | 25 |
| 2.94 µm · 300 ps (`2.94um_300ps_Innolas`) | 304 | 25 |

all with `--polarity negative`. For example:
```
python hdi_txt_to_imzml.py laser_comp_data/01_MSI_HDI_processed_data/2.94um_3ns_Montfort/2026_04_29_MB1_M2_33_Montfort_25um_5.txt 2026_04_29_MB1_M2_33_Montfort_25um_5_centroid --width 296 --pixel-size 25 --polarity negative
```

## Dependencies

NumPy throughout. `build_consensus_ion_list.py` also uses pandas, matchms and SciPy (for
the *m/z* clustering). `relative_ion_intensity.py` uses pandas and SciPy (peak finding, t distribution and
Benjamini–Hochberg q-values). `afm_crater_depth.py` uses igor2 to read the Igor binary
waves. `signal_per_crater_volume.py` and `crater_area_statistics.py`
use pandas to read CSV files. `hdi_txt_to_imzml.py` and
`imzml_to_intensity_cube.py` (and so `signal_per_crater_volume.py`) require pyimzML. See
`requirements.txt`.

## Licence

MIT (see `LICENSE`).

## Citation

Please cite the associated article (see `CITATION.cff`). Software DOI (all versions): https://doi.org/10.5281/zenodo.23028050.
