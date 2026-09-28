#!/usr/bin/env python3
"""Consensus ion list for brain homogenate (Supplementary Table S1; Fig. 4).

Builds the consensus ion list from the six homogenate replicates per laser and
the six neat-IPA background acquisitions, then overlays the single 150-peak
peak list of each intact-tissue image as an extra column per laser.

    python build_consensus_ion_list.py laser_comp_data

Method
------
1. Noise floor. For each background acquisition, the peaks are ranked by
   Avg Intensity and the 300 lowest-ranked of the 500 (ranks 201-500) are
   taken as the noise population. The floor for that acquisition is their
   mean + 3 s.d. (population s.d.). The highest of the six floors is used as a
   single global threshold.
2. IPA solvent background. All background peaks are pooled and clustered by
   m/z at 10 ppm (step 3). A cluster found above the noise floor in at least
   4 of the 6 background acquisitions is a solvent-background ion.
3. Clustering. All homogenate peaks (18 replicates) are pooled. Two m/z values
   match if |a - b| / mean(a, b) <= 10 ppm (matchms PrecursorMzMatch), and a
   cluster is a connected component of that match graph
   (scipy.sparse.csgraph.connected_components). The consensus m/z of a cluster
   is the mean of its member m/z values.
4. Detection. A laser detects a cluster if the cluster is above the noise
   floor in at least 4 of that laser's 6 replicates. The reported homogenate
   intensity is the mean Avg Intensity over the replicates above the floor.
5. A cluster is kept if at least one laser detects it and its consensus m/z
   is not within 10 ppm of a solvent-background ion.
6. Intact-tissue overlay. Peaks of each intact 150-peak list that are above
   the noise floor and not solvent background are matched to the nearest
   consensus m/z within 10 ppm. Unmatched intact peaks are counted only.

The calculation functions take NumPy arrays (or a DataFrame of the result) and
return values. File reading is confined to ``read_peak_list`` and
``read_consensus_inputs``.

Inputs (deposited data)
-----------------------
    03_brain_homogenate_MS/<laser folder>/<Montfort|Ivy|Innolas>_1..6.csv
    03_brain_homogenate_MS/background/Background_1..6.csv
    01_MSI_HDI_processed_data/intact_peaklists/<Montfort|Ivy|Innolas>_150.csv

Each is a Waters HDImaging peak-list export with columns "M/z", "Max
Intensity", "Analyte", "Min Intensity", "Avg Intensity", "Sum Intensity".

Demo outputs (written to --out-dir, default the working directory)
------------------------------------------------------------------
    master_ion_list_laser_comp.csv   one row per consensus ion: m/z, then per
                                     laser "<laser> (Homogenate)", "<laser>
                                     (Homogenate n)" and "<laser> (Intact)"
    ipa_background_ions.csv          solvent-background m/z values, with the
                                     global noise floor in a leading comment
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PPM_TOLERANCE = 10.0
NOISE_TAIL_START = 200      # rank 201 (0-based index 200) ...
NOISE_TAIL_END = 500        # ... to rank 500: the lowest 300 of 500 peaks
MIN_REPS_DETECTED = 4
N_REPLICATES = 6

# Laser label -> (data folder, file stem). Longest to shortest pulse duration.
LASERS = {
    "2.94 µm · 3 ns":   ("2.94um_3ns_Montfort",   "Montfort"),
    "2.72 µm · 2 ns":   ("2.72um_2ns_Ivy",        "Ivy"),
    "2.94 µm · 300 ps": ("2.94um_300ps_Innolas",  "Innolas"),
}


# ================================================================ readers ===

def read_peak_list(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read one HDImaging peak-list CSV. Returns (m/z, Avg Intensity)."""
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    return df["M/z"].to_numpy(float), df["Avg Intensity"].to_numpy(float)


def read_consensus_inputs(data_dir: Path) -> tuple[dict, list, dict]:
    """Read every input of the consensus list from the deposited data.

    Returns
    -------
    replicates : {laser label: [(mz, avg_intensity), ...] for replicates 1..6}
    backgrounds : [(mz, avg_intensity), ...] for Background_1..6
    intact : {laser label: (mz, avg_intensity)} from the 150-peak lists
    """
    homog = Path(data_dir) / "03_brain_homogenate_MS"
    intact_dir = Path(data_dir) / "01_MSI_HDI_processed_data" / "intact_peaklists"
    replicates = {
        label: [read_peak_list(homog / folder / f"{stem}_{i}.csv")
                for i in range(1, N_REPLICATES + 1)]
        for label, (folder, stem) in LASERS.items()
    }
    backgrounds = [read_peak_list(homog / "background" / f"Background_{i}.csv")
                   for i in range(1, N_REPLICATES + 1)]
    intact = {label: read_peak_list(intact_dir / f"{stem}_150.csv")
              for label, (_, stem) in LASERS.items()}
    return replicates, backgrounds, intact


# ==================================================== calculation functions ===

def per_file_noise_floor(avg_intensity: np.ndarray) -> float:
    """Noise floor of one acquisition: mean + 3 s.d. of peaks ranked 201-500.

    ``avg_intensity`` is the Avg Intensity column of one 500-peak list. The
    s.d. is the population value (ddof = 0).
    """
    avg_intensity = np.asarray(avg_intensity, dtype=float)
    if avg_intensity.size < NOISE_TAIL_END:
        raise ValueError(f"noise floor needs >= {NOISE_TAIL_END} peaks, "
                         f"got {avg_intensity.size}")
    tail = np.sort(avg_intensity)[::-1][NOISE_TAIL_START:NOISE_TAIL_END]
    return float(tail.mean() + 3 * tail.std(ddof=0))


def global_noise_floor(background_intensities: list[np.ndarray]) -> float:
    """Highest of the per-acquisition noise floors of the background scans."""
    return max(per_file_noise_floor(a) for a in background_intensities)


def cluster_mz(mz_values: np.ndarray, ppm_tolerance: float = PPM_TOLERANCE) -> np.ndarray:
    """Cluster m/z values at a ppm tolerance; returns one integer label per value.

    Two values match if |a - b| / mean(a, b) <= ppm_tolerance * 1e-6. Clusters
    are the connected components of the pairwise match graph.
    """
    from matchms import Spectrum
    from matchms.similarity import PrecursorMzMatch
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components

    spectra = [Spectrum(mz=np.array([]), intensities=np.array([]),
                        metadata={"precursor_mz": float(mz)})
               for mz in mz_values]
    similarity = PrecursorMzMatch(tolerance=ppm_tolerance, tolerance_type="ppm")
    match_matrix = similarity.matrix(spectra, spectra, is_symmetric=True)
    _, labels = connected_components(csr_matrix(match_matrix), directed=False)
    return labels


def is_background_ion(mz: float, exclusion_mz: np.ndarray,
                      ppm_tolerance: float = PPM_TOLERANCE) -> bool:
    """True if ``mz`` is within ``ppm_tolerance`` of any value in ``exclusion_mz``."""
    exclusion_mz = np.asarray(exclusion_mz, dtype=float)
    if exclusion_mz.size == 0:
        return False
    mean_mz = (mz + exclusion_mz) / 2
    return bool(np.any(np.abs(mz - exclusion_mz) / mean_mz <= ppm_tolerance * 1e-6))


def _pool(peak_lists: list[tuple[np.ndarray, np.ndarray]], group: str = "") -> pd.DataFrame:
    """Stack (mz, intensity) pairs into one table with a replicate index."""
    frames = [pd.DataFrame({"mz": mz, "avg_intensity": inten, "group": group,
                            "replicate": i})
              for i, (mz, inten) in enumerate(peak_lists)]
    return pd.concat(frames, ignore_index=True)


def _replicate_maxima(cluster_rows: pd.DataFrame, n_replicates: int) -> pd.Series:
    """Highest Avg Intensity of a cluster in each replicate (0 where absent)."""
    return (cluster_rows.groupby("replicate")["avg_intensity"].max()
            .reindex(range(n_replicates), fill_value=0.0))


def ipa_background_ions(backgrounds: list[tuple[np.ndarray, np.ndarray]],
                        noise_floor: float) -> np.ndarray:
    """Consensus m/z of clusters above the noise floor in >= 4 of 6 backgrounds."""
    pooled = _pool(backgrounds)
    pooled["cluster"] = cluster_mz(pooled["mz"].to_numpy())
    consensus = pooled.groupby("cluster")["mz"].mean()
    out = []
    for cluster_id, rows in pooled.groupby("cluster"):
        per_rep = _replicate_maxima(rows, len(backgrounds))
        if (per_rep >= noise_floor).sum() >= MIN_REPS_DETECTED:
            out.append(consensus.loc[cluster_id])
    return np.array(out)


def consensus_ion_list(replicates: dict[str, list[tuple[np.ndarray, np.ndarray]]],
                       noise_floor: float, exclusion_mz: np.ndarray) -> pd.DataFrame:
    """Steps 3-5: the consensus list from the homogenate replicates.

    Parameters
    ----------
    replicates : {laser label: [(mz, avg_intensity), ...]}, six replicates each
    noise_floor : global noise floor
    exclusion_mz : solvent-background m/z values

    Returns
    -------
    DataFrame sorted by m/z with "m/z" and, per laser, "<label> (Homogenate)"
    (mean intensity over replicates above the floor, NaN if not detected) and
    "<label> (Homogenate n)" ("k/6", replicates above the floor).
    """
    pooled = pd.concat([_pool(reps, label) for label, reps in replicates.items()],
                       ignore_index=True)
    pooled["cluster"] = cluster_mz(pooled["mz"].to_numpy())
    consensus = pooled.groupby("cluster")["mz"].mean()

    rows = []
    for cluster_id, group in pooled.groupby("cluster"):
        mz = consensus.loc[cluster_id]
        if is_background_ion(mz, exclusion_mz):
            continue
        row = {"m/z": mz}
        any_detected = False
        for label, reps in replicates.items():
            per_rep = _replicate_maxima(group[group["group"] == label], len(reps))
            above = per_rep >= noise_floor
            n_above = int(above.sum())
            row[f"{label} (Homogenate n)"] = f"{n_above}/{len(reps)}"
            if n_above >= MIN_REPS_DETECTED:
                row[f"{label} (Homogenate)"] = per_rep[above].mean()
                any_detected = True
            else:
                row[f"{label} (Homogenate)"] = np.nan
        if any_detected:
            rows.append(row)

    columns = ["m/z"] + [c for label in replicates
                         for c in (f"{label} (Homogenate)", f"{label} (Homogenate n)")]
    return pd.DataFrame(rows, columns=columns).sort_values("m/z").reset_index(drop=True)


def add_intact_overlay(master: pd.DataFrame,
                       intact: dict[str, tuple[np.ndarray, np.ndarray]],
                       noise_floor: float, exclusion_mz: np.ndarray
                       ) -> tuple[pd.DataFrame, dict]:
    """Step 6: add an "<label> (Intact)" column per laser from its 150-peak list.

    Returns the extended table and, per laser, the number of intact peaks above
    the floor (solvent background removed), matched to the list, and unmatched.
    """
    master = master.copy()
    ref = master["m/z"].to_numpy()
    stats = {}
    for label, (mz, inten) in intact.items():
        keep = [(m, a) for m, a in zip(mz, inten)
                if a >= noise_floor and not is_background_ion(m, exclusion_mz)]
        col = np.full(len(master), np.nan)
        n_matched = 0
        for m, a in keep:
            d = np.abs(ref - m) / ((ref + m) / 2)
            idx = int(np.argmin(d))
            if d[idx] <= PPM_TOLERANCE * 1e-6:
                if np.isnan(col[idx]) or a > col[idx]:
                    col[idx] = a
                n_matched += 1
        master[f"{label} (Intact)"] = col
        stats[label] = dict(n_above_floor=len(keep), n_matched=n_matched,
                            n_unmatched=len(keep) - n_matched)
    return master, stats


def detection_overlap(master: pd.DataFrame, labels: list[str]) -> dict:
    """Ion counts for each region of the three-laser overlap (Fig. 4).

    Keys are tuples of the labels detecting the ions in that region only,
    e.g. (label_a,) is "detected by label_a alone".
    """
    detected = np.column_stack([master[f"{lab} (Homogenate)"].notna().to_numpy()
                                for lab in labels])
    counts = {}
    for pattern in sorted({tuple(r) for r in detected}, key=lambda p: (-sum(p), p)):
        members = tuple(lab for lab, on in zip(labels, pattern) if on)
        counts[members] = int((detected == np.array(pattern)).all(axis=1).sum())
    return counts


# =================================================================== demo ===

def write_ipa_background(path: Path, exclusion_mz: np.ndarray, noise_floor: float) -> None:
    """Write the solvent-background m/z list with the noise floor as a comment."""
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write(f"# noise_floor: {float(noise_floor)!r}\n")
        fh.write("mz\n")
        for mz in np.sort(np.asarray(exclusion_mz, dtype=float)):
            fh.write(f"{float(mz)!r}\n")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("data_dir", type=Path, help="unzipped data record, e.g. laser_comp_data")
    p.add_argument("--out-dir", type=Path, default=Path("."),
                   help="folder for the two output CSVs (default: working directory)")
    args = p.parse_args()

    replicates, backgrounds, intact = read_consensus_inputs(args.data_dir)
    floor = global_noise_floor([a for _, a in backgrounds])
    exclusion = ipa_background_ions(backgrounds, floor)
    master = consensus_ion_list(replicates, floor, exclusion)
    master, intact_stats = add_intact_overlay(master, intact, floor, exclusion)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    master.to_csv(args.out_dir / "master_ion_list_laser_comp.csv", index=False)
    write_ipa_background(args.out_dir / "ipa_background_ions.csv", exclusion, floor)

    labels = list(LASERS)
    print(f"Global noise floor: {floor:.2f} counts")
    print(f"IPA solvent-background ions excluded: {len(exclusion)}")
    print(f"Consensus ions: {len(master)}")
    for label in labels:
        print(f"  {label:<18} detected in homogenate: "
              f"{int(master[f'{label} (Homogenate)'].notna().sum())}")
    print("Overlap of homogenate detections (ions in each region only):")
    for members, n in detection_overlap(master, labels).items():
        print(f"  {' + '.join(members):<52} {n}")
    print("Intact-tissue overlay (one 150-peak list per laser):")
    for label, s in intact_stats.items():
        print(f"  {label:<18} {s['n_above_floor']} peaks above floor, "
              f"{s['n_matched']} matched, {s['n_unmatched']} unmatched")
    print(f"\nWrote {args.out_dir / 'master_ion_list_laser_comp.csv'}")
    print(f"Wrote {args.out_dir / 'ipa_background_ions.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
