from pathlib import Path
import argparse
import sys
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from figures._shared import read_table, numeric, new_output, save_table, read_bed, overlap_flags


def composition(peaks, motifs):
    peaks = peaks.reset_index(drop=True).copy()
    for family, sites in motifs.items():
        peaks[family] = overlap_flags(peaks, sites)
    names = ["AP1", "TEAD", "HIF"]
    peaks["motif_pattern"] = ["/".join(name for name, flag in zip(names, flags) if flag) or "none"
                              for flags in peaks[names].itertuples(index=False, name=None)]
    counts = peaks.groupby(["stratum", "motif_pattern"]).size().rename("n").reset_index()
    counts["percent"] = counts.n / counts.groupby("stratum").n.transform("sum") * 100
    return peaks, counts


def nearest_distances(query, motifs):
    records = []
    for family, target in motifs.items():
        targets = {chrom: np.sort((part.start.to_numpy(float) + part.end.to_numpy(float)) / 2)
                   for chrom, part in target.groupby("chrom")}
        for region, row in enumerate(query.itertuples(index=False)):
            candidates = targets.get(row.chrom)
            distance = np.nan
            if candidates is not None and len(candidates):
                midpoint = (row.start + row.end) / 2
                index = np.searchsorted(candidates, midpoint)
                near = candidates[max(0, index - 1):min(len(candidates), index + 1)]
                distance = np.min(np.abs(near - midpoint))
            records.append(dict(region=region, chrom=row.chrom, start=row.start, end=row.end,
                                family=family, distance_bp=distance,
                                log10_distance_plus1=np.log10(distance + 1)))
    return pd.DataFrame(records)


def imaging(frame):
    numeric(frame, ["value"])
    if frame.duplicated(["panel", "group", "nucleus"]).any():
        raise ValueError("Duplicate nucleus")
    records = []
    for panel, part in frame.groupby("panel"):
        groups = part.group.unique()
        if len(groups) != 2 or "Ctrl" not in groups:
            raise ValueError("Ctrl and one treated group required")
        treated = next(x for x in groups if x != "Ctrl")
        a, b = [part.loc[part.group == group, "value"] for group in ["Ctrl", treated]]
        test = mannwhitneyu(a, b, alternative="two-sided")
        records.append(dict(panel=panel, treated=treated, n_ctrl=len(a), n_treated=len(b), pvalue=test.pvalue))
    return pd.DataFrame(records)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("analysis", choices=["composition", "distance", "volcano", "imaging"])
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--ap1", type=Path)
    p.add_argument("--tead", type=Path)
    p.add_argument("--hif", type=Path)
    p.add_argument("--alpha", type=float)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    if a.analysis in ["composition", "distance"] and not all([a.tead, a.hif]):
        p.error("--tead and --hif required")
    if a.analysis == "composition" and not a.ap1:
        p.error("--ap1 required")
    if a.analysis == "volcano" and a.alpha is None:
        p.error("--alpha required")
    out = new_output(a.output)
    if a.analysis == "composition":
        peaks = read_table(a.input, ["chrom", "start", "end", "stratum"])
        values, counts = composition(peaks, {name: read_bed(path) for name, path in [("AP1", a.ap1), ("TEAD", a.tead), ("HIF", a.hif)]})
        save_table(values, out / "peaks.tsv")
        save_table(counts, out / "composition.tsv")
    elif a.analysis == "distance":
        result = nearest_distances(read_bed(a.input), {"TEAD": read_bed(a.tead), "HIF": read_bed(a.hif)})
        save_table(result, out / "distances.tsv")
    elif a.analysis == "volcano":
        frame = read_table(a.input, ["motif", "differential_score", "padj"])
        numeric(frame, ["differential_score", "padj"])
        if not frame.padj.between(0, 1).all():
            raise ValueError("Adjusted P values outside [0, 1]")
        frame["minus_log10_padj"] = -np.log10(frame.padj.clip(lower=np.nextafter(0., 1.)))
        frame["status"] = np.select([(frame.padj < a.alpha) & (frame.differential_score < 0),
                                     (frame.padj < a.alpha) & (frame.differential_score > 0)], ["Ctrl", "KD"], "NS")
        save_table(frame, out / "volcano.tsv")
    else:
        frame = read_table(a.input, ["panel", "group", "nucleus", "value"])
        save_table(imaging(frame), out / "tests.tsv")
        save_table(frame, out / "values.tsv")


if __name__ == "__main__":
    main()
