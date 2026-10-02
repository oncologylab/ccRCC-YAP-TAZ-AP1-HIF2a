from pathlib import Path
import argparse
import sys
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from figures._shared import DATA, read_table, numeric, bh, new_output, save_table, read_bed, overlap_flags


def matched_peaks(frame):
    numeric(frame, ["log2FC"])
    if set(frame.panel) != set(range(1, 7)):
        raise ValueError("All six comparisons are required for the BH family")
    records = []
    for panel, part in frame.groupby("panel", sort=True):
        a = part.loc[part.peak_status == "bound", "log2FC"].to_numpy()
        b = part.loc[part.peak_status == "unbound", "log2FC"].to_numpy()
        if min(len(a), len(b)) == 0:
            raise ValueError("Bound and unbound peaks required")
        test = mannwhitneyu(a, b, alternative="two-sided")
        records.append(dict(panel=panel, n_bound=len(a), n_unbound=len(b),
                            median_bound=np.median(a), median_unbound=np.median(b), raw_p=test.pvalue))
    result = pd.DataFrame(records)
    result["FDR_p"] = bh(result.raw_p)
    return result


def growth(frame, log2_change, baseline_time):
    numeric(frame, ["time", "volume"])
    if frame.duplicated(["group", "mouse", "time"]).any() or (frame.volume <= 0).any():
        raise ValueError("Duplicate observations or nonpositive volumes")
    frame = frame.copy()
    if log2_change:
        baseline = frame.loc[frame.time == baseline_time, ["group", "mouse", "volume"]].rename(columns={"volume": "baseline"})
        frame = frame.merge(baseline, on=["group", "mouse"], how="left", validate="many_to_one")
        if frame.baseline.isna().any():
            raise ValueError("Missing mouse baseline")
        frame["value"] = np.log2(frame.volume / frame.baseline)
    else:
        frame["value"] = frame.volume
    summary = frame.groupby(["group", "time"]).value.agg(n="count", mean="mean", sd="std", sem="sem").reset_index()
    return frame, summary


def growth_workbook(path):
    frame = pd.read_excel(path, header=None)
    groups = frame.iloc[1].ffill()
    records = []
    for column in range(1, len(frame.columns)):
        group = str(groups.iloc[column])
        for row in range(2, len(frame)):
            time, value = frame.iloc[row, 0], frame.iloc[row, column]
            if pd.notna(time) and pd.notna(value):
                records.append(dict(group=group, mouse=f"column_{column}", time=float(time), value=float(value)))
    values = pd.DataFrame(records)
    if values.empty:
        raise ValueError("No growth measurements")
    numeric(values, ["time", "value"])
    summary = values.groupby(["group", "time"]).value.agg(n="count", mean="mean", sd="std", sem="sem").reset_index()
    return values, summary


def main():
    p = argparse.ArgumentParser()
    p.add_argument("analysis", choices=["matched", "growth", "growth-workbook", "overlap"])
    p.add_argument("--input", type=Path, default=DATA / "14_baseline_matched_effects_values.tsv.gz")
    p.add_argument("--other", type=Path)
    p.add_argument("--log2-change", action="store_true")
    p.add_argument("--baseline-time", type=float)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    if a.log2_change and a.baseline_time is None:
        p.error("--log2-change requires --baseline-time")
    if a.analysis == "overlap" and a.other is None:
        p.error("overlap requires --other")
    out = new_output(a.output)
    if a.analysis == "matched":
        frame = read_table(a.input, ["panel", "peak_status", "log2FC"])
        save_table(matched_peaks(frame), out / "tests.tsv")
        save_table(frame.loc[frame.panel.isin([1, 2])], out / "values.tsv")
    elif a.analysis == "growth":
        values, summary = growth(read_table(a.input, ["group", "mouse", "time", "volume"]), a.log2_change, a.baseline_time)
        save_table(values, out / "values.tsv")
        save_table(summary, out / "growth.tsv")
    elif a.analysis == "growth-workbook":
        values, summary = growth_workbook(a.input)
        save_table(values, out / "values.tsv")
        save_table(summary, out / "growth.tsv")
    else:
        first, second = read_bed(a.input), read_bed(a.other)
        first["overlap"] = overlap_flags(first, second)
        second["overlap"] = overlap_flags(second, first)
        save_table(first, out / "first_peaks.tsv")
        save_table(second, out / "second_peaks.tsv")
        save_table(pd.DataFrame([dict(peak_set=label, total=len(part), overlapping=int(part.overlap.sum()),
                                     nonoverlapping=int((~part.overlap).sum()))
                                 for label, part in [("first", first), ("second", second)]]), out / "overlap.tsv")


if __name__ == "__main__":
    main()
