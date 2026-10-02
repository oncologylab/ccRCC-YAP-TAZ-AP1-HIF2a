from pathlib import Path
import argparse
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from figures._shared import read_table, numeric, new_output, save_table, read_bed, overlap_flags


def occupancy(regions, peaks):
    frame = regions.copy()
    for name, sites in peaks.items():
        frame[name] = overlap_flags(frame, sites)
    names = list(peaks)
    frame["class"] = ["+".join(name for name, flag in zip(names, flags) if flag) or "None"
                      for flags in frame[names].itertuples(index=False, name=None)]
    counts = frame.groupby("class").size().rename("n_regions").reset_index()
    counts["percent"] = counts.n_regions / len(frame) * 100
    return frame, counts


def sort_signals(frame, signal_columns, control_columns):
    numeric(frame, signal_columns)
    if frame.region.duplicated().any():
        raise ValueError("Duplicate regions")
    frame = frame.copy()
    frame["control_mean"] = frame[control_columns].mean(axis=1)
    return frame.sort_values(["occupancy", "change", "control_mean"], ascending=[True, True, False], kind="stable")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("analysis", choices=["occupancy", "heatmap", "profiles"])
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--yap", type=Path)
    p.add_argument("--hif2a", type=Path)
    p.add_argument("--jun", type=Path)
    p.add_argument("--signals", nargs="+")
    p.add_argument("--control", nargs="+")
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    if a.analysis == "occupancy" and not all([a.yap, a.hif2a, a.jun]):
        p.error("occupancy requires --yap, --hif2a and --jun")
    if a.analysis == "heatmap" and (not a.signals or not a.control or not set(a.control).issubset(a.signals)):
        p.error("heatmap requires --signals and a subset --control")
    out = new_output(a.output)
    if a.analysis == "occupancy":
        regions, counts = occupancy(read_bed(a.input), {name: read_bed(path) for name, path in [("YAP", a.yap), ("HIF2a", a.hif2a), ("JUN", a.jun)]})
        save_table(regions, out / "regions.tsv")
        save_table(counts, out / "counts.tsv")
    elif a.analysis == "heatmap":
        frame = read_table(a.input, ["region", "occupancy", "change", *a.signals])
        sorted_frame = sort_signals(frame, a.signals, a.control)
        save_table(sorted_frame[["region", "occupancy", "change", "control_mean"]], out / "regions.tsv")
        sorted_frame.set_index("region")[a.signals].to_csv(out / "heatmap.tsv", sep="\t")
    else:
        frame = read_table(a.input, ["panel", "group", "region", "x", "value"])
        numeric(frame, ["x", "value"])
        if frame.duplicated(["panel", "group", "region", "x"]).any():
            raise ValueError("Duplicate region/bin")
        profile = frame.groupby(["panel", "group", "x"]).value.agg(value="mean", sd="std", n="count").reset_index()
        save_table(profile, out / "profiles.tsv")


if __name__ == "__main__":
    main()
