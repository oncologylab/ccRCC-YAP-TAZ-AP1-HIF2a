from pathlib import Path
import argparse
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from figures._shared import read_table, numeric, new_output, save_table


def pulldown(frame, reference):
    keys = ["panel", "target", "condition", "replicate"]
    if frame.duplicated(keys).any():
        raise ValueError("Duplicate band measurement")
    numeric(frame, ["ip", "ip_background", "input", "input_background"])
    frame = frame.copy()
    frame["corrected_ip"] = frame.ip - frame.ip_background
    frame["corrected_input"] = frame["input"] - frame.input_background
    if (frame.corrected_input <= 0).any() or (frame.corrected_ip < 0).any():
        raise ValueError("Invalid background-subtracted band intensities")
    frame["ip_input_ratio"] = frame.corrected_ip / frame.corrected_input
    baseline = frame.loc[frame.condition == reference, ["panel", "target", "replicate", "ip_input_ratio"]]
    baseline = baseline.rename(columns={"ip_input_ratio": "reference_ratio"})
    frame = frame.merge(baseline, on=["panel", "target", "replicate"], how="left", validate="many_to_one")
    if frame.reference_ratio.isna().any() or (frame.reference_ratio <= 0).any():
        raise ValueError("Positive matched reference bands required")
    frame["relative_binding"] = frame.ip_input_ratio / frame.reference_ratio
    summary = frame.groupby(["panel", "target", "condition"]).relative_binding.agg(n="count", mean="mean", sd="std", sem="sem").reset_index()
    return frame, summary


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--reference", required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    frame = read_table(a.input, ["panel", "target", "condition", "replicate", "ip", "ip_background", "input", "input_background"])
    values, summary = pulldown(frame, a.reference)
    out = new_output(a.output)
    save_table(values, out / "binding.tsv")
    save_table(summary, out / "summary.tsv")
