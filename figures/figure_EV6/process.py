from pathlib import Path
import argparse
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from figures._shared import read_table, numeric, new_output, save_table


def qpcr(frame, mode):
    frame = frame.copy()
    keys = ["panel", "antibody", "condition", "replicate"]
    if frame.duplicated(keys).any():
        raise ValueError("One measurement per biological replicate required")
    if mode == "ct":
        numeric(frame, ["ct_input", "ct_ip", "input_dilution"])
        if (frame.input_dilution < 1).any():
            raise ValueError("Input dilution must be at least one")
        frame["adjusted_input_ct"] = frame.ct_input - np.log2(frame.input_dilution)
        frame["percent_input"] = 100 * 2 ** (frame.adjusted_input_ct - frame.ct_ip)
    else:
        numeric(frame, ["percent_input"])
    if (frame.percent_input < 0).any():
        raise ValueError("Negative percent input")
    summary = frame.groupby(["panel", "antibody", "condition"]).percent_input.agg(n="count", mean="mean", sd="std", sem="sem").reset_index()
    return frame, summary


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--mode", choices=["ct", "percent-input"], required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    required = ["panel", "antibody", "condition", "replicate"]
    required += ["ct_input", "ct_ip", "input_dilution"] if a.mode == "ct" else ["percent_input"]
    values, summary = qpcr(read_table(a.input, required), a.mode)
    out = new_output(a.output)
    save_table(values, out / "qpcr.tsv")
    save_table(summary, out / "summary.tsv")
