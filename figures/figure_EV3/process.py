from pathlib import Path
import argparse
import sys
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from figures._shared import read_table, numeric, bh, new_output, save_table


def paired_signal(frame):
    numeric(frame, ["control", "treated"])
    if (frame[["control", "treated"]] < 0).any().any():
        raise ValueError("Negative chromatin signal")
    if frame.duplicated(["comparison", "region"]).any():
        raise ValueError("Duplicate paired region")
    values = frame.copy()
    values["log2FC"] = np.log2(values.treated + 1) - np.log2(values.control + 1)
    records = []
    for comparison, part in frame.groupby("comparison", sort=False):
        test = stats.wilcoxon(part.control, part.treated, alternative="two-sided")
        records.append(dict(comparison=comparison, n_regions=len(part), mean_control=part.control.mean(),
                            mean_treated=part.treated.mean(), sd_control=part.control.std(), sd_treated=part.treated.std(), raw_p=test.pvalue))
    result = pd.DataFrame(records)
    result["FDR_p"] = bh(result.raw_p)
    return values, result


def tf_regulation(frame, alpha, lfc, pcolumn, top_n):
    if frame.gene.duplicated().any():
        raise ValueError("Duplicate transcription factors")
    numeric(frame, ["abundance"])
    frame = frame.sort_values("abundance", ascending=False).copy()
    frame["status"] = np.select([(frame[pcolumn] < alpha) & (frame.log2FoldChange < -lfc),
                                (frame[pcolumn] < alpha) & (frame.log2FoldChange > lfc)], ["Down", "Up"], "NC")
    records = []
    for name, subset in [("All TFs", frame), (f"Top {top_n}", frame.head(top_n))]:
        for status in ["Down", "NC", "Up"]:
            n = int((subset.status == status).sum())
            records.append(dict(subset=name, status=status, n=n, total=len(subset), percent=100 * n / len(subset)))
    return frame, pd.DataFrame(records)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("analysis", choices=["signal", "transcription-factors"])
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--alpha", type=float)
    p.add_argument("--lfc", type=float)
    p.add_argument("--pcolumn", choices=["pvalue", "padj"])
    p.add_argument("--top-n", type=int, default=100)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    if a.analysis == "transcription-factors" and any(x is None for x in [a.alpha, a.lfc, a.pcolumn]):
        p.error("transcription-factors requires --alpha, --lfc and --pcolumn")
    if a.analysis == "signal":
        values, summary = paired_signal(read_table(a.input, ["comparison", "region", "control", "treated"]))
    else:
        values, summary = tf_regulation(read_table(a.input, ["gene", "abundance", "log2FoldChange", a.pcolumn]), a.alpha, a.lfc, a.pcolumn, a.top_n)
    out = new_output(a.output)
    save_table(values, out / "values.tsv")
    save_table(summary, out / "summary.tsv")


if __name__ == "__main__":
    main()
