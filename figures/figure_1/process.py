from pathlib import Path
import argparse
import sys
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from figures._shared import read_table, numeric, new_output, save_table
from figures.survival import kaplan_meier, logrank


def protein(workbook):
    values, tests = [], []
    for gene in ["SAV1", "LATS1"]:
        frame = pd.read_excel(workbook, sheet_name=gene)
        frame = frame.loc[frame.CCRCC == "Yes"].dropna(subset=[f"{gene}.CNV", f"{gene}.z.Protein"])
        cnv = frame[f"{gene}.CNV"]
        selected = frame.loc[(cnv < -.2) | ((cnv > -.2) & (cnv < .2))].copy()
        selected["group"] = np.where(selected[f"{gene}.CNV"] < -.2, "Loss", "Neutral")
        selected = selected.rename(columns={"SAMPLE": "sample", f"{gene}.z.Protein": "value"})
        selected["gene"] = gene
        neutral = selected.loc[selected.group == "Neutral", "value"].to_numpy(float)
        loss = selected.loc[selected.group == "Loss", "value"].to_numpy(float)
        numeric(selected, ["value"])
        test = stats.mannwhitneyu(neutral, loss, alternative="two-sided", method="exact")
        tests.append(dict(gene=gene, n_neutral=len(neutral), n_loss=len(loss), statistic=test.statistic, pvalue=test.pvalue))
        values.append(selected[["gene", "sample", "group", "value"]])
    return pd.concat(values, ignore_index=True), pd.DataFrame(tests)


def cna(frame):
    numeric(frame, ["cnv"])
    if frame.duplicated(["cohort", "gene", "sample"]).any():
        raise ValueError("Duplicate sample/gene")
    frame = frame.assign(status=np.select([frame.cnv < -.2, frame.cnv > .2], ["Loss", "Gain"], "Neutral"))
    counts = frame.groupby(["cohort", "gene", "status"]).size().rename("n").reset_index()
    counts["percent"] = counts.n / counts.groupby(["cohort", "gene"]).n.transform("sum") * 100
    return counts


def expression(frame):
    numeric(frame, ["cnv", "value"])
    frame = frame.loc[(frame.cnv < -.2) | ((frame.cnv > -.2) & (frame.cnv < .2))].copy()
    frame["group"] = np.where(frame.cnv < -.2, "Loss", "Neutral")
    records = []
    for (cohort, gene), part in frame.groupby(["cohort", "gene"]):
        a, b = [part.loc[part.group == g, "value"].to_numpy(float) for g in ["Neutral", "Loss"]]
        if min(len(a), len(b)) < 2:
            raise ValueError("At least two observations in each group")
        if cohort == "CPTAC" and gene in ["SAV1", "LATS1"]:
            equal = gene == "LATS1"
            test = stats.ttest_ind(a, b, equal_var=equal)
            method = "Student" if equal else "Welch"
        else:
            test = stats.mannwhitneyu(a, b, alternative="two-sided")
            method = "Mann-Whitney"
        records.append(dict(cohort=cohort, gene=gene, n_neutral=len(a), n_loss=len(b), method=method, pvalue=test.pvalue))
    return frame, pd.DataFrame(records)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("analysis", choices=["protein", "cna", "expression", "survival"])
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    out = new_output(a.output)
    if a.analysis == "protein":
        values, tests = protein(a.input)
        save_table(values, out / "values.tsv")
        save_table(tests, out / "tests.tsv")
    elif a.analysis == "cna":
        save_table(cna(read_table(a.input, ["cohort", "gene", "sample", "cnv"])), out / "cna.tsv")
    elif a.analysis == "expression":
        values, tests = expression(read_table(a.input, ["cohort", "gene", "sample", "cnv", "value"]))
        save_table(values, out / "values.tsv")
        save_table(tests, out / "tests.tsv")
    else:
        frame = read_table(a.input, ["sample", "group", "time", "event"])
        save_table(kaplan_meier(frame), out / "survival.tsv")
        save_table(logrank(frame), out / "tests.tsv")


if __name__ == "__main__":
    main()
