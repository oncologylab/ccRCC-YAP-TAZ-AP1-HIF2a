from pathlib import Path
import argparse
import sys
import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import linkage, leaves_list

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from figures._shared import read_table, numeric, new_output, save_table
from figures.survival import kaplan_meier, logrank


def pdx_signature(path):
    frame = pd.read_excel(path, sheet_name="HIF2-i RNASeq YAP-TAZ Signature")
    frame = frame.rename(columns={"VehicleTx-SensitiveTumors": "VS", "VehicleTx-ResistantTumors": "VR",
                                  "HIF2aiTx-SensitiveTumors": "PS", "HIF2aiTx-ResistantTumors": "PR"})
    frame = frame[["VS", "VR", "PS", "PR"]].dropna(how="all")
    values = numeric(frame, frame.columns).T
    if frame.shape != (31, 4):
        raise ValueError("Expected 31 genes and four groups")
    ranks = stats.rankdata(values.ravel())
    means = ranks.reshape(values.shape).mean(axis=1)
    variance = len(ranks) * (len(ranks) + 1) / 12 * stats.tiecorrect(ranks)
    tests = []
    omnibus = stats.kruskal(*values)
    tests.append(dict(comparison="All groups", method="Kruskal-Wallis", statistic=omnibus.statistic, pvalue=omnibus.pvalue))
    for i, j in [(0, 1), (0, 2), (2, 3), (1, 3)]:
        z = (means[i] - means[j]) / np.sqrt(variance * (1 / len(values[i]) + 1 / len(values[j])))
        tests.append(dict(comparison=f"{frame.columns[i]} vs {frame.columns[j]}", method="Dunn uncorrected", statistic=z, pvalue=2 * stats.norm.sf(abs(z))))
    long = frame.rename_axis("gene_index").reset_index().melt(id_vars="gene_index", var_name="group", value_name="value")
    return long, pd.DataFrame(tests)


def baseline(frame):
    genes = ["YAP1", "WWTR1", "ARNT", "EPAS1", "TEAD1", "TEAD2", "TEAD4"]
    frame = frame.loc[frame.gene.isin(genes)].copy()
    numeric(frame, ["value"])
    if frame.duplicated(["gene", "sample"]).any():
        raise ValueError("Duplicate gene/sample")
    records = []
    for gene, part in frame.groupby("gene", sort=False):
        a, b = [part.loc[part.group == group, "value"] for group in ["Sensitive", "Resistant"]]
        if min(len(a), len(b)) == 0:
            raise ValueError("Sensitive and Resistant groups required")
        test = stats.mannwhitneyu(a, b, alternative="two-sided")
        records.append(dict(gene=gene, n_sensitive=len(a), n_resistant=len(b), pvalue=test.pvalue))
    return frame, pd.DataFrame(records)


def baseline_workbook(path):
    records = []
    sheets = {"YAP1": "YAP1", "WWTR1 (TAZ)": "WWTR1", "ARNT (HIF1β)": "ARNT",
              "EPAS1 (HIF2⍺)": "EPAS1", "TEAD1": "TEAD1", "TEAD2": "TEAD2", "TEAD4": "TEAD4"}
    for sheet, gene in sheets.items():
        frame = pd.read_excel(path, sheet_name=sheet)
        for group in ["Sensitive", "Resistant"]:
            records.extend(dict(gene=gene, group=group, sample=f"{group}_{index + 1}", value=float(value))
                           for index, value in frame[group].dropna().items())
    return baseline(pd.DataFrame(records))


def signature_heatmap(frame, method, metric):
    numeric(frame, frame.columns)
    z = frame.sub(frame.mean(axis=1), axis=0).div(frame.std(axis=1, ddof=1).replace(0, np.nan), axis=0).dropna()
    if min(z.shape) < 2:
        raise ValueError("At least two variable genes and two samples required")
    genes = leaves_list(linkage(z.to_numpy(), method=method, metric=metric))
    samples = leaves_list(linkage(z.to_numpy().T, method=method, metric=metric))
    return z.iloc[genes, samples]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("analysis", choices=["pdx", "baseline", "heatmap", "survival"])
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--method", choices=["average", "complete", "single"])
    p.add_argument("--metric", choices=["euclidean", "correlation"])
    a = p.parse_args()
    if a.analysis == "heatmap" and (not a.method or not a.metric):
        p.error("heatmap requires --method and --metric")
    out = new_output(a.output)
    if a.analysis in ["pdx", "baseline"]:
        if a.analysis == "pdx":
            values, tests = pdx_signature(a.input)
        elif a.input.suffix == ".xlsx":
            values, tests = baseline_workbook(a.input)
        else:
            values, tests = baseline(read_table(a.input, ["gene", "sample", "group", "value"]))
        save_table(values, out / "values.tsv")
        save_table(tests, out / "tests.tsv")
    elif a.analysis == "heatmap":
        z = signature_heatmap(pd.read_csv(a.input, sep="\t", index_col=0), a.method, a.metric)
        z.to_csv(out / "heatmap.tsv", sep="\t")
    else:
        frame = read_table(a.input, ["sample", "group", "time", "event"])
        save_table(kaplan_meier(frame), out / "survival.tsv")
        save_table(logrank(frame), out / "tests.tsv")


if __name__ == "__main__":
    main()
