from pathlib import Path
import argparse
import sys
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from figures._shared import read_table, numeric, new_output, save_table


def clinical(frame):
    if frame["sample"].duplicated().any() or frame[["group", "category"]].isna().any().any():
        raise ValueError("Duplicate samples or missing annotations")
    counts = pd.crosstab(frame.group, frame.category)
    if len(counts) != 2:
        raise ValueError("Two signature groups required")
    percentages = counts.div(counts.sum(axis=1), axis=0) * 100
    records = []
    for category in counts.columns:
        table = np.column_stack([counts[category], counts.sum(axis=1) - counts[category]])
        result = fisher_exact(table)
        records.append(dict(category=category, odds_ratio=result.statistic, pvalue=result.pvalue))
    long = counts.rename_axis("group").reset_index().melt(id_vars="group", var_name="category", value_name="n")
    long["percent"] = [percentages.loc[g, c] for g, c in zip(long.group, long.category)]
    return long, pd.DataFrame(records)


def expression(counts, results, clusters, alpha, lfc, pcolumn):
    numeric(counts, counts.columns)
    if (counts < 0).any().any() or counts.index.duplicated().any():
        raise ValueError("Nonnegative counts and unique genes required")
    selected = set()
    for frame in results:
        selected.update(frame.loc[(frame[pcolumn] < alpha) & (frame.log2FoldChange.abs() > lfc), "ENS_ID"])
    missing = selected - set(counts.index)
    if missing:
        raise ValueError(f"{len(missing)} selected genes missing from counts")
    if clusters.gene.duplicated().any():
        raise ValueError("Duplicate cluster assignments")
    definitions = clusters.loc[clusters.gene.isin(selected)].sort_values(["cluster", "gene"])
    if set(definitions.gene) != selected:
        raise ValueError("Cluster assignments required for every selected gene")
    logged = np.log2(counts.loc[definitions.gene] + 1)
    z = logged.sub(logged.mean(axis=1), axis=0).div(logged.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    if z.isna().any().any():
        raise ValueError("Constant-expression gene")
    return z, definitions


def main():
    p = argparse.ArgumentParser()
    p.add_argument("analysis", choices=["clinical", "expression"])
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--results", type=Path, nargs="+")
    p.add_argument("--clusters", type=Path)
    p.add_argument("--alpha", type=float)
    p.add_argument("--lfc", type=float)
    p.add_argument("--pcolumn", choices=["pvalue", "padj"])
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    if a.analysis == "expression" and any(x is None for x in [a.results, a.clusters, a.alpha, a.lfc, a.pcolumn]):
        p.error("expression requires --results, --clusters, --alpha, --lfc and --pcolumn")
    out = new_output(a.output)
    if a.analysis == "clinical":
        values, tests = clinical(read_table(a.input, ["sample", "group", "category"]))
        save_table(values, out / "clinical.tsv")
        save_table(tests, out / "tests.tsv")
    else:
        z, definitions = expression(pd.read_csv(a.input, sep="\t", index_col=0),
                                    [read_table(x, ["ENS_ID", "log2FoldChange", a.pcolumn]) for x in a.results],
                                    read_table(a.clusters, ["gene", "cluster"]), a.alpha, a.lfc, a.pcolumn)
        z.to_csv(out / "heatmap.tsv", sep="\t")
        save_table(definitions, out / "clusters.tsv")


if __name__ == "__main__":
    main()
