from pathlib import Path
import argparse
import sys
import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import linkage, leaves_list

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from figures._shared import read_table, numeric, new_output, save_table


def methylation(frame):
    frame = frame.dropna(subset=["methylation_beta", "expression_z"]).copy()
    numeric(frame, ["methylation_beta", "expression_z"])
    if len(frame) < 3 or frame["sample"].duplicated().any():
        raise ValueError("At least three unique samples required")
    if not frame.methylation_beta.between(0, 1).all():
        raise ValueError("Methylation beta must be between 0 and 1")
    records = []
    for method, function in [("Pearson", stats.pearsonr), ("Spearman", stats.spearmanr)]:
        result = function(frame.methylation_beta, frame.expression_z)
        records.append(dict(method=method, n=len(frame), correlation=result.statistic, pvalue=result.pvalue))
    return frame, pd.DataFrame(records)


def select_signature(effects, ranks, maximum_rank, cutoff):
    if effects.gene.duplicated().any() or ranks.gene.duplicated().any():
        raise ValueError("Duplicate gene symbols")
    joined = effects.merge(ranks[["gene", "rank"]], on="gene", validate="one_to_one")
    numeric(joined, ["log2FC", "rank"])
    return joined.loc[(joined.log2FC < cutoff) & (joined["rank"] <= maximum_rank)].sort_values("rank")


def heatmap(matrix, genes, method, metric):
    if matrix.index.duplicated().any():
        raise ValueError("Duplicate expression genes")
    missing = set(genes) - set(matrix.index)
    if missing:
        raise ValueError(f"Missing signature genes: {sorted(missing)}")
    selected = matrix.loc[genes]
    numeric(selected, selected.columns)
    z = selected.sub(selected.mean(axis=1), axis=0).div(selected.std(axis=1, ddof=1).replace(0, np.nan), axis=0).dropna()
    row_order = leaves_list(linkage(z.to_numpy(), method=method, metric=metric))
    col_order = leaves_list(linkage(z.to_numpy().T, method=method, metric=metric))
    return z.iloc[row_order, col_order]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("analysis", choices=["methylation", "signature", "heatmap"])
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--ranks", type=Path)
    p.add_argument("--max-rank", type=int)
    p.add_argument("--log2fc-cutoff", type=float)
    p.add_argument("--genes", type=Path)
    p.add_argument("--method", choices=["average", "complete", "single"])
    p.add_argument("--metric", choices=["euclidean", "correlation"])
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    if a.analysis == "signature" and any(x is None for x in [a.ranks, a.max_rank, a.log2fc_cutoff]):
        p.error("signature requires --ranks, --max-rank and --log2fc-cutoff")
    if a.analysis == "heatmap" and not all([a.genes, a.method, a.metric]):
        p.error("heatmap requires --genes, --method and --metric")
    out = new_output(a.output)
    if a.analysis == "methylation":
        values, tests = methylation(read_table(a.input, ["sample", "methylation_beta", "expression_z"]))
        save_table(values, out / "values.tsv")
        save_table(tests, out / "correlations.tsv")
    elif a.analysis == "signature":
        selected = select_signature(read_table(a.input, ["gene", "log2FC"]), read_table(a.ranks, ["gene", "rank"]), a.max_rank, a.log2fc_cutoff)
        save_table(selected, out / "signature.tsv")
    else:
        matrix = pd.read_csv(a.input, sep="\t", index_col=0)
        genes = read_table(a.genes, ["gene"]).gene.drop_duplicates().to_list()
        heatmap(matrix, genes, a.method, a.metric).to_csv(out / "heatmap.tsv", sep="\t")


if __name__ == "__main__":
    main()
