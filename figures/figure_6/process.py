from pathlib import Path
import argparse
import sys
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from figures._shared import DATA, bh, new_output, save_table
from figures.figure_4.process import prism_values

GENES = ["VEGFA", "SLC2A1", "MTHFD2", "CCND1", "MYC"]
CONTRASTS = ["shYT", "shHIF2a", "shJUNs"]


def gene_scores(data=DATA):
    columns = [f"786O_Ctrl_{i}" for i in range(1, 4)] + [f"786O_shJUN_{i}" for i in range(1, 4)]
    counts = pd.read_csv(data / "figure6c_normalized_counts.tsv", sep="\t", index_col="ensembl_gene_id")[columns]
    lognorm = np.log2(counts + 1)
    z = lognorm.sub(lognorm.mean(axis=1), axis=0).div(lognorm.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    definitions = pd.read_csv(data / "figure6c_cluster_definitions.tsv", sep="\t")
    definitions["gene"] = definitions.rows.str.replace(r"\..*$", "", regex=True)
    records = []
    for revised, original in {1: 3, 2: 5, 3: 2, 4: 1, 5: 4}.items():
        genes = pd.Index(definitions.loc[definitions.number == original, "gene"].drop_duplicates()).intersection(z.dropna().index)
        for group, cols in [("Ctrl", columns[:3]), ("shJUN", columns[3:])]:
            records.extend(dict(cluster=revised, group=group, gene=gene, value=value)
                           for gene, value in z.loc[genes, cols].mean(axis=1).items())
    return pd.DataFrame(records)


def cluster_tests(scores):
    records = []
    for cluster, part in scores.groupby("cluster"):
        paired = part.pivot(index="gene", columns="group", values="value")
        if paired.isna().any().any():
            raise ValueError("Incomplete gene pairs")
        test = stats.wilcoxon(paired.Ctrl, paired.shJUN, alternative="two-sided")
        records.append(dict(cluster=cluster, paired_genes_n=len(paired), raw_p=test.pvalue))
    result = pd.DataFrame(records)
    result["FDR_p"] = bh(result.raw_p)
    return result


def verify_d(data=DATA):
    expected = pd.read_csv(data / "figure6d_displayed_values.tsv", sep="\t")
    normalized = pd.read_csv(data / "figure6d_normalized_counts.tsv", sep="\t").set_index("gene")
    records = []
    for contrast in CONTRASTS:
        frame = pd.read_csv(data / f"figure6d_DESeq2_{contrast}.tsv.gz", sep="\t")
        if len(frame) != 13833 or frame.ENS_ID.nunique() != 13833 or frame.pvalue.isna().any():
            raise ValueError("Expected complete 13,833-gene result family")
        np.testing.assert_allclose(2 * stats.norm.sf(abs(frame.stat)), frame.pvalue, rtol=1e-10, atol=1e-250)
        np.testing.assert_allclose(bh(frame.pvalue), frame.padj, rtol=2e-12, atol=1e-250)
        results = frame.set_index("ENS_ID")
        for gene in GENES:
            row = normalized.loc[gene]
            selected = expected.loc[(expected.gene == gene) & (expected.contrast == contrast)].sort_values("replicate")
            if selected.replicate.to_list() != [1, 2, 3]:
                raise ValueError("Expected three replicates")
            controls = row[[f"PAR.786O.C{i}" for i in range(1, 4)]].to_numpy(float)
            knockdown = row[[f"{contrast}.786O.C{i}" for i in range(1, 4)]].to_numpy(float)
            values = np.log2(knockdown / controls.mean())
            np.testing.assert_allclose(values, selected.displayed_log2FC, atol=5.1e-8, rtol=0)
            saved = results.loc[row.ENS_ID]
            if saved.V2 != gene:
                raise ValueError("Gene annotation mismatch")
            np.testing.assert_allclose(saved.padj, selected.adjusted_p, rtol=1e-12, atol=0)
            records.append(dict(gene=gene, contrast=contrast, n_per_condition=3, tested_genes=len(frame),
                                adjusted_p=float(saved.padj), mean_displayed_log2FC=values.mean(), displayed_SD=values.std(ddof=1)))
    return pd.DataFrame(records)


def imaging(path):
    ctrl, treated = prism_values(path, "Table11")
    test = stats.ttest_ind(ctrl, treated, equal_var=False)
    records = [dict(group=group, nucleus=i + 1, value=value)
               for group, values in [("Ctrl", ctrl), ("shJUNs", treated)] for i, value in enumerate(values)]
    return pd.DataFrame(records), pd.DataFrame([dict(n_ctrl=len(ctrl), n_shJUNs=len(treated), method="Welch", pvalue=test.pvalue)])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, default=DATA)
    p.add_argument("--prism", type=Path)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    scores = gene_scores(a.data)
    c, d = cluster_tests(scores), verify_d(a.data)
    out = new_output(a.output)
    save_table(scores, out / "figure6C_values.tsv")
    save_table(c, out / "figure6C_tests.tsv")
    save_table(d, out / "figure6D_tests.tsv")
    save_table(pd.read_csv(a.data / "figure6d_displayed_values.tsv", sep="\t"), out / "figure6D_values.tsv")
    if a.prism:
        values, tests = imaging(a.prism)
        save_table(values, out / "figure6B_values.tsv")
        save_table(tests, out / "figure6B_tests.tsv")


if __name__ == "__main__":
    main()
