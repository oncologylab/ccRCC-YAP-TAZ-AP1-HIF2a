from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from figures._shared import ROOT, DATA, bh, new_output, save_table
from figures.figure_3.process import matched_peaks
from figures.figure_6.process import gene_scores, cluster_tests


def chromatin():
    return matched_peaks(pd.read_csv(DATA / "14_baseline_matched_effects_values.tsv.gz", sep="\t"))


def rna():
    return cluster_tests(gene_scores())


def validate(actual, reference, id_col, counts):
    expected = pd.read_csv(reference, sep="\t").set_index(id_col)
    actual = actual.set_index(id_col).loc[expected.index]
    for column in counts:
        np.testing.assert_array_equal(actual[column], expected[column])
    for column in ["raw_p", "FDR_p"]:
        np.testing.assert_allclose(actual[column], expected[column], rtol=1e-8, atol=0)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    c, r = chromatin(), rna()
    validate(c, DATA / "14_baseline_matched_effects_stats.tsv", "panel", ["n_bound", "n_unbound"])
    validate(r, DATA / "13_shJUN_RNA_cluster_stats.tsv", "cluster", ["paired_genes_n"])
    out = new_output(a.output)
    save_table(c, out / "figure3E_tests.tsv")
    save_table(r, out / "figure6C_tests.tsv")
    print("PASS: 11 comparisons.")
