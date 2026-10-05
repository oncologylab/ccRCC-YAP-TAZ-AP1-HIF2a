from pathlib import Path
import importlib
import os
import subprocess
import sys
import tempfile
import unittest
import numpy as np
import pandas as pd

from figures._shared import DATA
from figures.figure_3.process import matched_peaks
from figures.figure_6.process import gene_scores, cluster_tests, verify_d
from figures.figure_EV7.process import downregulated

ROOT = Path(__file__).resolve().parents[1]


class FigureWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name)

    def cli(self, script, *args):
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", MPLCONFIGDIR=str(self.path / "mpl"))
        return subprocess.run([sys.executable, str(ROOT / script), *map(str, args)], check=True, capture_output=True, text=True, env=env)

    def test_all_figure_command_lines(self):
        for name in [str(i) for i in range(1, 7)] + [f"EV{i}" for i in range(1, 9)]:
            for operation in ["process", "plot"]:
                with self.subTest(figure=name, operation=operation):
                    self.cli(f"figures/figure_{name}/{operation}.py", "--help")

    def test_matched_peak_tests(self):
        frame = pd.read_csv(DATA / "14_baseline_matched_effects_values.tsv.gz", sep="\t")
        result = matched_peaks(frame).set_index("panel")
        expected = pd.read_csv(DATA / "14_baseline_matched_effects_stats.tsv", sep="\t").set_index("panel")
        np.testing.assert_array_equal(result.n_bound, expected.n_bound)
        np.testing.assert_array_equal(result.n_unbound, expected.n_unbound)
        np.testing.assert_allclose(result.raw_p, expected.raw_p, rtol=1e-8, atol=0)
        np.testing.assert_allclose(result.FDR_p, expected.FDR_p, rtol=1e-8, atol=0)

    def test_cluster_gene_counts_and_pairing(self):
        scores = gene_scores()
        self.assertFalse(scores.duplicated(["cluster", "group", "gene"]).any())
        result = cluster_tests(scores)
        self.assertEqual(result.paired_genes_n.to_list(), [702, 551, 983, 1365, 1095])
        expected = pd.read_csv(DATA / "13_shJUN_RNA_cluster_stats.tsv", sep="\t")
        np.testing.assert_allclose(result.FDR_p, expected.FDR_p, rtol=1e-8, atol=0)

    def test_deseq2_all_families(self):
        result = verify_d()
        self.assertEqual(len(result), 15)
        self.assertEqual(set(result.tested_genes), {13833})

    def test_gene_intersections_are_disjoint_and_complete(self):
        results = {name: pd.read_csv(DATA / f"figure6d_DESeq2_{name}.tsv.gz", sep="\t")
                   for name in ["shYT", "shHIF2a", "shJUNs"]}
        values, counts = downregulated(results, .05, 1, "padj")
        self.assertEqual(values.gene.nunique(), 13833)
        self.assertEqual(counts.n_genes.sum(), (values.membership != "None").sum())
        for name, frame in results.items():
            expected = set(frame.loc[(frame.padj < .05) & (frame.log2FoldChange < -1), "ENS_ID"])
            self.assertEqual(set(values.loc[values[name], "gene"]), expected)

    def test_figure3_cli_outputs_and_no_overwrite(self):
        out = self.path / "figure3"
        self.cli("figures/figure_3/process.py", "matched", "--output", out)
        self.cli("figures/figure_3/plot.py", "--kind", "matched", "--input", out / "values.tsv", "--output", out / "plot.pdf")
        self.assertGreater((out / "plot.pdf").stat().st_size, 1000)
        with self.assertRaises(subprocess.CalledProcessError):
            self.cli("figures/figure_3/process.py", "matched", "--output", out)

    def test_figure6_cli_and_plots(self):
        out = self.path / "figure6"
        self.cli("figures/figure_6/process.py", "--output", out)
        for panel in ["C", "D"]:
            self.cli("figures/figure_6/plot.py", "--panel", panel, "--input", out / f"figure6{panel}_values.tsv", "--output", out / f"{panel}.pdf")
            self.assertGreater((out / f"{panel}.pdf").stat().st_size, 1000)

    def test_missing_input_rejected_before_external_execution(self):
        out = self.path / "missing"
        with self.assertRaises(subprocess.CalledProcessError):
            self.cli("preprocessing/chromatin.py", "--execute", "peaks", "--bam", DATA / "absent.bam",
                     "--chrom-sizes", DATA / "absent.sizes", "--style", "factor", "--output", out)
        self.assertFalse(out.exists())


if __name__ == "__main__":
    unittest.main()
