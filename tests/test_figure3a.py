from pathlib import Path
import csv
import math
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/figure3a_tumor_growth.csv"
SCRIPT = ROOT / "figures/figure_3/tumor_growth.R"
EXPECTED = {
    "0_1": (0.0263552725310127, 204, 17, 15, 37, "*"),
    "0_4": (0.0003224836438281, 252, 21, 19, 37, "***"),
    "2_3": (0.0312009123945325, 212, 15, 13, 51, "*"),
    "2_5": (0.349428734490235, 213, 14, 12, 51, "ns"),
}


class Figure3ADataTests(unittest.TestCase):
    def test_complete_source_measurements(self):
        with DATA.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 662)
        self.assertEqual({int(row["group"]) for row in rows}, set(range(6)))
        self.assertEqual(len({(row["mouse"], row["day"]) for row in rows}), 662)
        self.assertTrue(all(math.isfinite(float(row["log2_fold_change"])) for row in rows))
        day37 = [row for row in rows if int(row["group"]) in (2, 5) and float(row["day"]) == 37]
        self.assertEqual(len(day37), 12)


@unittest.skipUnless(shutil.which("Rscript"), "Rscript unavailable")
class Figure3AAnalysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        dependency = subprocess.run(
            ["Rscript", "--vanilla", "-e", "quit(status=if(requireNamespace('nlme', quietly=TRUE)) 0 else 1)"],
            capture_output=True, text=True, check=False)
        if dependency.returncode:
            raise unittest.SkipTest("R package nlme unavailable")
        cls.workspace = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.workspace.cleanup)
        cls.output = Path(cls.workspace.name) / "results"
        subprocess.run(["Rscript", "--vanilla", str(SCRIPT), str(DATA), str(cls.output)],
                       check=True, capture_output=True, text=True)
        with (cls.output / "treatment_effects.tsv").open(newline="") as stream:
            cls.results = {row["comparison"]: row for row in csv.DictReader(stream, delimiter="\t")}

    def test_all_four_final_treatment_effects(self):
        self.assertEqual(set(self.results), set(EXPECTED))
        for comparison, (p, n, mice, df, day_max, stars) in EXPECTED.items():
            with self.subTest(comparison=comparison):
                row = self.results[comparison]
                self.assertTrue(math.isclose(float(row["p_value"]), p, rel_tol=1e-7, abs_tol=1e-12))
                self.assertEqual(int(row["n_observations"]), n)
                self.assertEqual(int(row["n_mice"]), mice)
                self.assertEqual(int(row["df_num"]), 1)
                self.assertEqual(int(row["df_den"]), df)
                self.assertEqual(float(row["day_min"]), 0)
                self.assertEqual(float(row["day_max"]), day_max)
                self.assertEqual(row["annotation"], stars)
                self.assertEqual(row["effect"], "Treatment main effect")
                self.assertEqual(row["adjustment"], "None (separate pairwise models)")
                method = "Two-way repeated-measures ANOVA" if comparison.startswith("0_") else "Additive REML mixed-effects model"
                self.assertEqual(row["method"], method)

    def test_models_and_software_versions_saved(self):
        self.assertTrue((self.output / "models.rds").stat().st_size > 0)
        session = (self.output / "sessionInfo.txt").read_text()
        self.assertIn("R version", session)
        self.assertIn("nlme", session)

    def test_existing_output_is_not_overwritten(self):
        before = (self.output / "treatment_effects.tsv").read_bytes()
        result = subprocess.run(["Rscript", "--vanilla", str(SCRIPT), str(DATA), str(self.output)],
                                capture_output=True, text=True, check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Choose a new output directory", result.stderr)
        self.assertEqual((self.output / "treatment_effects.tsv").read_bytes(), before)


if __name__ == "__main__":
    unittest.main()
