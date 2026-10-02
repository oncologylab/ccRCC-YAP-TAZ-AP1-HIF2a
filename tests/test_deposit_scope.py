from pathlib import Path
import ast
import unittest

ROOT = Path(__file__).resolve().parents[1]


class DepositScopeTests(unittest.TestCase):
    def test_readme(self):
        self.assertEqual((ROOT / "README.md").read_text().strip(),
                         "Code for “YAP/TAZ–AP-1 Cooperate with HIF2α to Drive Oncogenic Transcription in ccRCC.”")

    def test_no_archives_examples_or_review_notes(self):
        for name in ["historical_sources", "examples", "docs"]:
            self.assertFalse(any(path.is_file() for path in (ROOT / name).rglob("*")))
        self.assertFalse((ROOT / "figures/common.py").exists())

    def test_figure_files_have_local_analysis_and_plotting(self):
        for name in [str(i) for i in range(1, 7)] + [f"EV{i}" for i in range(1, 9)]:
            for operation in ["process", "plot"]:
                path = ROOT / "figures" / f"figure_{name}" / f"{operation}.py"
                tree = ast.parse(path.read_text())
                functions = [node for node in tree.body if isinstance(node, ast.FunctionDef)]
                self.assertTrue(any(node.end_lineno - node.lineno >= 15 for node in functions), str(path))
                self.assertNotIn("process_main(", path.read_text())
                self.assertNotIn("plot_main(", path.read_text())


if __name__ == "__main__":
    unittest.main()
