from pathlib import Path
import hashlib
import re
import shutil
import subprocess
import unittest

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "preprocessing/image_analysis"
FILES = {
    "Part1_Isolate_Nuclei_from_Confocal_Tiffs_and_Save_Combined_Images.ijm":
        "835f3552341782895ece894380cea85b613c81c644642b4b539f944ad0b910f2",
    "Part2_Extract_Puncta_Information_From_Two_Channel_Combined_Image.ijm":
        "e80f653d86995b4fa1274873a79ad2f6e1304dc1b67b0773cd87b5fc9ed593ee",
    "Part3_R_Measuring_2Ch_puncta_Colocalization_from_ImageJ_Output.Rmd":
        "f59fa44054f0d878e1aeff28e16c072e0f2cfb31c6ed9c36c67b333cab44ad37",
}


class ImageAnalysisTests(unittest.TestCase):
    def test_supplied_files_unchanged(self):
        self.assertEqual({p.name for p in SOURCE.iterdir()}, set(FILES))
        for name, expected in FILES.items():
            self.assertEqual(hashlib.sha256((SOURCE / name).read_bytes()).hexdigest(), expected, name)

    def test_original_attribution(self):
        notebook = SOURCE / next(name for name in FILES if name.endswith(".Rmd"))
        self.assertIn('author: "Alec McIntosh"', notebook.read_text())
        self.assertIn("preprocessing/image_analysis/: Alec McIntosh.", (ROOT / "NOTICE").read_text())

    def test_macro_delimiters(self):
        pattern = r'"(?:\\.|[^"\\])*"|//[^\n]*|/\*.*?\*/|[(){}\[\]]'
        pairs = {')': '(', '}': '{', ']': '['}
        for path in SOURCE.glob("*.ijm"):
            stack = []
            for token in re.findall(pattern, path.read_text(), flags=re.DOTALL):
                if token in ('(', '{', '['):
                    stack.append(token)
                elif token in pairs:
                    self.assertTrue(stack, path.name)
                    self.assertEqual(stack.pop(), pairs[token], path.name)
            self.assertEqual(stack, [], path.name)

    @unittest.skipUnless(shutil.which("Rscript"), "Rscript unavailable")
    def test_r_chunks_parse_without_execution(self):
        notebook = SOURCE / next(name for name in FILES if name.endswith(".Rmd"))
        code = r'''
lines <- readLines(commandArgs(trailingOnly = TRUE)[1])
active <- FALSE
chunks <- list()
code <- character()
for (line in lines) {
    if (grepl("^```\\{r", line)) {
        stopifnot(!active)
        active <- TRUE
        code <- character()
    } else if (active && grepl("^```", line)) {
        chunks[[length(chunks) + 1]] <- code
        active <- FALSE
    } else if (active) {
        code <- c(code, line)
    }
}
stopifnot(!active, length(chunks) == 3)
for (chunk in chunks) parse(text = chunk)
'''
        subprocess.run(["Rscript", "--vanilla", "-e", code, str(notebook)],
                       check=True, capture_output=True, text=True)


if __name__ == "__main__":
    unittest.main()
