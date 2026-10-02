from pathlib import Path
import re
import shutil
import subprocess
import unittest

ROOT = Path(__file__).resolve().parents[1]
EMAILS = {"yl814@georgetown.edu", "18198096+YaoxiangLi@users.noreply.github.com"}


class AuthorshipTests(unittest.TestCase):
    def test_sole_software_author(self):
        citation = (ROOT / "CITATION.cff").read_text()
        authors = citation.split("authors:\n", 1)[1].split("\nrepository-code:", 1)[0]
        self.assertEqual(len(re.findall(r"^  - ", authors, re.MULTILINE)), 1)
        self.assertIn('family-names: "Li"', authors)
        self.assertIn('given-names: "Yaoxiang"', authors)
        self.assertIn('url: "https://github.com/YaoxiangLi"', authors)

    def test_copyright_holder(self):
        license_text = (ROOT / "LICENSE").read_text()
        self.assertIn("Copyright (c) 2026 Yaoxiang Li", license_text)
        self.assertNotIn("Yi Laboratory contributors", license_text)

    @unittest.skipUnless((ROOT / ".git").exists() and shutil.which("git"), "Git history unavailable in source ZIP")
    def test_commit_identities_and_no_coauthors(self):
        output = subprocess.check_output(
            ["git", "log", "--all", "--format=%an%x09%ae%x09%cn%x09%ce"],
            cwd=ROOT, text=True,
        )
        for row in output.splitlines():
            author, author_email, committer, committer_email = row.split("\t")
            self.assertEqual(author, "Yaoxiang Li")
            self.assertEqual(committer, "Yaoxiang Li")
            self.assertIn(author_email, EMAILS)
            self.assertIn(committer_email, EMAILS)
        messages = subprocess.check_output(["git", "log", "--all", "--format=%B"], cwd=ROOT, text=True)
        self.assertIsNone(re.search(r"^Co-authored-by:", messages, re.IGNORECASE | re.MULTILINE))


if __name__ == "__main__":
    unittest.main()
