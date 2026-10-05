import unittest
import numpy as np
from reproduce import ROOT, bh, chromatin, rna, validate

class ReproductionTests(unittest.TestCase):
    def test_bh_known_values(self):
        np.testing.assert_allclose(bh([.01, .04, .03]), [.03, .04, .04])

    def test_chromatin_source(self):
        validate(chromatin(), ROOT / 'data/14_baseline_matched_effects_stats.tsv', 'panel', ['n_bound','n_unbound'])

    def test_rna_source(self):
        validate(rna(), ROOT / 'data/13_shJUN_RNA_cluster_stats.tsv', 'cluster', ['paired_genes_n'])

if __name__ == '__main__':
    unittest.main()
