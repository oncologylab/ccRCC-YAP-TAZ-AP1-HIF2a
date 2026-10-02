import unittest
import reproduce_figure6d


class Figure6DTests(unittest.TestCase):
    def test_complete_original_result_families_and_display_values(self):
        result = reproduce_figure6d.verify()
        self.assertEqual(len(result), 15)
        self.assertEqual(set(result.tested_genes), {13833})
        self.assertEqual(set(result.n_per_condition), {3})
        self.assertTrue((result.adjusted_p < .0001).all())


if __name__ == '__main__':
    unittest.main()
