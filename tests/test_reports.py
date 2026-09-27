"""Lightweight checks of the generated scientific artifacts."""
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]


class ReportTests(unittest.TestCase):
    def test_symbolic_report(self):
        d = json.loads((ROOT/'results/validation_report.json').read_text())
        self.assertEqual(sum(v == 'PASS' for v in d.values()), 49)
        for sector in ('quark', 'gluon'):
            cat = d[f'{sector}_catalogue']
            self.assertEqual(len(cat), 32)
            self.assertEqual([sum(row['K'] == k for row in cat) for k in range(4)], [2, 6, 10, 14])
            self.assertEqual(d[f'{sector} 64-by-32 tensor map has rank 32'], 'PASS')

    def test_born_report(self):
        d = json.loads((ROOT/'results/gluon_born_report.json').read_text())
        self.assertTrue(all(d['checks'].values()))
        self.assertTrue(all(v < 1e-11 for v in d['relative_Ward_residuals'].values()))
        self.assertTrue(all(v >= -1e-11 for v in d['B_eigenvalues']))
        br, bi = d['B_real'], d['B_imag']
        self.assertAlmostEqual(br[0][1], br[1][0], places=9)
        self.assertAlmostEqual(bi[0][1], -bi[1][0], places=9)

    def test_grid(self):
        d = json.loads((ROOT/'results/gluon_born_grid_report.json').read_text())
        self.assertEqual(d['points'], 36)
        self.assertEqual(len(d['results']), 36)
        self.assertTrue(d['all_pass'])
        self.assertTrue(all(all(row['checks'].values()) for row in d['results']))


if __name__ == '__main__':
    unittest.main()
