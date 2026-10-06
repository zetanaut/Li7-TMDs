"""Execute the legacy scientific calculations rather than trusting stored JSON."""
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'examples'))
from spin32_symbolic import run_checks  # noqa: E402
from validation_manifest import LEGACY_SYMBOLIC_LABELS  # noqa: E402
from gluon_born_response import evaluate  # noqa: E402
from reproduce_gluon_grid import evaluate_grid  # noqa: E402


class ReportTests(unittest.TestCase):
    def test_symbolic_report(self):
        report, ranks = run_checks()
        self.assertEqual([key for key, value in report.items() if value == 'PASS'],
                         list(LEGACY_SYMBOLIC_LABELS))
        self.assertEqual(ranks, {'quark': 32, 'gluon': 32})
        for sector in ('quark', 'gluon'):
            catalogue = report[f'{sector}_catalogue']
            self.assertEqual(len(catalogue), 32)
            self.assertEqual([sum(row['K'] == k for row in catalogue) for k in range(4)],
                             [2, 6, 10, 14])

    def test_born_report(self):
        report = evaluate()
        self.assertTrue(all(report['checks'].values()))
        self.assertTrue(all(value < 1e-11 for value in report['relative_Ward_residuals'].values()))
        self.assertTrue(all(value >= -1e-11 for value in report['B_eigenvalues']))
        real, imag = report['B_real'], report['B_imag']
        self.assertAlmostEqual(real[0][1], real[1][0], places=9)
        self.assertAlmostEqual(imag[0][1], -imag[1][0], places=9)

    def test_grid(self):
        report = evaluate_grid()
        self.assertEqual(report['points'], 36)
        self.assertEqual(len(report['results']), 36)
        self.assertTrue(report['all_pass'])
        self.assertTrue(all(all(row['checks'].values()) for row in report['results']))


if __name__ == '__main__':
    unittest.main()
