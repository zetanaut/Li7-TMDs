"""Regenerate and validate the 4 × 3 × 3 Born-response grid."""
import argparse
import itertools
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from gluon_born_response import evaluate  # noqa: E402


def evaluate_grid() -> dict[str, object]:
    """Evaluate the documented 4 × 3 × 3 physical grid without writing."""
    rows = []
    for theta, phi, helicity in itertools.product(
            [0.4, 0.8, 1.7, 2.7], [0.0, 0.4, 1.2], [-1.0, 0.0, 1.0]):
        report = evaluate(theta=theta, phi=phi, helicity=helicity)
        rows.append({key: report[key] for key in (
            'convention','convention_digest','reference_source_sha256',
            'inputs', 'Stokes', 'B_eigenvalues', 'relative_Ward_residuals',
            'relative_Hermiticity_residual', 'checks')})
    all_pass = len(rows) == 36 and all(all(row['checks'].values()) for row in rows)
    return {'convention':'born-current-v2-physical-h',
            'points': len(rows), 'all_pass': all_pass, 'results': rows}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'results/gluon_born_grid_report.json')
    args = parser.parse_args()
    report = evaluate_grid()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(f"{report['points']} points; all_pass={report['all_pass']}; report: {args.output}")
    return 0 if report['all_pass'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
