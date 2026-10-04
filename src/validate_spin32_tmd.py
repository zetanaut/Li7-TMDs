#!/usr/bin/env python3
"""Legacy CLI for exact spin-3/2, STF, and selected SIDIS checks."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
if __package__:
    from .spin32_symbolic import run_checks
else:
    from spin32_symbolic import run_checks


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('results/validation_report.json'),
                        help='Path for the JSON report (default: results/validation_report.json).')
    args = parser.parse_args()
    checks, _ranks = run_checks(emit=lambda line: print(line, flush=True))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(checks, indent=2) + '\n', encoding='utf-8')
    print(f'All {sum(v == "PASS" for v in checks.values())} checks passed. Report: {args.output}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
