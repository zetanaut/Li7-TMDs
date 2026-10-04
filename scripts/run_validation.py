#!/usr/bin/env python3
"""Execute a run-scoped baseline; full reports pending scientific suites."""
from __future__ import annotations
import argparse
import math
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'examples'))
from spin32_symbolic import run_checks  # noqa: E402
from gluon_born_response import evaluate  # noqa: E402
from reproduce_gluon_grid import evaluate_grid  # noqa: E402
from validation_manifest import (BORN_CHECKS, FULL_PENDING_SUITES, LEGACY_IDS,
                                 LEGACY_SYMBOLIC_LABELS, SOFTWARE_TESTS)  # noqa: E402
from validation_evidence import (GRID_AXES, REFERENCE_INPUTS, EvidenceError, finish_run,
                                 new_run, result, validate_run,
                                 write_results)  # noqa: E402


def execute(run_dir: Path, manifest: dict, reference_inputs: dict,
            results: list[dict], ranks: dict[str, int]) -> None:
    legacy, computed_ranks = run_checks()
    ranks.update(computed_ranks)
    if [key for key, value in legacy.items() if value == 'PASS'] != list(LEGACY_SYMBOLIC_LABELS):
        raise EvidenceError('legacy symbolic identities differ from required manifest')
    for label in LEGACY_SYMBOLIC_LABELS:
        payload = {'legacy_label': label}
        for kind in ('quark', 'gluon'):
            if label == f'{kind} 64-by-32 tensor map has rank 32':
                payload['computed_rank'] = ranks[kind]
                payload['matrix_shape'] = [64, len(legacy[f'{kind}_catalogue'])]
            if label == f'{kind} catalogue has 32 entries':
                catalogue = legacy[f'{kind}_catalogue']
                payload['catalogue'] = catalogue
                payload['catalogue_size'] = len(catalogue)
                payload['target_rank_counts'] = [sum(row['K'] == k for row in catalogue)
                                                 for k in range(4)]
        results.append(result(LEGACY_IDS[label], 'PASS',
                              'exact_rank' if 'computed_rank' in payload else 'exact_identity',
                              ['Exact finite-dimensional SymPy algebra; legacy convention'],
                              payload, manifest))
    write_results(run_dir, manifest, results)

    born = evaluate(**reference_inputs)
    if born['inputs'] != reference_inputs or set(born['checks']) != set(BORN_CHECKS):
        raise EvidenceError('Born report input or check mismatch')
    for name in BORN_CHECKS:
        results.append(result('born.' + name, 'PASS' if born['checks'][name] else 'FAIL',
                              'numerical_diagnostic', ['Specified coupling-stripped Born kinematics'],
                              {'inputs': born['inputs'], 'report': born}, manifest))
    write_results(run_dir, manifest, results)

    grid = evaluate_grid()
    cases = grid['results']
    expected_cases = math.prod(len(axis) for axis in GRID_AXES.values())
    grid_ok = (grid['all_pass'] and grid['points'] == expected_cases and
               len(cases) == expected_cases and
               all(set(row['checks']) == set(BORN_CHECKS) and
                   all(row['checks'].values()) for row in cases))
    results.append(result('grid.complete', 'PASS' if grid_ok else 'FAIL',
                          'numerical_parameter_cases', ['Specified 4 × 3 × 3 physical grid'],
                          {'number_of_cases': grid['points'], 'cases': cases}, manifest))
    write_results(run_dir, manifest, results)
    if not grid_ok:
        raise EvidenceError('Born grid failed or incomplete')

    proc = subprocess.run([sys.executable, '-m', 'unittest', 'discover', '-s', 'tests',
                           '-p', 'test_reports.py', '-v'], cwd=ROOT,
                          text=True, capture_output=True)
    (run_dir / 'software.log').write_text(proc.stdout + proc.stderr, encoding='utf-8')
    passed = set(re.findall(r'^([a-z_]+) \(test_reports\.ReportTests\.[a-z_]+\) \.\.\. ok$',
                            proc.stdout + proc.stderr, flags=re.MULTILINE))
    for name in SOFTWARE_TESTS:
        results.append(result('software.' + name,
                              'PASS' if proc.returncode == 0 and name in passed else 'FAIL',
                              'software_test', ['Executes current calculations in a subprocess'],
                              {'test_id': name, 'exit_code': proc.returncode}, manifest))
    write_results(run_dir, manifest, results)
    if proc.returncode or passed != set(SOFTWARE_TESTS):
        raise EvidenceError('software tests failed, omitted, or changed')


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile', choices=('baseline', 'full'), default='baseline')
    parser.add_argument('--output-root', type=Path, default=ROOT / 'validation_runs')
    parser.add_argument('--born-sqrt-s', type=float, default=REFERENCE_INPUTS['sqrt_s'],
                        help='Reference-point diagnostic override; a different point cannot certify baseline.')
    args = parser.parse_args()
    inputs = dict(REFERENCE_INPUTS, sqrt_s=args.born_sqrt_s)
    run_dir, manifest = new_run(args.output_root, args.profile, inputs)
    print(f'Run directory: {run_dir}', flush=True)
    results: list[dict] = []
    ranks: dict[str, int] = {}
    try:
        execute(run_dir, manifest, inputs, results, ranks)
        if args.profile == 'full':
            finish_run(run_dir, manifest, results, 'INCOMPLETE', 'MISSING', ranks)
            print('Full validation MISSING required suites: ' + ', '.join(FULL_PENDING_SUITES),
                  file=sys.stderr)
            return 2
        finish_run(run_dir, manifest, results, 'COMPLETE', 'PASS', ranks)
        validate_run(run_dir)
    except BaseException as exc:
        try:
            finish_run(run_dir, manifest, results, 'FAILED', 'FAIL', ranks)
        except BaseException as report_exc:
            print(f'Failure report could not be written: {report_exc}', file=sys.stderr)
        print(f'Validation FAILED: {type(exc).__name__}: {exc}', file=sys.stderr)
        return 1
    print(f'Baseline PASS: {len(results)} required results; '
          f'{len(LEGACY_SYMBOLIC_LABELS)} legacy symbolic outcomes; '
          f"{next(item['result_payload']['number_of_cases'] for item in results if item['check_id'] == 'grid.complete')} "
          'Born grid cases. Full manuscript coverage remains incomplete.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
