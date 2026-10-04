"""Write a summary only from one complete, current validation run."""
from __future__ import annotations
import argparse
import os
from pathlib import Path
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from validation_evidence import EvidenceError, validate_run  # noqa: E402
from validation_manifest import BORN_CHECKS, LEGACY_IDS  # noqa: E402


def render_summary(run_dir: Path) -> str:
    manifest, document = validate_run(run_dir)
    if manifest['profile'] not in ('baseline', 'foundations', 'quark-processes'):
        raise EvidenceError('only complete implemented profiles can be summarized')
    by_id = {item['check_id']: item for item in document['results']}
    counts = {}
    roles = {}
    for item in document['results']:
        kind = item['evidence_type']
        counts[kind] = counts.get(kind, 0) + 1
        role = item['claim_role']
        roles[role] = roles.get(role, 0) + 1
    grid = by_id['grid.complete']['result_payload']
    point = by_id['born.on_shell_and_conservation']['result_payload']['report']
    lines = [
        '# Generated result summary', '',
        f"This page describes the **{manifest['profile']} profile only**. The complete manuscript "
        'validation profile has required suites that are still missing.', '',
        'The summary was generated from a complete run whose required IDs, source and input '
        'digests, computed ranks, and result identities were checked together.', '',
        f"Source revision: `{manifest['source_revision']}`; dirty source: "
        f"`{str(manifest['dirty_state']).lower()}`; scientific source digest: "
        f"`{manifest['scientific_source_digest']}`.", '',
        '| Executed evidence | Count |', '|---|---:|',
    ]
    for kind, count in sorted(counts.items()):
        lines.append(f'| {kind.replace("_", " ")} | {count} |')
    lines += ['', '| Evidence role | Records |', '|---|---:|']
    for role, count in sorted(roles.items()):
        lines.append(f'| {role.replace("_", " ")} | {count} |')
    lines += ['', '| Computed symbolic quantity | Value |', '|---|---:|']
    for kind in ('quark', 'gluon'):
        catalogue_id = LEGACY_IDS[f'{kind} catalogue has 32 entries']
        rank_id = LEGACY_IDS[f'{kind} 64-by-32 tensor map has rank 32']
        lines.append(f'| {kind.title()} catalogue entries | '
                     f"{by_id[catalogue_id]['result_payload']['catalogue_size']} |")
        lines.append(f'| {kind.title()} tensor-map rank | '
                     f"{by_id[rank_id]['result_payload']['computed_rank']} |")
    lines += ['', '## Reference Born point', '', '| Input | Value |', '|---|---:|']
    for name, value in point['inputs'].items():
        lines.append(f'| `{name}` | {value:g} |')
    lines += ['', '| Observable | Value |', '|---|---:|']
    for name, value in point['Stokes'].items():
        lines.append(f'| `{name}` | {value:.9g} |')
    for index, value in enumerate(point['B_eigenvalues'], 1):
        lines.append(f'| B eigenvalue {index} | {value:.9g} |')
    for name, value in point['relative_Ward_residuals'].items():
        lines.append(f'| {name} Ward residual | {value:.4g} |')
    lines.append(f"| Hermiticity residual | {point['relative_Hermiticity_residual']:.4g} |")
    born_count = sum('born.' + name in by_id for name in BORN_CHECKS)
    lines += ['', f"All {born_count} hard-response diagnostics passed at the reference point and "
              f"all {grid['number_of_cases']} specified grid cases.", '',
              'The grid is a set of numerical cases, not additional independent scientific '
              'claims. These coupling-stripped hard responses are not lithium-7 cross-section '
              'predictions.', '']
    if manifest['profile'] in ('foundations','quark-processes'):
        lines += ['## Independent foundations', '',
                  'These exact checks cover spin and target-state algebra, transverse STF tensors, '
                  'the lower-spin gluon dictionary, quark/gluon covariants, and coefficient recovery. '
                  'They do not cover the remaining process-response suites.', '',
                  '| Certificate | Verified value |', '|---|---:|']
        for species in ('quark', 'gluon'):
            payload = by_id[f'{species}.rank_certificate']['result_payload']
            lines.append(f"| {species.title()} exact minor rank bound | {payload['rank_lower_bound']} |")
            lines.append(f"| {species.title()} invariant Hermitian dimension | "
                         f"{by_id[f'{species}.parity_bound']['result_payload']['real_invariant_dimension']} |")
            lines.append(f"| {species.title()} recovered coefficients | "
                         f"{by_id[f'{species}.independent_linear_recovery']['result_payload']['coefficients']} |")
        lines += ['', 'The seven-direction target-response determinant is '
                  f"`{by_id['spin.seven_direction_tomography']['result_payload']['determinant']}`. "
                  'This is target-response tomography, not fourteen-TMD separation.', '']
    if manifest['profile']=='quark-processes':
        from validation_manifest import PROCESS_ROW_IDS,PROCESS_INTEGRAL_IDS,REVERSAL_IDS
        sidis=[name for name in PROCESS_ROW_IDS if name.startswith('sidis.')]
        dy=[name for name in PROCESS_ROW_IDS if name.startswith('dy.')]
        lines += ['## Quark process evidence', '',
                  f'{len(sidis)} SIDIS and {len(dy)} octupole DY row expressions were '
                  'compared exactly with separately curated fixtures. '
                  f'{len(PROCESS_INTEGRAL_IDS)} row kernels were evaluated using '
                  'Cartesian quadrature, harmonic quadrature, and exact Gaussian moments.', '',
                  f'{len(REVERSAL_IDS)} quark/gluon coefficient sign records check finite '
                  'PT algebra conditional on the stated field/link transformation. '
                  'They do not establish QCD factorization or an evolution kernel.', '',
                  'The full profile still requires gluon response separation, an independent '
                  'heavy-pair Born program, and positivity/collinear/Fourier suites.', '']
        for process,ids in (('SIDIS',sidis),('DY',dy)):
            lines += [f'### {process} reviewed response rows', '',
                      '| Row ID | Target `(K,m)` | Projection | Orbital rank | Weight | Phase `(recoil,target,lepton)` | Sign | Max integral residual |',
                      '|---|---|---|---:|---|---|---:|---:|']
            for check_id in ids:
                row=by_id[check_id]['result_payload']
                label=row['label']
                integral=by_id['integral.'+check_id]['result_payload']
                residual=max(values[-1] for values in integral['absolute_residuals'].values())
                lines.append(f'| `{check_id}` | `({label["K"]},{label["m"]})` | '
                             f'`{label["channel"]}` | {label["orbital_rank"]} | '
                             f'`{row["weight"]}` | `{tuple(row["phase_coefficients"])}` | '
                             f'{row["signed_prefactor"]:+d} | {residual:.3g} |')
            lines.append('')
    return '\n'.join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', required=True, type=Path,
                        help='Directory of one complete validation run.')
    parser.add_argument('--output', type=Path, default=ROOT / 'docs/generated-results.md')
    args = parser.parse_args()
    try:
        content = render_summary(args.run_dir)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', dir=args.output.parent,
                                         prefix='.' + args.output.name + '.', delete=False) as stream:
            temp_path = Path(stream.name)
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temp_path, args.output)
    except (EvidenceError, OSError, KeyError, TypeError, ValueError) as exc:
        print(f'Summary rejected: {exc}', file=sys.stderr)
        return 1
    print(f"Updated {args.output} from validated {args.run_dir.name} evidence")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
