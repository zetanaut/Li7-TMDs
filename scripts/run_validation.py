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
                                 LEGACY_SYMBOLIC_LABELS, SOFTWARE_TESTS,
                                 FOUNDATION_REQUIRED, BASELINE_REQUIRED,
                                 FOUNDATION_NEGATIVE_TESTS, PROFILE_REQUIRED,
                                 PROCESS_NEGATIVE_TESTS, PROCESS_ROW_IDS,
                                 PROCESS_INTEGRAL_IDS, REVERSAL_IDS)  # noqa: E402
from validation_evidence import (GRID_AXES, REFERENCE_INPUTS, EvidenceError, finish_run,
                                 new_run, result, validate_run, atomic_json, file_digest,
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


def execute_foundations(run_dir: Path, manifest: dict, results: list[dict]) -> None:
    from spin_foundations import run_checks as spin_checks
    from spin_state_foundations import run_checks as state_checks
    from transverse_foundations import run_checks as stf_checks
    from correlator_foundations import (run_checks as covariant_checks, parity_bound,
                                        rotation_covariance, joint_expectation_check)
    from gluon_dictionary_foundations import run_checks as dictionary_checks
    from projector_foundations import run_checks as projector_checks
    from foundation_certificates import read_certificate, verify_certificate
    found = {}
    for compute in (spin_checks, state_checks, stf_checks, covariant_checks,
                    dictionary_checks, projector_checks):
        suite = compute()
        if set(found) & set(suite):
            raise EvidenceError('duplicate foundation check ID')
        found.update(suite)
    for species in ('quark','gluon'):
        found[f'{species}.parity_bound'] = parity_bound(species)
        found[f'{species}.rotation_covariance'] = rotation_covariance(species)
        found[f'{species}.joint_expectation'] = joint_expectation_check(species)
        certificate = read_certificate(ROOT / 'certificates' / f'{species}_rank.json')
        found[f'{species}.rank_certificate'] = verify_certificate(certificate)
    expected = set(FOUNDATION_REQUIRED) - set(BASELINE_REQUIRED) - {
        'software.foundation_negative.'+name for name in FOUNDATION_NEGATIVE_TESTS} - {
        'software.foundation_evidence_integrity'}
    if set(found) != expected:
        raise EvidenceError('foundation scientific ID mismatch: missing '+str(sorted(expected-set(found)))+
                            '; unexpected '+str(sorted(set(found)-expected)))
    for check_id in (item for item in FOUNDATION_REQUIRED if item in found):
        payload=found[check_id]
        role = ('rank_witness' if check_id.endswith('rank_certificate') else
                'independent_bound' if check_id.endswith('parity_bound') else
                'component_case' if check_id.startswith(('spin.direction_unit_', 'stf.rank_')) else
                'independent_comparison' if check_id.endswith(('comparison','symbolic_covariants',
                                      'independent_linear_recovery','dictionary_complete')) else
                'unique_claim')
        results.append(result(check_id,'PASS',
                              'exact_rank' if check_id.endswith(('rank_certificate','parity_bound')) else 'exact_identity',
                              ['Exact finite-dimensional algebra; M_A>0 where divided',
                               'Approved manuscript equations encoded in source; no runtime manuscript'],
                              payload,manifest,claim_role=role))
    write_results(run_dir,manifest,results)
    proc=subprocess.run([sys.executable,'-m','unittest','discover','-s','tests',
                         '-p','test_foundation_failures.py','-v'],cwd=ROOT,text=True,capture_output=True)
    (run_dir/'foundation_negative.log').write_text(proc.stdout+proc.stderr,encoding='utf-8')
    passed=set(re.findall(r'^(test_[a-z_]+) \(test_foundation_failures\.FoundationFailureTests\.[a-z_]+\) \.\.\. ok$',
                          proc.stdout+proc.stderr,flags=re.MULTILINE))
    for name in FOUNDATION_NEGATIVE_TESTS:
        results.append(result('software.foundation_negative.'+name,
                              'PASS' if proc.returncode==0 and name in passed else 'FAIL',
                              'software_test',['Injected local mutation reaches named scientific diagnostic'],
                              {'test_id':name,'exit_code':proc.returncode},manifest,
                              claim_role='negative_control'))
    if proc.returncode or passed!=set(FOUNDATION_NEGATIVE_TESTS):
        raise EvidenceError('foundation negative controls failed, missing or changed')
    atomic_json(run_dir/'certificate_attestation.json',{
        'run_id':manifest['run_id'],
        'scientific_source_digest':manifest['scientific_source_digest'],
        'certificates':{species:found[f'{species}.rank_certificate']['certificate_sha256']
                        for species in ('quark','gluon')},
    })


def check_foundation_evidence_integrity(run_dir: Path) -> dict:
    import copy
    import json
    import shutil
    import tempfile
    from validation_evidence import atomic_json
    validate_run(run_dir)
    with tempfile.TemporaryDirectory(prefix='li7-foundation-integrity-') as temp:
        root=Path(temp)
        missing=root/'missing';shutil.copytree(run_dir,missing)
        manifest=json.loads((missing/'manifest.json').read_text())
        document=json.loads((missing/'results.json').read_text())
        document['results']=[row for row in document['results'] if row['check_id']!='spin.seven_direction_tomography']
        manifest['executed_check_ids']=[row['check_id'] for row in document['results']]
        atomic_json(missing/'results.json',document)
        manifest['results_digest']=file_digest(missing/'results.json')
        atomic_json(missing/'manifest.json',manifest)
        try:validate_run(missing)
        except EvidenceError as exc:
            if 'missing or unexpected required check ID' not in str(exc):raise
        else:raise EvidenceError('missing foundation check was accepted')
        for mutation in ('missing-grid-point','duplicate-grid-point','failed-grid-diagnostic'):
            folder=root/mutation;shutil.copytree(run_dir,folder)
            changed_manifest=json.loads((folder/'manifest.json').read_text())
            changed_document=json.loads((folder/'results.json').read_text())
            grid=next(row for row in changed_document['results'] if row['check_id']=='grid.complete')['result_payload']
            cases=grid['cases']
            if mutation=='missing-grid-point':cases.pop()
            elif mutation=='duplicate-grid-point':cases[-1]=copy.deepcopy(cases[0])
            else:cases[-1]['checks']['photon_Ward']=False
            atomic_json(folder/'results.json',changed_document)
            changed_manifest['results_digest']=file_digest(folder/'results.json')
            atomic_json(folder/'manifest.json',changed_manifest)
            try:validate_run(folder)
            except EvidenceError as exc:
                if 'incomplete or failed grid cases' not in str(exc):raise
            else:raise EvidenceError(mutation+' was accepted')
        mixed=root/'mixed';shutil.copytree(run_dir,mixed)
        attestation=json.loads((mixed/'certificate_attestation.json').read_text())
        attestation['run_id']='another-run'
        atomic_json(mixed/'certificate_attestation.json',attestation)
        try:validate_run(mixed)
        except EvidenceError as exc:
            if 'certificate attestation run identity mismatch' not in str(exc):raise
        else:raise EvidenceError('mixed-run certificates were accepted')
    return {'negative_cases':['omitted tomography ID','different certificate run ID','missing grid point',
                              'duplicate grid point','failed grid diagnostic'],
            'method':'tampered copies of the completed run rejected by summary validator'}


def execute_processes(run_dir: Path, manifest: dict, results: list[dict]) -> None:
    from process_dirac import check_algebra
    from process_normalization import sidis_normalization,dy_normalization,qed_current_checks
    from process_observables import target_checks,flavor_checks,sidis_flavor_checks,spin_difference_checks
    from process_convolutions import check_row_integral,momentum_sign_control
    from process_responses import compare_row
    from process_reversal import check_reversal,conditional_evolution_check
    from response_fixtures import SIDIS_ROWS,DY_ROWS,row_id
    algebra=check_algebra()
    fixed={
        'process.dirac.algebra':{key:algebra[key] for key in ('dirac','chiral','basis_conversion')},
        'process.dirac.sidis_trace':algebra['SIDIS'],
        'process.dirac.dy_trace':algebra['DY'],
        'process.sidis.normalization':sidis_normalization(),
        'process.sidis.flavors':sidis_flavor_checks(),
        'process.dy.current':qed_current_checks(),
        'process.dy.normalization':dy_normalization(),
        'process.dy.flavors':flavor_checks(),
        'process.target.preparations':target_checks(),
        'process.spin_differences':spin_difference_checks(),
        'process.convolution.momentum_sign':momentum_sign_control(),
    }
    reversal=check_reversal()
    fixed['process.reversal.density_links']={key:value for key,value in reversal.items() if key!='tables'}
    fixed['process.reversal.conditional_evolution']=conditional_evolution_check()
    for check_id,payload in fixed.items():
        evidence_type=('conditional_algebra' if check_id.startswith('process.reversal.') else
                       'numerical_diagnostic' if check_id=='process.dy.current' else 'exact_identity')
        results.append(result(check_id,'PASS',evidence_type,
                              ['Four-dimensional leading-power Born conventions',
                               'Operator reversal conditional on eq:PTprojection where applicable'],
                              payload,manifest,claim_role='unique_claim'))
    fixture_sha=file_digest(ROOT/'src'/'response_fixtures.py')
    for process,rows in (('SIDIS',SIDIS_ROWS),('DY',DY_ROWS)):
        for row in rows:
            check_id=row_id(process,row)
            payload=compare_row(process,row)
            payload['fixture_sha256']=fixture_sha
            results.append(result(check_id,'PASS','exact_fixture',
                                  ['eq:SFrow target amplitude and beam factor outside convolution',
                                   'Reflection of radial TMDs about recoil axis'],payload,manifest,
                                  claim_role='independent_comparison'))
            integral=check_row_integral(process,row)
            integral['source_row_id']=check_id
            results.append(result('integral.'+check_id,'PASS','numerical_diagnostic',
                                  ['Synthetic radial Gaussian functions; not nuclear predictions',
                                   'Fourier angular factors evaluated at generic fixed phases'],
                                  integral,manifest,claim_role='numerical_corroboration'))
    for species in ('quark','gluon'):
        for item in reversal['tables'][species]['rows']:
            check_id=(f'process.reversal.{species}.{item["channel"]}.'
                      f'{item["K"]}{item["m"]}.{item["orbital_rank"]}')
            results.append(result(check_id,'PASS','conditional_algebra',
                                  ['eq:PTprojection field/link transformation is analytic input',
                                   'Momentum labels held fixed under combined PT'],
                                  item,manifest,claim_role='component_case'))
    expected=set(PROCESS_ROW_IDS)|set(PROCESS_INTEGRAL_IDS)|set(REVERSAL_IDS)
    present={item['check_id'] for item in results}
    if not expected<=present:
        raise EvidenceError('missing process leaf IDs '+str(sorted(expected-present)))
    proc=subprocess.run([sys.executable,'-m','unittest','discover','-s','tests',
                         '-p','test_process_failures.py','-v'],cwd=ROOT,text=True,capture_output=True)
    (run_dir/'process_negative.log').write_text(proc.stdout+proc.stderr,encoding='utf-8')
    passed=set(re.findall(r'^(test_[a-z_]+) \(test_process_failures\.ProcessFailureTests\.[a-z_]+\) \.\.\. ok$',
                          proc.stdout+proc.stderr,flags=re.MULTILINE))
    for name in PROCESS_NEGATIVE_TESTS:
        results.append(result('software.process_negative.'+name,
                              'PASS' if proc.returncode==0 and name in passed else 'FAIL',
                              'software_test',['Test-local convention mutation'],
                              {'test_id':name,'exit_code':proc.returncode},manifest,
                              claim_role='negative_control'))
    if proc.returncode or passed!=set(PROCESS_NEGATIVE_TESTS):
        raise EvidenceError('process negative controls failed, missing or changed')


def check_process_evidence_integrity(run_dir: Path) -> dict:
    import copy,json,shutil,tempfile
    validate_run(run_dir)
    cases=[]
    with tempfile.TemporaryDirectory(prefix='li7-process-integrity-') as temp:
        for mutation in ('missing-row','duplicate-row','replaced-row','fixture-digest','mixed-run','analytic-promotion'):
            folder=Path(temp)/mutation;shutil.copytree(run_dir,folder)
            m=json.loads((folder/'manifest.json').read_text())
            d=json.loads((folder/'results.json').read_text())
            row=next(item for item in d['results'] if item['check_id']==PROCESS_ROW_IDS[0])
            if mutation=='missing-row':d['results'].remove(row)
            elif mutation=='duplicate-row':d['results'].append(copy.deepcopy(row))
            elif mutation=='replaced-row':
                other=next(item for item in d['results'] if item['check_id']==PROCESS_ROW_IDS[1])
                row['result_payload']=copy.deepcopy(other['result_payload'])
            elif mutation=='fixture-digest':row['result_payload']['fixture_sha256']='0'*64
            elif mutation=='mixed-run':row['run_id']='another-run'
            else:row['evidence_type']='analytic_only'
            m['executed_check_ids']=[item['check_id'] for item in d['results']]
            atomic_json(folder/'results.json',d)
            m['results_digest']=file_digest(folder/'results.json')
            atomic_json(folder/'manifest.json',m)
            try:validate_run(folder)
            except EvidenceError as exc:
                cases.append({'mutation':mutation,'rejected_by':str(exc)})
            else:raise EvidenceError(mutation+' was accepted')
    return {'negative_cases':cases,'method':'tampered complete run copies rejected'}

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile', choices=('baseline', 'foundations', 'quark-processes', 'full'), default='baseline')
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
        if args.profile in ('foundations','quark-processes','full'):
            execute_foundations(run_dir,manifest,results)
            results.append(result('software.foundation_evidence_integrity','PASS','software_test',
                                  ['Run-specific complete evidence'],{'stage':'pending independent mutation'},manifest,
                                  claim_role='negative_control'))
            if args.profile != 'foundations':
                manifest['profile']='foundations'
                manifest['required_check_ids']=list(PROFILE_REQUIRED['foundations'])
            finish_run(run_dir,manifest,results,'COMPLETE','PASS',ranks)
            try:
                payload=check_foundation_evidence_integrity(run_dir)
            finally:
                if args.profile != 'foundations':
                    manifest['profile']=args.profile
                    manifest['required_check_ids']=list(PROFILE_REQUIRED[args.profile])
            results[-1]=result('software.foundation_evidence_integrity','PASS','software_test',
                               ['Tampered copies of this completed run'],payload,manifest,
                               claim_role='negative_control')
        if args.profile in ('quark-processes','full'):
            execute_processes(run_dir,manifest,results)
            results.append(result('software.process_evidence_integrity','PASS','software_test',
                                  ['Run-specific complete process evidence'],
                                  {'stage':'pending independent mutation'},manifest,
                                  claim_role='negative_control'))
            if args.profile=='full':
                manifest['profile']='quark-processes'
                manifest['required_check_ids']=list(PROFILE_REQUIRED['quark-processes'])
            finish_run(run_dir,manifest,results,'COMPLETE','PASS',ranks)
            try:
                payload=check_process_evidence_integrity(run_dir)
            finally:
                if args.profile=='full':
                    manifest['profile']='full'
                    manifest['required_check_ids']=list(PROFILE_REQUIRED['full'])
            results[-1]=result('software.process_evidence_integrity','PASS','software_test',
                               ['Tampered copies of this completed run'],payload,manifest,
                               claim_role='negative_control')
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
    print(f'{args.profile.capitalize()} PASS: {len(results)} required results; '
          f'{len(LEGACY_SYMBOLIC_LABELS)} legacy symbolic outcomes; '
          f"{next(item['result_payload']['number_of_cases'] for item in results if item['check_id'] == 'grid.complete')} "
          'Born grid cases. Full manuscript coverage remains incomplete.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
