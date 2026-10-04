#!/usr/bin/env python3
"""Execute cumulative run-scoped validation profiles and required leaves."""
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
                                 PROCESS_INTEGRAL_IDS, REVERSAL_IDS,
                                 GLUON_ROW_IDS,GLUON_OCT_IDS,GLUON_FIXED_IDS,
                                 GLUON_NEGATIVE_TESTS,COLLINEAR_ROW_IDS,
                                 COLLINEAR_FIXED_IDS,FOURIER_EXACT_IDS,
                                 FOURIER_NUMERIC_IDS,FOURIER_FIXED_IDS,
                                 LOCAL_RANK_IDS,LOCAL_PARITY_IDS,LOCAL_FIXED_IDS,
                                 POSITIVITY_IDS,LIMITS_NEGATIVE_TESTS)  # noqa: E402
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


def execute_gluon(run_dir: Path, manifest: dict, results: list[dict]) -> None:
    from gluon_stokes import exact_checks,numerical_checks
    from gluon_responses import verify_rows
    from gluon_octupole import verify_physical_rates,projection_checks
    from gluon_angular import certificate,check as check_angular
    from gluon_reconstruction import run_reconstruction
    from born_direct import spinor_and_ward_checks
    from born_validation import (precision_set,grid,dense_scan,broader_cases,
                                 normalization_checks,domain_checks,verify_scan_payload)
    import json
    fixed={
      'gluon.stokes.exact':exact_checks(),
      'gluon.stokes.complex':numerical_checks(),
      'gluon.angular.certificate':certificate(),
      'gluon.oct.physical_rates':{'four_rates_per_term':True,'state_construction':'rho=I/4+eta beta3 X3'},
      'gluon.oct.projections':projection_checks(),
      'gluon.oct.reconstruction':run_reconstruction(),
      'gluon.born.spinor_ward':spinor_and_ward_checks(),
      'gluon.born.precision':precision_set(),
      'gluon.born.grid':grid(),
      'gluon.born.dense_scan':dense_scan(),
      'gluon.born.broader':broader_cases(),
      'gluon.born.normalization':normalization_checks(),
      'gluon.born.domains':domain_checks(),
    }
    verify_scan_payload(fixed['gluon.born.dense_scan'])
    if set(fixed)!=set(GLUON_FIXED_IDS):raise EvidenceError('gluon fixed ID mismatch')
    angular_path=ROOT/'certificates'/'gluon_angular.json'
    check_angular(angular_path)
    atomic_json(run_dir/'gluon_certificate_attestation.json',{
        'run_id':manifest['run_id'],
        'scientific_source_digest':manifest['scientific_source_digest'],
        'certificate_sha256':file_digest(angular_path),
    })
    for check_id in GLUON_FIXED_IDS:
        results.append(result(check_id,'PASS',
            'exact_identity' if check_id.startswith(('gluon.stokes.exact','gluon.angular')) else
            'numerical_parameter_cases' if check_id.startswith('gluon.born.') else 'numerical_diagnostic',
            ['Approved source convention; physical spinor helicity h=-source lambda for matrix comparison',
             'Synthetic TMD coefficients are not lithium-7 predictions'],fixed[check_id],manifest,
            claim_role='independent_comparison'))
    rows=verify_rows()
    row_by_id={f'gluon.response.{K}{m}.{ch}.{n}':rows[f'gluon.{ch}.{K}{m}[{n}]']
               for K,m,ch,n,*_ in __import__('gluon_response_fixture').ROWS}
    if set(row_by_id)!=set(GLUON_ROW_IDS):raise EvidenceError('gluon row ID mismatch')
    for check_id in GLUON_ROW_IDS:
        payload=dict(row_by_id[check_id],fixture_sha256=file_digest(ROOT/'src'/'gluon_response_fixture.py'))
        results.append(result(check_id,'PASS','exact_fixture',
            ['Cartesian hard trace and independent complex-helicity route'],payload,manifest,
            claim_role='independent_comparison'))
    oct_rows=verify_physical_rates()
    if set('gluon.oct.'+name for name in oct_rows)!=set(GLUON_OCT_IDS):
        raise EvidenceError('octupole ID mismatch')
    for check_id in GLUON_OCT_IDS:
        name=check_id.removeprefix('gluon.oct.')
        results.append(result(check_id,'PASS','exact_fixture',
            ['Physical density matrices and four target/beam rates'],oct_rows[name],manifest,
            claim_role='independent_comparison'))
    write_results(run_dir,manifest,results)
    proc=subprocess.run([sys.executable,'-m','unittest','discover','-s','tests',
                         '-p','test_gluon_failures.py','-v'],cwd=ROOT,text=True,capture_output=True)
    (run_dir/'gluon_negative.log').write_text(proc.stdout+proc.stderr,encoding='utf-8')
    passed=set(re.findall(r'^(test_[a-z_]+) \(test_gluon_failures\.GluonFailureTests\.[a-z_]+\) \.\.\. ok$',
                          proc.stdout+proc.stderr,flags=re.MULTILINE))
    for name in GLUON_NEGATIVE_TESTS:
        results.append(result('software.gluon_negative.'+name,
            'PASS' if proc.returncode==0 and name in passed else 'FAIL','software_test',
            ['Test-local scientific mutation'],{'test_id':name,'exit_code':proc.returncode},manifest,
            claim_role='negative_control'))
    if proc.returncode or passed!=set(GLUON_NEGATIVE_TESTS):
        raise EvidenceError('gluon negative controls failed or omitted')


def check_gluon_evidence_integrity(run_dir: Path) -> dict:
    import copy,json,shutil,tempfile
    validate_run(run_dir)
    rejected=[]
    with tempfile.TemporaryDirectory(prefix='li7-gluon-integrity-') as temp:
        for mutation in ('missing-scan-case','duplicate-scan-case','wrong-helicity',
                         'missing-row','wrong-fixture','mixed-certificate','stale-precision'):
            folder=Path(temp)/mutation;shutil.copytree(run_dir,folder)
            m=json.loads((folder/'manifest.json').read_text())
            d=json.loads((folder/'results.json').read_text())
            by={x['check_id']:x for x in d['results']}
            scan=by['gluon.born.dense_scan']['result_payload']
            if mutation=='missing-scan-case':scan['cases'].pop()
            elif mutation=='duplicate-scan-case':scan['cases'][8]=copy.deepcopy(scan['cases'][7])
            elif mutation=='wrong-helicity':scan['cases'][8]['inputs']['helicity']=0.
            elif mutation=='missing-row':d['results'].remove(by[GLUON_ROW_IDS[0]])
            elif mutation=='wrong-fixture':by[GLUON_ROW_IDS[0]]['result_payload']['fixture_sha256']='0'*64
            elif mutation=='mixed-certificate':
                att=json.loads((folder/'gluon_certificate_attestation.json').read_text())
                att['run_id']='another-run';atomic_json(folder/'gluon_certificate_attestation.json',att)
            else:by['gluon.born.precision']['result_payload']['cases'][0]['dps80']['stokes'][1]='0'
            m['executed_check_ids']=[x['check_id'] for x in d['results']]
            atomic_json(folder/'results.json',d)
            m['results_digest']=file_digest(folder/'results.json')
            atomic_json(folder/'manifest.json',m)
            try:validate_run(folder)
            except EvidenceError as exc:rejected.append({'mutation':mutation,'reason':str(exc)})
            else:raise EvidenceError(f'{mutation} accepted')
    return {'rejected':rejected}


def execute_limits(run_dir: Path, manifest: dict, results: list[dict]) -> None:
    from collinear_limits import selection,operator_projection,operations_and_tail,inclusive_born_bookkeeping
    from fourier_limits import exact_gaussian,numeric_comparison
    from local_moments import check_rank_and_charge,nuclear_bookkeeping
    from conditional_positivity import (fixed_target_checks,joint_gram_checks,
        matrix_algorithm_counterexamples,collinear_blocks,collinear_witnesses,
        convention_reproducer)
    from collinear_block_certificate import check as check_blocks
    from source_joint_positivity import source_mapping_check,spectral_recovery
    by_species={name:selection(name) for name in ('quark','gluon')}
    rows={f'limits.angular.{row["id"]}':row
          for data in by_species.values() for row in data['rows']}
    if set(rows)!=set(COLLINEAR_ROW_IDS):raise EvidenceError('angular row identity mismatch')
    for check_id in COLLINEAR_ROW_IDS:
        results.append(result(check_id,'PASS','exact_identity',
            ['fixed-radius dphi/(2*pi); independent radial coefficient'],
            rows[check_id],manifest,claim_role='component_case'))
    fixed={
      'limits.collinear.selection.quark':by_species['quark'],
      'limits.collinear.selection.gluon':by_species['gluon'],
      'limits.collinear.operator_projection':operator_projection(),
      'limits.collinear.operations_uv':operations_and_tail(),
      'limits.collinear.inclusive_born':inclusive_born_bookkeeping(),
    }
    if set(fixed)!=set(COLLINEAR_FIXED_IDS):raise EvidenceError('collinear fixed ID mismatch')
    for check_id in COLLINEAR_FIXED_IDS:
        results.append(result(check_id,'PASS','conditional_algebra',
            ['Straight-link PT selection is an analytic operator input; no radial QCD matching'],
            fixed[check_id],manifest))
    exact=exact_gaussian()
    for n,check_id in enumerate(FOURIER_EXACT_IDS):
        results.append(result(check_id,'PASS','exact_identity',
            ['Lambda,M_A positive; convergent Gaussian polynomial moments'],
            {'rank':n,'component_count':1 if n==0 else 2**n,
             'source_labels':exact['source_labels']},manifest,claim_role='component_case'))
    numeric=numeric_comparison()
    mapped={f'limits.fourier.numeric.{case["kind"]}.rank_{case["rank"]}.branch_{"plus" if case["branch"]==1 else "minus"}':case
            for case in numeric['cases']}
    if set(mapped)!=set(FOURIER_NUMERIC_IDS):raise EvidenceError('Fourier numeric case identity mismatch')
    for check_id in FOURIER_NUMERIC_IDS:
        results.append(result(check_id,'PASS','numerical_parameter_cases',
            ['Convergent Gaussian or quartic radial fixture; declared finite quadrature domain'],
            mapped[check_id],manifest,claim_role='independent_comparison'))
    for check_id,payload in {
      FOURIER_FIXED_IDS[0]:{'inverse_normalization':exact['inverse_normalization'],
                             'scalar_transform':exact['scalar_transform'],
                             'inverse_cases':numeric['inverse_cases']},
      FOURIER_FIXED_IDS[1]:{'dimensions':numeric['dimensions'],
                             'max_absolute_error':numeric['max_absolute_error'],
                             'max_scaled_error':numeric['max_scaled_error'],
                             'conjugation':'F_real(b)* = F_real(-b); odd-rank phase retained'},
    }.items():
        results.append(result(check_id,'PASS','exact_identity',
            ['Convergent radial fixtures; no QCD matching claimed'],payload,manifest))
    local=check_rank_and_charge()
    rank_rows={f'limits.local.rotation.N{x["N"]}.{x["bilinear"]}':x
               for x in local['rotational_rows']}
    parity_rows={f'limits.local.antiquark.N{x["N"]}.{x["channel"]}':x
                 for x in local['parity_rows']}
    if set(rank_rows)!=set(LOCAL_RANK_IDS) or set(parity_rows)!=set(LOCAL_PARITY_IDS):
        raise EvidenceError('local moment row identity mismatch')
    for check_id in LOCAL_RANK_IDS:
        results.append(result(check_id,'PASS','exact_identity',
            ['Ambient SO(3) rotation content before Lorentz trace/mixed projection'],
            rank_rows[check_id],manifest,claim_role='component_case'))
    for check_id in LOCAL_PARITY_IDS:
        results.append(result(check_id,'PASS','conditional_algebra',
            ['Source positive-x antiquark and charge-conjugation convention'],
            parity_rows[check_id],manifest,claim_role='component_case'))
    for check_id,payload in {
      LOCAL_FIXED_IDS[0]:nuclear_bookkeeping(),
      LOCAL_FIXED_IDS[1]:{k:v for k,v in local.items() if k not in ('rotational_rows','parity_rows')},
    }.items():
        results.append(result(check_id,'PASS','conditional_algebra',
            ['Renormalized local-operator moment identification and convergence are analytic inputs'],
            payload,manifest))
    certificate_path=ROOT/'certificates'/'collinear_blocks.json'
    block_verdict=check_blocks(certificate_path)
    positivity={
      POSITIVITY_IDS[0]:fixed_target_checks(),
      POSITIVITY_IDS[1]:joint_gram_checks(),
      POSITIVITY_IDS[2]:matrix_algorithm_counterexamples(),
      POSITIVITY_IDS[3]:collinear_blocks('quark'),
      POSITIVITY_IDS[4]:collinear_blocks('gluon'),
      POSITIVITY_IDS[5]:collinear_witnesses(),
      POSITIVITY_IDS[6]:block_verdict,
      POSITIVITY_IDS[7]:spectral_recovery('quark'),
      POSITIVITY_IDS[8]:spectral_recovery('gluon'),
      POSITIVITY_IDS[9]:source_mapping_check(),
    }
    for check_id,payload in positivity.items():
        results.append(result(check_id,'PASS','exact_identity',
            ['Positive spectral/input prescription; source-indexed operator definition'],
            payload,manifest,claim_role='independent_comparison'))
    atomic_json(run_dir/'limits_certificate_attestation.json',{
        'run_id':manifest['run_id'],'scientific_source_digest':manifest['scientific_source_digest'],
        'certificate_sha256':file_digest(certificate_path)})
    discrepancy=convention_reproducer()
    atomic_json(run_dir/'source_joint_convention_mapping.json',discrepancy)
    write_results(run_dir,manifest,results)
    proc=subprocess.run([sys.executable,'-m','unittest','discover','-s','tests',
                         '-p','test_limits_failures.py','-v'],cwd=ROOT,text=True,capture_output=True)
    (run_dir/'limits_negative.log').write_text(proc.stdout+proc.stderr,encoding='utf-8')
    passed=set(re.findall(r'^(test_[a-z_]+) \(test_limits_failures\.LimitsFailureTests\.[a-z_]+\) \.\.\. ok$',
                          proc.stdout+proc.stderr,flags=re.MULTILINE))
    for name in LIMITS_NEGATIVE_TESTS:
        results.append(result('software.limits_negative.'+name,
            'PASS' if proc.returncode==0 and name in passed else 'FAIL','software_test',
            ['Test-local scientific mutation'],{'test_id':name,'exit_code':proc.returncode},
            manifest,claim_role='negative_control'))
    if proc.returncode or passed!=set(LIMITS_NEGATIVE_TESTS):
        raise EvidenceError('limits negative controls failed or omitted')


def check_limits_evidence_integrity(run_dir: Path) -> dict:
    import copy,json,shutil,tempfile
    validate_run(run_dir)
    rejected=[]
    with tempfile.TemporaryDirectory(prefix='li7-limits-integrity-') as temp:
        for mutation in ('missing-rank','wrong-fourier-domain','wrong-local-sign',
                         'mixed-certificate','mixed-run','analytic-promotion'):
            folder=Path(temp)/mutation;shutil.copytree(run_dir,folder)
            m=json.loads((folder/'manifest.json').read_text())
            d=json.loads((folder/'results.json').read_text())
            by={x['check_id']:x for x in d['results']}
            if mutation=='missing-rank':
                d['results'].remove(by[FOURIER_EXACT_IDS[-1]])
            elif mutation=='wrong-fourier-domain':
                by[FOURIER_NUMERIC_IDS[-1]]['result_payload']['M_A']=1.0
            elif mutation=='wrong-local-sign':
                by[LOCAL_PARITY_IDS[0]]['result_payload']['antiquark_sign']*=-1
            elif mutation=='mixed-certificate':
                att=json.loads((folder/'limits_certificate_attestation.json').read_text())
                att['run_id']='another-run';atomic_json(folder/'limits_certificate_attestation.json',att)
            elif mutation=='mixed-run':
                by[COLLINEAR_ROW_IDS[0]]['run_id']='another-run'
            else:
                by['limits.positivity.source_joint_convention']['evidence_type']='analytic_only'
            m['executed_check_ids']=[x['check_id'] for x in d['results']]
            atomic_json(folder/'results.json',d)
            m['results_digest']=file_digest(folder/'results.json')
            atomic_json(folder/'manifest.json',m)
            try:validate_run(folder)
            except EvidenceError as exc:rejected.append({'mutation':mutation,'reason':str(exc)})
            else:raise EvidenceError(f'limits evidence mutation accepted: {mutation}')
    return {'rejected':rejected,'scope':'tampered complete-run copies; source and certificate checks'}

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile', choices=('baseline', 'foundations', 'quark-processes',
                                               'gluon-processes','limits-positivity','full'), default='baseline')
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
        if args.profile in ('foundations','quark-processes','gluon-processes','limits-positivity','full'):
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
        if args.profile in ('quark-processes','gluon-processes','limits-positivity','full'):
            execute_processes(run_dir,manifest,results)
            results.append(result('software.process_evidence_integrity','PASS','software_test',
                                  ['Run-specific complete process evidence'],
                                  {'stage':'pending independent mutation'},manifest,
                                  claim_role='negative_control'))
            if args.profile in ('gluon-processes','limits-positivity','full'):
                manifest['profile']='quark-processes'
                manifest['required_check_ids']=list(PROFILE_REQUIRED['quark-processes'])
            finish_run(run_dir,manifest,results,'COMPLETE','PASS',ranks)
            try:
                payload=check_process_evidence_integrity(run_dir)
            finally:
                if args.profile in ('gluon-processes','limits-positivity','full'):
                    manifest['profile']=args.profile
                    manifest['required_check_ids']=list(PROFILE_REQUIRED[args.profile])
            results[-1]=result('software.process_evidence_integrity','PASS','software_test',
                               ['Tampered copies of this completed run'],payload,manifest,
                               claim_role='negative_control')
        if args.profile in ('gluon-processes','limits-positivity','full'):
            execute_gluon(run_dir,manifest,results)
            results.append(result('software.gluon_evidence_integrity','PASS','software_test',
                                  ['Run-specific complete gluon evidence'],
                                  {'stage':'pending independent mutation'},manifest,
                                  claim_role='negative_control'))
            if args.profile in ('limits-positivity','full'):
                manifest['profile']='gluon-processes'
                manifest['required_check_ids']=list(PROFILE_REQUIRED['gluon-processes'])
            finish_run(run_dir,manifest,results,'COMPLETE','PASS',ranks)
            try:payload=check_gluon_evidence_integrity(run_dir)
            finally:
                if args.profile in ('limits-positivity','full'):
                    manifest['profile']=args.profile
                    manifest['required_check_ids']=list(PROFILE_REQUIRED[args.profile])
            results[-1]=result('software.gluon_evidence_integrity','PASS','software_test',
                               ['Tampered copies of this completed run'],payload,manifest,
                               claim_role='negative_control')
        if args.profile in ('limits-positivity','full'):
            execute_limits(run_dir,manifest,results)
            if FULL_PENDING_SUITES:
                finish_run(run_dir, manifest, results, 'INCOMPLETE', 'MISSING', ranks)
                print(args.profile+' MISSING required scientific ID: ' + ', '.join(FULL_PENDING_SUITES),
                      file=sys.stderr)
                return 2
            results.append(result('software.limits_evidence_integrity','PASS','software_test',
                ['Run-specific complete evidence'],{'stage':'pending independent mutation'},
                manifest,claim_role='negative_control'))
            finish_run(run_dir,manifest,results,'COMPLETE','PASS',ranks)
            payload=check_limits_evidence_integrity(run_dir)
            results[-1]=result('software.limits_evidence_integrity','PASS','software_test',
                ['Tampered copies of this completed run'],payload,manifest,
                claim_role='negative_control')
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
          'Born grid cases. Analytical QCD assumptions remain outside computational evidence.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
