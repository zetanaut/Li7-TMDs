"""Run-scoped evidence and strict validation of the implemented baseline."""
from __future__ import annotations

import hashlib
import importlib.metadata
import itertools
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import uuid

from validation_manifest import (BORN_CHECKS, LEGACY_IDS, PROFILE_REQUIRED,
                                 SOFTWARE_TESTS, PROCESS_ROW_IDS,
                                 PROCESS_INTEGRAL_IDS, REVERSAL_IDS,
                                 GLUON_ROW_IDS,GLUON_OCT_IDS,GLUON_NEGATIVE_TESTS,
                                 COLLINEAR_ROW_IDS,FOURIER_EXACT_IDS,FOURIER_NUMERIC_IDS,
                                 LOCAL_RANK_IDS,LOCAL_PARITY_IDS,LIMITS_NEGATIVE_TESTS,
                                 CONVENTION_IDS,CONVENTION_NEGATIVE_TESTS)

SCHEMA_VERSION = 2
LEGACY_SOURCE_SHA256 = '3b5aaff51a77932ad561c1137a6d1bb5f0e4c60353add1d6f9034f2d7e2b892d'
REFERENCE_SOURCE_SHA256 = '7436c7b7dab536999f7caf42af54c5820fa474ca81ab4ff8e4cf15259c9a24ee'
CONVENTION_ID = 'born-current-v2-physical-h'
LEGACY_CONVENTION_ID = 'born-current-v1-legacy-source-label'
CONVENTION_REVIEW = {
    'computational_execution_status':'COMPLETE_CORRECTED_COMPARISONS',
    'source_formula_agreement_status':'CORRECTED_PHYSICAL_HELICITY_AGREES',
    'physical_helicity_anchor_status':'VERIFIED_EIGENVALUE',
    'current_index_contract_status':'APPROVED_CORRECTION_EXECUTED',
    'source_gram_index_status':'VERIFIED_CONDITIONAL_SPECTRAL_ORDER',
    'unresolved_scientific_issues':[],
    'publication_eligibility':'READY_FOR_PUBLICATION_REVIEW',
    'publication_authorized':False,
    'correction_id':'born-current-index-c7-author-approved',
    'legacy_source_sha256':LEGACY_SOURCE_SHA256,
    'corrected_source_sha256':REFERENCE_SOURCE_SHA256,
    'convention_id':CONVENTION_ID,
}
ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOTS = ('src', 'examples', 'scripts', 'tests')
REFERENCE_INPUTS = {
    'sqrt_s': 5.0, 'Q2': 4.0, 'mass': 1.5, 'theta': 0.8,
    'phi': 0.4, 'lepton_energy': 10.0, 'helicity': 1.0,
}
GRID_AXES = {
    'theta': [0.4, 0.8, 1.7, 2.7],
    'phi': [0.0, 0.4, 1.2],
    'helicity': [-1.0, 0.0, 1.0],
}
EVIDENCE_TYPES = {'exact_identity', 'exact_rank', 'numerical_diagnostic',
                  'numerical_parameter_cases', 'software_test',
                  'exact_fixture', 'conditional_algebra', 'analytic_only'}
CLAIM_ROLES = {'legacy_identity','legacy_rank','unique_claim','component_case',
               'independent_comparison','rank_witness','independent_bound',
               'numerical_corroboration','negative_control','software_regression'}


class EvidenceError(ValueError):
    """A run cannot support the requested claim of completion."""


def assert_convention_review(manifest: dict, current_payload: dict | None = None) -> None:
    """Require executed corrected agreement and explicit legacy mismatch."""
    if manifest.get('convention_review') != CONVENTION_REVIEW:
        raise EvidenceError('convention review altered')
    if current_payload is not None and current_payload.get('status') != {
        'diagnostic_execution':'PASS',
        'legacy_literal_physical_agreement':'EXPECTED_MISMATCH',
        'corrected_physical_source_agreement':'AGREES',
        'publication_eligibility':'READY_FOR_PUBLICATION_REVIEW',
    }:
        raise EvidenceError('corrected physical comparison absent')
    if current_payload is not None and (not current_payload.get('born_cases') or any(
            case.get('corrected_max_abs',float('inf'))>1e-8 or
            case.get('migration_max_abs',float('inf'))>1e-8 or
            case.get('literal_max_abs',0)<1
            for case in current_payload['born_cases'])):
        raise EvidenceError('corrected physical comparison failed')


def canonical_bytes(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(',', ':'),
                      ensure_ascii=True, allow_nan=False).encode('utf-8')


def digest(value: object) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def file_digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def scientific_source_digest() -> str:
    files = [path for folder in SOURCE_ROOTS for path in (ROOT / folder).rglob('*.py')]
    files += [ROOT / 'requirements.txt', ROOT / 'requirements-repro-py311.txt',
              ROOT / 'mkdocs.yml']
    files += list((ROOT / '.github/workflows').glob('*.yml'))
    return digest({str(path.relative_to(ROOT)): file_digest(path) for path in sorted(files)})


def input_spec(reference_inputs: dict | None = None) -> dict:
    return {'reference': reference_inputs or REFERENCE_INPUTS, 'grid_axes': GRID_AXES,
            'convention_id':CONVENTION_ID,'legacy_convention_id':LEGACY_CONVENTION_ID,
            'legacy_source_sha256':LEGACY_SOURCE_SHA256,
            'symbolic': 'exact finite algebra; no runtime TeX input'}


def input_digest(reference_inputs: dict | None = None) -> str:
    return digest(input_spec(reference_inputs))


def git_provenance() -> tuple[str, bool]:
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    dirty = bool(subprocess.check_output(
        ['git', 'status', '--porcelain', '--untracked-files=all'], cwd=ROOT, text=True).strip())
    return revision, dirty


def dependency_versions() -> dict[str, str]:
    packages = ('numpy', 'sympy', 'matplotlib', 'mkdocs', 'pymdown-extensions','mpmath')
    return {'python': sys.version.split()[0],
            **{name: importlib.metadata.version(name) for name in packages}}


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    data = json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n'
    with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', dir=path.parent,
                                     prefix='.' + path.name + '.', delete=False) as stream:
        temp_path = Path(stream.name)
        try:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        except BaseException:
            temp_path.unlink(missing_ok=True)
            raise
    os.replace(temp_path, path)


def new_run(output_root: Path, profile: str,
            reference_inputs: dict | None = None) -> tuple[Path, dict]:
    if profile not in PROFILE_REQUIRED:
        raise EvidenceError(f'unknown profile: {profile}')
    run_id = uuid.uuid4().hex
    run_dir = output_root / ('run-' + run_id)
    run_dir.mkdir(parents=True, exist_ok=False)
    revision, dirty = git_provenance()
    manifest = {
        'schema_version': SCHEMA_VERSION, 'run_id': run_id, 'profile': profile,
        'run_state': 'RUNNING', 'status': 'INCOMPLETE',
        'required_check_ids': list(PROFILE_REQUIRED[profile]),
        'executed_check_ids': [], 'missing_check_ids': [],
        'source_revision': revision, 'dirty_state': dirty,
        'reference_sha256': REFERENCE_SOURCE_SHA256,
        'scientific_source_digest': scientific_source_digest(),
        'input_spec': input_spec(reference_inputs),
        'input_digest': input_digest(reference_inputs),
        'dependency_versions': dependency_versions(),
        'results_digest': None, 'rank_metadata': {},
    }
    atomic_json(run_dir / 'manifest.json', manifest)
    return run_dir, manifest


def result(check_id: str, status: str, evidence_type: str, assumptions: list[str],
           payload: object, manifest: dict, *, claim_role: str | None = None) -> dict:
    if claim_role is None:
        claim_role = {'exact_identity':'legacy_identity','exact_rank':'legacy_rank',
                      'numerical_diagnostic':'numerical_corroboration',
                      'numerical_parameter_cases':'component_case',
                      'software_test':'software_regression',
                      'exact_fixture':'independent_comparison',
                      'conditional_algebra':'unique_claim',
                      'analytic_only':'unique_claim'}[evidence_type]
    if claim_role not in CLAIM_ROLES:
        raise EvidenceError('unknown claim role')
    return {
        'claim_role': claim_role, 'check_id': check_id, 'run_id': manifest['run_id'],
        'scientific_source_digest': manifest['scientific_source_digest'],
        'input_digest': manifest['input_digest'], 'status': status,
        'evidence_type': evidence_type, 'assumptions': assumptions,
        'result_payload': payload,
    }


def write_results(run_dir: Path, manifest: dict, results: list[dict]) -> None:
    atomic_json(run_dir / 'results.json', {
        'schema_version': SCHEMA_VERSION, 'run_id': manifest['run_id'],
        'scientific_source_digest': manifest['scientific_source_digest'],
        'input_digest': manifest['input_digest'], 'results': results,
    })


def finish_run(run_dir: Path, manifest: dict, results: list[dict],
               state: str, status: str, ranks: dict[str, int] | None = None) -> None:
    write_results(run_dir, manifest, results)
    manifest['executed_check_ids'] = [item['check_id'] for item in results]
    manifest['missing_check_ids'] = sorted(set(manifest['required_check_ids']) -
                                           set(manifest['executed_check_ids']))
    manifest['rank_metadata'] = ranks or {}
    manifest['results_digest'] = file_digest(run_dir / 'results.json')
    manifest['run_state'] = state
    manifest['status'] = status
    atomic_json(run_dir / 'manifest.json', manifest)


def _unique_pairs(pairs: list[tuple[str, object]]) -> dict:
    output = {}
    for key, value in pairs:
        if key in output:
            raise EvidenceError(f'duplicate JSON key: {key}')
        output[key] = value
    return output


def _read_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(encoding='utf-8'),
                           object_pairs_hook=_unique_pairs,
                           parse_constant=lambda value: (_ for _ in ()).throw(
                               EvidenceError(f'non-finite JSON constant: {value}')))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise EvidenceError(f'missing or corrupt report: {path.name}: {exc}') from exc
    if not isinstance(value, dict):
        raise EvidenceError(f'malformed report: {path.name}')
    return value


def _finite(value: object) -> bool:
    if isinstance(value, float):
        return math.isfinite(value)
    if isinstance(value, dict):
        return all(isinstance(k, str) and _finite(v) for k, v in value.items())
    if isinstance(value, list):
        return all(_finite(v) for v in value)
    return value is None or isinstance(value, (str, int, bool))


def validate_run(run_dir: Path, *, require_current: bool = True) -> tuple[dict, dict]:
    manifest = _read_json(run_dir / 'manifest.json')
    profile = manifest.get('profile')
    if manifest.get('schema_version') != SCHEMA_VERSION or profile not in PROFILE_REQUIRED:
        raise EvidenceError('unknown schema or profile')
    expected = PROFILE_REQUIRED[profile]
    if not expected or manifest.get('required_check_ids') != list(expected):
        raise EvidenceError('required-check manifest mismatch or empty collection')
    if manifest.get('run_state') != 'COMPLETE' or manifest.get('status') != 'PASS':
        raise EvidenceError(f"run is {manifest.get('run_state')}/{manifest.get('status')}")
    if (not isinstance(manifest.get('run_id'), str) or not manifest['run_id'] or
            not isinstance(manifest.get('dirty_state'), bool) or
            not isinstance(manifest.get('dependency_versions'), dict)):
        raise EvidenceError('malformed run provenance')
    if manifest.get('reference_sha256') != REFERENCE_SOURCE_SHA256:
        raise EvidenceError('reference version mismatch')
    if profile == 'full':assert_convention_review(manifest)
    if manifest.get('missing_check_ids'):
        raise EvidenceError('missing required suite or check')
    if manifest.get('input_spec') != input_spec() or manifest.get('input_digest') != digest(manifest['input_spec']):
        raise EvidenceError('input digest or baseline input mismatch')
    if require_current:
        revision, _dirty = git_provenance()
        if manifest.get('source_revision') != revision:
            raise EvidenceError('source revision mismatch')
        if manifest.get('scientific_source_digest') != scientific_source_digest():
            raise EvidenceError('scientific source digest mismatch')
        if manifest['dependency_versions'] != dependency_versions():
            raise EvidenceError('dependency version mismatch')
    results_path = run_dir / 'results.json'
    if not results_path.is_file():
        raise EvidenceError('missing results file')
    if manifest.get('results_digest') != file_digest(results_path):
        raise EvidenceError('results digest mismatch')
    document = _read_json(results_path)
    for field in ('run_id', 'scientific_source_digest', 'input_digest'):
        if document.get(field) != manifest.get(field):
            raise EvidenceError(f'results {field} mismatch')
    if document.get('schema_version') != SCHEMA_VERSION:
        raise EvidenceError('results schema mismatch')
    results = document.get('results')
    if not isinstance(results, list) or not results:
        raise EvidenceError('empty or malformed result collection')
    if not all(isinstance(item, dict) and isinstance(item.get('check_id'), str)
               for item in results):
        raise EvidenceError('malformed check result or ID')
    ids = [item.get('check_id') if isinstance(item, dict) else None for item in results]
    if len(ids) != len(set(ids)):
        raise EvidenceError('duplicate check ID')
    executed = manifest.get('executed_check_ids')
    if not isinstance(executed, list) or not all(isinstance(value, str) for value in executed):
        raise EvidenceError('malformed executed-check IDs')
    if set(ids) != set(expected) or set(executed) != set(expected):
        raise EvidenceError('missing or unexpected required check ID')
    if manifest.get('executed_check_ids') != ids:
        raise EvidenceError('executed-check order mismatch')
    by_id = {}
    for item in results:
        for field in ('run_id', 'scientific_source_digest', 'input_digest'):
            if item.get(field) != manifest.get(field):
                raise EvidenceError(f'check {field} mismatch')
        if item.get('status') != 'PASS':
            raise EvidenceError(f"failed or inconclusive check: {item.get('check_id')}")
        if item.get('evidence_type') not in EVIDENCE_TYPES:
            raise EvidenceError('unknown evidence type')
        if item.get('evidence_type') == 'analytic_only':
            raise EvidenceError('analytic-only assertion cannot certify executed check')
        if item.get('claim_role') not in CLAIM_ROLES:
            raise EvidenceError('unknown claim role')
        if (not isinstance(item.get('assumptions'), list) or
                not isinstance(item.get('result_payload'), dict) or
                not _finite(item['result_payload'])):
            raise EvidenceError('malformed or non-finite result payload')
        by_id[item['check_id']] = item
    if not isinstance(manifest.get('rank_metadata'), dict):
        raise EvidenceError('malformed rank metadata')
    for kind in ('quark', 'gluon'):
        check_id = LEGACY_IDS[f'{kind} 64-by-32 tensor map has rank 32']
        actual = by_id[check_id]['result_payload'].get('computed_rank')
        if (not isinstance(actual, int) or isinstance(actual, bool) or actual != 32 or
                manifest.get('rank_metadata', {}).get(kind) != actual):
            raise EvidenceError(f'{kind} rank metadata disagrees with computed output')
        catalogue_id = LEGACY_IDS[f'{kind} catalogue has 32 entries']
        catalogue = by_id[catalogue_id]['result_payload']
        rows = catalogue.get('catalogue')
        if (not isinstance(rows, list) or
                not all(isinstance(row, dict) and isinstance(row.get('K'), int)
                        for row in rows) or
                catalogue.get('catalogue_size') != len(rows) or len(rows) != 32 or
                catalogue.get('target_rank_counts') !=
                [sum(row['K'] == k for row in rows) for k in range(4)] or
                catalogue.get('target_rank_counts') != [2, 6, 10, 14]):
            raise EvidenceError(f'{kind} catalogue payload mismatch')
    if profile in ('foundations', 'quark-processes', 'gluon-processes', 'limits-positivity','full'):
        from foundation_certificates import read_certificate, verify_certificate
        attestation = _read_json(run_dir / 'certificate_attestation.json')
        if attestation.get('run_id') != manifest['run_id']:
            raise EvidenceError('certificate attestation run identity mismatch')
        if attestation.get('scientific_source_digest') != manifest['scientific_source_digest']:
            raise EvidenceError('certificate attestation source mismatch')
        expected_certificates = {}
        for species in ('quark', 'gluon'):
            certificate = read_certificate(ROOT / 'certificates' / f'{species}_rank.json')
            try:
                verdict = verify_certificate(certificate, require_current=require_current)
            except (OSError, ValueError) as exc:
                raise EvidenceError(f'{species} rank certificate invalid: {exc}') from exc
            expected_certificates[species] = verdict['certificate_sha256']
            payload = by_id[f'{species}.rank_certificate']['result_payload']
            if payload.get('certificate_sha256') != verdict['certificate_sha256'] or payload.get('determinant') != verdict['determinant']:
                raise EvidenceError(f'{species} certificate result mismatch')
        if attestation.get('certificates') != expected_certificates:
            raise EvidenceError('certificate attestation digest mismatch')
    if profile in ('quark-processes','gluon-processes','limits-positivity','full'):
        from response_fixtures import SIDIS_ROWS,DY_ROWS,row_id
        fixture_digest=file_digest(ROOT/'src'/'response_fixtures.py')
        fixture_rows={row_id(process,row):row for process,rows in
                      (('SIDIS',SIDIS_ROWS),('DY',DY_ROWS)) for row in rows}
        for check_id in PROCESS_ROW_IDS:
            row=by_id[check_id]
            K,m,channel,n,weight,phase,sign=fixture_rows[check_id]
            payload=row['result_payload']
            if row['evidence_type']!='exact_fixture' or row['claim_role']!='independent_comparison' or \
                    payload.get('fixture_sha256')!=fixture_digest or \
                    payload.get('label')!={'K':K,'m':m,'channel':channel,'orbital_rank':n} or \
                    payload.get('weight')!=weight or payload.get('phase_coefficients')!=list(phase) or \
                    payload.get('signed_prefactor')!=sign:
                raise EvidenceError('response fixture digest or evidence type mismatch')
        for check_id in PROCESS_INTEGRAL_IDS:
            row=by_id[check_id]
            payload=row['result_payload']
            if row['evidence_type']!='numerical_diagnostic' or \
                    row['claim_role']!='numerical_corroboration' or \
                    len(payload.get('orders',()))<3 or \
                    payload.get('source_row_id')!=check_id.removeprefix('integral.') or \
                    set(payload.get('absolute_residuals',{}))!={'cartesian','harmonic'} or \
                    not isinstance(payload.get('relative_scale'),(int,float)) or \
                    payload['relative_scale']<=0 or \
                    any(not isinstance(values,list) or len(values)!=len(payload['orders']) or
                        values[-1]>2e-11*payload['relative_scale']
                        for values in payload['absolute_residuals'].values()):
                raise EvidenceError('missing independent integral evidence')
        for check_id in REVERSAL_IDS:
            row=by_id[check_id]
            if row['evidence_type']!='conditional_algebra':
                raise EvidenceError('operator reversal classified incorrectly')
    if profile in ('gluon-processes','limits-positivity','full'):
        from gluon_angular import certificate,check as check_angular
        from gluon_response_fixture import ROWS
        from born_validation import verify_scan_payload,precision_set,BASE
        attestation=_read_json(run_dir/'gluon_certificate_attestation.json')
        cert_path=ROOT/'certificates'/'gluon_angular.json'
        try:check_angular(cert_path)
        except (OSError,ValueError,AssertionError) as exc:
            raise EvidenceError(f'angular certificate invalid: {exc}') from exc
        if (attestation.get('run_id')!=manifest['run_id'] or
                attestation.get('scientific_source_digest')!=manifest['scientific_source_digest'] or
                attestation.get('certificate_sha256')!=file_digest(cert_path) or
                by_id['gluon.angular.certificate']['result_payload']!=certificate()):
            raise EvidenceError('angular certificate attestation mismatch')
        fixture_sha=file_digest(ROOT/'src'/'gluon_response_fixture.py')
        for row,check_id in zip(ROWS,GLUON_ROW_IDS):
            K,m,ch,n,*_=row
            payload=by_id[check_id]['result_payload']
            if (check_id!=f'gluon.response.{K}{m}.{ch}.{n}' or
                    payload.get('fixture_sha256')!=fixture_sha or
                    (payload.get('K'),payload.get('m'),payload.get('channel'),payload.get('n'))!=(K,m,ch,n) or
                    by_id[check_id]['evidence_type']!='exact_fixture'):
                raise EvidenceError('gluon response fixture mismatch')
        if [x.removeprefix('gluon.oct.') for x in GLUON_OCT_IDS]!=[x[0] for x in certificate()['modes']]:
            raise EvidenceError('octupole angular row identity mismatch')
        for check_id in GLUON_OCT_IDS:
            if (by_id[check_id]['result_payload'].get('four_spin_rates')!=4 or
                    by_id[check_id]['evidence_type']!='exact_fixture'):
                raise EvidenceError('octupole physical-rate evidence missing')
        try:verify_scan_payload(by_id['gluon.born.dense_scan']['result_payload'])
        except (ValueError,TypeError,KeyError) as exc:
            raise EvidenceError(f'dense scan evidence invalid: {exc}') from exc
        grid=by_id['gluon.born.grid']['result_payload']
        expected_grid=[dict(BASE,theta=theta,phi=phi,helicity=helicity)
                       for theta in (.4,.8,1.7,2.7) for phi in (0.,.4,1.2)
                       for helicity in (-1.,0.,1.)]
        if (grid.get('convention')!=CONVENTION_ID or
                grid.get('case_count')!=36 or grid.get('independent_direct_count')!=36 or
                [x.get('inputs') for x in grid.get('cases',[])]!=expected_grid or
                any(x.get('convention')!=CONVENTION_ID or
                    x.get('convention_digest')!=digest(CONVENTION_ID) or
                    x.get('reference_source_sha256')!=REFERENCE_SOURCE_SHA256
                    for x in grid['cases']) or
                any(x.get('direct_residual_abs',float('inf'))>1e-8 for x in grid['cases'])):
            raise EvidenceError('independent Born grid incomplete')
        if by_id['gluon.born.precision']['result_payload']!=precision_set():
            raise EvidenceError('stale or mismatched high-precision reference')
        for name in GLUON_NEGATIVE_TESTS:
            payload=by_id['software.gluon_negative.'+name]['result_payload']
            if payload.get('test_id')!=name or payload.get('exit_code')!=0:
                raise EvidenceError('gluon negative-control evidence missing')
    if profile in ('limits-positivity','full'):
        same=lambda left,right: canonical_bytes(left)==canonical_bytes(right)
        from collinear_limits import selection,operator_projection,operations_and_tail,inclusive_born_bookkeeping
        from fourier_limits import exact_gaussian,numeric_comparison
        from local_moments import check_rank_and_charge,nuclear_bookkeeping
        from conditional_positivity import (fixed_target_checks,joint_gram_checks,
            matrix_algorithm_counterexamples,collinear_blocks,collinear_witnesses)
        from collinear_block_certificate import check as check_blocks
        from source_joint_positivity import source_mapping_check,spectral_recovery
        selections={kind:selection(kind) for kind in ('quark','gluon')}
        for key in COLLINEAR_ROW_IDS+FOURIER_EXACT_IDS+LOCAL_RANK_IDS:
            if by_id[key]['evidence_type']!='exact_identity':
                raise EvidenceError('exact limits claim downgraded')
        for key in FOURIER_NUMERIC_IDS:
            if by_id[key]['evidence_type']!='numerical_parameter_cases':
                raise EvidenceError('Fourier independent comparison type mismatch')
        for key in LOCAL_PARITY_IDS:
            if by_id[key]['evidence_type']!='conditional_algebra':
                raise EvidenceError('antiquark convention is conditional algebra')
        angular={f'limits.angular.{row["id"]}':row
                 for data in selections.values() for row in data['rows']}
        if set(angular)!=set(COLLINEAR_ROW_IDS) or any(
            not same(by_id[key]['result_payload'],row) for key,row in angular.items()):
            raise EvidenceError('collinear angular row payload mismatch')
        fixed={
          'limits.collinear.selection.quark':selections['quark'],
          'limits.collinear.selection.gluon':selections['gluon'],
          'limits.collinear.operator_projection':operator_projection(),
          'limits.collinear.operations_uv':operations_and_tail(),
          'limits.collinear.inclusive_born':inclusive_born_bookkeeping(),
        }
        if any(not same(by_id[key]['result_payload'],value) for key,value in fixed.items()):
            raise EvidenceError('collinear fixed payload mismatch')
        exact=exact_gaussian()
        for n,key in enumerate(FOURIER_EXACT_IDS):
            if not same(by_id[key]['result_payload'],{'rank':n,'component_count':1 if n==0 else 2**n,
                                               'source_labels':exact['source_labels']}):
                raise EvidenceError('Fourier exact rank payload mismatch')
        numeric=numeric_comparison()
        numerical={f'limits.fourier.numeric.{case["kind"]}.rank_{case["rank"]}.branch_{"plus" if case["branch"]==1 else "minus"}':case
                   for case in numeric['cases']}
        if set(numerical)!=set(FOURIER_NUMERIC_IDS) or any(
            not same(by_id[key]['result_payload'],value) for key,value in numerical.items()):
            raise EvidenceError('independent Fourier route/domain payload mismatch')
        if (not same(by_id['limits.fourier.inverse_normalization']['result_payload'],
                     {'inverse_normalization':exact['inverse_normalization'],
                      'scalar_transform':exact['scalar_transform'],
                      'inverse_cases':numeric['inverse_cases']}) or
            not same(by_id['limits.fourier.dimensions_conjugation']['result_payload'],
                     {'dimensions':numeric['dimensions'],
                      'max_absolute_error':numeric['max_absolute_error'],
                      'max_scaled_error':numeric['max_scaled_error'],
                      'conjugation':'F_real(b)* = F_real(-b); odd-rank phase retained'})):
            raise EvidenceError('Fourier inverse/dimension payload mismatch')
        local=check_rank_and_charge()
        local_rows={f'limits.local.rotation.N{x["N"]}.{x["bilinear"]}':x
                    for x in local['rotational_rows']}
        local_rows.update({f'limits.local.antiquark.N{x["N"]}.{x["channel"]}':x
                           for x in local['parity_rows']})
        if set(local_rows)!=set(LOCAL_RANK_IDS+LOCAL_PARITY_IDS) or any(
            not same(by_id[key]['result_payload'],value) for key,value in local_rows.items()):
            raise EvidenceError('local rank or antiquark payload mismatch')
        if not same(by_id['limits.local.nuclear_bookkeeping']['result_payload'],nuclear_bookkeeping()):
            raise EvidenceError('nuclear bookkeeping payload mismatch')
        if not same(by_id['limits.local.selection_and_moments']['result_payload'],
                    {k:v for k,v in local.items() if k not in ('rotational_rows','parity_rows')}):
            raise EvidenceError('local selection payload mismatch')
        positivity={
          'limits.positivity.fixed_target':fixed_target_checks(),
          'limits.positivity.joint_gram':joint_gram_checks(),
          'limits.positivity.counterexamples':matrix_algorithm_counterexamples(),
          'limits.positivity.collinear_blocks.quark':collinear_blocks('quark'),
          'limits.positivity.collinear_blocks.gluon':collinear_blocks('gluon'),
          'limits.positivity.collinear_witnesses':collinear_witnesses(),
          'limits.positivity.finite_k_gram.quark':spectral_recovery('quark'),
          'limits.positivity.finite_k_gram.gluon':spectral_recovery('gluon'),
          'limits.positivity.source_joint_convention':source_mapping_check(),
        }
        if any(not same(by_id[key]['result_payload'],value) for key,value in positivity.items()):
            raise EvidenceError('positivity payload mismatch')
        cert_path=ROOT/'certificates'/'collinear_blocks.json'
        verdict=check_blocks(cert_path)
        att=_read_json(run_dir/'limits_certificate_attestation.json')
        if (by_id['limits.positivity.block_certificate']['result_payload']!=verdict or
            att.get('run_id')!=manifest['run_id'] or
            att.get('scientific_source_digest')!=manifest['scientific_source_digest'] or
            att.get('certificate_sha256')!=verdict['certificate_sha256']):
            raise EvidenceError('collinear block certificate attestation mismatch')
        for name in LIMITS_NEGATIVE_TESTS:
            row=by_id['software.limits_negative.'+name]
            if row['result_payload']!={'test_id':name,'exit_code':0}:
                raise EvidenceError('limits negative-control evidence mismatch')
    if profile == 'full':
        from current_index_conventions import run as current_run
        from source_joint_positivity import spectral_index_check,source_mapping_check
        expected_conventions={
            'convention.current_order':current_run(),
            'convention.sidis_order':current_run()['sidis'],
            'convention.source_spectral':spectral_index_check(),
            'convention.auxiliary_map':source_mapping_check(),
        }
        if any(canonical_bytes(by_id[key]['result_payload'])!=canonical_bytes(value)
               for key,value in expected_conventions.items()):
            raise EvidenceError('convention diagnostic payload mismatch')
        assert_convention_review(manifest,by_id['convention.current_order']['result_payload'])
        for name in CONVENTION_NEGATIVE_TESTS:
            row=by_id['software.convention_negative.'+name]
            if row['result_payload']!={'test_id':name,'exit_code':0}:
                raise EvidenceError('convention negative-control evidence mismatch')
    for label, check_id in LEGACY_IDS.items():
        if by_id[check_id]['result_payload'].get('legacy_label') != label:
            raise EvidenceError('legacy ID-to-check mapping mismatch')
    for name in SOFTWARE_TESTS:
        payload = by_id['software.' + name]['result_payload']
        if payload.get('test_id') != name or payload.get('exit_code') != 0:
            raise EvidenceError('software test payload mismatch')
    born = by_id['born.on_shell_and_conservation']['result_payload']
    if (born.get('inputs') != REFERENCE_INPUTS or
            born['report'].get('convention') != CONVENTION_ID or
            born['report'].get('convention_digest') != digest(CONVENTION_ID) or
            born['report'].get('reference_source_sha256') != REFERENCE_SOURCE_SHA256):
        raise EvidenceError('Born reference input mismatch')
    reference_report = born.get('report')
    if (not isinstance(reference_report, dict) or
            not isinstance(reference_report.get('checks'), dict) or
            reference_report.get('inputs') != REFERENCE_INPUTS):
        raise EvidenceError('Born report payload mismatch')
    for name in BORN_CHECKS:
        payload = by_id['born.' + name]['result_payload']
        if (payload.get('report') != reference_report or
                reference_report.get('checks', {}).get(name) is not True):
            raise EvidenceError('Born diagnostic payload mismatch')
    grid = by_id['grid.complete']['result_payload']
    expected_cases = math.prod(len(axis) for axis in GRID_AXES.values())
    expected_points = list(itertools.product(GRID_AXES['theta'], GRID_AXES['phi'],
                                             GRID_AXES['helicity']))
    cases = grid.get('cases')
    if (not isinstance(cases, list) or
            grid.get('number_of_cases') != expected_cases or
            len(cases) != expected_cases or
            not all(isinstance(case, dict) and
                    case.get('convention') == CONVENTION_ID and
                    case.get('convention_digest') == digest(CONVENTION_ID) and
                    case.get('reference_source_sha256') == REFERENCE_SOURCE_SHA256 and
                    case.get('checks') == {name: True for name in BORN_CHECKS} and
                    case.get('inputs') == dict(REFERENCE_INPUTS, theta=theta, phi=phi,
                                               helicity=helicity)
                    for case, (theta, phi, helicity) in zip(cases, expected_points))):
        raise EvidenceError('incomplete or failed grid cases')
    return manifest, document
