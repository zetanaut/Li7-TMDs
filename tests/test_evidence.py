"""Negative controls for run identity, completeness, and stale evidence."""
from __future__ import annotations
import copy
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
import uuid

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'scripts'))
from validation_evidence import (EvidenceError, atomic_json, file_digest, new_run,
                                 validate_run)  # noqa: E402
from update_report_summary import render_summary  # noqa: E402
from check_public_scope import forbidden, forbidden_bytes  # noqa: E402
from validation_manifest import PROFILE_REQUIRED  # noqa: E402


class EvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix='li7-evidence-')
        cls.root = Path(cls.temp.name)
        proc = cls.run_cli('--profile', 'baseline')
        assert proc.returncode == 0, proc.stdout + proc.stderr
        cls.valid = cls.single_new_run()
        validate_run(cls.valid)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    @classmethod
    def run_cli(cls, *args):
        return subprocess.run([sys.executable, str(ROOT / 'scripts/run_validation.py'),
                               '--output-root', str(cls.root / 'attempts'), *args],
                              cwd=ROOT, text=True, capture_output=True)

    @classmethod
    def single_new_run(cls):
        return max((cls.root / 'attempts').iterdir(), key=lambda path: path.stat().st_mtime_ns)

    def clone(self):
        target = self.root / ('edited-' + uuid.uuid4().hex)
        shutil.copytree(self.valid, target)
        return target

    def mutate(self, folder, *, manifest=None, results=None):
        manifest_path = folder / 'manifest.json'
        results_path = folder / 'results.json'
        if results is not None:
            atomic_json(results_path, results)
        if manifest is not None:
            if results is not None:
                manifest['results_digest'] = file_digest(results_path)
            atomic_json(manifest_path, manifest)

    def assert_rejected(self, folder, text):
        with self.assertRaisesRegex(EvidenceError, re.escape(text)):
            validate_run(folder)
        with self.assertRaises(EvidenceError):
            render_summary(folder)

    def test_imports_are_inert_with_unrelated_arguments(self):
        directory = self.root / 'import-probe'
        directory.mkdir(exist_ok=True)
        code = ('import sys; sys.argv=["host", "--unrelated-option"]; '
                'import spin32_symbolic, validate_spin32_tmd, gluon_born_response, '
                'reproduce_gluon_grid, reproduce_symbolic_validation, reproduce_gluon_point')
        env = dict(os.environ, PYTHONPATH=os.pathsep.join([str(ROOT / 'src'),
                    str(ROOT / 'examples')]), PYTHONDONTWRITEBYTECODE='1')
        proc = subprocess.run([sys.executable, '-c', code], cwd=directory, env=env,
                              text=True, capture_output=True)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(proc.stdout, '')
        self.assertEqual(list(directory.iterdir()), [])

    def test_invalid_born_cannot_reuse_old_pass(self):
        old = self.root / 'old-barn-report.json'
        old.write_text('{"checks":{"old":"PASS"}}')
        before = old.read_bytes()
        proc = subprocess.run([sys.executable, str(ROOT / 'src/gluon_born_response.py'),
                               '--sqrt-s', '2', '--output', str(old)],
                              text=True, capture_output=True)
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn('Require finite inputs', proc.stderr)
        self.assertEqual(old.read_bytes(), before)
        failed = self.run_cli('--born-sqrt-s', '2')
        self.assertEqual(failed.returncode, 1)
        folder = self.single_new_run()
        self.assertEqual(json.loads((folder / 'manifest.json').read_text())['run_state'], 'FAILED')
        self.assert_rejected(folder, 'run is FAILED/FAIL')

    def test_missing_partial_corrupt_failed_and_incomplete(self):
        missing = self.clone()
        (missing / 'results.json').unlink()
        self.assert_rejected(missing, 'missing results file')
        partial = self.clone()
        manifest = json.loads((partial / 'manifest.json').read_text())
        manifest['run_state'] = 'RUNNING'
        self.mutate(partial, manifest=manifest)
        self.assert_rejected(partial, 'run is RUNNING/PASS')
        corrupt = self.clone()
        (corrupt / 'results.json').write_text('{broken')
        manifest = json.loads((corrupt / 'manifest.json').read_text())
        manifest['results_digest'] = file_digest(corrupt / 'results.json')
        self.mutate(corrupt, manifest=manifest)
        self.assert_rejected(corrupt, 'missing or corrupt report')
        failed = self.clone()
        manifest = json.loads((failed / 'manifest.json').read_text())
        manifest.update(run_state='FAILED', status='FAIL')
        self.mutate(failed, manifest=manifest)
        self.assert_rejected(failed, 'run is FAILED/FAIL')
        failed_check = self.clone()
        manifest = json.loads((failed_check / 'manifest.json').read_text())
        document = json.loads((failed_check / 'results.json').read_text())
        document['results'][0]['status'] = 'FAIL'
        self.mutate(failed_check, manifest=manifest, results=document)
        self.assert_rejected(failed_check, 'failed or inconclusive check')
        incomplete = self.clone()
        manifest = json.loads((incomplete / 'manifest.json').read_text())
        manifest.update(run_state='INCOMPLETE', status='MISSING')
        self.mutate(incomplete, manifest=manifest)
        self.assert_rejected(incomplete, 'run is INCOMPLETE/MISSING')

    def test_source_input_and_run_identity_mismatch(self):
        folder = self.clone()
        manifest = json.loads((folder / 'manifest.json').read_text())
        manifest['source_revision'] = '07aee1d96de4c4db7a6d021c3eb4a52cc7c2bfe6'
        self.mutate(folder, manifest=manifest)
        self.assert_rejected(folder, 'source revision mismatch')
        for field, message in [('scientific_source_digest', 'scientific source digest mismatch'),
                               ('input_digest', 'input digest or baseline input mismatch')]:
            with self.subTest(field=field):
                folder = self.clone()
                manifest = json.loads((folder / 'manifest.json').read_text())
                manifest[field] = '0' * 64
                self.mutate(folder, manifest=manifest)
                self.assert_rejected(folder, message)
        folder = self.clone()
        manifest = json.loads((folder / 'manifest.json').read_text())
        manifest['reference_sha256'] = '0' * 64
        self.mutate(folder, manifest=manifest)
        self.assert_rejected(folder, 'reference version mismatch')
        for value in (None,'0'*64):
            folder = self.clone()
            manifest = json.loads((folder / 'manifest.json').read_text())
            document = json.loads((folder / 'results.json').read_text())
            for item in document['results']:
                if item['check_id'].startswith('born.'):
                    item['result_payload']['report']['convention_digest'] = value
            self.mutate(folder, manifest=manifest, results=document)
            self.assert_rejected(folder, 'Born reference input mismatch')
        folder = self.clone()
        manifest = json.loads((folder / 'manifest.json').read_text())
        manifest['input_spec']['convention_id'] = 'born-current-v1-legacy-source-label'
        self.mutate(folder, manifest=manifest)
        self.assert_rejected(folder, 'input digest or baseline input mismatch')
        folder = self.clone()
        manifest = json.loads((folder / 'manifest.json').read_text())
        manifest['run_id'] = 'different'
        self.mutate(folder, manifest=manifest)
        self.assert_rejected(folder, 'results run_id mismatch')

    def test_missing_duplicate_unexpected_and_empty_ids(self):
        for case, message in [('missing', 'missing or unexpected required check ID'),
                              ('duplicate', 'duplicate check ID'),
                              ('unexpected', 'missing or unexpected required check ID'),
                              ('empty', 'empty or malformed result collection')]:
            with self.subTest(case=case):
                folder = self.clone()
                manifest = json.loads((folder / 'manifest.json').read_text())
                document = json.loads((folder / 'results.json').read_text())
                if case == 'missing':
                    document['results'].pop()
                    manifest['executed_check_ids'] = [x['check_id'] for x in document['results']]
                elif case == 'duplicate':
                    document['results'].append(copy.deepcopy(document['results'][0]))
                elif case == 'unexpected':
                    document['results'][-1]['check_id'] = 'unregistered'
                    manifest['executed_check_ids'] = [x['check_id'] for x in document['results']]
                else:
                    document['results'] = []
                    manifest['executed_check_ids'] = []
                self.mutate(folder, manifest=manifest, results=document)
                self.assert_rejected(folder, message)

    def test_interrupted_run_rejected(self):
        folder, _ = new_run(self.root / 'interrupted', 'baseline')
        self.assert_rejected(folder, 'run is RUNNING/INCOMPLETE')

    def test_rank_metadata_and_nonfinite_payload_rejected(self):
        folder = self.clone()
        manifest = json.loads((folder / 'manifest.json').read_text())
        manifest['rank_metadata']['quark'] = 31
        self.mutate(folder, manifest=manifest)
        self.assert_rejected(folder, 'quark rank metadata disagrees')
        folder = self.clone()
        manifest = json.loads((folder / 'manifest.json').read_text())
        document = json.loads((folder / 'results.json').read_text())
        document['results'][0]['result_payload']['bad_float'] = float('inf')
        # Force invalid JSON while updating the enclosing file digest.
        (folder / 'results.json').write_text(json.dumps(document))
        manifest['results_digest'] = file_digest(folder / 'results.json')
        self.mutate(folder, manifest=manifest)
        self.assert_rejected(folder, 'non-finite JSON constant')

    def test_mixed_successful_runs_rejected(self):
        proc = self.run_cli('--profile', 'baseline')
        self.assertEqual(proc.returncode, 0, proc.stderr)
        other = self.single_new_run()
        folder = self.clone()
        shutil.copy2(other / 'results.json', folder / 'results.json')
        manifest = json.loads((folder / 'manifest.json').read_text())
        manifest['results_digest'] = file_digest(folder / 'results.json')
        self.mutate(folder, manifest=manifest)
        self.assert_rejected(folder, 'results run_id mismatch')

    def test_grid_missing_duplicate_and_failed_point_rejected(self):
        for mutation in ('missing', 'duplicate', 'failed', 'legacy-convention'):
            with self.subTest(mutation=mutation):
                folder = self.clone()
                manifest = json.loads((folder / 'manifest.json').read_text())
                document = json.loads((folder / 'results.json').read_text())
                grid = next(item for item in document['results']
                            if item['check_id'] == 'grid.complete')['result_payload']
                cases = grid['cases']
                if mutation == 'missing':
                    cases.pop()
                elif mutation == 'duplicate':
                    cases[-1] = copy.deepcopy(cases[0])
                else:
                    if mutation == 'failed':
                        cases[-1]['checks']['photon_Ward'] = False
                    else:
                        cases[-1]['convention'] = 'born-current-v1-legacy-source-label'
                self.mutate(folder, manifest=manifest, results=document)
                self.assert_rejected(folder, 'incomplete or failed grid cases')

    def test_incomplete_full_profile_rejected(self):
        folder = self.clone()
        manifest = json.loads((folder / 'manifest.json').read_text())
        manifest['profile']='full'
        manifest['required_check_ids']=list(PROFILE_REQUIRED['full'])
        manifest['missing_check_ids']=sorted(set(PROFILE_REQUIRED['full'])-
                                             set(manifest['executed_check_ids']))
        manifest['run_state']='INCOMPLETE'
        manifest['status']='MISSING'
        self.mutate(folder,manifest=manifest)
        self.assertEqual(manifest['run_state'], 'INCOMPLETE')
        self.assertTrue(manifest['missing_check_ids'])
        self.assert_rejected(folder, 'run is INCOMPLETE/MISSING')

    def test_public_scope_filters_handoff_artifacts(self):
        for path in ('reference/spin32_tmd_prd.tex', 'assets/manuscript/figure.png',
                     'docs/handoff.zip', 'site-copy.tar.gz', 'notes/CODEX_prompt.txt'):
            self.assertTrue(forbidden(path), path)
        self.assertTrue(forbidden('sitemap.xml.gz'))
        self.assertFalse(forbidden('sitemap.xml.gz', site=True))
        self.assertFalse(forbidden('figures/b_U.png', site=True))
        self.assertTrue(forbidden_bytes(b'%PDF-1.7'))
        self.assertTrue(forbidden_bytes(b'PK\x03\x04hidden archive'))
        self.assertTrue(forbidden_bytes(b'\\documentclass{article}'))
        self.assertFalse(forbidden_bytes(b'# Scientific validation\n'))


if __name__ == '__main__':
    unittest.main()
