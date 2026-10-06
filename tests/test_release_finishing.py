"""Release gate regressions against calculations and optional complete run evidence."""
from __future__ import annotations

import ast
import copy
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'scripts'))
from born_validation import dense_scan, precision_set, verify_scan_payload  # noqa: E402
from fourier_limits import exact_gaussian  # noqa: E402
from validation_evidence import EvidenceError, atomic_json, file_digest, validate_run  # noqa: E402
from validation_manifest import BORN_CHECKS  # noqa: E402
from check_publication_eligibility import require_commit_authorization  # noqa: E402


class ReleaseFinishingTests(unittest.TestCase):
    def test_commit_scoped_approval_rejects_absent_and_wrong_sha(self):
        sha = 'a' * 40
        with self.assertRaisesRegex(ValueError, 'approval is absent'):
            require_commit_authorization('false', sha, sha, sha)
        with self.assertRaisesRegex(ValueError, 'full approved commit SHA'):
            require_commit_authorization('true', 'main', sha, sha)
        with self.assertRaisesRegex(ValueError, 'revisions disagree'):
            require_commit_authorization('true', 'b' * 40, sha, sha)
        with self.assertRaisesRegex(ValueError, 'revisions disagree'):
            require_commit_authorization('true', sha, sha, 'b' * 40)
        require_commit_authorization('true', sha, sha, sha)

    def test_authorized_gate_rejects_bad_evidence_and_dirty_tree(self):
        location = os.environ.get('LI7_RELEASE_FULL_RUN_DIR')
        if not location:
            self.skipTest('set LI7_RELEASE_FULL_RUN_DIR to test a complete final-revision run')
        with tempfile.TemporaryDirectory(prefix='li7-authorized-gate-') as temp:
            clone = Path(temp) / 'checkout'
            subprocess.run(['git', 'clone', '--quiet', '--shared', str(ROOT), str(clone)],
                           check=True)
            revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'],
                                               cwd=clone, text=True).strip()
            base = [sys.executable, str(clone / 'scripts/check_publication_eligibility.py'),
                    '--run-dir', location, '--require-publication-authorization']
            def invoke(*options):
                return subprocess.run(base + list(options), cwd=clone,
                                      text=True, capture_output=True)
            absent = invoke('--approved-revision', revision)
            self.assertEqual(absent.returncode, 2)
            self.assertIn('approval is absent', absent.stderr)
            wrong = invoke('--publication-approval', 'true', '--approved-revision', 'b'*40)
            self.assertEqual(wrong.returncode, 2)
            self.assertIn('revisions disagree', wrong.stderr)
            approved = ['--publication-approval', 'true', '--approved-revision', revision]
            valid = invoke(*approved)
            self.assertEqual(valid.returncode, 0, valid.stderr)
            original = Path(location)
            for name, field, value, diagnostic in (
                ('failed', 'run_state', 'FAILED', 'run is FAILED'),
                ('stale', 'source_revision', '0'*40, 'source revision mismatch'),
            ):
                with self.subTest(name=name):
                    target = Path(temp) / name
                    shutil.copytree(original, target)
                    manifest = json.loads((target / 'manifest.json').read_text())
                    manifest[field] = value
                    atomic_json(target / 'manifest.json', manifest)
                    bad = invoke('--run-dir', str(target), *approved)
                    self.assertEqual(bad.returncode, 2)
                    self.assertIn(diagnostic, bad.stderr)
            with (clone / 'README.md').open('a') as stream:
                stream.write('\n<!-- disposable post-evidence edit -->\n')
            dirty = invoke(*approved)
            self.assertEqual(dirty.returncode, 2)
            self.assertIn('current working tree is not clean', dirty.stderr)

    def test_scan_requires_each_named_boolean_diagnostic(self):
        payload = dense_scan()
        verify_scan_payload(payload)
        self.assertEqual(set(payload['cases'][0]['checks']), set(BORN_CHECKS))
        mutations = {
            'missing': lambda row: row.pop('checks'),
            'empty': lambda row: row.__setitem__('checks', {}),
            'removed': lambda row: row['checks'].pop(next(iter(BORN_CHECKS))),
            'replacement': lambda row: row.__setitem__('checks', {'unknown': True}),
            'integer': lambda row: row['checks'].__setitem__(next(iter(BORN_CHECKS)), 1),
        }
        for name, mutate in mutations.items():
            with self.subTest(name=name):
                bad = copy.deepcopy(payload)
                mutate(bad['cases'][0])
                with self.assertRaisesRegex(ValueError, 'dense scan required diagnostics'):
                    verify_scan_payload(bad)

    def test_precision_identifier_and_no_duplicate_dictionary_key(self):
        value = precision_set()
        self.assertEqual(value['convention'], 'born-current-v2-physical-h')
        self.assertEqual(value['convention_description'],
                         'physical electron helicity h=+1 in spinor and trace')
        tree = ast.parse((ROOT / 'src/born_validation.py').read_text())
        fn = next(node for node in tree.body
                  if isinstance(node, ast.FunctionDef) and node.name == 'precision_set')
        for node in ast.walk(fn):
            if isinstance(node, ast.Dict):
                keys = [key.value for key in node.keys if isinstance(key, ast.Constant)
                        and isinstance(key.value, str)]
                self.assertEqual(len(keys), len(set(keys)))

    def test_mass_conversion_description(self):
        self.assertEqual(exact_gaussian()['mass_rescaling'],
                         'c(M_0)=(M_0/M_A)^n c(M_A)')

    def test_complete_run_rejects_missing_scan_diagnostics(self):
        location = os.environ.get('LI7_RELEASE_FULL_RUN_DIR')
        if not location:
            self.skipTest('set LI7_RELEASE_FULL_RUN_DIR to test a complete final-revision run')
        original = Path(location)
        validate_run(original)
        with tempfile.TemporaryDirectory(prefix='li7-scan-evidence-') as temp:
            target = Path(temp) / 'run'
            shutil.copytree(original, target)
            manifest = json.loads((target / 'manifest.json').read_text())
            results = json.loads((target / 'results.json').read_text())
            scan = next(row['result_payload'] for row in results['results']
                        if row['check_id'] == 'gluon.born.dense_scan')
            scan['cases'][0]['checks'].pop(next(iter(BORN_CHECKS)))
            atomic_json(target / 'results.json', results)
            manifest['results_digest'] = file_digest(target / 'results.json')
            atomic_json(target / 'manifest.json', manifest)
            with self.assertRaisesRegex(EvidenceError, 'dense scan evidence invalid'):
                validate_run(target)

    def test_complete_run_rejects_precision_convention_change(self):
        location = os.environ.get('LI7_RELEASE_FULL_RUN_DIR')
        if not location:
            self.skipTest('set LI7_RELEASE_FULL_RUN_DIR to test a complete final-revision run')
        original = Path(location)
        validate_run(original)
        with tempfile.TemporaryDirectory(prefix='li7-precision-evidence-') as temp:
            target = Path(temp) / 'run'
            shutil.copytree(original, target)
            manifest = json.loads((target / 'manifest.json').read_text())
            results = json.loads((target / 'results.json').read_text())
            precision = next(row['result_payload'] for row in results['results']
                             if row['check_id'] == 'gluon.born.precision')
            precision['convention'] = 'born-current-v1-legacy-source-label'
            atomic_json(target / 'results.json', results)
            manifest['results_digest'] = file_digest(target / 'results.json')
            atomic_json(target / 'manifest.json', manifest)
            with self.assertRaisesRegex(EvidenceError, 'high-precision reference'):
                validate_run(target)

    def test_gate_checks_live_dirty_tree(self):
        location = os.environ.get('LI7_RELEASE_FULL_RUN_DIR')
        if not location:
            self.skipTest('set LI7_RELEASE_FULL_RUN_DIR to test a complete final-revision run')
        with tempfile.TemporaryDirectory(prefix='li7-dirty-gate-') as temp:
            clone = Path(temp) / 'checkout'
            subprocess.run(['git', 'clone', '--quiet', '--shared', str(ROOT), str(clone)],
                           check=True)
            command = [sys.executable, str(clone / 'scripts/check_publication_eligibility.py'),
                       '--run-dir', location]
            clean = subprocess.run(command, cwd=clone, text=True, capture_output=True)
            self.assertEqual(clean.returncode, 0, clean.stderr)
            with (clone / 'README.md').open('a') as stream:
                stream.write('\n<!-- disposable dirty-tree regression -->\n')
            dirty = subprocess.run(command, cwd=clone, text=True, capture_output=True)
            self.assertEqual(dirty.returncode, 2)
            self.assertIn('current working tree is not clean', dirty.stderr)


if __name__ == '__main__':
    unittest.main()
