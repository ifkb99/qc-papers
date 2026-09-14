"""Mechanical workflow regressions; no scientific experiment or live board mutation.

Run: uv run python tools/test_evidence_lint.py
"""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import evidence_lint as lint


GOOD = ("from lab.harness import Experiment\ne=Experiment('fixture')\n"
        "e.predict('P1','p')\ne.must_fail('C1','c')\ne.check('P1',True)\n"
        "e.fail_check('C1',True)\ne.finish()\n")
BAD = GOOD.replace("e.finish()", "print('ALL TESTS PASSED')")


class EvidenceLintTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.findings = lint.Findings()

    def codes(self):
        return {x['code'] for x in self.findings.errors}

    def entry(self, name, text):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        blob = path.read_bytes()
        sha = hashlib.sha256(blob).hexdigest()
        archive = self.root / 'blobs' / sha
        archive.parent.mkdir(exist_ok=True)
        archive.write_bytes(blob)
        return dict(source=name, archive=str(archive.relative_to(self.root)),
                    sha256=sha, bytes=len(blob))

    def submission(self, entries, **fields):
        return dict(manifest=dict(evidence=entries, limitations='mechanical test', **fields))

    def cli(self, *args):
        return subprocess.run([sys.executable, '-B', lint.__file__, '--project',
                               str(self.root), '--json', *args], capture_output=True, text=True)

    def test_good_harness_and_bad_harness_discriminate(self):
        lint.lint_text('good.py', GOOD, self.findings)
        self.assertFalse(self.findings.errors)
        lint.lint_text('bad.py', BAD, self.findings)
        self.assertTrue({'HARNESS-NO-FINISH', 'UNCONDITIONAL-SUCCESS'} <= self.codes())

    def test_unreadable_submission_is_not_a_pass(self):
        lint.lint_submission({'_error': {'code': 'candidate_changed', 'message': 'changed'}},
                             self.root, self.findings)
        self.assertIn('UNREADABLE:candidate_changed', self.codes())

    def test_json_failure_with_success_exit(self):
        lint.lint_text('report.json', json.dumps({'ok': False, 'checks': [
            {'id': 'P1', 'ok': False, 'kind': 'check'}]}), self.findings, 0)
        self.assertIn('REPORT-CONTRADICTS-EXIT', self.codes())

    def test_report_verdict_and_resolution(self):
        report = dict(ok=True, checks=[dict(id='P1', ok=False, kind='check')],
                      predictions={'P1': 'p', 'P2': 'p'}, controls={'C1': 'c'})
        lint.lint_text('report.json', json.dumps(report), self.findings)
        self.assertTrue({'REPORT-VERDICT-MISMATCH', 'REPORT-UNRESOLVED'} <= self.codes())

    def test_filename_and_header_do_not_suppress_active_errors(self):
        lint.lint_text('failed_reference.py', '# PRESERVED FAILURE\n' + BAD, self.findings)
        self.assertIn('HARNESS-NO-FINISH', self.codes())

    def test_explicit_preservation_requires_matching_hash(self):
        item = self.entry('reference_v1.py', BAD)
        manifest = self.root / 'roles.json'
        record = dict(source=item['source'], sha256=item['sha256'], reason='superseded by corrected v2')
        manifest.write_text(json.dumps(dict(schema_version=1, preserved=[record])))
        result = self.cli('reference_v1.py', '--preserved-manifest', 'roles.json')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn('PRESERVED:HARNESS-NO-FINISH', result.stdout)
        (self.root / 'reference_v1.py').write_text(BAD + '\n# changed\n')
        result = self.cli('reference_v1.py', '--preserved-manifest', 'roles.json')
        self.assertEqual(result.returncode, 1, result.stdout)
        self.assertIn('PRESERVATION-HASH-MISMATCH', result.stdout)

    def test_historical_failure_cannot_excuse_successful_check(self):
        lint.lint_text('run.log', 'Traceback (most recent call last):\n',
                       self.findings, 0, preserved=True)
        self.assertIn('EXIT-CONTRADICTS-LOG', self.codes())

    def test_metadata_for_review_must_itself_be_frozen(self):
        item = self.entry('reference_v1.py', BAD)
        lint.lint_submission(self.submission([item]), self.root, self.findings, 'roles.json')
        self.assertIn('PRESERVATION-NOT-FROZEN', self.codes())

    def test_frozen_preservation_uses_archived_metadata(self):
        item = self.entry('reference_v1.py', BAD)
        meta = self.entry('roles.json', json.dumps(dict(schema_version=1, preserved=[
            dict(source=item['source'], sha256=item['sha256'], reason='superseded') ])))
        (self.root / 'roles.json').write_text('{}')
        lint.lint_submission(self.submission([item, meta]), self.root, self.findings, 'roles.json')
        self.assertFalse(self.findings.errors, self.findings.items)
        self.assertTrue(any(x['code'].startswith('PRESERVED:') for x in self.findings.items))

    def test_check_path_cannot_choose_last_archive(self):
        old = self.entry('reused.log', 'Traceback (most recent call last):\n')
        new = self.entry('reused.log', '=== x: 1/1 checks pass ===\n')
        lint.lint_submission(self.submission([old, new], checks=[
            dict(command='x', exit_code=0, log='reused.log')]), self.root, self.findings)
        self.assertIn('CHECK-LOG-AMBIGUOUS', self.codes())

    def test_unique_logs_and_frozen_content(self):
        old = self.entry('run_v1.log', 'Traceback (most recent call last):\n')
        new = self.entry('run_v2.log', '=== x: 1/1 checks pass ===\n')
        (self.root / 'run_v2.log').write_text('working copy differs')
        lint.lint_submission(self.submission([old, new], checks=[
            dict(command='v2', exit_code=0, log='run_v2.log')], runs=[
            dict(id='R1', state='finished', exit_code=1, log=old),
            dict(id='R2', state='finished', exit_code=0, log=new)]), self.root, self.findings)
        self.assertFalse(self.findings.errors, self.findings.items)
        self.assertTrue(any(x['code'] == 'WORKING-COPY-CHANGED' for x in self.findings.items))

    def test_changed_canonical_candidate_is_reported(self):
        item = self.entry('claims/C1.md', 'claim\n')
        result = self.submission([item])
        result['candidate_changed'] = ['claims/C1.md']
        lint.lint_submission(result, self.root, self.findings)
        self.assertIn(('WARN', 'CANDIDATE-CHANGED'), {(x['level'], x['code']) for x in self.findings.items})
        self.assertFalse(self.findings.errors)

    def test_archive_corruption_is_an_error(self):
        item = self.entry('reference.py', GOOD)
        (self.root / item['archive']).write_text('corrupted')
        lint.lint_submission(self.submission([item]), self.root, self.findings)
        self.assertIn('ARCHIVE-TAMPERED', self.codes())

    def test_missing_archive_hash_cannot_waive_errors(self):
        item = self.entry('reference.py', BAD)
        del item['sha256']
        lint.lint_submission(self.submission([item]), self.root, self.findings)
        self.assertIn('ARCHIVE-INVALID', self.codes())

    def test_empty_manifest_is_not_a_pass(self):
        lint.lint_submission(self.submission([]), self.root, self.findings)
        self.assertIn('MANIFEST-INVALID', self.codes())

    def test_run_only_json_failure_with_success_exit(self):
        script = self.entry('reference.py', GOOD)
        log = self.entry('run.json', json.dumps(dict(ok=False, checks=[
            dict(id='P1', ok=False, kind='check')])))
        lint.lint_submission(self.submission([script], runs=[
            dict(id='R1', state='finished', exit_code=0, log=log)]), self.root, self.findings)
        self.assertIn('REPORT-CONTRADICTS-EXIT', self.codes())

    def test_same_bytes_different_run_path_still_checks_exit(self):
        other = self.entry('other.log', 'Traceback (most recent call last):\n')
        log = self.entry('actual.log', 'Traceback (most recent call last):\n')
        lint.lint_submission(self.submission([other], runs=[
            dict(id='R1', state='finished', exit_code=0, log=log)]), self.root, self.findings)
        self.assertIn('EXIT-CONTRADICTS-LOG', self.codes())

    def test_empty_invocation_cannot_pass(self):
        self.assertEqual(self.cli().returncode, 2)

    def test_test_data_does_not_claim_to_be_harness(self):
        lint.lint_text('unit_test.py', 'GOOD = ' + repr(GOOD), self.findings)
        self.assertFalse(self.findings.errors)


if __name__ == '__main__':
    unittest.main()
