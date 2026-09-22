"""Workflow-helper regressions; disposable boards only, no science and no live board mutation.

Run: uv run python tools/test_swarm.py
Board-backed tests need the arb CLI on PATH and are skipped without it.
"""
import io
import json
from contextlib import redirect_stdout
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import swarm


def arb(root, actor, *args, **fields):
    pairs = [f"{k}:={json.dumps(v)}" for k, v in fields.items()]
    out = subprocess.run(["arb", "--project", str(root), "--actor", actor, *args, *pairs],
                         capture_output=True, text=True, stdin=subprocess.DEVNULL, timeout=60)
    data = json.loads(out.stdout)
    assert data["ok"], data
    return data["result"]


class OfflineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        for folder in ("claims", "notes", "todo/open", "todo/done"):
            (self.root / folder).mkdir(parents=True)

    def test_citations_skip_generated_aggregates_and_hash_lookalikes(self):
        (self.root / "notes/XX-a.md").write_text("review V0269acf7cd364806 and task T083e119b28fe40eb\nsha d8c19b61ccb842826d4e\n")
        (self.root / "notes/INDEX.md").write_text("V1111111111111111\n")
        (self.root / "NOTES.md").write_text("V2222222222222222\n")
        (self.root / "HANDOFF.md").write_text("see XT083e119b28fe40eb and A31e305c5e8c44063.\n")
        found = swarm.citations(self.root)
        self.assertEqual(sorted(found), ["A31e305c5e8c44063", "T083e119b28fe40eb", "V0269acf7cd364806"])
        self.assertEqual(found["V0269acf7cd364806"], ["notes/XX-a.md:1"])

    def test_find_todo_by_key_open_done_and_ambiguity(self):
        (self.root / "todo/open/50-carry.md").write_text("x")
        (self.root / "todo/done/12d-peak.md").write_text("x")
        self.assertEqual(swarm.find_todo(self.root, "50").name, "50-carry.md")
        self.assertEqual(swarm.find_todo(self.root, "12d").parent.name, "done")
        (self.root / "todo/done/50-old.md").write_text("x")
        with self.assertRaises(SystemExit):
            swarm.find_todo(self.root, "50")
        with self.assertRaises(SystemExit):
            swarm.find_todo(self.root, "99")

    def test_check_records_real_exit_and_refuses_logs_outside_project(self):
        buffer = io.StringIO()
        with redirect_stdout(buffer):
            code = swarm.main(["--project", str(self.root), "check", "out/c.log", "--", "sh", "-c", "echo visible; exit 3"])
        check = json.loads(buffer.getvalue())
        self.assertEqual((code, check["exit_code"], check["log"]), (3, 3, "out/c.log"))
        self.assertIn("visible", (self.root / "out/c.log").read_text())
        (self.root / "claims/C1.md").write_text("claim\n")
        for bad in ("../escape.log", "claims/C1.md", "claims/new.log", "out/agent-board/state/board.sqlite3",
                    "out/agent-board/artifacts/x.log", "out/c.log", "out/x.txt"):
            with self.assertRaises(SystemExit, msg=bad):
                swarm.main(["--project", str(self.root), "check", bad, "--", "true"])
        self.assertEqual((self.root / "claims/C1.md").read_text(), "claim\n")
        with redirect_stdout(io.StringIO()):
            for log, cmd in (("out/a.log", "true"), ("out/b.log", "false"), ("out/d.log", "no-such-program-xyz")):
                swarm.main(["--project", str(self.root), "check", "--append", "out/checks.json", log, "--", cmd])
        collected = json.loads((self.root / "out/checks.json").read_text())
        self.assertEqual([(c["log"], c["exit_code"]) for c in collected], [("out/a.log", 0), ("out/b.log", 1), ("out/d.log", 127)])
        self.assertIn("[exit 127]", (self.root / "out/d.log").read_text())
        with redirect_stdout(io.StringIO()):
            swarm.main(["--project", str(self.root), "check", "out/e.log", "--append", "out/checks.json", "--", "true"])
        self.assertEqual(json.loads((self.root / "out/checks.json").read_text())[-1]["log"], "out/e.log")
        (self.root / "out/report.json").write_text('[{"value": 1}]')
        for bad in ("out/report.json",):
            with self.assertRaises(SystemExit):
                swarm.main(["--project", str(self.root), "check", "--append", bad, "out/f.log", "--", "true"])
        self.assertFalse((self.root / "out/f.log").exists())  # refused before running
        self.assertEqual((self.root / "out/report.json").read_text(), '[{"value": 1}]')
        with self.assertRaises(SystemExit):
            swarm.main(["--project", str(self.root), "check", "out/g.log", "--bogus", "--", "true"])
        with redirect_stdout(io.StringIO()):
            code = swarm.main(["--project", str(self.root), "check", "--append", "out/checks.json", "out/h.log", "--",
                               "sh", "-c", 'printf "%s|" "$@"', "argv0", "--", "x"])
        self.assertEqual(code, 0)
        self.assertIn("--|x|", (self.root / "out/h.log").read_text())
        with self.assertRaises(SystemExit):
            swarm.main(["--project", str(self.root), "check", "out/i.log", "true"])
        with redirect_stdout(io.StringIO()):
            code = swarm.main(["--project", str(self.root), "check", "out/k.log", "--", "sh", "-c", "echo ran", "-h"])
        self.assertEqual(code, 0)
        self.assertIn("ran", (self.root / "out/k.log").read_text())
        (self.root / "todo/open/check-x.md").write_text("x")
        with redirect_stdout(io.StringIO()):  # "check" as a value elsewhere is not the check subcommand
            swarm.main(["--project", str(self.root), "check", "--append", "out/j.json", "out/j.log", "--", "true"])
        self.assertEqual(swarm.find_todo(self.root, "check").name, "check-x.md")
        seen = []
        original = swarm.run_check
        swarm.run_check = lambda project, log, command: seen.append(log) or {"command": "x", "exit_code": 0, "log": log}
        try:
            with redirect_stdout(io.StringIO()):
                swarm.main(["--project", str(self.root), "doc-gate", "--prefix", "check", "out/gate2"])
        finally:
            swarm.run_check = original
        self.assertEqual(seen, ["out/gate2/checkreindex.log", "out/gate2/checkcheck.log"])

    def test_doc_gate_commands_and_existing_log_refusal(self):
        seen = []
        original = swarm.run_check
        swarm.run_check = lambda project, log, command: seen.append(command) or {"command": " ".join(command), "exit_code": 0, "log": log}
        try:
            with redirect_stdout(io.StringIO()):
                swarm.main(["--project", str(self.root), "doc-gate", "out/gate"])
                swarm.main(["--project", str(self.root), "doc-gate", "--reindex", "out/gate"])
        finally:
            swarm.run_check = original
        self.assertEqual(seen[0][-1], "--check")
        self.assertEqual(seen[2][-1], "tools/reindex.py")
        self.assertEqual(seen[1][-1], "tools/check.py")
        (self.root / "out/gate").mkdir(parents=True)
        (self.root / "out/gate/check.log").write_text("old")
        with self.assertRaises(SystemExit):
            swarm.main(["--project", str(self.root), "doc-gate", "out/gate"])

    def test_log_created_during_setup_is_preserved_and_command_not_run(self):
        log = self.root / "out/race.log"
        mkdir = Path.mkdir

        def concurrent_mkdir(path, *args, **kwargs):
            mkdir(path, *args, **kwargs)
            if path == log.parent:
                log.write_text("another writer's evidence\n")

        with patch.object(Path, "mkdir", concurrent_mkdir), patch.object(swarm.subprocess, "run") as run:
            with self.assertRaises(SystemExit):
                swarm.run_check(self.root, "out/race.log", ["true"])
            run.assert_not_called()
        self.assertEqual(log.read_text(), "another writer's evidence\n")

    def test_existing_symlink_is_not_followed_when_opening_a_log(self):
        (self.root / "out").mkdir()
        target = self.root / "out/missing.log"
        (self.root / "out/link.log").symlink_to(target)
        with patch.object(swarm.subprocess, "run") as run:
            with self.assertRaises(SystemExit):
                swarm.run_check(self.root, "out/link.log", ["true"])
            run.assert_not_called()
        self.assertFalse(target.exists())

    def test_review_start_requires_readable_consistent_lint_results(self):
        submission = dict(task="T1", author="worker", manifest=dict(
            summary="bounded", limitations="one case", outcome="result", evidence=[], checks=[], runs=[]))
        task = dict(task=dict(state="submitted", submission="S1", spec=dict(
            mode="read", kind="research", title="Review", reviewer="referee", independent_review=True,
            source_ref="todo/open/50-carry.md", acceptance=["Evidence"])))
        error = dict(level="ERROR", code="BROKEN", where="report", message="bad evidence")
        info = dict(level="INFO", code="NO-HARNESS-SUMMARY", where="lint.log", message="mechanical check")
        advisory = dict(level="ADVISE", code="NOVELTY-ASSERTION", where="report", message="check scope")
        cases = (
            ("not JSON", 0, 2),
            ("null", 0, 2),
            ('{"error": "could not read evidence"}', 0, 2),
            ('{"errors": 0}', 0, 2),
            ('{"errors": 0, "findings": [null]}', 0, 2),
            (json.dumps(dict(errors=0, findings=[error])), 0, 2),
            (json.dumps(dict(errors=1, findings=[error])), 0, 2),
            (json.dumps(dict(errors=0, findings=[])), 0, 0),
            (json.dumps(dict(errors=0, findings=[info, advisory])), 0, 0),
            (json.dumps(dict(errors=1, findings=[error])), 1, 1),
            (json.dumps(dict(errors=1, findings=[error, info])), 1, 1),
            ("not JSON", 7, 7),
        )
        for stdout, lint_exit, expected in cases:
            for presentation in ([], ["--verbose"], ["--json"]):
                with self.subTest(stdout=stdout, lint_exit=lint_exit, presentation=presentation):
                    with patch.object(swarm, "arb", side_effect=[{"record": submission}, {"record": task}]), \
                            patch.object(swarm.subprocess, "run", return_value=subprocess.CompletedProcess(
                                args=[], returncode=lint_exit, stdout=stdout, stderr="")), \
                            redirect_stdout(io.StringIO()) as output:
                        code = swarm.main(["--project", str(self.root), "review-start", "S1"]
                                          + presentation)
                    self.assertEqual(code, expected, output.getvalue())
                    if presentation == ["--json"]:
                        packet = json.loads(output.getvalue())
                        self.assertEqual(packet["lint_exit"], lint_exit)
                        if expected == 2:
                            self.assertIn("lint_error", packet)
                    elif expected == 2:
                        self.assertIn("lint result invalid", output.getvalue())


@unittest.skipUnless(shutil.which("arb"), "arb CLI not installed")
class BoardTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / "arb.toml").write_text('schema_version = 1\n[records]\nclaims = "claims/*.md"\n[closure]\nrequire_records = true\n')
        for folder in ("claims", "notes", "todo/open", "tools"):
            (self.root / folder).mkdir(parents=True)
        (self.root / "todo/open/50-carry.md").write_text('---\ntitle: "Carry question"\noutcome: "open"\n---\n')
        (self.root / "claims/C1.md").write_text("claim\n")
        (self.root / "report.md").write_text("bounded report\n")
        shutil.copy(Path(__file__).with_name("evidence_lint.py"), self.root / "tools/evidence_lint.py")
        subprocess.run(["arb", "--project", str(self.root), "init"], capture_output=True, check=True, stdin=subprocess.DEVNULL)

    def run_swarm(self, *args):
        buffer = io.StringIO()
        with redirect_stdout(buffer):
            code = swarm.main(["--project", str(self.root), *args])
        return code, buffer.getvalue()

    def closed_task(self):
        task = arb(self.root, "coordinator", "call", "task.create", title="Audit carry", objective="Bound it",
                   source_ref="todo/open/50-carry.md", acceptance=["Evidence"], reviewer="referee-1")["id"]
        attempt = arb(self.root, "worker-1", "call", "task.claim", task=task)["attempt"]["id"]
        sub = arb(self.root, "worker-1", "call", "submission.create", task=task, attempt=attempt, summary="Bounded",
                  evidence=["report.md"], limitations="One case")["submission"]
        return task, sub

    def test_sweep_reports_missing_citations_and_attention(self):
        task, sub = self.closed_task()
        (self.root / "notes/XX-a.md").write_text(f"found in {sub}; lost review V0000000000000000\n")
        code, text = self.run_swarm("sweep", "--strict")
        self.assertEqual(code, 1, text)  # the unresolved citation fails strict mode
        self.assertIn("0 need attention", text)  # a fresh submission is not yet waiting past the threshold
        self.assertIn("1 unresolved", text)
        self.assertIn("V0000000000000000: notes/XX-a.md:1", text)
        code, text = self.run_swarm("sweep", "--json")
        self.assertEqual(list(json.loads(text)["missing_citations"]), ["V0000000000000000"])

    def test_read_commands_leave_the_board_unchanged(self):
        import sqlite3
        task, sub = self.closed_task()
        def counts():
            db = sqlite3.connect(self.root / ".agent-board/state/board.sqlite3")
            try:
                return [db.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0] for t in ("objects", "events", "requests", "deliveries", "cursors")]
            finally:
                db.close()
        before = counts()
        for args in (("sweep",), ("todo", "50", "--records"), ("review-start", sub)):
            self.run_swarm(*args)
        self.assertEqual(before, counts())

    def test_review_start_then_todo_lists_closure_records(self):
        task, sub = self.closed_task()
        code, text = self.run_swarm("review-start", sub)
        self.assertIn("Acceptance:\n  - Evidence", text)
        self.assertIn("current candidate: yes", text)
        self.assertIn("report.md", text)
        arb(self.root, "referee-1", "call", "review.create", submission=sub, disposition="accept", body="Checked archive")
        arb(self.root, "coordinator", "call", "task.close", task=task, summary="Recorded", records=["claims/C1.md"])
        code, text = self.run_swarm("todo", "50", "--records")
        self.assertIn("1 board tasks: 1 closed", text)
        self.assertIn("records: claims/C1.md", text)


if __name__ == "__main__":
    unittest.main()
