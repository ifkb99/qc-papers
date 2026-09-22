#!/usr/bin/env python3
"""Mechanical pre-review lint for qsim swarm evidence.

Catches mechanical defect classes from SWARM_REVIEW_2026-09-12.md: harness scripts that
never call finish(), print success unconditionally, or leave predictions and
must-fail controls unresolved; recorded check exit codes that contradict their
logs; archived evidence whose hash no longer matches; quantitative memory or
novelty assertions that need measured/primary-source backing.

It never judges mathematics, relevance of a control, or scope. A clean lint is
a precondition for review, not evidence of correctness.

Usage (from the research root):
  evidence_lint.py --submission S...      lint a frozen submission via arb
  evidence_lint.py FILE [FILE ...]        lint working files before submitting
  --json                                   machine-readable output
  --preserved-manifest FILE                hash-bound historical artifacts, for a file
                                           list; a submission's own frozen metadata
                                           is found without it, and naming it pins
Exit 1 when any ERROR is found, else 0.
No scientific code is executed. AST checks establish presence, not reachability;
dynamic checks and mathematical validity still require a reviewer.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

SUCCESS_WORDS = re.compile(r"\b(all\s+(tests|checks)\s+pass(ed)?|pass(ed)?|success(ful)?|ok|verified)\b", re.I)
MEMORY_CLAIM = re.compile(r"\b\d+(\.\d+)?\s*(KiB|MiB|GiB|KB|MB|GB)\b|\bpeak (memory|allocation)\b|\bmemory (cap|bound)\b", re.I)
NOVELTY_CLAIM = re.compile(r"\b(novel|first (to|known)|breakthrough|state[- ]of[- ]the[- ]art|no prior (work|art))\b", re.I)
LOG_FAILURE = re.compile(r"Traceback \(most recent call last\)|\bFAILED\b|^\s*FAIL\s|CONTROL DID NOT FAIL|Fatal Python error", re.M)
LOG_WARNING = re.compile(r"^\s*WARNING: .*", re.M)


class Findings:
    def __init__(self):
        self.items: list[dict] = []

    def add(self, level: str, where: str, code: str, message: str):
        self.items.append(dict(level=level, where=where, code=code, message=message))

    @property
    def errors(self):
        return [i for i in self.items if i["level"] == "ERROR"]


def _call_name(node: ast.Call) -> str:
    f = node.func
    if isinstance(f, ast.Attribute):
        return f.attr
    if isinstance(f, ast.Name):
        return f.id
    return ""


def _first_str(node: ast.Call) -> str | None:
    if node.args and isinstance(node.args[0], ast.Constant) and isinstance(node.args[0].value, str):
        return node.args[0].value
    return None


def _joined_text(node: ast.AST) -> str:
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    if isinstance(node, ast.JoinedStr):
        return "".join(v.value for v in node.values if isinstance(v, ast.Constant) and isinstance(v.value, str))
    return ""


def lint_python(path: str, text: str, f: Findings):
    try:
        tree = ast.parse(text)
    except SyntaxError as e:
        f.add("ERROR", path, "PY-SYNTAX", f"does not parse: {e}")
        return
    calls = [n for n in ast.walk(tree) if isinstance(n, ast.Call)]
    uses_harness = any(_call_name(c) == "Experiment" for c in calls) or any(
        isinstance(n, ast.ImportFrom) and n.module and n.module.split('.')[0] == 'lab'
        and any(a.name == 'Experiment' for a in n.names) for n in ast.walk(tree))
    predicted = {s for c in calls if _call_name(c) == "predict" and (s := _first_str(c))}
    controls = {s for c in calls if _call_name(c) == "must_fail" and (s := _first_str(c))}
    checked = {s for c in calls if _call_name(c) == "check" and (s := _first_str(c))}
    fail_checked = {s for c in calls if _call_name(c) == "fail_check" and (s := _first_str(c))}
    # A must_fail() whose id is built at run time (e.g. in a loop) still registers a
    # control; only a script with no must_fail() call at all lacks one.
    registers_control = any(_call_name(c) == "must_fail" for c in calls)
    dynamic_ids = any(_call_name(c) in {"check", "fail_check"} and _first_str(c) is None for c in calls)

    if uses_harness:
        if not any(_call_name(c) == "finish" for c in calls):
            f.add("ERROR", path, "HARNESS-NO-FINISH",
                  "uses Experiment but never calls finish(): no verdict, no warnings, exit code is not the result")
        if not registers_control:
            f.add("ERROR", path, "HARNESS-NO-CONTROL",
                  "no must_fail() control registered; a test only every case passes proves nothing")
        if not dynamic_ids:
            for pid in sorted(predicted - checked - fail_checked):
                f.add("ERROR", path, "HARNESS-UNCHECKED-PREDICTION", f"prediction {pid!r} is declared but never checked")
            for pid in sorted(controls - fail_checked - checked):
                f.add("ERROR", path, "HARNESS-UNRESOLVED-CONTROL", f"control {pid!r} is never resolved")
            for pid in sorted((controls - fail_checked) & checked):
                f.add("WARN", path, "CONTROL-VIA-CHECK",
                      f"control {pid!r} resolved with check(), not fail_check(): confirm the asserted condition is the control failing")
        else:
            f.add("INFO", path, "HARNESS-DYNAMIC-IDS", "check ids built dynamically; prediction/control resolution not statically verified")

    # unconditional success prints at module level (not inside if/for/def/try)
    for node in tree.body:
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call) and _call_name(node.value) == "print":
            msg = " ".join(_joined_text(a) for a in node.value.args)
            if SUCCESS_WORDS.search(msg):
                f.add("ERROR", f"{path}:{node.lineno}", "UNCONDITIONAL-SUCCESS",
                      f"module-level print({msg.strip()[:60]!r}) runs regardless of results")
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr in {"exit", "_exit"}:
            if node.args and isinstance(node.args[0], ast.Constant) and node.args[0].value == 0:
                f.add("WARN", f"{path}:{node.lineno}", "FORCED-EXIT-0", "explicit exit(0) can mask a failing verdict")


def lint_log(path: str, text: str, f: Findings, claimed_exit: int | None = None):
    # swarm.py check writes the child's actual exit as the last line. Earlier
    # markers can be child output (including nested checks), so only bind the
    # terminal marker. Plain logs need not use this wrapper format.
    marker = re.search(r"^\[exit (-?\d+)\][ \t]*(?:\r?\n)?\Z", text, re.M)
    if marker and claimed_exit is not None and int(marker.group(1)) != claimed_exit:
        f.add("ERROR", path, "EXIT-MARKER-MISMATCH",
              f"recorded exit_code {claimed_exit} but terminal log marker records {marker.group(1)}")
    failure = LOG_FAILURE.search(text)
    warnings = LOG_WARNING.findall(text)
    summary = re.findall(r"=== .*?: (\d+)/(\d+) checks pass(, (\d+) FAILED)? ===", text)
    if claimed_exit == 0 and failure:
        f.add("ERROR", path, "EXIT-CONTRADICTS-LOG",
              f"recorded exit_code 0 but log contains {failure.group(0)!r}")
    if claimed_exit not in (None, 0) and summary and all(s[0] == s[1] for s in summary):
        f.add("WARN", path, "NONZERO-EXIT-ALL-PASS",
              "nonzero exit but every harness summary passes: crash after verdict, or the wrong log")
    if claimed_exit is not None and not summary and not failure:
        f.add("INFO", path, "NO-HARNESS-SUMMARY", "log has no harness verdict line; confirm it is the log of the recorded command")
    for w in warnings:
        f.add("WARN", path, "HARNESS-WARNING", w.strip())


def lint_report(path: str, text: str, f: Findings, claimed_exit: int | None = None):
    try:
        rep = json.loads(text)
    except json.JSONDecodeError as error:
        f.add("ERROR", path, "JSON-SYNTAX", str(error))
        return
    if not isinstance(rep, dict) or "checks" not in rep:
        return
    if rep.get("ok") is False and claimed_exit == 0:
        f.add("ERROR", path, "REPORT-CONTRADICTS-EXIT", "harness report ok=false but a check recorded exit 0")
    checks = rep["checks"]
    if not isinstance(checks, list) or any(not isinstance(c, dict) or
            type(c.get("ok")) is not bool or not isinstance(c.get("id"), str) for c in checks):
        f.add("ERROR", path, "REPORT-INVALID-CHECKS", "expected check objects with string id and boolean ok")
        return
    if type(rep.get("ok")) is not bool or rep["ok"] != all(c["ok"] for c in checks):
        f.add("ERROR", path, "REPORT-VERDICT-MISMATCH", "verdict disagrees with individual checks")
    resolved = {c["id"] for c in checks}
    failed_controls = {c["id"] for c in checks if c.get("kind") == "must-fail"}
    for field, done in (("predictions", resolved), ("controls", failed_controls)):
        if field in rep:
            if not isinstance(rep[field], dict):
                f.add("ERROR", path, "REPORT-INVALID-DECLARATIONS", f"{field} must be an object")
                continue
            for pid in sorted(set(rep[field]) - done):
                f.add("ERROR", path, "REPORT-UNRESOLVED", f"{field} id {pid!r} has no matching check")
    if "controls" in rep and not rep["controls"]:
        f.add("ERROR", path, "REPORT-NO-CONTROL", "no must-fail control declared")
    for w in rep.get("warnings") or []:
        f.add("WARN", path, "REPORT-WARNING", str(w))


def lint_prose(path: str, text: str, f: Findings):
    for n, line in enumerate(text.splitlines(), 1):
        if MEMORY_CLAIM.search(line):
            f.add("ADVISE", f"{path}:{n}", "MEMORY-ASSERTION",
                  "quantitative memory statement: reviewer must see a measured peak (tracemalloc/RSS) or a proof, not source reading")
        if NOVELTY_CLAIM.search(line):
            f.add("ADVISE", f"{path}:{n}", "NOVELTY-ASSERTION",
                  "novelty/priority wording: needs a primary-source prior-art check, not an abstract")


def lint_text(path: str, text: str, f: Findings, claimed_exit: int | None = None,
              *, preserved: bool = False):
    """Only explicit, validated provenance marks an artifact historical.

    A historical artifact cannot excuse a successful check citing its failure.
    Archive, provenance and IO errors are generated outside this function and
    are never downgraded.
    """
    if preserved:
        inner = Findings()
        _lint_text(path, text, inner, claimed_exit)
        for i in inner.items:
            if claimed_exit == 0 and i["level"] == "ERROR":
                f.add(i["level"], i["where"], i["code"], i["message"])
            else:
                f.add("INFO", i["where"], "PRESERVED:" + i["code"], i["message"])
        return
    _lint_text(path, text, f, claimed_exit)


def _lint_text(path: str, text: str, f: Findings, claimed_exit: int | None = None):
    suffix = Path(path).suffix
    if suffix == ".py":
        lint_python(path, text, f)
    elif suffix == ".json":
        lint_report(path, text, f, claimed_exit)
    elif suffix == ".log" or claimed_exit is not None:
        lint_log(path, text, f, claimed_exit)
    if suffix in {".md", ".txt"}:
        lint_prose(path, text, f)


def arb_get(submission: str, project: str | None) -> dict:
    cmd = ["arb"] + (["--project", project] if project else []) + ["call", "submission.get", "--data", json.dumps({"submission": submission})]
    out = subprocess.run(cmd, capture_output=True, text=True)
    data = json.loads(out.stdout or "{}")
    if not data.get("ok"):
        return {"_error": data.get("error") or {"code": "arb_failed", "message": out.stderr.strip()}}
    return data["result"]


def _blob(root: Path, entry: dict, f: Findings) -> str | None:
    if not isinstance(entry.get("sha256"), str) or not re.fullmatch(r"[0-9a-f]{64}", entry["sha256"]):
        f.add("ERROR", entry["source"], "ARCHIVE-INVALID", "archive requires its full sha256")
        return None
    archive = root / entry.get("archive", "")
    try:
        blob = archive.read_bytes()
    except OSError:
        f.add("ERROR", entry["source"], "ARCHIVE-MISSING", f"archived blob {archive} not readable")
        return None
    if entry.get("sha256") and hashlib.sha256(blob).hexdigest() != entry["sha256"]:
        f.add("ERROR", entry["source"], "ARCHIVE-TAMPERED", "archive content does not match its recorded sha256")
    current = root / entry["source"]
    if current.is_file() and hashlib.sha256(current.read_bytes()).hexdigest() != entry.get("sha256"):
        f.add("INFO", entry["source"], "WORKING-COPY-CHANGED", "working file differs from the frozen archive; review the archive")
    return blob.decode("utf-8", "replace")


def looks_preservation(text: str) -> bool:
    """Whether an artifact is offering itself as preservation metadata.

    Any JSON object carrying a `preserved` key counts, so a malformed manifest
    reaches `preservation_records` and is reported rather than silently ignored.
    """
    try:
        data = json.loads(text)
    except ValueError:
        return False
    return isinstance(data, dict) and "preserved" in data


def preservation_records(text: str, where: str, root: Path, f: Findings) -> dict:
    """Return source/sha identities with an explicit reason, without trusting names."""
    try:
        data = json.loads(text)
        if data.get("schema_version") != 1 or not isinstance(data.get("preserved"), list):
            raise ValueError("expected schema_version=1 and preserved list")
        result = {}
        for entry in data["preserved"]:
            if not isinstance(entry, dict) or not isinstance(entry.get("source"), str) or not entry["source"]:
                raise ValueError("each artifact needs a source")
            sha = entry.get("sha256", "")
            if not isinstance(sha, str) or not re.fullmatch(r"[0-9a-f]{64}", sha):
                raise ValueError("each artifact needs its full sha256")
            if not isinstance(entry.get("reason"), str) or not entry["reason"].strip():
                raise ValueError("each artifact needs a historical reason")
            key = str((root / entry["source"]).resolve())
            if key in result:
                raise ValueError("duplicate preservation source")
            result[key] = sha
        return result
    except (ValueError, TypeError, AttributeError) as error:
        f.add("ERROR", where, "PRESERVATION-INVALID", str(error))
        return {}


def lint_submission(result: dict, root: Path, f: Findings, preserved_manifest: str | None = None):
    if "_error" in result:
        e = result["_error"]
        f.add("ERROR", "submission", "UNREADABLE:" + str(e.get("code")),
              f"{e.get('message')} (evidence was not inspected; this cannot be a clean lint)")
        return
    for path in result.get("candidate_changed") or []:
        f.add("WARN", path, "CANDIDATE-CHANGED", "canonical file differs from this frozen submission; "
              "review.create and task.close will refuse it until a new version is submitted")
    sub = result.get("manifest", result)
    if not isinstance(sub, dict) or not isinstance(sub.get("evidence"), list) or not sub["evidence"]:
        f.add("ERROR", "submission", "MANIFEST-INVALID", "expected a nonempty frozen evidence manifest")
        return
    if any(not isinstance(e, dict) or not isinstance(e.get("source"), str) or not e["source"]
           for e in sub["evidence"]):
        f.add("ERROR", "submission", "MANIFEST-INVALID", "evidence entries need source paths")
        return
    evidence = [e for e in sub.get("evidence") or [] if isinstance(e, dict) and "source" in e]
    by_source = {}
    for e in evidence:
        by_source.setdefault(e["source"], {})[e.get("sha256")] = e
    exit_by_ref = {}
    def bind_exit(log, code):
        key = (log["source"], log.get("sha256"))
        if key in exit_by_ref and exit_by_ref[key] != code:
            f.add("ERROR", log["source"], "LOG-EXIT-CONFLICT", "the same log version has contradictory exits")
        exit_by_ref[key] = code

    for c in sub.get("checks") or []:
        log = c.get("log")
        if isinstance(log, dict) and log.get("sha256"):
            bind_exit(log, c.get("exit_code"))
        elif isinstance(log, str):
            if log in by_source:
                if len(by_source[log]) != 1:
                    f.add("ERROR", log, "CHECK-LOG-AMBIGUOUS", "check path names multiple archived versions; no last-version inference")
                else:
                    bind_exit(next(iter(by_source[log].values())), c.get("exit_code"))
            else:
                f.add("ERROR", log, "CHECK-LOG-NOT-EVIDENCE", "a recorded check's log is not among the frozen evidence")
    run_logs = []
    for r in sub.get("runs") or []:
        if r.get("state") != "finished":
            f.add("ERROR", r.get("id", "run"), "RUN-NOT-FINISHED", "submitted run record was never finished")
        log = r.get("log")
        if isinstance(log, dict) and log.get("sha256"):
            bind_exit(log, r.get("exit_code"))
            if (log["source"], log["sha256"]) not in {(e["source"], e.get("sha256")) for e in evidence}:
                run_logs.append((log, r.get("exit_code")))
    versions: dict[str, set] = {}
    for e in evidence:
        versions.setdefault(e["source"], set()).add(e.get("sha256"))
    for src, shas in versions.items():
        if len(shas) > 1:
            f.add("WARN", src, "LOG-PATH-REUSED",
                  f"{len(shas)} different archived versions share this path; give each run its own log name so checks cannot point at the wrong one")
    if not (sub.get("limitations") or "").strip():
        f.add("ERROR", "submission", "NO-LIMITATIONS", "empty limitations: every submission states what was not checked")
    if not sub.get("runs") and any(e["source"].endswith(".py") for e in evidence):
        f.add("WARN", "submission", "NO-RUN-RECORDS", "script evidence without run records: execution order is unrecorded")
    # Read each archived blob once: _blob reports tampering and working-copy
    # drift, and those findings must not be duplicated by a second pass.
    texts = [(e, _blob(root, e, f)) for e in evidence]
    # Preservation metadata is frozen evidence (SWARM.md), so it is discovered
    # here rather than named on the command line: the submission already
    # carries the author's intent, and a name supplied later can disagree
    # with it. Discovery is reported, never silent.
    preserved: dict[str, str] = {}
    claimed_by: dict[str, str] = {}
    found: list[str] = []
    for e, text in texts:
        if text is None or not looks_preservation(text):
            continue
        found.append(e["source"])
        for src, sha in preservation_records(text, e["source"], root, f).items():
            if src in preserved and preserved[src] != sha:
                f.add("ERROR", src, "PRESERVATION-CONFLICT",
                      f"{claimed_by[src]} and {e['source']} preserve this source at different hashes")
                continue
            preserved[src] = sha
            claimed_by[src] = e["source"]
    if preserved_manifest:      # explicit pin: this artifact, and no other, supplies preservation
        wanted = str((root / preserved_manifest).resolve())
        named = [s for s in found if str((root / s).resolve()) == wanted]
        if not named:
            f.add("ERROR", preserved_manifest, "PRESERVATION-NOT-FROZEN",
                  "named preservation metadata is not frozen preservation evidence of this submission")
        for other in sorted(set(found) - set(named)):
            f.add("ERROR", other, "PRESERVATION-UNNAMED",
                  f"also supplies preservation, but --preserved-manifest named {preserved_manifest}")
    for src in sorted(found):
        covered = sum(1 for k, v in claimed_by.items() if v == src)
        f.add("INFO", src, "PRESERVATION-APPLIED",
              f"frozen preservation metadata; downgrades findings for {covered} artifact(s)")
    all_entries = evidence + [log for log, _ in run_logs]
    for src, sha in preserved.items():
        if not any(str((root / e["source"]).resolve()) == src and e.get("sha256") == sha for e in all_entries):
            f.add("ERROR", src, "PRESERVATION-HASH-MISMATCH", "historical identity is absent from frozen evidence")
    for e, text in texts:
        if text is not None:
            historical = bool(e.get("sha256")) and preserved.get(str((root / e["source"]).resolve())) == e["sha256"]
            lint_text(e["source"], text, f, exit_by_ref.get((e["source"], e.get("sha256"))), preserved=historical)
    for log, code in run_logs:
        text = _blob(root, log, f)
        if text is not None:
            lint_text(log["source"], text, f, code)
    for field in ("summary", "outcome"):
        if sub.get(field):
            lint_prose(f"submission.{field}", sub[field], f)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="*")
    ap.add_argument("--submission")
    ap.add_argument("--project")
    ap.add_argument("--preserved-manifest", help="hash-bound historical artifact metadata for a file list; "
                                                 "a submission's frozen metadata is discovered, and naming it pins that one artifact")
    ap.add_argument("--exit-code", type=int, help="claimed exit code for a single log file")
    ap.add_argument("--json", action="store_true")
    ns = ap.parse_args(argv)
    if not ns.submission and not ns.files:
        ap.error("provide --submission or at least one file")
    if ns.exit_code is not None and (ns.submission or len(ns.files) != 1):
        ap.error("--exit-code requires exactly one file and no submission")
    f = Findings()
    root = Path(ns.project or ".").resolve()
    if ns.submission:
        lint_submission(arb_get(ns.submission, str(root)), root, f, ns.preserved_manifest)
    preserved = {}
    if ns.files and ns.preserved_manifest:
        try:
            p = root / ns.preserved_manifest
            preserved = preservation_records(p.read_text(), str(p), root, f)
        except OSError as error:
            f.add("ERROR", ns.preserved_manifest, "PRESERVATION-UNREADABLE", str(error))
        submitted = {str((root / name).resolve()) for name in ns.files}
        for src, sha in preserved.items():
            p = Path(src)
            if src not in submitted or not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest() != sha:
                f.add("ERROR", src, "PRESERVATION-HASH-MISMATCH", "historical file is missing, changed, or absent from the file list")
    for name in ns.files:
        try:
            p = root / name
            blob = p.read_bytes()
            text = blob.decode("utf-8", "replace")
        except OSError as e:
            f.add("ERROR", name, "UNREADABLE", str(e))
            continue
        historical = preserved.get(str(p.resolve())) == hashlib.sha256(blob).hexdigest()
        lint_text(name, text, f, ns.exit_code, preserved=historical)
    if ns.json:
        print(json.dumps(dict(errors=len(f.errors), findings=f.items), indent=2))
    else:
        for i in f.items:
            print(f"{i['level']:6} {i['code']:28} {i['where']}: {i['message']}")
        print(f"-- {len(f.errors)} error(s), {len(f.items)} finding(s)")
    return 1 if f.errors else 0


if __name__ == "__main__":
    sys.exit(main())
