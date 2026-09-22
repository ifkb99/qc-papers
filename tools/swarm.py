#!/usr/bin/env python3
"""Shortcuts for the swarm workflow in SWARM.md; each replaces several arb/shell steps.

Run from the research root with `uv run python tools/swarm.py COMMAND`:

  sweep                  coordinator round start: board attention, open questions,
                         held resources, and board IDs cited in the ledger that the
                         board cannot resolve (read-only; no inbox receipts)
  todo 50                every board assignment against one TODO, with closure records
  review-start S_ID      referee packet: task acceptance, frozen manifest and evidence lint
  check LOG -- CMD ...   run CMD, write its output to LOG, print the arb check object
                         (--append FILE adds it to a JSON list for `checks:=@FILE`)
  doc-gate OUTDIR        documentation gate logs plus the arb checks list
                         (reindex --check by default; --reindex regenerates, coordinator only)

Nothing here creates, claims, reviews or closes board records, and nothing runs
science. arb remains the record; these commands only read it or format checks.
"""
from __future__ import annotations

import argparse
import json
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BOARD_ID = re.compile(r"(?<![A-Za-z0-9])[TASVCMRL][0-9a-f]{16}(?![A-Za-z0-9])")
# Hand-written homes only: generated aggregates repeat the notes they are built from.
CITATION_SOURCES = ["claims", "notes", "todo", "HANDOFF.md", "METHOD.md", "SWARM.md", "SWARM_REVIEW_*.md", "PAPER_*.md", "README.md"]
GENERATED = {"CLAIMS.md", "NOTES.md", "INDEX.md"}


class ArbError(RuntimeError):
    pass


def arb(project: Path, *args: str, timeout: int = 60) -> dict:
    exe = shutil.which("arb")
    if not exe:
        raise ArbError("arb is not on PATH")
    proc = subprocess.run([exe, "--project", str(project), *args], capture_output=True, text=True,
                          stdin=subprocess.DEVNULL, timeout=timeout)
    try:
        data = json.loads(proc.stdout)
    except json.JSONDecodeError:
        raise ArbError(f"arb {' '.join(args)} returned no JSON: {proc.stderr.strip()[:300]}") from None
    if not data.get("ok"):
        error = data.get("error", {})
        raise ArbError(f"arb {' '.join(args[:2])}: {error.get('code')}: {error.get('message')}"
                       + (f" ({error['hint']})" if error.get("hint") else ""))
    return data["result"]


def citations(root: Path, sources=CITATION_SOURCES) -> dict[str, list[str]]:
    """Board ID -> ['path:line', ...] for IDs cited in hand-written project documents."""
    found: dict[str, list[str]] = {}
    for source in sources:
        base = root / source
        files = sorted(base.rglob("*.md")) if base.is_dir() else sorted(root.glob(source))
        for path in files:
            if path.name in GENERATED:
                continue
            for number, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
                for match in BOARD_ID.findall(line):
                    found.setdefault(match, []).append(f"{path.relative_to(root)}:{number}")
    return found


def find_todo(root: Path, key: str) -> Path:
    matches = [p for folder in ("todo/open", "todo/done") for p in sorted((root / folder).glob(f"{key}-*.md"))]
    if not matches:
        raise SystemExit(f"no todo/open or todo/done file starts with {key}-")
    if len(matches) > 1:
        raise SystemExit("ambiguous TODO key; candidates: " + ", ".join(str(p.relative_to(root)) for p in matches))
    return matches[0]


def all_pages(project: Path, op: str, **fields) -> list[dict]:
    items, offset = [], 0
    while offset is not None:
        pairs = [f"{k}={v}" for k, v in fields.items()]
        page = arb(project, "call", op, *pairs, f"offset={offset}", "limit=100")
        items += page["items"]
        offset = page["next_offset"]
    return items


# ------------------------------------------------------------------ commands --

def cmd_sweep(ns) -> int:
    board = arb(ns.project, "call", "board.read", "limit=100")
    cited = citations(ns.project)
    checked = arb(ns.project, "call", "record.check", "ids:=" + json.dumps(sorted(cited))) if cited else {"found": {}, "missing": []}
    missing = {key: cited[key] for key in checked["missing"]}
    if ns.json:
        print(json.dumps({"attention": board["attention"], "questions": board["questions"], "resources": board["resources"],
                          "open_tasks": board["task_count"], "cited_ids": len(cited), "missing_citations": missing}, indent=2))
    else:
        print(f"Board: {board['task_count']} open tasks, {board['attention_count']} need attention, "
              f"{board['question_count']} open questions, {len(board['resources'])} held resources")
        for item in board["attention"]:
            label = item.get("task") or item.get("resource")
            print(f"  {label}  {item.get('state', 'resource')}: {'; '.join(item['reasons'])}")
            if item.get("title"):
                print(f"      {item['title'][:100]}")
        for question in board["questions"]:
            print(f"  question {question['id']} [{question['topic']}]: {question['preview'][:120]}")
        for resource in board["resources"]:
            print(f"  resource {resource['name']} held by {resource['owner']} ({resource['id']})")
        print(f"Ledger citations: {len(cited)} board IDs cited, {len(missing)} unresolved on this board")
        for key, where in sorted(missing.items()):
            print(f"  {key}: {', '.join(where[:3])}" + (f" (+{len(where)-3} more)" if len(where) > 3 else ""))
        if board["attention"]:
            print("Clear attention items before new work: review, reclaim (after confirming the worker stopped) or cancel.")
    return 1 if ns.strict and (board["attention"] or missing) else 0


def cmd_todo(ns) -> int:
    path = find_todo(ns.project, ns.key)
    rel = str(path.relative_to(ns.project))
    text = path.read_text()
    title = re.search(r'^title:\s*"?(.*?)"?\s*$', text, re.M)
    outcome = re.search(r'^outcome:\s*"?(.*?)"?\s*$', text, re.M)
    tasks = all_pages(ns.project, "task.list", source_ref=f"todo/{path.parent.name}/{ns.key}-")
    if path.parent.name == "done":  # tasks created while the TODO was still open keep that source_ref
        tasks += all_pages(ns.project, "task.list", source_ref=f"todo/open/{ns.key}-")
    rows = []
    for task in tasks:
        row = {k: task[k] for k in ("id", "state", "owner", "title", "mode", "submission")}
        if ns.records and task["state"] == "closed":
            history = arb(ns.project, "call", "task.history", f"task={task['id']}")
            recorded = [c for c in history["closures"] if c.get("records") is not None]
            row["records"] = sorted({r for c in recorded for r in c["records"]}) if recorded else None
            row["records_note"] = next((c["records_note"] for c in recorded if c.get("records_note")), None)
        rows.append(row)
    if ns.json:
        print(json.dumps({"todo": rel, "title": title and title.group(1), "tasks": rows}, indent=2))
        return 0
    print(f"{rel}: {title.group(1) if title else ''}")
    if outcome:
        print(f"  outcome: {outcome.group(1)}")
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["state"]] = counts.get(row["state"], 0) + 1
    print(f"  {len(rows)} board tasks: " + ", ".join(f"{n} {s}" for s, n in sorted(counts.items())))
    for row in rows:
        extra = ""
        if "records" in row:
            extra = ("  records: n/a (closed before the records policy)" if row["records"] is None else
                     f"  records: {', '.join(row['records']) or 'none — ' + (row['records_note'] or '')[:60]}")
        print(f"  {row['id']}  {row['state']:9} {row['mode']:8} {(row['owner'] or '-')[:28]:28} {row['title'][:90]}{extra}")
    print(f"  more: arb search <terms> --kind message --kind review; arb show <ID>")
    return 0


def lint_result_error(result, exit_code: int) -> str | None:
    """Require the evidence linter's complete result before reporting success."""
    if not isinstance(result, dict) or type(result.get("errors")) is not int or not isinstance(result.get("findings"), list):
        return "expected an object with integer errors and a findings list"
    items = result["findings"]
    if any(not isinstance(item, dict) or
           any(not isinstance(item.get(key), str) for key in ("level", "code", "where", "message")) or
           item["level"] not in {"ERROR", "WARNING", "WARN", "ADVISE", "INFO"} for item in items):
        return "malformed lint finding"
    errors = sum(item["level"] == "ERROR" for item in items)
    if result["errors"] != errors:
        return "error count disagrees with lint findings"
    if (exit_code == 0) != (errors == 0):
        return "lint exit code disagrees with lint findings"
    return None


def cmd_review_start(ns) -> int:
    submission = arb(ns.project, "show", ns.submission)["record"]
    manifest = submission["manifest"]
    task = arb(ns.project, "show", submission["task"])["record"]
    spec = task["task"]["spec"]
    lint = subprocess.run([sys.executable, str(ns.project / "tools" / "evidence_lint.py"), "--project", str(ns.project),
                           "--submission", ns.submission, "--json"], capture_output=True, text=True, timeout=300)
    try:
        findings = json.loads(lint.stdout)
        lint_error = lint_result_error(findings, lint.returncode)
    except json.JSONDecodeError:
        findings = None
        lint_error = "unreadable JSON: " + (lint.stdout + lint.stderr).strip()[-500:]
    exit_code = lint.returncode or (2 if lint_error else 0)
    if ns.json:
        packet = {"submission": submission, "task": task["task"], "lint": findings, "lint_exit": lint.returncode}
        if lint_error:
            packet["lint_error"] = lint_error
        print(json.dumps(packet, indent=2))
        return exit_code
    current = task["task"]["submission"] == ns.submission and task["task"]["state"] == "submitted"
    print(f"Submission {ns.submission} by {submission['author']} for {submission['task']} ({spec['mode']}/{spec['kind']})")
    print(f"  title: {spec['title']}")
    print(f"  reviewer: {spec['reviewer']}  independent_review: {spec['independent_review']}  "
          f"current candidate: {'yes' if current else 'NO — state ' + task['task']['state']}")
    print(f"  source_ref: {spec['source_ref']}")
    print("Acceptance:")
    for item in spec["acceptance"]:
        print(f"  - {item}")
    print(f"Summary: {manifest['summary']}")
    print(f"Limitations: {manifest['limitations']}")
    print(f"Outcome: {manifest['outcome']}")
    print("Evidence (read the archive, not the working file):")
    for ref in manifest["evidence"]:
        print(f"  {ref['source']}  {ref['bytes']} B  sha {ref['sha256'][:12]}  archive {ref['archive']}")
    for check in manifest["checks"]:
        print(f"  check exit {check['exit_code']}: {check['command']}  log {check['log']}")
    for run in manifest["runs"]:
        print(f"  run {run['id']} exit {run.get('exit_code')}: {run.get('control_result', '')[:120]}")
    if submission.get("candidate_changed"):
        print(f"  WARNING canonical candidates changed since submission: {', '.join(submission['candidate_changed'])}")
    print(f"Evidence lint exit {lint.returncode}:")
    if lint_error:
        print(f"  lint result invalid: {lint_error}")
    else:
        items = findings["findings"]
        grouped: dict[tuple[str, str], list[dict]] = {}
        for item in items:
            if ns.verbose or item["level"] in ("ERROR", "WARNING", "WARN"):
                print(f"  {item['level']} {item['code']} {item['where']}: {item['message']}")
            else:
                grouped.setdefault((item["level"], item["code"]), []).append(item)
        for (level, code), group in sorted(grouped.items()):
            files = sorted({g["where"].split(":")[0] for g in group})
            print(f"  {level} {code} x{len(group)} in {', '.join(files[:4])}" + (f" (+{len(files)-4} files)" if len(files) > 4 else "")
                  + f": {group[0]['message']}")
        if grouped:
            print("  (grouped; --verbose lists each line)")
        if not items:
            print("  no findings")
    return exit_code


def output_path(project: Path, value: str, suffix: str) -> Path:
    """A new-or-owned output file under out/, never board storage or a canonical record."""
    project = project.resolve()
    path = (project / value).resolve()
    storage = [project / "out/agent-board/state", project / "out/agent-board/artifacts", project / ".agent-board"]
    if not path.is_relative_to(project / "out") or any(path.is_relative_to(s) for s in storage):
        raise SystemExit(f"{value} must be under out/ (outside board state and artifacts) so arb can archive it")
    if path.suffix != suffix:
        raise SystemExit(f"{value} must end in {suffix}")
    return path


def run_check(project: Path, log: str, command: list[str]) -> dict:
    output_path(project, log, ".log")  # validate the resolved destination
    requested = project / log
    # Resolve parents, but retain the leaf so exclusive creation also refuses
    # a dangling symlink instead of following it to a new file.
    path = requested.parent.resolve() / requested.name
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        stream = path.open("x")
    except FileExistsError:
        raise SystemExit(f"{log} already exists; give each run its own log (SWARM.md) rather than overwriting evidence") from None
    with stream:
        stream.write(f"$ {shlex.join(command)}\n")
        stream.flush()
        try:
            code = subprocess.run(command, cwd=project, stdout=stream, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL).returncode
        except OSError as error:
            stream.write(f"{type(error).__name__}: {error}\n")
            code = 127
        stream.write(f"\n[exit {code}]\n")
    return {"command": shlex.join(command), "exit_code": code, "log": str(path.relative_to(project.resolve()))}


def existing_checks(project: Path, name: str) -> list[dict]:
    """The current list in a checks file; refuses any other JSON so a report cannot be overwritten."""
    path = output_path(project, name, ".json")
    if not path.is_file():
        return []
    try:
        data = json.loads(path.read_text())
    except json.JSONDecodeError:
        data = None
    if not isinstance(data, list) or not all(isinstance(c, dict) and set(c) == {"command", "exit_code", "log"} for c in data):
        raise SystemExit(f"{name} exists and is not a checks list; choose a new --append file")
    return data


def append_checks(project: Path, name: str | None, checks: list[dict]) -> None:
    if not name:
        return
    path = output_path(project, name, ".json")
    existing = existing_checks(project, name)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(existing + checks, indent=1) + "\n")


def cmd_check(ns) -> int:
    command = ns.command
    if not command:
        raise SystemExit("give the command after --")
    if ns.append:
        existing_checks(ns.project, ns.append)  # validate before running anything
    check = run_check(ns.project, ns.log, command)
    append_checks(ns.project, ns.append, [check])
    print(json.dumps(check))
    return check["exit_code"]


def cmd_doc_gate(ns) -> int:
    reindex = ["uv", "run", "python", "tools/reindex.py"] + ([] if ns.reindex else ["--check"])
    if ns.append:
        existing_checks(ns.project, ns.append)
    for name in ("reindex", "check"):  # refuse before running either command
        if output_path(ns.project, f"{ns.outdir}/{ns.prefix}{name}.log", ".log").exists():
            raise SystemExit(f"{ns.outdir}/{ns.prefix}{name}.log already exists; choose a new --prefix")
    checks = [run_check(ns.project, f"{ns.outdir}/{ns.prefix}reindex.log", reindex),
              run_check(ns.project, f"{ns.outdir}/{ns.prefix}check.log", ["uv", "run", "python", "tools/check.py"])]
    append_checks(ns.project, ns.append, checks)
    print(json.dumps(checks))
    return max(c["exit_code"] for c in checks)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--project", type=Path, default=ROOT)
    commands = parser.add_subparsers(dest="name", required=True)
    sweep = commands.add_parser("sweep", help="round-start board and citation sweep")
    sweep.add_argument("--json", action="store_true")
    sweep.add_argument("--strict", action="store_true", help="exit 1 when anything needs attention")
    todo = commands.add_parser("todo", help="board assignments for one TODO")
    todo.add_argument("key", help="TODO number or key, e.g. 50 or 12f")
    todo.add_argument("--records", action="store_true", help="also list records named by each closure")
    todo.add_argument("--json", action="store_true")
    review = commands.add_parser("review-start", help="referee packet for one submission")
    review.add_argument("submission")
    review.add_argument("--json", action="store_true")
    review.add_argument("--verbose", action="store_true", help="include INFO lint findings")
    check = commands.add_parser("check", help="run a command into a log and print its arb check object")
    check.add_argument("--append", help="JSON list file collecting checks for checks:=@FILE")
    check.add_argument("log")
    gate = commands.add_parser("doc-gate", help="documentation gate logs and arb checks list")
    gate.add_argument("--append", help="JSON list file collecting checks for checks:=@FILE")
    gate.add_argument("outdir")
    gate.add_argument("--reindex", action="store_true", help="regenerate indexes (coordinator integration only)")
    gate.add_argument("--prefix", default="", help="log name prefix, e.g. C98_")
    argv = list(sys.argv[1:] if argv is None else argv)
    command = None
    at = 0
    while at < len(argv) and argv[at].startswith("-"):  # global options before the subcommand
        at += 2 if argv[at] == "--project" else 1
    before = argv[at:argv.index("--", at)] if "--" in argv[at:] else argv[at:]
    if at < len(argv) and argv[at] == "check" and not ({"-h", "--help"} & set(before)):  # -h after -- belongs to the command
        # Everything after the first -- following `check` is the command, verbatim,
        # including any -- of its own; options and LOG come before it in any order.
        if "--" not in argv[at:]:
            parser.error("check needs the command after --: check [--append FILE] LOG -- COMMAND ...")
        split = argv.index("--", at)
        argv, command = argv[:split], argv[split + 1:]
    ns = parser.parse_args(argv)
    if ns.name == "check":
        ns.command = command
    ns.project = ns.project.resolve()
    handler = {"sweep": cmd_sweep, "todo": cmd_todo, "review-start": cmd_review_start,
               "check": cmd_check, "doc-gate": cmd_doc_gate}[ns.name]
    try:
        return handler(ns)
    except ArbError as error:
        print(f"swarm.py: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
