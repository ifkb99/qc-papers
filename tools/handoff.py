#!/usr/bin/env python3
"""Generate HANDOFF.md's state block from the systems that own the facts.

HANDOFF.md is for *where things stand*. The board owns live assignments and
their leases, git owns what is uncommitted, and `todo/open` owns what is next,
so restating any of them by hand makes HANDOFF a mirror -- the one failure mode
this repository has actually suffered (CLAUDE.md). This script writes that part;
everything outside the markers stays hand-written and carries what cannot be
derived: the current direction, why it was chosen, and the caveats that change
how a headline reads.

    uv run python tools/handoff.py            # rewrite the generated block
    uv run python tools/handoff.py --print    # print it without touching the file
    uv run python tools/handoff.py --check    # exit 1 if the markers or budgets fail

Freshness is deliberately NOT gated: the block reports live state, so it is
stale the moment anything moves, and a check that failed on that would fail on
every commit. It carries its own timestamp instead, and is regenerated when a
round opens or closes and when a session resumes. What the gate does enforce is
size, which is the failure that actually happened: HANDOFF reached 1034 lines,
96% of it history that claims, notes, todos and the board already owned.
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
from swarm import ArbError, arb  # noqa: E402  (same directory, one board client)

START = "<!-- generated:state -->"
END = "<!-- /generated:state -->"
HAND_BUDGET = 60    # lines outside the block
FILE_BUDGET = 200   # lines in total
MAX_TASKS = 12
MAX_TODOS = 20


def git(*args: str) -> str:
    out = subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, text=True, timeout=60)
    return out.stdout if out.returncode == 0 else ""


def mark(task: dict, now: float) -> str:
    """A lease only says something while the work is still in a worker's hands.

    Once a submission exists the attempt is frozen and its lease keeps ticking
    down with nothing depending on it, so reading an expired lease there as
    abandoned work is wrong -- it was, on 2026-09-17, for all four open tasks.
    """
    if task.get("submission"):
        return f"{task['state']}; lease n/a (submitted)"
    lease = task.get("lease_until")
    if not lease:
        return "unclaimed"
    if lease > now:
        return f"lease live {int((lease - now) // 60)}m"
    return f"LEASE EXPIRED {int((now - lease) // 60)}m ago"


def board_lines() -> list[str]:
    """Open assignments with lease state, attention items, questions, resources."""
    try:
        board = arb(ROOT, "call", "board.read", "limit=100")
    except (ArbError, OSError, subprocess.SubprocessError) as exc:
        return [f"- board unreadable: {exc}. Check the daemon before trusting this block."]

    now = time.time()
    # board.read lists the open tasks. task.list pages from the oldest, so its first
    # page missed every task after the hundredth and reported an open board as empty.
    open_tasks = [t for t in board["tasks"] if t["state"] not in ("closed", "cancelled")]
    lines = [f"- {len(open_tasks)} open assignment(s); "
             f"{board['attention_count']} need attention; {board['question_count']} open question(s)."]
    for task in open_tasks[:MAX_TASKS]:
        lines.append(f"  - `{task['id']}` {task['state']:10} {(task['owner'] or '-')[:24]:24} "
                     f"[{mark(task, now)}] {task['title'][:64]}")
    if len(open_tasks) > MAX_TASKS or board.get("tasks_truncated"):
        more = f"+{len(open_tasks) - MAX_TASKS}" + (" or more, board list truncated" if board.get("tasks_truncated") else "")
        lines.append(f"  - ({more}; `arb --project ROOT board`)")
    working = sorted({t["owner"] for t in open_tasks
                      if not t.get("submission") and (t.get("lease_until") or 0) > now and t.get("owner")})
    waiting = [t for t in open_tasks if t.get("submission")]
    lines.append(f"- Actors holding work with a live lease (a session may be running): "
                 f"{', '.join(working) if working else 'none'}.")
    if waiting:
        lines.append(f"- {len(waiting)} submitted result(s) waiting on the coordinator, not on a worker. "
                     f"Integrate what belongs in the ledger, then close with records.")
    for item in board["attention"][:6]:
        label = item.get("task") or item.get("resource")
        lines.append(f"- attention: {label} {item.get('state', 'resource')}: {'; '.join(item.get('reasons', []))[:90]}")
    for question in board["questions"][:4]:
        lines.append(f"- open question {question['id']} [{question['topic']}]: {question['preview'][:80]}")
    for resource in board["resources"][:4]:
        lines.append(f"- resource {resource['name']} held by {resource['owner']} ({resource['id']})")
    return lines


def tree_lines() -> list[str]:
    status = [ln for ln in git("status", "--porcelain").splitlines() if ln.strip()]
    if not status:
        return ["- Working tree clean."]
    modified = [ln[3:] for ln in status if not ln.startswith("??")]
    untracked = [ln[3:] for ln in status if ln.startswith("??")]
    branch = git("rev-parse", "--abbrev-ref", "HEAD").strip()
    lines = [f"- Branch `{branch}`: {len(modified)} modified, {len(untracked)} untracked. Nothing here is committed."]
    for label, prefix in (("claims", "claims/"), ("notes", "notes/"), ("transfers", "transfers/"),
                          ("experiments", "experiments/"),
                          ("lab modules", "lab/"), ("tools", "tools/")):
        new = [p for p in untracked if p.startswith(prefix)]
        if new:
            shown = ", ".join(f"`{Path(p).name}`" for p in new[:6])
            lines.append(f"  - untracked {label}: {shown}" + (f" (+{len(new) - 6})" if len(new) > 6 else ""))
    lines.append("  - untracked evidence under `out/` is gitignored: nothing in git reproduces its numbers.")
    return lines


def todo_lines() -> list[str]:
    rows = []
    for path in sorted((ROOT / "todo" / "open").glob("*.md")):
        text = path.read_text(errors="replace")
        title = re.search(r'^title:\s*"?(.*?)"?\s*$', text, re.M)
        key = path.name.split("-")[0]
        rows.append(f"  - {key}: {title.group(1)[:88] if title else path.stem}")
    head = [f"- {len(rows)} open item(s) in `todo/open` (titles own the question; `todo/INDEX.md` is generated):"]
    return head + rows[:MAX_TODOS] + ([f"  - (+{len(rows) - MAX_TODOS} more)"] if len(rows) > MAX_TODOS else [])


def block() -> str:
    stamp = time.strftime("%Y-%m-%d %H:%M", time.localtime())
    return "\n".join([
        START,
        f"<!-- GENERATED by tools/handoff.py at {stamp}. Do not edit inside these markers:",
        "     rerun the tool. Live state, so it is only as current as its timestamp. -->",
        "",
        "## Board", *board_lines(), "",
        "## Working tree", *tree_lines(), "",
        "## Next", *todo_lines(), "",
        END,
    ])


def check(text: str) -> list[str]:
    if START not in text or END not in text:
        # Since the 2026-09-18 migration a missing marker is a regression: without
        # this, deleting the markers would make the size budgets pass vacuously.
        return [f"the {START} / {END} markers are missing; restore them and rerun tools/handoff.py"]
    bad = []
    hand = len(re.sub(re.escape(START) + r".*?" + re.escape(END), "", text, flags=re.S).splitlines())
    total = len(text.splitlines())
    if hand > HAND_BUDGET:
        bad.append(f"hand-written part is {hand} lines, budget {HAND_BUDGET}: it is accumulating history again")
    if total > FILE_BUDGET:
        bad.append(f"HANDOFF.md is {total} lines, budget {FILE_BUDGET}")
    return bad


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="verify markers and size budgets")
    ap.add_argument("--print", dest="show", action="store_true", help="print the block, write nothing")
    ns = ap.parse_args()
    path = ROOT / "HANDOFF.md"
    text = path.read_text()

    if ns.check:
        problems = check(text)
        for problem in problems:
            print(f"HANDOFF.md: {problem}")
        return 1 if problems else 0

    if ns.show:
        print(block())
        return 0

    if START not in text or END not in text:
        print(f"HANDOFF.md has no generated markers. Restore the {START} / {END} pair\n"
              "(see todo/done/61) before this tool can write. Use --print meanwhile.", file=sys.stderr)
        return 1
    new = re.sub(re.escape(START) + r".*?" + re.escape(END), lambda _: block(), text, flags=re.S)
    path.write_text(new)
    print(f"rewrote the state block in {path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
