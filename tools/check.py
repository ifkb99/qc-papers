"""Documentation consistency gate -- suite 10.

Every check here exists because a real defect got through without it. The
2026-08-08 audit found sixteen, essentially all of one shape: a fact updated in
one place and not its mirror. Cleaning that up once was worth little; making it
fail a test is the point.

    uv run python tools/check.py          # all checks
    uv run python tools/check.py -v       # list every check as it runs

Exits nonzero on any failure so it cannot be ignored.
"""
from __future__ import annotations

import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
PAPERS = ["PAPER_A.md", "PAPER_B.md"]

# Sections of *other people's* papers we cite by number; not our headings.
EXTERNAL_SECTION_REFS = {"3.5", "5.2"}

failures: list[str] = []
checks_run = 0
VERBOSE = "-v" in sys.argv


def check(name: str):
    def deco(fn):
        global checks_run
        checks_run += 1
        try:
            problems = fn() or []
        except Exception as exc:  # a broken check is a failure, not a pass
            problems = [f"check raised {type(exc).__name__}: {exc}"]
        if problems:
            failures.append(name)
            print(f"FAIL  {name}")
            for p in problems:
                print(f"        {p}")
        elif VERBOSE:
            print(f"ok    {name}")
        return fn

    return deco


def read(name: str) -> str:
    return (ROOT / name).read_text()


def body_of(paper: str) -> str:
    """Paper text excluding its appendix claim map."""
    return read(paper).split("## Appendix A")[0]


def headings(text: str) -> set[str]:
    return {m.group(1) for m in re.finditer(r"^#{2,4}\s+([0-9]+(?:\.[0-9]+)?)\.?\s", text, re.M)}


def claim_ids() -> tuple[set[str], set[str]]:
    """(live, retracted) claim ids from claims/, or from CLAIMS.md pre-migration."""
    cdir = ROOT / "claims"
    if cdir.exists():
        live, dead = set(), set()
        for p in cdir.rglob("*.md"):
            if p.name == "INDEX.md":
                continue
            m = re.search(r"^id:\s*(\S+)", p.read_text(), re.M)
            if not m:
                continue
            (dead if "retracted" in p.parts else live).add(m.group(1))
        return live, dead
    text = read("CLAIMS.md")
    head, _, tail = text.partition("## Retracted / dead")
    grab = lambda s: set(re.findall(r"^\|\s*\*?\*?([CFH]\d+)\*?\*?\s*\|", s, re.M))
    return grab(head), grab(tail.split("## F-number")[0])


# ---------------------------------------------------------------- checks --


@check("1. every claim id cited in a paper appendix exists")
def _():
    live, dead = claim_ids()
    known = live | dead
    bad = []
    for paper in PAPERS:
        app = read(paper).split("## Appendix A")
        if len(app) < 2:
            bad.append(f"{paper}: no Appendix A claim map")
            continue
        for cid in sorted(set(re.findall(r"\b[CF]\d+\b", app[1]))):
            if cid not in known:
                bad.append(f"{paper}: cites {cid}, which is in no claim file")
    return bad


@check("2. no paper cites a retracted claim as support")
def _():
    live, dead = claim_ids()
    bad = []
    for paper in PAPERS:
        parts = read(paper).split("## Appendix A")
        if len(parts) < 2:
            continue
        for line in parts[1].split("\n"):
            if not line.startswith("|"):
                continue
            # a row naming the retractions section may legitimately cite dead ids
            if "etract" in line:
                continue
            for cid in re.findall(r"\b[CF]\d+\b", line):
                if cid in dead:
                    bad.append(f"{paper}: {cid} is retracted but cited in {line.split('|')[1].strip()!r}")
    return bad


@check("3. appendix claim-map rows match real section headings")
def _():
    bad = []
    for paper in PAPERS:
        parts = read(paper).split("## Appendix A")
        if len(parts) < 2:
            continue
        heads = headings(parts[0])
        for m in re.finditer(r"^\|\s*([0-9]+(?:\.[0-9]+)?)\.?\s", parts[1], re.M):
            if m.group(1) not in heads:
                bad.append(f"{paper}: appendix row for §{m.group(1)}, which is not a heading")
    return bad


@check("4. every section cross-reference resolves")
def _():
    bad = []
    for paper in PAPERS:
        text = read(paper)
        heads = headings(text) | {"A"}
        refs = {m.group(1) for m in re.finditer(r"§\s*([0-9]+(?:\.[0-9]+)?|[A-Z])", text)}
        for r in sorted(refs - heads - EXTERNAL_SECTION_REFS):
            bad.append(f"{paper}: §{r} does not resolve to a heading")
    return bad


@check("5. no lab notation in paper prose")
def _():
    bad = []
    for paper in PAPERS:
        body = body_of(paper)
        for m in re.finditer(r"§[A-Z][A-Z0-9]*\b", body):
            bad.append(f"{paper}: internal notes reference {m.group(0)!r}")
        for cid in re.findall(r"\b[CF]\d+\b", body):
            if cid != "C0":
                bad.append(f"{paper}: claim id {cid!r} in prose (belongs in Appendix A)")
        for m in re.finditer(r"\(TODO [0-9]", body):
            bad.append(f"{paper}: TODO reference {m.group(0)!r} in prose")
    return bad


@check("6. numbered statements are unique and sequential")
def _():
    # A *declaration* opens a bold span and is followed immediately by "." or a
    # parenthetical: "**Corollary 3.**", "**Lemma 1 (Paper A, Theorem 1).**".
    # Bold prose that merely names one ("**Corollary 3 is in fact stronger...**")
    # has words after the number and must not count.
    decl = re.compile(r"\*\*(Theorem|Proposition|Corollary|Lemma) (\d+)(?=\.|\s*\()")
    bad = []
    for paper in PAPERS:
        found = decl.findall(read(paper))
        nums = [int(n) for _, n in found]
        seen: dict[int, str] = {}
        for kind, n in found:
            n = int(n)
            if n in seen:
                bad.append(f"{paper}: {kind} {n} collides with {seen[n]} {n}")
            seen[n] = kind
        if nums and sorted(nums) != list(range(min(nums), min(nums) + len(nums))):
            bad.append(f"{paper}: statement numbers not contiguous: {nums}")
    return bad


@check("7. reproducibility commands name files that exist")
def _():
    bad = []
    for paper in PAPERS:
        text = read(paper)
        mods = re.findall(r"experiments\.(\w+)", text)
        if not mods:
            bad.append(f"{paper}: no reproducibility commands at all")
        for m in mods:
            if not (ROOT / "experiments" / f"{m}.py").exists():
                bad.append(f"{paper}: experiments/{m}.py does not exist")
        for t in re.findall(r"uv run python (test_\w+)\.py", text):
            if not (ROOT / f"{t}.py").exists():
                bad.append(f"{paper}: {t}.py does not exist")
    return bad


@check("8. generated indexes are up to date")
def _():
    r = subprocess.run([sys.executable, str(ROOT / "tools" / "reindex.py"), "--check"],
                       capture_output=True, text=True)
    return [] if r.returncode == 0 else [l for l in r.stdout.strip().split("\n") if l]


@check("9. prose suite count matches ls test_*.py")
def _():
    n = len(list(ROOT.glob("test_*.py")))
    words = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six",
             7: "seven", 8: "eight", 9: "nine", 10: "ten", 11: "eleven", 12: "twelve"}
    want = words.get(n, str(n))
    # Only count number-words qualifying "suite"; "regression suite" is not a count.
    pat = re.compile(r"\b(" + "|".join(words.values()) + r"|\d+)[- ]suite\b", re.I)
    bad = []
    for name in PAPERS + ["README.md", "CLAUDE.md"]:
        if not (ROOT / name).exists():
            continue
        for m in pat.finditer(read(name)):
            if m.group(1).lower() != want:
                bad.append(f"{name}: says {m.group(0)!r} but there are {n} suites ({want})")
    return bad


@check("10. a number in two documents must have a home in claims/")
def _():
    """The direct answer to the duplication problem: 0.716 lived in seven files.
    A figure repeated across documents must also appear in a claim file, so
    there is one authoritative copy to correct."""
    cdir = ROOT / "claims"
    if not cdir.exists():
        return []  # pre-migration; nothing to check against
    claim_text = "".join(p.read_text() for p in cdir.rglob("*.md") if p.name != "INDEX.md")
    claim_text = re.sub(r"(?<=\d),(?=\d\d\d)", "", claim_text)

    docs = PAPERS + ["HANDOFF.md", "README.md"]
    seen: dict[str, set[str]] = {}
    for name in docs:
        if not (ROOT / name).exists():
            continue
        text = re.sub(r"(?<=\d),(?=\d\d\d)", "", read(name))
        text = re.sub(r"```.*?```", "", text, flags=re.S)          # code blocks
        text = re.sub(r"(?:arXiv:|quant-ph/)[\d.]+", "", text)     # citations
        for m in re.finditer(r"(?<![\w.:/-])(\d{4,}|0\.\d{3,})(?![\w.:/-])", text):
            seen.setdefault(m.group(1), set()).add(name)

    bad = []
    for num, files in sorted(seen.items()):
        if len(files) >= 2 and num not in claim_text:
            bad.append(f"{num} appears in {sorted(files)} but in no claim file")
    return bad


# ------------------------------------------------------------------ main --

if __name__ == "__main__":
    print(f"tools/check.py -- {checks_run} checks")
    if failures:
        print(f"\n{len(failures)} FAILED: {', '.join(failures)}")
        sys.exit(1)
    print(f"all {checks_run} checks pass")
