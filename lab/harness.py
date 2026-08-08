"""The working protocol, made executable.

Encodes the rules from METHOD.md / the computational-research skill that have
each already caught a real error in this project:

  * predictions are DECLARED BEFORE MEASURING (derive-then-test) -- declaring
    one after a check has run earns a warning in the final report;
  * every experiment carries a MUST-FAIL control -- a place the effect must
    vanish; its absence earns a warning (the step-9 vacuous pass survived
    exactly until its must-fail control failed to fail);
  * every prediction must be resolved by a check -- unresolved predictions
    are reported, so quiet dropping is impossible;
  * finish() exits nonzero on any violated check, so a failing experiment
    cannot look like a passing one in a log.

Usage:

    exp = Experiment("myexp", doc=__doc__)
    exp.predict("P1", "support is constant across n_exp")
    exp.must_fail("C1", "the beta>1 control must grow")
    ...measure...
    exp.check("P1", sp == base, f"{sp} vs {base}")
    exp.fail_check("C1", grew, f"{sp2} vs {base2}")
    exp.finish()
"""
from __future__ import annotations
import sys


class Experiment:
    def __init__(self, name: str, doc: str | None = None,
                 exit_on_fail: bool = True):
        self.name = name
        self.exit_on_fail = exit_on_fail
        self._predictions: dict[str, str] = {}
        self._controls: dict[str, str] = {}
        self._results: list[tuple[str, str, bool, str]] = []
        self._resolved: set[str] = set()
        self._measured = False
        self.warnings: list[str] = []
        print(f"=== {name} ===")
        if doc:
            head = doc.strip().splitlines()[0]
            print(f"    {head}")

    # -- declaration (before measurement) -----------------------------------
    def predict(self, pid: str, text: str):
        if self._measured:
            self.warnings.append(
                f"prediction {pid!r} declared AFTER measurement began "
                f"(derive-then-test violated)")
        self._predictions[pid] = text
        print(f"  PREDICT {pid}: {text}")

    def must_fail(self, pid: str, text: str = ""):
        """Register a control that must fail (the effect must vanish here)."""
        self._controls[pid] = text
        print(f"  CONTROL {pid} (must fail): {text}")

    # -- reporting ----------------------------------------------------------
    def section(self, title: str):
        print(f"\n-- {title}")

    def log(self, *args):
        print(" ", *args)

    # -- resolution ---------------------------------------------------------
    def check(self, pid: str, ok: bool, detail: str = ""):
        self._measured = True
        self._resolved.add(pid)
        self._results.append((pid, "check", bool(ok), detail))
        print(f"  {'PASS' if ok else 'FAIL'}  {pid}  {detail}")
        return bool(ok)

    def fail_check(self, pid: str, failed: bool, detail: str = ""):
        """Resolve a must-fail control: `failed` must be True (it DID fail)."""
        self._measured = True
        self._resolved.add(pid)
        self._results.append((pid, "must-fail", bool(failed), detail))
        tag = "PASS (control failed as required)" if failed \
            else "FAIL (CONTROL DID NOT FAIL -- test may be vacuous)"
        print(f"  {tag}  {pid}  {detail}")
        return bool(failed)

    # -- verdict ------------------------------------------------------------
    def finish(self) -> bool:
        for pid in self._predictions:
            if pid not in self._resolved:
                self.warnings.append(f"prediction {pid!r} was never checked")
        for pid in self._controls:
            if pid not in self._resolved:
                self.warnings.append(f"must-fail control {pid!r} was never "
                                     f"resolved with fail_check()")
        if not self._controls:
            self.warnings.append("no must-fail control was registered "
                                 "(a test only every case passes proves nothing)")

        bad = [r for r in self._results if not r[2]]
        print(f"\n=== {self.name}: "
              f"{len(self._results) - len(bad)}/{len(self._results)} checks pass"
              + (f", {len(bad)} FAILED" if bad else "") + " ===")
        for w in self.warnings:
            print(f"  WARNING: {w}")
        ok = not bad
        if not ok and self.exit_on_fail:
            sys.exit(1)
        return ok
