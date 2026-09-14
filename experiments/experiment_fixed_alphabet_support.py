"""Bounded exact support audit for a fixed reflection-label alphabet.

Question answered: does a fixed alphabet of d route labels give a polynomial
finite-time backward support cover, even when its static subgroup is all of
Z_M?  This is an integer support experiment only; it does not propagate fine
work amplitudes or claim a simulation speedup.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  The d=3 word q=(0,1,101) repeated cyclically through k=32 obeys the
      independent normal-form L1 cover at every reverse intermediate prefix.
  P2  The sharper d<=2 bounds min(M,2^ell,1+2*ell), and the d=1 bound,
      survive collisions, fixed points, even M and noninvertible differences.
  P3  Exhaustive selected-history enumeration agrees with the recurrence for
      the small k<=8 rows, while the growing-alphabet pair word has support
      beyond both a universal two-label and a fixed-quadratic bound.

  C1  A universal two-label bound applied to the growing alphabet must fail.
  C2  A fixed-quadratic (d<=3) bound applied to the growing alphabet must
      fail at a finite, nonvacuous row.
  C3  A single gamma=0 support is not a global support envelope when another
      sector attains the disjoint E/O envelope.

The large-k rows use only set recurrence and coefficient covers; no 2^k path
enumeration is performed there.  This is a known generalized-dihedral
normal-form check, not a novelty or hardness claim.
"""
from __future__ import annotations

import itertools
import math
import time
from datetime import datetime, timezone
from pathlib import Path

from lab import Experiment
from lab.coherent_reverse import route_support_bound as production_route_support_bound


MAX_SET_ENTRIES = 65_536
MAX_OPS = 2_000_000
WORD_PROBE_MAX_OPS = 100_000
FIXED_M = 1_000_003
FIXED_ALPHABET = (0, 1, 101)
MAX_K = 32
SMALL_EXHAUSTIVE_K = 8

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__, exit_on_fail=False)
exp.predict("P1", "fixed d=3 cyclic rows stay inside the independent L1 normal-form cover")
exp.predict("P2", "d<=2 sharper bounds survive wrap, collisions, fixed points and noninvertible differences")
exp.predict("P3", "small exhaustive histories agree and growing alphabets exceed the stated controls")
exp.predict("P4", "production support bounds match the independent tight L1 bound at every included suffix")
exp.predict("P5", "independent E/O envelopes match every reverse prefix and product-condition rows are globally tight")
exp.must_fail("C1", "the universal two-label bound must fail for a growing alphabet")
exp.must_fail("C2", "the fixed-quadratic d<=3 bound must fail for a growing alphabet")
exp.must_fail("C3", "gamma=0 is not a global support envelope for the cyclic 0,1,2 word")


class Budget:
    """Guarded candidate/update-operation count, not a CPU instruction count."""

    def __init__(self, limit=MAX_OPS):
        self.limit = int(limit)
        self.operations = 0

    def add(self, amount=1):
        self.operations += int(amount)
        if self.operations > self.limit:
            raise MemoryError(f"integer support operation budget exceeded: {self.operations}")


def report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"fixed_alphabet_support_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    serial = 0
    while path.exists():
        serial += 1
        path = Path("out") / f"fixed_alphabet_support_{stamp}_{serial}.json"
    return path


def guard_entries(entries, label):
    if type(entries) is not int or entries < 0 or entries > MAX_SET_ENTRIES:
        raise MemoryError(f"{label} integer-set allocation exceeds {MAX_SET_ENTRIES} entries")


def reverse_prefixes(modulus, gamma, chronological, budget):
    """Return actual reverse supports, using no history enumeration."""
    reached = {int(gamma) % int(modulus)}
    guard_entries(len(reached), "initial support")
    prefixes = [reached]
    retained_entries = len(reached)
    for q in reversed(tuple(int(q) % modulus for q in chronological)):
        budget.add(len(reached))
        guard_entries(retained_entries + min(modulus, 2*len(reached)),
                      "retained prefix supports before update")
        nxt = set(reached)
        guard_entries(len(nxt) + len(reached), "support update")
        for alpha in reached:
            nxt.add((-alpha - q) % modulus)
        guard_entries(len(nxt), "support")
        reached = nxt
        prefixes.append(reached)
        retained_entries += len(reached)
        guard_entries(retained_entries, "retained prefix supports")
    return prefixes


def eo_prefixes(modulus, chronological, budget):
    """Independent E/O normal form, retaining every reverse prefix.

    The updates use the old E and O simultaneously.  This is deliberately a
    second recurrence rather than a transformation of ``reverse_prefixes``.
    """
    modulus = int(modulus)
    even = {0}
    odd = set()
    rows = [(even, odd)]
    retained_entries = 1
    for q in reversed(tuple(int(q) % modulus for q in chronological)):
        old_even, old_odd = even, odd
        candidate_entries = 2 * (len(old_even) + len(old_odd))
        guard_entries(candidate_entries, "E/O update before growth")
        guard_entries(retained_entries+candidate_entries, "retained E/O prefixes before growth")
        budget.add(len(old_even) + len(old_odd))
        even = set(old_even)
        even.update((-value - q) % modulus for value in old_odd)
        odd = set(old_odd)
        odd.update((-value - q) % modulus for value in old_even)
        guard_entries(len(even) + len(odd), "E/O prefix")
        retained_entries += len(even)+len(odd)
        rows.append((even, odd))
    return rows


def eo_reached(modulus, gamma, pair):
    """Reachable labels represented by one E/O prefix."""
    even, odd = pair
    return ({(int(gamma) + value) % modulus for value in even}
            | {(-int(gamma) + value) % modulus for value in odd})


def eo_global_stats(modulus, pair, budget, compute_witness=True):
    """Return the E/O global envelope and a product-condition witness."""
    even, odd = pair
    total = len(even) + len(odd)
    gcd_factor = math.gcd(2, int(modulus))
    product = len(even) * len(odd)
    condition = product < int(modulus) // gcd_factor
    witness = None
    difference_size = None
    if compute_witness and condition:
        # The product condition guarantees a missing doubled-sector value;
        # guard the Cartesian materialization before allocating it.
        guard_entries(product, "E/O difference set before growth")
        budget.add(product)
        differences = {(o - e) % modulus for o in odd for e in even}
        difference_size = len(differences)
        guard_entries(difference_size, "E/O difference set")
        for gamma in range(min(int(modulus) // gcd_factor,
                              difference_size + 1)):
            budget.add(1)
            if (2 * gamma) % modulus not in differences:
                witness = gamma
                break
    return {
        "even_size": len(even), "odd_size": len(odd),
        "envelope": min(int(modulus), total), "product": product,
        "gcd_2_M": gcd_factor, "product_condition": condition,
        "difference_size": difference_size, "tight_witness_gamma": witness,
        "tight_by_product_condition": bool(condition and witness is not None),
    }


def word_support_audit(name, modulus, chronological, gammas, budget,
                       compute_global=True):
    """Compare E/O and actual supports at every reverse prefix."""
    chronological = tuple(int(q) % int(modulus) for q in chronological)
    eo_rows = eo_prefixes(modulus, chronological, budget)
    actual_rows = {
        int(gamma): reverse_prefixes(modulus, gamma, chronological, budget)
        for gamma in gammas
    }
    prefix_mismatches = 0
    max_actual = 0
    prefix_rows = []
    for depth, pair in enumerate(eo_rows):
        expected_sizes = {}
        actual_sizes = {}
        for gamma, actual in actual_rows.items():
            expected = eo_reached(modulus, gamma, pair)
            reached = actual[depth]
            prefix_mismatches += expected != reached
            expected_sizes[str(gamma)] = len(expected)
            actual_sizes[str(gamma)] = len(reached)
            max_actual = max(max_actual, len(reached))
        prefix_rows.append({"depth": depth, "E_size": len(pair[0]),
                            "O_size": len(pair[1]),
                            "expected_sizes": expected_sizes,
                            "actual_sizes": actual_sizes})
    final_stats = eo_global_stats(modulus, eo_rows[-1], budget,
                                  compute_witness=compute_global)
    witness = final_stats["tight_witness_gamma"]
    final_stats["witness_actual_support"] = None
    if witness is not None:
        actual = reverse_prefixes(modulus,witness,chronological,budget)[-1]
        final_stats["witness_actual_support"] = len(actual)
        prefix_mismatches += actual != eo_reached(modulus,witness,eo_rows[-1])
        final_stats["tight_by_product_condition"] &= len(actual)==final_stats["envelope"]
    final_stats.update({"name": name, "modulus": modulus,
                        "chronological": chronological,
                        "gammas": tuple(gammas), "prefix_mismatches": prefix_mismatches,
                        "max_actual_support": max_actual,
                        "prefix_rows": prefix_rows})
    return final_stats


def selected_history_support(modulus, gamma, chronological, budget):
    """Independent small-k support by enumerating selected reflection subsets."""
    word = tuple(int(q) % modulus for q in chronological)
    if len(word) > SMALL_EXHAUSTIVE_K:
        raise ValueError("selected-history reference is capped at k<=8")
    reached = set()
    for mask in range(1 << len(word)):
        budget.add(len(word) + 1)
        alpha = int(gamma) % modulus
        for index in range(len(word) - 1, -1, -1):
            if mask & (1 << index):
                alpha = (-alpha - word[index]) % modulus
        reached.add(alpha)
    guard_entries(len(reached), "selected-history support")
    return reached


def fixed_point_roots(modulus, q, budget):
    """Solve 2*alpha == -q (mod modulus), without scanning the modulus."""
    modulus, q = int(modulus), int(q) % int(modulus)
    g = math.gcd(2, modulus)
    rhs = -q
    budget.add(4)
    if rhs % g:
        return []
    reduced_modulus = modulus // g
    reduced_rhs = (rhs // g) % reduced_modulus
    inverse = pow(2 // g, -1, reduced_modulus)
    base = (reduced_rhs * inverse) % reduced_modulus
    roots = [(base + j * reduced_modulus) % modulus for j in range(g)]
    budget.add(g)
    return roots


def normal_form_cover(modulus, gamma, alphabet, depth, budget):
    """Independent L1 coefficient cover for all words of length <= depth."""
    alphabet = tuple(dict.fromkeys(int(q) % modulus for q in alphabet))
    if not alphabet:
        return {int(gamma) % modulus}
    q0 = alphabet[0]
    deltas = tuple(q - q0 for q in alphabet[1:])
    if depth == 0:
        return {int(gamma) % modulus}
    cover = set()
    # This is a coefficient cover, not a word/path enumeration.  The number
    # of vectors is at most (2*depth+1)^(d-1), guarded before insertion.
    coefficient_count = (2 * depth + 1) ** len(deltas)
    guard_entries(2 * coefficient_count, "normal-form coefficient cover")
    for coeffs in itertools.product(range(-depth, depth + 1), repeat=len(deltas)):
        budget.add(1)
        if sum(abs(n) for n in coeffs) > depth:
            continue
        shift = sum(n * delta for n, delta in zip(coeffs, deltas))
        cover.add((int(gamma) + shift) % modulus)
        cover.add((-int(gamma) - q0 + shift) % modulus)
    guard_entries(len(cover), "normal-form cover")
    return cover


def sharp_bound(modulus, depth, alphabet):
    d = len(set(int(q) % modulus for q in alphabet))
    if d == 0:
        return 1
    if d == 1:
        return min(modulus, 1 << depth, 1 if depth == 0 else 2)
    if d == 2:
        return min(modulus, 1 << depth, 1 + 2 * depth)
    return min(modulus, 1 << depth, 2 * (2 * depth + 1) ** (d - 1))


def l1_ball_count(dimension, radius):
    """Number of integer vectors in an L1 ball of the given radius."""
    return sum(2 ** j * math.comb(dimension, j) * math.comb(radius, j)
               for j in range(min(dimension, radius) + 1))


def tight_support_bound(modulus, depth, alphabet):
    """Independent bound, retaining the established sharper d<=2 cases."""
    d = len(set(int(q) % modulus for q in alphabet))
    if d == 0:
        return 1
    if d == 1:
        return min(modulus, 1 << depth, 1 if depth == 0 else 2)
    if d == 2:
        return min(modulus, 1 << depth, 1 + 2 * depth)
    return min(modulus, 1 << depth, 2 * l1_ball_count(d - 1, depth))


def cyclic_word(alphabet, length):
    alphabet = tuple(alphabet)
    return tuple(alphabet[i % len(alphabet)] for i in range(length))


def fixed_d3_audit(budget):
    rows = []
    gammas = (0, 1, FIXED_M - 1, FIXED_M // 2)
    for k in range(MAX_K + 1):
        word = cyclic_word(FIXED_ALPHABET, k)
        max_leaf = max_peak = 0
        normal_violations = exhaustive_mismatches = 0
        production_violations = production_mismatches = 0
        for gamma in gammas:
            prefixes = reverse_prefixes(FIXED_M, gamma, word, budget)
            max_leaf = max(max_leaf, len(prefixes[-1]))
            for depth, reached in enumerate(prefixes):
                max_peak = max(max_peak, len(reached))
                cover = normal_form_cover(FIXED_M, gamma, word[:k], depth, budget)
                if not reached <= cover:
                    normal_violations += 1
                suffix = word[k - depth:] if depth else ()
                expected = tight_support_bound(FIXED_M, depth, suffix)
                actual_production = production_route_support_bound(FIXED_M, suffix)
                production_mismatches += actual_production != expected
                production_violations += len(reached) > actual_production
            if k <= SMALL_EXHAUSTIVE_K:
                exact = selected_history_support(FIXED_M, gamma, word, budget)
                if exact != prefixes[-1]:
                    exhaustive_mismatches += 1
        rows.append({
            "k": k, "word": word, "gamma_count": len(gammas),
            "max_leaf_support": max_leaf, "max_peak_support": max_peak,
            "normal_form_bound": sharp_bound(FIXED_M, k, FIXED_ALPHABET),
            "tight_l1_bound": tight_support_bound(FIXED_M, k, FIXED_ALPHABET),
            "normal_form_violations": normal_violations,
            "exhaustive_selected_history_mismatches": exhaustive_mismatches,
            "production_bound_violations": production_violations,
            "production_bound_mismatches": production_mismatches,
        })
    return rows


def two_label_audit(budget):
    specs = [
        ("odd_wrap_pair01", 1_000_003, (0, 1), (0, 1, 1_000_002)),
        ("even_noninvertible_pair26", 16, (2, 6), (0, 1, 7, 8, 15)),
        ("even_noninvertible_pair04", 8, (0, 4), (0, 1, 4, 7)),
        ("coincident_fixed_pair22", 16, (2, 2), (0, 7, 15)),
        ("coincident_pair00", 7, (0, 0), (0, 1, 6)),
        ("even_d3_noninvertible048", 12, (0, 4, 8), (0, 1, 6, 11)),
        ("odd_d3_noninvertible036", 15, (0, 3, 6), (0, 1, 7, 14)),
    ]
    rows = []
    for name, modulus, alphabet, gammas in specs:
        max_violation = 0
        production_violations = production_mismatches = 0
        max_peak = 0
        fixed_points = fixed_point_roots(modulus, alphabet[0], budget)
        for k in range(MAX_K + 1):
            word = cyclic_word(alphabet, k)
            for gamma in gammas:
                prefixes = reverse_prefixes(modulus, gamma, word, budget)
                sizes = [len(support) for support in prefixes]
                max_peak = max(max_peak, max(sizes))
                for depth, size in enumerate(sizes):
                    max_violation = max(max_violation,
                                        size - sharp_bound(modulus, depth, alphabet))
                    suffix = word[k - depth:] if depth else ()
                    expected = tight_support_bound(modulus, depth, suffix)
                    actual_production = production_route_support_bound(modulus, suffix)
                    production_mismatches += actual_production != expected
                    production_violations += size > actual_production
        rows.append({"name": name, "modulus": modulus, "alphabet": alphabet,
                     "gammas": gammas, "fixed_points": fixed_points,
                     "max_bound_excess": max_violation, "max_peak_support": max_peak,
                     "production_bound_violations": production_violations,
                     "production_bound_mismatches": production_mismatches})
    return rows


def growing_alphabet_audit(budget):
    rows = []
    for m in range(1, 9):
        word = tuple(value for j in range(m) for value in (0, 3 ** j))
        modulus = 3 ** (m + 1)
        prefixes = reverse_prefixes(modulus, 0, word, budget)
        support = len(prefixes[-1])
        two_label = min(modulus, 1 << len(word), 1 + 2 * len(word))
        quadratic = min(modulus, 1 << len(word), 2 * (2 * len(word) + 1) ** 2)
        exhaustive_match = None
        if len(word) <= SMALL_EXHAUSTIVE_K:
            exhaustive_match = selected_history_support(modulus, 0, word, budget) == prefixes[-1]
        production_mismatches = production_violations = 0
        for depth, reached in enumerate(prefixes):
            suffix = word[len(word) - depth:] if depth else ()
            expected = tight_support_bound(modulus, depth, suffix)
            actual_production = production_route_support_bound(modulus, suffix)
            production_mismatches += actual_production != expected
            production_violations += len(reached) > actual_production
        rows.append({"m": m, "k": len(word), "alphabet_size": len(set(word)),
                     "modulus": modulus, "support": support,
                     "two_label_bound": two_label,
                     "fixed_quadratic_bound": quadratic,
                     "exceeds_two_label": support > two_label,
                     "exceeds_fixed_quadratic": support > quadratic,
                     "exhaustive_selected_history_match": exhaustive_match,
                     "production_bound_violations": production_violations,
                     "production_bound_mismatches": production_mismatches})
    return rows


def word_probe_audit(budget):
    """Bounded independent audit of the word-specific E/O envelope."""
    rows = []
    fixed_word = cyclic_word(FIXED_ALPHABET, MAX_K)
    # This long fixed-alphabet row checks every reverse prefix without
    # repeatedly retaining all k<=32 E/O histories.
    rows.append(word_support_audit(
        "fixed_d3_k32", FIXED_M, fixed_word,
        (0, 1, FIXED_M - 1), budget, compute_global=False))

    separated_modulus = (1 << 40) - 1
    separated_alphabet = (0, 679_535_556_937, 314_159_265_359)
    separated_rows = []
    for k in range(17):
        row = word_support_audit(
            f"separated_k{k}", separated_modulus,
            cyclic_word(separated_alphabet, k), (0, 1), budget,
            compute_global=True)
        separated_rows.append(row)
    rows.extend(separated_rows)

    counterexample = word_support_audit(
        "M101_cyclic012_k8", 101, cyclic_word((0, 1, 2), 8),
        (0, 5), budget, compute_global=True)
    gamma0_size = counterexample["prefix_rows"][-1]["actual_sizes"]["0"]
    gamma5_size = counterexample["prefix_rows"][-1]["actual_sizes"]["5"]
    counterexample["gamma0_size"] = gamma0_size
    counterexample["gamma5_size"] = gamma5_size
    counterexample["gamma0_is_global_envelope"] = (
        gamma0_size >= counterexample["envelope"])
    rows.append(counterexample)
    rows.append(word_support_audit("even_M1000_cyclic012_k8",1000,
        cyclic_word((0,1,2),8),(0,5),budget,compute_global=True))

    product_rows = [row for row in rows if row["product_condition"]
                    and row["difference_size"] is not None]
    p5 = (len(product_rows)==len(separated_rows)+2
          and all(row["prefix_mismatches"] == 0 for row in rows)
          and all(row["tight_by_product_condition"] for row in product_rows))
    c3 = (not counterexample["gamma0_is_global_envelope"]
          and gamma0_size < counterexample["envelope"]
          and gamma5_size == counterexample["envelope"])
    return {"rows": rows, "counterexample": counterexample,
            "product_rows_checked": len(product_rows),
            "p5": p5, "c3": c3}


def json_safe(value):
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [json_safe(v) for v in value]
    return value


def main():
    started = time.perf_counter()
    report = {"status": "PASS", "fixed_d3": [], "two_label": [],
              "growing_alphabet": [], "word_probe": {},
              "operations": 0, "word_probe_operations": 0}
    p1 = p2 = p3 = p4 = p5 = True
    try:
        budget = Budget()
        report["fixed_d3"] = fixed_d3_audit(budget)
        report["two_label"] = two_label_audit(budget)
        report["growing_alphabet"] = growing_alphabet_audit(budget)
        p1 = all(row["normal_form_violations"] == 0
                 and row["exhaustive_selected_history_mismatches"] == 0
                 for row in report["fixed_d3"])
        p2 = all(row["max_bound_excess"] <= 0 for row in report["two_label"])
        p3 = all(row["exhaustive_selected_history_match"] is not False
                 for row in report["growing_alphabet"])
        all_rows = (report["fixed_d3"] + report["two_label"]
                    + report["growing_alphabet"])
        p4 = all(row["production_bound_violations"] == 0
                 and row["production_bound_mismatches"] == 0
                 for row in all_rows)
        c1 = any(row["exceeds_two_label"] for row in report["growing_alphabet"])
        c2 = any(row["exceeds_fixed_quadratic"] for row in report["growing_alphabet"])
        report["operations"] = budget.operations
        word_budget = Budget(WORD_PROBE_MAX_OPS)
        report["word_probe"] = word_probe_audit(word_budget)
        report["word_probe_operations"] = word_budget.operations
        p5 = report["word_probe"]["p5"]
        c3 = report["word_probe"]["c3"]
        report["status"] = ("PASS" if p1 and p2 and p3 and p4 and p5
                             and c1 and c2 and c3 else "FAIL")
        exp.check("P1", p1, "d=3 cyclic prefixes lie in the independent L1 cover")
        exp.check("P2", p2, "d<=2 sharp bounds cover all edge cases")
        exp.check("P3", p3, "small exact selected-history checks agree")
        exp.check("P4", p4, "production suffix bounds match independent L1 bounds")
        exp.check("P5", p5, "E/O word envelopes match all audited reverse prefixes and tightness witnesses")
        exp.fail_check("C1", c1, "growing alphabet exceeds universal two-label bound")
        exp.fail_check("C2", c2, "growing alphabet exceeds fixed quadratic bound")
        exp.fail_check("C3", c3, "gamma=0 fails as a global envelope in the M=101 control")
    except Exception as exc:
        report["status"] = "FAIL"
        report["exception"] = repr(exc)
        report["traceback"] = __import__("traceback").format_exc()
        exp.log("EXCEPTION", repr(exc))
        for name in ("P1", "P2", "P3", "P4", "P5"):
            exp.check(name, False, "exception before completion")
        for name in ("C1", "C2", "C3"):
            exp.fail_check(name, False, "exception before control")
    report["elapsed_seconds"] = time.perf_counter() - started
    path = report_path()
    ok = exp.finish(report_path=path, rows=[json_safe(report)],
                    metadata={"max_k": MAX_K, "fixed_alphabet": FIXED_ALPHABET,
                              "fixed_modulus": FIXED_M,
                              "max_set_entries": MAX_SET_ENTRIES,
                              "max_operations": MAX_OPS,
                              "word_probe_max_operations": WORD_PROBE_MAX_OPS,
                              "operation_scope": "enumerated candidate/update operations, not CPU instructions",
                              "detail_report_is_first_row": True})
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
