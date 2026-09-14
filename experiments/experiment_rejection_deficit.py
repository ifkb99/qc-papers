"""Exact Fraction audit of conservative rejection acceptance intervals.

For valid endpoints PL<=P<=PU and DL<=D<=DU, the conservative dyadic
acceptance is

    a = floor(2^L * max(0,PL) / DU) / 2^L,  if DU>0, else 0.

The finite diagnostic checks 0<=a<=1 and

  0 <= P-D*a <= (PU-PL)+(DU-DL)+D/2^L,

including P=D=0, cancellation P=0<D, and D as small as 2^-1024.  A fixed
two-output accepted-law fixture then checks the C60 normalized-law bound while
sweeping only L.  All arithmetic remains exact Fraction arithmetic.

P1: the conservative floor acceptance and residual bound hold on fixed edge
    cases, including tiny denominators and exact cancellation.
P2: the normalized accepted fixture law obeys the stated C60 TV bound for a
    fixed q-perturbation while only the random-word precision L varies.
C1: exact ideal acceptance ratios do not certify a wrong proposal law.
C2: silently deleting accepted mass changes the normalized output law.

This is a mathematical acceptance/normalization diagnostic, not a production
sampler or a numerical floating-point certificate.  No array or RSS estimate
is used; only bounded operation counts and exact rational outputs are logged.
"""
from __future__ import annotations

import json
import time
import traceback
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

from lab import Experiment


MAX_J = 1024
MAX_L = 16


def report_path(prefix: str = "rejection_deficit") -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"{prefix}_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def conservative_acceptance(PL: Fraction, PU: Fraction, DL: Fraction,
                            DU: Fraction, bits: int) -> Fraction:
    if (not all(isinstance(x, Fraction) for x in (PL, PU, DL, DU))
            or PL > PU or PU < 0 or DL > DU or DU < 0 or PL > DU
            or type(bits) is not int or not 0 <= bits <= MAX_L):
        raise ValueError("invalid exact acceptance endpoints")
    if DU == 0:
        return Fraction(0)
    # The denominator belongs inside the floor: floor(2^L*PL/DU)/2^L.
    scaled = max(Fraction(0), PL) * (1 << bits) / DU
    a = Fraction(scaled.numerator // scaled.denominator, 1 << bits)
    # This is conservative even when the lower endpoint is zero.
    if a < 0 or a > 1:
        raise ArithmeticError("conservative acceptance escaped [0,1]")
    return a


def residual_bound(P: Fraction, D: Fraction, PL: Fraction, PU: Fraction,
                   DL: Fraction, DU: Fraction, bits: int) -> tuple[Fraction, Fraction]:
    a = conservative_acceptance(PL, PU, DL, DU, bits)
    residual = P - D * a
    bound = (PU - PL) + (DU - DL) + D / (1 << bits)
    return residual, bound


def total_variation(left: list[Fraction], right: list[Fraction]) -> Fraction:
    return sum(abs(a - b) for a, b in zip(left, right)) / 2


def main() -> None:
    exp = Experiment("rejection_deficit", doc=__doc__)
    exp.predict("P1", "conservative acceptance residual bound holds on exact edge cases")
    exp.predict("P2", "normalized accepted-law fixture obeys C60 TV bound across L")
    exp.must_fail("C1", "zero acceptance-ratio error alone cannot certify a wrong proposal")
    exp.must_fail("C2", "uncharged accepted-mass deletion changes the normalized law")
    started = time.perf_counter()
    rows: list[dict] = []

    exp.section("P1 endpoint and residual inequalities")
    edge_cases = [
        # P=D=0 and DU=0 branch.
        (Fraction(0), Fraction(0), Fraction(0), Fraction(0)),
        # Exact cancellation P=0<D.
        (Fraction(0), Fraction(0), Fraction(1, 3), Fraction(1, 3)),
        # Tiny denominators through j=1024, with exact target ratio 1/3.
        *[(Fraction(1, 3 * (1 << j)), Fraction(1, 3 * (1 << j)),
           Fraction(1, 1 << j), Fraction(1, 1 << j))
          for j in (0, 1, 16, 128, 512, 1024)],
        # Nonzero endpoint widths and a fractional floor ratio.
        (Fraction(1, 5), Fraction(21, 100), Fraction(2, 5), Fraction(41, 100)),
    ]
    edge_failures = []
    worst_edge_ratio = Fraction(0)
    worst_edge_case: tuple[int, int] | None = None
    # Nonnegative quantities may have negative outward LOWER endpoints.
    # They must not be rejected or midpoint-interpreted as negative targets.
    targets = [((PL+PU)/2, (DL+DU)/2) for PL,PU,DL,DU in edge_cases]
    edge_cases += [(Fraction(-1,100),Fraction(1,100),Fraction(-1,100),Fraction(1,100)),
                   (Fraction(-1,100),Fraction(1,100),Fraction(1,3),Fraction(1,3))]
    targets += [(Fraction(0),Fraction(0)), (Fraction(0),Fraction(1,3))]
    from lab.verified_rejection import conservative_threshold
    for index, ((PL, PU, DL, DU), (P,D)) in enumerate(zip(edge_cases, targets)):
        assert 0 <= P <= D and PL <= P <= PU and DL <= D <= DU
        for bits in (0, 1, 4, 8, 16):
            residual, bound = residual_bound(P, D, PL, PU, DL, DU, bits)
            a = conservative_acceptance(PL, PU, DL, DU, bits)
            ratio = residual / bound if bound else Fraction(0)
            if ratio > worst_edge_ratio:
                worst_edge_ratio = ratio
                worst_edge_case = (index, bits)
            production = Fraction(conservative_threshold(PL,DU,bits),1 << bits)
            if not (Fraction(0) <= a <= 1 and 0 <= residual <= bound and production == a):
                edge_failures.append((index, bits, str(a), str(residual), str(bound)))
    exp.check("P1", not edge_failures,
              f"cases={len(edge_cases)}, precision_points={len(edge_cases)*5}, "
              f"failures={edge_failures}")
    edge_endpoints = [dict(PL=str(PL), PU=str(PU), DL=str(DL), DU=str(DU), P=str(P), D=str(D))
                      for (PL, PU, DL, DU),(P,D) in zip(edge_cases,targets)]
    rows.append(dict(series="edge_cases", case_count=len(edge_cases),
                     precision_points=len(edge_cases)*5,
                     tiny_denominator_max_j=MAX_J,
                     endpoints=edge_endpoints,
                     worst_residual_to_bound=str(worst_edge_ratio),
                     worst_residual_case=worst_edge_case,
                     failures=edge_failures))

    exp.section("P2 fixed accepted-law fixture; sweep only L")
    C = Fraction(2)
    ideal_p = [Fraction(1, 2), Fraction(1, 2)]
    ideal_q = [Fraction(1, 4), Fraction(3, 4)]
    ideal_P = ideal_p[:]  # accepted target; ideal_P <= ideal_D componentwise
    ideal_D = [C * q for q in ideal_q]
    eta = Fraction(1, 128)
    q_tilde = [ideal_q[0] + eta, ideal_q[1] - eta]
    PL = ideal_P
    PU = ideal_P
    # The proposal perturbation is enclosed around the ideal denominator.
    DL = [D - C * eta for D in ideal_D]
    DU = [D + C * eta for D in ideal_D]
    fixture_failures = []
    worst_ratio = Fraction(0)
    for bits in range(2, MAX_L + 1):
        acceptance = [conservative_acceptance(PL[i], PU[i], DL[i], DU[i], bits)
                      for i in range(2)]
        accepted = [q_tilde[i] * acceptance[i] for i in range(2)]
        mass = sum(accepted)
        if mass == 0:
            fixture_failures.append((bits, "zero accepted mass"))
            continue
        normalized = [value / mass for value in accepted]
        tv = total_variation(normalized, ideal_p)
        widths = sum((PU[i] - PL[i]) + (DU[i] - DL[i]) for i in range(2))
        bound = 2 * C * eta + widths + C / (1 << bits)
        ratio = tv / bound if bound else Fraction(0)
        worst_ratio = max(worst_ratio, ratio)
        if tv > bound or sum(normalized) != 1:
            fixture_failures.append((bits, str(tv), str(bound)))
        rows.append(dict(series="fixture", bits=bits,
                         endpoints={"PL": [str(x) for x in PL],
                                    "PU": [str(x) for x in PU],
                                    "DL": [str(x) for x in DL],
                                    "DU": [str(x) for x in DU]},
                         acceptance=[str(x) for x in acceptance],
                         accepted_mass=str(mass),
                         ideal_accepted_mass=str(Fraction(1, 2)),
                         mass_loss=str(Fraction(1, 2) - mass),
                         normalized=[str(x) for x in normalized],
                         tv=str(tv), bound=str(bound), ratio=str(ratio)))
    exp.check("P2", not fixture_failures,
              f"L=2..{MAX_L}, eta={eta}, worst_tv_to_bound={worst_ratio}, "
              f"failures={fixture_failures}")

    exp.section("must-fail proposal and mass controls")
    exact_acceptance = [ideal_P[i] / ideal_D[i] for i in range(2)]
    wrong_proposal = [Fraction(1), Fraction(0)]
    wrong_accepted = [wrong_proposal[i] * exact_acceptance[i] for i in range(2)]
    wrong_mass = sum(wrong_accepted)
    wrong_law = [value / wrong_mass for value in wrong_accepted]
    wrong_tv = total_variation(wrong_law, ideal_p)
    exp.fail_check("C1", wrong_tv > 0,
                   f"exact_acceptance={exact_acceptance}, wrong_proposal={wrong_proposal}, "
                   f"normalized_law={wrong_law}, TV={wrong_tv}")

    deleted = [Fraction(0), ideal_P[1]]
    deleted_mass = sum(deleted)
    deleted_law = [value / deleted_mass for value in deleted]
    deleted_tv = total_variation(deleted_law, ideal_p)
    exp.fail_check("C2", deleted_tv > 0,
                   f"deleted_accepted_mass={deleted}, normalized_law={deleted_law}, "
                   f"TV={deleted_tv}")
    rows.append(dict(series="controls", exact_acceptance=[str(x) for x in exact_acceptance],
                     wrong_proposal=[str(x) for x in wrong_proposal],
                     wrong_tv=str(wrong_tv), deleted_tv=str(deleted_tv)))

    report = report_path()
    exp.finish(report_path=report, rows=rows,
               metadata=dict(max_j=MAX_J, max_precision_bits=MAX_L,
                             fixture_C=str(C), fixture_eta=str(eta),
                             arithmetic="exact Fraction only",
                             allocation_bound="bounded operation counts; no RSS claim",
                             elapsed_seconds=time.perf_counter() - started))
    print(f"report: {report}")


if __name__ == "__main__":
    try:
        main()
    except BaseException as exc:
        failure = report_path("rejection_deficit_failure")
        failure.write_text(json.dumps(
            dict(ok=False, error=repr(exc), traceback=traceback.format_exc()),
            indent=2) + "\n")
        raise
