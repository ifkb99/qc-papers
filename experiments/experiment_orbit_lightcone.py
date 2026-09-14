"""Bound the pre-kick support created by a few periodic orbit mixers.

DERIVED BEFORE MEASUREMENT.  Use an indexed orbit of length r, block size
b=3, and the ascending controlled shifts U^(2^i).  Just before a localized
kick after s controls, the early exponent l is uniform in [0,L), L=2^s.
Each of m previous repeated b-by-b mixers changes an orbit label by q-p in
[-(b-1), b-1].  Therefore every branch started at l is supported in the
circular radius R=m*(b-1) around l mod r, independently of interference and
of the matrices supplied at the different insertion times.

For a supplied kick support S, the pre-kick probability F is bounded by the
fraction of l whose residue lies in the R-expanded circular neighborhood of
S.  The exact finite count is at most

    min(1, d*(2*R+1)*(1/r + 1/L)),  d=|S|.

P1: finite branch products obey the R light-cone support bound for fixed
    dense b=3 blocks, including wrapped supports and several mixer counts.
P2: brute residue coverage agrees with an independent quotient/remainder
    count and never exceeds the stated scalar upper bound.
P3: m=0 reproduces the bare uniform residue count exactly; R>=r/2 gives
    full coverage, with both L<r and L>r represented.

C1: incorrectly setting R=0 must fail on a fixed nonzero-mixer case.
C2: using the bare unexpanded residue count must underbound actual F on that
    same fixed case.  If either control does not fail, this run exits failed.

This is a bounded support diagnostic, not a generic propagator or a claim
about post-kick output sampling.  The branch state is formed only by small
finite matrix products (r<=18, t<=7); all dense allocations are guarded at
16 MiB before allocation.  The fixed W is reused across every sweep row.

Run: OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with
     'numpy<2.5' python -m experiments.experiment_orbit_lightcone
"""
from __future__ import annotations

import math
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from lab import Experiment


MAX_BYTES = 16 * 1024 * 1024
MAX_R = 18
MAX_T = 7
B = 3
TOL = 2e-11


def guard(shape, dtype, label: str) -> None:
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_BYTES:
        raise MemoryError(f"{label} allocation {payload} bytes exceeds 16 MiB")


def alloc(shape, dtype, label: str) -> np.ndarray:
    guard(shape, dtype, label)
    return np.empty(shape, dtype=dtype)


def report_path() -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    root = Path("out")
    root.mkdir(parents=True, exist_ok=True)
    path = root / f"orbit_lightcone_{stamp}.json"
    serial = 0
    while path.exists():
        serial += 1
        path = root / f"orbit_lightcone_{stamp}_{serial}.json"
    return path


def fixed_dense_block() -> np.ndarray:
    """One fixed dense unitary, held unchanged over every sweep row."""
    # The 3-point Fourier unitary gives all allowed within-block transitions
    # nonzero, making the support diagnostic sensitive to both signs of delta.
    guard((B, B), np.complex128, "fixed b=3 W")
    p = np.arange(B)
    block = np.exp(2j * np.pi * p[:, None] * p[None, :] / B) / np.sqrt(B)
    return block


def periodic_block(r: int, block: np.ndarray) -> np.ndarray:
    if not 1 <= r <= MAX_R or r % B:
        raise ValueError("r must be a positive multiple of b=3")
    if not isinstance(block, np.ndarray) or block.shape != (B, B):
        raise ValueError("invalid block")
    block = np.asarray(block, dtype=complex)
    if not np.all(np.isfinite(block)):
        raise ValueError("nonfinite block")
    if np.max(np.abs(block.conj().T @ block - np.eye(B))) > 2e-12:
        raise ValueError("block must be unitary")
    guard((r, r), np.complex128, "periodic block")
    result = np.zeros((r, r), dtype=complex)
    for start in range(0, r, B):
        result[start:start + B, start:start + B] = block
    return result


def shift_branch(vec: np.ndarray, power: int) -> np.ndarray:
    """Apply the indexed orbit translation U^power."""
    return np.roll(vec, int(power))


def branch_before_kick(r: int, s: int, l: int, mixer_count: int,
                       block: np.ndarray) -> np.ndarray:
    """Finite algebraic branch product for one early exponent l.

    Mixers are inserted after controls 0,...,m-1.  This is enough to test the
    displacement lemma; no generic circuit simulator or post-kick evolution
    is introduced.
    """
    if not 1 <= r <= MAX_R or r % B or not 0 <= s <= MAX_T or not 0 <= l < 1 << s:
        raise ValueError("invalid branch dimensions")
    if not 0 <= mixer_count <= s:
        raise ValueError("mixer count must be at most the number of controls")
    v = periodic_block(r, block)
    state = np.zeros(r, dtype=complex)
    state[0] = 1.0
    for i in range(s):
        if (l >> i) & 1:
            state = shift_branch(state, 1 << i)
        if i < mixer_count:
            state = v @ state
    return state


def circular_interval(center: int, radius: int, r: int) -> set[int]:
    return {(center + delta) % r for delta in range(-radius, radius + 1)}


def expanded_support(S: set[int], radius: int, r: int) -> set[int]:
    expanded: set[int] = set()
    for p in S:
        expanded.update(circular_interval(p, radius, r))
    return expanded


def brute_cover_count(S: set[int], radius: int, r: int, L: int) -> int:
    expanded = expanded_support(S, radius, r)
    return sum((l % r) in expanded for l in range(L))


def quotient_cover_count(S: set[int], radius: int, r: int, L: int) -> int:
    """Independent count using residue multiplicities, not l enumeration."""
    expanded = expanded_support(S, radius, r)
    quotient, remainder = divmod(L, r)
    return quotient * len(expanded) + sum(x < remainder for x in expanded)


def bare_count(S: set[int], r: int, L: int) -> int:
    return sum((l % r) in S for l in range(L))


def branch_summary(r: int, s: int, mixer_count: int, block: np.ndarray,
                   S: set[int]) -> dict:
    L = 1 << s
    radius = mixer_count * (B - 1)
    expanded = expanded_support(S, radius, r)
    total_mass = 0.0
    max_norm_error = 0.0
    max_outside = 0.0
    max_outside_r0 = 0.0
    support_violation = False
    support_violation_r0 = False
    for l in range(L):
        state = branch_before_kick(r, s, l, mixer_count, block)
        norm_error = abs(float(np.vdot(state, state).real) - 1.0)
        max_norm_error = max(max_norm_error, norm_error)
        total_mass += float(np.sum(np.abs(state[list(S)]) ** 2))
        allowed = circular_interval(l % r, radius, r)
        forbidden = [j for j in range(r) if j not in allowed]
        outside = float(np.max(np.abs(state[forbidden]))) if forbidden else 0.0
        max_outside = max(max_outside, outside)
        support_violation |= outside > TOL
        allowed_r0 = {l % r}
        forbidden_r0 = [j for j in range(r) if j not in allowed_r0]
        outside_r0 = (float(np.max(np.abs(state[forbidden_r0])))
                      if forbidden_r0 else 0.0)
        max_outside_r0 = max(max_outside_r0, outside_r0)
        support_violation_r0 |= outside_r0 > TOL
        # This branch-level implication is the pointwise version of the
        # coverage bound: S can receive mass only for l mod r in expanded.
        if np.sum(np.abs(state[list(S)]) ** 2) > TOL and (l % r) not in expanded:
            raise AssertionError("mass reached S outside expanded preimage")
    actual_f = total_mass / L
    cover_count = brute_cover_count(S, radius, r, L)
    quotient_count = quotient_cover_count(S, radius, r, L)
    bare = bare_count(S, r, L)
    scalar_bound = min(1.0, len(S) * (2 * radius + 1) * (1 / r + 1 / L))
    return dict(r=r, s=s, L=L, mixer_count=mixer_count, radius=radius,
                support_size=len(S), support=sorted(S),
                regime="L<r" if L < r else ("L=r" if L == r else "L>r"),
                actual_F=actual_f, cover_count=cover_count,
                quotient_cover_count=quotient_count,
                F_cover=cover_count / L, bare_count=bare,
                bare_F=bare / L, scalar_bound=scalar_bound,
                full_coverage=(radius >= r / 2),
                max_norm_error=max_norm_error,
                max_outside_radius=max_outside,
                max_outside_R0=max_outside_r0,
                R0_support_violation=support_violation_r0,
                support_violation=support_violation,
                zero_mixer=(mixer_count == 0))


def main() -> None:
    exp = Experiment("orbit_lightcone", doc=__doc__)
    exp.predict("P1", "all finite branches stay within radius R=m*(b-1)")
    exp.predict("P2", "brute and quotient coverage counts agree and obey scalar bound")
    exp.predict("P3", "zero-mixer baseline and full-coverage/L<r/L>r regimes behave exactly")
    exp.must_fail("C1", "setting R=0 misses a fixed nonzero-mixer branch support")
    exp.must_fail("C2", "bare unexpanded residue count underbounds actual pre-kick F")

    start = time.perf_counter()
    block = fixed_dense_block()
    rows: list[dict] = []
    max_support_violation = 0.0
    max_bound_excess = 0.0
    max_count_disagreement = 0
    zero_rows = 0
    full_rows = 0
    small_rows = 0
    large_rows = 0

    exp.section("P1/P2/P3 finite light cones and residue coverage")
    radii_cases = 0
    for r in (6, 9, 12, 15, 18):
        for s in (1, 2, 3, 4, 5, 6):
            # Hold one W fixed while changing only the number of prior mixer
            # insertions; include the first full-coverage m when available.
            m_full = math.ceil(r / (2 * (B - 1)))
            counts = sorted({0, min(1, s), min(2, s), min(s, m_full), s})
            supports = (
                {0},
                {r - 1, 0, 1},
                {(r // 2 - 1) % r, r // 2, (r // 2 + 1) % r},
            )
            for mixer_count in counts:
                for S in supports:
                    row = branch_summary(r, s, mixer_count, block, set(S))
                    rows.append(row)
                    radii_cases += 1
                    max_support_violation = max(max_support_violation,
                                                row["max_outside_radius"])
                    max_bound_excess = max(max_bound_excess,
                                           row["actual_F"] - row["F_cover"])
                    max_count_disagreement = max(
                        max_count_disagreement,
                        abs(row["cover_count"] - row["quotient_cover_count"]))
                    if row["regime"] == "L<r":
                        small_rows += 1
                    if row["regime"] == "L>r":
                        large_rows += 1
                    if row["full_coverage"]:
                        full_rows += 1
                    if row["zero_mixer"]:
                        zero_rows += 1
                    exp.check("P1", not row["support_violation"] and row["max_norm_error"] < 3e-12,
                              f"r={r},s={s},m={mixer_count},S={row['support']}: "
                              f"outside-R={row['max_outside_radius']:.2e}")
                    exp.check("P2", row["cover_count"] == row["quotient_cover_count"]
                              and row["actual_F"] <= row["F_cover"] + 3e-12
                              and row["F_cover"] <= row["scalar_bound"] + 3e-15,
                              f"r={r},s={s},m={mixer_count},d={len(S)}: "
                              f"F={row['actual_F']:.6g}, cover={row['F_cover']:.6g}, "
                              f"bound={row['scalar_bound']:.6g}")
                    if mixer_count == 0:
                        exp.check("P3", abs(row["actual_F"] - row["bare_F"]) < 3e-12
                                  and row["cover_count"] == row["bare_count"],
                                  f"zero mixer r={r},s={s},S={row['support']}: "
                                  f"F={row['actual_F']:.6g}, bare={row['bare_F']:.6g}")
                    if row["full_coverage"]:
                        exp.check("P3", row["cover_count"] == row["L"],
                                  f"full coverage r={r},s={s},m={mixer_count}: "
                                  f"{row['cover_count']}/{row['L']}")

    # The two controls use a single fixed row and are deliberately not inferred
    # from the sweep.  A same-row failure of both is required evidence that an
    # R=0/bare count shortcut is not entitled to this nonzero mixer.
    exp.section("must-fail controls")
    # Wrapped support is intentional: this fixed row has both a genuine
    # R=0 support violation and actual F above the bare residue count.
    control = branch_summary(6, 2, 1, block, {4, 5, 0})
    r0_underbound = (control["R0_support_violation"]
                     or control["actual_F"] > control["bare_F"] + 1e-10)
    bare_underbound = control["actual_F"] > control["bare_F"] + 1e-10
    exp.fail_check("C1", r0_underbound,
                   f"fixed r=6,s=2,m=1,S={{4,5,0}}: "
                   f"R0 support violation={control['R0_support_violation']}, "
                   f"F={control['actual_F']:.6g}, bare={control['bare_F']:.6g}")
    exp.fail_check("C2", bare_underbound,
                   f"fixed r=6,s=2,m=1,S={{4,5,0}}: "
                   f"actual F={control['actual_F']:.6g} vs bare={control['bare_F']:.6g}")
    rows.append(dict(series="controls", **control))

    exp.finish(
        report_path=report_path(),
        rows=rows,
        metadata=dict(
            block="fixed 3-point Fourier unitary, reused unchanged across all rows",
            b=B, max_r=MAX_R, max_t=MAX_T, max_dense_bytes=MAX_BYTES,
            tested_case_count=radii_cases, zero_mixer_rows=zero_rows,
            full_coverage_rows=full_rows, L_less_than_r_rows=small_rows,
            L_greater_than_r_rows=large_rows,
            max_support_outside_radius=max_support_violation,
            max_actual_minus_cover=max_bound_excess,
            max_brute_vs_quotient_count_disagreement=max_count_disagreement,
            fixed_control=control,
            references="bounded finite matrix products only; no generic post-kick propagator",
            interpretation="pre-kick work mass on supplied label support, conditioned on uniform early l",
            elapsed_seconds=time.perf_counter() - start,
        ),
    )


if __name__ == "__main__":
    main()
