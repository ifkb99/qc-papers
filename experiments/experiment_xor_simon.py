"""Does Contract-B Grover output on the factor predicate hand over p xor q, and does that factor N?

Context: draft C122 (submission Sc377b4ddf8aa4d5c, section 6) left open whether
exact Contract-B sampling (Grover from uniform, then H on every bit) for C119's
two-orientation semiprime predicate f_N(p,q)=[pq=N, 2<=p,q<2^m] is hard. Its
derived law: P(0)=cos^2(2t theta), and conditioned on y!=0, y=(y1,y2) is
uniform on {y!=0 : (y1 xor y2).s = 0}, s = p xor q (Simon-like). With the
canonical ordering p<=q (M=1) the conditional law is uniform on y!=0.
About m samples give s by GF(2) linear algebra. LSB branch-and-prune from
(N,s) is a known method (github.com/sliedes/xor_factor, after a Math.SE answer;
no complexity claim there). Heuristic derivation: at bit k the product bit is
p_k+q_k+c(lower bits) = s_k+c, so each node keeps two xor-consistent children
and the next check prunes about half: a critical branching process, whose
tracked set is expected to grow linearly in m. Heuristic, not a theorem.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  Dense reference (vector updates + walsh.wht; the derived law is not used
      to compute it), m=4..7, t=round(pi/(4 theta)): |P(0)-cos^2(2t theta)|<1e-9
      and the y!=0 conditional equals the uniform-complement law to 1e-9
      relative; canonical M=1 conditional is uniform on y!=0 to the same bound.
  P2  For m in {32,64,128,256}, 20 semiprimes each, m+8 samples of the derived
      law give a GF(2) kernel equal to exactly {0,s} in every instance.
  P3  Branch-and-prune with the true s returns (p,q) in every instance.
  P4  (heuristic) max tracked-set size <= 4m in every instance; the ratio
      mean(max tracked)/m at m=256 is at most twice that at m=32.
  C0  must fail: omitting the oracle leaves no y!=0 mass (P1 is not vacuous).
  C1  must fail: branch-and-prune WITHOUT the hint keeps max tracked <= 4m (m=16).
  C2  must fail: a random wrong hint s' (bit 0 = 0) still factors N (m=64).

Run:  uv run python -m experiments.experiment_xor_simon     (from research/)

Reviewed as S3ad6a279f9f2462c (script sha256 82ac0380...12e5, run Rf0b35178211640b5,
review Va27665ae402a41b6). Path-only edits for experiments/: ROOT, OUT and this Run line.
C0 shows only that y!=0 mass needs the oracle; C2 is guaranteed by the xor invariant and
does not discriminate (P3 is the check that catches a disabled filter). See C123, which
supersedes this header's heuristic context paragraph (Heninger-Shacham attribution, scope).
"""
from __future__ import annotations
import json
import math
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402
from lab.harness import Experiment  # noqa: E402
from lab.gf2 import rank_kernel  # noqa: E402
from walsh import wht  # noqa: E402

OUT = ROOT / "out" / "experiment_xor_simon"
rng = random.Random(20260921)

exp = Experiment("experiment_xor_simon", doc=__doc__)
exp.predict("P1", "dense Contract-B law = cos^2 coin + uniform complement of (s,s); M=1 uniform")
exp.predict("P2", "m+8 derived-law samples give kernel exactly {0,s}")
exp.predict("P3", "branch-and-prune with true s factors N every time")
exp.predict("P4", "heuristic: max tracked <= 4m; mean(max)/m ratio growth <= 2x from m=32 to 256")
exp.must_fail("C0", "omitted oracle leaves nonzero y!=0 mass")
exp.must_fail("C1", "no-hint branch-and-prune stays <= 4m at m=16")
exp.must_fail("C2", "random wrong hint still factors N at m=64")


def is_prime(n):
    if n < 2:
        return False
    for sp in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % sp == 0:
            return n == sp
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for _ in range(32):
        a = rng.randrange(2, n - 1)
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def random_prime(bits):
    while True:
        c = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if is_prime(c):
            return c


def semiprime(bits):
    while True:
        p, q = random_prime(bits), random_prime(bits)
        if p != q:
            return p, q


def dense_contract_b(m, N, canonical, oracle=True):
    """Dense route: returns (theta, t, post-H probabilities). Registers: x = p + (q << m)."""
    n = 2 * m
    x = np.arange(1 << n, dtype=np.uint64)
    P, Q = x & np.uint64((1 << m) - 1), x >> np.uint64(m)
    sel = (P * Q == np.uint64(N)) & (P >= 2) & (Q >= 2)
    if canonical:
        sel &= P <= Q
    M = int(sel.sum())
    theta = math.asin(math.sqrt(M / (1 << n)))
    t = round(math.pi / (4 * theta))
    state = np.ones(1 << n)
    for _ in range(t):
        if oracle:
            state[sel] *= -1
        state = 2 * state.mean() - state
    return theta, t, (wht(state) / (1 << n)) ** 2


def parity(v):
    return v.bit_count() & 1


def sample_derived(m, s):
    while True:
        y1, y2 = rng.getrandbits(m), rng.getrandbits(m)
        if (y1 or y2) and parity((y1 ^ y2) & s) == 0:
            return y1 ^ y2


def branch_and_prune(N, m, s=None, cap=200_000):
    """LSB branch-and-prune on sorted (p mod 2^k, q mod 2^k); s=None means no hint."""
    tracked = {(1, 1)}
    sizes = [1]
    for k in range(2, m + 1):
        bit, mask = 1 << (k - 1), (1 << k) - 1
        nk = N & mask
        new = set()
        for a, b in tracked:
            for pa in (a, a | bit):
                for qb in (b, b | bit):
                    if s is not None and ((pa ^ qb) & bit) != (s & bit):
                        continue
                    if (pa * qb) & mask == nk:
                        new.add((min(pa, qb), max(pa, qb)))
        tracked = new
        sizes.append(len(tracked))
        if len(tracked) > cap:
            return None, sizes
    for a, b in tracked:
        if a * b == N and a > 1:
            return (a, b), sizes
    return None, sizes


rows = {"P1": [], "P2P4": [], "C1": None, "C2": []}

exp.section("P1  dense reference vs derived Contract-B law, m=4..7")
p1_ok, worst = True, 0.0
for m in range(4, 8):
    primes = [v for v in range(3, 1 << m) if is_prime(v)]
    for _ in range(3):
        p, q = rng.sample(primes, 2)
        N, s = p * q, p ^ q
        n = 2 * m
        y = np.arange(1 << n, dtype=np.uint64)
        z = (y & np.uint64((1 << m) - 1)) ^ (y >> np.uint64(m))
        zpar = np.array([parity(int(v) & s) for v in z])
        for canonical in (False, True):
            theta, t, prob = dense_contract_b(m, N, canonical)
            coin_err = abs(prob[0] - math.cos(2 * t * theta) ** 2)
            cond = prob[1:] / prob[1:].sum()
            law = np.ones((1 << n) - 1) if canonical else (zpar[1:] == 0).astype(float)
            law /= law.sum()
            rel = float(np.abs(cond - law).max() / law.max())
            worst = max(worst, coin_err, rel)
            p1_ok &= coin_err < 1e-9 and rel < 1e-9
            rows["P1"].append(dict(m=m, p=p, q=q, canonical=canonical, t=t,
                                   coin_err=coin_err, cond_rel_err=rel, p_nonzero=float(prob[1:].sum())))
exp.check("P1", p1_ok, f"worst discrepancy {worst:.3e} over {len(rows['P1'])} cases")

exp.section("C0  omitted oracle")
_, _, prob0 = dense_contract_b(6, 7 * 11, False, oracle=False)
exp.fail_check("C0", float(prob0[1:].sum()) < 1e-12, f"y!=0 mass {float(prob0[1:].sum()):.3e}")

exp.section("P2-P4  s recovery and branch-and-prune, m in {32,64,128,256}")
p2_ok = p3_ok = p4_bound = True
ratio = {}
for m in (32, 64, 128, 256):
    maxes = []
    for _ in range(20):
        p, q = semiprime(m)
        N, s = p * q, p ^ q
        zs = [sample_derived(m, s) for _ in range(m + 8)]
        rank, kernel = rank_kernel(zs, m)
        rec = kernel[0] if len(kernel) == 1 else None
        p2_ok &= rec == s
        found, sizes = branch_and_prune(N, m, rec if rec is not None else s)
        p3_ok &= found == (min(p, q), max(p, q))
        mx = max(sizes)
        maxes.append(mx)
        p4_bound &= mx <= 4 * m
        rows["P2P4"].append(dict(m=m, rank=rank, recovered=rec == s, factored=found is not None,
                                 max_tracked=mx, final_tracked=sizes[-1]))
    ratio[m] = sum(maxes) / len(maxes) / m
    exp.log(f"m={m}: mean max tracked {sum(maxes)/len(maxes):.1f}, max {max(maxes)}, mean/m {ratio[m]:.3f}")
exp.check("P2", p2_ok, "kernel == {0,s} in all 80 instances" if p2_ok else "recovery failed")
exp.check("P3", p3_ok, "factored in all 80 instances" if p3_ok else "factor missing")
exp.check("P4", p4_bound and ratio[256] <= 2 * ratio[32],
          f"bound ok={p4_bound}; mean(max)/m: " + ", ".join(f"{k}:{v:.3f}" for k, v in ratio.items()))

exp.section("C1  no hint, m=16")
p, q = semiprime(16)
_, sizes = branch_and_prune(p * q, 16, None)
rows["C1"] = dict(m=16, sizes=sizes)
exp.fail_check("C1", max(sizes) > 4 * 16, f"max tracked {max(sizes)} (sizes {sizes})")

exp.section("C2  random wrong hint, m=64")
succ = 0
for _ in range(10):
    p, q = semiprime(64)
    s = p ^ q
    while True:
        wrong = rng.getrandbits(64) & ~1
        if wrong != s:
            break
    found, sizes = branch_and_prune(p * q, 64, wrong)
    succ += found is not None
    rows["C2"].append(dict(found=found is not None, max_tracked=max(sizes)))
exp.fail_check("C2", succ == 0, f"{succ}/10 wrong-hint runs factored")

exp.finish(report_path=OUT / "report_v1.json", rows=rows,
           metadata={"seed": 20260921, "script": "experiment_xor_simon.py"})
