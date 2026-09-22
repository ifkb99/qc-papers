"""Does a Pauli-LIMDD of the post-modexp state escape the exponential in exponent-first orders?

Context: register row TX29 lists Pauli-LIMDDs as an open escape with no lower bound
known; C107/C108 cover ordered *linear* representations only.  The object here is
the state a decision-diagram simulator holds after modular exponentiation,
|Psi> = sum_{e < 2^t} |e>|a^e mod N>|0...0>, produced by ToffoliModExp.build() on
the clean code.  Derivation (note LM): in any order that reads all exponent qubits
before the x qubits (scratch anywhere), let j be the last exponent qubit read and
A = a^(2^j) mod N.  The sub-states at the node level of e_j are the pair states
|0>|c> + |1>|Ac mod N>, c = a^e' (e' with bit j clear), and

  L1 (exact)  two of them are Pauli-LIM equivalent iff phi(c) = phi(c'),
              phi(c) = c XOR (Ac mod N);
  L2 (fiber)  if A - 1 and A + 1 are units mod N, every fiber of phi on [0, N) has
              at most 2^min(wt d, n - wt d) <= 2^floor(n/2) points;
  L3          D_j = #{a^e' : e' < 2^t, bit_j(e') = 0} >= ceil(r/2) when 2^t >= r.

So the LIMDD has >= ceil(r/2) / 2^floor(n/2) nodes at that level (> r/(2 sqrt(2N))).

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P0  build() maps |e>|0> to |e>|a^e mod N> with all scratch returned to 0, for every
      e < 2^t (ties the object to the construction), n <= 5.
  P1  L1: brute-force LIM class counts (orbit canonicalisation over all XOR shifts,
      no use of phi) equal #distinct phi, at every (N, a, j) tested, n <= 7.
      (First written as n <= 9; the code brute-forces n <= 7 and uses phi alone at
      n = 8, 9. Corrected after run 1, which checked n <= 7 only.)
  P2  L2: max fiber of phi over [0, N) is <= 2^min(wt d, n - wt d) whenever A -+ 1 are
      units, every odd N with 5 <= N < 2^9 and every A in [2, N) with A - 1 and A + 1
      units mod N (A itself need not be a unit; "every unit A" until v2).
  P3  L3 and the bound: D_j >= ceil(r/2) and classes >= ceil(D_j / 2^floor(n/2)) at
      every tested case satisfying the hypothesis.
  P4  L1 including phases and scalars: at n <= 3 an explicit search over all
      lambda * P (P in {I,X,Y,Z}^(n+1)) on the dense vectors gives the same classes.
  P5  OPEN, no commitment: the ratio classes / D_j.  H-random: close to 1 (phi behaves
      like a random map into ~N values); H-collapse: small.  Reported, not graded.
  P6  POST-HOC (added after run 1 and two scratch scans, so NOT predictions; logged
      only): (a) classes / D_j for 40 random semiprimes per n = 10, 12, 14, 16 (seed 1);
      (b) x-first order, level after every x qubit: XOR-shift classes of the
      sub-states 1[e = l mod r] on [0, 2^t), against ceil(beta/2), r = beta * 2^alpha.
      v2 (after review Va10c9887fa3f4116, B1): run 2's P6a took the prefix set as
      {a^e : e < r} and divided by r.  That is |phi(<a>)| / r, which differs from
      classes / D_j at j = 0 for even r (C_0 = <a^2>).  P6a now uses prefix_values and
      logs min / median / max of classes / D_j, and the set of D_j / r.

Scope of P0-P5 (v2 note, review B2): cases() takes the smallest unit of order > 2,
which is 2 for every odd N >= 5 (ord_N(2) <= 2 only when N | 3).  Every P0-P5 case
therefore uses base a = 2, and the sweep runs over odd N, not primes.  The bases used are
logged.
  C1  must fail: N = 17, a = 3 (r = 16, a power of 2), t = 6, j = 5 gives A = 1,
      outside the hypothesis; the bound's conclusion (>= 16/4 = 4 classes) must fail.
  C2  must fail: the mutant invariant c + (Ac mod N) (integer sum, not XOR) must
      disagree with the brute-force class count somewhere, else P1 is not
      discriminating.

Run:  uv run python -m experiments.experiment_limdd_modexp     (from research/)
"""
from __future__ import annotations

import itertools
import math

import numpy as np

from lab import Experiment
from lab.modarith import order
from toffoli_arith import ToffoliModExp
from walsh import classical_images

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__)

exp.predict("P0", "build() on |e>|0> gives |e>|a^e mod N>, scratch clean, n <= 5")
exp.predict("P1", "brute-force LIM classes == #distinct phi, n <= 7")
exp.predict("P2", "fiber of phi <= 2^min(wt d, n - wt d) when A -+ 1 are units")
exp.predict("P3", "D_j >= ceil(r/2) and classes >= ceil(D_j / 2^floor(n/2)) under the hypothesis")
exp.predict("P4", "dense Pauli-LIM search (with lambda, Z, Y) agrees with support-shift classes, n <= 3")
exp.must_fail("C1", "N=17, a=3, t=6, j=5 (A = 1) violates the bound's conclusion")
exp.must_fail("C2", "integer-sum mutant invariant disagrees with brute force somewhere")


def primes_of(N):
    out, m, p = [], N, 2
    while p * p <= m:
        if m % p == 0:
            out.append(p)
            while m % p == 0:
                m //= p
        p += 1
    if m > 1:
        out.append(m)
    return out


def hypothesis(N, A):
    return all(A % p not in (1, p - 1) for p in primes_of(N))


def prefix_values(N, a, t, j):
    """c = a^e' mod N for e' < 2^t with bit j clear (distinct values, sorted)."""
    es = [e for e in range(1 << t) if not (e >> j) & 1]
    return sorted({pow(a, e, N) for e in es})


def phi(c, A, N):
    return c ^ (A * c % N)


def brute_classes(cs, A, N, n):
    """LIM classes of pair supports {(0,c),(1,Ac)} by orbit canonicalisation."""
    canon = set()
    for c in cs:
        pts = ((0, c), (1, A * c % N))
        best = None
        for se in (0, 1):
            for sx in range(1 << n):
                key = tuple(sorted((b ^ se, y ^ sx) for b, y in pts))
                if best is None or key < best:
                    best = key
        canon.add(best)
    return len(canon)


def cases(nmax, t_of):
    for N in range(5, 1 << nmax, 2):
        n = N.bit_length()
        units = [a for a in range(2, N) if math.gcd(a, N) == 1]
        # one parameter varies (N); the base is the smallest unit with r > 2
        for a in units:
            if order(a, N) > 2:
                yield N, n, a, t_of(n)
                break


# ---------------------------------------------------------------- P0
exp.section("P0  the object is the circuit's clean-code output")
bad0 = 0
checked0 = 0
for N, n, a, t in cases(5, lambda n: n + 1):
    me = ToffoliModExp(N, a, n_exp=t)
    qc = me.build()
    ins = np.array([sum(((e >> i) & 1) << q for i, q in enumerate(me.exp))
                    for e in range(1 << t)], dtype=np.int64)
    outs = classical_images(qc, ins)
    for e, y in zip(range(1 << t), outs):
        xval = sum(((int(y) >> q) & 1) << i for i, q in enumerate(me.x))
        rest = int(y) & ~sum(1 << q for q in me.x) & ~sum(1 << q for q in me.exp)
        eval_ = sum(((int(y) >> q) & 1) << i for i, q in enumerate(me.exp))
        checked0 += 1
        if xval != pow(a, e, N) or rest != 0 or eval_ != e:
            bad0 += 1
exp.check("P0", checked0 > 0 and bad0 == 0, f"{checked0} basis inputs, {bad0} mismatches")

# ---------------------------------------------------------------- P1, P3, P5
exp.section("P1/P3/P5  classes at the last exponent level, every j")
rows = []
bad1 = bad3 = 0
n_hyp = 0
min_ratio = (2.0, None)
for N, n, a, t in cases(9, lambda n: n + 2):
    r = order(a, N)
    if (1 << t) < r:
        continue
    for j in range(t):
        A = pow(a, 1 << j, N)
        cs = prefix_values(N, a, t, j)
        k_phi = len({phi(c, A, N) for c in cs})
        k_bf = brute_classes(cs, A, N, n) if n <= 7 else k_phi
        if n <= 7 and k_bf != k_phi:
            bad1 += 1
        if hypothesis(N, A):
            n_hyp += 1
            ok = (len(cs) >= math.ceil(r / 2)
                  and k_phi >= math.ceil(len(cs) / 2 ** (n // 2)))
            bad3 += not ok
            ratio = k_phi / len(cs)
            if ratio < min_ratio[0]:
                min_ratio = (ratio, (N, a, j, k_phi, len(cs)))
        rows.append(dict(N=N, n=n, a=a, r=r, t=t, j=j, A=A, D=len(cs),
                         classes=k_phi, hyp=hypothesis(N, A)))
exp.check("P1", bad1 == 0, f"{sum(1 for x in rows if x['n'] <= 7)} brute-force cases, {bad1} disagreements")
exp.check("P3", n_hyp > 0 and bad3 == 0, f"{n_hyp} cases under the hypothesis, {bad3} violations")
exp.log(f"bases used in P1/P3/P5: {sorted({x['a'] for x in rows})}")
exp.log(f"P5 min classes/D_j under the hypothesis: {min_ratio[0]:.3f} at (N, a, j, classes, D) = {min_ratio[1]}")
big = [x for x in rows if x['hyp'] and x['n'] == 9]
if big:
    worst = min(big, key=lambda x: x['classes'] / x['D'])
    exp.log(f"P5 at n = 9: {len(big)} cases, worst {worst}")

# ---------------------------------------------------------------- P2
exp.section("P2  fiber bound for phi")
bad2 = 0
nfib = 0
for N in range(5, 1 << 9, 2):
    n = N.bit_length()
    for A in range(2, N):
        if math.gcd(A - 1, N) != 1 or math.gcd(A + 1, N) != 1:
            continue
        vals = np.arange(N, dtype=np.int64)
        ph = vals ^ (A * vals % N)
        cnt = np.bincount(ph, minlength=1 << n)
        nfib += 1
        for d in np.nonzero(cnt)[0]:
            w = bin(int(d)).count("1")
            if cnt[d] > 2 ** min(w, n - w):
                bad2 += 1
                break
exp.check("P2", nfib > 0 and bad2 == 0, f"{nfib} (N, A) pairs, {bad2} violating")

# ---------------------------------------------------------------- P4
exp.section("P4  dense Pauli-LIM search, n <= 3")
I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]])
Z = np.diag([1, -1])
Y = 1j * X @ Z
PAULIS = [I2, X, Y, Z]


def lim_equiv(u, v, q):
    for ps in itertools.product(PAULIS, repeat=q):
        M = ps[0]
        for p in ps[1:]:
            M = np.kron(M, p)   # qubit order irrelevant: all tensor products tried
        w = M @ u
        k = np.flatnonzero(np.abs(w) > 1e-12)
        if len(k) and np.array_equal(k, np.flatnonzero(np.abs(v) > 1e-12)):
            lam = v[k[0]] / w[k[0]]
            if np.allclose(lam * w, v):
                return True
    return False


bad4 = 0
n4 = 0
for N, n, a, t in cases(3, lambda n: n + 2):
    for j in range(t):
        A = pow(a, 1 << j, N)
        cs = prefix_values(N, a, t, j)
        q = n + 1
        vecs = []
        for c in cs:
            v = np.zeros(1 << q, dtype=complex)
            v[(0 << n) | c] = 1          # e_j is the top qubit
            v[(1 << n) | (A * c % N)] = 1
            vecs.append(v)
        reps = []
        for v in vecs:
            if not any(lim_equiv(v, u, q) for u in reps):
                reps.append(v)
        n4 += 1
        bad4 += len(reps) != len({phi(c, A, N) for c in cs})
exp.check("P4", n4 > 0 and bad4 == 0, f"{n4} cases, {bad4} disagreements")

# ---------------------------------------------------------------- P6 (post-hoc, logged)
exp.section("P6a  POST-HOC: random semiprimes, classes / D_j (logged, not graded)")
import random
random.seed(1)
for n in (10, 12, 14, 16):
    rat, worst, per_r, dfrac = [], None, [], set()
    for _ in range(40):
        while True:
            N = random.randrange(2 ** (n - 1) + 1, 2 ** n, 2)
            ps = primes_of(N)
            if len(ps) == 2 and ps[0] * ps[1] == N and min(ps) > 7:
                break
        while True:
            a = random.randrange(2, N)
            if math.gcd(a, N) == 1:
                break
        r = order(a, N)
        t = max(1, math.ceil(math.log2(r))) + 1
        for j in (0, t - 1):
            A = pow(a, 1 << j, N)
            if not hypothesis(N, A):
                continue
            Cj = np.array(prefix_values(N, a, t, j), dtype=np.int64)
            # detector for review B1: C_j must be the bit-j-clear prefix set, not <a>
            es = np.array([e for e in range(1 << t) if not (e >> j) & 1])
            assert len(Cj) == len(np.unique(es % r)), "prefix set is not C_j"
            k = len(np.unique(Cj ^ (A * Cj % N)))
            q = k / len(Cj)
            rat.append(q)
            dfrac.add(round(len(Cj) / r, 3))
            per_r.append(k / r)
            if worst is None or q < worst[0]:
                worst = (round(q, 3), N, a, r, j, len(Cj))
    exp.log(f"n={n}: {len(rat)} cases, classes/D_j min {min(rat):.3f}, median {np.median(rat):.3f}, "
            f"max {max(rat):.3f}; classes/r min {min(per_r):.3f} max {max(per_r):.3f}; "
            f"D_j/r values {sorted(dfrac)}; worst (ratio, N, a, r, j, D_j) = {worst}; "
            f"proved floor 2^-floor(n/2) = {2 ** -(n // 2):.4f}")

exp.section("P6b  POST-HOC: x-first order, XOR-shift classes of 1[e = l mod r] (logged)")


def ap_classes(r, t):
    M = 1 << t
    sh = np.arange(M)[:, None]
    parent = list(range(r))

    def root(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for l in range(r):
        ap = np.arange(l, M, r)[None, :]
        res = (ap ^ sh) % r
        same = (res == res[:, :1]).all(axis=1)
        for l2 in np.unique(res[same, 0]):
            if len(range(int(l2), M, r)) == ap.shape[1]:
                parent[root(l)] = root(int(l2))
    return len({root(x) for x in range(r)})


for r in (12, 15, 23, 31, 45, 63, 100, 127):
    beta = r
    while beta % 2 == 0:
        beta //= 2
    tb = math.ceil(math.log2(r))
    exp.log(f"r={r} beta={beta}: classes at t = {tb + 1}..{2 * tb + 2}: "
            f"{[ap_classes(r, t) for t in range(tb + 1, 2 * tb + 3)]}; ceil(beta/2) = {(beta + 1) // 2}")

# ---------------------------------------------------------------- C1, C2
exp.section("C1  outside the hypothesis: A = 1")
N, a, t, j = 17, 3, 6, 5
A = pow(a, 1 << j, N)
cs = prefix_values(N, a, t, j)
k = brute_classes(cs, A, N, N.bit_length())
exp.log(f"A = {A}, D = {len(cs)}, classes = {k}, bound would claim >= {math.ceil(len(cs) / 2 ** (N.bit_length() // 2))}")
exp.fail_check("C1", k < math.ceil(len(cs) / 2 ** (N.bit_length() // 2)), "bound conclusion fails at A = 1")

exp.section("C2  mutant invariant")
dis = 0
for x in rows:
    if x['n'] > 7:
        continue
    cs = prefix_values(x['N'], x['a'], x['t'], x['j'])
    mut = len({c + x['A'] * c % x['N'] for c in cs})
    dis += mut != brute_classes(cs, x['A'], x['N'], x['n'])
exp.fail_check("C2", dis > 0, f"integer-sum invariant disagrees in {dis} cases")

import os
exp.finish(report_path=os.environ.get("LIMDD_REPORT", "out/limdd_modexp/report.json"), rows=rows)
