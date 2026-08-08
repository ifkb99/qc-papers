"""TODO 12d -- why is the peak ratio exactly `rot = 2*perm - 2` for modexp?

Paper A SS5's factor-of-two peak reduction is the only quantitative claim in
that paper that is measured rather than derived. It is 2.000000 for ripple-carry
adders (128/64, 512/256, 2048/1024) and 1.9997 for modular exponentiation, where
every instance instead satisfies

    rot = 2*perm - 2        13666 = 2*6834 - 2      (N=5,  a=2)
                            16386 = 2*8194 - 2      (N=7,  a=3)
                           128138 = 2*64070 - 2     (N=15, a=7)

exactly, 3 of 3, across wildly different magnitudes. A constant additive deficit
of 2 that does not scale with the instance is not a coincidence.

READING THE CONSTRUCTION FIRST (this is where the predictions come from).

Under Heisenberg propagation a Z-type string commutes with every Z-rotation, so
it cannot branch until an H turns a Z into an X. In the Toffoli gadget the only
H acts on the target c, so ONLY strings carrying Z_c ever leave the diagonal.
Once such a string has become X_c-type, the four T gates on c rotate it inside
the two-dimensional space spanned by {X_c, Y_c}:

    X_c -> cos X_c - sin Y_c        Y_c -> cos Y_c + sin X_c

which is CLOSED. So each Z_c-carrying string contributes exactly 2 Paulis, not
2^4, and each string missing Z_c contributes exactly 1, untouched. Writing S for
the live Z-set and B = {z in S : z_c = 0},

    rot_peak = 2|S| - |B|

**so the deficit IS |B|, the number of live strings that miss the Toffoli's own
target qubit.** "rot = 2*perm - 2" says exactly two live strings miss it;
"rot = 2*perm" for adders says none do. That turns an unexplained constant into
a set we can print.

FIRST PASS (P1, P2) WAS REFUTED, AND THE REFUTATION IS THE INTERESTING PART.
The identity above is right, but the first attempt paired the peak against the
Z-set at the *gadget boundary*, reconstructed by re-emitting the circuit and
mapping gate indices to logical ops. That gave 2056 where the peak needs 6834,
and P1/P2 failed 0/9 while P4/P5 -- which only counted B -- passed 6/6. The
diagnosis, from dumping the peak set: the peak sits deep inside the gadget,
*after* the Toffoli's own 4-way expansion has run at rotation level, so no
gadget boundary carries the right set. The correct partner is the perm-level
peak itself, and the gadget target c does not need reconstructing at all --
it can be read straight off the peak set, which is the second finding:

**at the peak every Pauli has x-mask either 0 or exactly (1 << c).** Two Paulis
are pure Z-type; all 13,664 others share the single x-mask {t3}. The
gate-boundary machinery was solving a problem that did not exist.

PREDICTIONS FOR THIS PASS, WRITTEN BEFORE MEASURING.
(P1-P5 are kept verbatim below in `REFUTED` so the record shows the order.)

  P6  At the rot-level peak, every Pauli has x-mask in {0, 1<<c} for a single
      qubit c. Nothing branches into a second X coordinate, because only one
      qubit ever meets an H while carrying a Z.

  P7  Fold each peak Pauli (x, z) to the Z-string  z | (1<<c) if x != 0 else z.
      The result is exactly the perm-level PEAK set S -- as a set -- with
      multiplicity 2 on {z in S : z_c = 1} and 1 on {z in S : z_c = 0}. Hence
      rot_max = 2*perm_max - |B| exactly, with B = {z in S : z_c = 0}.

  P8  |B| = 2 for every modexp instance and 0 for every adder, and in every
      modexp instance the two members of B are the SAME shape: Z on the
      observed x-register bit alone, and that bit together with one exponent
      qubit.

  P9  Those two are not propagation trivia. SS W3/W4 found the modexp spectrum
      dominated by exactly two coefficients of magnitude 1/2, supported on the
      measured x bit plus one exponent qubit. **Predict: B is exactly the set
      of |c| = 1/2 coefficients of the final spectrum.** If it holds, the
      deficit is the dominant Fourier modes, and it is 2 because there are two
      of them.

  C1  MUST FAIL -- ripple adders, where the ratio is exactly 2.000000. There
      |B| must be 0: every live string at the peak carries Z_c.

  C2  MUST FAIL, and it is the sharper control -- **the deficit must not be a
      property of "modexp".** If the mechanism is "strings missing the peak
      gadget's target do not double", then keeping the circuit fixed and moving
      the OBSERVABLE to a qubit that sits inside the scratch register must
      change the deficit away from 2. If it stays 2 for every observable, the
      explanation above is decoration and the constant is something else.

Method note: no new propagator is written. Heisenberg propagation of the last m
gates IS the state after m steps, so both peak states come from running the
existing verified propagators on a gate/logical suffix.

OUTCOME (2026-08-08): 6/6, and TODO 12d is closed. rot_max = 2*perm_max - |B|
holds exactly on 9/9 instances, with the peak confined to {I, X_c} for a single
c and the fold recovering the perm peak set for set. B = {Z_x0, Z_x0 Z_e0} in
6/6 modexp instances, and those are exactly the |coefficient| = 1/2 Walsh terms
-- the same pair SS W3/W4 found from the final spectrum by an unrelated route.
C2 fires hard: moving the observable gives deficits 4004 (b0), 4014 (anc) and 0
(t0, ratio exactly 2), so the constant belongs to the OBSERVABLE and not to
modular exponentiation. Written up in NOTES.md SS PK; claim C44; Paper A SS5
rewritten from "an observation" to Proposition 2.

Still NOT proved: |B| = 2 itself. The identity is derived, the membership of B
is measured. Do not upgrade it without a characterisation of which peak-time
strings avoid the scratch register.

Run:  uv run python -m experiments.experiment_c17_deficit     (~4 min)

REFUTED (first pass, kept for the record):
  P1  rot_max == 2|A| + |B| at the enclosing GADGET BOUNDARY.   FAILED 0/9 --
      right identity, wrong boundary; the peak is mid-gadget.
  P2  the peak Pauli set is that BOUNDARY set, A doubled.       FAILED 0/9 --
      same cause. The fold itself was correct and is reused as P7.
  P3  the trajectory plateaus rather than doubling per T gate.  PASSED.
  P4  deficit == 2 beyond the three logged instances.           PASSED 6/6.
  P5  is B the same recognisable pair every time?               YES, see P8.
"""
from __future__ import annotations

import time

import numpy as np

from lab import Experiment
from circuits import ripple_adder
from toffoli_arith import ToffoliModExp
import pps
import perm_pps
import walsh

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__)

exp.predict("P6", "at the peak every Pauli has x-mask in {0, 1<<c}, one c")
exp.predict("P7", "folding the peak set gives the perm PEAK set, "
                  "multiplicity 2 on z_c=1 and 1 on z_c=0")
exp.predict("P8", "|B| = 2 for modexp, 0 for adders, same shape every time")
exp.predict("P9", "B is exactly the set of |coefficient| = 1/2 Walsh terms")
exp.must_fail("C1", "adders must have |B| = 0")
exp.must_fail("C2", "moving the observable into the scratch register must "
                    "change the deficit away from 2")


class _GateSuffix:
    """The last `m` rotations of a circuit -- i.e. the state after m steps."""

    def __init__(self, circuit, m):
        self.n = circuit.n
        self.gates = circuit.gates[len(circuit.gates) - m:]


class _LogicalSuffix:
    """The last `j` logical ops, for perm_pps."""

    def __init__(self, circuit, j):
        self.n = circuit.n
        self.logical = circuit.logical[len(circuit.logical) - j:]
        self.gates = []

    def is_classical(self):
        return True


def analyse(label, circuit, target_qubit):
    t0 = time.time()
    rot = pps.propagate(circuit, {(0, 1 << target_qubit): 1.0}, delta=0.0)
    perm = perm_pps.propagate_perm(circuit, 1 << target_qubit)

    m_star = int(np.argmax(rot.n_terms)) + 1
    peak = pps.propagate(_GateSuffix(circuit, m_star),
                         {(0, 1 << target_qubit): 1.0}, delta=0.0).final_terms

    # P6: read the gadget target straight off the peak set instead of
    # reconstructing gate boundaries (which is what broke the first pass)
    # P6 says the x-masks lie IN {0, 1<<c}; it does not say both occur. When
    # |B| = 0 nothing stays pure Z, so the identity sector is simply absent --
    # which is the same fact seen from the other side, not a counterexample.
    xmasks = sorted({x for (x, _) in peak})
    nonzero = [x for x in xmasks if x]
    single_c = (len(nonzero) == 1 and nonzero[0].bit_count() == 1
                and set(xmasks) <= {0, nonzero[0]})
    c = nonzero[0].bit_length() - 1 if single_c else None

    # P7: fold X_c/Y_c partners back onto their parent Z-string
    fold: dict[int, int] = {}
    if c is not None:
        for (x, z) in peak:
            key = (z | (1 << c)) if x else z
            fold[key] = fold.get(key, 0) + 1

    j_star = int(np.argmax(perm.n_terms)) + 1
    S = set(perm_pps.propagate_perm(_LogicalSuffix(circuit, j_star),
                                    1 << target_qubit).final_terms)
    B = sorted(z for z in S if c is not None and not ((z >> c) & 1))

    return {"label": label, "q": circuit.n, "rot_max": rot.n_max,
            "perm_max": perm.n_max, "deficit": 2 * perm.n_max - rot.n_max,
            "xmasks": xmasks, "single_c": single_c, "c": c,
            "fold": fold, "S": S, "B": B,
            "traj": rot.n_terms, "m_star": m_star, "secs": time.time() - t0}


def header():
    print(f"    {'instance':>18} {'q':>3} {'obs':>5} {'perm':>8} {'rot':>9} "
          f"{'deficit':>8} {'peak x-masks':>16} {'c':>4} {'|B|':>5} {'secs':>7}")


def show(r, obs=""):
    xm = "{" + ",".join(str(x) for x in r["xmasks"]) + "}"
    print(f"    {r['label']:>18} {r['q']:3d} {obs:>5} {r['perm_max']:8,d} "
          f"{r['rot_max']:9,d} {r['deficit']:8d} {xm:>16} "
          f"{str(r['c']):>4} {len(r['B']):5d} {r['secs']:7.1f}", flush=True)


def qubit_names(me):
    names = {}
    for i, q in enumerate(me.b):
        names[q] = f"b{i}"
    for i, q in enumerate(me.t):
        names[q] = f"t{i}"
    for i, q in enumerate(me.x):
        names[q] = f"x{i}"
    names[me.c0] = "c0"
    names[me.anc] = "anc"
    for i, q in enumerate(me.exp):
        names[q] = f"e{i}"
    return names


def render(z, names, n):
    if z == 0:
        return "IDENTITY"
    return "{" + ",".join(names.get(k, f"q{k}")
                          for k in range(n) if (z >> k) & 1) + "}"


# ---------------------------------------------------------------------------
exp.section("modexp -- three logged instances plus three never used for this")
header()
MODEXP = [("N=5 a=2", 5, 2), ("N=7 a=3", 7, 3), ("N=15 a=7", 15, 7),
          ("N=5 a=4", 5, 4), ("N=7 a=6", 7, 6), ("N=15 a=2", 15, 2)]
mres = []
for label, N, a in MODEXP:
    me = ToffoliModExp(N=N, a=a, n_exp=1)
    r = analyse(label, me.build(), me.x[0])
    r["me"] = me
    mres.append(r)
    show(r, "x0")

exp.section("C1 control -- ripple-carry adders, ratio exactly 2")
header()
ares = []
for nb in (3, 4, 5):
    add, lay = ripple_adder(nb)
    ares.append(analyse(f"{nb}-bit adder", add, lay["b"][0]))
    show(ares[-1], "b0")

# ---------------------------------------------------------------------------
exp.section("P6  is the peak confined to a single X coordinate?")
for r in mres + ares:
    exp.log(f"{r['label']:>18}: x-masks {r['xmasks']} -> "
            + (f"single c = {r['c']}"
               + ("" if 0 in r["xmasks"] else "  (no identity sector: |B|=0)")
               if r["single_c"] else "NOT single"))
exp.check("P6", all(r["single_c"] for r in mres + ares),
          f"{sum(r['single_c'] for r in mres + ares)}/{len(mres) + len(ares)} "
          f"peaks live in {{I, X_c}} for one qubit c")

exp.section("P7  does the fold land exactly on the perm PEAK set?")
p7 = []
for r in mres + ares:
    want = {z: (2 if (z >> r["c"]) & 1 else 1) for z in r["S"]}
    ok = r["fold"] == want
    p7.append(ok)
    exp.log(f"{r['label']:>18}: fold has {len(r['fold']):,} keys, perm peak "
            f"{len(r['S']):,}; multiplicities {'match' if ok else 'DIFFER'}; "
            f"2|S|-|B| = {2 * len(r['S']) - len(r['B']):,} vs rot_max "
            f"{r['rot_max']:,}")
exp.check("P7", all(p7) and all(r["rot_max"] == 2 * len(r["S"]) - len(r["B"])
                                for r in mres + ares),
          f"{sum(p7)}/{len(p7)} instances: the peak IS the perm peak with the "
          f"z_c=1 half doubled, so the deficit is exactly |B|")

exp.section("P8  what are the two strings?")
for r in mres:
    names = qubit_names(r["me"])
    for z in r["B"]:
        exp.log(f"{r['label']:>18}: z = {z:>7} = {render(z, names, r['q'])}"
                f"   [gadget target c = {names.get(r['c'], r['c'])}]")
shapes = []
for r in mres:
    me = r["me"]
    obs, e0 = 1 << me.x[0], 1 << me.exp[0]
    shapes.append(sorted(r["B"]) == sorted([obs, obs | e0]))
exp.check("P8", all(shapes) and all(len(r["B"]) == 0 for r in ares),
          f"{sum(shapes)}/{len(shapes)} modexp instances have "
          f"B = {{Z_x0, Z_x0 Z_e0}} exactly; adders have B empty")

exp.section("P9  is B the set of |coefficient| = 1/2 Walsh terms?")
p9 = []
for r in mres:
    me = r["me"]
    coeffs = walsh.pullback_coefficients(me.build(), me.x[0])
    half = set(np.nonzero(np.abs(np.abs(coeffs) - 0.5) < 1e-12)[0].tolist())
    ok = half == set(r["B"])
    p9.append(ok)
    top = np.argsort(-np.abs(coeffs))[:4]
    names = qubit_names(me)
    exp.log(f"{r['label']:>18}: |c|=1/2 terms {sorted(half)} vs B {r['B']} "
            f"-> {'SAME' if ok else 'DIFFER'};  top-4 |c| = "
            + ", ".join(f"{render(int(z), names, r['q'])}:{abs(coeffs[z]):.4f}"
                        for z in top))
    del coeffs
exp.check("P9", all(p9),
          f"{sum(p9)}/{len(p9)} instances: the strings that fail to double are "
          f"exactly the dominant Fourier modes")

# ---------------------------------------------------------------------------
exp.section("C1  control -- adders")
for r in ares:
    exp.log(f"{r['label']:>18}: perm {r['perm_max']:,}, rot {r['rot_max']:,}, "
            f"deficit {r['deficit']}, |B| = {len(r['B'])}")
exp.fail_check("C1", all(r["deficit"] == 0 and len(r["B"]) == 0 for r in ares),
               f"adders: deficits {[r['deficit'] for r in ares]}, "
               f"|B| {[len(r['B']) for r in ares]}")

exp.section("C2  control -- same circuit, observable moved into the scratch")
header()
c2 = []
me = ToffoliModExp(N=5, a=2, n_exp=1)
qc = me.build()
names = qubit_names(me)
for q in (me.b[0], me.t[0], me.anc):
    r = analyse("N=5 a=2", qc, q)
    r["me"] = me
    c2.append(r)
    show(r, names[q])
    if r["c"] is not None:
        head = ", ".join(render(z, names, r["q"]) for z in r["B"][:4])
        exp.log(f"     |B| = {len(r['B']):,}"
                + (f"; first few: {head}{' ...' if len(r['B']) > 4 else ''}"
                   if r["B"] else " (empty)"))
exp.fail_check("C2", any(r["deficit"] != 2 for r in c2),
               f"deficits by observable: "
               + ", ".join(f"{names[q]}:{r['deficit']}"
                           for q, r in zip((me.b[0], me.t[0], me.anc), c2))
               + " -- if these were all 2 the |B| story would be decoration")

exp.finish()
