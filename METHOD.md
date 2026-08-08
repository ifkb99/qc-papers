# Method notes

## Is this different from the scientific method?

No. It *is* the scientific method — conjecture, predict, test, revise. Nothing
here is a new epistemology and it would be silly to dress it up as one.

What *is* worth writing down is that computational/mathematical research inverts
one prior badly, and that inversion drives most of the technique below.

**In empirical science, a surprising result might be a discovery. In
computational research on your own code, a surprising result is almost always
your own bug.** Nature is not being subtle; you made a mistake. Every anomaly
starts guilty.

Three further differences from lab science shape everything:

1. **Experiments are exact and deterministic.** No noise, no statistics, no
   p-values. `other = 0 in 7/7` is categorical. Falsification is sharp — and a
   single counterexample is fatal rather than "anomalous".
2. **Ground truth is usually independently computable.** The same quantity can
   often be obtained by a completely unrelated route. That is far stronger than
   more tests of the same code.
3. **Empirical regularities can be promoted to theorems.** This is the real
   difference from lab science. If you see an invariance, try to *prove* it —
   and the proof usually comes from reading the construction, not from more
   measurements.

---

## The loop, as actually run here

```
   observe something odd
        |
        v
   [BUG CHECK]  <-- assume it is your fault; discharge this FIRST
        |
        v
   conjecture a mechanism
        |
        v
   derive a consequence you have NOT yet measured
        |
        v
   test that consequence
        |
    +---+---+
    |       |
  fails   holds
    |       |
    v       v
 retract  try to promote to a proof (read the construction)
 + log      |
            v
     state scope limits, then move on
```

The two steps people skip are the bug check and *deriving before measuring*.
Both earned their place here expensively.

---

## Techniques that paid for themselves

**Prior art first.** Cheap, decisive, and it can invalidate everything
downstream. Ran it as step 1 before any writing. It narrowed the novelty claim
substantially and cost hours instead of a referee report.

**Derive-then-test.** Predict a consequence *before* looking. Used twice:
branch symmetry `|ĥ|=|ĝ|`, and support confinement `z_I ∈ {0,1_I}`. Both
confirmed with zero violations. Evidence from a prediction is worth far more
than an explanation fitted to data already seen — and it is the only way to
avoid dressing up a curve fit as a mechanism.

**Every experiment carries a control that must FAIL.** The one rule it is most
tempting to skip once a hypothesis is going well, and it has now caught two
vacuous measurements of mine that had already produced confident-looking
numbers:

- *step 9* — synthetic blocks acted only on b-qubits while the observable was
  `Z_x0`, so both block types came out constant for a trivial reason;
- *TODO 12e* — the tail-confinement test `z_I ∈ {0, 1_I}` is vacuous at
  |I| = 1, because one bit *is* all-zeros or all-ones. The β>1 control passed
  when it had to fail, which is the only reason it was noticed — and three
  rows of the positive result were being read as evidence too.

In both cases nothing measured was wrong. What was wrong was how much the
measurement was entitled to say, and only the control could tell the
difference. `lab.harness` warns when no must-fail control is registered.

(The project's other two self-caught errors came from different rules, not
this one: "vary exactly one parameter" caught a sweep that drew a new random
table per t, and null-model hygiene caught random tables degenerating to
constants at small r. Worth keeping straight — the rules are not
interchangeable.)

**Precision sweep to separate bug from float error.** Re-run in `longdouble`.
If the error is *identical*, it is a logic bug, full stop. This caught the θ=π
bug that six test suites had missed.

**An independent exact reference beats more tests.** The Walsh transform
computed the same quantity by a completely different route and immediately
exposed a propagator bug that every existing suite passed. Build the second
route early.

**Vary exactly one parameter.** Fix the modulus, vary the base. An
ancilla-count confound (10 qubits vs 15) had previously killed an entire thesis
that looked well-supported.

**Step by 1, not 2.** The `n_exp = v₂(r)+1` threshold was invisible while
stepping n_exp in twos — and the original sweep happened to *start* exactly on
the threshold, so the claim looked stronger than it was. Pure luck, in the
dangerous direction.

**Randomised tests only cover what they sample.** The random-circuit suite drew
from `{h,t,cnot,rx,rz}` and so never emitted an `X`. Enumerate the gate *set*;
do not sample it and call it coverage.

**Check that the measured quantity is the quantity that matters.** "Exact cost
model" predicted `N_final` while PPS memory is `N_max`. The identity was right;
the claim was about the wrong number.

**Read the construction when experiments stall.** Two experimental routes to the
C15 proof failed. The proof came from reading four lines of the circuit
definition and noticing `u_a(·,1) = A⁻¹SA` — a conjugate of an involution. No
measurement would have produced that.

---

## Failure modes seen here, with names

| failure | instance |
|---|---|
| **Right conclusion, wrong mechanism** | adder collapses to 1 term — because permutation (wrong) vs because affine (right). Still counts as a failure. |
| **Toy that isn't a proxy** | H → adder → QFT "Shor sandwich": the observable barely met the arithmetic. |
| **Confound in a matched comparison** | compilation A/B with mismatched qubit counts. |
| **Bug agreeing with a hypothesis** | the θ=π bug produced numbers consistent with the compilation thesis, which is why it survived. |
| **Lucky sampling hiding a threshold** | every sweep started at n_exp = α+1. |
| **Trusting a summary over source** | web summaries claimed Qiskit `pauli-prop` supports Toffoli; the source rejects it. |
| **Reading an abstract, not the body** | Dang et al.'s abstract says "factors of r"; §4 is explicitly 2-adic. Nearly caused a wrong retraction. |
| **A test with only one possible answer** | `z_I ∈ {0, 1_I}` at \|I\| = 1. Caught by the must-fail control, which passed. See §OS3. |

---

## Bookkeeping that mattered

- **A retraction log is worth more than a results list.** Half of what was
  believed at various points is now marked dead. Without the log, a future
  session re-derives dead ends.
- **Write for a reader with no context** — including yourself in a later
  session. `NOTES.md` is context restoration, not a report.
- **Mark stale files in the file itself**, not only in the index. Retracted
  experiments carry warning headers because an index entry will be missed.
- **Record *why* something was believed**, not just that it was wrong. "Killed
  by the θ=π bug" is more useful than "wrong".
