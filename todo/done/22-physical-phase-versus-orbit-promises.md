---
id: 22
state: done
title: "Does a physical work-local phase escape the orbit-index shortcuts?"
outcome: "bounded physical test complete: pi-sector routing and period-six regrouping preserve cheap shortcuts; generic coarse dephasing and naive exponent transfer fail on the fixed fixture"
claims: [C53, C56, C57, C58]
---
# A natural gate with extensive orbit support

Completed 2026-09-11 UTC at bounded first-pass scope. C58 owns the exact
control-independent routing criterion and small-orbit identities; §SR records
independent physical tests on both moduli, failures and main audit. This does
not solve arbitrary physical phases or establish hardness. The original
question below is retained; TODO 23 owns the next research step.

C57's localized kick genuinely breaks coarse-sector dephasing but has an
omission guarantee once its reachable support is small. That exact sampler
question remains open; the complete-output boundary identity does not settle
its low-bit marginals. Keep that open question explicit rather than calling
the finite tests a proof of hardness or quietly labeling the approximation
as an exact simulator.

The next discriminating family changes the support promise in a physically
simple way: replace the one-pair kick by Rz(theta) on one actual modular work
qubit, at the SAME insertion inside the same periodic background. A diagonal
physical phase cannot leak out of the clean orbit, yet K-I typically acts on
the entire orbit. Thus few-physical-qubit locality does not imply C55/C57's
finite orbit-basis support. Its orbit displacement is zero, but that alone
does not make the extensive-support omission bound useful.

## Smallest tests and stronger baselines

Start with the fixed N=7,a=3,t=5, W01@s1,W12@s3,W01@s4 fixture, phase on work
bit 1 after s=2. Vary only theta, including zero and both signs. Use the
actual Circuit/statevec reference and existing full-r `sequential_path`, with
known order/index and construction costs explicit. Ask whether wrong initial
coarse dephasing changes measured probabilities; noncommutation alone is not
evidence. Controls: zero angle, end insertion, a genuinely periodic orbit
diagonal phase, and deleting a detected coherent sector contribution.

Before scaling, derive the ACTUAL orbit phase sequence

    g(j)=exp[-i*theta*(-1)^bit(a^j mod N)/2]

and its action on the coarse sectors. On a tiny orbit, a complement symmetry
or short Fourier expansion may make this another small exact shortcut. Check
that possibility before celebrating the failure of the old mixture. C53's
input-specific exponent-phase transfer also remains a required comparator;
with earlier mixers it cannot simply be assumed valid. Preserve original gate
order throughout. Compare the unperturbed sampler at the same output accuracy.

If the N=7 fixture has a special collapse, identify its algebra explicitly and
test the predicted scope on one independently simulable next modulus, with a
preflight state-memory and gate/runtime budget. Do not turn two different
moduli into an unsupported scaling law. Do not hide orbit enumeration in
constructing an apparently small transition matrix. Initial testing should
again go to lower-cost agents while main audits the algebra and the strongest
known baselines.

## Decision point

The useful question is which description controls prediction cost: physical
gate locality, translation symmetry, Fourier coupling bandwidth, or a more
task-specific sufficient statistic. A finite observable discrepancy is not a
simulation lower bound. A short sector-coupling formula is not yet a sampler;
it must give normalized conditional outputs without an uncharged all-sector
sum. If another exact small shortcut wins, record its promise and avoid a
larger toy benchmark. If the structure becomes genuinely different, consult
current primary sources for that specific family before choosing a framework
or asserting novelty. Precision certification remains a separate unresolved
requirement for every existing floating-point sampler.
