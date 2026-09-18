---
id: 59
state: open
title: "Promote the DD-native PPS pilot and its proposed module to the tracked tree so note DN's numbers are reproducible from git"
claims: [C104]
---
# Promote the DD-native PPS module (the DD analogue of TODO 58)

Opened 2026-09-17 from the DN integration review (`V1383c571ed7f471f`). The accepted
promotion round (`T163d291c2cdd4cd0`, submission `S4279919efc104d37`, review
`V0c95e8735e50477a`) produced a proposed `lab` module and experiment, and they were
never integrated: the coordinator held the integration (`M889892e48cad42b9`) because the
round's registered P5 trend clause is contradicted by the accepted six-point series in
DN. So `lab/dd_pps.py` and `experiments/experiment_dd_pps.py` do not exist, the frozen
pilot and every reading live under gitignored `out/`, and nothing in git reproduces DN.

What this needs, if the direction is ever revived: the accepted module and experiment
integrated with a reconciled prediction (a v5 registering the real trend, not P5's), the
print-only validator replaced by one that exits nonzero (`M841349bedd124e52`), and DN
updated to cite the tracked script. The same gap for K2 is TODO 58.

Priority is low: DN closes this direction negatively for the peak-memory goal, and C104
supersedes it under a streamed output. Do not revive it without a reason that survives
DN's byte figures.
