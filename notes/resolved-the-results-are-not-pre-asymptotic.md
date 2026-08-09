---
code: resolved-the-results-are-not-pre-asymptotic
title: "C7 RESOLVED — the results are NOT pre-asymptotic"
outcome: solved
claims: [C7]
todo: []
---
# C7 RESOLVED — the results are NOT pre-asymptotic

F9 removed the need for a faster simulator: since PPS term count *equals* Walsh
sparsity, exact PPS cost is computable in O(2ⁿn) without running PPS. That
reached **24 qubits** (dim 1.7e7) where PPS itself stalls around 15–17.

`experiment_c7.py` / `/tmp/c7a.log`, Toffoli modexp, observable Z_x0:

```
   N   n   r  qubits        2^q   sparsity  density  r pow2    time
   5   3   4      15      32768      15493  0.47281     yes    0.1s
   7   3   6      15      32768      15539  0.47421      no    0.0s
  15   4   4      18     262144     127936  0.48804     yes    0.8s
  21   5   6      21    2097152    1037322  0.49463      no   16.5s
  33   6  10      24   16777216    8347241  0.49753      no  258.2s
  35   6  12      24   16777216    8346759  0.49751      no  262.1s

  log2(sparsity) grows 1.008 bits/qubit  (1.000 = exactly Theta(2^n))
  density: 0.473 -> 0.474 -> 0.488 -> 0.495 -> 0.498 -> 0.498
```

**Density converges monotonically to 1/2 and the growth slope is 1.008
bits/qubit.** So PPS on modular exponentiation costs Θ(2ⁿ) — asymptotically no
better than a state-vector simulation. The small-n numbers were already in the
asymptotic regime; C7 is answered, and negatively for the threat.

For contrast, PPS at 17 qubits ran >20 min without finishing while the Walsh
computation of the same quantity took **0.01s**.
