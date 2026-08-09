---
id: 8
state: done
title: Decode the dominant coefficients
outcome: DONE inside step 2
claims: []
---
# Decode the dominant coefficients

For modexp N=5 the four dominant coefficients are `{x8,e13}` and `{x8}` at
|c| = 0.5 (weights 2 and 1), plus a pair at |c| = 0.1465 (weights 8 and 9). They
sum to exactly ⟨O⟩; the other 3082 sum to exactly 0. See `NOTES.md` §W3–W4.

Still open: the two large coefficients are supported almost entirely on the
x-register bit being measured plus one exponent qubit. Whether that generalises
across N, a and observable is untested and would be cheap to check.

---
