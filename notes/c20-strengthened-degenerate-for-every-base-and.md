---
code: c20-strengthened-degenerate-for-every-base-and
title: "C20 STRENGTHENED — N=15 is degenerate for EVERY base, and necessarily so"
outcome: record
claims: [C20]
todo: []
---
# C20 STRENGTHENED — N=15 is degenerate for EVERY base, and necessarily so

Every order divides the Carmichael function λ(N), so if λ(N) is a power of two
then **every** base has β=1 and the whole modulus is in the free branch.

```
    N   factors  lambda(N)  pow2?          orders present
    5         5          4    YES               [1, 2, 4]
   15       3x5          4    YES               [1, 2, 4]
   21       3x7          6     no            [1, 2, 3, 6]
   33      3x11         10     no           [1, 2, 5, 10]
   51      3x17         16    YES        [1, 2, 4, 8, 16]
   85      5x17         16    YES        [1, 2, 4, 8, 16]
  143     11x13         60     no   [1,2,3,4,5,6,10,12,15,20,30,60]
```

```
  N=15: 7 usable bases, 7 with beta=1  -> 100%
  N=21: 11 usable bases, 3 with beta=1 ->  27%
```

λ(N) is a power of two exactly when **N = 2^a × (product of distinct Fermat
primes)** — since p−1 must be a power of two for each odd prime p, and p^k needs
k=1. Known Fermat primes: 3, 5, 17, 257, 65537.

**So the odd semiprimes in the free branch are exactly p·q with both p and q
Fermat primes: 15 = 3×5, 51 = 3×17, 85 = 5×17, … and 15 is the smallest.**

That is a sharper version of C20 than "N=15, a=7 happens to have r=4". The
canonical demonstration instance is degenerate **for every base**, and it is
degenerate *because* it is the smallest product of two Fermat primes — the same
property that makes it the natural smallest demo. The degeneracy is forced by
the choice of N, not by the choice of a.
