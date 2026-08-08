"""TODO step 4 -- third modulus for C15, plus a structural question it raises.

Two things.

(1) N=5 was never swept across n_exp. It is worth more than "a third data
    point": N=5 with a=2 has r=4, i.e. alpha=2, whereas every constancy sweep so
    far used r=2 (alpha=1). If the invariance is really about beta=1 and not
    about alpha, both should be constant -- at different constant values.

(2) While setting it up: N=5 and N=15 have NO base with an odd-factor order.
    lambda(5)=4 and lambda(15)=lcm(2,4)=4, so every order divides 4 and is a
    power of two. That would mean those moduli are in the free branch for EVERY
    base, not just for the usual demo base -- a much stronger statement than C20
    currently makes. Check it, and characterise which N have the property.
"""
from __future__ import annotations
import numpy as np, time, math
from toffoli_arith import ToffoliModExp
import walsh


def order(a, N):
    r, v = 1, a % N
    while v != 1:
        v = (v * a) % N
        r += 1
    return r


def v2(r):
    k = 0
    while r % 2 == 0:
        r //= 2
        k += 1
    return k


def carmichael(N):
    """lambda(N) via lcm of orders (N is small here)."""
    units = [a for a in range(1, N) if math.gcd(a, N) == 1]
    lam = 1
    for a in units:
        lam = lam * order(a, N) // math.gcd(lam, order(a, N))
    return lam


print("=" * 78)
print("1. THIRD MODULUS: N=5 swept across n_exp  (r=4, i.e. alpha=2)")
print("=" * 78)
print("  Prior sweeps used r=2 (alpha=1). If beta=1 is what matters, r=4 should")
print("  also be constant -- at a different constant.\n")
print(f"  {'a':>3} {'r':>3} {'alpha':>6} {'beta':>5} {'n_exp':>6} {'q':>4} "
      f"{'|support|':>10} {'density':>9} {'vs first':>9} {'t':>7}")
for a in (2, 4):
    r = order(a, 5)
    base = None
    for ne in (2, 4, 6, 8):
        me = ToffoliModExp(N=5, a=a, n_exp=ne)
        qc = me.build()
        t0 = time.time()
        c = walsh.pullback_coefficients(qc, me.x[0])
        el = time.time() - t0
        sp = int((np.abs(c) > 1e-12).sum())
        dim = 1 << me.n_qubits
        if base is None:
            base = sp
        tag = "SAME" if sp == base else f"{sp - base:+d}"
        print(f"  {a:3d} {r:3d} {v2(r):6d} {r >> v2(r):5d} {ne:6d} {me.n_qubits:4d} "
              f"{sp:10d} {sp/dim:9.6f} {tag:>9} {el:6.1f}s", flush=True)
        del c
    print()

print("=" * 78)
print("2. WHICH MODULI ARE FREE FOR EVERY BASE?")
print("=" * 78)
print("  lambda(N) is the Carmichael function; every order divides it. If")
print("  lambda(N) is a power of two then EVERY base has beta=1.\n")
print(f"  {'N':>5} {'factors':>12} {'lambda(N)':>10} {'pow2?':>6} "
      f"{'orders present':>28}")
for N in (5, 7, 15, 21, 33, 35, 51, 85, 91, 143, 255):
    lam = carmichael(N)
    orders = sorted({order(a, N) for a in range(1, N) if math.gcd(a, N) == 1})
    f = []
    m = N
    for p in range(2, N + 1):
        while m % p == 0:
            f.append(p); m //= p
    ispow2 = (lam & (lam - 1)) == 0
    print(f"  {N:5d} {'x'.join(map(str, f)):>12} {lam:10d} "
          f"{'YES' if ispow2 else 'no':>6} {str(orders):>28}")

print("""
  lambda(N) is a power of two exactly when every prime power in N contributes a
  power of two, i.e. N = 2^a * (product of DISTINCT FERMAT PRIMES).
  Known Fermat primes: 3, 5, 17, 257, 65537.
  So the odd semiprimes in the free branch are exactly p*q with both p and q
  Fermat primes: 15=3x5, 51=3x17, 85=5x17, ... and 15 is the smallest.""")

print()
print("=" * 78)
print("3. CONSEQUENCE FOR THE STANDARD DEMO INSTANCE")
print("=" * 78)
for N in (15, 21):
    units = [a for a in range(2, N) if math.gcd(a, N) == 1]
    ords = {a: order(a, N) for a in units}
    free = [a for a, r in ords.items() if (r >> v2(r)) == 1]
    print(f"  N={N}: {len(units)} usable bases, "
          f"{len(free)} with beta=1 (free branch) -> {100*len(free)/len(units):.0f}%")
    print(f"        orders: {dict(sorted(ords.items()))}")
