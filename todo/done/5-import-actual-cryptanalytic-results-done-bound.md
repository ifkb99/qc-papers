---
id: 5
state: done
title: Import actual cryptanalytic results
outcome: DONE, bound is real but weak
claims: [C12, C25, C26]
---
# Import actual cryptanalytic results

See `NOTES.md` §X; file `experiment_crypto.py`.

- **C25.** Derived and verified **S ≥ (1 − NL/2ⁿ⁻¹)⁻²** from Parseval. Converts
  any published nonlinearity into a PPS cost lower bound for *every* circuit
  computing the function, with no simulation. Tight at both extremes.
- **Validation.** AES S-box nonlinearity comes out at exactly **112** over all
  255 nonzero linear combinations — the published constant. An external check on
  the whole Walsh pipeline.
- **End-to-end.** Built reversible circuits computing the inner product (bent)
  and propagated them: support exactly 2^(2m) (64/256/1024), matching the bent
  prediction. No truncation is available when every coefficient has the same
  magnitude — the clean worst-case statement.
- **C26, the honest limit.** Loose away from the extremes: AES bound 64 vs 239
  actual, modexp 4 vs 3086. NL uses only `max|c|` and throws away the rest of the
  spectrum. Follow-up worth doing: for crypto families whose **full** Walsh
  value/multiplicity distribution is published (the AES inverse among them), S is
  determined *exactly* rather than bounded — a much stronger import.
- **C12 correction.** The "0.74 of the bent bound" figure is stable across N at
  fixed n_exp, **not** across widths (0.50 at n_exp=1). Fixed in the abstract.
