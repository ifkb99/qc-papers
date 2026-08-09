---
code: caveats-that-decide-whether-real
title: "Caveats that decide whether F3 is real"
outcome: record
claims: []
todo: []
---
# Caveats that decide whether F3 is real

- **Scale gap is severe.** Mine: 6–14 qubits, 31–129 T gates. Paper: 127 qubits,
  5000–8000 gates. Paper says the power law only emerges after ~1/3 of the
  circuit. My circuits may simply be **pre-asymptotic**, which would make F3 an
  artifact. This is the single biggest threat to the result.
- Python dict-based PPS dies around ~15 qubits. Paper's **Appendix B** gives the
  bit-packed representation (ν_P vectors in uint64 arrays) — that's the fix if
  scale is needed.
- `<O>=0` degeneracy above means F3's cliff was measured on a trivial observable.

---
