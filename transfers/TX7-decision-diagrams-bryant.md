---
id: TX7
field: "Binary decision diagrams (Bryant)"
status: imported
effect: order-dependent
one_line: "Adder carry chains have linear ROBDDs in interleaved order; per-block modexp ROBDDs measured far below Walsh"
source: "Bryant, ACM Computing Surveys 24(3) 1992, §1.3–1.4, Table 1, Fig. 4 (read in C102)"
claims: [C102, C103, C106]
notes: [HD, HB, DN]
todo: [13, 59, 68]
---
# TX7 — Binary decision diagrams (Bryant)

## Dictionary

ROBDD of the pulled-back bit ↔ an exact compressed PPS output.

## Hypotheses

Order-sensitive; C103's bound holds in exponent-first order only.

## Consequence for the goal

Compresses the final object; DN shows propagation-native diagrams lose in bytes.
Asymptotics at growing n are bounded by TX8, which must be checked first.
