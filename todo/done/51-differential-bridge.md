---
id: 51
state: done
title: "Does differential cryptanalysis count off-diagonal Pauli paths the way linear cryptanalysis counts diagonal ones?"
outcome: "C99/DB prove the DDT identity, bounds and APN extremality for X/Y-type pullbacks; independent dense, gate-level PPS and 14-qubit GPU checks pass, including the out-of-sample F10 counts; no simulation speedup follows"
claims: [C8, C12, C13, C25, C99]
---
# The differential bridge for off-diagonal observables

A Claude-originated question, raised on 2026-09-14 at the user's invitation
to propose an idea. C13 and F10 leave X/Y-type observables outside the Walsh
identity; F10 reports only a cap hit. The question was whether the
cryptanalysis dictionary of C12/C25 has a differential half covering exactly
those observables, with a derived count, bound and extremal family.

Completed at bounded scope. C99 owns the identity, counting theorem, bounds
and scope limits. DB owns the derivation history, fixtures, controls and the
out-of-sample F10 check. No follow-up TODO was opened; DB lists candidate
directions for the user to decide on.
