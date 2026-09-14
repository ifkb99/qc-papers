---
id: 34
state: open
title: "Resolve recurring native crashes before trusting long sampling comparisons"
outcome: "Two word-envelope sweeps aborted at different locations; old processor microcode is a concrete host-stability lead, not a proven root cause"
claims: [C71, C72, C73, C74]
---
# Reliability before another long sweep

CF records a further terminal native failure during the full legacy
permutation suite in the CNOT-frame investigation. The raw status and empty
buffered log are retained there; that gate remains incomplete. Short core,
law and allocation passes in CF do not resolve host reliability or replace
the failed full suite.

WE owns the failed-run evidence and partial checkpoint; PC owns the independent
exact-arithmetic follow-up. TODO 32's full serial timing comparison remains
unfinished. Do not combine selected successful rows from aborted runs into a
claimed complete benchmark, or infer a scientific failure from a native abort.

## Read-only diagnostic findings, 2026-09-11

The host reports Intel Core i9-14900K, microcode 0x11d, ASUS ProArt
Z790-CREATOR WIFI, BIOS 1501 dated 10/06/2023. Intel's current
[stability guidance](https://www.intel.com/content/www/us/en/support/articles/000102331/processors.html)
recommends Intel Default Settings and BIOS microcode 0x12F or later. This
mismatch merits attention but does NOT prove CPU degradation or explain a
specific NumPy/Python fault. Software, memory and other causes remain open.
Several earlier kernel fault entries name logical CPU 10; no matched affinity
experiment or hardware diagnosis has been completed.

Main asked the user about other crashes and BIOS history. Do not flash firmware,
change voltage/power/overclock settings, reboot, install system microcode or
run a prolonged stress test without explicit direction. No such changes have
been made. Do not claim the existing Python-3.14 warning explains these
separate Python-3.12 crashes.

## Next bounded checks

1. Preserve each crash trace, exact versions, native build metadata and kernel
   log evidence. A passing isolated row is not a fix. The comparison now saves
   per-row incomplete JSON checkpoints outside its measured intervals.
2. The short system-Python core/law/preflight checks now pass; WE records them
   as a separate environment, not pooled timing. Change one runtime component at
   a time if further isolation is useful; no uncontrolled full-sweep retries.
3. With the user's hardware context, choose between a bounded software/CPU-
   affinity reproduction, a second host, or user-managed firmware/support
   investigation. Treat affinity as a diagnostic, not a correctness certificate.
4. Only after a defensible stable-run plan, repeat TODO 32 serially with frozen
   inputs, caps, counters and all failures retained. Maintain one complete
   environment per timing comparison. Numerical certification remains TODO 24.
