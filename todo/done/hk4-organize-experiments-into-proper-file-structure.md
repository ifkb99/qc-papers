---
id: hk4
state: done
title: Organize experiments into a proper file structure, with a general
claims: []
---
# Organize experiments into a proper file structure, with a general

Organize experiments into a proper file structure, with a general
engine — **DONE 2026-08-08.** `lab/` package (harness with enforced
predictions/controls, cached measurement, GF(2) structure tools, modarith,
wrap-registry variants, null models); experiments moved to `experiments/`
unchanged (records, not code to DRY); `experiment3.py` to `archive/`; two
new gate suites (`test_lab.py` pins the engine to logged numbers,
`test_claims.py` re-verifies the headline claims); `experiments/TEMPLATE.py`
for new work. Deliberately NOT a `src/` package — the unit of value here is
the readable record, and root-level instruments keep every historical
import working.
