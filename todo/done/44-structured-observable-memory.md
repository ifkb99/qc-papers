---
id: 44
state: done
title: "Can exact Boolean-observable structure reduce CNOT-frame storage beyond sparse Walsh enumeration?"
outcome: "C83/QC prove and implement exact quadratic-cell closure; the stronger recognized formula wins and broad novelty is excluded"
claims: [C8, C17, C45, C82, C83]
---
# Compress the represented observable, not only its coordinate updates

Completed at mathematical and bounded exact-query scope. C83 owns closure,
resource bounds and the implementation contract; QC owns independent evidence,
failed verifiers, the stronger-baseline loss and validation. TODO45 owns the
remaining phase-sensitive stabilizer comparison. The original brief follows.

C82/CF resolve the allocation-only discriminator. Do not repeat its depth
sweep or call another lazy-frame/active-span implementation a breakthrough.
Clifft and classical Clifford/quadratic-form methods already cover broad
dynamic-coordinate and low-active-dimension ideas. A new result needs a
specified prediction task and a smaller sufficient representation that survives
the actual nonlinear circuit, with preprocessing and output extraction charged.

First bounded candidate: pullback signs (-1)^f(x) with f a quadratic Boolean
polynomial. A product of disjoint nonlinear target updates can have exponential
Walsh support while admitting compact quadratic data. That fact is known
stabilizer/quadratic-form mathematics. Audit the exact closure boundary under
Toffoli and CNOT, and compare it with a stabilizer-state representation of the
normalized sign vector and a direct GF(2) quadratic Gauss-sum reference.

An unmeasured derivation to audit is substitution at Toffoli target t:
f(x+t*a*b)=f(x)+a*b*derivative_t(f). For a quadratic f, genuinely cubic
terms arise from target neighbors outside the two control variables; repeated
Boolean variables must be reduced using x^2=x. Determine the exact necessary
and sufficient closure condition, including linear terms, before coding.
Preserve a cubic witness that must FAIL the quadratic promise. Do not silently
drop higher-degree terms or confuse final quadratic degree with intermediate
degree. Non-Clifford gate count alone is not the resource parameter.

Freeze n and vary one violating gate or structural parameter. Requested
outputs can be a selected physical Walsh coefficient, an expectation on a
specified product input, or a streamed/sampled coefficient law with justified
normalization; producing an explicit exponential list still costs its length.
State the actual useful task before measurement. Existing sparse PPS and the
strongest compact classical reference must receive the same inputs. If the
candidate only rediscovers standard stabilizer evaluation, retain the exact
boundary as negative guidance and redirect to an extension that makes an
additional testable prediction. Primary prior art and a bounded independent
control come before a simulator or a scaling claim.

The missing legacy full-suite validation remains with TODO34; short bounded
science can continue under its existing constraints, with core first. Use the
board swarm with one coordinator source writer and independent mathematical
and technical review. Existing TODO42 remains open on its separate sampler
cost-selection question.
