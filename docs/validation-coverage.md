# Validation coverage

The `baseline` profile executes the original finite-algebra and Born calculations with run-scoped evidence. It is **not** the full validation of *Octupole parton structure of spin-3/2 targets*. The `full` profile returns `MISSING` until its required scientific suites are implemented. Source labels below refer to the current submission-review definitions. A `PASS` applies only to the stated test scope.

| Claim ID | Source label | Baseline status | Scope and next evidence |
|---|---|---|---|
| `spin.inverse` | `eq:Q`, `eq:QO`, `eq:rho` | PASS, partial | The legacy exact density reconstruction runs on 16 Hermitian coordinates. Independent ladder/CG constructions, trace metrics, rotations, positivity, and physical preparation tests are MISSING. |
| `stf.branches` | `eq:Kn`, `eq:Kcomplex`, `eq:rankrule` | PASS, partial | Complex-component STF and branch identities run. Independent Cartesian construction and symbolic mass normalization are MISSING. |
| `gluon.dyadic` | `eq:referenceSTF`, `eq:dyadic_not_STF` | MISSING | Existing related null tensors do not test the source's STF-versus-dyadic distinction directly. |
| `basis.independence` | `eq:Ffull`, `eq:Gfull`, `eq:rankrule` | PASS, partial | The implemented quark/gluon maps each have computed rank 32 at one momentum. Independent covariant review, rank witnesses, spanning, and inverse projectors are MISSING. |
| `basis.name_mapping` | `eq:Ffull`, `eq:Gfull` | BLOCKED | The legacy `f/g/h` catalogue lacks flavor, antiquark, link/color, and hadron fields. A one-to-one mapping to every named manuscript coefficient requires the complete correlator review. |
| `sidis.response` | `eq:SIDISmaster` | PASS, partial | Twelve hard-trace pairs and a reduced contraction run. The 32-row catalogue and independent convolution are MISSING. |
| `dy.response` | `eq:DYconv` | MISSING | No DY projection, antiquark, angular, or convolution suite exists. |
| `link.reversal` | `eq:PTspin`, `eq:octSIDISDY` | MISSING | No operator-derived sign or Wilson-line descriptor suite exists. |
| `gluon.born` | `eq:BornB` | PASS, partial | Five diagnostics run at the specified point and 36 grid cases. Independent spinor/high-precision checks and the dense scan are MISSING. |
| `gluon.separation` | `eq:gkernelrule`, `eq:sevendirections` | MISSING | The target response and preparation reconstruction are not implemented. |
| `qcd.factorization` | process-dependent factorization definitions | ANALYTIC_ONLY | Finite algebra and Born diagnostics cannot prove factorization, evolution, or universality of measured asymmetries. |
| `nuclear.magnitudes` | nuclear interpretation | ANALYTIC_ONLY | No lithium-7 wave function, distribution fit, rate, or sensitivity prediction is computed. |

The [convention map](conventions.md) explains which source definitions the existing code uses. The [reproduction instructions](reproducing-results.md) explain how to obtain and validate a fresh baseline run. Expected missing suites are explicit in `validation_manifest.FULL_PENDING_SUITES`; they are not fabricated passing results.
