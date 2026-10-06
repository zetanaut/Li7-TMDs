# Li7-TMDs

This site explains the executable checks of spin-3/2 target-polarization algebra, leading-twist quark and gluon tensor bases, a selected electromagnetic SIDIS trace, and a Born-level gluon hard response. The motivating target is polarized lithium-7. The algebraic target-spin basis is general to spin 3/2.

The run-scoped `baseline` profile verifies the historical calculations. The cumulative profiles add independent finite-dimensional, process, and limit checks. The `full` profile executes the declared finite computational program and the [convention closure review](convention-closure.md). The author-approved Born index correction is included. Deployment requires separate exact-commit approval after full validation.

The independent coefficient count is **32 quark and 32 gluon coefficients** per specified flavor and gauge-link/color class. In either sector the ranks K = 0, 1, 2, 3 contribute 2, 6, 10, and 14 coefficients. The new spin-3/2 sector is the rank-3 octupole.

## Recommended reading/running order

1. [Physics overview](physics-overview.md)
2. [Spin-3/2 polarization](spin32-polarization.md)
3. [Two-dimensional transverse STF tensors](transverse-stf.md)
4. Run `python src/validate_spin32_tmd.py`
5. [Quark basis](quark-basis.md)
6. [Gluon basis](gluon-basis.md)
7. [SIDIS trace validation](sidis-validation.md)
8. Run `python src/gluon_born_response.py`
9. [Gluon Born response](gluon-hard-response.md)
10. Inspect the [plots and JSON reports](reproducing-results.md)

A **spin-density multipole** is a part of the target's 4 × 4 spin-density matrix that transforms with a definite rotational rank. An **STF tensor** is symmetric under exchange of its indices and has zero contraction (trace) over any index pair. A **TMD** is a scalar coefficient in a correlator that retains partonic transverse momentum.
