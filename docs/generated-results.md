# Generated result summary

This page describes the **quark-processes profile only**. The complete manuscript validation profile has required suites that are still missing.

The summary was generated from a complete run whose required IDs, source and input digests, computed ranks, and result identities were checked together.

Source revision: `fac802c1323562849f4524406b6d90436ce1205f`; dirty source: `false`; scientific source digest: `fa847dd28c6a426fe13badef5f73b9d970986edfa069fe2cd157c50ed236359f`.

| Executed evidence | Count |
|---|---:|
| conditional algebra | 66 |
| exact fixture | 46 |
| exact identity | 178 |
| exact rank | 6 |
| numerical diagnostic | 52 |
| numerical parameter cases | 1 |
| software test | 22 |

| Evidence role | Records |
|---|---:|
| component case | 122 |
| independent bound | 2 |
| independent comparison | 53 |
| legacy identity | 47 |
| legacy rank | 2 |
| negative control | 19 |
| numerical corroboration | 51 |
| rank witness | 2 |
| software regression | 3 |
| unique claim | 70 |

| Computed symbolic quantity | Value |
|---|---:|
| Quark catalogue entries | 32 |
| Quark tensor-map rank | 32 |
| Gluon catalogue entries | 32 |
| Gluon tensor-map rank | 32 |

## Reference Born point

| Input | Value |
|---|---:|
| `Q2` | 4 |
| `helicity` | 1 |
| `lepton_energy` | 10 |
| `mass` | 1.5 |
| `phi` | 0.4 |
| `sqrt_s` | 5 |
| `theta` | 0.8 |

| Observable | Value |
|---|---:|
| `b_C` | 110.166177 |
| `b_G` | -123.644617 |
| `b_S` | 5.97692567 |
| `b_U` | 1684.01119 |
| B eigenvalue 1 | 1518.29969 |
| B eigenvalue 2 | 1849.72269 |
| gluon Ward residual | 9.051e-17 |
| photon Ward residual | 1.113e-16 |
| Hermiticity residual | 6.785e-17 |

All 5 hard-response diagnostics passed at the reference point and all 36 specified grid cases.

The grid is a set of numerical cases, not additional independent scientific claims. These coupling-stripped hard responses are not lithium-7 cross-section predictions.

## Independent foundations

These exact checks cover spin and target-state algebra, transverse STF tensors, the lower-spin gluon dictionary, quark/gluon covariants, and coefficient recovery. They do not cover the remaining process-response suites.

| Certificate | Verified value |
|---|---:|
| Quark exact minor rank bound | 32 |
| Quark invariant Hermitian dimension | 32 |
| Quark recovered coefficients | 32 |
| Gluon exact minor rank bound | 32 |
| Gluon invariant Hermitian dimension | 32 |
| Gluon recovered coefficients | 32 |

The seven-direction target-response determinant is `-1400/35937`. This is target-response tomography, not fourteen-TMD separation.

## Quark process evidence

32 SIDIS and 14 octupole DY row expressions were compared exactly with separately curated fixtures. 46 row kernels were evaluated using Cartesian quadrature, harmonic quadrature, and exact Gaussian moments.

64 quark/gluon coefficient sign records check finite PT algebra conditional on the stated field/link transformation. They do not establish QCD factorization or an evolution kernel.

The full profile still requires gluon response separation, an independent heavy-pair Born program, and positivity/collinear/Fourier suites.

### SIDIS reviewed response rows

| Row ID | Target `(K,m)` | Projection | Orbital rank | Weight | Phase `(recoil,target,lepton)` | Sign | Max integral residual |
|---|---|---|---:|---|---|---:|---:|
| `sidis.row.f.00.0` | `(0,0)` | `f` | 0 | `w0` | `(0, 0, 0)` | +1 | 3.33e-16 |
| `sidis.row.h.00.1` | `(0,0)` | `h` | 1 | `b0` | `(2, 0, 0)` | -1 | 2.78e-17 |
| `sidis.row.g.10.0` | `(1,0)` | `g` | 0 | `w0` | `(0, 0, 0)` | +1 | 3.33e-16 |
| `sidis.row.h.10.1` | `(1,0)` | `h` | 1 | `b0` | `(2, 0, 0)` | +1 | 5.55e-17 |
| `sidis.row.f.11.1` | `(1,1)` | `f` | 1 | `w1` | `(1, -1, 0)` | +1 | 6.94e-18 |
| `sidis.row.g.11.1` | `(1,1)` | `g` | 1 | `w1` | `(1, -1, 0)` | +1 | 5.55e-17 |
| `sidis.row.h.11.0` | `(1,1)` | `h` | 0 | `l1` | `(1, 1, 0)` | +1 | 2.22e-16 |
| `sidis.row.h.11.2` | `(1,1)` | `h` | 2 | `r1` | `(3, -1, 0)` | +1 | 6.94e-18 |
| `sidis.row.f.20.0` | `(2,0)` | `f` | 0 | `w0` | `(0, 0, 0)` | +1 | 3.33e-16 |
| `sidis.row.h.20.1` | `(2,0)` | `h` | 1 | `b0` | `(2, 0, 0)` | -1 | 2.78e-17 |
| `sidis.row.f.21.1` | `(2,1)` | `f` | 1 | `w1` | `(1, -1, 0)` | +1 | 5.55e-17 |
| `sidis.row.g.21.1` | `(2,1)` | `g` | 1 | `w1` | `(1, -1, 0)` | +1 | 6.94e-18 |
| `sidis.row.h.21.0` | `(2,1)` | `h` | 0 | `l1` | `(1, 1, 0)` | -1 | 2.22e-16 |
| `sidis.row.h.21.2` | `(2,1)` | `h` | 2 | `r1` | `(3, -1, 0)` | -1 | 2.49e-18 |
| `sidis.row.f.22.2` | `(2,2)` | `f` | 2 | `w2` | `(2, -2, 0)` | +1 | 2.08e-17 |
| `sidis.row.g.22.2` | `(2,2)` | `g` | 2 | `w2` | `(2, -2, 0)` | +1 | 6.94e-18 |
| `sidis.row.h.22.1` | `(2,2)` | `h` | 1 | `l2` | `(0, 2, 0)` | -1 | 2.26e-17 |
| `sidis.row.h.22.3` | `(2,2)` | `h` | 3 | `r2` | `(4, -2, 0)` | -1 | 4.34e-19 |
| `sidis.row.g.30.0` | `(3,0)` | `g` | 0 | `w0` | `(0, 0, 0)` | +1 | 3.33e-16 |
| `sidis.row.h.30.1` | `(3,0)` | `h` | 1 | `b0` | `(2, 0, 0)` | +1 | 5.55e-17 |
| `sidis.row.f.31.1` | `(3,1)` | `f` | 1 | `w1` | `(1, -1, 0)` | +1 | 6.94e-18 |
| `sidis.row.g.31.1` | `(3,1)` | `g` | 1 | `w1` | `(1, -1, 0)` | +1 | 5.55e-17 |
| `sidis.row.h.31.0` | `(3,1)` | `h` | 0 | `l1` | `(1, 1, 0)` | +1 | 2.22e-16 |
| `sidis.row.h.31.2` | `(3,1)` | `h` | 2 | `r1` | `(3, -1, 0)` | +1 | 6.94e-18 |
| `sidis.row.f.32.2` | `(3,2)` | `f` | 2 | `w2` | `(2, -2, 0)` | +1 | 6.94e-18 |
| `sidis.row.g.32.2` | `(3,2)` | `g` | 2 | `w2` | `(2, -2, 0)` | +1 | 2.08e-17 |
| `sidis.row.h.32.1` | `(3,2)` | `h` | 1 | `l2` | `(0, 2, 0)` | +1 | 2.08e-17 |
| `sidis.row.h.32.3` | `(3,2)` | `h` | 3 | `r2` | `(4, -2, 0)` | +1 | 7.81e-18 |
| `sidis.row.f.33.3` | `(3,3)` | `f` | 3 | `w3` | `(3, -3, 0)` | +1 | 1.73e-18 |
| `sidis.row.g.33.3` | `(3,3)` | `g` | 3 | `w3` | `(3, -3, 0)` | +1 | 1.73e-18 |
| `sidis.row.h.33.2` | `(3,3)` | `h` | 2 | `l3` | `(-1, 3, 0)` | +1 | 1.73e-17 |
| `sidis.row.h.33.4` | `(3,3)` | `h` | 4 | `r3` | `(5, -3, 0)` | +1 | 2.6e-18 |

### DY reviewed response rows

| Row ID | Target `(K,m)` | Projection | Orbital rank | Weight | Phase `(recoil,target,lepton)` | Sign | Max integral residual |
|---|---|---|---:|---|---|---:|---:|
| `dy.row.f.31.1` | `(3,1)` | `f` | 1 | `w1` | `(1, -1, 0)` | +1 | 4.16e-17 |
| `dy.row.f.32.2` | `(3,2)` | `f` | 2 | `w2` | `(2, -2, 0)` | +1 | 2.08e-17 |
| `dy.row.f.33.3` | `(3,3)` | `f` | 3 | `w3` | `(3, -3, 0)` | +1 | 8.67e-18 |
| `dy.row.g.30.0` | `(3,0)` | `g` | 0 | `w0` | `(0, 0, 0)` | -1 | 2.22e-16 |
| `dy.row.g.31.1` | `(3,1)` | `g` | 1 | `w1` | `(1, -1, 0)` | -1 | 5.55e-17 |
| `dy.row.g.32.2` | `(3,2)` | `g` | 2 | `w2` | `(2, -2, 0)` | -1 | 1.39e-17 |
| `dy.row.g.33.3` | `(3,3)` | `g` | 3 | `w3` | `(3, -3, 0)` | -1 | 6.94e-18 |
| `dy.row.h.30.1` | `(3,0)` | `h` | 1 | `b0` | `(0, 0, 2)` | +1 | 2.78e-17 |
| `dy.row.h.31.0` | `(3,1)` | `h` | 0 | `l1` | `(1, -1, 2)` | +1 | 2.22e-16 |
| `dy.row.h.31.2` | `(3,1)` | `h` | 2 | `r1` | `(-1, 1, 2)` | +1 | 4.34e-18 |
| `dy.row.h.32.1` | `(3,2)` | `h` | 1 | `l2` | `(2, -2, 2)` | +1 | 2.78e-17 |
| `dy.row.h.32.3` | `(3,2)` | `h` | 3 | `r2` | `(-2, 2, 2)` | +1 | 3.25e-19 |
| `dy.row.h.33.2` | `(3,3)` | `h` | 2 | `l3` | `(3, -3, 2)` | +1 | 5.55e-17 |
| `dy.row.h.33.4` | `(3,3)` | `h` | 4 | `r3` | `(-3, 3, 2)` | +1 | 3.25e-19 |
