# Generated result summary

This page describes the **full profile**. The declared executable program passed; analytical QCD inputs and excluded dynamic calculations remain separate.

The summary was generated from a complete run whose required IDs, source and input digests, computed ranks, and result identities were checked together.

Source revision: `3d250407e061fb14a0f36b73f341ae68c4fffa92`; dirty source: `false`; scientific source digest: `633e28117ac78c47c95067bab3a4e17fa92f03c1f4fa6d6045f64769695cfde7`.

| Executed evidence | Count |
|---|---:|
| conditional algebra | 85 |
| exact fixture | 92 |
| exact identity | 274 |
| exact rank | 6 |
| numerical diagnostic | 56 |
| numerical parameter cases | 28 |
| software test | 38 |

| Evidence role | Records |
|---|---:|
| component case | 216 |
| independent bound | 2 |
| independent comparison | 142 |
| legacy identity | 49 |
| legacy rank | 2 |
| negative control | 35 |
| numerical corroboration | 51 |
| rank witness | 2 |
| software regression | 3 |
| unique claim | 77 |

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

These exact checks cover spin and target-state algebra, transverse STF tensors, the lower-spin gluon dictionary, quark/gluon covariants, and coefficient recovery. Their scope is exact finite-dimensional algebra.

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

## Gluon responses and independent Born checks

32 gluon rows were compared through Cartesian-trace and complex-helicity routes. 14 octupole rows were obtained from four physical spin rates each.

The exact angular Gram determinant is `1/8192`. The independent-Born synthetic response has rank 14 and maximum coefficient recovery error 4.66e-14.

The 36-case grid compared 36 direct amplitudes; the dense scan evaluated and independently compared 483 cases. Worst matrix component residual: 2.73e-12 GeV².

4 independently calculated cases used 50 and 80 decimal digits from string inputs. The physical spinor helicity is the negative of the source lambda label in the matrix comparison; this sign-label discrepancy awaits author review.

| Gluon response ID | Target rank | Channel | Orbital rank |
|---|---:|---|---:|
| `gluon.response.00.f.0` | 0 | `f` | 0 |
| `gluon.response.00.h.2` | 0 | `h` | 2 |
| `gluon.response.10.g.0` | 1 | `g` | 0 |
| `gluon.response.10.h.2` | 1 | `h` | 2 |
| `gluon.response.11.f.1` | 1 | `f` | 1 |
| `gluon.response.11.g.1` | 1 | `g` | 1 |
| `gluon.response.11.h.1` | 1 | `h` | 1 |
| `gluon.response.11.h.3` | 1 | `h` | 3 |
| `gluon.response.20.f.0` | 2 | `f` | 0 |
| `gluon.response.20.h.2` | 2 | `h` | 2 |
| `gluon.response.21.f.1` | 2 | `f` | 1 |
| `gluon.response.21.g.1` | 2 | `g` | 1 |
| `gluon.response.21.h.1` | 2 | `h` | 1 |
| `gluon.response.21.h.3` | 2 | `h` | 3 |
| `gluon.response.22.f.2` | 2 | `f` | 2 |
| `gluon.response.22.g.2` | 2 | `g` | 2 |
| `gluon.response.22.h.0` | 2 | `h` | 0 |
| `gluon.response.22.h.4` | 2 | `h` | 4 |
| `gluon.response.30.g.0` | 3 | `g` | 0 |
| `gluon.response.30.h.2` | 3 | `h` | 2 |
| `gluon.response.31.f.1` | 3 | `f` | 1 |
| `gluon.response.31.g.1` | 3 | `g` | 1 |
| `gluon.response.31.h.1` | 3 | `h` | 1 |
| `gluon.response.31.h.3` | 3 | `h` | 3 |
| `gluon.response.32.f.2` | 3 | `f` | 2 |
| `gluon.response.32.g.2` | 3 | `g` | 2 |
| `gluon.response.32.h.0` | 3 | `h` | 0 |
| `gluon.response.32.h.4` | 3 | `h` | 4 |
| `gluon.response.33.f.3` | 3 | `f` | 3 |
| `gluon.response.33.g.3` | 3 | `g` | 3 |
| `gluon.response.33.h.1` | 3 | `h` | 1 |
| `gluon.response.33.h.5` | 3 | `h` | 5 |

| Octupole response | Signed cone factor | Orbital rank | Source factor |
|---|---:|---:|---:|
| `gluon.oct.g.30` | -9/25 | 0 | 1 |
| `gluon.oct.h.30.2` | -9/25 | 2 | 1/2 |
| `gluon.oct.f.31` | 8/25 | 1 | 1 |
| `gluon.oct.g.31` | 8/25 | 1 | 1 |
| `gluon.oct.h.31.1` | 8/25 | 1 | 1/2 |
| `gluon.oct.h.31.3` | 8/25 | 3 | 1/4 |
| `gluon.oct.f.32` | 48/25 | 2 | 1 |
| `gluon.oct.g.32` | 48/25 | 2 | 1 |
| `gluon.oct.h.32.0` | 48/25 | 0 | 1 |
| `gluon.oct.h.32.4` | 48/25 | 4 | 1/4 |
| `gluon.oct.f.33` | 32/25 | 3 | 1 |
| `gluon.oct.g.33` | 32/25 | 3 | 1 |
| `gluon.oct.h.33.1` | 32/25 | 1 | 1 |
| `gluon.oct.h.33.5` | 32/25 | 5 | 1/4 |

## Collinear, Fourier, local moments, and conditional positivity

The full Cartesian angular projection evaluated 64 covariants. Quarks: 7 angular candidates and 6 straight-link survivors. Gluons: 6 and 5.

The cutoff-tail fixture integrates to `-log(Lambda**2) + log(Lambda**2 + R**2)`; its radial limit diverges logarithmically. Angular selection is not a renormalized PDF integral.

Exact Fourier differentiation covered 6 ranks and 20 independently integrated Gaussian/quartic rank-branch cases, plus 6 Gaussian inverse Hankel cases. The worst scaled Cartesian/Bessel error was 1.1e-12 for gaussian rank 5 branch -1.

Rotational coupling checked 12 finite N/bilinear cases; 12 charge-conjugation rows retain independent antiquark terms. No numerical nuclear moment or QCD sum-rule value is inferred.

The independent complex Gram example has rank 3; the parity-averaged source-index Gram has rank 6 and recovers 32 quark coefficients. The collinear block checker verified 3 quark and 2 gluon coupled blocks. 11 exact collinear cases include interiors, boundaries, and violations.

These positivity statements assume a positive spectral/input prescription. They do not impose pointwise positivity on arbitrary subtracted TMDs.
