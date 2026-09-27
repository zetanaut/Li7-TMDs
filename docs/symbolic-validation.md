# Exact symbolic validation

`src/validate_spin32_tmd.py` is a SymPy program for exact finite-dimensional algebra. It does not fit data, choose numerical TMD models, perform Monte Carlo inference, or calculate nuclear wave functions. Run it with

```bash
python src/validate_spin32_tmd.py --output results/validation_report.json
```

A failed `check()` raises an assertion and exits nonzero. The JSON report has a `PASS` value for each successful check plus the two complete catalogues. The committed report was regenerated from the source in this repository: **all 49 exact symbolic checks pass**.

## What the 49 checks cover

- Spin-3/2 quadrupole and octupole trace removal; inversion of the density-matrix expansion on all 16 Hermitian matrix coordinates; longitudinal and helicity-changing operator matrices. The inverse-map check establishes completeness of K = 0, 1, 2, 3 spin multipoles.
- Two-dimensional momentum STF tensors through rank five, three vanishing candidate gluon tensors, low and high helicity branches for target transverse ranks 1–3, and the sine/dual relation.
- Quark and gluon catalogues of 32 entries each, target-rank counts 2, 6, 10, 14, and exact rank 32 for each 64 × 32 tensor map at nonzero transverse momentum `(kx,ky)=(2,1)`.
- Duality of quark and fragmentation reconstruction matrices with their Dirac projection matrices. The Gram matrices are both the 4 × 4 identity. This matters because the scalar coefficients are defined through operator projections and must reconstruct the same correlator.
- Twelve channel-specific electromagnetic hard traces, the scalar–helicity–Collins master contraction, and three ideal spin-transition moment identities.

The leading-twist quark projections are Φ[γ⁺], Φ[γ⁺γ₅], and Φ[iσⁱ⁺γ₅]. The script creates Dirac gamma matrices through SymPy, light-front γ⁺ and γ⁻, and the reconstruction/projection matrix pairs. Its trace checks use these actual matrices.

## Generated report summary

The [generated result summary](generated-results.md) reads the current JSON reports and displays the check and catalogue counts.

The rank test is essential: 32 formal terms can still be dependent. At generic nonzero transverse momentum, the independent columns show that the proposed covariants span a 32-dimensional space. The script uses exact rational and symbolic arithmetic for these identities.
