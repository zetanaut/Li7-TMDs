# Physics overview

```text
spin-3/2 density matrix
        ↓
irreducible multipoles K = 0,1,2,3
        ↓
light-front polarization tensors
        ↓
quark/gluon operator correlators
        ↓
2D STF covariants
        ↓
independent scalar TMD coefficients
        ↓
hard-process projection
        ↓
observable angular response
```

The symbolic code tests the spin algebra, transverse tensor construction, completeness, independence, Dirac projections, and selected hard-trace contractions. The numerical code tests a specific Born hard response. Neither code calculates the unknown nonperturbative scalar TMD functions or their lithium-7 magnitudes.

The target-spin space has dimension four, so a Hermitian density matrix has 16 real operator coordinates: one normalization coordinate and 15 polarization coordinates. These organize as 1 + 3 + 5 + 7 for rotational ranks 0 through 3. The transverse correlator count of 2 + 6 + 10 + 14 is a different count: it includes partonic channels and allowed transverse momentum harmonics.
