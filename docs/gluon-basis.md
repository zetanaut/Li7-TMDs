# Gluon basis

The transverse gluon correlator has trace/unpolarized, antisymmetric/circular, and symmetric-traceless/linear polarization sectors. Their transverse channel ranks are r = 0, 0, and 2. The code labels them `f`, `g`, and `h`.

For the linear sector, the allowed momentum ranks are:

| Target transverse rank m | Independent orbital ranks n |
|---:|---|
| 0 | 2 |
| 1 | 1, 3 |
| 2 | 0, 4 |
| 3 | 1, 5 |

Only the sum and difference helicity harmonics exist in two physical transverse dimensions. In particular, m = 2 has no independent n = 2 middle branch. The exact rank-two and rank-three STF null tests on the [transverse STF page](transverse-stf.md) show why the tempting extra tensor expressions do not add coefficients.

The scalar channels follow the same K and m loop as for quarks. Each target multipole rank contributes 2 + 4K entries, hence 2, 6, 10, and 14 for K = 0, 1, 2, 3. The script checks both the 32-entry catalogue and the exact rank 32 of its 64 × 32 tensor map. This count is per specified gauge-link/color class.
