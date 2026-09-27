# Spin-3/2 polarization

The validator constructs the 4 × 4 spin generators Jx, Jy, Jz in the Jz eigenbasis (3/2, 1/2, −1/2, −3/2). It constructs the irreducible quadrupole and octupole operators

\[
Q_{ij}=\tfrac12\{J_i,J_j\}-\tfrac54\delta_{ij}I,
\qquad
O_{ijk}=J_{(i}J_jJ_{k)}-\tfrac{41}{60}(\delta_{ij}J_k+\delta_{ik}J_j+\delta_{jk}J_i).
\]

The symmetrized cubic product averages all six index permutations. The subtraction makes the octupole traceless. The code checks the traces of Q and O exactly, including all three octupole trace directions. Because the generators are Hermitian, the symmetrized operators are Hermitian; the exact reconstruction on all 16 Hermitian basis matrices also tests the operator expansion.

For any 4 × 4 Hermitian matrix ρ, the script checks the inverse expansion

\[
\rho=\frac{\operatorname{tr}\rho}{4}I
+\frac15\sum_i\operatorname{tr}(\rho J_i)J_i
+\frac16\sum_{ij}\operatorname{tr}(\rho Q_{ij})Q_{ij}
+\frac29\sum_{ijk}\operatorname{tr}(\rho O_{ijk})O_{ijk}.
\]

Testing this identity on 16 real Hermitian basis matrices establishes completeness and correct normalization of the rank-0, 1, 2, and 3 operator sectors. The program checks Ozzz = diag(3/10, −9/10, 9/10, −3/10), a helicity-one octupole matrix, and a helicity-two quadrupole matrix. These check the longitudinal spectrum and off-diagonal spin coherences, respectively. The final three symbolic checks also verify ideal resolved-transition expressions for vector, quadrupole, and octupole longitudinal moments.
