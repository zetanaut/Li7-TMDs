# Two-dimensional transverse STF tensors

For transverse rank n > 0, a real symmetric traceless tensor has just two independent real components. The code stores them as `u + i v`; for example rank two has the matrix `[[u,v],[v,-u]]`. A transverse rotation acts on this complex component with harmonics \(e^{+in\phi}\) and \(e^{-in\phi}\). Momentum tensors are built from `(kx + i ky)^n / 2^(n−1)` through rank five and checked to be traceless.

A target tensor with transverse rank m > 0 and a partonic channel of rank r contains only helicities ±m and ±r. Their product gives sum and difference harmonics, so the orbital ranks are

\[
n=|m-r|\quad\text{or}\quad n=m+r.
\]

For m = 0, n = r. This is **SO(2) transverse-helicity coupling**, not the full SO(3) addition rule with all intermediate angular momenta. In two dimensions an STF tensor has only two harmonics; there are no intermediate transverse irreducible states to support the extra ranks.

The linearly polarized gluon channel has r = 2. For m = 2 its branches are n = 0 and 4, not n = 0, 2, and 4. An explicit code identity makes the missing middle branch visible. With

\[
A=\begin{pmatrix}a&b\\b&-a\end{pmatrix},\quad
B=\begin{pmatrix}c&d\\d&-c\end{pmatrix},
\]

`stf2(A*B)` is identically zero: the symmetric part of the product is `(ac+bd)I`, whose trace subtraction vanishes. The code also verifies `stf2(E*A*B)=0` for E = [[0,1],[-1,0]], and a corresponding dual rank-three double-contraction identity. These three exact null tests remove the candidate intermediate gluon structures that would inflate the count to 35. A pure trace can be absorbed into a trace-sector coefficient; it is not a new independent linear-polarization tensor.
