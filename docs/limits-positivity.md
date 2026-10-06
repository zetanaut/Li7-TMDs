# Collinear limits, Fourier tensors, local moments, and positivity

The `limits-positivity` and `full` profiles execute the finite calculations described here. The source is the approved Le–Keller manuscript identified by SHA-256 `3b5aaff51a77932ad561c1137a6d1bb5f0e4c60353add1d6f9034f2d7e2b892d`. Validation runs use the compact source definitions encoded in this repository; they do not read the external manuscript.

## Angular projection and the straight link

The code integrates each Cartesian covariant over azimuth at fixed nonzero $|k_T|$ using $d\phi/(2\pi)$ and independently computes its harmonic form. It does so for all 32 quark and 32 gluon coefficients. Seven quark candidates survive angular averaging: `f_00`, `f_20`, `g_10`, `g_30`, `h_11[0]`, `h_21[0]`, and `h_31[0]`. The six gluon candidates are `f_00`, `f_20`, `g_10`, `g_30`, `h_22[0]`, and `h_32[0]`. The straight-link transformation then excludes quark `h_21[0]` and gluon `h_32[0]` from collinear PDFs. They remain legitimate finite-$k_T$ coefficients; in particular, the gluon LTT coherence remains present in the TMD basis.

For a smooth radial fixture, setting $k_T=0$ leaves the same algebraic zero-rank terms, but the inverse projector can be singular. Neither this forward evaluation nor the fixed-radius angular average is the renormalized straight-link operator. The synthetic tail $1/[\pi(k_T^2+\Lambda^2)]$ has cutoff integral $\log(1+R^2/\Lambda^2)$, which diverges as $R\to\infty$. This illustrates why angular selection does not establish radial convergence or small-$b_T$ matching.

Population operators give $U_{3/2}=f+f_Q$, $U_{1/2}=f-f_Q$, $G_{3/2}=3g/2+3g_O/10$, and $G_{1/2}=g/2-9g_O/10$. The quark raising block yields outer coherence $\sqrt3(h+2h_O/5)$ and middle coherence $2h-6h_O/5$. The inclusive Born check uses the existing transverse Dirac currents and $F_1=\frac12\sum_q e_q^2(q+\bar q)$ for each allowed target multipole. Deriving an inclusive observable from identified-hadron SIDIS additionally needs an energy-weighted hadron sum, the fragmentation momentum sum rule, and the variable-change Jacobians.

## Fourier phases and normalization

With the positive forward exponential, the Gaussian $c(k_T^2)=e^{-k_T^2/\Lambda^2}/(\pi\Lambda^2)$ transforms to $e^{-\Lambda^2b_T^2/4}$. The rank-one harmonic transforms to $i\Lambda^2(b_x+i b_y)e^{-\Lambda^2b_T^2/4}/(2M_A)$; the rank-five harmonic carries $i^5\Lambda^{10}(b_x+i b_y)^5e^{-\Lambda^2b_T^2/4}/(2M_A)^5$. The source STF complex component has an extra $2^{1-n}$ for $n>0$, while rank zero is exactly one. The transform of the conjugated momentum harmonic is evaluated separately. A real odd-rank Cartesian component has an imaginary transform and satisfies $\widetilde F(b)^*=\widetilde F(-b)$.

All independent Cartesian components through rank five are differentiated exactly. Cartesian two-dimensional quadrature and radial $J_n$ integration are separately coded for both harmonic branches, a normalized Gaussian, and a smooth quartic exponential. The inverse Gaussian Hankel calculation checks ranks zero, one, and five at two nonzero momenta with the $(2\pi)^{-2}$ inverse normalization. The calculations require integrable polynomial-weighted radial fixtures and $M_A,\Lambda>0$; they are Fourier identities, not perturbative QCD matching or a claim about an unregulated TMD tail.

## Local rotation content and antiquarks

A vector or axial bilinear restricts to spin $0\oplus1$, the antisymmetric tensor to $1\oplus1$, and each of the $N-1$ derivative indices to $0\oplus1$. Two independent coupling methods agree through $N=4$, with an induction check through $N=7$: the ambient maximum spin is $N$. Trace and mixed-symmetry projection cannot create an absent higher spin. Thus a rank-three target moment is absent for $N=1,2$; a reduced rank-three matrix element is first allowed at $N=3$.

An explicit Dirac charge-conjugation matrix gives bilinear parities $(-,+,-)$ for vector, axial, and transversity. With $N-1$ derivatives, the positive-$x$ combinations are $f^q+(-1)^Nf^{\bar q}$, $g^q+(-1)^{N-1}g^{\bar q}$, and $h^q+(-1)^Nh^{\bar q}$. The first rank-two vector and rank-three axial/transversity combinations vanish as *local reduced matrix elements*; the rank-three second moments also vanish. Ordinary-integral interpretations need small-$x$ convergence and distributional care. The calculation does not infer a gluon-helicity first-moment sum rule.

The whole-nucleus variables satisfy $x_N=x/y$, $k_{\mathrm{rel},T}=k_T-(x/y)p_{NT}$, and $0<x\le y\le1$ for the positive-$x$ term. A Mellin variable change supplies $y^{N-1}$. Proton/neutron number normalizations of 3 and 4 alone do not fix total momentum. The total tensor-polarized momentum constraint is conditional on the full QCD momentum generator and does not impose a separate zero on every flavor or on gluons.

## Three positivity questions

Target density positivity, the independently constructed Born hard-matrix Gram positivity, and parton–target correlator positivity have distinct assumptions. The last is tested for a specified positive spectral/input prescription. It is not imposed pointwise on arbitrary UV/rapidity-subtracted TMDs.

For a fixed target state, the quark and gluon $2\times2$ characteristic polynomial is $t^2-Ft+(F^2-G^2-X^2-Y^2)/4$. Positive semidefiniteness requires $F\ge\sqrt{G^2+X^2+Y^2}$, including a nonnegative trace. Exact interior, boundary, negative-trace, and indefinite examples are checked.

The joint source-index matrices are built from the defining spin operators with the manuscript's upper off-diagonal `+i` convention. Exact rational-complex spectral amplitudes give a positive Gram matrix; parity averaging preserves positivity, and the source-index 32-column map recovers its coefficients at an admissible nonzero $k_T$. The earlier foundation auxiliary map uses the opposite Pauli-$\sigma_2$ orientation. Its relation to the source matrix is an explicitly checked parton-index transpose and factor of two. A Bell-state counterexample demonstrates that partial transpose itself is **not** positivity preserving; the spectral example is constructed directly in source indices.

The exact collinear block certificate retains the common $1/2$ matrix factor. In the source conventions, the quark outer and middle conditions are $3(h+2h_O/5)^2\le(U_{3/2}+G_{3/2})(U_{1/2}-G_{1/2})$ and $|2h-6h_O/5|\le U_{1/2}+G_{1/2}$. The gluon two-unit flip gives $12h_{1TT}^{g,2}\le(U_{3/2}+G_{3/2})(U_{1/2}+G_{1/2})$, together with the diagonal bounds $U_\Lambda\ge|G_\Lambda|$. A separate checker reconstructs every block entry and determinant from the curated source amplitudes. Exact interior, saturated, zero-diagonal, and violating matrices are evaluated.

Pairwise minors are insufficient: the exact matrix with diagonal 1 and off-diagonal $-3/4$ has all pairwise principal minors $7/16$ but spectrum $(-1/2,7/4,7/4)$. A singular example with all leading principal minors zero has a negative nonleading minor. A singular positive matrix is accepted as a boundary.

## Analytical scope

The finite checks condition on the stated straight-link transformation, renormalized local-operator identification, suitable moment convergence, and a positive spectral prescription where positivity is invoked. They do not calculate QCD matching or evolution kernels, Wilson-line matrix elements, a factorization theorem, nuclear distributions, or experimental sensitivity. A passing computational profile means the declared executable claims ran; it does not settle these analytical inputs or authorize publication.

See the [generated result summary](generated-results.md) for the actual run-specific counts, errors, certificate outcomes, and provenance.
