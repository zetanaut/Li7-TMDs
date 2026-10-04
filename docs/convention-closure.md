# Ordered currents and joint-matrix review

This page records a computational diagnosis against the approved Le–Keller source. The source manuscript and its reference values are unchanged. The historical 579-leaf checkpoint-5 full run remains a record of execution under its implemented conventions. A new `full` run adds convention diagnostics and keeps publication blocked pending an author decision.

## Physical electron and lepton indices

The direct massless spinor `u_h(p)` is an eigenvector of `diag(sigma·p_hat,sigma·p_hat)` with eigenvalue `h=±1` (twice physical helicity), obeys `gamma5 u_h=h u_h`, has norm `2p0`, and satisfies `u_h ubar_h=(1+h gamma5) slash(p)/2`. These identities are checked after arbitrary azimuthal rotations. The conventions are `g=(+---)`, `epsilon^0123=+1`, `epsilon_0123=-1`, `gamma5=i gamma0 gamma1 gamma2 gamma3`, and `Tr(gamma0 gamma1 gamma2 gamma3 gamma5)=-4i`.

The separate annihilation helper's `electron_u` branch has spinor helicity and `gamma5` eigenvalue equal to its input label. Its `antiparticle_v` branch has spinor eigenvalues opposite its input label; the barred projectors are checked independently. Interpreting the antiparticle *state* helicity requires its creation-operator convention. This branch is not used to infer the incoming-electron sign in heavy-pair production.

For `j_mu=ubar(l') gamma_mu u_h(l)`, an explicit final-spin sum and an independent Dirac trace give

```
J_mu,nu(h) = sum_h' j_mu j_nu* =
  2 [l_mu l'_nu + l_nu l'_mu - g_mu,nu (l·l')
     - i h epsilon_mu,nu,alpha,beta l^alpha l'^beta].
K_mu,nu(h) = sum_h' j_mu* j_nu = J_nu,mu(h).
```

The approved source's `eq:leptonic` has the opposite antisymmetric sign. Consequently `J(h)=L_source(h)^T=L_source(-h)`. At `l=(5,0,0,5)`, `l'=(5,3,0,4)`, physical `h=+1`, one finds `J_xy=-10i` and literal `L_source,xy(+1)=+10i` GeV². The spinor label is physically correct; the opposite source-label adapter is an index-order conversion in this comparison, not evidence of literal agreement.

Writing the full heavy-pair amplitude as `M_i=sum_mu j_mu Q_i^mu`, its square has `H_ij^mu,nu=sum Q_i^mu Q_j^nu*` and `B_ij=sum_mu,nu J_mu,nu H_ij^mu,nu`. The approved `eq:BornB` uses this amplitude-first hard trace but contracts it with `L_source,mu,nu`; the candidate contraction for physical `h` is `L_source,nu,mu H_ij^mu,nu`. This candidate is evaluated separately. For the reference point at physical `h=+1`, the literal source tuple `(b_U,b_G,b_C,b_S)` is approximately `(1684.011188,-123.644617,110.166177,5.976926)` GeV², while the amplitude-first candidate has `b_G=+123.644617` with the other three unchanged. The literal full-matrix difference is about `247.289234` GeV²; the candidate agrees with direct amplitudes at roundoff scale. Circular and elliptic polarization rates expose the difference.

The SIDIS trace is ordered differently: `Tr(Phi gamma^i Delta gamma^j)=A^j A^{i*}` for explicit spinor projectors, while the heavy Born trace is `A^i A^{j*}`. Thus the source lepton tensor can contract the SIDIS trace consistently at physical `h` without a global helicity flip. The SIDIS response, DY symmetric final-lepton tensor, beam-axis choices, target gluon circular convention, and link reversal require separate definitions. The present diagnosis changes none of them.

The proposed author action is to resolve the ordered `eq:BornB` contraction and its prose. If the amplitude-first hard trace is retained, interchange the two indices on `L_source` in that equation; keep `eq:leptonic`, `eq:Lnorm`, and `eq:SIDISmaster` under their separately verified SIDIS ordering. Recompute and relabel the physical-helicity Born `b_G/b_U` plot, reference and scan labels, beam-odd octupole moments, and longitudinal double differences only after approval. The unpolarized and linear-polarization coefficients, hard eigenvalues, angular rank, coefficient count, and operator link classifications are invariant under this antisymmetric lepton-index interchange. No manuscript correction or fixture replacement has been made here.

The existing `gluon_reconstruction.analyzer()` takes `compare_source_label()` and then uses its `b_G` with a source-labelled beam sign. Its fourteen-column linear reconstruction, angular Gram rank, and known-acceptance conditioning are unaffected by a consistent beam-label conversion, but the interpreted sign of beam-odd synthetic rates and the TTT beam-odd observable changes. Existing `b_G` figures and high-precision source-label comparisons remain historical source-labelled outputs. No candidate figure replaces them in the public site.

## Source and auxiliary parton indices

The physical source basis is ordered `(parton i, target Lambda)` with target order `3/2,1/2,-1/2,-3/2`. Direct construction from the source's `eq:qD`, `eq:gD`, and `eq:jointM` gives, for every one of the 32 covariant columns of each species at exact nonzero momentum,

```
M_source(i Lambda,j Lambda') =
  M_aux(j Lambda,i Lambda') / 2.
```

The operation is parton partial transpose and a factor of two. It is an invertible permutation and scaling of the 64 real-linear matrix coordinates, so the coefficient-map rank 32 and its inverse recovery survive the corresponding row change. A selected minor determinant depends on that row choice and normalization. Partial transpose is not a unitary similarity and does not preserve the matrix rank, eigenvalues, or PSD of an individual joint matrix. The exact Bell control has Gram rank one and partial-transpose eigenvalue `-1/2`; a nontrivial parity-averaged source-positive finite-`k_T` Gram also becomes indefinite after partial transpose.

In a positive spectral/input prescription, let `A_(X,i,Lambda)` be removal amplitudes. The source-index matrix is `M_(i Lambda,j Lambda')=sum_X A*_(X,i,Lambda) A_(X,j,Lambda')`. For a pure target `c`, `D_ij(c)=sum_Lambda,Lambda' c*_Lambda M_(i Lambda,j Lambda') c_Lambda'`. For mixed `rho`, the factor is `rho_(Lambda',Lambda)`. The code checks all 64 coefficient columns against the operator substitution with imaginary target coherence and checks a complex amplitude rate. These physical checks construct `M_source` directly; they do not treat the auxiliary matrix as PSD.

The auxiliary map is consumed by the foundation parity, rotation, coefficient-rank and expectation-coordinate checks, and by diagnostic comparisons. Those are linear coordinates. Fixed-target spectra, exact collinear blocks, finite-`k_T` spectral Grams, and source response rates use source-index matrices. The auxiliary certificate is retained as a coordinate-rank witness and is not a physical positivity certificate.

| Consumer | Matrix use | Physical-rate or PSD claim |
|---|---|---|
| `correlator_foundations.joint_covariant`, `parity_bound`, `rotation_covariance`, `joint_expectation_check` | Auxiliary linear coordinates, symmetry and coefficient checks | None; its positive target density only probes complex coordinate expectations. |
| `foundation_certificates`, analytical projectors and independent inversion | 64-by-32 coefficient map and selected minor | Rank and recovery only; selected determinant depends on row normalization. |
| `source_joint_positivity.source_mapping_check`, `conditional_positivity.convention_reproducer` | Exact source/auxiliary comparison | No PSD transfer. |
| `source_joint_positivity.spectral_recovery`, `conditional_positivity` source matrices and block checker | Direct source-index operator and spectral Gram | Conditional PSD, fixed-state spectra, collinear bounds. |
| `gluon_reconstruction`, `gluon_stokes`, response fixtures and event rates | Gluon 2-by-2 source correlator and Born analyzer | No call to the auxiliary 8-by-8 map; Born beam-label review remains separate. |

## Evidence status

A convention diagnostic may pass because it correctly detects a mismatch. The `full` manifest separately records `LITERAL_PHYSICAL_HELICITY_MISMATCH`, `AUTHOR_CORRECTION_REQUIRED`, and `BLOCKED_AUTHOR_REVIEW`. The Pages publication guard fails closed until a later author-approved same-revision integration resolves the physical assertion. Runtime validation needs neither the external manuscript nor network access.
