# Ordered currents and joint-matrix closure

The author-approved checkpoint-7 correction replaces the Born lepton-index contraction while retaining the source's plus-i tensor and the physical spinor label. The immutable pre-correction source and the corrected manuscript copy are held outside this repository. The old literal mismatch remains an expected negative regression.

## Physical electron and lepton indices

For a positive-energy electron, `h=±1` is the eigenvalue of `diag(sigma·p_hat,sigma·p_hat)`. The independent spinor also obeys `gamma5 u_h=h u_h`. With `j_mu=ubar(l') gamma_mu u_h(l)`, the retained source tensor is `L_source[mu,nu]=sum j_mu* j_nu`, with its plus-i epsilon term for `epsilon^0123=+1`. The amplitude-first tensor is `J[mu,nu]=sum j_mu j_nu*=L_source[nu,mu]`.

The heavy trace is `H_ij[mu,nu]=sum Q_i^mu Q_j^{nu*}`. Squaring the scalar amplitude requires `B_ij=N_B sum L_source[nu,mu] H_ij[mu,nu]`. The pre-correction contraction used `L_source[mu,nu]` with the same heavy trace. For the reference event at physical `h=+1`, the corrected reduced coefficients are approximately `(1684.011188,+123.644617,110.166177,5.976926)` GeV². The old wrong-order calculation gives negative `b_G` at the same numeric label and differs in the full matrix by about `247.289234` GeV². It is available only as `evaluate_legacy_source_label` or `--legacy-source-label`.

The independent scalar-amplitude spinors and the corrected trace agree at both physical helicities. A separate SIDIS trace has reverse hard-current ordering, so its current contraction and helicity signs remain unchanged. The target gluon Gram matrix, Wilson-link relations, and TMD definitions remain unchanged. Beam-odd Born analyzer and synthetic rates are regenerated with the corrected physical label. No measured asymmetry or experimental fit is claimed.

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

The `full` profile requires current and trace agreement, the legacy mismatch negative control, and the source/auxiliary partial-transpose checks. A complete current run can establish `READY_FOR_PUBLICATION_REVIEW`; its scientific record grants no standing deployment permission. Publication requires a separate exact-commit workflow approval. Validation is offline and does not read the external manuscript.
