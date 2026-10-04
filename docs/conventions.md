# Conventions

- **Metrics:** Transverse Cartesian tensor algebra uses positive Euclidean δᵢⱼ. Lorentz transverse components satisfy gTⁱʲ = −δⁱʲ. The Born code uses the Minkowski metric diag(+1,−1,−1,−1).
- **Epsilon:** The transverse dualization matrix is E = [[0,1],[-1,0]], so εTˣʸ = +1. For a transverse vector represented by `z = vx + i vy`, E maps `(vx,vy)` to `(vy,−vx)`, hence `z → −i z`. This is the code's dual convention. The Born lepton tensor uses ε⁰¹²³ = +1 and therefore ε₀₁₂₃ = −1.
- **Complex STF:** A rank-n STF tensor is stored using two real numbers `(u,v)` as `u+i v`. The momentum rank-n tensor uses `(kx+i ky)^n / 2^(n−1)` for n>0. Its conjugate is the opposite helicity harmonic.
- **Target axis:** The natural z axis of the target correlator is the target-momentum/light-front axis. Moving to a virtual-photon axis requires rotating target-polarization tensors; relabeling azimuths alone is insufficient.
- **Mass and fraction:** Target-level transverse-rank normalization uses the whole target/nuclear mass MA where applicable, and target-level x is the partonic plus-momentum fraction relative to the whole nucleus, `x = k⁺/PA⁺`. The legacy symbolic script uses dimensionless polynomial momentum variables. `transverse_foundations.py` verifies the symbolic relation to tensors normalized by positive M_A and the rank-dependent mass-conversion rule. The Born code's `mass` is the heavy-quark mass mQ, not MA.
- **Rank-indexed labels:** A label such as h⁽ⁿ⁾Km records target multipole rank K, transverse target rank m, and orbital momentum rank n. The generated catalogue stores these as JSON fields. This disambiguates higher-spin structures better than familiar lower-spin names.
- **Gluon indices:** The correlator convention places transverse index i on the conjugate-amplitude field at zero and j on the field at ξ. The hard matrix order is Bᵢⱼ with `B_ij = M_i M_j*` after summing unobserved states. The response convention is `Tr(B D_g)` with `D_g = Γ/x`. The gluon spin average is in `D_g = I/2` for an unpolarized incoming gluon, not another factor in B.

## Source-to-code convention map

The current submission-review source is *Octupole parton structure of spin-3/2 targets* by Brandon B. Le and Dustin Keller. Equation labels below identify definitions; the normal test run does not read the manuscript.

| Source definition | Implemented mapping and current limit |
|---|---|
| `eq:J`, `eq:Q`, `eq:QO`, `eq:rho`, `eq:trace_metrics` | `spin32_symbolic.run_checks` uses descending helicity order `(3/2,1/2,-1/2,-3/2)`, the displayed Cartesian operators, and equivalent inverse coefficients `1/5`, `1/6`, `2/9` after distributing the factor `1/4`. The legacy run tests density reconstruction, but it does not independently derive the trace metrics or physical-state bounds. |
| `eq:euclidean`, `eq:Kn`, `eq:Kcomplex`, `eq:dual` | The code uses Euclidean transverse components, `E=[[0,1],[-1,0]]`, and the complex phase `-i` under dualization. Its `ktensor(n,x,y)` is a *dimensionful numerator*: divide by the positive target mass `M_A**n` to obtain the source's dimensionless `mathsf k_n`. The current checks omit a symbolic mass. |
| `eq:referenceSTF`, `eq:dyadic_not_STF` | The source's two-index momentum is STF, `k_i k_j - delta_ij |k|²/2`, and differs from a dyadic. Existing code checks three projected null tensors; the explicit STF-versus-dyadic regression and complete lower-spin conversion remain pending. |
| `eq:componentdefinitions`, `eq:rotatedrho`, `eq:gkernelrule` | Target rank `K` and orbital rank are distinct. The generated catalogue records `K`, `m`, and orbital rank; the independent target-state module additionally checks preparation strength `beta_3` and signed rotated components. Full gluon response kernels remain later work. |
| `eq:pairkin`, `eq:chi`, `eq:hardcontraction`, `eq:Bhard` | The Born code uses the heavy-quark mass for `mass` and returns a coupling-stripped `B` with `b_G=Im(B_xy)` and conjugate-amplitude-first `Tr(B D_g)`. Its input is one specified partonic point, not a calculation of `x_g` or a nuclear TMD. |
| `eq:SIDISdelta`, `eq:DYconv` | The existing SIDIS code checks selected Dirac traces only. Its kinematics do not implement the complete fragmentation convolution. DY, with a **sum** of incoming transverse momenta and separate antiquark labels, is absent. |
| `eq:PTspin`, `eq:octSIDISDY` | Link reversal and flavor/antiquark bookkeeping are not implemented. Any universality or factorization claim beyond finite algebra remains an analytic QCD assumption. |

Display labels such as `f`, `g`, and `h` in the legacy catalogue identify channels only. They do not specify quark flavor, antiquark status, detected hadron, Wilson-line direction, or gluon color class. The foundation catalogue now carries separate semantic fields with explicitly unspecified values where a flavor, hadron, link, or color class has not been selected.
