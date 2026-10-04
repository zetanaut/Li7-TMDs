# Conventions

- **Metrics:** Transverse Cartesian tensor algebra uses positive Euclidean δᵢⱼ. Lorentz transverse components satisfy gTⁱʲ = −δⁱʲ. The Born code uses the Minkowski metric diag(+1,−1,−1,−1).
- **Epsilon:** The transverse dualization matrix is E = [[0,1],[-1,0]], so εTˣʸ = +1. For a transverse vector represented by `z = vx + i vy`, E maps `(vx,vy)` to `(vy,−vx)`, hence `z → −i z`. This is the code's dual convention. The Born lepton tensor uses ε⁰¹²³ = +1 and therefore ε₀₁₂₃ = −1.
- **Complex STF:** A rank-n STF tensor is stored using two real numbers `(u,v)` as `u+i v`. The momentum rank-n tensor uses `(kx+i ky)^n / 2^(n−1)` for n>0. Its conjugate is the opposite helicity harmonic.
- **Target axis:** The natural z axis of the target correlator is the target-momentum/light-front axis. Moving to a virtual-photon axis requires rotating target-polarization tensors; relabeling azimuths alone is insufficient.
- **Mass and fraction:** Target-level transverse-rank normalization uses the whole target/nuclear mass MA where applicable, and target-level x is the partonic plus-momentum fraction relative to the whole nucleus, `x = k⁺/PA⁺`. The legacy symbolic script uses dimensionless polynomial momentum variables. `transverse_foundations.py` verifies the symbolic relation to tensors normalized by positive M_A and the rank-dependent mass-conversion rule. The Born code's `mass` is the heavy-quark mass mQ, not MA.
- **Rank-indexed labels:** A label such as h⁽ⁿ⁾Km records target multipole rank K, transverse target rank m, and orbital momentum rank n. The generated catalogue stores these as JSON fields. This disambiguates higher-spin structures better than familiar lower-spin names.
- **Gluon indices:** The correlator convention places transverse index i on the conjugate-amplitude field at zero and j on the field at ξ. The hard matrix order is Bᵢⱼ with `B_ij = M_i M_j*` after summing unobserved states. The response convention is `Tr(B D_g)` with `D_g = Γ/x`. The gluon spin average is in `D_g = I/2` for an unpolarized incoming gluon, not another factor in B.
- **Electron current indices:** The physical incoming-electron helicity eigenvalue is `h`. For `j_mu=ubar(l') gamma_mu u_h(l)`, the amplitude-first tensor is `J_mu,nu=sum j_mu j_nu*`, while the historical source lepton tensor is `L_source(h)=J(h)^T`. The source Born hard trace is amplitude-first, so its literal same-label contraction is under [author review](convention-closure.md). This does not change the separately ordered SIDIS trace.
- **Joint parton–target indices:** The source spectral matrix uses `(parton i,target Lambda)` and `M=A* A`. The older foundation matrix is an auxiliary coefficient map; `M_source=T_parton(M_aux)/2`. The row transform preserves coefficient-map rank but not individual-matrix positivity.

## Source-to-code convention map

The current submission-review source is *Octupole parton structure of spin-3/2 targets* by Brandon B. Le and Dustin Keller. Equation labels below identify definitions; the normal test run does not read the manuscript.

| Source definition | Implemented mapping and current limit |
|---|---|
| `eq:J`, `eq:Q`, `eq:QO`, `eq:rho`, `eq:trace_metrics` | The legacy route uses descending helicity order `(3/2,1/2,-1/2,-3/2)`; the independent ladder and Clebsch–Gordan routes derive trace metrics, density inversion, and physical-state bounds. |
| `eq:euclidean`, `eq:Kn`, `eq:Kcomplex`, `eq:dual` | The code uses Euclidean transverse components and `E=[[0,1],[-1,0]]`. The legacy `ktensor(n,x,y)` is a dimensionful numerator; the independent rank-zero-through-five route uses positive symbolic `M_A` and checks the normalized tensors and mass conversion. |
| `eq:referenceSTF`, `eq:dyadic_not_STF` | The source's two-index momentum is STF, `k_i k_j - delta_ij |k|²/2`, and differs from a dyadic. The explicit STF-versus-dyadic regression and lower-spin conversion now run in `foundations`. |
| `eq:componentdefinitions`, `eq:rotatedrho`, `eq:gkernelrule` | Target rank `K` and orbital rank are distinct. The catalogue, independent target-state preparations, and all 32 gluon response kernels now have finite checks. |
| `eq:pairkin`, `eq:chi`, `eq:leptonic`, `eq:BornB`, `eq:hardcontraction`, `eq:Bhard` | The Born code returns a coupling-stripped `B` with `b_G=Im(B_xy)` and `Tr(B D_g)`. Its source-labelled lepton contraction disagrees with physical same-label electron helicity; the candidate index-interchanged trace is segregated for author review. |
| `eq:SIDISdelta`, `eq:DYconv` | All 32 SIDIS and 14 octupole DY finite response rows and independent Gaussian convolution routes run. The SIDIS hadronic trace is checked with its own reversed hard-current order. |
| `eq:PTspin`, `eq:octSIDISDY` | Finite link-reversal and flavor/antiquark bookkeeping run. Wilson-line matrix elements and factorization remain analytical QCD assumptions. |

Display labels such as `f`, `g`, and `h` in the legacy catalogue identify channels only. They do not specify quark flavor, antiquark status, detected hadron, Wilson-line direction, or gluon color class. The foundation catalogue now carries separate semantic fields with explicitly unspecified values where a flavor, hadron, link, or color class has not been selected.
