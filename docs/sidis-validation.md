# SIDIS trace validation

The symbolic validator checks selected leading electromagnetic semi-inclusive deep-inelastic scattering (SIDIS) kernels. It forms the leading transverse hard tensor from actual gamma-matrix traces of quark reconstruction matrices, photon vertices γⁱ, and fragmentation reconstruction matrices. It checks twelve channel/output combinations against their expected 2 × 2 transverse matrices.

It then contracts a symbolic quark/fragmentation hard tensor with a symbolic leptonic transverse tensor. The tested master expression has the scalar term `F D`, helicity term `λ dep G D`, and Collins-like transverse-spin term `ε H (Tx py + Ty px)`. The test checks the azimuthal and spin signs that follow from the adopted gamma and epsilon conventions, instead of assigning those kernels by analogy with proton expressions.

This is an exact tree-level Dirac-trace identity in the specified leading-power setup. It is not an all-orders SIDIS factorization proof, nor does it calculate unknown TMD coefficient values.
