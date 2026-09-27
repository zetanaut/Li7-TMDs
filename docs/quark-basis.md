# Quark basis

At leading twist the quark projections are the unpolarized channel Φ[γ⁺], helicity channel Φ[γ⁺γ₅], and transverse-spin channel Φ[iσⁱ⁺γ₅]. Their transverse operator ranks are r = 0, 0, and 1. The script labels their catalogue channels `f`, `g`, and `h`.

For each target multipole K, the code loops over target transverse ranks m = 0,…,K. Parity selects one scalar channel at m = 0 and both scalar channels for m > 0. The transverse-spin channel has one momentum branch for m = 0 and two branches for m > 0: n = |m−1| and m+1. Thus each K contributes 2 + 4K entries, giving the checked count below.

| Target multipole rank K | Independent coefficients |
|---|---:|
| 0 | 2 |
| 1 | 6 |
| 2 | 10 |
| 3 | 14 |
| **Total** | **32** |

As examples, the unpolarized K=0 target has an `f` scalar and an `h` rank-one momentum structure. A K=3, m=3 octupole contributes scalar harmonics of rank three and transverse-spin branches of orbital ranks two and four. The labels in the JSON record channel, K, m, and n; they identify the implemented covariants without relying on potentially ambiguous conventional TMD names.

The complete catalogue is also converted into a 64 × 32 linear map over the 16 target-density coordinates and four quark-channel output coordinates. Its exact generic rank is 32. Merely listing 32 labels would not establish independence.
