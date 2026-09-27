# Spin-3/2 TMD computational companion

[![Validation](https://github.com/zetanaut/Li7-TMDs/actions/workflows/validation.yml/badge.svg)](https://github.com/zetanaut/Li7-TMDs/actions/workflows/validation.yml)
[![Pages](https://github.com/zetanaut/Li7-TMDs/actions/workflows/pages.yml/badge.svg)](https://github.com/zetanaut/Li7-TMDs/actions/workflows/pages.yml)

## What this repository is

Computational checks and examples for leading-twist quark and gluon transverse-momentum-dependent (TMD) structures of a polarized spin-3/2 target. Polarized lithium-7 motivates the study, while the implemented spin algebra is general to spin 3/2. The repository contains code, generated validation reports, diagnostic plots, and explanations of the code.

## Main result checked

**32 independent quark structures and 32 independent gluon structures**, per specified flavor and gauge-link/color class. Each sector has **2 + 6 + 10 + 14 = 32** coefficients at target multipole ranks K = 0, 1, 2, 3.

## Why spin 3/2 is different

A spin-3/2 density matrix has rank-0 unpolarized, rank-1 vector, rank-2 quadrupole, and rank-3 octupole parts. The octupole creates the new target-polarization structures. Two-dimensional symmetric traceless (STF) tensor counting is essential, especially for linearly polarized gluons.

## What the symbolic script checks

`src/validate_spin32_tmd.py` checks spin multipoles, density-matrix completeness and inversion, STF identities, elimination of spurious gluon structures, 32-column basis independence, quark Dirac projection identities, and selected electromagnetic SIDIS gamma-matrix trace identities. **All 49 exact symbolic checks pass** in the generated [validation report](results/validation_report.json).

## What the Born script checks

`src/gluon_born_response.py` calculates a specified, coupling-stripped Born response for γ* g → Q Q̄. It checks physical kinematics, photon and gluon Ward identities, hard-matrix Hermiticity, positivity, and its four Stokes coefficients. The [reference point](results/gluon_born_report.json) and [36-point grid](results/gluon_born_grid_report.json) pass their required checks.

## What the repository does not do

It does not calculate nonperturbative lithium-7 TMDs, fit data, contain a nuclear wave-function model, or prove arbitrary-process all-orders TMD factorization. The Born coefficients are hard-response diagnostics, not lithium-7 cross-section predictions.

## Install and run

```bash
git clone https://github.com/zetanaut/Li7-TMDs.git
cd Li7-TMDs
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python src/validate_spin32_tmd.py
python src/gluon_born_response.py
```

Full reproduction:

```bash
python examples/reproduce_symbolic_validation.py
python examples/reproduce_gluon_point.py
python examples/reproduce_gluon_grid.py
python scripts/make_plots.py
python scripts/update_report_summary.py
python -m unittest discover -s tests -v
```

See the [documentation site](https://zetanaut.github.io/Li7-TMDs/) for conventions, derivations, checks, and a recommended learning path.

Dustin Keller · University of Virginia · MIT license.
