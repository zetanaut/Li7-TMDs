# Reproducing results

Use Python 3.10 or newer. Install the listed dependencies in a separate environment. A fresh run needs no manuscript, bibliography, or network access after installation.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python scripts/run_validation.py --profile baseline
```

`requirements-repro-py311.txt` records the exact package versions used for the supported Python 3.11.15 macOS arm64 checkpoint run. Install that file instead of `requirements.txt` when reproducing this environment; the recorded versions have been executed and passed `pip check` locally.

The runner prints a unique `Run directory: ...` path. Use that exact path for subsequent commands:

```bash
python scripts/update_report_summary.py --run-dir /path/printed/by/runner
python scripts/make_plots.py --run-dir /path/printed/by/runner
python -m unittest discover -s tests -v
mkdocs build --strict
python scripts/check_public_scope.py --site site
```

The baseline manifest records 49 legacy exact outcomes, five numerical diagnostics at the reference Born point, completeness of the 36-case grid, and three software tests that execute the calculations. The grid cases are parameter evaluations, not 36 new independent scientific identities. `validate_run` requires a complete run, the full required ID set, one run identity, matching source and input digests, computed ranks, and finite results. It rejects a failed or interrupted attempt even if an older passing run remains on disk. The generated summary says explicitly that full manuscript coverage is incomplete.

```bash
python scripts/run_validation.py --profile full
```

The full profile currently exits nonzero and lists its unimplemented scientific suites. Do not interpret a passing baseline as a full validation result.

## Legacy commands

The original entry points and output options remain available:

```bash
python src/validate_spin32_tmd.py --output results/validation_report.json
python src/gluon_born_response.py --output results/gluon_born_report.json
python examples/reproduce_gluon_grid.py --output results/gluon_born_grid_report.json
```

These commands reproduce the historical JSON format. They do not create a run manifest, so a saved JSON file by itself is not current validation evidence. The example wrappers and plotter also retain their original paths and default arguments. Numerical last digits can vary by platform; compare physical values with tolerances rather than matching the last digit of Ward residuals.
