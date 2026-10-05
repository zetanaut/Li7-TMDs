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

The baseline manifest records 49 legacy exact outcomes, five numerical diagnostics at the reference Born point, completeness of the 36-case grid, and three software tests that execute the calculations. The grid cases are parameter evaluations, not 36 new independent scientific identities. `validate_run` requires a complete run, the full required ID set, one run identity, matching source and input digests, computed ranks, and finite results. It rejects a failed or interrupted attempt even if an older passing run remains on disk.

```bash
python scripts/run_validation.py --profile full
```

The cumulative foundations profile and standalone exact certificates can be run with:

```sh
python scripts/run_validation.py --profile foundations --output-root validation_runs
python scripts/check_foundation_certificates.py
```

Pass the printed foundations run directory to `scripts/update_report_summary.py --run-dir RUN --output /tmp/li7-foundations-summary.md`. A passing foundations profile certifies only its listed finite checks.

The cumulative quark-process profile adds 32 SIDIS and 14 octupole DY exact response rows, an independent Gaussian integral for every row, massless QED current checks, and finite link-reversal algebra conditional on the field rule:

```sh
python scripts/run_validation.py --profile quark-processes --output-root /tmp/li7-runs
python scripts/update_report_summary.py --run-dir /tmp/li7-runs/PRINTED-RUN-DIRECTORY --output /tmp/li7-summary.md
```

The cumulative gluon profile adds the 32 gluon response rows, fourteen physical octupole modes, the angular certificate, independent direct Born amplitudes, high-precision cases, and the actual 483-case dense scan:

```sh
python scripts/run_validation.py --profile gluon-processes --output-root /tmp/li7-runs
python scripts/check_gluon_angular_certificate.py
```

The limit/moment/positivity profile contains all earlier required checks and the new finite calculations. The full profile additionally requires the physical-helicity and matrix-index convention-closure checks and must finish `COMPLETE/PASS` before a full-success summary is generated. Analytical assumptions remain listed separately.

```sh
python scripts/run_validation.py --profile limits-positivity --output-root /tmp/li7-runs
python scripts/check_collinear_block_certificate.py
python scripts/run_validation.py --profile full --output-root /tmp/li7-runs
python scripts/update_report_summary.py --run-dir /tmp/li7-runs/PRINTED-FULL-RUN --output docs/generated-results.md
```

Use the printed run directory for `make_plots.py`; it reads the validated dense-scan payload. A deliberately incomplete or failed run is rejected by the summary generator. The external manuscript and network are not runtime inputs.

## Legacy commands

The original entry points and output options remain available:

```bash
python src/validate_spin32_tmd.py --output results/validation_report.json
python src/gluon_born_response.py --output results/gluon_born_report.json
python examples/reproduce_gluon_grid.py --output results/gluon_born_grid_report.json
```

These commands retain the historical JSON structure with corrected physical-helicity semantics. Use `--legacy-source-label` on the Born CLI only to reproduce the pre-correction indexing. They do not create a run manifest, so a saved JSON file by itself is not current validation evidence. The example wrappers and plotter retain their paths and default arguments; old unlabeled data are rejected for physical-helicity plots. Numerical last digits can vary by platform; compare physical values with tolerances rather than matching the last digit of Ward residuals.
