# Reproducing results

Use Python 3.10 or newer. From a fresh checkout:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python examples/reproduce_symbolic_validation.py
python examples/reproduce_gluon_point.py
python examples/reproduce_gluon_grid.py
python scripts/make_plots.py
python scripts/update_report_summary.py
python -m unittest discover -s tests -v
mkdocs build --strict
```

The commands regenerate the three JSON files in `results/` and six figures in both PNG and SVG under `docs/figures/`. The output paths are configurable with `--output` on the validation/response runners and `--input`/`--output-dir` on the plotter. The main scripts can also be run directly:

```bash
python src/validate_spin32_tmd.py --output results/validation_report.json
python src/gluon_born_response.py --output results/gluon_born_report.json
python examples/reproduce_gluon_grid.py --output results/gluon_born_grid_report.json
```

The symbolic validator exits nonzero on a failed identity. The Born evaluator raises an error if a required kinematic, Ward, Hermiticity, or positivity test fails; the grid runner exits nonzero if any of its 36 evaluations fails. The JSON is indented, deterministic for the same Python/NumPy/SymPy environment, and contains no random sampling. Last-digit floating-point differences may occur across platforms.

The validation report has 49 `PASS` entries and two catalogues. The reference hard-response report has five passing checks. The grid JSON records 36 points with `all_pass=true`. The plots display only calculated hard coefficients and residuals, without a phenomenological TMD model.
