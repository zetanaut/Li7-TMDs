"""Run the exact symbolic validator from any working directory."""
from pathlib import Path
import subprocess
import sys
ROOT = Path(__file__).resolve().parents[1]
raise SystemExit(subprocess.call([sys.executable, str(ROOT/'src/validate_spin32_tmd.py'),
                                  '--output', str(ROOT/'results/validation_report.json')]))
