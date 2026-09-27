"""Regenerate the specified Born hard-response point."""
from pathlib import Path
import subprocess
import sys
ROOT = Path(__file__).resolve().parents[1]
raise SystemExit(subprocess.call([sys.executable, str(ROOT/'src/gluon_born_response.py'),
                                  '--output', str(ROOT/'results/gluon_born_report.json')]))
