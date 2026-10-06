"""Regenerate the specified Born hard-response point."""
from pathlib import Path
import subprocess
import sys
ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    return subprocess.call([sys.executable, str(ROOT/'src/gluon_born_response.py'),
                            '--output', str(ROOT/'results/gluon_born_report.json')])


if __name__ == '__main__':
    raise SystemExit(main())
