"""Standalone exact rank-certificate checker; no network or manuscript needed."""
from __future__ import annotations
import argparse
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from foundation_certificates import read_certificate,verify_certificate


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--certificate-dir',type=Path,default=ROOT/'certificates')
    args=parser.parse_args()
    for species in ('quark','gluon'):
        path=args.certificate_dir/f'{species}_rank.json'
        verdict=verify_certificate(read_certificate(path))
        print(f"{species}: rank >= {verdict['rank_lower_bound']}, determinant {verdict['determinant']}, SHA-256 {verdict['certificate_sha256']}")
    return 0

if __name__=='__main__':raise SystemExit(main())
