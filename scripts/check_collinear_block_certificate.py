#!/usr/bin/env python3
"""Check the source-index block witness from independent source amplitudes."""
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from collinear_block_certificate import check

if __name__=='__main__':
    path=Path(sys.argv[1]) if len(sys.argv)>1 else ROOT/'certificates'/'collinear_blocks.json'
    verdict=check(path)
    print(f"Collinear block certificate PASS: {verdict['quark_blocks']} quark and {verdict['gluon_blocks']} gluon coupled blocks")
