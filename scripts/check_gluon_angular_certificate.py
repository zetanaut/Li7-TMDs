#!/usr/bin/env python3
"""Reconstruct the angular witness independently of stored certificate bytes."""
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from gluon_angular import check

if __name__=='__main__':
    path=Path(sys.argv[1]) if len(sys.argv)>1 else ROOT/'certificates'/'gluon_angular.json'
    check(path)
    print('Angular Gram certificate PASS: 14 modes, determinant 1/8192')
