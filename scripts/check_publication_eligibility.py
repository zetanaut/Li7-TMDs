#!/usr/bin/env python3
"""Fail closed until same-revision full evidence and author resolution exist."""
from __future__ import annotations
import argparse
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from validation_evidence import CONVENTION_REVIEW,validate_run,assert_convention_review  # noqa: E402

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir',type=Path)
    parser.add_argument('--require-publication-authorization',action='store_true')
    args=parser.parse_args()
    if args.run_dir is None:
        print('Publication blocked: no same-revision full convention evidence.',file=sys.stderr)
        return 2
    manifest,document=validate_run(args.run_dir)
    if manifest['profile']!='full':
        print('Publication blocked: full profile required.',file=sys.stderr)
        return 2
    review=manifest.get('convention_review')
    current=next(x for x in document['results'] if x['check_id']=='convention.current_order')
    try:
        assert_convention_review(manifest,current['result_payload'])
    except ValueError as exc:
        print(f'Publication blocked: {exc}',file=sys.stderr)
        return 2
    if review!=CONVENTION_REVIEW or review['publication_eligibility']!='READY_FOR_PUBLICATION_REVIEW':
        print('Publication blocked: corrected science review incomplete.',
              file=sys.stderr)
        return 2
    if args.require_publication_authorization and not review['publication_authorized']:
        print('Publication blocked: a separate remote-publication authorization is required.',file=sys.stderr)
        return 2
    print('Same-revision full science gate passed; remote publication remains separately authorized.')
    return 0

if __name__=='__main__':raise SystemExit(main())
