#!/usr/bin/env python3
"""Fail closed until same-revision full evidence and author resolution exist."""
from __future__ import annotations
import argparse
from pathlib import Path
import re
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from validation_evidence import (CONVENTION_REVIEW,EvidenceError,validate_run,
                                 assert_convention_review,git_provenance)  # noqa: E402


def require_commit_authorization(approval, approved_revision, current_revision, evidence_revision):
    """Check one invocation's approval without changing the scientific record."""
    if approval != 'true':
        raise ValueError('explicit publication approval is absent')
    if not approved_revision or re.fullmatch(r'[0-9a-f]{40}',approved_revision) is None:
        raise ValueError('a full approved commit SHA is required')
    if approved_revision != current_revision or approved_revision != evidence_revision:
        raise ValueError('approved, checked-out, and evidence revisions disagree')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir',type=Path)
    parser.add_argument('--require-publication-authorization',action='store_true')
    parser.add_argument('--publication-approval',choices=('true','false'),default='false')
    parser.add_argument('--approved-revision')
    args=parser.parse_args()
    if args.run_dir is None:
        print('Publication blocked: no same-revision full convention evidence.',file=sys.stderr)
        return 2
    revision,current_dirty=git_provenance()
    if current_dirty:
        print('Publication blocked: the current working tree is not clean.',file=sys.stderr)
        return 2
    if args.require_publication_authorization:
        try:
            require_commit_authorization(args.publication_approval,args.approved_revision,
                                         revision,revision)
        except ValueError as exc:
            print(f'Publication blocked: {exc}',file=sys.stderr)
            return 2
    try:
        manifest,document=validate_run(args.run_dir)
    except (EvidenceError,OSError,ValueError) as exc:
        print(f'Publication blocked: {exc}',file=sys.stderr)
        return 2
    if manifest['profile']!='full' or manifest['dirty_state']:
        print('Publication blocked: a clean same-revision full profile is required.',file=sys.stderr)
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
    if args.require_publication_authorization:
        try:
            require_commit_authorization(args.publication_approval,args.approved_revision,
                                         revision,manifest['source_revision'])
        except ValueError as exc:
            print(f'Publication blocked: {exc}',file=sys.stderr)
            return 2
        print('Same-revision full science gate and commit-scoped publication approval passed.')
        return 0
    print('Same-revision full science gate passed; deployment requires separate approval.')
    return 0

if __name__=='__main__':raise SystemExit(main())
