"""Reject manuscript, handoff, and agent artifacts in Git or a built site."""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
FORBIDDEN_SUFFIXES = ('.tex', '.bib', '.pdf', '.zip', '.tar', '.tgz', '.gz', '.7z')
FORBIDDEN_PARTS = ('reference', 'instructions', 'handoff', 'manuscript', 'transcript')
FORBIDDEN_NAMES = ('codex_', 'prompt', 'chat_export')


def forbidden(path: str, *, site: bool = False) -> bool:
    parts = Path(path).parts
    low = path.lower()
    if site and low == 'sitemap.xml.gz':
        return False
    return (low.endswith(FORBIDDEN_SUFFIXES) or
            any(part.lower() in FORBIDDEN_PARTS for part in parts) or
            any(token in low for token in FORBIDDEN_NAMES))


def forbidden_bytes(data: bytes) -> bool:
    """Catch renamed manuscript containers and obvious TeX/BibTeX sources."""
    if (data.startswith((b'%PDF-', b'PK\x03\x04', b'\x1f\x8b', b'7z\xbc\xaf\x27\x1c')) or
            data[257:262] == b'ustar'):
        return True
    head = data[:4096].lstrip()
    return head.startswith((b'\\documentclass', b'@article{', b'@book{', b'@misc{'))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--site', type=Path, help='Inspect an actual MkDocs output tree.')
    args = parser.parse_args()
    tracked = subprocess.check_output(['git', 'ls-files', '--cached', '-z'], cwd=ROOT)
    paths = [value.decode('utf-8', 'surrogateescape') for value in tracked.split(b'\0') if value]
    offenders = [f'Git: {path}' for path in paths if forbidden(path)]
    staged = subprocess.check_output(['git', 'ls-files', '--stage', '-z'], cwd=ROOT)
    for entry in staged.split(b'\0'):
        if not entry:
            continue
        metadata, raw_path = entry.split(b'\t', 1)
        mode = metadata.split(b' ', 1)[0]
        path = raw_path.decode('utf-8', 'surrogateescape')
        if mode == b'120000':
            offenders.append(f'Git symlink: {path}')
        elif mode in (b'100644', b'100755'):
            content = subprocess.check_output(['git', 'show', ':' + path], cwd=ROOT)
            if forbidden_bytes(content):
                offenders.append(f'Git content: {path}')
    if args.site:
        if not args.site.is_dir():
            offenders.append('Site output is missing')
        else:
            for path in args.site.rglob('*'):
                if path.is_symlink():
                    offenders.append(f'Site symlink: {path.relative_to(args.site)}')
                elif forbidden(str(path.relative_to(args.site)), site=True):
                    offenders.append(f'Site: {path.relative_to(args.site)}')
                elif path.is_file() and path.name != 'sitemap.xml.gz' and forbidden_bytes(path.read_bytes()):
                    offenders.append(f'Site content: {path.relative_to(args.site)}')
    if offenders:
        print('Public scope rejected:\n' + '\n'.join(offenders), file=sys.stderr)
        return 1
    print(f'Public scope PASS: {len(paths)} cached Git paths' +
          (f'; site {args.site}' if args.site else ''))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
