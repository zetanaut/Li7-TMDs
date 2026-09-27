"""Plot Born hard-response coefficients and validation diagnostics from the grid JSON."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]


def save(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    for suffix in ('.png', '.svg'):
        fig.savefig(path.with_suffix(suffix), dpi=180, bbox_inches='tight')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=ROOT/'results/gluon_born_grid_report.json')
    parser.add_argument('--output-dir', type=Path, default=ROOT/'docs/figures')
    args = parser.parse_args()
    report = json.loads(args.input.read_text())
    if report['points'] != 36 or not report['all_pass']:
        raise ValueError('Expected a validated 36-point grid')
    rows = report['results']
    thetas = sorted({r['inputs']['theta'] for r in rows})
    phis = sorted({r['inputs']['phi'] for r in rows})
    for key in ('b_U', 'b_G', 'b_C', 'b_S'):
        fig, ax = plt.subplots(figsize=(6, 4))
        for phi in phis:
            selected = [next(r for r in rows if r['inputs']['theta'] == t
                             and r['inputs']['phi'] == phi
                             and r['inputs']['helicity'] == 1.0) for t in thetas]
            ax.plot(thetas, [r['Stokes'][key] for r in selected], marker='o', label=f'φ={phi:g}')
        ax.set(xlabel='θ (rad)', ylabel=key, title=f'Born-level hard response: {key}')
        ax.legend(); ax.grid(alpha=.25)
        save(fig, args.output_dir/key)
    fig, ax = plt.subplots(figsize=(6, 4))
    for j in (0, 1):
        selected = [next(r for r in rows if r['inputs']['theta'] == t
                         and r['inputs']['phi'] == 0.4
                         and r['inputs']['helicity'] == 1.0) for t in thetas]
        ax.plot(thetas, [r['B_eigenvalues'][j] for r in selected], marker='o', label=f'eigenvalue {j+1}')
    ax.set(xlabel='θ (rad)', ylabel='eigenvalue', title='Born-level hard-matrix eigenvalues (φ=0.4, λ=1)')
    ax.legend(); ax.grid(alpha=.25)
    save(fig, args.output_dir/'eigenvalues')
    fig, ax = plt.subplots(figsize=(7, 4))
    for key in ('photon', 'gluon'):
        ax.semilogy(range(len(rows)), [max(r['relative_Ward_residuals'][key], 1e-18) for r in rows],
                    marker='.', linestyle='none', label=key)
    ax.axhline(1e-11, color='black', linestyle='--', label='acceptance tolerance')
    ax.set(xlabel='grid point index', ylabel='relative residual', title='Born-amplitude Ward validation (36 points)')
    ax.legend(); ax.grid(alpha=.25)
    save(fig, args.output_dir/'ward_residuals')
    print(f'Wrote six figures as PNG and SVG to {args.output_dir}')


if __name__ == '__main__':
    main()
