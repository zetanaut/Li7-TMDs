"""Plot Born hard-response coefficients and validation diagnostics from the grid JSON."""
import argparse
import json
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from validation_evidence import validate_run  # noqa: E402


def save(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    for suffix in ('.png', '.svg'):
        fig.savefig(path.with_suffix(suffix), dpi=180, bbox_inches='tight')
    svg=path.with_suffix('.svg')
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=ROOT/'results/gluon_born_grid_report.json')
    parser.add_argument('--run-dir', type=Path,
                        help='Use a complete current run instead of a historical grid JSON.')
    parser.add_argument('--output-dir', type=Path, default=ROOT/'docs/figures')
    args = parser.parse_args()
    if args.run_dir:
        _manifest, document = validate_run(args.run_dir)
        by_id={item['check_id']:item['result_payload'] for item in document['results']}
        grid = next(item['result_payload'] for item in document['results']
                    if item['check_id'] == 'grid.complete')
        report = {'points': grid['number_of_cases'], 'all_pass': True,
                  'results': grid['cases']}
    else:
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
    if args.run_dir and 'gluon.born.dense_scan' in by_id:
        dense=by_id['gluon.born.dense_scan']
        selected=[x for x in dense['cases'] if x['inputs']['helicity']==1.]
        angles=dense['angle_array']
        fig,ax=plt.subplots(figsize=(7,4))
        for name in ('b_G','b_C','b_S'):
            ax.plot(angles,[x['stokes'][name]/x['stokes']['b_U'] for x in selected],label=name+'/b_U')
        ax.set(xlabel='θ (rad)',ylabel='dimensionless analyzing ratio',
               title='Born hard analyzer, 161 calculated angles')
        ax.grid(alpha=.25);ax.legend();save(fig,args.output_dir/'dense_analyzing')
        fig,ax=plt.subplots(figsize=(7,4))
        for j in (0,1):
            ax.plot(angles,[x['eigenvalues'][j] for x in selected],label=f'eigenvalue {j+1}')
        ax.set(xlabel='θ (rad)',ylabel='reduced B eigenvalue (GeV²)',
               title='Born hard matrix, 161 calculated angles')
        ax.grid(alpha=.25);ax.legend();save(fig,args.output_dir/'dense_eigenvalues')
        fig,ax=plt.subplots(figsize=(7,4))
        floor=1e-18
        ax.semilogy(angles,[max(x['direct_residual_relative'],floor) for x in selected],label='direct/trace')
        for key in ('photon','gluon'):
            ax.semilogy(angles,[max(x['relative_wards'][key],floor) for x in selected],label=key+' Ward')
        ax.set(xlabel='θ (rad)',ylabel=f'relative residual (display floor {floor:g})',
               title='Born independent and Ward diagnostics')
        ax.grid(alpha=.25);ax.legend();save(fig,args.output_dir/'dense_residuals')
        recon=by_id['gluon.oct.reconstruction']
        fig,ax=plt.subplots(figsize=(5,4))
        keys=('condition_raw','condition_column_normalized')
        ax.bar(('physical columns','unit columns'),[recon[k] for k in keys])
        ax.set_yscale('log');ax.set(ylabel='2-norm condition number',
               title='Synthetic octupole response design')
        save(fig,args.output_dir/'octupole_condition')
        print(f'Wrote ten validated figures as PNG and SVG to {args.output_dir}')
    else:
        print(f'Wrote six figures as PNG and SVG to {args.output_dir}')


if __name__ == '__main__':
    main()
