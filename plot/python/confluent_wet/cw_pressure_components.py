"""Plot native sigma_bulk and pressure_lb PDFs from completed pressure-scan frames."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, ScalarFormatter


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def extract(task):
    a, rows, raw, old, manifest = task
    arrays = {'sigmaB': [], 'pLB': []}
    proofs = {}
    for row in rows:
        root = raw/row['case']
        previous = json.loads((old/row['case']/'pressure_result.json').read_text())
        if previous['run_dat_sha256'] != row['run_dat_sha256']:
            raise ValueError('Runcard provenance changed')
        for step in range(manifest['steady_start'], manifest['nsteps']+1, manifest['ninfo']):
            path = root/f'frame{step}.json'
            blob = path.read_bytes()
            sha = hashlib.sha256(blob).hexdigest()
            if sha != previous['frame_sha256'][path.name]:
                raise ValueError(f'Frame changed: {path}')
            data = json.loads(blob)['data']
            for key, field in [('sigmaB', 'sigma_bulk'), ('pLB', 'pressure_lb')]:
                arr = np.asarray(data[field]['value'], dtype=float).ravel()
                if arr.size != 256**2 or not np.isfinite(arr).all():
                    raise ValueError(f'Invalid {field}: {path}')
                arrays[key].append(arr)
            proofs[f'{row["case"]}/{path.name}'] = sha
    result = {'activity': a, 'source_sha256': proofs}
    for key in arrays:
        values = np.concatenate(arrays[key])
        if key == 'sigmaB' and values.min() < 0:
            raise ValueError('Negative sigma_bulk')
        density, edges = np.histogram(values, bins=240, density=True)
        if not np.isclose(np.sum(density*np.diff(edges)), 1.):
            raise ValueError('PDF normalization failed')
        result[key] = dict(edges=edges.tolist(), density=density.tolist(),
                           mean=float(values.mean()), sigma=float(values.std()), samples=len(values))
    return result


def plot(rows, out):
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False,
                         'axes.spines.right': False, 'pdf.fonttype': 42})
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    for row, color in zip(rows, plt.colormaps['viridis'](np.linspace(.08, .88, len(rows)))):
        for ax, key in zip(axes, ['sigmaB', 'sigmaB', 'pLB']):
            edges = np.asarray(row[key]['edges'])
            ax.plot(.5*(edges[1:]+edges[:-1]), row[key]['density'],
                    color=color, lw=1.7, label=f'{row["activity"]:g}')
    for ax, title, xlabel in zip(axes, [r'$\sigma_B$', r'$\sigma_B$ (log density)', r'$p_{\mathrm{LB}}$'],
                                 [r'$\sigma_B$', r'$\sigma_B$', r'$p_{\mathrm{LB}}$']):
        ax.set(title=title, xlabel=xlabel, ylabel='Probability density')
        ax.xaxis.set_major_locator(MaxNLocator(5))
        fmt=ScalarFormatter(useMathText=True); fmt.set_powerlimits((-3,3))
        ax.xaxis.set_major_formatter(fmt); ax.margins(x=.025)
    axes[1].set_yscale('log')
    for ax in [axes[0], axes[2]]:
        ax.set_ylim(bottom=0)
        ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,3), useMathText=True)
    h,l=axes[0].get_legend_handles_labels()
    fig.legend(h,l,title=r'$a=\zeta_{\mathrm{eff}}/\zeta$',ncol=6,loc='upper center',
               bbox_to_anchor=(.5,1.02),frameon=False,columnspacing=2.)
    fig.tight_layout(rect=(0,0,1,.84),w_pad=2.)
    for ext in ['png','pdf']:
        fig.savefig(out/f'sigmaB_pLB_distributions.{ext}',dpi=220,bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('raw',type=Path);p.add_argument('old',type=Path);p.add_argument('out',type=Path)
    p.add_argument('--workers',type=int,default=6)
    args=p.parse_args()
    manifest_path=Path(__file__).resolve().parents[3]/'cases/20260915/cw_pressure_compare/manifest.json'
    m=json.loads(manifest_path.read_text())
    activities=sorted({r['activity'] for r in m['cases']})
    tasks=[(a,[r for r in m['cases'] if r['activity']==a],args.raw,args.old,m) for a in activities]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        rows=list(pool.map(extract,tasks))
    args.out.mkdir(parents=True,exist_ok=True)
    result={'manifest_sha256':digest(manifest_path),'script_sha256':digest(Path(__file__)),
            'source':'Direct serialized sigma_bulk and pressure_lb; no pressure reconstruction.',
            'window_steps':[m['steady_start'],m['nsteps']], 'frames_per_seed':41,'seeds_per_activity':3,
            'histogram':'240 full-range bins per activity and field; equal site/time/seed weights; no smoothing or centering.',
            'activities':rows}
    (args.out/'component_distributions.json').write_text(json.dumps(result,indent=2)+'\n')
    plot(rows,args.out)
    print(json.dumps({'activities':len(rows),'verified_frames':sum(len(r['source_sha256']) for r in rows),
                      'output':str(args.out)}))
