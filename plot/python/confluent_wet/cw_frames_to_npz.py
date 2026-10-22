#!/usr/bin/env python3
"""Pack the dense full frames of a confluent-wet run into one compressed npz (for defect tracking).

  cw_frames_to_npz.py <run_dir> <out.npz> [--compare <original_run_dir>] [--fields chi m ...]

All frame*.json in <run_dir>, sorted by step, are stored as float16 arrays of shape (T, LX, LY)
(axis 0 = x, as plot_hd reshapes (LX, LY)) under the frame field names, plus `t` (steps).
float16 keeps ~3 significant digits: enough for chi, m in [0, 1], Q, u ~ 1e-2 and P ~ 2e-2
(spacing ~1.5e-5 at the threshold .0168). With --compare, the last frame is checked against the
frame with the same step in the original run (max |difference| per field, printed and stored in
`compare_maxabs`), to confirm that a run differing only in output settings is the same trajectory.
"""
import argparse
import json
from pathlib import Path
import re

import numpy as np


def load(path):
    d = json.loads(path.read_text())['data']
    return {k: np.asarray(v['value'], float) for k, v in d.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run', type=Path)
    ap.add_argument('out', type=Path)
    ap.add_argument('--compare', type=Path)
    ap.add_argument('--fields', nargs='*', help='keep only these frame fields (default: all)')
    a = ap.parse_args()
    p = json.loads((a.run/'parameters.json').read_text()); p = p.get('data', p)
    get = lambda k: p[k]['value'] if isinstance(p.get(k), dict) else p.get(k)
    LX, LY = int(get('LX')), int(get('LY'))
    files = sorted(a.run.glob('frame*.json'), key=lambda f: int(re.search(r'frame(\d+)', f.name).group(1)))
    steps = np.array([int(re.search(r'frame(\d+)', f.name).group(1)) for f in files])
    first = load(files[0])
    if a.fields:
        first = {k: first[k] for k in a.fields}
    data = {k: np.empty((len(files), LX, LY), np.float16) for k in first}
    for i, f in enumerate(files):
        fr = first if i == 0 else load(f)
        for k in data:
            data[k][i] = fr[k].reshape(LX, LY).astype(np.float16)
    out = dict(data, t=steps)
    if a.compare:
        ref = a.compare/files[-1].name
        last = load(files[-1]); orig = load(ref)
        diff = {k: float(np.max(np.abs(last[k] - orig[k]))) for k in last}
        out['compare_maxabs'] = np.array([diff[k] for k in sorted(diff)])
        out['compare_fields'] = np.array(sorted(diff))
        print(json.dumps({'compare_frame': str(ref), 'max_abs_diff': diff}))
    a.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(a.out, **out)
    print(json.dumps({'run': str(a.run), 'frames': len(files), 'first_step': int(steps[0]), 'last_step': int(steps[-1]),
                      'out': str(a.out), 'MB': round(a.out.stat().st_size/2**20, 1)}))


if __name__ == '__main__':
    main()
