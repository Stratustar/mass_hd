#!/usr/bin/env python3
"""Generate and inspect small pulse-integrator and output-invariance checks."""
import argparse
import json
from pathlib import Path
import numpy as np

NAMES = ('analytic', 'disabled', 'noop', 'diagnostics_off', 'diagnostics_on', 'production_pulse')


def generate(out):
    repo = Path(__file__).resolve().parents[1]
    p = {}
    baseline = repo/'cases/20260913/cw_pmem_pulse/controls/L256/tm1_initpatch_rep1/run.dat'
    for line in baseline.read_text().splitlines():
        if line.strip() and not line.startswith('#'):
            k,v = line.split('=',1); p[k.strip()] = v.strip()
    for name in NAMES:
        q = dict(p)
        q.update(LX=16, LY=16, nsteps=100, ninfo=100, seed=7351, noise=0.,
            **{'chi-config':'uniform','chi0':.5,'m0':.2,'tau-m':10.,'tau-chi':3.,
            'chi-freeze-steps':0,'chi-length':0.,'pmem-pulse-start':17,'pmem-pulse-steps':0,
            'pmem-pulse-value':-.01,'nresponse':7,'nvideo':11,'video-stride':1,
            'video-start':5,'nvideo-dense':3,'video-dense-start':5,'video-dense-end':40,
            'frame-light':0})
        if name == 'analytic':q['pmem-pulse-steps']=11
        if name == 'noop':q.update({'pmem-pulse-steps':11,'pmem-pulse-value':.016838})
        if name in ('diagnostics_off','diagnostics_on','production_pulse'):
            q.update(LX=32,LY=32,nsteps=1000,ninfo=1000,noise=.05)
            q.update({'tau-m':100.,'tau-chi':30.,'pmem-pulse-start':317,
                'nresponse':13 if name!='diagnostics_off' else 0,'nvideo':0,'nvideo-dense':0,
                'video-start':0,'video-dense-start':0,'video-dense-end':0})
        if name == 'production_pulse':q.update({'pmem-pulse-value':0.,'pmem-pulse-steps':202})
        path=out/name/'run.dat';path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text('# Local validation; not a production physics result.\n'+
                       ''.join(f'{k} = {v}\n' for k,v in q.items()))
    print(json.dumps({'generated':list(NAMES),'out':str(out)}))


def frame(root):
    path=sorted(root.glob('frame*.json'))[-1]
    j=json.loads(path.read_text())['data']
    return {k:np.asarray(v['value']) for k,v in j.items()}


def check(root):
    checks={}
    a=np.genfromtxt(root/'analytic/response.csv',delimiter=',',names=True)
    t=a['step'].astype(int)
    expected_t=sorted(set(range(0,101,7))|{0,17,28,100})
    assert t.tolist()==expected_t
    active=(t>=17)&(t<28)
    assert np.array_equal(a['pulse_on'],active)
    assert np.array_equal(a['pmem_effective'],np.where(active,-.01,.016838))
    m=np.empty(101);m[0]=.2;decay=1-.1+.5*.1**2
    for k in range(100):
        source=float(17<=k<28);m[k+1]=source+(m[k]-source)*decay
    error=float(np.max(np.abs(a['m_mean']-m[t])))
    assert error<2e-12,error
    assert np.max(abs(a['P_mean']))<1e-12 and np.max(abs(a['u_rms']))<1e-12
    checks['heun_memory_max_error']=error
    checks['exact_event_steps']=t.tolist()
    for left,right in [('disabled','noop'),('diagnostics_off','diagnostics_on')]:
        x,y=frame(root/left),frame(root/right)
        fields=('QQxx','QQyx','chi','m','ff')
        errors={k:float(np.max(np.abs(x[k]-y[k]))) for k in fields}
        assert all(v==0 for v in errors.values()),errors
        checks[left+'_vs_'+right]=errors
    v=np.genfromtxt(root/'analytic/video_meta.csv',delimiter=',',names=True,comments='#',skip_header=3)
    expected_video=[]
    for step in range(101):
        if step<5:continue
        if step in (5,17,28,100):due=True
        elif 5<=step<=40:due=(step-5)%3==0 or step==40
        else:due=(step-5)%11==0
        if due:expected_video.append(step)
    assert v['t'].astype(int).tolist()==expected_video
    for field in ('u','P','m','chi'):
        assert (root/'analytic'/f'video_{field}.u8').stat().st_size==len(v)*16*16
    checks['video_steps']=expected_video
    p=np.genfromtxt(root/'production_pulse/response.csv',delimiter=',',names=True)
    assert p['step'][-1]==1000 and {317,519}<=set(p['step'])
    assert np.all(p['pmem_effective'][(p['step']>=317)&(p['step']<519)]==0)
    assert np.all(p['pmem_effective'][(p['step']<317)|(p['step']>=519)]==.016838)
    checks['production_threshold_zero']=True
    out=root/'validation.json';out.write_text(json.dumps(checks,indent=2)+'\n')
    print(json.dumps(checks))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['generate','check']);p.add_argument('root',type=Path)
    a=p.parse_args();(generate if a.mode=='generate' else check)(a.root)
