#!/usr/bin/env python3
"""Compare moment-matched Gaussian and Beta closures on unquantised m snapshots.

Spatial column-block bootstrap retains ALL selected times in each resampled block.
Intervals are approximate spatial-ensemble uncertainty conditional on the saved
window, not independent-site KS p-values and not evidence of absolute model truth.
Main effect: D_G-D_B, positive favours Beta. Refit both moments on every bootstrap.
Cross-validation: five contiguous x strips, 25-cell periodic training buffer;
held-out binned log score, Beta minus Gaussian, without Gaussian renormalisation.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import ndtr, betainc, betaincc

TC = 674.3290333006435
EDGE = np.unique(np.r_[0., np.geomspace(1e-8,.01,41), np.linspace(0,1,1025),
                       1-np.geomspace(1e-8,.01,41), .21, 1.])
SCORE_EDGE = np.unique(np.r_[0., np.geomspace(1e-8,.01,21), np.linspace(0,1,257),
                             1-np.geomspace(1e-8,.01,21), 1.])


def params(mom):
    mu=mom[...,0]
    var=np.maximum(mom[...,1]-mu*mu,1e-16)
    concentration=mu*(1-mu)/var-1
    if np.any(concentration<=0):
        raise ValueError('Moment-matched continuous Beta undefined at Bernoulli limit')
    return mu,np.sqrt(var),mu*concentration,(1-mu)*concentration


def cdfs(mom,edges=EDGE):
    mu,sd,a,b=params(mom)
    return ndtr((edges-mu[...,None])/sd[...,None]), betainc(a[...,None],b[...,None],edges)


def log_bins(mom):
    mu,sd,a,b=params(mom)
    z=(SCORE_EDGE-mu[...,None])/sd[...,None]
    gc,gs=ndtr(z),ndtr(-z)
    bc=betainc(a[...,None],b[...,None],SCORE_EDGE)
    bs=betaincc(a[...,None],b[...,None],SCORE_EDGE)
    def diff(c,s):
        p=np.where(c[...,1:]<.5,np.diff(c,axis=-1),-np.diff(s,axis=-1))
        return np.log(np.maximum(p,1e-300))
    return diff(gc,gs),diff(bc,bs)


def dists(H,mom):
    G,B=cdfs(mom)
    # Include 0-left and 1-left, so point masses and Gaussian leakage are not hidden.
    def dist(F):
        return np.maximum.reduce([np.max(abs(H-F),axis=-1), F[...,0],
                                  abs(1-mom[...,5]-F[...,-1])])
    return dist(G),dist(B),G,B


def aggregate(a,side):
    # a starts as a 20x20 grid of 25-cell blocks, with trailing features.
    k=side//25
    n=20//k
    shape=(n,k,n,k)+a.shape[2:]
    return a.reshape(shape).mean(axis=(1,3)).reshape((n*n,)+a.shape[2:])


def interval(values,point):
    lo,hi=np.quantile(values,[.025,.975])
    # Centred bootstrap test of zero difference; approximate, not a KS null p-value.
    p=(1+np.count_nonzero(abs(values-point)>=abs(point)))/(len(values)+1)
    return {'estimate':float(point),'ci95':[float(lo),float(hi)],'p_approx':float(p)}


def bootstrap(H,M,reps,rng):
    n=len(M)
    DG,DB,_,_=dists(H.mean(0),M.mean(0))
    mc_i=int(np.flatnonzero(EDGE==.21)[0])
    values=[]; mc_values=[]
    for start in range(0,reps,100):
        q=min(100,reps-start)
        w=rng.multinomial(n,np.full(n,1/n),size=q)/n
        mm=w@M; hh=w@H
        dg,db,g,b=dists(hh,mm)
        values.extend(dg-db)
        mc_values.extend(abs(g[:,mc_i]-mm[:,6])-abs(b[:,mc_i]-mm[:,6]))
    g,b=cdfs(M.mean(0)); empirical=M.mean(0)[6]
    out=interval(np.asarray(values),DG-DB)
    out['blocks']=n
    out['mc_abs_error_difference']=interval(np.asarray(mc_values),
            abs(g[mc_i]-empirical)-abs(b[mc_i]-empirical))
    return out


def crossfit(M,S,reps,rng):
    # Use 25-cell base blocks. Five 100-cell wide x strips, buffered by one base block.
    x=np.repeat(np.arange(20),20)
    sets=[]
    for f in range(5):
        hold=(x//4==f)
        excluded=np.array([(j%20) for j in range(4*f-1,4*f+5)])
        train=~np.isin(x,excluded)
        sets.append((hold,train))
    def score(w):
        delta=np.zeros(len(w)); total=np.zeros(len(w))
        for hold,train in sets:
            tm=w[:,train]@M[train]/w[:,train].sum(1)[:,None]
            lg,lb=log_bins(tm)
            hist=w[:,hold]@S[hold]
            delta+=np.sum(hist*(lb-lg),axis=1)
            total+=w[:,hold].sum(1)
        return delta/total
    point=score(np.ones((1,400)))[0]
    # Resample 50-cell column blocks, propagate weights down to the 25-cell base grid.
    vals=[]
    for start in range(0,reps,100):
        q=min(100,reps-start)
        counts=rng.multinomial(100,np.full(100,.01),size=q).reshape(q,10,10)
        w=np.repeat(np.repeat(counts,2,axis=1),2,axis=2).reshape(q,400)
        vals.extend(score(w))
    return interval(np.asarray(vals),point)


def one_case(src,out,reps=1999):
    import cw_common as cw
    from archive.archive import loadarchive
    from cw_s3_sk import radial_corr,length_1e
    par=cw.read_params(str(src)); oa=loadarchive(str(src))
    L=int(par['LX']); assert L==500 and int(par['LY'])==500
    # A final off-cadence frame may also be saved. Counting filenames would turn
    # it into a nonexistent regular frame; retain only the common ninfo grid.
    idx=[i for i in range(cw.ph.frame_count(oa))
         if cw.ph.frame_time(oa,i)>=par['nsteps']-100*TC
         and (src/f'frame{oa.nstart+i*oa.ninfo}.json').exists()]
    assert len(idx)>=6
    rng=np.random.default_rng(84621+int(par.get('seed',0)))
    ox,oy=rng.integers(0,25,size=2)
    H=np.zeros((20,20,len(EDGE))); S=np.zeros((20,20,len(SCORE_EDGE)-1))
    M=np.zeros((20,20,7)); per=[]; hs=[]; power=np.zeros((L,L))
    for i in idx:
        frame=oa.read_frame(i)
        m=np.asarray(cw.ph.grid(frame,'m'),float)
        chi=np.asarray(cw.ph.grid(frame,'chi'),float)
        assert np.isfinite(m).all() and m.min()>=0 and m.max()<=1
        power+=abs(np.fft.fft2(m-m.mean()))**2/(L*L)**2
        # CDF uses <= with searchsorted, including any atoms exactly on grid edges.
        hh=np.bincount(np.searchsorted(EDGE,m.ravel(),side='left'),minlength=len(EDGE)+1)[:-1].cumsum()/m.size
        hs.append(hh)
        per.append({'t':cw.ph.frame_time(oa,i)/TC,'mean':float(m.mean()),'std':float(m.std()),
                    'cdf_mc':float((m<.21).mean()),'chi':float(chi.mean())})
        rolled=np.roll(m,(int(ox),int(oy)),axis=(0,1))
        for ix in range(20):
            for iy in range(20):
                v=rolled[ix*25:(ix+1)*25,iy*25:(iy+1)*25].ravel()
                count=np.bincount(np.searchsorted(EDGE,v,side='left'),minlength=len(EDGE)+1)[:-1]
                H[ix,iy]+=count.cumsum()/len(v)
                count=np.bincount(np.searchsorted(SCORE_EDGE,v,side='left'),minlength=len(SCORE_EDGE)+1)[:-1]
                count[1]+=count[0] # first [0,e1] bin includes the zero atom; no density clipping
                S[ix,iy]+=count[1:]/len(v)
                M[ix,iy]+=[np.mean(v),np.mean(v*v),np.mean(v**3),np.mean(v**4),
                            np.mean(v==0),np.mean(v==1),np.mean(v<.21)]
    H/=len(idx); S/=len(idx); M/=len(idx)
    mm=M.mean((0,1)); empirical=H.mean((0,1)); dg,db,g,b=dists(empirical,mm)
    mu,sd,alpha,beta=params(mm)
    half=len(per)//2
    early=np.mean([[p['mean'],p['std'],p['chi']] for p in per[:half]],axis=0)
    late=np.mean([[p['mean'],p['std'],p['chi']] for p in per[half:]],axis=0)
    drift=(late-early)
    ks_time=float(np.max(abs(np.mean(hs[:half],axis=0)-np.mean(hs[half:],axis=0))))
    flags=[]
    if abs(drift[0])/sd>.1: flags.append('mean drift > 0.1 pooled std')
    if abs(drift[1])/sd>.1: flags.append('std drift > 0.1 pooled std')
    if abs(drift[2])>.02: flags.append('chi drift > 0.02')
    if ks_time>.05: flags.append('early/late CDF distance > 0.05')
    corr=radial_corr(power/len(idx),L)
    mi=int(np.flatnonzero(EDGE==.21)[0])
    out.mkdir(parents=True,exist_ok=True)
    info={'case':src.name,'study':src.parent.name,'tau_m':float(par['tau_m'])/TC,
          'chi0':float(par['chi0']),'mc':.21,'window':[per[0]['t'],per[-1]['t']],
          'frames':len(idx),'sites_per_frame':L*L,'moments':mm.tolist(),
          'mean':float(mu),'std':float(sd),'alpha':float(alpha),'beta':float(beta),
          'skewness':float((mm[2]-3*mu*mm[1]+2*mu**3)/sd**3),
          'excess_kurtosis':float((mm[3]-4*mu*mm[2]+6*mu**2*mm[1]-3*mu**4)/sd**4-3),
          'ks_gaussian':float(dg),'ks_beta':float(db),'cdf_mc':float(mm[6]),
          'cdf_mc_gaussian':float(g[mi]),'cdf_mc_beta':float(b[mi]),
          'gaussian_outside_support':float(ndtr(-mu/sd)+ndtr((mu-1)/sd)),
          'edge_zero':float(mm[4]),'edge_one':float(mm[5]),
          'temporal_ks':ks_time,'drift_mean_over_sd':float(drift[0]/sd),
          'drift_chi':float(drift[2]),'nonstationary_flags':flags,
          'm_corr_1e':length_1e(corr),'m_corr_50':float(corr[50]),
          'per_frame':per,'edges':EDGE.tolist(),'cdf':empirical.tolist(),
          'gaussian_cdf':g.tolist(),'beta_cdf':b.tolist(),'bootstrap':{},
          'bootstrap_reps':reps,'block_offset':[int(ox),int(oy)]}
    for size in (25,50,100):
        info['bootstrap'][str(size)]=bootstrap(aggregate(H,size),aggregate(M,size),reps,rng)
    info['crossfit_logscore']=crossfit(M.reshape(400,7),S.reshape(400,-1),reps,rng)
    # Finite-window descriptors per frame expose mixing of evolving distributions.
    np.savez_compressed(out/'blocks.npz',H=H,M=M,S=S,per_cdf=np.array(hs))
    (out/'distribution.json').write_text(json.dumps(info))
    print(f'{src.name}: mean={mu:.5f} sd={sd:.5f} D(G/B)={dg:.4f}/{db:.4f} '
          f'CI50={info["bootstrap"]["50"]["ci95"]} flags={flags}',flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('inputdir',type=Path);ap.add_argument('outdir',type=Path)
    ap.add_argument('--reps',type=int,default=1999)
    ap.add_argument('--old-controls',action='store_true')
    args=ap.parse_args()
    if args.old_controls:
        for tag in ('0p3','3p43','9p68'):
            for start in ('chi0','chi1'):
                c=f'tm{tag}_{start}'
                one_case(args.inputdir/c,args.outdir/c,args.reps)
    else:
        one_case(args.inputdir,args.outdir,args.reps)
if __name__=='__main__':main()
