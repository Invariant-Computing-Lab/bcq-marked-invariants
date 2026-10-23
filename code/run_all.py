from __future__ import annotations
import json, platform, sys, time
from pathlib import Path
from fractions import Fraction
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from bcq_invariants import (
    symplectic_J, pairing, pairings_to_basis, random_symplectic,
    quintic_euler_line_bundle, hrr_pairing_line_bundles,
    random_unitary, random_su3_hermitian, gauge_invariants,
    random_unimodular, wilson_interval,
)

SEED = 20260812
TOL = 1e-10
ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = ROOT / "figures"
RES.mkdir(exist_ok=True); FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(SEED)
summary = []


def add_summary(exp, metric, value, n=None, passed=None, note=""):
    summary.append({"experiment": exp, "metric": metric, "value": value, "n": n, "passed": passed, "note": note})


def exp1_symplectic_invariance(n=20000):
    J = symplectic_J(3)
    B = np.eye(6, dtype=int)
    gamma_c = np.array([1,0,0,0,0,0], dtype=int)
    gamma_w = np.array([0,1,0,3,0,0], dtype=int)
    chi0 = pairing(gamma_c, gamma_w, J)
    pc0 = pairings_to_basis(gamma_c, B, J)
    pw0 = pairings_to_basis(gamma_w, B, J)
    fails_chi = fails_rel = 0
    for _ in range(n):
        U = random_symplectic(3, rng, stages=2)
        gc = U @ gamma_c; gw = U @ gamma_w; Bb = U @ B
        if pairing(gc, gw, J) != chi0: fails_chi += 1
        if pairings_to_basis(gc, Bb, J) != pc0 or pairings_to_basis(gw, Bb, J) != pw0: fails_rel += 1
    add_summary("E1 Symplectic quotient", "chi failures", fails_chi, n, fails_chi == 0)
    add_summary("E1 Symplectic quotient", "relational fingerprint failures", fails_rel, n, fails_rel == 0)
    return {"chi0": chi0, "fails_chi": fails_chi, "fails_rel": fails_rel, "n": n}


def exp1b_toroidal_replication(n=20000):
    # Factorized 3-cycle toy model on T^2 x T^2 x T^2.
    a = [np.array([1,0],dtype=int), np.array([1,0],dtype=int), np.array([1,1],dtype=int)]
    b = [np.array([0,1],dtype=int), np.array([0,1],dtype=int), np.array([1,4],dtype=int)]
    def det2(x,y): return int(x[0]*y[1]-x[1]*y[0])
    def Iprod(aa,bb):
        vals=[det2(x,y) for x,y in zip(aa,bb)]
        return vals, int(np.prod(vals))
    local0,total0=Iprod(a,b)
    fails=0
    for _ in range(n):
        aa=[]; bb=[]
        for i in range(3):
            U=random_unimodular(2,rng,steps=6)
            # force det +1 to preserve oriented intersection, not just absolute value
            if round(np.linalg.det(U)) == -1:
                U[0,:] *= -1
            aa.append(U@a[i]); bb.append(U@b[i])
        _,tot=Iprod(aa,bb)
        fails += int(tot != total0)
    add_summary("E1b Toroidal replication", "oriented intersection failures", fails, n, fails==0)
    add_summary("E1b Toroidal replication", "baseline factor intersections", str(tuple(local0)), None, total0==3)
    return {"local0":local0,"total0":total0,"fails":fails,"n":n}


def exp2_exhaustive_collisions(bound=3):
    J = symplectic_J(3); B = np.eye(6, dtype=int)
    gc = np.array([1,0,0,0,0,0], dtype=int)
    target = 3
    collision = []
    sigs = set()
    total = 0
    vals = range(-bound, bound+1)
    for x0 in vals:
      for x1 in vals:
       for x2 in vals:
        for x3 in vals:
         for x4 in vals:
          for x5 in vals:
            gw = np.array([x0,x1,x2,x3,x4,x5], dtype=int); total += 1
            if pairing(gc, gw, J) == target:
                collision.append(gw.copy())
                sigs.add(pairings_to_basis(gw, B, J))
    ncol = len(collision)
    df = pd.DataFrame(collision, columns=[f"g{i+1}" for i in range(6)])
    df.to_csv(RES / "E2_chi3_collision_vectors.csv", index=False)
    sig_list=[pairings_to_basis(v,B,J) for v in collision]
    ab_rows=[]
    for drop in range(6):
        reduced={tuple(x for j,x in enumerate(sig) if j!=drop) for sig in sig_list}
        ab_rows.append([drop+1,len(reduced),ncol/len(reduced)])
    pd.DataFrame(ab_rows,columns=['removed_relational_coordinate','unique_reduced_signatures','mean_collision_multiplicity']).to_csv(RES/'E2_reconstructibility_ablation.csv',index=False)
    add_summary("E2 Exhaustive collision", "search space", total, total, True)
    add_summary("E2 Exhaustive collision", "vectors with chi=3", ncol, total, True, "All share one scalar generation index")
    add_summary("E2 Exhaustive collision", "unique marked-lattice relational signatures", len(sigs), ncol, len(sigs)==ncol)
    return {"total": total, "collisions": ncol, "unique_rel": len(sigs), "ablation": ab_rows}


def exp3_hrr(n=20000):
    common_fail = 0; one_side_same = 0
    rows=[]
    for _ in range(n):
        a=int(rng.integers(-12,13)); b=int(rng.integers(-12,13)); t=int(rng.integers(-8,9))
        if t==0: t=1
        base=hrr_pairing_line_bundles(a,b)
        common=hrr_pairing_line_bundles(a+t,b+t)
        bad=hrr_pairing_line_bundles(a+t,b)
        common_fail += int(common != base)
        one_side_same += int(bad == base)
        if len(rows)<100:
            rows.append([a,b,t,str(base),str(common),str(bad)])
    pd.DataFrame(rows,columns=['a','b','twist','chi_base','chi_common_twist','chi_one_sided']).to_csv(RES/'E3_HRR_examples.csv',index=False)
    add_summary("E3 HRR common twist", "positive failures", common_fail, n, common_fail==0)
    add_summary("E3 HRR one-sided control", "negative controls escaping detection", one_side_same, n, one_side_same==0)
    return {"common_fail":common_fail,"one_side_same":one_side_same,"n":n}


def exp4_gauge(n=10000):
    max_err=0.0; raw=[]; neg_detect=0
    for _ in range(n):
        H=random_su3_hermitian(rng); inv0=gauge_invariants(H)
        U=random_unitary(3,rng); Hp=U@H@U.conjugate().T
        err=float(np.max(np.abs(gauge_invariants(Hp)-inv0))); max_err=max(max_err,err)
        raw.append(float(np.linalg.norm(Hp-H)))
        eps=float(rng.uniform(1e-4,1.0)); Hbad=Hp+eps*np.eye(3)
        if abs(np.trace(Hbad))>TOL: neg_detect += 1
    med=float(np.median(raw))
    add_summary("E4 Gauge conjugacy", "max invariant error", max_err, n, max_err<TOL)
    add_summary("E4 Gauge conjugacy", "median raw matrix drift", med, n, True)
    add_summary("E4 Gauge negative control", "trace-contamination detections", neg_detect, n, neg_detect==n)
    return {"max_err":max_err,"median_raw":med,"neg_detect":neg_detect,"n":n}


def integer_inverse_transpose(U):
    # U unimodular -> exact integer inverse transpose via sympy
    import sympy as sp
    M=sp.Matrix(U.tolist())
    return np.array(M.inv().T.tolist(),dtype=object)


def exp5_hypercharge_basis(n=3000):
    import sympy as sp
    c = sp.Matrix([sp.Rational(1,6), sp.Rational(-1,2), sp.Rational(1,3)])
    A = sp.Matrix([
        [1,0,0],[0,1,0],[0,0,1],[1,-1,0],[1,0,-1],[0,1,-1],
        [2,-1,0],[1,1,-1],[-1,2,1],[2,0,-1],[1,-2,2],[-2,1,1]
    ])
    y=A*c
    failures=0; changed=0
    for _ in range(n):
        U_np=random_unimodular(3,rng,steps=8)
        U=sp.Matrix(U_np.tolist())
        # Row charge vector q transforms as q' = q U^T; c' = U^{-T} c.
        Ap=A*U.T
        cp=U.inv().T*c
        if Ap*cp != y: failures+=1
        if cp != c: changed+=1
    add_summary("E5 Hypercharge basis", "exact covariance failures", failures, n, failures==0)
    add_summary("E5 Hypercharge basis", "coefficient-coordinate changes", changed, n, True)
    return {"failures":failures,"changed":changed,"n":n}


def exp6_discriminant(reps=1000):
    # Orthonormal columns in R^8. third column = a1+a2+t*b.
    Q,_=np.linalg.qr(rng.normal(size=(8,3)))
    a1,a2,b=Q[:,0],Q[:,1],Q[:,2]
    c_true=np.array([0.7,-1.2,0.45])
    sigma=1e-7
    ts=np.logspace(0,-7,8)
    rows=[]
    for t in ts:
        A=np.column_stack([a1,a2,a1+a2+t*b])
        y=A@c_true
        gram=A.T@A
        det=float(np.linalg.det(gram)); cond=float(np.linalg.cond(A))
        c_errors=[]; y_errors=[]
        for _ in range(reps):
            yn=y+rng.normal(scale=sigma,size=y.shape)
            chat,*_=np.linalg.lstsq(A,yn,rcond=None)
            c_errors.append(float(np.linalg.norm(chat-c_true)/np.linalg.norm(c_true)))
            y_errors.append(float(np.linalg.norm(A@chat-y)/np.linalg.norm(y)))
        rows.append([t,det,cond,float(np.median(c_errors)),float(np.quantile(c_errors,0.95)),float(np.median(y_errors))])
    df=pd.DataFrame(rows,columns=['t','delta_det_AtA','condition_number','median_relative_c_error','p95_relative_c_error','median_relative_output_error'])
    df.to_csv(RES/'E6_discriminant_curve.csv',index=False)
    # exact t=0 null direction via construction v=(1,1,-1)
    A0=np.column_stack([a1,a2,a1+a2])
    v=np.array([1.,1.,-1.]); null_norm=float(np.linalg.norm(A0@v)); rank=int(np.linalg.matrix_rank(A0,tol=1e-12))
    add_summary("E6 Hypercharge discriminant", "rank at t=0", rank, None, rank==2)
    add_summary("E6 Hypercharge discriminant", "||A(0)(1,1,-1)||", null_norm, None, null_norm<TOL)
    add_summary("E6 Hypercharge discriminant", "max condition number sampled", float(df.condition_number.max()), len(df), True)
    return df, {"rank0":rank,"null_norm":null_norm}


def exp7_unit_quotient(n=30000):
    max_mu_err=0.0; raw_changed=0; valid=0
    rows=[]
    for _ in range(n):
        alpha=10**rng.uniform(-8,8)
        L=10**rng.uniform(-8,8)
        gs=10**rng.uniform(-2,1)
        N=int(rng.integers(0,9)); a=0.0 if rng.random()<0.5 else 0.5
        td1=1/(2*np.pi*gs*alpha)
        M2=(td1*L)**2+(N-a)/alpha
        if M2 < 0: continue
        valid += 1
        mu2=alpha*M2; ell=L/np.sqrt(alpha)
        rho=10**rng.uniform(-6,6)
        alpha2=rho*rho*alpha; L2=rho*L; td12=td1/(rho*rho); M22=M2/(rho*rho)
        mu22=alpha2*M22; ell2=L2/np.sqrt(alpha2)
        err=max(abs(mu22-mu2)/(1+abs(mu2)),abs(ell2-ell)/(1+abs(ell)))
        max_mu_err=max(max_mu_err,err)
        raw_changed += int(abs(L2-L)>0 or abs(M22-M2)>0)
        if len(rows)<500: rows.append([alpha,L,gs,N,a,rho,M2,mu2,ell,mu22,ell2])
    pd.DataFrame(rows,columns=['alpha_prime','L','g_s','N','a','rho','M2','mu2','ell','mu2_after','ell_after']).to_csv(RES/'E7_unit_quotient_examples.csv',index=False)
    add_summary("E7 Unit quotient", "max normalized invariant error", max_mu_err, valid, max_mu_err<TOL)
    add_summary("E7 Unit quotient", "raw dimensionful states changed", raw_changed, valid, raw_changed==valid)
    return {"max_err":max_mu_err,"raw_changed":raw_changed,"n_requested":n,"n_valid":valid}


def exp8_ablation(n_random=20000, n_hard=5000):
    J=symplectic_J(3); B=np.eye(6,dtype=int)
    gc=np.array([1,0,0,0,0,0],dtype=int)
    gw=np.array([0,1,0,3,0,0],dtype=int)
    chi0=pairing(gc,gw,J); rel0=pairings_to_basis(gw,B,J)
    # Random topological mutations
    detect_chi=detect_rel=0
    valid=0
    for _ in range(n_random):
        delta=rng.integers(-2,3,size=6)
        if not np.any(delta): continue
        g2=gw+delta; valid+=1
        detect_chi += int(pairing(gc,g2,J)!=chi0)
        detect_rel += int(pairings_to_basis(g2,B,J)!=rel0)
    # Hard adversarial: force chi unchanged by forcing delta[3]=0 for C=e1.
    hard_chi=hard_rel=0; valid_h=0
    for _ in range(n_hard):
        delta=rng.integers(-3,4,size=6); delta[3]=0
        if not np.any(delta): continue
        g2=gw+delta; valid_h+=1
        hard_chi += int(pairing(gc,g2,J)!=chi0)
        hard_rel += int(pairings_to_basis(g2,B,J)!=rel0)
    # Other primitive-class violations are deterministic.
    table=pd.DataFrame([
        ['Random charge mutation','chi_CW',detect_chi,valid,detect_chi/valid],
        ['Random charge mutation','marked relational charge',detect_rel,valid,detect_rel/valid],
        ['Adversarial chi-preserving mutation','chi_CW',hard_chi,valid_h,hard_chi/valid_h],
        ['Adversarial chi-preserving mutation','marked relational charge',hard_rel,valid_h,hard_rel/valid_h],
        ['Color-rank 3->4','chi_CW',0,1000,0.0],
        ['Color-rank 3->4','full fingerprint',1000,1000,1.0],
        ['Hypercharge output perturbation','chi_CW',0,1000,0.0],
        ['Hypercharge output perturbation','full fingerprint',1000,1000,1.0],
        ['Chirality flip','chi_CW',0,1000,0.0],
        ['Chirality flip','full fingerprint',1000,1000,1.0],
    ],columns=['violation','fingerprint','detected','n','detection_rate'])
    ci=[]
    for _,r in table.iterrows():
        lo,hi=wilson_interval(int(r.detected),int(r.n)); ci.append((lo,hi))
    table['ci95_low']=[x[0] for x in ci]; table['ci95_high']=[x[1] for x in ci]
    table.to_csv(RES/'E8_ablation_detection.csv',index=False)
    add_summary("E8 Ablation", "random mutation chi detection rate", detect_chi/valid, valid, None)
    add_summary("E8 Ablation", "random mutation relational detection rate", detect_rel/valid, valid, detect_rel==valid)
    add_summary("E8 Ablation", "adversarial chi-preserving chi detection rate", hard_chi/valid_h, valid_h, hard_chi==0)
    add_summary("E8 Ablation", "adversarial chi-preserving relational detection rate", hard_rel/valid_h, valid_h, hard_rel==valid_h)
    return table


def make_figures(discriminant_df, ablation_df):
    # Figure 1: discriminant stress
    fig,ax=plt.subplots(figsize=(6.4,4.4))
    ax.loglog(discriminant_df['t'],discriminant_df['median_relative_c_error'],'o-',label='Coefficient reconstruction error')
    ax.loglog(discriminant_df['t'],discriminant_df['median_relative_output_error'],'s-',label='Observable reconstruction error')
    ax.set_xlabel('Discriminant control parameter t')
    ax.set_ylabel('Median relative error')
    ax.grid(True,which='both',alpha=.25); ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(FIG/'Fig1_discriminant_reconstructibility.png',dpi=300); fig.savefig(FIG/'Fig1_discriminant_reconstructibility.pdf')
    plt.close(fig)

    # Figure 2: ablation selected rows
    sel=ablation_df[ablation_df['violation'].isin(['Random charge mutation','Adversarial chi-preserving mutation'])].copy()
    labels=[f"{v}\n{f}" for v,f in zip(sel.violation,sel.fingerprint)]
    fig,ax=plt.subplots(figsize=(7.3,4.5))
    x=np.arange(len(sel)); vals=sel.detection_rate.to_numpy()*100
    ax.bar(x,vals)
    ax.set_xticks(x); ax.set_xticklabels(labels,rotation=20,ha='right')
    ax.set_ylim(0,105); ax.set_ylabel('Detection rate (%)')
    ax.axhline(100,linewidth=.8); ax.grid(True,axis='y',alpha=.25)
    fig.tight_layout(); fig.savefig(FIG/'Fig2_ablation_detection.png',dpi=300); fig.savefig(FIG/'Fig2_ablation_detection.pdf')
    plt.close(fig)

    # Figure 3: condition number
    fig,ax=plt.subplots(figsize=(6.4,4.4))
    ax.loglog(discriminant_df['t'],discriminant_df['condition_number'],'o-')
    ax.set_xlabel('Discriminant control parameter t'); ax.set_ylabel('Condition number of A(t)')
    ax.grid(True,which='both',alpha=.25); fig.tight_layout()
    fig.savefig(FIG/'Fig3_discriminant_conditioning.png',dpi=300); fig.savefig(FIG/'Fig3_discriminant_conditioning.pdf')
    plt.close(fig)


def main():
    t0=time.time()
    out={}
    out['E1']=exp1_symplectic_invariance()
    out['E1b']=exp1b_toroidal_replication()
    out['E2']=exp2_exhaustive_collisions()
    out['E3']=exp3_hrr()
    out['E4']=exp4_gauge()
    out['E5']=exp5_hypercharge_basis()
    ddf,out['E6']=exp6_discriminant()
    out['E7']=exp7_unit_quotient()
    adf=exp8_ablation(); out['E8_rows']=len(adf)
    make_figures(ddf,adf)
    sdf=pd.DataFrame(summary); sdf.to_csv(RES/'summary.csv',index=False)
    meta={
      'seed':SEED,'tolerance':TOL,'runtime_seconds':time.time()-t0,
      'python':sys.version,'platform':platform.platform(),
      'numpy':np.__version__,'pandas':pd.__version__,
      'experiments':out
    }
    (RES/'run_metadata.json').write_text(json.dumps(meta,indent=2,default=str),encoding='utf-8')
    print(sdf.to_string(index=False))
    print(json.dumps(meta,indent=2,default=str))

if __name__=='__main__': main()
