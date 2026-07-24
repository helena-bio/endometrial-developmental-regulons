#!/usr/bin/env python3
"""Frozen TASK B Phase 2 execution. ASCII-only; no adaptive repair or tuning."""
import argparse, csv, hashlib, json, math, os, platform, sys, tempfile
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

TARGETS = ["GATA2", "SOX9", "HOXA9", "WT1", "PAX8", "LHX1"]
SUBTYPES = ["POLE", "MMRd", "NSMP", "p53abn"]
MASTER_SEED = 20260713
EXPECTED = {
"spec":"14343c81ed4a9211e26117d5489ef24d3d158a028121b4c78e633e956a24ef35",
"engine":"4d7bfb997c9dc188daf9e584a791ff490dfcfdd7fa365c20bca8a2c066f42065",
"design":"682099a03923778d170ff4bd46e8b475405419fe1fadb3259ce4c2e75b43f2c1",
"scores":"956b7a968f35d1da57719e3bc7392decbdd6e292ac31d4dec100b6632518b0cb",
"covariates":"cc46495c0abc6877182c9bff8d4ba44f11599d8bc5bc7a8cb6f8331aa1ec3c9a",
"order":"e494200b1a9d1c026679d039181bf34cfa3be62892e3d55bca2aeec0dab5e432",
"rules":"aa1eaea70991d9e915f499243a11964e01b7242c2f190cac783661b896e1eb42",
"input_lock":"bdaa15a40800cb5384abda946ef2c388f97fbb4fdcc16c4c901d9d00452d9039",
"decision":"a9b01e817d9e84180dc4804811dd2d0e33c245925a11bc8d5996363c42e7bc94",
"harm_spec":"c4252fd1d17a07212852a894631afd7a8df3b331ccf9aed314869b0e5512a362",
"harm_map":"c255355de4b6c4b366a626d7e79b26f474d55ead67387d1fde519743ac7e316b",
"ledger":"c6087dc8e685b305ccd3ff9d1da02846ff8a306018cc1fb0820cf3825fc8e823",
"missing":"02e08cb20fe60056b6a23b4bd9837d07ba80a6c729f1fcf366e28d50d97dc69e",
"raw_counts":"09fe226376379a333715c3dcfce4d1282e882a368e92c3eb3f6a9ece9abbec80",
"proposed":"91726649b4e910d5f08318304f052a044fba786e581d0efca0fee39b4737e0dc",
"critic":"1b83bae3d188db5e844c1786f3d5773c0e2c183d8a050d5c0f575b8c6eeadbcf"}

def sha(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(1<<20),b""): h.update(b)
    return h.hexdigest()
def subseed(step): return int.from_bytes(hashlib.sha256(f"{MASTER_SEED}|TASKB_PHASE2|{step}".encode("ascii")).digest()[:8],"big")
def fnum(x):
    if x is None or not np.isfinite(x): return "NA"
    return format(float(x),".17g")
def write_tsv(path, rows, fields):
    with open(path,"w",newline="",encoding="ascii") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n",extrasaction="ignore"); w.writeheader()
        for r in rows: w.writerow({k:(fnum(v) if isinstance(v,(float,np.floating)) else v) for k,v in r.items()})
def json_ready(obj):
    """Normalize JSON containers/scalars without changing numeric values."""
    if isinstance(obj, dict):
        return {str(k):json_ready(v) for k,v in obj.items()}
    if isinstance(obj, (list,tuple)):
        return [json_ready(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return json_ready(obj.tolist())
    if isinstance(obj, np.generic):
        return json_ready(obj.item())
    if obj is None or isinstance(obj,(str,int,float,bool)):
        return obj
    raise TypeError(f"Unsupported JSON type: {type(obj).__module__}.{type(obj).__qualname__}")
def dump_json(path,obj):
    with open(path,"w",encoding="ascii",newline="\n") as f: json.dump(json_ready(obj),f,sort_keys=True,indent=2,allow_nan=False); f.write("\n")

def design(subt, cov, clinical=None):
    cols=[np.ones(len(subt))]+[(subt==s).astype(float) for s in SUBTYPES[1:]]
    names=["intercept","subtype_MMRd","subtype_NSMP","subtype_p53abn"]
    for k,v in cov.items(): cols.append(np.asarray(v,float)); names.append(k)
    if clinical is not None: cols.append(np.asarray(clinical,float)); names.append("clinical_binary")
    return np.column_stack(cols),names

def fit(y,X):
    n,p=X.shape; b,_,rank,s=np.linalg.lstsq(X,y,rcond=None)
    tol=np.finfo(float).eps*max(n,p)*s[0]; rank_exact=int(np.sum(s>tol))
    if rank_exact<p: raise np.linalg.LinAlgError("rank deficient")
    e=y-X@b; df=n-p; sse=float(e@e); sig=math.sqrt(sse/df); inv=np.linalg.inv(X.T@X)
    L=np.zeros(p); L[1]=.5; L[2]=-1
    beta=float(L@b); se=math.sqrt(float(sig*sig*(L@inv@L))); t=beta/se; pval=float(2*stats.t.sf(abs(t),df)); crit=float(stats.t.ppf(.975,df))
    return {"coef":b,"resid":e,"sigma":sig,"df":df,"inv":inv,"beta":beta,"se":se,"t":t,"p":pval,"ci_lo":beta-crit*se,"ci_hi":beta+crit*se,"d":beta/sig,"rank":rank_exact,"tol":tol,"singular":s}

def diagnostics(y,X,names,subt,clinical,node,base_d):
    n,p=X.shape; out={}; reasons=[]
    try: z=(X[:,1:]-X[:,1:].mean(0))/X[:,1:].std(0,ddof=0); sv=np.linalg.svd(z,compute_uv=False); cond=float(sv[0]/sv[-1])
    except Exception: cond=float("inf")
    R=np.corrcoef(X[:,1:],rowvar=False); invR=np.linalg.pinv(R); vifs=np.diag(invR); maxv=float(np.max(vifs))
    # subtype adjusted GVIF; determinant formula on predictor correlation matrix
    A=[0,1,2]; B=list(range(3,R.shape[0])); detR=float(np.linalg.det(R)); gv=float("inf") if detR<=0 else float((np.linalg.det(R[np.ix_(A,A)])*np.linalg.det(R[np.ix_(B,B)])/detR)**(1/6)) if B else 1.0
    H=np.sum(X*(X@np.linalg.inv(X.T@X)),axis=1); frac=float(np.mean(H>3*p/n)); maxhat=float(H.max())
    e=base_d["resid"]; sig=base_d["sigma"]; df=base_d["df"]; h=np.clip(H,0,1-1e-14)
    stud=e/(sig*np.sqrt(1-h)); ext=stud*np.sqrt((df-1)/np.maximum(1e-15,df-stud*stud)); cook=(e*e/(p*sig*sig))*h/(1-h)**2
    dfb=(np.linalg.inv(X.T@X)@X.T).T*e[:,None]/(sig*np.sqrt(np.diag(np.linalg.inv(X.T@X)))[None,:]*(1-h)[:,None])
    flagged=np.where((cook>4/df)|(np.max(np.abs(dfb),axis=1)>2/math.sqrt(n))|(np.abs(ext)>3))[0]
    max_delta=0.; direction_change=False
    for i in flagged:
        mask=np.arange(n)!=i
        try: loo=fit(y[mask],X[mask]); max_delta=max(max_delta,abs(loo["d"]-base_d["d"])); direction_change |= np.sign(loo["d"])!=np.sign(base_d["d"])
        except Exception: max_delta=float("inf"); direction_change=True
    if n<max(80,10*p): reasons.append("support_n")
    if df<30: reasons.append("residual_df")
    if set(subt)!=set(SUBTYPES): reasons.append("missing_subtype")
    counts={s:int(np.sum(subt==s)) for s in SUBTYPES}
    if node=="histology_matched" and min(counts.values())<5: reasons.append("subtype_support")
    if node=="endometrioid_grade_matched" and min(counts.values())<10: reasons.append("subtype_support")
    cram=None; dummycorr=None; clinical_counts={}
    if clinical is not None:
        clinical_counts={str(k):int(np.sum(clinical==k)) for k in [0.,1.]}
        if min(clinical_counts.values())<20: reasons.append("clinical_overall_support")
        tab=np.array([[np.sum((subt==s)&(clinical==k)) for k in [0.,1.]] for s in SUBTYPES])
        if np.sum(np.sum(tab>=5,axis=0)>=2)<2: reasons.append("clinical_cross_subtype_support")
        for j in range(2):
            if tab[:,j].sum() and tab[:,j].max()/tab[:,j].sum()>=.9: reasons.append("clinical_level_concentration")
        chi=stats.chi2_contingency(tab,correction=False)[0]; cram=float(math.sqrt(chi/(n*min(tab.shape[0]-1,tab.shape[1]-1))))
        dummycorr=float(np.nanmax(np.abs(R[:3,-1])))
        fr=max(clinical_counts.values())/max(1,min(clinical_counts.values())); pu=200/n
        if fr>19 and pu<=10: reasons.append("near_zero_variance")
        if cram>.9: reasons.append("cramers_v")
        if dummycorr>=.95: reasons.append("dummy_correlation")
    if base_d["rank"]<p: reasons.append("rank")
    if cond>100: reasons.append("condition_number_not_estimable")
    elif cond>30: reasons.append("condition_number_unstable")
    if maxv>10 or gv>10: reasons.append("vif_fail")
    elif maxv>5 or gv>5: reasons.append("vif_unstable")
    if maxhat>=.5 or frac>.05: reasons.append("leverage")
    if direction_change: reasons.append("single_case_direction_change")
    if max_delta>=.2: reasons.append("single_case_delta_d")
    status="PASS" if not reasons else ("NOT_ESTIMABLE" if any(x in reasons for x in ["rank","condition_number_not_estimable","vif_fail"]) else "UNSTABLE_NOT_INTERPRETED")
    out.update(n=n,p=p,rank=base_d["rank"],residual_df=df,condition_number=cond,max_vif=maxv,subtype_gvif_adjusted=gv,max_hat=maxhat,fraction_hat_gt_3p_n=frac,n_hat_gt_2p_n=int(np.sum(H>2*p/n)),n_cooks_flags=int(np.sum(cook>4/df)),n_dfbeta_flags=int(np.sum(np.max(np.abs(dfb),axis=1)>2/math.sqrt(n))),n_studentized_flags=int(np.sum(np.abs(ext)>3)),n_influence_union=int(len(flagged)),max_single_case_abs_d_change=max_delta,single_case_direction_change=direction_change,cramers_v=cram,max_clinical_dummy_abs_correlation=dummycorr,status=status,reasons=";".join(sorted(set(reasons))) or "none",subtype_counts=json.dumps(counts,sort_keys=True,separators=(",",":")),clinical_counts=json.dumps(clinical_counts,sort_keys=True,separators=(",",":")))
    infl=[{"row_index":int(i),"hat":H[i],"cooks_d":cook[i],"max_abs_dfbeta":float(np.max(np.abs(dfb[i]))),"externally_studentized":ext[i]} for i in flagged]
    return out,infl

def boot_single(y,subt,cov,clinical,seed,B):
    rng=np.random.default_rng(seed); vals=[]; n=len(y)
    for _ in range(B):
        ix=rng.integers(0,n,n)
        if len(set(subt[ix]))<4: continue
        try: q=fit(y[ix],design(subt[ix],{k:v[ix] for k,v in cov.items()},None if clinical is None else clinical[ix])[0]); vals.append((q["beta"],q["d"]))
        except (np.linalg.LinAlgError,FloatingPointError): pass
    a=np.asarray(vals,float)
    return (len(vals),np.percentile(a[:,0],[2.5,97.5]).tolist(),np.percentile(a[:,1],[2.5,97.5]).tolist()) if len(vals) else (0,[None,None],[None,None])

def perm_p(y,subt,cov,clinical,obs,seed,P):
    rng=np.random.default_rng(seed); ge=0
    for _ in range(P):
        sp=rng.permutation(subt)
        try: q=fit(y,design(sp,cov,clinical)[0]); ge += abs(q["d"])>=abs(obs)-1e-12
        except np.linalg.LinAlgError: pass
    return (1+ge)/(P+1)

def paired_boot(y,subt,cov,clinical,seed,B):
    rng=np.random.default_rng(seed); vals=[]; n=len(y)
    for _ in range(B):
        ix=rng.integers(0,n,n)
        if len(set(subt[ix]))<4: continue
        try:
            cb={k:v[ix] for k,v in cov.items()}; b=fit(y[ix],design(subt[ix],cb)[0]); a=fit(y[ix],design(subt[ix],cb,clinical[ix])[0])
            sd=a["d"]-b["d"]; ma=abs(b["d"])-abs(a["d"]); rs=b["beta"]*(1/a["sigma"]-1/b["sigma"]); cc=(a["beta"]-b["beta"])/a["sigma"]
            vals.append([sd,abs(sd),ma,(100*ma/abs(b["d"]) if abs(b["d"])>=.05 and np.sign(a["d"])==np.sign(b["d"]) else np.nan),a["beta"]-b["beta"],rs,cc])
        except (np.linalg.LinAlgError,FloatingPointError): pass
    a=np.asarray(vals,float); return len(vals),np.nanpercentile(a,[2.5,97.5],axis=0)

def taxonomy(b,a,unstable):
    if unstable:return "unstable"
    if abs(b)<.05:return "base_near_null_unclassifiable"
    if np.sign(b)!=np.sign(a):return "sign_reversal"
    if abs(a)>abs(b)+1e-12:return "amplified"
    if abs(a)<.05:return "near_null_or_disappears_descriptively"
    att=abs(b)-abs(a); pct=100*att/abs(b)
    return "same_direction_materially_attenuated" if att>=.10 and pct>=20 else "largely_retained"

def synthetic_test():
    rng=np.random.default_rng(17); n=120; order=np.array([f"P{i:03d}" for i in range(n)],object); subt=np.array(SUBTYPES*(n//4),object)
    cov={"M4_prolif":rng.normal(size=n),"purity_CPE":rng.normal(size=n),"composition":rng.normal(size=n)}; y=rng.normal(size=n)+.2*(subt=="MMRd")-.3*(subt=="NSMP")
    with tempfile.TemporaryDirectory(prefix="taskb_cycle6_synthetic_") as td:
        td=Path(td); order_path=td/"synthetic_order.json"; score_path=td/"synthetic_scores.npz"; json_path=td/"synthetic_fixture.json"
        with open(order_path,"w") as f: json.dump(order.tolist(),f)
        np.savez(score_path,patient_order=order,**{f"M3primary__{t}":y for t in TARGETS})
        z=np.load(score_path,allow_pickle=True); assert z["patient_order"].shape==(n,) and all(z[f"M3primary__{t}"].shape==(n,) for t in TARGETS); z.close()
        X,names=design(subt,cov); q=fit(y,X); assert np.isfinite(q["d"])
        dg,_=diagnostics(y,X,names,subt,None,"full_base",q)
        b=boot_single(y,subt,cov,None,subseed("TCGA-UCEC|primary_cpe|full_base|GATA2|C2|bootstrap"),20); assert b[0]>=19
        p=perm_p(y,subt,cov,None,q["d"],subseed("TCGA-UCEC|primary_cpe|full_base|GATA2|C2|permutation"),20); assert 0<p<=1
        fixture={
            "python_scalars":[None,True,False,0,1,-1,1.25,"ASCII"],
            "numpy_scalars":[np.bool_(True),np.bool_(False),np.int8(-2),np.int64(7),np.uint64(9),np.float32(1.5),np.float64(-2.25),np.str_("text")],
            "numpy_arrays":[np.array(True),np.array([True,False],dtype=np.bool_),np.array([1,2],dtype=np.int64),np.array([1.25,-3.5],dtype=np.float64)],
            "nested":{"tuple":(np.bool_(True),{"integer":np.int32(4),"array":np.array([[1.0],[2.0]])}),"diagnostics":dg},
        }
        dump_json(json_path,fixture)
        loaded=json.load(open(json_path,encoding="ascii"))
        assert loaded["numpy_scalars"]==[True,False,-2,7,9,1.5,-2.25,"text"]
        assert loaded["numpy_arrays"]==[True,[True,False],[1,2],[1.25,-3.5]]
        assert isinstance(loaded["nested"]["diagnostics"]["single_case_direction_change"],bool)
        try: json_ready(complex(1,2))
        except TypeError: pass
        else: raise AssertionError("unsupported complex scalar was silently serialized")

def execute(args):
    out=Path(args.output); out.mkdir(parents=True,exist_ok=False)
    paths=json.load(open(args.config,encoding="ascii"))["paths"]
    audit=[]
    for k,p in paths.items():
        got=sha(p); exp=EXPECTED.get(k); audit.append({"key":k,"path":p,"sha256":got,"expected_sha256":exp,"status":"MATCH" if got==exp else "MISMATCH"})
        if exp and got!=exp: raise RuntimeError(f"hash mismatch: {k}")
    # Object loading only after hash verification, only pinned string metadata.
    order=json.load(open(paths["order"],encoding="ascii")); cov=pd.read_csv(paths["covariates"],sep="\t",dtype={"patient_barcode":str})
    if cov.patient_barcode.tolist()!=order: raise RuntimeError("covariate order mismatch")
    z=np.load(paths["scores"],allow_pickle=True); npz_order=z["patient_order"].tolist()
    if npz_order!=order: raise RuntimeError("NPZ patient_order mismatch")
    outcomes={t:np.asarray(z[f"M3primary__{t}"],dtype=np.float64) for t in TARGETS}; z.close()
    led=pd.read_csv(paths["ledger"],sep="\t",dtype=str,keep_default_na=False)
    tc=led[led["cohort"]=="TCGA"].copy(); idcol="analytic_patient_id" if "analytic_patient_id" in tc else "patient_id"
    tc=tc.drop_duplicates(idcol).set_index(idcol); tc=tc.reindex(order)
    if tc.index.hasnans or len(tc)!=len(order): raise RuntimeError("ledger alignment failed")
    hcol="histology_proposed"; gcol="grade_proposed"
    hist=tc[hcol].to_numpy(str); grade=tc[gcol].to_numpy(str)
    subtype=cov.subtype.to_numpy(str); assert set(subtype)==set(SUBTYPES)
    base_specs={"primary_cpe":(["M4_prolif","purity_CPE","composition"],cov.cpe_complete_case.astype(str).str.lower().isin(["true","1"]).to_numpy()),"no_purity":(["M4_prolif","composition"],np.ones(len(cov),bool))}
    results=[]; diags=[]; influences=[]; decomp=[]; fitted={}
    for bs,(cnames,roster) in base_specs.items():
      for node in ["full_base","histology_matched","endometrioid_grade_matched"]:
        mask=roster.copy(); clinical_all=None; variants=["base"]
        if node=="histology_matched": mask &= np.isin(hist,["endometrioid","non_endometrioid"]); clinical_all=(hist=="non_endometrioid").astype(float); variants=["base_refit","base_plus_binary_histology"]
        if node=="endometrioid_grade_matched": mask &= (hist=="endometrioid") & np.isin(grade,["low_grade","high_grade"]); clinical_all=(grade=="low_grade").astype(float); variants=["base_refit","base_plus_binary_grade"]
        idx=np.where(mask)[0]; ss=subtype[idx]; cv={c:cov[c].to_numpy(float)[idx] for c in cnames}; cl=None if clinical_all is None else clinical_all[idx]
        for target in TARGETS:
          yy=outcomes[target][idx]
          for variant in variants:
            usecl=cl if variant.startswith("base_plus") else None; X,names=design(ss,cv,usecl); q=fit(yy,X)
            step="|".join(["TCGA-UCEC",bs,node,target,"C2"])
            nb,ci,di=boot_single(yy,ss,cv,usecl,subseed(step+"|bootstrap|"+variant),args.bootstrap)
            pp=perm_p(yy,ss,cv,usecl,q["d"],subseed(step+"|permutation|"+variant),args.permutations)
            dg,inf=diagnostics(yy,X,names,ss,usecl,node,q); unstable=dg["status"]!="PASS" or nb<1900
            if nb<1900: dg["status"]="UNSTABLE_NOT_INTERPRETED"; dg["reasons"]=(dg["reasons"]+";bootstrap_valid_lt_1900").strip(";")
            mid=f"{bs}|{node}|{variant}|{target}"; fitted[mid]=(q,dg,nb,ci,di)
            row={"cohort":"TCGA-UCEC","base_specification":bs,"node":node,"model":variant,"target":target,"contrast":"C2","n":len(idx),"p_design":X.shape[1],"coefficient":q["beta"],"se":q["se"],"t":q["t"],"residual_df":q["df"],"ci_lo":q["ci_lo"],"ci_hi":q["ci_hi"],"raw_p":q["p"],"direction":"positive" if q["beta"]>0 else "negative" if q["beta"]<0 else "zero","residual_sd":q["sigma"],"d":q["d"],"d_ci_lo":di[0],"d_ci_hi":di[1],"coefficient_boot_ci_lo":ci[0],"coefficient_boot_ci_hi":ci[1],"bootstrap_attempts":args.bootstrap,"bootstrap_valid":nb,"diagnostic_permutations":args.permutations,"permutation_p":pp,"interpretability":dg["status"]}; results.append(row)
            diags.append({"model_id":mid,**dg}); influences.extend({"model_id":mid,"patient_barcode":order[idx[x["row_index"]]],**x} for x in inf)
        if node!="full_base":
          for target in TARGETS:
            b,dbg,*_=fitted[f"{bs}|{node}|base_refit|{target}"]; a,dga,*_=fitted[f"{bs}|{node}|{variants[1]}|{target}"]
            sd=a["d"]-b["d"]; ma=abs(b["d"])-abs(a["d"]); rs=b["beta"]*(1/a["sigma"]-1/b["sigma"]); cc=(a["beta"]-b["beta"])/a["sigma"]
            if abs((rs+cc)-sd)>1e-12: raise RuntimeError("decomposition tolerance failure")
            step="|".join(["TCGA-UCEC",bs,node,target,"C2","paired_bootstrap"]); nv,pci=paired_boot(outcomes[target][idx],ss,cv,cl,subseed(step),args.bootstrap)
            unstable=dbg["status"]!="PASS" or dga["status"]!="PASS" or nv<1900
            pct=None if abs(b["d"])<.05 or np.sign(b["d"])!=np.sign(a["d"]) else 100*ma/abs(b["d"])
            decomp.append({"cohort":"TCGA-UCEC","base_specification":bs,"node":node,"target":target,"n":len(idx),"base_beta":b["beta"],"adjusted_beta":a["beta"],"base_residual_sd":b["sigma"],"adjusted_residual_sd":a["sigma"],"base_d":b["d"],"adjusted_d":a["d"],"signed_delta_d":sd,"absolute_delta_d":abs(sd),"magnitude_attenuation":ma,"percent_attenuation":pct,"beta_change":a["beta"]-b["beta"],"residual_scale_contribution":rs,"coefficient_contribution":cc,"decomposition_error":abs(rs+cc-sd),"paired_bootstrap_attempts":args.bootstrap,"paired_bootstrap_valid":nv,"signed_delta_d_ci_lo":pci[0,0],"signed_delta_d_ci_hi":pci[1,0],"taxonomy":taxonomy(b["d"],a["d"],unstable),"interpretability":"UNSTABLE_NOT_INTERPRETED" if unstable else "PASS"})
    rf=list(results[0]); write_tsv(out/"MODEL_RESULTS.tsv",results,rf); write_tsv(out/"SIX_TARGET_APPENDIX.tsv",results,rf)
    write_tsv(out/"MODEL_DIAGNOSTICS.tsv",diags,list(diags[0])); write_tsv(out/"INFLUENCE_RECORDS.tsv",influences,list(influences[0]) if influences else ["model_id","patient_barcode","row_index","hat","cooks_d","max_abs_dfbeta","externally_studentized"])
    write_tsv(out/"MATCHED_DECOMPOSITIONS.tsv",decomp,list(decomp[0]))
    counts=[]
    for bs,(_,rost) in base_specs.items():
      for node,mask in [("full_base",rost),("histology_matched",rost&np.isin(hist,["endometrioid","non_endometrioid"])),("endometrioid_grade_matched",rost&(hist=="endometrioid")&np.isin(grade,["low_grade","high_grade"]))]:
        for s in SUBTYPES: counts.append({"base_specification":bs,"node":node,"stratum":"subtype","level":s,"n":int(np.sum(mask&(subtype==s)))})
    write_tsv(out/"COHORT_COUNTS.tsv",counts,list(counts[0]));
    miss=[{"field":"histology_proposed","missing_n":int(np.sum(~np.isin(hist,["endometrioid","non_endometrioid"])))},{"field":"grade_proposed","missing_n":int(np.sum(~np.isin(grade,["low_grade","high_grade"])))},{"field":"purity_CPE","missing_n":int(np.sum(~np.isfinite(cov.purity_CPE.to_numpy(float))))}]
    write_tsv(out/"MISSINGNESS.tsv",miss,list(miss[0])); dist=[]
    for field,arr in [("histology_proposed",hist),("grade_proposed",grade)]:
      for val,n in sorted(zip(*np.unique(arr,return_counts=True))): dist.append({"field":field,"level":val or "MISSING","n":int(n)})
    write_tsv(out/"CLINICAL_DISTRIBUTIONS.tsv",dist,list(dist[0]));
    summary=[d for d in decomp if d["target"] in ["GATA2","SOX9"]]; write_tsv(out/"MANUSCRIPT_READY_SENSITIVITY_TABLE.tsv",summary,list(summary[0]))
    dump_json(out/"SCIENTIFIC_RESULTS.json",{"classification":"post_hoc_explanatory_sensitivity","models":results,"decompositions":decomp,"diagnostics":diags,"no_q_values":True,"no_verdict":True})
    text="# Proposed wording (proposal only; no manuscript edited)\n\n## Results\n\nIn the post-hoc TCGA-UCEC explanatory sensitivity, the frozen C2 contrasts were refit on exactly matched rows before binary histology adjustment and, within endometrioid carcinoma, binary FIGO grade adjustment. Results are descriptive and are provided in MANUSCRIPT_READY_SENSITIVITY_TABLE.tsv; raw p values and diagnostic permutation p values carry no confirmatory or multiplicity credit.\n\n## Methods\n\nWe used frozen TASK-028 signed regulon outcomes for six prespecified targets and OLS models containing the four-level molecular subtype factor, M4 proliferation, ESTIMATE-derived composition, and, in the primary model, CPE purity. C2 was 0.5*POLE + 0.5*MMRd - NSMP with p53abn retained at weight zero. Histology and endometrioid-only grade additions were treatment coded under the frozen Phase-1 map. We used 2,000 attempted patient bootstraps and 2,000 diagnostic subtype-label permutations with deterministic SHA-256-derived sub-seeds.\n\n## Limitations\n\nThis post-hoc analysis is not a unique causal estimand. Grade and histology may be confounders, mediators, consequences of subtype biology, pathology proxies, or composition markers. CPE and ESTIMATE are estimated proxies and do not establish purity independence or complete cell-type adjustment. TCGA-only internal bootstrap is not external validation. The analysis supports no individual-patient biomarker, treatment response, target category, q value, verdict revision, or manuscript change.\n"
    (out/"PROPOSED_WORDING.md").write_text(text,encoding="ascii",newline="\n")
    report=f"# Analytical report\n\nProducer execution status: COMPLETE.\n\nThis package contains {len(results)} model rows, {len(decomp)} matched-row decompositions, {len(diags)} diagnostic rows, and {len(influences)} flagged influence records. It is a post-hoc explanatory TCGA-only sensitivity. No scientific verdict is made. Every model used {args.bootstrap} attempted bootstraps and {args.permutations} diagnostic permutations.\n"
    (out/"ANALYTICAL_REPORT.md").write_text(report,encoding="ascii",newline="\n")
    write_tsv(out/"INPUT_ACCESS_INVENTORY.tsv",audit,list(audit[0])); return {"models":len(results),"decompositions":len(decomp),"influence_records":len(influences),"audit_inputs":len(audit)}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--synthetic-test",action="store_true"); ap.add_argument("--config"); ap.add_argument("--output"); ap.add_argument("--bootstrap",type=int,default=2000); ap.add_argument("--permutations",type=int,default=2000); a=ap.parse_args()
    if a.synthetic_test: synthetic_test(); print("SYNTHETIC_TEST_PASS"); return
    print(json.dumps(execute(a),sort_keys=True))
if __name__=="__main__": main()
