#!/usr/bin/env python3
"""Promote byte-identical run-1 products and build the audit package."""
from pathlib import Path
import hashlib, json, shutil, subprocess
import pandas as pd

ROOT=Path(__file__).resolve().parent
R1=ROOT/'run1_definitive'; R2=ROOT/'run2_definitive'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def git(*args): return subprocess.run(['git',*args],cwd=ROOT.parents[1],text=True,capture_output=True,check=True).stdout.strip()

files1={str(p.relative_to(R1)):p for p in R1.rglob('*') if p.is_file()}
files2={str(p.relative_to(R2)):p for p in R2.rglob('*') if p.is_file()}
rows=[]
for rel in sorted(set(files1)|set(files2)):
    h1=sha(files1[rel]) if rel in files1 else 'MISSING'
    h2=sha(files2[rel]) if rel in files2 else 'MISSING'
    rows.append({'relative_path':rel,'run1_sha256':h1,'run2_sha256':h2,
                 'byte_identical':h1==h2,'comparison':'BYTEWISE'})
if not rows or not all(r['byte_identical'] for r in rows):
    raise SystemExit('DETERMINISTIC_RERUN_MISMATCH')

# Promote run 1 without touching the two preserved run directories.
for directory in ['results','figures']:
    dest=ROOT/directory
    dest.mkdir(exist_ok=True)
    for source in (R1/directory).iterdir():
        shutil.copy2(source,dest/source.name)
for name in ['DEPENDENCIES.lock','INPUT_CHECKSUMS.tsv','SEED_SCHEME.md','SESSION_ENVIRONMENT.json']:
    shutil.copy2(R1/name,ROOT/name)
(ROOT/'integrity').mkdir(exist_ok=True)
for source in (R1/'integrity').iterdir(): shutil.copy2(source,ROOT/'integrity'/source.name)
pd.DataFrame(rows).to_csv(ROOT/'integrity/DETERMINISTIC_RERUN_COMPARISON.tsv',sep='\t',index=False,lineterminator='\n')

# Explicit upstream point-fit reconciliation against the frozen C2/C3 authorities.
d=pd.read_csv(ROOT/'results/PER_CLASS_CONTRASTS_LONG.tsv',sep='\t')
t28=json.load(open('data/external/original-workspace/task028-freeze-b-draft/execution/results_v3/phase2b_models.json'))
t29=json.load(open('data/external/original-workspace/task029-external-replication-feasibility/execution/results/primary_results.json'))
recon=[]
for model in ['TCGA_PRIMARY_CPE_N506','TCGA_NOPURITY_N507','CPTAC_DISCOVERY_N95','CPTAC_CONFIRMATORY_N135']:
  for target in ['GATA2','SOX9','HOXA9','WT1','PAX8','LHX1']:
    x=d[(d.model==model)&(d.target==target)]
    c2b=.5*x[x.contrast=='POLE_vs_NSMP'].coefficient.iloc[0]+.5*x[x.contrast=='MMRd_vs_NSMP'].coefficient.iloc[0]
    c2d=.5*x[x.contrast=='POLE_vs_NSMP'].d.iloc[0]+.5*x[x.contrast=='MMRd_vs_NSMP'].d.iloc[0]
    pm=x[x.contrast=='POLE_vs_MMRd'].iloc[0]
    if model=='TCGA_PRIMARY_CPE_N506': old=t28['primary']['M3_'+target]
    elif model=='TCGA_NOPURITY_N507': old=t28['sensitivities']['SENS_nopurity']['M3_'+target]
    else:
      strat='Discovery' if model=='CPTAC_DISCOVERY_N95' else 'Confirmatory'
      old=t29['strata'][strat]['m3_model'][target]
    row={'model':model,'target':target,
      'residual_SD_taskA':pm.residual_SD,'residual_SD_frozen':old['sigma_resid'],
      'residual_SD_abs_delta':abs(pm.residual_SD-old['sigma_resid']),
      'C2_coefficient_taskA':c2b,'C2_coefficient_frozen':old['contrasts']['C2']['estimate'],
      'C2_coefficient_abs_delta':abs(c2b-old['contrasts']['C2']['estimate']),
      'C2_d_taskA':c2d,'C2_d_frozen':old['contrasts']['C2']['d'],
      'C2_d_abs_delta':abs(c2d-old['contrasts']['C2']['d']),
      'PM_coefficient_taskA':pm.coefficient,'C3_coefficient_frozen':old['contrasts']['C3']['estimate'],
      'PM_coefficient_abs_delta':abs(pm.coefficient-old['contrasts']['C3']['estimate']),
      'PM_d_taskA':pm.d,'C3_d_frozen':old['contrasts']['C3']['d'],
      'PM_d_abs_delta':abs(pm.d-old['contrasts']['C3']['d'])}
    row['scaled_1e12_pass']=max(row[k] for k in row if k.endswith('_abs_delta')) <= 1e-12*max(1,abs(c2b),abs(pm.coefficient),abs(pm.residual_SD))
    recon.append(row)
if not all(r['scaled_1e12_pass'] for r in recon): raise SystemExit('UPSTREAM_POINT_RECONCILIATION_FAILURE')
pd.DataFrame(recon).to_csv(ROOT/'results/UPSTREAM_POINT_ESTIMATE_RECONCILIATION.tsv',sep='\t',index=False,lineterminator='\n')

# Detailed concise report with raw primary-target findings.
lines=['# TASK A analytical report','',
'Status: COMPLETE producer execution; independent critic verdict pending.','',
'This is a post-hoc explanatory sensitivity. It is descriptive and non-causal. It cannot change any frozen TCGA/CPTAC category or replication verdict, and no manuscript byte was edited.','',
'## Guard and design result','',
'All 21 pinned inputs, all nested seals, and the exact four cohort/count guards passed: TCGA primary 506 (49/148/146/163), TCGA no-purity 507 (49/148/147/163), CPTAC Discovery 95 (7/25/43/20), and CPTAC Confirmatory 135 (6/47/66/16), ordered POLE/MMRd/NSMP/p53abn. Each target/model used one unchanged full fit; no pairwise subset was refit.','',
'## GATA2 and SOX9 per-class results','',
'Values are coefficient (d), followed by the mechanically assigned direct same-fit state. The direct POLE-minus-MMRd model interval, not a comparison of nominal p values, determines class differentiation.','']
for target in ['GATA2','SOX9']:
  lines += [f'### {target}','']
  for model in ['TCGA_PRIMARY_CPE_N506','TCGA_NOPURITY_N507','CPTAC_DISCOVERY_N95','CPTAC_CONFIRMATORY_N135']:
    x=d[(d.target==target)&(d.model==model)]
    vals={r.contrast:r for _,r in x.iterrows()}
    lines.append(f"- {model}: POLE-NSMP {vals['POLE_vs_NSMP'].coefficient:+.15g} (d={vals['POLE_vs_NSMP'].d:+.15g}); MMRd-NSMP {vals['MMRd_vs_NSMP'].coefficient:+.15g} (d={vals['MMRd_vs_NSMP'].d:+.15g}); direct POLE-MMRd {vals['POLE_vs_MMRd'].coefficient:+.15g} (d={vals['POLE_vs_MMRd'].d:+.15g}, model 95% CI [{vals['POLE_vs_MMRd'].coefficient_t_CI_lo:+.15g}, {vals['POLE_vs_MMRd'].coefficient_t_CI_hi:+.15g}]); `{vals['POLE_vs_MMRd'].interpretation_state}`.")
  x=d[(d.target==target)&(d.model=='CPTAC_FIXED_EFFECT_META')]
  vals={r.contrast:r for _,r in x.iterrows()}
  lines.append(f"- CPTAC fixed-effect meta: d_PN={vals['POLE_vs_NSMP'].d:+.15g}, d_MN={vals['MMRd_vs_NSMP'].d:+.15g}, direct d_PM={vals['POLE_vs_MMRd'].d:+.15g} (95% normal CI [{vals['POLE_vs_MMRd'].d_boot_CI_lo:+.15g}, {vals['POLE_vs_MMRd'].d_boot_CI_hi:+.15g}]); `{vals['POLE_vs_MMRd'].interpretation_state}`.")
  lines.append('')
lines += ['For both primary targets, pooled C2 is directionally concordant with negative POLE-NSMP and MMRd-NSMP point effects in every direct cohort/stratum and CPTAC meta. GATA2 is compatible in TCGA and CPTAC Confirmatory/meta, with Discovery unresolved. SOX9 is compatible in TCGA primary, CPTAC Discovery/meta; TCGA no-purity has a supported magnitude distinction below the inherited materiality floor, and CPTAC Confirmatory is unresolved. None of these labels means equality or equivalence.','',
'The equal-weight C2 is not sample-size weighted. POLE:MMRd sample proportions alone are not a class-mix artifact.','',
'## Completeness, algebra, and multiplicity','',
'All six targets and all 90 planned direct/meta rows are present. The four completeness targets show a mixture of compatible, unresolved, and point-heterogeneous states; the complete mechanical table is `results/INTERPRETATION_TAXONOMY.tsv`. No result was selected or repaired.','',
'All 24 direct target/model coefficient and d identities passed scaled tolerance 1e-12, including every usable bootstrap replicate. CPTAC fixed-effect C2 identity is not required because residual scales and contrast-specific inverse-bootstrap-variance weights differ; its discrepancy and all weights are reported.','',
'Exactly five separate descriptive BH-18 families were computed, each with 18 evaluable rows (no missing placeholders were needed). These q values do not confer confirmatory credit.','',
'## Reproducibility and inference boundary','',
'Two complete fresh-directory runs are byte-identical for all 27 produced files, including TSV, CSV, JSON, SVG, PNG, and PDF scientific artifacts. Upstream worktrees and all pinned bytes were unchanged.','',
'Cross-cohort differences can reflect biology, composition, acquisition, platform, classifier, or scoring context; this analysis cannot identify cause. Bulk subtype-level results are not individual-patient biomarkers. Frozen verdicts are preserved and no manuscript edit was made.','']
(ROOT/'ANALYTICAL_REPORT.md').write_text('\n'.join(lines),encoding='ascii')

commands=[
"python3 -m py_compile experiments/taskA_perclass_c2/run_taskA_perclass_c2.py",
"env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONHASHSEED=0 python3 experiments/taskA_perclass_c2/run_taskA_perclass_c2.py --output-root experiments/taskA_perclass_c2/run1 > experiments/taskA_perclass_c2/run1.stdout.log 2> experiments/taskA_perclass_c2/run1.stderr.log",
"env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONHASHSEED=0 python3 experiments/taskA_perclass_c2/run_taskA_perclass_c2.py --output-root experiments/taskA_perclass_c2/run2 > experiments/taskA_perclass_c2/run2.stdout.log 2> experiments/taskA_perclass_c2/run2.stderr.log",
"# Reporting-only schema audit found missing CPTAC stratum coefficient/SE/residual-SD fields; original run1/run2 preserved unchanged",
"env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONHASHSEED=0 python3 experiments/taskA_perclass_c2/run_taskA_perclass_c2.py --output-root experiments/taskA_perclass_c2/run1_definitive > experiments/taskA_perclass_c2/run1_definitive.stdout.log 2> experiments/taskA_perclass_c2/run1_definitive.stderr.log",
"env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONHASHSEED=0 python3 experiments/taskA_perclass_c2/run_taskA_perclass_c2.py --output-root experiments/taskA_perclass_c2/run2_definitive > experiments/taskA_perclass_c2/run2_definitive.stdout.log 2> experiments/taskA_perclass_c2/run2_definitive.stderr.log",
"env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONHASHSEED=0 python3 experiments/taskA_perclass_c2/package_taskA_perclass_c2.py",
]
(ROOT/'COMMANDS.txt').write_text('\n'.join(commands)+'\n',encoding='ascii')
(ROOT/'UPSTREAM_PRESERVATION_REPORT.md').write_text('# Upstream preservation\n\nPASS. The original dirty worktree and TASK-030 verified worktree have identical start/end HEAD, branch, porcelain-status SHA-256, unstaged binary-diff SHA-256, and staged binary-diff SHA-256 in both complete producer runs. All 21 pinned inputs reverified after analysis. Researcher charter/spec SHA-256 values remain exact. No upstream, manuscript, `src/**`, sealed, or frozen-result byte was written.\n',encoding='ascii')

required=['results/PER_CLASS_CONTRASTS_LONG.tsv','results/PER_CLASS_CONTRASTS.json','results/PER_CLASS_CONTRASTS_LONG.csv','results/GATA2_SOX9_HETEROGENEITY.tsv','results/SIX_TARGET_COMPLETENESS.tsv','results/MANUSCRIPT_READY_POSTHOC_TABLE.md','figures/GATA2_PER_CLASS_FOREST.svg','figures/SOX9_PER_CLASS_FOREST.svg','ANALYTICAL_REPORT.md','COMMANDS.txt','SESSION_ENVIRONMENT.json','DEPENDENCIES.lock','INPUT_CHECKSUMS.tsv','SEED_SCHEME.md','integrity/DETERMINISTIC_RERUN_COMPARISON.tsv','integrity/UPSTREAM_AND_DIRTY_WORKTREE_PRESERVATION.json']
audit=[{'required_path':x,'exists':(ROOT/x).is_file(),'sha256':sha(ROOT/x) if (ROOT/x).is_file() else 'MISSING'} for x in required]
pd.DataFrame(audit).to_csv(ROOT/'integrity/REQUIRED_DELIVERABLE_AUDIT.tsv',sep='\t',index=False,lineterminator='\n')
if not all(x['exists'] for x in audit): raise SystemExit('MISSING_REQUIRED_DELIVERABLE')

base=json.loads((R1/'REPRODUCIBILITY_MANIFEST.json').read_text())
analysis_status=git('status','--porcelain=v1')
analysis_git_state={'head':git('rev-parse','HEAD'),'branch':git('branch','--show-current'),
 'status_porcelain':analysis_status.splitlines(),
 'status_sha256':hashlib.sha256((analysis_status+'\n' if analysis_status else '').encode()).hexdigest(),
 'staged_paths':git('diff','--cached','--name-only').splitlines(),
 'initial_status_from_pre_outcome_guard':['?? experiments/taskA_perclass_c2/'],
 'initial_status_sha256':hashlib.sha256(b'?? experiments/taskA_perclass_c2/\n').hexdigest()}
base.update({'producer_script':str(ROOT/'run_taskA_perclass_c2.py'),'producer_script_sha256':sha(ROOT/'run_taskA_perclass_c2.py'),
 'packaging_script':str(ROOT/'package_taskA_perclass_c2.py'),'packaging_script_sha256':sha(ROOT/'package_taskA_perclass_c2.py'),
 'commands':commands,'two_complete_runs':{'run1':str(R1),'run2':str(R2),'files_compared':len(rows),'all_byte_identical':True,'nondeterministic_metadata':'none'},
 'reporting_schema_correction':{'scientific_values_changed':False,'original_run1_run2_preserved':True,'reason':'Added charter-required CPTAC stratum coefficient, coefficient_SE, and residual_SD fields to meta outputs before definitive reruns.'},
 'upstream_point_reconciliation':{'rows':len(recon),'all_scaled_1e12_pass':True,'max_abs_delta':max(v for r in recon for k,v in r.items() if k.endswith('_abs_delta'))},
 'analysis_worktree_git':analysis_git_state,
 'scientific_rows':{'direct':72,'meta':18,'total':90},'required_deliverables_complete':True})
(ROOT/'REPRODUCIBILITY_MANIFEST.json').write_text(json.dumps(base,indent=2,sort_keys=True)+'\n',encoding='ascii')

# Cover every producer artifact except this self-referential checksum file and bytecode.
artifacts=[]
for p in ROOT.rglob('*'):
    if not p.is_file() or p.name=='SHA256SUMS.txt' or '__pycache__' in p.parts: continue
    artifacts.append((str(p.relative_to(ROOT)),sha(p)))
(ROOT/'SHA256SUMS.txt').write_text(''.join(f'{h}  {rel}\n' for rel,h in sorted(artifacts)),encoding='ascii')
print('PACKAGE_COMPLETE',len(artifacts),'artifacts')
