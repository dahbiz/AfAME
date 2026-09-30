#!/usr/bin/env python3
"""Reproduce the added trade-offs and eleven-qutrit certificate.
No third-party Python packages. Requires C++17. Run from the AfAME checkout:
python3 data/tradeoffs/run_checks.py --output NEW_DIRECTORY
Classification coverage is inherited from the cited published classifications.
"""
import argparse, collections, csv, hashlib, itertools, json, math, subprocess, sys
from pathlib import Path
HERE=Path(__file__).resolve().parent

def save(p,x): p.write_text(json.dumps(x,indent=2)+'\n')
def run(cmd,log):
 r=subprocess.run(list(map(str,cmd)),capture_output=True,text=True)
 log.write_text(r.stdout+r.stderr)
 if r.returncode: raise RuntimeError(f'{cmd}: {r.stderr}')
 return r.stdout

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--repo',type=Path,default=HERE.parents[1]);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--cxx',default='c++');args=ap.parse_args()
 repo=args.repo.resolve();exact=HERE.parent/'eight_qutrit_classification';out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
 sys.path.insert(0,str(exact));from finite_field_certificate import Field,verify
 from analyze_classification import read_records
 for name in ['binary_tradeoff','binary_tradeoff_span','support_check']:
  run([args.cxx,'-O3','-std=c++17','-Wall','-Wextra','-Wpedantic',HERE/(name+'.cpp'),'-o',out/name],out/(name+'_build.log'))
 for name in ['binary_tradeoff','binary_tradeoff_span']:
  run([out/name,exact/'binary_graph_atlas_rows.txt'],out/(name+'.log'))
 print('Binary cost/count Pareto boundary independently verified',flush=True)
 records,_=read_records(exact/'source_records_length8.tsv')
 run([args.cxx,'-O3','-std=c++17',exact/'verify_stabilizer_support.cpp','-o',out/'ternary_support'],out/'ternary_build.log')
 run([out/'ternary_support',exact/'source_records_length8.tsv',out/'ternary_support.csv'],out/'ternary_support.log')
 independent={(int(r['class']),int(r['mask'])):int(r['deficit']) for r in csv.DictReader((out/'ternary_support.csv').open())}
 cuts=[(k,s,[i for i in range(8) if i not in s],sum(1<<i for i in s)) for k in range(1,5) for s in itertools.combinations(range(8),k)]
 f=Field(3,1,[]);profiles=collections.defaultdict(list);points=[]
 for rec in records:
  h=collections.Counter();cost=0
  for k,s,t,m in cuts:
   d=k-f.rank([[rec['matrix'][i][j] for j in t] for i in s]);assert d==independent[rec['index'],m];cost+=d*d
   if k==4:h[d]+=1
  tail=tuple(sum(h[j] for j in range(d+1)) for d in range(4))
  profiles[tail].append(rec['index']);points.append({'class':rec['index'],'cost':cost,'n':list(h.get(d,0) for d in range(5)),'tails':tail})
 frontier={p:ids for p,ids in profiles.items() if not any(q!=p and all(x>=y for x,y in zip(q,p)) for q in profiles)}
 assert set(frontier)=={(54,70,70,70),(60,68,70,70),(66,66,70,70)}
 assert all(p['tails'][0]+3*p['tails'][1]<=264 for p in points)
 assert sum(len(x) for x in frontier.values())==10
 save(out/'ternary_all_profiles.json',points);save(out/'ternary_frontier.json',[{'tails':p,'classes':ids} for p,ids in sorted(frontier.items())])
 print('All 817 ternary classes: complete tail frontier checked twice',flush=True)
 def check_matrix(a,q,name):
  n=len(a);field=Field(q,1,[]);path=out/(name+'_matrix.txt');path.write_text('\n'.join(' '.join(map(str,r)) for r in a)+'\n')
  log=run([out/'support_check',n,q,path,out/(name+'_support.csv')],out/(name+'_support.log'))
  independent={int(r['mask']):int(r['deficit']) for r in csv.DictReader((out/(name+'_support.csv')).open())}
  ds={};hist=collections.Counter()
  for k in range(1,n//2+1):
   for s in itertools.combinations(range(n),k):
    t=[i for i in range(n) if i not in s];m=sum(1<<i for i in s);d=k-field.rank([[a[i][j] for j in t] for i in s]);assert d==independent[m];ds[m]=d;hist[k,d]+=1
  assert len(ds)==len(independent)
  return ds,dict((f'{k},{d}',v) for (k,d),v in sorted(hist.items())),log
 def certificate(matrices,primes,ds,name):
  n=len(matrices[0]);cost=sum(d*d for x in ds for d in x.values());bad=sum(any(x[m] for x in ds) for m in ds[0])
  c={'schema':1,'N':n,'d':math.prod(primes),'ame':bad==0,'cuts':len(ds[0]),'failing_cuts':bad,'verified_cost':cost,'fields':[{'p':q,'m':1,'q':q,'modulus':[],'matrix':a} for a,q in zip(matrices,primes)]}
  verify(c);save(out/(name+'_certificate.json'),c);return c
 a=json.loads((HERE/'danielsen11_matrix.json').read_text());d,h,weights=check_matrix(a,3,'danielsen11')
 c=certificate([a],[3],[d],'danielsen11');assert c['verified_cost']==6 and h=={'1,0':11,'2,0':55,'3,0':165,'4,0':330,'5,0':456,'5,1':6}
 assert list(map(int,weights.split('=')[1].split()))==[1,0,0,0,0,12,888,3960,14970,42500,66240,48576]
 bad=[m for m,v in d.items() if v];inc=[sum(bool(m>>i&1) for m in bad) for i in range(11)];assert inc==[6,6,2,2,2,2,2,2,2,2,2]
 save(out/'eleven_summary.json',{'cost':6,'histogram':h,'bad_subsets_one_based':[[i+1 for i in range(11) if m>>i&1] for m in bad],'bad_four_erasure_patterns_by_reference':inc,'four_erasure_patterns_per_reference':210,'mean_leakage_in_log3_units':'1/77','uniformity_forcing_bound':sum(math.comb(7,j) for j in range(2)),'classification_input':'Danielsen 2012 Eq. W11,alpha: alpha in {0,6,...,50,54,60}; unique alpha=0 class. Not regenerated.'})
 print('Published eleven-qutrit matrix: 1,023 cut ranks and all 177,147 stabilizer labels agree',flush=True)
 products=json.loads((HERE/'aligned_products.json').read_text());result={}
 for name,x in products.items():
  d2,_,_=check_matrix(x['binary_matrix'],2,name+'_binary');d3,_,_=check_matrix(x['ternary_matrix'],3,name+'_ternary');c=certificate([x['binary_matrix'],x['ternary_matrix']],[2,3],[d2,d3],name)
  full=sum(d2[m]==d3[m]==0 for m in d2 if m.bit_count()==4);assert full==x['fully_mixed_balanced']
  result[name]={'cost':c['verified_cost'],'full_balanced':full,'histogram':x['balanced_joint_histogram']}
 assert [(result[x]['cost'],result[x]['full_balanced']) for x in ['cost_optimal','T4']]==[(44,48),(48,56)]
 save(out/'product_summary.json',result)
 hashes={str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [exact/'source_records_length8.tsv',exact/'binary_graph_atlas_rows.txt']}
 save(out/'run_summary.json',{'binary_pareto':[[28,48],[32,56]],'qutrit_tail_pareto':list(sorted(frontier)),'product_pareto':[[44,48],[48,56]],'eleven_cost':6,'ternary_cut_ranks_compared':len(independent),'classification_inputs_sha256':hashes,'coverage':'Inherited from published classifications; eleven-qutrit lower bound additionally uses the manuscript proof and published weight-enumerator list.'})
 print('Product witnesses independently verified; all checks passed',flush=True)
if __name__=='__main__':main()
