#!/usr/bin/env python3
# AfAME - parallel-tempering search and algebraic certification of
# absolutely maximally entangled (AME) states.
# Copyright (C) 2026 Zakaria Dahbi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#
"""Regenerate the 150 search runs and independently verify their certificates.
Requires Python 3.10+ and a C++17 compiler. Never overwrites an existing run.
Timing is machine dependent and is not compared with the frozen CSV.
"""
import argparse,csv,json,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'eight_qutrit_classification'))
from finite_field_certificate import verify

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--cxx',default='c++');a=p.parse_args()
 a.output.mkdir(parents=True,exist_ok=False)
 exe=a.output.resolve()/'recombination_search';baseline=a.output.resolve()/'baseline_search'
 for src,dest in [(ROOT/'sector_recombination_benchmark.cpp',exe),(ROOT/'eight_qutrit_classification/afame_search_verifier.cpp',baseline)]:
  subprocess.run([a.cxx,'-O3','-std=c++17',str(src),'-o',str(dest)],check=True)
 rows=list(csv.DictReader((ROOT/'recombination_benchmark.csv').open()));assert len(rows)==150
 checked=[]
 for row in rows:
  n,seed=int(row['N']),int(row['seed']);name=f'N{n}_seed{seed}';dest=a.output/name
  args=['-N',str(n),'-d','6','--seed',str(seed),'--replicas','4','--steps','500','--restarts','1','--tmin','0.25','--tmax','50','--guide','0.85','--swap-interval','10','--log-interval','500']
  run=subprocess.run([str(exe),*args,'--output',str(dest)],capture_output=True,text=True)
  if run.returncode not in (0,2):raise RuntimeError(run.stderr+run.stdout)
  cert=json.loads((dest/'certificate.json').read_text());verify(cert)
  expected={'joint_cost':cert['best_joint_cost'],'recombined_cost':cert['verified_cost'],'binary_cost':cert['best_sector_costs'][0],'ternary_cost':cert['best_sector_costs'][1],'proposals':cert['proposals'],'gain':cert['best_joint_cost']-cert['verified_cost']}
  for k,v in expected.items():assert int(row[k])==v,(name,k,row[k],v)
  if seed in (1000,1001,1002):
   base=a.output/(name+'_baseline');run=subprocess.run([str(baseline),*args,'--output',str(base)],capture_output=True,text=True)
   if run.returncode not in (0,2):raise RuntimeError(run.stderr+run.stdout)
   bc=json.loads((base/'certificate.json').read_text());verify(bc);assert bc['verified_cost']==cert['best_joint_cost']
  checked.append({'N':n,'seed':seed,**expected})
  if len(checked)%10==0:print(f'{len(checked)}/150 runs reproduced and independently verified',flush=True)
 (a.output/'verification_summary.json').write_text(json.dumps({'runs':len(checked),'baseline_matches':9,'frozen_integer_results_match':True,'checks':checked},indent=2)+'\n')
if __name__=='__main__':main()
