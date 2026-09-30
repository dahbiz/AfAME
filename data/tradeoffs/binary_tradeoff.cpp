// AfAME - parallel-tempering search and algebraic certification of
// absolutely maximally entangled (AME) states.
// Copyright (C) 2026 Zakaria Dahbi
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.
//
#include <array>
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <climits>
using namespace std;
int rank2(const array<int,8>&a,int s){
 int basis[8]={},r=0;
 for(int i=0;i<8;i++)if(s>>i&1){int v=a[i]&(~s&255);while(v){int j=31-__builtin_clz((unsigned)v);if(basis[j])v^=basis[j];else{basis[j]=v;r++;break;}}}
 return r;
}
int main(int argc,char**argv){
 if(argc!=2)return 2;ifstream in(argv[1]);if(!in)return 3;
 array<int,8> b{},a{};int n7=0;long cases=0;map<int,int> maxfull;map<int,array<int,8>>witness;set<array<int,5>> optimal_profiles;
 while(in>>b[0]){for(int i=1;i<7;i++)if(!(in>>b[i]))return 4;n7++;
 for(int n=0;n<128;n++){a=b;a[7]=n;for(int i=0;i<7;i++)if(n>>i&1)a[i]|=128;
 int c=0,f=0;array<int,5>h{};
 for(int s=1;s<256;s++){int k=__builtin_popcount((unsigned)s);if(k>4)continue;int d=k-rank2(a,s);c+=d*d;if(k==4){h[d]++;f+=d==0;}}
 if(!maxfull.count(c)||f>maxfull[c]){maxfull[c]=f;witness[c]=a;}
 if(c==28)optimal_profiles.insert(h);cases++;
 }}
 cout<<"atlas7="<<n7<<" extensions="<<cases<<"\n";
 int record=-1;for(auto [c,f]:maxfull)if(f>record){cout<<"pareto cost="<<c<<" full="<<f<<" rows=";for(int x:witness[c])cout<<x<<' ';cout<<'\n';record=f;}
 for(auto h:optimal_profiles){cout<<"cost28_profile=";for(int x:h)cout<<x<<' ';cout<<'\n';}
 if(n7!=1044||cases!=133632||maxfull.begin()->first!=28||maxfull[28]!=48||maxfull[32]!=56)return 5;
 for(auto [c,f]:maxfull)if((c<32&&f>48)||f>56)return 6;
}
