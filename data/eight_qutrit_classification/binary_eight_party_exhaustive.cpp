#include <array>
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
#include <climits>
using namespace std;
int rank_gf2(const array<uint8_t,8>& a, int mask){
 int rows[4],nr=0,cm=((1<<8)-1)^mask;
 for(int i=0;i<8;i++)if(mask&(1<<i))rows[nr++]=a[i]&cm;
 int rank=0;
 for(int col=0;col<8;col++){
  if(!(cm&(1<<col)))continue;
  int p=rank;while(p<nr&&!(rows[p]&(1<<col)))p++;
  if(p==nr)continue;
  swap(rows[p],rows[rank]);
  for(int j=0;j<nr;j++)if(j!=rank&&(rows[j]&(1<<col)))rows[j]^=rows[rank];
  rank++;
  if(rank==nr)break;
 }
 return rank;
}
int main(int argc,char**argv){
 if(argc!=2){cerr<<"usage: binary atlas7_rows.txt\n";return 1;}
 ifstream in(argv[1]);array<uint8_t,8> base{},a{},best{};int bestcost=INT_MAX,seven=0;long long cases=0,winners=0;vector<int> cuts;
 for(int mask=1;mask<256;mask++)if(__builtin_popcount((unsigned)mask)<=4)cuts.push_back(mask);
 while(true){
  int x;if(!(in>>x))break;base[0]=x;
  for(int i=1;i<7;i++){in>>x;base[i]=x;}
  seven++;
  for(int n=0;n<128;n++){
   a=base;a[7]=n;
   for(int i=0;i<7;i++)if(n&(1<<i))a[i]|=1<<7;
   int cost=0;
   for(int mask:cuts){int k=__builtin_popcount((unsigned)mask);int d=k-rank_gf2(a,mask);cost+=d*d;
    if(cost>bestcost)break;}
   if(cost<bestcost){bestcost=cost;best=a;winners=1;}
   else if(cost==bestcost)winners++;
   cases++;
  }
 }
 cout<<"atlas7="<<seven<<" extensions="<<cases<<" cuts="<<cuts.size()<<" min_cost="<<bestcost<<" winning_extensions="<<winners<<"\n";
 cout<<"witness_rows=";for(int i=0;i<8;i++)cout<<(i?" ":"")<<int(best[i]);cout<<"\n";
}
