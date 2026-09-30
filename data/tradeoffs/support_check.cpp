// Independently enumerate all stabilizer labels and sum their support counts.
#include <fstream>
#include <iostream>
#include <vector>
#include <stdexcept>
using namespace std;
int main(int argc,char**argv){try{
 if(argc!=5)throw runtime_error("usage: support_check n q matrix.txt output.csv");
 int n=stoi(argv[1]),q=stoi(argv[2]);if(n<2||n>11||(q!=2&&q!=3))throw runtime_error("unsupported dimensions");
 vector<vector<int>>a(n,vector<int>(n));ifstream in(argv[3]);for(auto&r:a)for(int&x:r)if(!(in>>x)||x<0||x>=q)throw runtime_error("bad matrix");
 for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(a[i][j]!=a[j][i])throw runtime_error("not symmetric");
 int total=1;for(int i=0;i<n;i++)total*=q;
 vector<int>counts(1<<n),weights(n+1);
 for(int label=0;label<total;label++){vector<int>x(n);int v=label,s=0;for(int&i:x){i=v%q;v/=q;}
 for(int i=0;i<n;i++){int y=0;for(int j=0;j<n;j++)y+=a[i][j]*x[j];if(x[i]||y%q)s|=1<<i;}
 counts[s]++;weights[__builtin_popcount((unsigned)s)]++;}
 for(int i=0;i<n;i++)for(int s=0;s<(1<<n);s++)if(s>>i&1)counts[s]+=counts[s^(1<<i)];
 ofstream out(argv[4]);if(!out)throw runtime_error("output unavailable");out<<"mask,size,rank,deficit\n";
 for(int s=1;s<(1<<n);s++){int k=__builtin_popcount((unsigned)s);if(k>n/2)continue;int d=0,v=counts[s];while(v>1&&v%q==0){v/=q;d++;}if(v!=1)throw runtime_error("support count not prime power");out<<s<<','<<k<<','<<k-d<<','<<d<<'\n';}
 cout<<"weights=";for(int w:weights)cout<<w<<' ';cout<<'\n';
 }catch(const exception&e){cerr<<e.what()<<'\n';return 1;}}
