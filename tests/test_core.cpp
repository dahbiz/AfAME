#define AME_NO_MAIN
#include "../AME7DIA_PT.cpp"
#include <cassert>
using namespace ame;
int main() {
    require(!irreducible({2,1,1},11),"Old GF121 polynomial not rejected");
    require(!irreducible({1,1,1,1},7),"Old GF343 polynomial not rejected");
    require(!irreducible({2,1,0,0,1},5),"Old GF625 polynomial not rejected");
    RNG rng(417);
    for(auto pm:std::vector<std::pair<int,int>>{{2,1},{73,1},{137,1},{257,1},{2,2},{3,2},{11,2},{7,3},{5,4},{2,9},{2,10}}) {
        Field f(pm.first,pm.second);
        for(int a=1;a<f.q;++a)require(f.times(a,f.inv[a])==1,"Inverse identity");
        for(int i=0;i<400;++i) {
            int a=int(rng.bounded(f.q)),b=int(rng.bounded(f.q)),c=int(rng.bounded(f.q));
            require(f.times(a,f.minus(b,c))==f.minus(f.times(a,b),f.times(a,c)),"Distributivity");
            std::vector<Element> mat{Element(a),Element(b),Element(c),Element(rng.bounded(f.q))};
            int det=f.minus(f.times(mat[0],mat[3]),f.times(mat[1],mat[2]));
            int expect=det?2:((mat[0]||mat[1]||mat[2]||mat[3])?1:0);
            require(rank(mat,2,2,f)==expect,"2x2 determinant/rank mismatch");
        }
        if(f.q==257)require(rank({256},1,1,f)==1,"F257 regression");
    }
    int moves=0;
    for(int dimension:{2,4,6,9,12,121,257,10001}) {
        std::vector<Field> fs;for(auto pm:factors(dimension))fs.emplace_back(pm.first,pm.second);
        Geometry g(5,1,int(fs.size()),512);Replica r(512,.8,0);r.init(g,fs,true);
        for(int i=0;i<400;++i) {
            r.step(g,fs,.85);Cache full;full.rebuild(r.state,g,fs);
            require(full.total==r.cache.total&&full.costs==r.cache.costs&&full.ranks==r.cache.ranks&&full.fail==r.cache.fail&&full.fail_sum==r.cache.fail_sum,"Cache/guidance mismatch");
            for(int x:r.cache.fail)require(x>=0,"Negative guide count");++moves;
        }
    }
    // Exact detailed-balance comparison for all neighboring binary three-party graphs.
    std::vector<Field> fs;fs.emplace_back(2,1);Geometry g(3,1,1,512);
    for(int label=0;label<8;++label) {
        State s(1,std::vector<Element>(9));
        for(int e=0;e<3;++e){auto [u,v]=g.edges[e];s[0][u*3+v]=s[0][v*3+u]=(label>>e)&1;}
        Cache a;a.rebuild(s,g,fs);
        for(int e=0;e<3;++e) {
            auto [u,v]=g.edges[e];State t=s;t[0][u*3+v]^=1;t[0][v*3+u]^=1;Cache b;b.rebuild(t,g,fs);
            double qa=edge_probability(a,e,.85),qb=edge_probability(b,e,.85),temp=.7;
            double lp=-double(b.total-a.total)/temp+std::log(qb/qa);
            double lhs=std::exp(-double(a.total)/temp)*qa*std::exp(std::min(0.0,lp));
            double rhs=std::exp(-double(b.total)/temp)*qb*std::exp(std::min(0.0,-lp));
            require(std::abs(lhs-rhs)<1e-12,"Guided detailed balance");
        }
    }
    require(std::abs(swap_log(1,2,0,10)+5)<1e-12,"Swap acceptance sign");
    Workers workers(3);std::vector<int> counts(3);
    for(int i=0;i<50;++i)workers.run([&](size_t r){++counts[r];});
    for(auto n:counts)require(n==50,"Worker dispatch");
    bool caught=false;try{workers.run([](size_t r){if(r==1)throw std::runtime_error("test");});}catch(const std::runtime_error&){caught=true;}
    require(caught,"Worker error propagation");workers.run([](size_t){});
    std::cout<<"PASS: field regression tests, "<<moves<<" cache/guide checks, detailed balance, swap sign, worker completion/error recovery\n";
}
