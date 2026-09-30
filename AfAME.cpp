// AfAME - Developed by Dr. Zakaria Dahbi.
// Parallel tempering for quadratic-phase absolutely maximally entangled states.
// Refactoring and testing assisted by Codex; see AUTHORS.md.
// See README.md for provenance, phase conventions, validation, and limitations.
// C++17 standard-library threads; no MPI/PETSc dependency.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace ame {
using Clock = std::chrono::steady_clock;
using Element = uint16_t;
using Cost = int64_t;
using Poly = std::vector<int>; // coefficients in increasing degree
constexpr int MAX_FIELD = 1024;
constexpr const char* VERSION = "AME7DIA_PT-1.0";
void require(bool ok, const std::string& what) { if (!ok) throw std::runtime_error(what); }
int mod(int x,int p) { x%=p; return x<0?x+p:x; }
bool prime(int p) { if(p<2)return false;for(int i=2;int64_t(i)*i<=p;++i)if(p%i==0)return false;return true; }
int power_mod(int a,int e,int p) { int r=1;for(;e;e>>=1,a=a*a%p)if(e&1)r=r*a%p;return r; }
void trim(Poly& a) { while(!a.empty() && a.back()==0)a.pop_back(); }
Poly rem(Poly a,const Poly& b,int p) {
    require(!b.empty(),"Polynomial division by zero");
    int inv=power_mod(b.back(),p-2,p);
    while(a.size()>=b.size()) {
        int shift=int(a.size()-b.size()), c=a.back()*inv%p;
        for(size_t j=0;j<b.size();++j)a[shift+j]=mod(a[shift+j]-c*b[j],p);
        trim(a);
    }
    return a;
}
Poly pmul(const Poly& a,const Poly& b,const Poly& f,int p) {
    if(a.empty()||b.empty())return {};
    Poly v(a.size()+b.size()-1);
    for(size_t i=0;i<a.size();++i)for(size_t j=0;j<b.size();++j)v[i+j]=(v[i+j]+a[i]*b[j])%p;
    trim(v);return rem(v,f,p);
}
Poly ppow(Poly a,int e,const Poly& f,int p) {
    Poly r{1};for(;e;e>>=1,a=pmul(a,a,f,p))if(e&1)r=pmul(r,a,f,p);return r;
}
Poly gcd(Poly a,Poly b,int p) {while(!b.empty()){auto r=rem(a,b,p);a=std::move(b);b=std::move(r);}return a;}
bool irreducible(const Poly& f,int p) {
    int m=int(f.size())-1;if(m<1||f.back()!=1||!prime(p))return false;
    if(m==1)return true;
    Poly x{0,1}, h=x;
    for(int i=1;i<=m;++i) {
        h=ppow(h,p,f,p);Poly d=h;if(d.size()<2)d.resize(2);
        d[1]=mod(d[1]-1,p);trim(d);
        if(i<=m/2 && gcd(f,d,p).size()>1)return false;
        if(i==m && !d.empty())return false;
    }
    return true;
}
Poly modulus_for(int p,int m,int q) {
    for(int label=1;label<q;++label) {
        Poly f(m+1);f[m]=1;int v=label;
        for(int i=0;i<m;++i){f[i]=v%p;v/=p;}
        if(irreducible(f,p))return f;
    }
    throw std::runtime_error("No irreducible polynomial found");
}
struct Field {
    int p,m,q;Poly modulus;
    std::vector<Element> mul,sub,inv;
    Field(int pp,int mm):p(pp),m(mm),q(1) {
        require(prime(p)&&m>0,"Invalid field specification");
        for(int i=0;i<m;++i){require(q<=MAX_FIELD/p,"Field size exceeds 1024");q*=p;}
        if(m>1)modulus=modulus_for(p,m,q);
        mul.resize(size_t(q)*q);sub.resize(size_t(q)*q);inv.resize(q);
        for(int a=0;a<q;++a)for(int b=0;b<q;++b) {
            mul[size_t(a)*q+b]=Element(mul_raw(a,b));
            int aa=a,bb=b,r=0,pw=1;
            for(int i=0;i<m;++i){r+=mod(aa%p-bb%p,p)*pw;aa/=p;bb/=p;pw*=p;}
            sub[size_t(a)*q+b]=Element(r);
        }
        for(int a=1;a<q;++a) {
            int r=1,b=a,e=q-2;for(;e;e>>=1,b=times(b,b))if(e&1)r=times(r,b);
            inv[a]=Element(r);require(times(a,r)==1,"Invalid inverse in field construction");
        }
    }
    int mul_raw(int a,int b) const {
        if(m==1)return a*b%p;
        Poly aa(m),bb(m);for(int i=0;i<m;++i){aa[i]=a%p;bb[i]=b%p;a/=p;b/=p;}
        trim(aa);trim(bb);auto r=pmul(aa,bb,modulus,p);int v=0,pw=1;
        for(int c:r){v+=c*pw;pw*=p;}return v;
    }
    int times(int a,int b) const { return mul[size_t(a)*q+b]; }
    int minus(int a,int b) const { return sub[size_t(a)*q+b]; }
};
std::vector<std::pair<int,int>> factors(int d) {
    require(d>=2,"d must be >= 2");std::vector<std::pair<int,int>> f;
    for(int p=2;int64_t(p)*p<=d;++p)if(d%p==0){int m=0;do{d/=p;++m;}while(d%p==0);f.emplace_back(p,m);}
    if(d>1)f.emplace_back(d,1);return f;
}
int rank(std::vector<Element> a,int rows,int cols,const Field& f) {
    int r=0;
    for(int c=0;c<cols && r<rows;++c) {
        int pivot=r;while(pivot<rows && !a[pivot*cols+c])++pivot;
        if(pivot==rows)continue;
        for(int j=0;j<cols;++j)std::swap(a[r*cols+j],a[pivot*cols+j]);
        int v=f.inv[a[r*cols+c]];
        for(int j=c;j<cols;++j)a[r*cols+j]=Element(f.times(a[r*cols+j],v));
        for(int i=r+1;i<rows;++i)if(a[i*cols+c]) {
            int v2=a[i*cols+c];for(int j=c;j<cols;++j)a[i*cols+j]=Element(f.minus(a[i*cols+j],f.times(v2,a[r*cols+j])));
        }
        ++r;
    }
    return r;
}
uint64_t choose(int n,int k) {uint64_t r=1;for(int i=1;i<=k;++i)r=r*(n-k+i)/i;return r;}
struct Cut { uint64_t mask;int k;std::vector<uint16_t> edges; };
struct Geometry {
    int n;std::vector<std::pair<int,int>> edges;
    std::vector<Cut> cuts;std::vector<std::vector<uint32_t>> affected;
    Geometry(int nn,int replicas,int nf,double limit_mb):n(nn) {
        require(n>=2 && n<=24,"N must be in [2,24]; exhaustive cut enumeration is exponential");
        uint64_t count=0,incidences=0;
        for(int k=1;k<=n/2;++k){auto c=choose(n,k);count+=c;incidences+=c*k*(n-k);}
        long double bytes=incidences*12.0L+count*(128.0L+replicas*(24.0L+4*nf));
        require(bytes<=limit_mb*1024*1024,"Estimated cut/cache memory exceeds --memory-mb; reduce N/replicas or raise the explicit limit");
        for(int u=0;u<n;++u)for(int v=u+1;v<n;++v)edges.emplace_back(u,v);
        affected.resize(edges.size());cuts.reserve(count);
        std::function<void(int,int,uint64_t,int)> enumerate=[&](int start,int left,uint64_t mask,int k) {
            if(left){for(int i=start;i<=n-left;++i)enumerate(i+1,left-1,mask|(1ULL<<i),k);return;}
            uint32_t idx=uint32_t(cuts.size());Cut c{mask,k,{}};c.edges.reserve(k*(n-k));
            for(size_t e=0;e<edges.size();++e) {
                auto [u,v]=edges[e];if(((mask>>u)^(mask>>v))&1){c.edges.push_back(uint16_t(e));affected[e].push_back(idx);}
            }
            cuts.push_back(std::move(c));
        };
        for(int k=1;k<=n/2;++k)enumerate(0,k,0,k);
    }
};
using State = std::vector<std::vector<Element>>;
int cut_rank(const State& s,size_t fidx,const Cut& c,const Geometry& g,const Field& f) {
    std::vector<Element> a;a.reserve(c.k*(g.n-c.k));
    for(int i=0;i<g.n;++i)if((c.mask>>i)&1)
        for(int j=0;j<g.n;++j)if(!((c.mask>>j)&1))a.push_back(s[fidx][i*g.n+j]);
    return rank(std::move(a),c.k,g.n-c.k,f);
}
struct Cache {
    std::vector<int> ranks,costs,fail;
    Cost total=0;int64_t fail_sum=0;
    void rebuild(const State& s,const Geometry& g,const std::vector<Field>& fs) {
        ranks.assign(g.cuts.size()*fs.size(),0);costs.assign(g.cuts.size(),0);fail.assign(g.edges.size(),0);total=0;fail_sum=0;
        for(size_t b=0;b<g.cuts.size();++b) {
            int c=0;for(size_t f=0;f<fs.size();++f){int r=cut_rank(s,f,g.cuts[b],g,fs[f]);ranks[b*fs.size()+f]=r;int d=g.cuts[b].k-r;c+=d*d;}
            costs[b]=c;total+=c;if(c)for(int e:g.cuts[b].edges){++fail[e];++fail_sum;}
        }
    }
};
struct Change {uint32_t cut;int old_rank,new_rank,old_cost,new_cost;};
std::vector<Change> changes(const State& s,const Cache& c,int edge,int fi,const Geometry& g,const std::vector<Field>& fs) {
    std::vector<Change> out;out.reserve(g.affected[edge].size());
    for(auto b:g.affected[edge]) {
        int old=c.ranks[b*fs.size()+fi],r=cut_rank(s,fi,g.cuts[b],g,fs[fi]),k=g.cuts[b].k;
        if(r!=old)out.push_back({b,old,r,c.costs[b],c.costs[b]+(k-r)*(k-r)-(k-old)*(k-old)});
    }
    return out;
}
void apply(Cache& c,const std::vector<Change>& updates,int fi,const Geometry& g,size_t nf,bool forward=true) {
    for(auto v:updates) {
        int old=forward?v.old_cost:v.new_cost,nw=forward?v.new_cost:v.old_cost;
        c.total+=nw-old;c.costs[v.cut]=nw;c.ranks[v.cut*nf+fi]=forward?v.new_rank:v.old_rank;
        int d=(nw>0)-(old>0);
        if(d)for(auto e:g.cuts[v.cut].edges){c.fail[e]+=d;c.fail_sum+=d;}
    }
}
struct RNG {
    std::mt19937_64 engine;
    explicit RNG(uint64_t seed):engine(seed){}
    uint64_t bounded(uint64_t n) {
        require(n>0,"Empty random range");uint64_t threshold=(uint64_t(0)-n)%n;
        for(;;){uint64_t v=engine();if(v>=threshold)return v%n;}
    }
    double unit() {return double(engine()>>11)*0x1.0p-53;} // [0,1)
};
uint64_t mix(uint64_t x) {x+=0x9e3779b97f4a7c15ULL;x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;x=(x^(x>>27))*0x94d049bb133111ebULL;return x^(x>>31);}
double edge_probability(const Cache& c,int e,double guide) {
    double uniform=1.0/c.fail.size();
    return c.fail_sum ? (1-guide)*uniform+guide*double(c.fail[e])/c.fail_sum : uniform;
}
int select_edge(const Cache& c,double guide,RNG& rng) {
    if(!c.fail_sum || rng.unit()>=guide)return int(rng.bounded(c.fail.size()));
    int64_t target=int64_t(rng.bounded(uint64_t(c.fail_sum)));
    for(size_t e=0;e<c.fail.size();++e){target-=c.fail[e];if(target<0)return int(e);}
    throw std::runtime_error("Inconsistent edge guidance");
}
bool accept_log(double lp,RNG& rng) { return lp>=0 || std::log(rng.unit())<lp; }
double swap_log(double ti,double tj,Cost ei,Cost ej) {return (1/ti-1/tj)*double(ei-ej);}
struct Replica {
    State state,best;Cache cache;RNG rng;double temperature;
    Cost best_cost=std::numeric_limits<Cost>::max();uint64_t seed,proposed=0,accepted=0,uphill=0,best_step=0;
    int walker;
    Replica(uint64_t sd,double t,int id):rng(sd),temperature(t),seed(sd),walker(id){}
    void remember() {if(cache.total<best_cost){best_cost=cache.total;best=state;best_step=proposed;}}
    void init(const Geometry& g,const std::vector<Field>& fs,bool diagonal) {
        state.assign(fs.size(),std::vector<Element>(g.n*g.n));
        for(size_t f=0;f<fs.size();++f)for(int u=0;u<g.n;++u)for(int v=diagonal?u:u+1;v<g.n;++v)
            state[f][u*g.n+v]=state[f][v*g.n+u]=Element(rng.bounded(fs[f].q));
        cache.rebuild(state,g,fs);remember();
    }
    void step(const Geometry& g,const std::vector<Field>& fs,double guide) {
        int e=select_edge(cache,guide,rng),fi=int(rng.bounded(fs.size()));auto [u,v]=g.edges[e];
        auto& mat=state[fi];Element old=mat[u*g.n+v];
        double forward=edge_probability(cache,e,guide);Cost old_cost=cache.total;
        mat[u*g.n+v]=mat[v*g.n+u]=Element((old+1+rng.bounded(fs[fi].q-1))%fs[fi].q);
        auto update=changes(state,cache,e,fi,g,fs);apply(cache,update,fi,g,fs.size());
        double reverse=edge_probability(cache,e,guide);
        double lp=-double(cache.total-old_cost)/temperature+std::log(reverse/forward);
        ++proposed;
        if(accept_log(lp,rng)){++accepted;if(cache.total>old_cost)++uphill;remember();}
        else {apply(cache,update,fi,g,fs.size(),false);mat[u*g.n+v]=mat[v*g.n+u]=old;}
    }
};

// Persistent threads; dispatch finishes before replica configurations are exchanged.
class Workers {
    std::mutex mutex;std::condition_variable start,done;
    bool stopping=false;size_t generation=0,remaining=0;
    std::function<void(size_t)> task;std::exception_ptr error;std::vector<std::thread> threads;
public:
    explicit Workers(size_t n) {
        try {for(size_t id=0;id<n;++id)threads.emplace_back([this,id]{
            size_t seen=0;std::unique_lock<std::mutex> lock(mutex);
            for(;;){start.wait(lock,[&]{return stopping||generation!=seen;});if(stopping)return;
                seen=generation;auto job=task;lock.unlock();std::exception_ptr caught;
                try{job(id);}catch(...){caught=std::current_exception();}
                lock.lock();if(caught&&!error)error=caught;if(--remaining==0)done.notify_one();
            }
        });}catch(...){shutdown();throw;}
    }
    void run(std::function<void(size_t)> job) {
        std::unique_lock<std::mutex> lock(mutex);task=std::move(job);error=nullptr;remaining=threads.size();++generation;start.notify_all();
        done.wait(lock,[&]{return remaining==0;});if(error)std::rethrow_exception(error);
    }
    void shutdown(){ {std::lock_guard<std::mutex> lock(mutex);stopping=true;}start.notify_all();for(auto& t:threads)if(t.joinable())t.join();}
    ~Workers(){shutdown();}
};
struct Options {
    int n=7,d=2,replicas=4,steps=5000,restarts=20,seed_tries=1,swap_interval=25,log_interval=200;
    double tmin=.25,tmax=50,guide=.85,memory_mb=512;
    uint64_t seed=1;bool diagonal=false,verify_only=false,keep_going=false,debug=false,print_entropy=false;
    std::string output="",input="";
};
Options parse(int argc,char**argv) {
    Options o;
    auto integer=[](const std::string& s){size_t n;long long v=std::stoll(s,&n);require(n==s.size()&&v>=0&&v<=std::numeric_limits<int>::max(),"Invalid integer: "+s);return int(v);};
    auto real=[](const std::string& s){size_t n;double v=std::stod(s,&n);require(n==s.size()&&std::isfinite(v),"Invalid number: "+s);return v;};
    for(int i=1;i<argc;++i) {
        std::string a=argv[i];auto value=[&](){require(i+1<argc,"Missing value for "+a);return std::string(argv[++i]);};
        if(a=="-N"||a=="--N")o.n=integer(value());
        else if(a=="-d"||a=="--d")o.d=integer(value());
        else if(a=="-replicas"||a=="--replicas")o.replicas=integer(value());
        else if(a=="-steps"||a=="--steps")o.steps=integer(value());
        else if(a=="-restarts"||a=="--restarts")o.restarts=integer(value());
        else if(a=="-seed_tries"||a=="--seed-tries")o.seed_tries=integer(value());
        else if(a=="-swap_interval"||a=="--swap-interval")o.swap_interval=integer(value());
        else if(a=="-log_interval"||a=="--log-interval")o.log_interval=integer(value());
        else if(a=="--tmin")o.tmin=real(value());
        else if(a=="--tmax")o.tmax=real(value());
        else if(a=="--guide")o.guide=real(value());
        else if(a=="--memory-mb")o.memory_mb=real(value());
        else if(a=="-seed"||a=="--seed") {auto s=value();size_t end;require(!s.empty()&&s[0]!='-',"Invalid seed");o.seed=std::stoull(s,&end);require(end==s.size(),"Invalid seed");}
        else if(a=="--output"||a=="-outprefix")o.output=value();
        else if(a=="--input")o.input=value();
        else if(a=="--verify-only")o.verify_only=true;
        else if(a=="--keep-going")o.keep_going=true;
        else if(a=="--debug-checks")o.debug=true;
        else if(a=="--diagonal")o.diagonal=true;
        else if(a=="-print_search_entropy")o.print_entropy=true;
        else if(a=="-tensor_phases") {if(i+1<argc && std::string(argv[i+1])=="1")++i;}
        else throw std::runtime_error("Unknown option: "+a+" (use --help)");
    }
    require(o.n>=2&&o.n<=24&&o.d>=2,"Require 2 <= N <= 24 and d >= 2");
    require(o.replicas>=1&&o.replicas<=64,"Require 1 <= replicas <= 64");
    require(o.steps>0&&o.restarts>0&&o.seed_tries>0&&int64_t(o.restarts)*o.seed_tries<=1000000,"Positive steps/restarts/seed tries required; at most 1,000,000 restarts");
    require(o.swap_interval>0&&o.log_interval>0,"Intervals must be positive");
    require(o.tmin>=1e-6&&o.tmax>=o.tmin&&o.tmax<=1e9,"Require 1e-6 <= tmin <= tmax <= 1e9");
    require(o.replicas==1||o.tmax>o.tmin,"Multiple replicas require tmax > tmin");
    require(o.guide>=0&&o.guide<1,"Guide must be in [0,1)");require(o.memory_mb>0,"Memory limit must be positive");
    require(!o.verify_only||!o.input.empty(),"--verify-only requires --input");
    require(o.verify_only||o.input.empty(),"--input currently requires --verify-only");
    return o;
}
std::ofstream output_file(const std::filesystem::path& p) {std::ofstream f(p);require(bool(f),"Cannot write "+p.string());f<<std::setprecision(17);f.exceptions(std::ios::badbit|std::ios::failbit);return f;}
State read_state(const Options& o,const std::vector<Field>& fs) {
    std::ifstream in(o.input);require(bool(in),"Cannot read input matrix");std::vector<int> values;std::string line;
    while(std::getline(in,line)){auto pos=line.find('#');if(pos!=std::string::npos)line.resize(pos);std::istringstream ss(line);int v;while(ss>>v)values.push_back(v);require(ss.eof(),"Invalid matrix token");}
    require(values.size()==size_t(o.n*o.n),"Input requires exactly N*N integer entries");
    if(fs.size()>1)for(const auto& f:fs)require(f.m==1,"CRT input is supported only for square-free d; use the JSON verifier for extension-field tensor certificates");
    State s(fs.size(),std::vector<Element>(values.size()));
    for(int i=0;i<o.n;++i)for(int j=0;j<o.n;++j){int v=values[i*o.n+j];require(v>=0&&v<o.d,"Input element outside [0,d)");require(v==values[j*o.n+i],"Input must be symmetric");for(size_t f=0;f<fs.size();++f)s[f][i*o.n+j]=Element(v%fs[f].q);}
    return s;
}
void matrix_json(std::ostream& out,const std::vector<Element>& p,int n) {
    out<<'[';for(int i=0;i<n;++i){if(i)out<<',';out<<'[';for(int j=0;j<n;++j){if(j)out<<',';out<<p[i*n+j];}out<<']';}out<<']';
}
struct Timings {double setup=0,search=0,verification=0,export_data=0,total=0;};
void export_data(const Options& o,const Geometry& g,const std::vector<Field>& fs,const State& s,const Cache& verified,const std::filesystem::path& dir) {
    auto entropy=output_file(dir/"entropy_bips.tsv"),hist=output_file(dir/"rank_histogram.tsv"),per=output_file(dir/"entropy_per_k.tsv"),heat=output_file(dir/"hot_edges.tsv");
    entropy<<"cut\tk\tmask\trank_sum\tmax_rank_sum\tentropy_nats\tentropy_bits\tmax_entropy_bits\tcost\n";
    hist<<"k\tfactor\trank\tcount\tfraction\n";per<<"k\tcuts\tmean_bits\tstd_bits\tmin_bits\tmax_bits\tpassing\n";
    std::vector<double> bits(g.cuts.size());
    for(size_t b=0;b<g.cuts.size();++b){double en=0;int ranks=0;for(size_t f=0;f<fs.size();++f){int r=verified.ranks[b*fs.size()+f];en+=r*std::log(fs[f].q);ranks+=r;}bits[b]=en/std::log(2.0);entropy<<b<<'\t'<<g.cuts[b].k<<'\t'<<g.cuts[b].mask<<'\t'<<ranks<<'\t'<<g.cuts[b].k*fs.size()<<'\t'<<en<<'\t'<<bits[b]<<'\t'<<g.cuts[b].k*std::log2(o.d)<<'\t'<<verified.costs[b]<<'\n';}
    for(int k=1;k<=o.n/2;++k){size_t count=0,pass=0;double mean=0,m2=0,lo=1e100,hi=0;std::vector<std::vector<int>> h(fs.size(),std::vector<int>(k+1));
        for(size_t b=0;b<g.cuts.size();++b)if(g.cuts[b].k==k){++count;pass+=verified.costs[b]==0;double delta=bits[b]-mean;mean+=delta/count;m2+=delta*(bits[b]-mean);lo=std::min(lo,bits[b]);hi=std::max(hi,bits[b]);for(size_t f=0;f<fs.size();++f)++h[f][verified.ranks[b*fs.size()+f]];}
        per<<k<<'\t'<<count<<'\t'<<mean<<'\t'<<std::sqrt(std::max(0.0,count>1?m2/(count-1):0))<<'\t'<<lo<<'\t'<<hi<<'\t'<<pass<<'\n';
        for(size_t f=0;f<fs.size();++f)for(int r=0;r<=k;++r)if(h[f][r])hist<<k<<'\t'<<f<<'\t'<<r<<'\t'<<h[f][r]<<'\t'<<double(h[f][r])/count<<'\n';
    }
    heat<<"i\tj\tfailing_cuts\n";for(size_t e=0;e<g.edges.size();++e)heat<<g.edges[e].first<<'\t'<<g.edges[e].second<<'\t'<<verified.fail[e]<<'\n';
    for(size_t f=0;f<fs.size();++f){auto out=output_file(dir/("factor_"+std::to_string(f)+".txt"));out<<"# GF("<<fs[f].p<<'^'<<fs[f].m<<"), base-p polynomial labels; modulus (low to high):";for(int c:fs[f].modulus)out<<' '<<c;out<<'\n';for(int i=0;i<o.n;++i){for(int j=0;j<o.n;++j)out<<s[f][i*o.n+j]<<' ';out<<'\n';}}
    bool squarefree=true;for(auto& f:fs)squarefree&=f.m==1;
    if(squarefree){auto out=output_file(dir/"crt_matrix.txt");out<<"# CRT residue matrix; same cut ranks. Sector character weights may differ.\n";
        for(int i=0;i<o.n;++i){for(int j=0;j<o.n;++j){int64_t value=0;for(size_t f=0;f<fs.size();++f){int quotient=o.d/fs[f].p;int inv=power_mod(quotient%fs[f].p,fs[f].p-2,fs[f].p);value=(value+int64_t(s[f][i*o.n+j])*quotient*inv)%o.d;}out<<value<<' ';}out<<'\n';}}
}
int run(const Options& o) {
    auto start=Clock::now();std::filesystem::path dir=o.output;
    if(dir.empty())dir="AME_N"+std::to_string(o.n)+"_d"+std::to_string(o.d)+"_seed"+std::to_string(o.seed)+"_"+std::to_string(std::chrono::system_clock::now().time_since_epoch().count());
    require(!std::filesystem::exists(dir),"Output directory already exists; select a new --output to preserve earlier runs");
    std::vector<Field> fs;for(auto [p,m]:factors(o.d))fs.emplace_back(p,m);
    Geometry g(o.n,o.verify_only?1:o.replicas,int(fs.size()),o.memory_mb);
    std::filesystem::create_directories(dir);
    std::cout<<"AfAME - Developed by Dr. Zakaria Dahbi\n";
    std::cout<<VERSION<<" N="<<o.n<<" d="<<o.d<<" cuts="<<g.cuts.size()<<" replicas="<<(o.verify_only?0:o.replicas)<<"\n";
    for(auto& f:fs){std::cout<<"GF("<<f.q<<") modulus coefficients:";for(int c:f.modulus)std::cout<<' '<<c;std::cout<<'\n';}
    auto conv=output_file(dir/"convergence.tsv"),swaps=output_file(dir/"swaps.tsv"),restarts=output_file(dir/"restart_stats.tsv");
    conv<<"restart\tstep_per_replica\tslot\ttemperature\twalker\tcost\tbest_cost\tproposed\taccepted\tuphill\n";
    swaps<<"restart\tround\tstep\ti\tj\tenergy_i\tenergy_j\tlog_acceptance\taccepted\twalker_i_before\twalker_j_before\n";
    restarts<<"restart\tslot\tseed\tproposed\taccepted\tuphill\tbest_cost\tbest_step\tfinal_cost\n";
    auto setup_end=Clock::now();State best;Cost best_cost=std::numeric_limits<Cost>::max();uint64_t winning_seed=o.seed,total_moves=0,swap_attempts=0,swap_accepts=0;int winner_restart=-1,winner_slot=-1;uint64_t winner_step=0;
    if(o.verify_only){best=read_state(o,fs);}
    else {
        Workers pool(o.replicas);RNG exchange(mix(o.seed^0x729acd));
        std::vector<double> temps(o.replicas);for(int r=0;r<o.replicas;++r)temps[r]=o.replicas==1?o.tmin:o.tmin*std::pow(o.tmax/o.tmin,double(r)/(o.replicas-1));
        for(int restart=0;restart<o.restarts*o.seed_tries;++restart) {
            std::vector<std::unique_ptr<Replica>> rs;for(int r=0;r<o.replicas;++r)rs.emplace_back(new Replica(mix(o.seed+uint64_t(restart)*o.replicas+r),temps[r],r));
            pool.run([&](size_t r){rs[r]->init(g,fs,o.diagonal);});
            auto collect=[&](){for(int r=0;r<o.replicas;++r)if(rs[r]->best_cost<best_cost){best_cost=rs[r]->best_cost;best=rs[r]->best;winning_seed=rs[r]->seed;winner_restart=restart;winner_slot=r;winner_step=rs[r]->best_step;}};
            auto log=[&](int step){for(int r=0;r<o.replicas;++r){auto& v=*rs[r];conv<<restart<<'\t'<<step<<'\t'<<r<<'\t'<<v.temperature<<'\t'<<v.walker<<'\t'<<v.cache.total<<'\t'<<best_cost<<'\t'<<v.proposed<<'\t'<<v.accepted<<'\t'<<v.uphill<<'\n';}};
            collect();log(0);int step=0,round=0;
            while(step<o.steps && (best_cost!=0||o.keep_going)) {
                int chunk=std::min(o.swap_interval-step%o.swap_interval,o.steps-step);
                chunk=std::min(chunk,o.log_interval-step%o.log_interval);
                pool.run([&](size_t r){for(int j=0;j<chunk;++j){if(!o.keep_going && rs[r]->best_cost==0)break;rs[r]->step(g,fs,o.guide);}
                    if(o.debug){Cache check;check.rebuild(rs[r]->state,g,fs);require(check.total==rs[r]->cache.total&&check.ranks==rs[r]->cache.ranks&&check.fail==rs[r]->cache.fail&&check.fail_sum==rs[r]->cache.fail_sum,"Incremental cache mismatch");}
                });step+=chunk;collect();
                if(step%o.swap_interval==0 && (best_cost!=0||o.keep_going)) {
                    int first=o.replicas==2?0:round%2;
                    for(int i=first;i+1<o.replicas;i+=2){auto& a=*rs[i];auto& b=*rs[i+1];double lp=swap_log(a.temperature,b.temperature,a.cache.total,b.cache.total);bool yes=accept_log(lp,exchange);++swap_attempts;swap_accepts+=yes;
                        swaps<<restart<<'\t'<<round<<'\t'<<step<<'\t'<<i<<'\t'<<i+1<<'\t'<<a.cache.total<<'\t'<<b.cache.total<<'\t'<<lp<<'\t'<<yes<<'\t'<<a.walker<<'\t'<<b.walker<<'\n';
                        if(yes){std::swap(a.state,b.state);std::swap(a.cache,b.cache);std::swap(a.walker,b.walker);a.remember();b.remember();}
                    }++round;
                }
                if(step%o.log_interval==0||step==o.steps||(!o.keep_going&&best_cost==0)){log(step);std::cout<<"restart="<<restart<<" step="<<step<<" best_cost="<<best_cost<<" swaps="<<swap_accepts<<'/'<<swap_attempts<<'\n';}
            }
            // Include success at initialization and during the last restart.
            for(int r=0;r<o.replicas;++r){auto& v=*rs[r];total_moves+=v.proposed;restarts<<restart<<'\t'<<r<<'\t'<<v.seed<<'\t'<<v.proposed<<'\t'<<v.accepted<<'\t'<<v.uphill<<'\t'<<v.best_cost<<'\t'<<v.best_step<<'\t'<<v.cache.total<<'\n';}
            if(best_cost==0&&!o.keep_going)break;
        }
    }
    conv.close();swaps.close();restarts.close();auto search_end=Clock::now();
    if(o.verify_only)setup_end=search_end; // Input loading belongs to setup; no search occurred.
    Cache verified;verified.rebuild(best,g,fs);
    require(o.verify_only||verified.total==best_cost,"Final verification disagrees with search; no certificate issued");
    size_t failing=std::count_if(verified.costs.begin(),verified.costs.end(),[](int c){return c!=0;});bool ame=failing==0;
    auto verify_end=Clock::now();export_data(o,g,fs,best,verified,dir);auto export_end=Clock::now();
    auto seconds=[](auto a,auto b){return std::chrono::duration<double>(b-a).count();};
    Timings t{seconds(start,setup_end),seconds(setup_end,search_end),seconds(search_end,verify_end),seconds(verify_end,export_end),seconds(start,export_end)};
    auto cert=output_file(dir/"certificate.json");
    cert<<"{\n\"schema\":1,\"version\":\""<<VERSION<<"\",\"N\":"<<o.n<<",\"d\":"<<o.d<<",\"ame\":"<<(ame?"true":"false")<<",\"verified_cost\":"<<verified.total<<",\"cuts\":"<<g.cuts.size()<<",\"failing_cuts\":"<<failing<<",\n";
    cert<<"\"phase_convention\":\"tensor of fields; exp(2*pi*i*Tr(sum_{i<=j} Pij*xi*xj)/p); labels in polynomial basis\",\n\"fields\":[";
    for(size_t f=0;f<fs.size();++f){if(f)cert<<',';cert<<"{\"p\":"<<fs[f].p<<",\"m\":"<<fs[f].m<<",\"q\":"<<fs[f].q<<",\"modulus\":[";for(size_t j=0;j<fs[f].modulus.size();++j){if(j)cert<<',';cert<<fs[f].modulus[j];}cert<<"],\"matrix\":";matrix_json(cert,best[f],o.n);cert<<'}';}
    cert<<"],\n\"parameters\":{\"seed\":"<<o.seed<<",\"replicas\":"<<o.replicas<<",\"steps\":"<<o.steps<<",\"restarts\":"<<o.restarts<<",\"seed_tries\":"<<o.seed_tries<<",\"swap_interval\":"<<o.swap_interval<<",\"log_interval\":"<<o.log_interval<<",\"tmin\":"<<o.tmin<<",\"tmax\":"<<o.tmax<<",\"guide\":"<<o.guide<<",\"diagonal\":"<<(o.diagonal?"true":"false")<<",\"keep_going\":"<<(o.keep_going?"true":"false")<<",\"verify_only\":"<<(o.verify_only?"true":"false")<<"},\n";
    cert<<"\"winning_slot_seed\":"<<winning_seed<<",\"winning_restart\":"<<winner_restart<<",\"winning_slot\":"<<winner_slot<<",\"winning_step\":"<<winner_step<<",\"proposals\":"<<total_moves<<",\"swap_attempts\":"<<swap_attempts<<",\"swap_accepts\":"<<swap_accepts<<",\n";
    cert<<"\"timings_seconds\":{\"setup\":"<<t.setup<<",\"search\":"<<t.search<<",\"verification\":"<<t.verification<<",\"data_export\":"<<t.export_data<<",\"total_through_data_export\":"<<t.total<<"}\n}\n";cert.close();
    auto summary=output_file(dir/"run_summary.txt");summary<<VERSION<<"\n"<<(ame?"CONFIRMED AME":"NOT AME: search best only; no optimality or nonexistence claim")<<"\nN="<<o.n<<" d="<<o.d<<" cuts="<<g.cuts.size()<<" failing="<<failing<<" cost="<<verified.total<<"\nseed="<<o.seed<<" proposals="<<total_moves<<" swap_attempts="<<swap_attempts<<" swap_accepts="<<swap_accepts<<"\nsetup_seconds="<<t.setup<<"\nsearch_seconds="<<t.search<<"\nverification_seconds="<<t.verification<<"\ndata_export_seconds="<<t.export_data<<"\ntotal_through_data_export_seconds="<<t.total<<"\nTiming excludes writing this summary and certificate metadata.\n";
    std::cout<<(ame?"CONFIRMED AME":"NOT AME")<<" failing="<<failing<<'/'<<g.cuts.size()<<" total_seconds="<<t.total<<"\nOutput: "<<dir.string()<<'\n';
    if(o.print_entropy)std::cout<<"Per-cut and per-size entropies: entropy_bips.tsv and entropy_per_k.tsv\n";
    return ame?0:2;
}
} // namespace ame
#ifndef AME_NO_MAIN
int main(int argc,char**argv) {
    try {
        for(int i=1;i<argc;++i)if(std::string(argv[i])=="--help"||std::string(argv[i])=="-h") {
            std::cout<<"AfAME - Developed by Dr. Zakaria Dahbi\n"
             "AME7DIA_PT (C++17 CPU threads)\n"
             "-N 7 -d 2 --replicas 4 --steps 5000 --restarts 20 --seed 1\n"
             "--tmin 0.25 --tmax 50 --guide 0.85 --swap-interval 25 --log-interval 200\n"
             "--output NEW_DIRECTORY  --memory-mb 512  --seed-tries 1\n"
             "--input MATRIX.txt --verify-only   --debug-checks   --keep-going\n"
             "--diagonal: random local diagonal phases, fixed during search\n"
             "Legacy -tensor_phases [1] is accepted; dimension factorization is automatic.\n"
             "Exit codes: 0 certified AME; 2 verified non-AME best; 1 error.\n"
             "Fields q<=1024. N<=24 and a preallocation memory check.\n";return 0;
        }
        return ame::run(ame::parse(argc,argv));
    }catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}
}
#endif
