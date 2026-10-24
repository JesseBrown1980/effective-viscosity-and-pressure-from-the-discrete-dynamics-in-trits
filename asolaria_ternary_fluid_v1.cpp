// Asolaria-informed discrete-fluid experiment (NOT a reproduction of ASI-OS).
// Integer-only, 3D, 27 balanced-ternary velocity directions; each velocity
// population is an occupancy trit (0, 1 or 2). No floating point.
// Local collisions are deterministic involutions within exact mass/momentum/
// kinetic-energy classes. Streaming is a permutation of site/channel slots.
// This is a lattice gas, NOT a molecular validation or proof of Navier-Stokes.
#include <array>
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <vector>

using std::array;
using i64 = std::int64_t;
using u64 = std::uint64_t;
constexpr int NX=24, NY=18, NZ=12, K=27, CELLS=NX*NY*NZ;
constexpr int FINAL=288;
using State = std::vector<std::uint8_t>;
struct Velocity {int x,y,z,e;};
struct Totals {i64 m=0, px=0, py=0, pz=0, energy2=0;};
struct RowStats {i64 contrast=0, roughness=0;};
struct Collision {
    int a,b,c,d;
    array<std::uint8_t,81> lut{};
};

static array<Velocity,K> velocities() {
    array<Velocity,K> out{};
    int k=0;
    for (int z=-1;z<=1;++z) for (int y=-1;y<=1;++y)
      for (int x=-1;x<=1;++x) out[k++]={x,y,z,x*x+y*y+z*z};
    return out;
}
static int index(int x,int y,int z) {return (z*NY+y)*NX+x;}
static int mod(int a,int n) {int r=a%n;return r<0?r+n:r;}
static u64 next64(u64& seed) {
    seed+=0x9e3779b97f4a7c15ull;
    u64 z=seed;
    z=(z^(z>>30))*0xbf58476d1ce4e5b9ull;
    z=(z^(z>>27))*0x94d049bb133111ebull;
    return z^(z>>31);
}
static int encode(int a,int b,int c,int d){return a+3*b+9*c+27*d;}
static array<int,4> decode(int code) {
    array<int,4> out{};
    for (int i=0;i<4;++i) {out[i]=code%3;code/=3;}
    return out;
}
static std::vector<Collision> make_collisions(const array<Velocity,K>& v) {
    // All distinct four-channel configurations admitting a physically exact
    // 2->2 collision with equal vector momentum and equal squared speed sum.
    using Key=std::tuple<int,int,int,int>;
    std::map<Key,std::vector<array<int,2>>> pairs;
    for(int a=0;a<K;++a) for(int b=a+1;b<K;++b) {
        pairs[{v[a].x+v[b].x,v[a].y+v[b].y,
               v[a].z+v[b].z,v[a].e+v[b].e}].push_back({a,b});
    }
    std::set<array<int,4>> quads;
    for(const auto& entry:pairs) {
        const auto& p=entry.second;
        for(std::size_t i=0;i<p.size();++i)
          for(std::size_t j=i+1;j<p.size();++j) {
              auto q=array<int,4>{p[i][0],p[i][1],p[j][0],p[j][1]};
              std::sort(q.begin(),q.end());
              if(std::unique(q.begin(),q.end())==q.end()) quads.insert(q);
          }
    }
    std::vector<Collision> out;
    for(const auto& q:quads) {
        Collision c{q[0],q[1],q[2],q[3],{}};
        std::map<std::array<int,5>,std::vector<int>> cls;
        for(int mask=0;mask<81;++mask) {
            auto d=decode(mask);
            int m=0,px=0,py=0,pz=0,e=0;
            for(int j=0;j<4;++j) {
                auto vel=v[q[j]];
                m+=d[j];px+=d[j]*vel.x;py+=d[j]*vel.y;
                pz+=d[j]*vel.z;e+=d[j]*vel.e;
            }
            cls[{m,px,py,pz,e}].push_back(mask);
            c.lut[mask]=static_cast<std::uint8_t>(mask);
        }
        for(const auto& entry:cls) {
            const auto& codes=entry.second;
            for(std::size_t j=0;j+1<codes.size();j+=2) {
                c.lut[codes[j]]=static_cast<std::uint8_t>(codes[j+1]);
                c.lut[codes[j+1]]=static_cast<std::uint8_t>(codes[j]);
            }
        }
        int swaps=0;
        for(int mask=0;mask<81;++mask) {
            int dest=c.lut[mask];
            if(c.lut[dest]!=mask)throw std::runtime_error("collision not involutive");
            auto invariant=[&](int code){
                auto d=decode(code);
                int m=0,px=0,py=0,pz=0,e=0;
                for(int j=0;j<4;++j){
                   m+=d[j];px+=d[j]*v[q[j]].x;py+=d[j]*v[q[j]].y;
                   pz+=d[j]*v[q[j]].z;e+=d[j]*v[q[j]].e;
                }
                return std::array<int,5>{m,px,py,pz,e};
            };
            if(invariant(mask)!=invariant(dest))throw std::runtime_error("bad conservation LUT");
            swaps+=(dest!=mask);
        }
        if(swaps)out.push_back(c);
    }
    return out;
}
static State initial(const array<Velocity,K>& v) {
    State f(CELLS*K,0);
    u64 s=0xA501A21B00000003ull;
    for (int z=0;z<NZ;++z) for(int y=0;y<NY;++y) for(int x=0;x<NX;++x) {
        int dir=y<NY/2?1:-1;
        for (int k=0;k<K;++k) {
            int r=static_cast<int>(next64(s)%12);
            // Bias only the x-direction, leaving y and z motion populated.
            int occup;
            if (v[k].x*dir==1) occup=r<2?0:(r<6?1:2);
            else if(v[k].x*dir==-1) occup=r<6?0:(r<10?1:2);
            else occup=r<4?0:(r<8?1:2);
            f[index(x,y,z)*K+k]=static_cast<std::uint8_t>(occup);
        }
    }
    return f;
}
static Totals totals(const State& f,const array<Velocity,K>& v){
    Totals t;
    for(int i=0;i<CELLS;++i) for(int k=0;k<K;++k){
        int n=f[i*K+k];
        t.m+=n;t.px+=n*v[k].x;t.py+=n*v[k].y;t.pz+=n*v[k].z;
        t.energy2+=n*v[k].e;
    }
    return t;
}
static RowStats row_stats(const State& f,const array<Velocity,K>& v) {
    i64 rows[NY]={};
    i64 global=0;
    for(int z=0;z<NZ;++z)for(int y=0;y<NY;++y)for(int x=0;x<NX;++x){
        int cell=index(x,y,z);
        for(int k=0;k<K;++k)rows[y]+=f[cell*K+k]*v[k].x;
    }
    RowStats st;
    for(int y=0;y<NY;++y){
        global+=rows[y];
        st.contrast+=(y<NY/2?1:-1)*rows[y];
    }
    for(int y=0;y<NY;++y){
        i64 delta=rows[y]*NY-global;
        st.roughness+=delta<0?-delta:delta;
    }
    return st;
}
static void apply_collisions(State& f,const std::vector<Collision>& c,int step){
    // A stack of local microreversible ternary collision gates.
    const int count=32;
    const std::size_t total=c.size();
    for(int cell=0;cell<CELLS;++cell) {
        int base=cell*K;
        for(int j=0;j<count;++j) {
            std::size_t which=(static_cast<std::size_t>(step)*count+
                static_cast<std::size_t>(cell)*17+static_cast<std::size_t>(j)*37)%total;
            const auto& col=c[which];
            int code=encode(f[base+col.a],f[base+col.b],f[base+col.c],f[base+col.d]);
            auto dest=decode(col.lut[code]);
            f[base+col.a]=static_cast<std::uint8_t>(dest[0]);
            f[base+col.b]=static_cast<std::uint8_t>(dest[1]);
            f[base+col.c]=static_cast<std::uint8_t>(dest[2]);
            f[base+col.d]=static_cast<std::uint8_t>(dest[3]);
        }
    }
}
static void stream(const State&,State&,const array<Velocity,K>&);
static void undo_collisions(State& f,const std::vector<Collision>& c,int step){
    // Inverse of the time-varying stack: each gate is an exact involution.
    const int count=32;
    const std::size_t total=c.size();
    for(int cell=0;cell<CELLS;++cell){
        const int base=cell*K;
        for(int j=count-1;j>=0;--j){
            std::size_t which=(static_cast<std::size_t>(step)*count+
                static_cast<std::size_t>(cell)*17+static_cast<std::size_t>(j)*37)%total;
            const auto& col=c[which];
            int code=encode(f[base+col.a],f[base+col.b],f[base+col.c],f[base+col.d]);
            auto dest=decode(col.lut[code]);
            f[base+col.a]=static_cast<std::uint8_t>(dest[0]);
            f[base+col.b]=static_cast<std::uint8_t>(dest[1]);
            f[base+col.c]=static_cast<std::uint8_t>(dest[2]);
            f[base+col.d]=static_cast<std::uint8_t>(dest[3]);
        }
    }
}
static void unstream(const State& after,State& before,const array<Velocity,K>& v){
    for(int z=0;z<NZ;++z)for(int y=0;y<NY;++y)for(int x=0;x<NX;++x){
        int dst=index(x,y,z)*K;
        for(int k=0;k<K;++k){
            int src=index(mod(x+v[k].x,NX),mod(y+v[k].y,NY),
                mod(z+v[k].z,NZ))*K+k;
            before[dst+k]=after[src];
        }
    }
}
static void verify_reverse(const State& start,const array<Velocity,K>& v,
                           const std::vector<Collision>& c){
    const int count=72;
    State cur=start,next(CELLS*K,0);
    for(int t=0;t<count;++t){
        apply_collisions(cur,c,t);
        stream(cur,next,v);
        cur.swap(next);
    }
    if(cur==start)throw std::runtime_error("unexpected complete state recurrence");
    for(int t=count-1;t>=0;--t){
        unstream(cur,next,v);
        cur.swap(next);
        undo_collisions(cur,c,t);
    }
    if(cur!=start)throw std::runtime_error("backward replay was not exact");
    std::cout<<"REVERSE|steps="<<count<<"|every_occupancy_trit_restored=1|pass=1\n";
}
static void stream(const State& current,State& next,const array<Velocity,K>& v) {
    for(int z=0;z<NZ;++z)for(int y=0;y<NY;++y)for(int x=0;x<NX;++x){
        const int src=index(x,y,z)*K;
        for(int k=0;k<K;++k) {
            int dst=index(mod(x+v[k].x,NX),mod(y+v[k].y,NY),mod(z+v[k].z,NZ))*K+k;
            next[dst]=current[src+k];
        }
    }
}
static void simulate(bool collision_on,const State& start,
    const array<Velocity,K>& v,const std::vector<Collision>& c,
    std::ofstream& csv) {
    State cur=start,next(CELLS*K,0);
    const Totals baseline=totals(cur,v);
    const int checkpoints[]={0,6,12,18,24,36,48,60,72,96,120,144,180,216,288};
    int checkpoint=0;
    for(int t=0;t<=FINAL;++t){
        if (checkpoint<15 && t==checkpoints[checkpoint]) {
            const Totals now=totals(cur,v);
            if(now.m!=baseline.m||now.px!=baseline.px||now.py!=baseline.py||
                now.pz!=baseline.pz||now.energy2!=baseline.energy2){
                throw std::runtime_error("global physical conservation law violated");
            }
            const auto s=row_stats(cur,v);
            const bool original=cur==start;
            csv<<(collision_on?"collisions":"stream_only")<<","<<t<<","<<s.contrast
              <<","<<s.roughness<<","<<now.m<<","<<now.px<<","<<now.py<<","<<now.pz
              <<","<<now.energy2<<","<<(original?1:0)<<"\n";
            std::cout<<(collision_on?"COLL":"FREE")<<"|t="<<t
              <<"|shear="<<s.contrast<<"|rough="<<s.roughness<<"|return="<<original<<"\n";
            checkpoint++;
        }
        if(t==FINAL)break;
        if(collision_on)apply_collisions(cur,c,t);
        stream(cur,next,v);
        cur.swap(next);
    }
}
int main(int argc,char**argv){
    try {
        const auto v=velocities();
        const auto c=make_collisions(v);
        State start=initial(v);
        const auto t=totals(start,v);
        const auto s=row_stats(start,v);
        const char* path=(argc>1?argv[1]:"asolaria_d3q27_integer_fluid_results.csv");
        std::ofstream csv(path);
        if(!csv)throw std::runtime_error("cannot create CSV");
        csv<<"mode,step,shear_contrast,row_momentum_roughness,mass,px,py,pz,kinetic_energy_x2,full_state_returned\n";
        std::cout<<"D3Q27|dimensions="<<NX<<"x"<<NY<<"x"<<NZ
            <<"|cells="<<CELLS<<"|occupancy_trits="<<CELLS*K
            <<"|collision_rules="<<c.size()<<"|floating_point=0\n";
        std::cout<<"GLOBAL|mass="<<t.m<<"|px="<<t.px<<"|py="<<t.py
            <<"|pz="<<t.pz<<"|kinetic_energy_x2="<<t.energy2
            <<"|initial_shear="<<s.contrast<<"\n";
        simulate(false,start,v,c,csv);
        simulate(true,start,v,c,csv);
        verify_reverse(start,v,c);
        std::cout<<"RESULT|global_mass_momentum_kinetic_energy_exact_conservation=PASS\n";
        return 0;
    }catch(const std::exception& e){std::cerr<<"FAIL|"<<e.what()<<"\n";return 1;}
}
