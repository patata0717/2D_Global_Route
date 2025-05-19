// PJ4 .hpp
#ifndef DATASTRUCTURE_HPP
#define DATASTRUCTURE_HPP
#include <vector>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <iterator>
#include <limits.h>
#include <tuple>
#include <random>
#include <numeric>
using namespace std;

extern int M;
extern int N;


class GridMap;
class Net;
struct Coord;


struct Coord {
    int x, y;
    bool operator==(const Coord& o) const { return x==o.x && y==o.y; }
    bool operator<(const Coord& o) const  {                 // NEW
        return x < o.x || (x == o.x && y < o.y);
    }
};


struct hash_pair {
    size_t operator()(const Coord& c) const {
        return hash<int>()(c.x) ^ (hash<int>()(c.y) << 1);
    }
};

class GridMap {
public:
    vector<vector<bool>> horizontal;
    vector<vector<bool>> vertical;

    GridMap() : horizontal(M, vector<bool>(N - 1, false)), vertical(M - 1, vector<bool>(N, false)) {}
};

class Net {
public:
    Net() : id(-1) {}
    Net(int id) : id(id) {}

    int id;
    vector<Coord> Gcells;
    GridMap best_steiner;


    unordered_set<Coord, hash_pair> steiner_nodes;
    unordered_map<Coord, vector<Coord>, hash_pair> adjacency;
};

struct DSU {
    std::vector<int> p, sz;
    DSU(size_t n = 0){ reset(n); }
    void reset(size_t n){
        p.resize(n); sz.assign(n,1);
        std::iota(p.begin(), p.end(), 0);
    }
    int find(int x){ return p[x]==x?x:p[x]=find(p[x]); }
    void unite(int a,int b){
        a=find(a); b=find(b);
        if(a==b) return;
        if(sz[a]<sz[b]) std::swap(a,b);
        p[b]=a; sz[a]+=sz[b];
    }
    int num_sets() const{
        int c=0; for(size_t i=0;i<p.size();++i) if(p[i]==(int)i) ++c;
        return c;
    }
    int first_set() const{
        for(size_t i=0;i<p.size();++i) if(p[i]==(int)i) return i;
        return -1;
    }
};

#endif