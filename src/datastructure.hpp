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
using namespace std;

extern int M;
extern int N;


class GridMap;
class Net;
struct Coord;


struct Coord {
    int x;
    int y;

    bool operator==(const Coord& other) const {
        return x == other.x && y == other.y;
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
    GridMap zigzag;
    GridMap best_steiner;

    unordered_set<Coord, hash_pair> steiner_nodes;
    unordered_map<Coord, vector<Coord>, hash_pair> adjacency;
};

#endif