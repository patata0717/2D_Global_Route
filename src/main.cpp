// PJ4 main
#include <iostream>
#include <set>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <stack>
#include <array>
#include <queue>
#include <cstdlib>
#include <ctime>
#include <sys/time.h>
#include <time.h>
#include <cmath>
#include <random>
#include <algorithm>
#include "datastructure.hpp"
#include <chrono>
using namespace std;

enum Direction {
    UP,
    RIGHT,
    DOWN,
    LEFT
};

// Global variable
int M, N;
int H_limit, V_limit;

double GetTime(void);
bool GetLine(istream &input_file, stringstream &line);
void Zigzag_route(Net& n);
void Best_steiner_route(Net& n);
void Best_steiner_obstacle_route(GridMap obstacle, Net n);
void Rip_up_reroute(GridMap gridmap, int seed);
Coord Search_nearest(Coord src, Direction dir, const unordered_set<Coord, hash_pair>& steiner_nodes);
vector<pair<Coord, Coord>> find_cycle(const unordered_map<Coord, vector<Coord>, hash_pair>& adjacency);
void MapAdjacencyToGridMap(Net &net);
void PrintGridMap(const GridMap &grid);

// void Tree2Route(adjlist steinertree, GridMap gridmap) {

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <input_file> <output_file> <from> <to>\n";
        return 1;
    }

    double t0;
    t0 = GetTime(); // start time

/*---Read input---*/
    string input_filepath = argv[1];  // ../testcase/sample.txt
    ifstream input_file(input_filepath);
    stringstream line;

    // Read Number of Hardblocks
    // Gcell_grid 4 5
    if (GetLine(input_file, line)) {
        string dummy_Gcell_grid;
        line >> dummy_Gcell_grid >> M >> N;
    }
    cout << "M: " << M << ", N: " << N << endl;
    // HTrack 4
    if (GetLine(input_file, line)) {
        string HTrack;
        line >> HTrack >> H_limit;
    }
    // VTrack 3
    if (GetLine(input_file, line)) {
        string VTrack;
        line >> VTrack >> V_limit;
    }

    // nets
    // N1  (0 0) (1 0) (2 2)
    // ...
    // end_of_nets
    vector<Net> nets(16 + 1);
    while (GetLine(input_file, line)) {
        string keyword;
        line >> keyword;
        if (keyword == "nets") break;
    }
    while (GetLine(input_file, line)) {
        string first_word;
        line >> first_word;
        if (first_word == "end_of_nets") break;

        int id = stoi(first_word.substr(1));
        nets[id] = Net(id);

        int x, y;
        char Lbracket, Rbracket;
        while (line >> Lbracket >> x >> y >> Rbracket) {
            nets[id].Gcells.push_back(Coord{x, y});
        }
    }

    // Print M N
    // Print HTrack and VTrack limits
    cout << "M: " << M << ", N: " << N << endl;
    cout << "HTrack limit: " << H_limit << ", VTrack limit: " << V_limit << endl;

    // Optional: Print for verification
    for (const auto& net : nets) {
        cout << "Net " << net.id << ": ";
        for (const auto& c : net.Gcells) {
            cout << "(" << c.x << "," << c.y << ") ";
        }
        cout << endl;
    }
        
    input_file.close();

/*---Init Route---*/

    cout << M << " " << N << endl;
    vector<vector<int>> super_horizontal(M, vector<int>(N - 1, 0));
    vector<vector<int>> super_vertical(M - 1, vector<int>(N, 0));

    for (int i = stoi(argv[3]); i <= stoi(argv[4]); i++) {
        cout << "Net " << i << ":\n";
        Zigzag_route(nets[i]);
        Best_steiner_route(nets[i]);
        // print adjacency
        cout << "Adjacency List:\n";
        for (const auto& entry : nets[i].adjacency) {
            const Coord& from = entry.first;
            const vector<Coord>& neighbors = entry.second;
            cout << "(" << from.x << "," << from.y << "): ";
            for (const auto& to : neighbors) {
                cout << "(" << to.x << "," << to.y << ") ";
            }
            cout << endl;
        }
        MapAdjacencyToGridMap(nets[i]);
        // cout << "\nMapped Grid for Net " << nets[i].id << ":\n";
        // Accumulate usage count
        for (int y = 0; y < M; ++y)
            for (int x = 0; x < N - 1; ++x)
                if (nets[i].best_steiner.horizontal[y][x])
                    ++super_horizontal[y][x];
        for (int y = 0; y < M - 1; ++y)
            for (int x = 0; x < N; ++x)
                if (nets[i].best_steiner.vertical[y][x])
                    ++super_vertical[y][x];
    }

    // Print superimposed net count
    cout << "\nSuperimposed GridMap (Edge Usage Count):\n";
    cout << "Horizontal Matrix:\n";
    for (int y = 0; y < M; ++y) {
        for (int x = 0; x < N - 1; ++x) {
            cout << super_horizontal[y][x] << " ";
        }
        cout << "\n";
    }
    cout << "\n";
    cout << "Vertical Matrix:\n";
    for (int y = 0; y < M - 1; ++y) {
        for (int x = 0; x < N; ++x) {
            cout << super_vertical[y][x] << " ";
        }
        cout << "\n";
    }




/*---Write result to output file*/
    string output_filepath = argv[2];  // ../output/sample.out
    cout << "Writing result to " << output_filepath << endl;
    ofstream output_file(output_filepath);
    output_file << "M " << M << endl;
    output_file << "N " << N << endl;
    output_file << endl;
    output_file << "Horizontal Matrix:\n";
    for (int y = 0; y < M; ++y) {
        for (int x = 0; x < N - 1; ++x) {
            output_file << super_horizontal[y][x] << " ";
        }
        output_file << "\n";
    }
    output_file << "\n";

    output_file << "Vertical Matrix:\n";
    for (int y = 0; y < M - 1; ++y) {
        for (int x = 0; x < N; ++x) {
            output_file << super_vertical[y][x] << " ";
        }
        output_file << "\n";
    }
    output_file.close();
    cout << "Time: " << GetTime() - t0 << endl;

    return 0;
}


double GetTime(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}

bool GetLine(istream &input_file, stringstream &line) {
    string input;
    while (getline(input_file, input)) {
        if (input[0] != '#') {
            line.clear();
            line.str(input);
            return true;  // Found a valid (non-comment) line
        }
    }
    return false;  // No more valid lines
}

void Zigzag_route(Net& n) {
    n.adjacency.clear();

    // Group by x-coordinate (column-wise)
    unordered_map<int, vector<Coord>> cols;
    for (size_t i = 0; i < n.Gcells.size(); ++i) {
        const Coord& c = n.Gcells[i];
        cols[c.x].push_back(c);
    }

    vector<Coord> top_nodes;

    // Vertical connections (within columns)
    for (unordered_map<int, vector<Coord>>::iterator it = cols.begin(); it != cols.end(); ++it) {
        vector<Coord>& col = it->second;


        sort(col.begin(), col.end(), [](const Coord& a, const Coord& b) {
            return a.y < b.y;
        });

        // Always add the top-most node
        top_nodes.push_back(col.back());

        if (col.size() <= 1)
            continue;

        for (size_t i = 1; i < col.size(); ++i) {
            Coord a = col[i - 1];
            Coord b = col[i];
            n.adjacency[a].push_back(b);
            n.adjacency[b].push_back(a);
        }
    }

    // Horizontal connection across top-most nodes
    sort(top_nodes.begin(), top_nodes.end(), [](const Coord& a, const Coord& b) {
        return a.x < b.x;
    });

    for (size_t i = 1; i < top_nodes.size(); ++i) {
        Coord a = top_nodes[i - 1];
        Coord b = top_nodes[i];

        n.adjacency[a].push_back(b);
        n.adjacency[b].push_back(a);
    }
}


inline int manhattan(const Coord &a, const Coord &b) {
    return std::abs(a.x - b.x) + std::abs(a.y - b.y);
}

static void add_edge(std::unordered_map<Coord, std::vector<Coord>, hash_pair> &adj,
                     const Coord &u, const Coord &v) {
    adj[u].push_back(v);
    adj[v].push_back(u);
}

static void erase_edge(std::unordered_map<Coord, std::vector<Coord>, hash_pair> &adj,
                       const Coord &u, const Coord &v) {
    std::vector<Coord> &au = adj[u];
    au.erase(std::remove(au.begin(), au.end(), v), au.end());
    std::vector<Coord> &av = adj[v];
    av.erase(std::remove(av.begin(), av.end(), u), av.end());
}

// ───────── algorithm ───────────────────────────────────────

void Best_steiner_route(Net &net) {
    // ─── Hard‑reset any pre‑existing Steiner points so each run
    //      starts with *only* the true terminals. ────────────────
    int   best_gain  = std::numeric_limits<int>::min();
    Coord best_coord = { -1, -1 };

    const std::array<Direction, 4> dirs = { UP, RIGHT, DOWN, LEFT };

    // convert Gcells once for O(1) lookup
    std::unordered_set<Coord, hash_pair> fixed(net.Gcells.begin(), net.Gcells.end());

    bool do_flag;
    do {
        for (int y = 0; y < M; ++y) {
            for (int x = 0; x < N; ++x) {
                Coord cur{ x, y };
                if (fixed.count(cur) || net.steiner_nodes.count(cur))
                    continue;

//                std::cout << "\n[Evaluating] Coord: (" << cur.x << "," << cur.y << ")\n";

                // local working copies
                std::unordered_map<Coord, std::vector<Coord>, hash_pair> local_adj = net.adjacency;
                std::unordered_set<Coord, hash_pair> local_nodes = fixed;
                local_nodes.insert(net.steiner_nodes.begin(), net.steiner_nodes.end());
                local_nodes.insert(cur);

                int total_gain = 0;

                for (size_t k = 0; k < dirs.size(); ++k) {
                    Direction dir = dirs[k];
//                    std::cout << "  Direction " << static_cast<int>(dir) << " → ";

                    //–– robust neighbour search ––
                    Coord nei = Search_nearest(cur, dir, local_nodes);
                    if (nei.x == -1 || !local_nodes.count(nei)) {
//                        std::cout << "No Steiner neighbor found." << std::endl;
                        continue;
                    }

//                    std::cout << "Found neighbor: (" << nei.x << "," << nei.y << ")" << std::endl;

                    add_edge(local_adj, cur, nei);
                    int dist = manhattan(cur, nei);
                    total_gain -= dist;
//                    std::cout << "    Added edge: (" << cur.x << "," << cur.y << ") - ("
//                              << nei.x << "," << nei.y << "), dist = " << dist
//                              << ", gain -= " << dist << std::endl;

                    std::vector<std::pair<Coord, Coord> > cycle_edges = find_cycle(local_adj);
                    if (!cycle_edges.empty()) {
                        std::pair<Coord, Coord> largest = *std::max_element(
                            cycle_edges.begin(), cycle_edges.end(),
                            [](const std::pair<Coord, Coord> &a,
                               const std::pair<Coord, Coord> &b) {
                                return manhattan(a.first, a.second) < manhattan(b.first, b.second);
                            });
                        int longest = manhattan(largest.first, largest.second);
                        erase_edge(local_adj, largest.first, largest.second);
                        total_gain += longest;
//                        std::cout << "    🔁 Loop detected! Removing largest edge: ("
//                                  << largest.first.x << "," << largest.first.y << ") - ("
//                                  << largest.second.x << "," << largest.second.y << "), gain += "
//                                  << longest << std::endl;
                    }

//                    std::cout << "    → Intermediate total gain: " << total_gain << std::endl;
                }

//                std::cout << "Total gain for (" << cur.x << "," << cur.y << ") = "
//                          << total_gain << std::endl;

                if (total_gain > best_gain) {
                    best_gain  = total_gain;
                    best_coord = cur;
                }
            }
        }

        do_flag = false;
        if (best_gain >= 0) {
            do_flag = true;
            // ─── commit winner ────────────────────────────────────────────
            net.steiner_nodes.insert(best_coord);

            // Candidate set for neighbour search = pins ∪ current Steiner nodes
            std::unordered_set<Coord, hash_pair> connect_set = fixed;
            connect_set.insert(net.steiner_nodes.begin(), net.steiner_nodes.end());

            for (size_t k = 0; k < dirs.size(); ++k) {
                Direction dir = dirs[k];
                Coord nei = Search_nearest(best_coord, dir, connect_set);
                if (nei.x == -1) continue;

                add_edge(net.adjacency, best_coord, nei);
                std::vector<std::pair<Coord, Coord> > cycle_edges = find_cycle(net.adjacency);
                if (!cycle_edges.empty()) {
                    std::pair<Coord, Coord> largest = *std::max_element(
                        cycle_edges.begin(), cycle_edges.end(),
                        [](const std::pair<Coord, Coord> &a,
                           const std::pair<Coord, Coord> &b) {
                            return manhattan(a.first, a.second) < manhattan(b.first, b.second);
                        });
                    erase_edge(net.adjacency, largest.first, largest.second);
                }
            }
        }

//        std::cout << "\nBest gain = " << best_gain
//                  << " at (" << best_coord.x << ", " << best_coord.y << ")" << std::endl;
        // Reset best gain
        best_gain = std::numeric_limits<int>::min();
        best_coord = { -1, -1 };
    } while (do_flag);
}

void Best_steiner_obstacle_route(GridMap obstacle, Net n) {
    GridMap best_steiner;
    // Implement best Steiner tree logic with obstacles here
}

void Rip_up_reroute(GridMap gridmap, int seed) {
    // Implement rip up and reroute logic here
}

Coord Search_nearest(Coord src, Direction dir, const unordered_set<Coord, hash_pair>& steiner_nodes) {
    Coord invalid = {-1, -1};
    unordered_set<Coord, hash_pair> visited;
    queue<Coord> q;

    auto in_bounds = [&](Coord c) {
        return 0 <= c.x && c.x < N && 0 <= c.y && c.y < M;
    };

    auto in_cone = [&](Coord from, Coord to) -> bool {
        int dx = to.x - from.x;
        int dy = to.y - from.y;

        switch (dir) {
            case UP:
                return dy > 0 && dx > -dy && dx <= dy;
            case DOWN:
                return dy < 0 && dx >= dy && dx < -dy;
            case LEFT:
                return dx < 0 && dy > dx && dy <= -dx;
            case RIGHT:
                return dx > 0 && dy >= -dx && dy < dx;
        }
        return false;
    };

    vector<Coord> directions = {
        {0, 1},   // UP
        {-1, 0},  // LEFT
        {0, -1},  // DOWN
        {1, 0}    // RIGHT
    };

    visited.insert(src);
    q.push(src);


    while (!q.empty()) {
        Coord cur = q.front(); q.pop();

        for (const auto& d : directions) {
            Coord next = {cur.x + d.x, cur.y + d.y};
            // Out of M by N or visited
            if (!in_bounds(next) || visited.count(next)) continue;
            // Not in <dir> region
            if (!in_cone(src, next)) continue;
            // Not match
            if (steiner_nodes.count(next)) return next;

            visited.insert(next);
            q.push(next);
        }
    }

    return invalid;  // not found
}

std::vector<std::pair<Coord, Coord> > find_cycle(const std::unordered_map<Coord, std::vector<Coord>, hash_pair> &adj)
{
    std::unordered_map<Coord, Coord, hash_pair> parent;
    std::unordered_set<Coord, hash_pair>        visited;
    std::stack<Coord>                           st;

    for (const auto &kv : adj) {
        Coord start = kv.first;
        if (visited.count(start)) continue;

        parent[start] = Coord{ -1, -1 };
        st.push(start);

        while (!st.empty()) {
            Coord u = st.top(); st.pop();
            visited.insert(u);

            for (const Coord &v : adj.at(u)) {
                if (v == parent[u]) continue; // skip the tree‑edge back to parent

                if (!visited.count(v)) {
                    parent[v] = u;
                    st.push(v);
                } else {
                    // ── Found a back‑edge u‑v → cycle. Collect *all* edges. ──
                    std::vector<std::pair<Coord, Coord> > path_u, path_v;
                    Coord a = u, b = v;

                    // climb from a to root, store edges
                    while (a.x != -1) {
                        Coord p = parent[a];
                        if (p.x == -1) break;
                        path_u.push_back(std::make_pair(a, p));
                        if (p == v) break; // met other side early
                        a = p;
                    }
                    // climb from b to root, store edges
                    while (b.x != -1) {
                        Coord p = parent[b];
                        if (p.x == -1) break;
                        path_v.push_back(std::make_pair(b, p));
                        if (p == u) break;
                        b = p;
                    }
                    // The cycle consists of u‑v plus both paths
                    std::vector<std::pair<Coord, Coord> > cycle_edges;
                    cycle_edges.push_back(std::make_pair(u, v));
                    cycle_edges.insert(cycle_edges.end(), path_u.begin(), path_u.end());
                    cycle_edges.insert(cycle_edges.end(), path_v.begin(), path_v.end());

                    // Remove duplicate edge directions (v,u) vs (u,v)
                    std::unordered_set<long long> seen;
                    std::vector<std::pair<Coord, Coord> > unique_edges;
                    for (auto &e : cycle_edges) {
                        long long key1 = (static_cast<long long>(e.first.x) << 48) ^ (static_cast<long long>(e.first.y) << 32) ^ (static_cast<long long>(e.second.x) << 16) ^ e.second.y;
                        long long key2 = (static_cast<long long>(e.second.x) << 48) ^ (static_cast<long long>(e.second.y) << 32) ^ (static_cast<long long>(e.first.x) << 16) ^ e.first.y;
                        if (!seen.count(key1) && !seen.count(key2)) {
                            seen.insert(key1);
                            unique_edges.push_back(e);
                        }
                    }
                    return unique_edges;
                }
            }
        }
    }
    return {}; // no cycle
}

void MapAdjacencyToGridMap(Net &net) {
    net.best_steiner = GridMap(); // Reset the grid

    for (const auto &entry : net.adjacency) {
        const Coord &from = entry.first;
        const vector<Coord> &neighbors = entry.second;

        for (const Coord &to : neighbors) {
            int dx = to.x - from.x;
            int dy = to.y - from.y;
            // Only consider unit Manhattan distance (adjacent)
            if (abs(dx) + abs(dy) != 1)
                continue;

            if (dx == 1) {
                net.best_steiner.horizontal[from.y][from.x] = true;
            } else if (dx == -1) {
                net.best_steiner.horizontal[from.y][to.x] = true;
            } else if (dy == 1) {
                net.best_steiner.vertical[from.y][from.x] = true;
            } else if (dy == -1) {
                net.best_steiner.vertical[to.y][from.x] = true;
            }
        }
    }
}


void PrintGridMap(const GridMap &grid) {
    cout << "Horizontal Matrix:\n";
    for (int y = 0; y < M; ++y) {
        for (int x = 0; x < N - 1; ++x) {
            cout << grid.horizontal[y][x] << " ";
        }
        cout << "\n";
    }

    cout << "Vertical Matrix:\n";
    for (int y = 0; y < M - 1; ++y) {
        for (int x = 0; x < N; ++x) {
            cout << grid.vertical[y][x] << " ";
        }
        cout << "\n";
    }
}
