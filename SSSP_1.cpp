#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <limits>
#include <set>
#include <map>

using namespace std;

typedef long long ll;
const ll INF = 1e18;

struct Edge {
    int to;
    ll weight;
};

int n, m;
vector<vector<Edge>> adj, radj;
vector<ll> db;
vector<int> pred;
int k, t, max_level;

bool tie_breaker(int u, int v) {
    return u < v;
}

class DataStructureD {
private:
    struct Block {
        vector<pair<int, ll>> elements;
        ll separator_lower;
        ll separator_upper;
        
        Block(ll lower = INF, ll upper = -INF) : separator_lower(lower), separator_upper(upper) {}
        
        bool empty() const { return elements.empty(); }
        size_t size() const { return elements.size(); }
        
        void add(int vertex, ll dist) {
            elements.push_back({vertex, dist});
        }
        
        ll min_value() const { 
            if (empty()) return INF;
            ll min_val = INF;
            for (const auto& [v, d] : elements) {
                min_val = min(min_val, d);
            }
            return min_val;
        }
        
        ll max_value() const { 
            if (empty()) return -INF;
            ll max_val = -INF;
            for (const auto& [v, d] : elements) {
                max_val = max(max_val, d);
            }
            return max_val;
        }
        
        void update_bounds() {
            if (empty()) {
                separator_lower = INF;
                separator_upper = -INF;
            } else {
                separator_lower = min_value();
                separator_upper = max_value();
            }
        }
    };
    
    vector<Block> d0_sequence;
    vector<Block> d1_sequence;
    unordered_map<int, pair<ll, pair<int, int>>> vertex_index;
    int M;
    ll B;

public:
    void initialize(int capacity, ll bound) {
        M = capacity;
        B = bound;
        d0_sequence.clear();
        d1_sequence.clear();
        vertex_index.clear();
        
        d1_sequence.emplace_back(-INF, B);
    }
    
    void insert(int vertex, ll distance) {
        if (distance >= B) return;
        
        if (vertex_index.count(vertex)) {
            auto [old_dist, loc] = vertex_index[vertex];
            if (distance >= old_dist) return;
            remove_vertex(vertex);
        }
        
        int target_block = find_d1_block(distance);
        
        while (target_block >= (int)d1_sequence.size()) {
            d1_sequence.emplace_back(B, B);
        }
        
        d1_sequence[target_block].add(vertex, distance);
        vertex_index[vertex] = {distance, {1, target_block}};
        
        if ((int)d1_sequence[target_block].size() > M) {
            split_d1_block(target_block);
        }
    }
    
    void batch_prepend(const vector<pair<int, ll>>& elements) {
        if (elements.empty()) return;
        
        ll d1_min = B;
        for (const auto& block : d1_sequence) {
            if (!block.empty()) {
                d1_min = min(d1_min, block.min_value());
            }
        }
        
        vector<pair<int, ll>> filtered;
        for (const auto& [vertex, distance] : elements) {
            if (distance < B && distance < d1_min && 
                (!vertex_index.count(vertex) || vertex_index[vertex].first > distance)) {
                if (vertex_index.count(vertex)) {
                    remove_vertex(vertex);
                }
                filtered.push_back({vertex, distance});
            }
        }
        
        if (filtered.empty()) return;
        
        sort(filtered.begin(), filtered.end(),
             [](const pair<int, ll>& a, const pair<int, ll>& b) {
                 if (a.second != b.second) return a.second < b.second;
                 return a.first < b.first;
             });
        
        vector<Block> new_blocks;
        int L = filtered.size();
        int block_size = max(1, M/2);
        
        for (int start = 0; start < L; start += block_size) {
            int end = min(L, start + block_size);
            Block block;
            
            for (int i = start; i < end; i++) {
                block.add(filtered[i].first, filtered[i].second);
            }
            block.update_bounds();
            new_blocks.push_back(block);
        }
        
        int offset = new_blocks.size();
        for (auto& [vertex, data] : vertex_index) {
            if (data.second.first == 0) {
                data.second.second += offset;
            }
        }
        
        for (int i = 0; i < (int)new_blocks.size(); i++) {
            for (const auto& [vertex, distance] : new_blocks[i].elements) {
                vertex_index[vertex] = {distance, {0, i}};
            }
        }
        
        d0_sequence.insert(d0_sequence.begin(), new_blocks.begin(), new_blocks.end());
    }
    
    pair<ll, vector<int>> pull() {
        vector<pair<int, ll>> candidates;
        
        for (const auto& block : d0_sequence) {
            for (const auto& elem : block.elements) {
                candidates.push_back(elem);
                if ((int)candidates.size() >= M) break;
            }
            if ((int)candidates.size() >= M) break;
        }
        
        if ((int)candidates.size() < M) {
            for (const auto& block : d1_sequence) {
                for (const auto& elem : block.elements) {
                    candidates.push_back(elem);
                    if ((int)candidates.size() >= M) break;
                }
                if ((int)candidates.size() >= M) break;
            }
        }
        
        sort(candidates.begin(), candidates.end(),
             [](const pair<int, ll>& a, const pair<int, ll>& b) {
                 if (a.second != b.second) return a.second < b.second;
                 return a.first < b.first;
             });
        
        int extract_count = min(M, (int)candidates.size());
        
        vector<int> result;
        for (int i = 0; i < extract_count; i++) {
            result.push_back(candidates[i].first);
            remove_vertex(candidates[i].first);
        }
        
        ll separator = B;
        if (extract_count < (int)candidates.size()) {
            separator = candidates[extract_count].second;
        } else {
            for (const auto& block : d0_sequence) {
                if (!block.empty()) {
                    separator = min(separator, block.min_value());
                }
            }
            for (const auto& block : d1_sequence) {
                if (!block.empty()) {
                    separator = min(separator, block.min_value());
                }
            }
        }
        
        cleanup_structure();
        return {separator, result};
    }
    
    bool is_empty() const {
        for (const auto& block : d0_sequence) {
            if (!block.empty()) return false;
        }
        for (const auto& block : d1_sequence) {
            if (!block.empty()) return false;
        }
        return true;
    }

private:
    int find_d1_block(ll distance) {
        for (int i = 0; i < (int)d1_sequence.size(); i++) {
            if (distance <= d1_sequence[i].separator_upper || 
                (i == (int)d1_sequence.size() - 1 && distance < B)) {
                return i;
            }
        }
        return max(0, (int)d1_sequence.size() - 1);
    }
    
    void split_d1_block(int block_idx) {
        if (block_idx >= (int)d1_sequence.size()) return;
        
        auto& block = d1_sequence[block_idx];
        if ((int)block.size() <= M) return;
        
        vector<pair<int, ll>> elements = block.elements;
        sort(elements.begin(), elements.end(),
             [](const pair<int, ll>& a, const pair<int, ll>& b) {
                 if (a.second != b.second) return a.second < b.second;
                 return a.first < b.first;
             });
        
        int split_idx = elements.size() / 2;
        ll split_value = elements[split_idx].second;
        
        Block left_block(block.separator_lower, split_value);
        Block right_block(split_value, block.separator_upper);
        
        for (int i = 0; i < split_idx; i++) {
            left_block.add(elements[i].first, elements[i].second);
        }
        for (int i = split_idx; i < (int)elements.size(); i++) {
            right_block.add(elements[i].first, elements[i].second);
        }
        
        left_block.update_bounds();
        right_block.update_bounds();
        
        d1_sequence[block_idx] = left_block;
        d1_sequence.insert(d1_sequence.begin() + block_idx + 1, right_block);
        
        for (auto& [vertex, data] : vertex_index) {
            if (data.second.first == 1 && data.second.second > block_idx) {
                data.second.second++;
            }
        }
        
        if ((int)d1_sequence[block_idx + 1].size() > M) {
            split_d1_block(block_idx + 1);
        }
    }
    
    void remove_vertex(int vertex) {
        if (!vertex_index.count(vertex)) return;
        
        auto [dist, loc] = vertex_index[vertex];
        int seq = loc.first;
        int block_idx = loc.second;
        
        if (seq == 0 && block_idx < (int)d0_sequence.size()) {
            remove_from_block(d0_sequence[block_idx], vertex);
        } else if (seq == 1 && block_idx < (int)d1_sequence.size()) {
            remove_from_block(d1_sequence[block_idx], vertex);
        }
        
        vertex_index.erase(vertex);
    }
    
    void remove_from_block(Block& block, int vertex) {
        block.elements.erase(
            remove_if(block.elements.begin(), block.elements.end(),
                      [vertex](const pair<int, ll>& p) { return p.first == vertex; }),
            block.elements.end());
        block.update_bounds();
    }
    
    void cleanup_structure() {
        d0_sequence.erase(
            remove_if(d0_sequence.begin(), d0_sequence.end(),
                      [](const Block& b) { return b.empty(); }),
            d0_sequence.end());
        
        d1_sequence.erase(
            remove_if(d1_sequence.begin(), d1_sequence.end(),
                      [](const Block& b) { return b.empty(); }),
            d1_sequence.end());
        
        for (auto& block : d0_sequence) {
            if ((int)block.size() > M) {
                split_d0_block(&block - &d0_sequence[0]);
            }
        }
        for (auto& block : d1_sequence) {
            if ((int)block.size() > M) {
                split_d1_block(&block - &d1_sequence[0]);
            }
        }
    }
    
    void split_d0_block(int block_idx) {
        if (block_idx >= (int)d0_sequence.size()) return;
        
        auto& block = d0_sequence[block_idx];
        if ((int)block.size() <= M) return;
        
        vector<pair<int, ll>> elements = block.elements;
        sort(elements.begin(), elements.end(),
             [](const pair<int, ll>& a, const pair<int, ll>& b) {
                 if (a.second != b.second) return a.second < b.second;
                 return a.first < b.first;
             });
        
        int split_idx = elements.size() / 2;
        ll split_value = elements[split_idx].second;
        
        Block left_block(block.separator_lower, split_value);
        Block right_block(split_value, block.separator_upper);
        
        for (int i = 0; i < split_idx; i++) {
            left_block.add(elements[i].first, elements[i].second);
        }
        for (int i = split_idx; i < (int)elements.size(); i++) {
            right_block.add(elements[i].first, elements[i].second);
        }
        
        left_block.update_bounds();
        right_block.update_bounds();
        
        d0_sequence[block_idx] = left_block;
        d0_sequence.insert(d0_sequence.begin() + block_idx + 1, right_block);
        
        for (auto& [vertex, data] : vertex_index) {
            if (data.second.first == 0 && data.second.second > block_idx) {
                data.second.second++;
            }
        }
    }
};

void exact_degree_reduction() {
    vector<vector<Edge>> original_adj = adj;
    int original_n = n;
    
    radj.resize(n + 1);
    for (int u = 1; u <= n; u++) {
        for (const Edge& e : adj[u]) {
            radj[e.to].push_back({u, e.weight});
        }
    }
    
    int vertex_counter = n + 1;
    int max_nodes = original_n + 6 * m;
    
    adj.clear();
    radj.clear();
    adj.resize(max_nodes + 1);
    radj.resize(max_nodes + 1);
    
    for (int v = 1; v <= original_n; v++) {
        const auto& out_edges = original_adj[v];
        const auto& in_edges = radj[v];
        
        if (out_edges.size() + in_edges.size() <= 2) {
            adj[v] = out_edges;
            for (const Edge& e : in_edges) {
                radj[v].push_back(e);
            }
            continue;
        }
        
        vector<int> cycle_nodes;
        int total_degree = out_edges.size() + in_edges.size();
        
        for (int i = 0; i < total_degree; i++) {
            cycle_nodes.push_back(vertex_counter++);
        }
        
        for (int i = 0; i < total_degree; i++) {
            int next = (i + 1) % total_degree;
            adj[cycle_nodes[i]].push_back({cycle_nodes[next], 0});
            radj[cycle_nodes[next]].push_back({cycle_nodes[i], 0});
        }
        
        adj[v].push_back({cycle_nodes[0], 0});
        radj[cycle_nodes[0]].push_back({v, 0});
        
        int idx = 1;
        for (const Edge& e : out_edges) {
            adj[cycle_nodes[idx]].push_back({e.to, e.weight});
            idx++;
        }
        
        for (const Edge& e : in_edges) {
            radj[cycle_nodes[idx]].push_back({e.to, e.weight});
            idx++;
        }
    }
    
    n = vertex_counter - 1;
    adj.resize(n + 1);
    radj.resize(n + 1);
}

bool relax_paper_exact(int u, int v, ll weight) {
    ll candidate = db[u] + weight;
    if (candidate < db[v]) {
        db[v] = candidate;
        pred[v] = u;
        return true;
    } else if (candidate == db[v] && pred[v] != -1 && tie_breaker(u, pred[v])) {
        pred[v] = u;
        return false;
    }
    return false;
}

pair<vector<int>, vector<int>> find_pivots_paper_spec(ll B, const vector<int>& S) {
    vector<int> W = S;
    unordered_set<int> W_set(S.begin(), S.end());
    vector<int> frontier = S;
    
    for (int step = 1; step <= k && (int)W.size() <= k * (int)S.size(); step++) {
        vector<int> next_frontier;
        
        for (int u : frontier) {
            for (const Edge& e : adj[u]) {
                int v = e.to;
                ll new_dist = db[u] + e.weight;
                
                if (new_dist < B && new_dist <= db[v]) {
                    bool relaxed = relax_paper_exact(u, v, e.weight);
                    
                    if ((relaxed || new_dist == db[v]) && !W_set.count(v)) {
                        W.push_back(v);
                        next_frontier.push_back(v);
                        W_set.insert(v);
                    }
                }
            }
        }
        
        frontier = next_frontier;
        if (frontier.empty()) break;
    }
    
    unordered_map<int, vector<int>> children;
    unordered_set<int> has_parent;
    
    for (int u : W) {
        for (const Edge& e : adj[u]) {
            int v = e.to;
            if (W_set.count(v) && db[v] == db[u] + e.weight && pred[v] == u) {
                children[u].push_back(v);
                has_parent.insert(v);
            }
        }
    }
    
    function<int(int)> compute_size = [&](int root) -> int {
        int size = 1;
        for (int child : children[root]) {
            size += compute_size(child);
        }
        return size;
    };
    
    vector<int> P;
    for (int s : S) {
        if (W_set.count(s) && !has_parent.count(s) && compute_size(s) >= k) {
            P.push_back(s);
        }
    }
    
    return {P, W};
}

pair<ll, vector<int>> base_case_paper_spec(ll B, const vector<int>& S) {
    if (S.empty()) {
        return {B, {}};
    }
    
    vector<int> result;
    priority_queue<pair<ll, int>, vector<pair<ll, int>>, greater<pair<ll, int>>> pq;
    vector<bool> processed(n + 1, false);
    
    for (int source : S) {
        if (db[source] < B) {
            pq.push({db[source], source});
        }
    }
    
    while (!pq.empty() && (int)result.size() < k) {
        auto [d, u] = pq.top();
        pq.pop();
        
        if (processed[u] || d > db[u] || d >= B) continue;
        processed[u] = true;
        
        result.push_back(u);
        
        for (const Edge& e : adj[u]) {
            int v = e.to;
            ll new_dist = d + e.weight;
            
            if (new_dist < B && !processed[v]) {
                if (relax_paper_exact(u, v, e.weight)) {
                    pq.push({db[v], v});
                }
            }
        }
    }
    
    ll boundary = ((int)result.size() == k && !result.empty()) ? db[result[k-1]] : B;
    return {boundary, result};
}

pair<ll, vector<int>> bounded_mssp_paper_spec(int level, ll B, const vector<int>& S) {
    if (S.empty()) return {B, {}};
    if (level == 0) return base_case_paper_spec(B, S);
    
    auto [P, W] = find_pivots_paper_spec(B, S);
    
    DataStructureD D;
    
    long long target_ll = k;
    int safe_exp = min(level * t, min(60, (int)log2(n) + 5));
    for (int i = 0; i < safe_exp && target_ll <= (1LL << 30); i++) {
        target_ll *= 2;
    }
    int target_size = (int)min(target_ll, (long long)n);
    
    D.initialize(max(1, target_size), B);
    
    for (int p : P) {
        D.insert(p, db[p]);
    }
    
    ll B_prime = B;
    if (!P.empty()) {
        B_prime = INF;
        for (int p : P) {
            B_prime = min(B_prime, db[p]);
        }
    }
    
    vector<int> U;
    ll current_bound = B_prime;
    
    while ((int)U.size() < target_size && !D.is_empty()) {
        auto [Bi, Si] = D.pull();
        if (Si.empty()) break;
        
        auto [Bi_prime, Ui] = bounded_mssp_paper_spec(level - 1, Bi, Si);
        
        for (int u : Ui) {
            U.push_back(u);
        }
        current_bound = Bi_prime;
        
        vector<pair<int, ll>> insert_list, batch_list;
        unordered_set<int> processed_vertices;
        
        for (int u : Ui) {
            for (const Edge& e : adj[u]) {
                int v = e.to;
                ll new_dist = db[u] + e.weight;
                
                if (new_dist < B && new_dist <= db[v] && !processed_vertices.count(v)) {
                    ll old_dist = db[v];
                    bool relaxed = relax_paper_exact(u, v, e.weight);
                    
                    if (relaxed || new_dist == old_dist) {
                        processed_vertices.insert(v);
                        ll curr_dist = db[v];
                        
                        if (curr_dist >= Bi && curr_dist < B) {
                            insert_list.push_back({v, curr_dist});
                        } else if (curr_dist >= Bi_prime && curr_dist < Bi) {
                            batch_list.push_back({v, curr_dist});
                        }
                    }
                }
            }
        }
        
        for (const auto& [v, d] : insert_list) {
            D.insert(v, d);
        }
        if (!batch_list.empty()) {
            D.batch_prepend(batch_list);
        }
        
        if ((int)U.size() >= target_size) break;
    }
    
    for (int w : W) {
        if (db[w] < current_bound && find(U.begin(), U.end(), w) == U.end()) {
            U.push_back(w);
        }
    }
    
    return {min(current_bound, B), U};
}

vector<ll> duan_sssp_exact_implementation(int source) {
    int original_n = n;
    
    exact_degree_reduction();
    
    double log_n = log2(max(2.0, (double)n));
    k = max(1, (int)floor(pow(log_n, 1.0/3.0)));
    t = max(1, (int)floor(pow(log_n, 2.0/3.0)));
    max_level = max(1, (int)ceil(log_n / max(1.0, (double)t)));
    
    db.assign(n + 1, INF);
    pred.assign(n + 1, -1);
    db[source] = 0;
    
    vector<int> initial = {source};
    auto [bound, vertices] = bounded_mssp_paper_spec(max_level, INF, initial);
    
    vector<ll> result;
    for (int i = 1; i <= original_n; i++) {
        result.push_back(db[i]);
    }
    return result;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);
    
    cin >> n >> m;
    adj.resize(n + 1);
    
    for (int i = 0; i < m; i++) {
        int u, v;
        ll w;
        cin >> u >> v >> w;
        adj[u].push_back({v, w});
    }
    
    vector<ll> distances = duan_sssp_exact_implementation(1);
    
    for (size_t i = 0; i < distances.size(); i++) {
        cout << distances[i];
        if (i < distances.size() - 1) cout << " ";
    }
    cout << "\n";
    
    return 0;
}
