#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
const ll INF = 1e18;

struct Edge {
    int to;
    ll weight;
};

int n, m;
vector<vector<Edge>> adj;
vector<ll> dist;
vector<int> pred;
vector<bool> complete;
int k, t, lmax;

mt19937 rng(123456); 

struct DataStructureD {
    struct Block {
        ll min_dist, max_dist;
        vector<int> vertices;
        ll range_bound; 
        
        Block(ll B, int k_param) {
            min_dist = INF;
            max_dist = -INF;
            range_bound = max(1LL, B / max(1LL, (1LL * k_param * k_param)));
        }
        
        void add_vertex(int v, ll d) {
            if (vertices.empty()) {
                min_dist = max_dist = d;
            } else {
                min_dist = min(min_dist, d);
                max_dist = max(max_dist, d);
            }
            vertices.push_back(v);
        }
        
        bool can_fit(ll new_dist) const {
            if (vertices.empty()) return true;
            ll new_min = min(min_dist, new_dist);
            ll new_max = max(max_dist, new_dist);
            return (new_max - new_min) <= range_bound;
        }
        
        bool is_full(int capacity) const {
            return (int)vertices.size() >= capacity;
        }
    };
    
    vector<Block> small_blocks; 
    vector<pair<ll, Block>> large_blocks; 
    
    unordered_map<int, pair<ll, int>> vertex_info; 
    unordered_set<int> all_vertices; 
    
    uint64_t capacity;
    ll bound;
    int k_param;
    int block_size;
    
    DataStructureD(ll B, uint64_t M, int k) : bound(B), capacity(M), k_param(k) {
        block_size = max(1, (int)sqrt(M));
        small_blocks.reserve(block_size * 2);
        large_blocks.reserve(block_size);
    }
    
    bool empty() const {
        return small_blocks.empty() && large_blocks.empty();
    }
    
    void insert(int vertex, ll distance) {
        if (distance >= bound) return;
        
        if (all_vertices.count(vertex)) {
            auto existing = vertex_info[vertex];
            if (distance >= existing.first) return;
            remove_vertex(vertex);
        }
        
        bool inserted = false;
        for (int i = 0; i < (int)small_blocks.size(); i++) {
            auto& block = small_blocks[i];
            if (!block.is_full(block_size) && block.can_fit(distance)) {
                block.add_vertex(vertex, distance);
                vertex_info[vertex] = {distance, i};
                all_vertices.insert(vertex);
                inserted = true;
                break;
            }
        }
        
        if (!inserted) {
            Block new_block(bound, k_param);
            new_block.add_vertex(vertex, distance);
            small_blocks.push_back(new_block);
            vertex_info[vertex] = {distance, (int)small_blocks.size() - 1};
            all_vertices.insert(vertex);
        }
        
        promote_full_blocks();
    }
    
    void batch_prepend(vector<pair<int, ll>>& pairs) {
        if (pairs.empty()) return;

        sort(pairs.begin(), pairs.end(), 
             [](const pair<int,ll>& a, const pair<int,ll>& b) { 
                 return a.second < b.second; 
             });
        
        vector<pair<int, ll>> clean_pairs;
        for (const auto& [vertex, distance] : pairs) {
            if (distance < bound && !all_vertices.count(vertex)) {
                clean_pairs.push_back({vertex, distance});
                all_vertices.insert(vertex);
            }
        }
        
        vector<Block> new_blocks;
        for (int start = clean_pairs.size() - 1; start >= 0; ) {
            Block block(bound, k_param);
            int end = max(0, start - block_size + 1);
            
            for (int i = start; i >= end; i--) {
                block.add_vertex(clean_pairs[i].first, clean_pairs[i].second);
            }
            
            new_blocks.push_back(block);
            start = end - 1;
        }
        
        small_blocks.insert(small_blocks.begin(), new_blocks.begin(), new_blocks.end());
        reindex_blocks();
        promote_full_blocks();
    }
    
    pair<vector<int>, ll> pull() {
        vector<int> result;
        uint64_t count = 0;
        
        while (!small_blocks.empty() && count < capacity) {
            auto& block = small_blocks.front();
            
            while (!block.vertices.empty() && count < capacity) {
                int vertex = block.vertices.back();
                block.vertices.pop_back();
                vertex_info.erase(vertex);
                all_vertices.erase(vertex);
                result.push_back(vertex);
                count++;
            }
            
            if (block.vertices.empty()) {
                small_blocks.erase(small_blocks.begin());
                reindex_blocks();
            }
        }
        
        while (!large_blocks.empty() && count < capacity) {
            auto& [min_dist, block] = large_blocks.front();
            
            while (!block.vertices.empty() && count < capacity) {
                int vertex = block.vertices.back();
                block.vertices.pop_back();
                vertex_info.erase(vertex);
                all_vertices.erase(vertex);
                result.push_back(vertex);
                count++;
            }
            
            if (block.vertices.empty()) {
                large_blocks.erase(large_blocks.begin());
            }
        }
        
        ll separator = bound;
        if (!small_blocks.empty()) {
            separator = small_blocks.front().min_dist;
        } else if (!large_blocks.empty()) {
            separator = large_blocks.front().second.min_dist;
        }
        
        return {result, separator};
    }
    
private:
    void promote_full_blocks() {
        auto it = small_blocks.begin();
        while (it != small_blocks.end()) {
            if (it->is_full(block_size)) {
                large_blocks.push_back({it->min_dist, move(*it)});
                it = small_blocks.erase(it);
            } else {
                ++it;
            }
        }
        
        sort(large_blocks.begin(), large_blocks.end(),
             [](const pair<ll, Block>& a, const pair<ll, Block>& b) {
                 return a.first < b.first;
             });
        reindex_blocks();
    }
    
    void reindex_blocks() {
        for (int i = 0; i < (int)small_blocks.size(); i++) {
            for (int vertex : small_blocks[i].vertices) {
                if (vertex_info.count(vertex)) {
                    vertex_info[vertex].second = i;
                }
            }
        }
    }
    
    void remove_vertex(int vertex) {
        if (!all_vertices.count(vertex)) return;
        
        auto [distance, block_idx] = vertex_info[vertex];
        
        if (block_idx >= 0 && block_idx < (int)small_blocks.size()) {
            auto& vertices = small_blocks[block_idx].vertices;
            auto pos = find(vertices.begin(), vertices.end(), vertex);
            if (pos != vertices.end()) {
                swap(*pos, vertices.back());
                vertices.pop_back();
            }
        } else {
            for (auto& [min_dist, block] : large_blocks) {
                auto& vertices = block.vertices;
                auto pos = find(vertices.begin(), vertices.end(), vertex);
                if (pos != vertices.end()) {
                    swap(*pos, vertices.back());
                    vertices.pop_back();
                    break;
                }
            }
        }
        
        vertex_info.erase(vertex);
        all_vertices.erase(vertex);
    }
};

struct DegreeTransform {
    int original_n;
    vector<vector<Edge>> new_adj;
    int new_n;
    
    DegreeTransform(int n, const vector<vector<Edge>>& adj) : original_n(n) {
        transform(adj);
    }
    
private:
    void transform(const vector<vector<Edge>>& adj) {
        new_n = original_n;
        
        int total_edges = 0;
        for (int u = 1; u <= original_n; u++) {
            total_edges += adj[u].size();
        }
        
        new_adj.resize(original_n + 2 * total_edges + 5);
        
        for (int u = 1; u <= original_n; u++) {
            const auto& edges = adj[u];
            
            if (edges.size() <= 2) {
                new_adj[u] = edges;
            } else {
                vector<int> cycle;
                for (size_t i = 0; i < edges.size(); i++) {
                    cycle.push_back(++new_n);
                }
                
                for (size_t i = 0; i < cycle.size(); i++) {
                    int curr = cycle[i];
                    int next = cycle[(i + 1) % cycle.size()];
                    new_adj[curr].push_back({next, 0});
                }
                
                for (size_t i = 0; i < edges.size(); i++) {
                    new_adj[cycle[i]].push_back(edges[i]);
                }
                
                new_adj[u].clear();
                new_adj[u].push_back({cycle[0], 0});
            }
        }
        
        new_adj.resize(new_n + 1);
    }
};

void init_params(int graph_size) {
    double log_n = max(1.0, log2(graph_size));
    k = max(1, (int)floor(pow(log_n, 1.0/3.0)));
    t = max(1, (int)floor(pow(log_n, 2.0/3.0)));
    lmax = max(1, min(15, (int)ceil(log_n / t))); 
}

pair<vector<int>, vector<int>> find_separators(ll B, const vector<int>& sources) {
    vector<ll> temp_dist = dist;
    unordered_set<int> reachable(sources.begin(), sources.end());
    queue<int> q;
    
    for (int s : sources) {
        q.push(s);
    }
    
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        
        vector<int> edge_order(adj[u].size());
        iota(edge_order.begin(), edge_order.end(), 0);
        shuffle(edge_order.begin(), edge_order.end(), rng);
        
        for (int idx : edge_order) {
            const Edge& e = adj[u][idx];
            int v = e.to;
            ll new_dist = temp_dist[u] + e.weight;
            
            if (new_dist < B && new_dist < temp_dist[v]) {
                temp_dist[v] = new_dist;
                
                if (reachable.find(v) == reachable.end()) {
                    reachable.insert(v);
                    q.push(v);
                    
                    if (reachable.size() > 8ULL * k * sources.size()) {
                        goto done;
                    }
                }
            }
        }
    }
    
    done:
    if (reachable.size() > 8ULL * k * sources.size()) {
        return {sources, vector<int>(reachable.begin(), reachable.end())};
    }
    
    vector<int> reachable_vec(reachable.begin(), reachable.end());
    unordered_map<int, int> parent;
    unordered_map<int, vector<int>> children;
    
    for (int u : reachable_vec) {
        for (const Edge& e : adj[u]) {
            int v = e.to;
            if (reachable.count(v) && temp_dist[u] + e.weight == temp_dist[v]) {
                if (parent.find(v) == parent.end() || 
                    temp_dist[u] < temp_dist[parent[v]]) {
                    parent[v] = u;
                }
            }
        }
    }
    
    for (const auto& [child, par] : parent) {
        children[par].push_back(child);
    }
    
    function<int(int)> get_subtree_size = [&](int u) -> int {
        int size = 1;
        for (int child : children[u]) {
            size += get_subtree_size(child);
        }
        return size;
    };
    
    vector<int> separators;
    double c = 2.0;
    double log_k = max(1.0, log2(k));
    
    for (int attempt = 0; attempt < 3 && separators.empty(); attempt++) {
        for (int s : sources) {
            if (reachable.count(s) && parent.find(s) == parent.end()) {
                int size = get_subtree_size(s);
                
                if (size >= k) {
                    double log_size = max(1.0, log2(size));
                    double prob = min(1.0, c * log_size / log_k * (1.0 + 0.1 * attempt));
                    
                    uniform_real_distribution<double> coin(0.0, 1.0);
                    if (coin(rng) < prob) {
                        separators.push_back(s);
                    }
                }
            }
        }
    }
    
    if (separators.empty()) {
        for (int s : sources) {
            if (reachable.count(s) && parent.find(s) == parent.end()) {
                if (get_subtree_size(s) >= k) {
                    separators.push_back(s);
                    break;
                }
            }
        }
    }
    
    return {separators, reachable_vec};
}

pair<ll, vector<int>> base_case(ll B, const vector<int>& sources) {
    if (sources.size() != 1) {
        throw runtime_error("Base case needs exactly one source");
    }
    
    int source = sources[0];
    vector<int> result;
    
    priority_queue<pair<ll, int>, vector<pair<ll, int>>, greater<pair<ll, int>>> pq;
    vector<bool> processed(n + 1, false);
    
    pq.push({dist[source], source});
    
    while (!pq.empty() && (int)result.size() < k) {
        auto [d, u] = pq.top();
        pq.pop();
        
        if (processed[u] || d > dist[u] || d >= B) continue;
        
        processed[u] = true;
        result.push_back(u);
        complete[u] = true;
        
        for (const Edge& e : adj[u]) {
            int v = e.to;
            ll new_dist = d + e.weight;
            
            if (new_dist < dist[v] && new_dist < B) {
                dist[v] = new_dist;
                pred[v] = u;
                
                if (!processed[v]) {
                    pq.push({new_dist, v});
                }
            }
        }
    }
    
    if ((int)result.size() >= k) {
        ll kth_dist = dist[result[k-1]];
        while ((int)result.size() > k && dist[result.back()] > kth_dist) {
            result.pop_back();
        }
        return {kth_dist, result};
    }
    
    return {B, result};
}

pair<ll, vector<int>> solve_recursive(int level, ll B, const vector<int>& sources) {
    if (sources.empty()) return {B, {}};
    
    if (level == 0) {
        return base_case(B, sources);
    }
    
    auto [separators, reachable] = find_separators(B, sources);
    
    uint64_t capacity = 1;
    int exp = (level - 1) * t;
    if (exp <= 60) {
        capacity = 1ULL << exp;
    } else {
        capacity = UINT64_MAX / 64;
    }
    
    DataStructureD data_structure(B, capacity, k);
    
    for (int sep : separators) {
        data_structure.insert(sep, dist[sep]);
    }
    
    vector<int> all_results;
    unordered_set<int> result_set;
    ll final_bound = B;
    
    while (!data_structure.empty()) {
        auto [current_sources, bound_i] = data_structure.pull();
        
        if (current_sources.empty()) continue;
        
        auto [bound_i_prime, partial_results] = solve_recursive(level - 1, bound_i, current_sources);
        
        for (int u : partial_results) {
            if (result_set.find(u) == result_set.end()) {
                all_results.push_back(u);
                result_set.insert(u);
            }
        }
        
        final_bound = min(final_bound, bound_i_prime);
        
        vector<pair<int, ll>> batch_items;
        unordered_set<int> batch_set;
        
        for (int u : partial_results) {
            for (const Edge& e : adj[u]) {
                int v = e.to;
                ll new_dist = dist[u] + e.weight;
                
                if (new_dist < dist[v] || (new_dist == dist[v] && u < pred[v])) {
                    dist[v] = new_dist;
                    pred[v] = u;
                }
                
                if (new_dist < B && batch_set.find(v) == batch_set.end()) {
                    if (new_dist > bound_i_prime) {
                        data_structure.insert(v, new_dist);
                    } else if (new_dist <= bound_i) {
                        batch_items.push_back({v, new_dist});
                        batch_set.insert(v);
                    }
                }
            }
        }
        
        if (!batch_items.empty()) {
            data_structure.batch_prepend(batch_items);
        }
    }
    
    for (int w : reachable) {
        if (dist[w] < final_bound && result_set.find(w) == result_set.end()) {
            all_results.push_back(w);
            result_set.insert(w);
            complete[w] = true;
        }
    }
    
    return {final_bound, all_results};
}

struct RadixHeap {
    static const int BUCKET_COUNT = 64; 
    vector<vector<pair<ll, int>>> buckets;
    vector<ll> bucket_mins;
    ll global_min;
    
    RadixHeap() : buckets(BUCKET_COUNT), bucket_mins(BUCKET_COUNT, INF), global_min(0) {}
    
    bool empty() const {
        for (int i = 0; i < BUCKET_COUNT; i++) {
            if (!buckets[i].empty()) return false;
        }
        return true;
    }
    
    void push(ll distance, int vertex) {
        if (distance < global_min) return; 
        
        int bucket = get_bucket(distance);
        buckets[bucket].push_back({distance, vertex});
        bucket_mins[bucket] = min(bucket_mins[bucket], distance);
    }
    
    pair<ll, int> pop() {
        int bucket = 0;
        while (bucket < BUCKET_COUNT && buckets[bucket].empty()) {
            bucket++;
        }
        
        if (bucket == BUCKET_COUNT) {
            return {INF, -1};
        }
        
        if (bucket == 0 || buckets[bucket].size() == 1) {
            auto result = buckets[bucket].back();
            buckets[bucket].pop_back();
            
            if (buckets[bucket].empty()) {
                bucket_mins[bucket] = INF;
            }
            
            global_min = max(global_min, result.first);
            return result;
        }
        
        redistribute_bucket(bucket);
        return pop(); 
    }
    
private:
    int get_bucket(ll distance) const {
        if (distance == global_min) return 0;
        ll diff = distance - global_min;
        return min(BUCKET_COUNT - 1, 63 - __builtin_clzll(diff));
    }
    
    void redistribute_bucket(int bucket) {
        vector<pair<ll, int>> items = move(buckets[bucket]);
        buckets[bucket].clear();
        bucket_mins[bucket] = INF;
        
        ll actual_min = INF;
        for (const auto& [dist, vertex] : items) {
            actual_min = min(actual_min, dist);
        }
        
        global_min = actual_min;
        
        for (const auto& [dist, vertex] : items) {
            int new_bucket = get_bucket(dist);
            buckets[new_bucket].push_back({dist, vertex});
            bucket_mins[new_bucket] = min(bucket_mins[new_bucket], dist);
        }
    }
};

void final_cleanup() {
    RadixHeap heap;
    vector<bool> in_heap(n + 1, false);
    
    for (int i = 1; i <= n; i++) {
        if (dist[i] < INF) {
            heap.push(dist[i], i);
            in_heap[i] = true;
        }
    }
    
    while (!heap.empty()) {
        auto [d, u] = heap.pop();
        in_heap[u] = false;
        
        if (d > dist[u]) continue;
        
        for (const Edge& e : adj[u]) {
            int v = e.to;
            ll new_dist = d + e.weight;
            
            if (new_dist < dist[v]) {
                dist[v] = new_dist;
                pred[v] = u;
                
                if (!in_heap[v]) {
                    heap.push(new_dist, v);
                    in_heap[v] = true;
                }
            }
        }
    }
}

vector<ll> fast_sssp(int source) {
    DegreeTransform transform(n, adj);
    n = transform.new_n;
    adj = transform.new_adj;
    
    init_params(n);
    
    dist.assign(n + 1, INF);
    pred.assign(n + 1, -1);
    complete.assign(n + 1, false);
    
    dist[source] = 0;
    complete[source] = true;
    
    vector<int> initial = {source};
    solve_recursive(lmax, INF, initial);
    
    final_cleanup();
    
    vector<ll> result;
    for (int i = 1; i <= transform.original_n; i++) {
        result.push_back(dist[i]);
    }
    
    return result;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);
    
    cin >> n >> m;
    adj.resize(n + 1);
    
    for (int i = 0; i < m; i++) {
        int a, b;
        ll c;
        cin >> a >> b >> c;
        adj[a].push_back({b, c});
    }
    
    vector<ll> distances = fast_sssp(1);
    
    for (size_t i = 0; i < distances.size(); i++) {
        cout << distances[i];
        if (i < distances.size() - 1) cout << " ";
    }
    cout << "\n";
    
    return 0;
}