struct Centroid {
    Centroid(const vector<vector<int>>& tree_): n(tree_.size()), tree(tree_), mark(n, false), sizes(n) {
        build(0);
    }

    fft::Poly get_d() {
        return result;
    }

private:
    void build(int v) {
        int K = compute_sizes(v);

        while (true) {
            int next_v = v;
            for (int u: tree[v]) {
                if (u == v || mark[u] || sizes[u] > sizes[v]) {
                    continue;
                }
                if (2 * sizes[u] >= K) {
                    next_v = u;
                    break;
                }
            }
            if (next_v == v) {
                break;
            }
            v = next_v;
        }

        process(v);
        mark[v] = true;
        for (int u: tree[v]) {
            if (!mark[u]) {
                build(u);
            }
        }
    }

    void process(int v) {
        fft::Poly total = {1};
        fft::Poly result_this;

        for (int u : tree[v]) {
            if (mark[u]) {
                continue;
            }

            fft::Poly local{};

            std::function<void(int, int, int)> go (
                [&local, &total, this, &go](int vert, int par, int dist) {
                    local.add_monom(dist, 1);
                    total.add_monom(dist, 1);

                    for (int neigh: tree[vert]) {
                        if (mark[neigh] || neigh == par) {
                            continue;
                        }

                        go(neigh, vert, dist + 1);   
                    }
                }
            );

            go(u, v, 1);
            result_this.sub(fft::square(local));
        }
        result_this.add(fft::square(total));
        result_this.add_monom(0, 1);
        for (int& x: result_this) {
            LASSERT(x % 2 == 0);
            x /= 2;
        }
        result.add(result_this);
    }

    int compute_sizes(int v, int p = -1) {
        sizes[v] = 1;
        for (int u: tree[v]) {
            if (mark[u] || u == p) {
                continue;
            }
            sizes[v] += compute_sizes(u, v);
        }
        return sizes[v];
    }

    int n;
    const vector<vector<int>>& tree;
    vector<int> mark;
    vector<int> sizes;

    fft::Poly result;
};
