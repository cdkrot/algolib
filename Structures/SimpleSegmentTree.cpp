struct MinF {
    template<typename T>
    T operator()(T a, T b) const {
        return (a < b ? a : b);
    }

    T neutral() const {
        return std::numeric_limits<T>::min();
    }
};

struct MaxF {
    template<typename T>
    T operator()(T a, T b) const {
        return (a > b ? a : b);
    }

    T neutral() const {
        return std::numeric_limits<T>::max();
    }
};

template <typename T, typename FunctionT>
struct SegmentTree {
    static inline FunctionT function{};
    
    segment_tree(const vector<T>& in): n (in.size()), tree(2 * in.size()) {
        for (int i = n - 1; i >= 1; --i) {
            tree[i] = function(tree[2 * i], tree[2 * i + 1]);
        }
    };

    T get(int l, int r) {
        l += n;
        r += n + 1;
        T res = function.neutral();
        while (l < r) {
            if (l % 2 == 1) {
                res = function(tree[l++], res);
            }
            if (r % 2 == 1) {
                res = function(tree[--r], res);
            }
            l /= 2;
            r /= 2;
        }
    }
        
    int n;
    vector<T> tree;
};
