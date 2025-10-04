template <typename T>
struct vec {
    T x{};
    T y{};

    bool operator==(vec<T> other) const {
        return x == other.x && y == other.y;
    }

    bool operator!=(vec<T> other) const {
        return !(*this == other);
    }

    vec<T> operator-(vec<T> other) const {
        return vec<T> {x - other.x, y - other.y};
    }

    template <typename Stream>
    friend Stream& operator>>(Stream& s, vec<T>& self) {
        return s >> self.x >> self.y;
    }
};

int64_t dot(vec<int> a, vec<int> b) {
    return a.x * int64_t(b.x) + a.y * int64_t(b.y);
}

int64_t cross(vec<int> a, vec<int> b) {
    return a.x * int64_t(b.y) - a.y * int64_t(b.x);
}
