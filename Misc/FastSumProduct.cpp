constexpr const uint64_t mod28 = 8 * uint64_t(mod) * mod;

int fastmult(const int*__restrict a, const int* __restrict b, int m) {
    uint64_t ans = 0;
    int forb = 0;
    for (forb = 0; forb + 8 <= m;) {
        for (int t = 0; t < 8; ++t, ++forb) {
            ans += a[max_n - 1 - m + forb] * uint64_t(1) * b[forb + 1];
        }
        if (ans >= mod28) {
            ans -= mod28;
        }
    }

    for (; forb <= m; ++forb) {
        ans += a[max_n - 1 - m + forb] * uint64_t(1) * b[forb + 1];
    }
    ans %= mod;

    return ans;
}
