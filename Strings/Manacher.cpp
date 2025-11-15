// mana[i] == length of max even palindrome centered at i-0.5.
// containerT can be vector or string.
template <typename ContainerT> vector<int> get_mana(const ContainerT& s) {
    vector<int> mana(SZ(s), -1);

    int L = -1, R = -1;
    for (int i = 1; i < SZ(s); ++i) {
        if (i <= R)
            mana[i] = min(R - i + 1, mana[L + R - i + 1]);
            
        if (s[i - 1] == s[i]) {
            while (i + mana[i] < SZ(s) and i >= 1 + mana[i] and s[i + mana[i]] == s[i - 1 - mana[i]])
                ++mana[i];

            if (i + mana[i] - 1 > R)
                L = i - mana[i], R = i + mana[i] - 1;
        }
    }
    return mana;
}
