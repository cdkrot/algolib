namespace fft {
// partially adapted from cryptozoology notebook
const int mod = 998244353;
const int L = 22;
const int N = 1 << L;

class Poly: public std::vector<int> {
public:
    using std::vector<int>::vector;

    void trim() {
        while (size() > 0 && back() == 0) {
            pop_back();
        }
    }

    void add_monom(int i, int val) {
        if (size_t(i) >= size()) {
            resize(i + 1, 0);
        }
        int& x = (*this)[i];
        if ((x += val) >= mod) {
            x -= mod;
        }
    }

    void add(const std::vector<int>& other) {
        resize(max(size(), other.size()), 0);
        
        for (size_t i = 0; i < other.size(); ++i) {
            int& x = (*this)[i];
            if ((x += other[i]) >= mod) {
                x -= mod;
            }
        }
        trim();
    }

    void sub(const std::vector<int>& other) {
        resize(max(size(), other.size()), 0);
        
        for (size_t i = 0; i < other.size(); ++i) {
            int& x = (*this)[i];
            if ((x -= other[i]) < 0) {
                x += mod;
            }
        }
        trim();
    }
};


long long pw(long long a, long long b) {
	long long res = 1;
	while (b) {
		if (b & 1ll) {
			res = res * a % mod;
		}
		b >>= 1;
		a = a * a % mod;
	}
	return res;
}

int getRoot() {
	int root = 1;
	while (pw(root, 1 << L) != 1 || pw(root, 1 << (L - 1)) == 1) {
		++root;
	}
	return root;
}

const int root = getRoot();

long long angles[N + 1];
int bitrev[N];

inline int revBit(int x, int len) {
	return bitrev[x] >> (L - len);
}

void fft(vector<int>& a, bool inverse = false) {
	int n = a.size();
	assert(!(n & (n - 1)));	// work only with powers of two
	int l = __builtin_ctz(n);

	for (int i = 0; i < n; ++i) {
		int j = revBit(i, l);
		if (i < j) {
			swap(a[i], a[j]);
		}
	}

	for (int len = 1; len < n; len *= 2) {
		for (int start = 0; start < n; start += 2 * len) {
			for (int i = 0; i < len; ++i) {
				long long x = a[start + i], y = a[start + len + i];
				int idx = N / 2 / len * i;
				long long w = angles[inverse ? N - idx : idx];
				w = w * y % mod;
				a[start + i] = x + w;
				if (a[start + i] >= mod) {
					a[start + i] -= mod;
				}
				a[start + len + i] = x - w;
				if (a[start + len + i] < 0) {
					a[start + len + i] += mod;
				}
			}
		}
	}

	if (inverse) {
		int rev_deg = 1;
		for (int i = 0; i < l; ++i) {
			rev_deg = (rev_deg % 2) ? ((rev_deg + mod) / 2) : (rev_deg / 2);
		}
		for (auto& x : a) {
			x = x * int64_t(rev_deg) % mod;
		}
	}
}

Poly multiply(Poly a, Poly b) {
	int n = 1;
	while (n < (int)a.size() || n < (int)b.size()) {
		n *= 2;
	}
	a.resize(n + n);
	b.resize(n + n);
	fft(a);
	fft(b);
	for (int i = 0; i < n + n; ++i) {
		a[i] = a[i] * int64_t(b[i]) % mod;
	}
	fft(a, true);
    a.trim();
	return a;
}

Poly square(Poly a) {
    int n = 1;
	while (n < (int)a.size()) {
		n *= 2;
	}
	a.resize(n + n);
	fft(a);
	for (int i = 0; i < n + n; ++i) {
		a[i] = a[i] * int64_t(a[i]) % mod;
	}
	fft(a, true);
    a.trim();
	return a;
}

auto fft_init = []() {
	angles[0] = 1;
	for (int i = 1; i <= N; ++i) {
		angles[i] = angles[i - 1] * root % mod;
	}

	for (int i = 0; i < N; ++i) {
		int x = i;
		for (int j = 0; j < L; ++j) {
			bitrev[i] = (bitrev[i] << 1) | (x & 1);
			x >>= 1;
		}
	}
    return std::monostate{};
}();
}
