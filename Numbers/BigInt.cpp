struct BigInt {
    BigInt(): BigInt(0) {}
    BigInt(int val) {
        assert(val >= 0);
        while (val != 0) {
            digits.push_back(val % 10);
            val /= 10;
        }
    }

    BigInt(const std::string& s): digits(s.rbegin(), s.rend()) {
        for (char& ch: digits) {
            ch -= '0';
        }
        trim();
    }

    std::strong_ordering operator<=>(const BigInt& other) const {
        if (num_digits() < other.num_digits()) {
            return std::strong_ordering::less;
        } else if (num_digits() > other.num_digits()) {
            return std::strong_ordering::greater;
        }

        for (int i = 0; i < num_digits(); ++i) {
            if (nth_highest_digit(i) < other.nth_highest_digit(i)) {
                return std::strong_ordering::less;
            } else if (nth_highest_digit(i) < other.nth_highest_digit(i)) {
                return std::strong_ordering::greater;
            }
        }
        return std::strong_ordering::equal;
    }

    int operator%(int mod) {
        int64_t result = 0;
        for (auto it = digits.rbegin(); it != digits.rend(); ++it) {
            result = (10 * result + *it) % mod; 
        }
        return int(result);
    }
    BigInt& operator-=(const BigInt& other) {
        assert(num_digits() >= other.num_digits());
        for (int i = 0; i < other.num_digits(); ++i) {
            digits[i] -= other.digits[i];
            if (digits[i] < 0) {
                digits[i] += 10;
                assert(i + 1 < num_digits());
                digits[i + 1] -= 1;
            }
        }
        return *this;
    }
    BigInt operator-(const BigInt& other) const {
        BigInt res = *this;
        res -= other;
        return res;
    }

    int nth_highest_digit(int i) const {
        return nth_lowest_digit(num_digits() - 1 - i);
    }
    int nth_lowest_digit(int i) const {
        return digits[i];
    }
    int num_digits() const {
        return digits.size(); 
    }

private:
    void trim() {
        while (digits.size() > 0 && digits.back() == 0) {
            digits.pop_back();
        }
    }

    std::vector<char> digits;
};
