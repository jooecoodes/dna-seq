#include <vector>
#include <string>
#include <algorithm>
#include <cstring>
#include <unordered_map>

using namespace std;

// KMP Algorithm
vector<int> kmp_search(const string& text, const string& pattern) {
    int n = text.length();
    int m = pattern.length();
    vector<int> matches;
    if (m == 0) return matches;

    vector<int> lps(m, 0);
    int len = 0;
    for (int i = 1; i < m; ) {
        if (pattern[i] == pattern[len]) {
            lps[i++] = ++len;
        } else if (len != 0) {
            len = lps[len - 1];
        } else {
            lps[i++] = 0;
        }
    }

    int i = 0, j = 0;
    while (i < n) {
        if (pattern[j] == text[i]) {
            i++;
            j++;
        }
        if (j == m) {
            matches.push_back(i - j);
            j = lps[j - 1];
        } else if (i < n && pattern[j] != text[i]) {
            if (j != 0) j = lps[j - 1];
            else i++;
        }
    }
    return matches;
}

// Boyer-Moore Algorithm (Bad Character Heuristic)
vector<int> boyer_moore_search(const string& text, const string& pattern) {
    int n = text.length();
    int m = pattern.length();
    vector<int> matches;
    if (m == 0) return matches;

    unordered_map<char, int> bad_char_shift;
    for (int i = 0; i < m; i++) {
        bad_char_shift[pattern[i]] = m - i - 1;
    }

    int i = 0;
    while (i <= n - m) {
        int j = m - 1;
        while (j >= 0 && pattern[j] == text[i + j]) j--;
        if (j < 0) {
            matches.push_back(i);
            i += m;
        } else {
            i += max(1, bad_char_shift.count(text[i + j]) ? bad_char_shift[text[i + j]] : m);
        }
    }
    return matches;
}

// Bit-Parallel Algorithm (Shift-Or for exact matching)
vector<int> bit_parallel_search(const string& text, const string& pattern) {
    int m = pattern.length();
    int n = text.length();
    vector<int> matches;
    if (m == 0 || m > 64) return matches; // Limited to 64 bits

    unsigned long long pattern_mask[256] = {~0ULL};
    unsigned long long R = ~1ULL;

    for (int i = 0; i < m; i++) {
        pattern_mask[pattern[i]] &= ~(1ULL << i);
    }

    for (int i = 0; i < n; i++) {
        R |= pattern_mask[text[i]];
        R <<= 1;
        if ((R & (1ULL << m)) == 0) {
            matches.push_back(i - m + 1);
        }
    }
    return matches;
}

// Extern C functions for each algorithm
extern "C" {
    int* kmp_search_c(const char* text, const char* pattern, int* result_size) {
        string s_text(text);
        string s_pattern(pattern);
        vector<int> matches = kmp_search(s_text, s_pattern);
        *result_size = matches.size();
        int* result = new int[matches.size()];
        copy(matches.begin(), matches.end(), result);
        return result;
    }

    int* boyer_moore_search_c(const char* text, const char* pattern, int* result_size) {
        string s_text(text);
        string s_pattern(pattern);
        vector<int> matches = boyer_moore_search(s_text, s_pattern);
        *result_size = matches.size();
        int* result = new int[matches.size()];
        copy(matches.begin(), matches.end(), result);
        return result;
    }

    int* bit_parallel_search_c(const char* text, const char* pattern, int* result_size) {
        string s_text(text);
        string s_pattern(pattern);
        vector<int> matches = bit_parallel_search(s_text, s_pattern);
        *result_size = matches.size();
        int* result = new int[matches.size()];
        copy(matches.begin(), matches.end(), result);
        return result;
    }

    void free_result(int* result) {
        delete[] result;
    }
}