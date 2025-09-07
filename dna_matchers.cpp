#include <vector>
#include <string>
#include <algorithm>
#include <cstring>
#include <unordered_map>
#include <cmath>
#include <array>
#include <bitset>

using namespace std;

// Ultra-optimized KMP with precomputation optimization
vector<int> kmp_search(const string& text, const string& pattern) {
    int n = text.length();
    int m = pattern.length();
    vector<int> matches;
    if (m == 0) return matches;

    // Precompute LPS array with single pass
    vector<int> lps(m, 0);
    int len = 0;
    int i = 1;
    
    while (i < m) {
        if (pattern[i] == pattern[len]) {
            len++;
            lps[i] = len;
            i++;
        } else {
            if (len != 0) {
                len = lps[len - 1];
            } else {
                lps[i] = 0;
                i++;
            }
        }
    }

    // Optimized search with minimal branching
    i = 0;
    int j = 0;
    while (i < n) {
        while (j < m && i < n && pattern[j] == text[i]) {
            j++;
            i++;
        }
        
        if (j == m) {
            matches.push_back(i - j);
            j = lps[j - 1];
        } else if (i < n) {
            if (j != 0) {
                j = lps[j - 1];
            } else {
                i++;
            }
        }
    }
    
    return matches;
}

// Highly optimized Boyer-Moore with combined heuristics
vector<int> boyer_moore_search(const string& text, const string& pattern) {
    int n = text.length();
    int m = pattern.length();
    vector<int> matches;
    if (m == 0) return matches;

    // Bad character heuristic with direct array access
    const int ALPHABET_SIZE = 256;
    array<int, ALPHABET_SIZE> badChar;
    fill(badChar.begin(), badChar.end(), m);
    
    for (int i = 0; i < m; i++) {
        badChar[(unsigned char)pattern[i]] = m - i - 1;
    }

    // Good suffix heuristic (simplified)
    vector<int> goodSuffix(m + 1, m);
    vector<int> borderPos(m + 1, 0);
    
    int i = m, j = m + 1;
    borderPos[i] = j;
    
    while (i > 0) {
        while (j <= m && pattern[i - 1] != pattern[j - 1]) {
            if (goodSuffix[j] == m) {
                goodSuffix[j] = j - i;
            }
            j = borderPos[j];
        }
        i--;
        j--;
        borderPos[i] = j;
    }
    
    j = borderPos[0];
    for (i = 0; i <= m; i++) {
        if (goodSuffix[i] == m) {
            goodSuffix[i] = j;
        }
        if (i == j) {
            j = borderPos[j];
        }
    }

    // Search with combined heuristics
    int s = 0;
    while (s <= n - m) {
        int j = m - 1;
        
        // Fast match from right to left
        while (j >= 0 && pattern[j] == text[s + j]) {
            j--;
        }
        
        if (j < 0) {
            matches.push_back(s);
            s += goodSuffix[0];
        } else {
            s += max(goodSuffix[j + 1], j - badChar[text[s + j]]);
        }
    }
    
    return matches;
}

// Optimized Bit-Parallel without register keyword
vector<int> bit_parallel_search(const string& text, const string& pattern) {
    int m = pattern.length();
    int n = text.length();
    vector<int> matches;
    if (m == 0 || m > 64) return matches;

    // Use array for pattern mask
    array<unsigned long long, 256> patternMask;
    for (int i = 0; i < 256; i++) {
        patternMask[i] = ~0ULL;
    }
    
    for (int i = 0; i < m; i++) {
        patternMask[(unsigned char)pattern[i]] &= ~(1ULL << i);
    }

    unsigned long long R = ~1ULL;
    unsigned long long mask = (1ULL << m);
    
    for (int i = 0; i < n; i++) {
        R = (R << 1) | patternMask[(unsigned char)text[i]];
        
        if ((R & mask) == 0) {
            matches.push_back(i - m + 1);
        }
    }
    
    return matches;
}

// Lightweight pattern analysis with minimal computation
struct PatternAnalysis {
    int length;
    float gcContent;
    bool isRepetitive;
    int maxRepeat;
};

PatternAnalysis analyze_pattern_fast(const string& pattern) {
    int length = pattern.length();
    PatternAnalysis result;
    result.length = length;
    
    if (length == 0) {
        result.gcContent = 0;
        result.isRepetitive = false;
        result.maxRepeat = 0;
        return result;
    }
    
    // Calculate GC content with single pass
    int gc_count = 0;
    int max_repeat = 1;
    int current_repeat = 1;
    char prev_char = '\0';
    
    for (int i = 0; i < length; i++) {
        char c = pattern[i];
        if (c == 'G' || c == 'C') {
            gc_count++;
        }
        
        if (c == prev_char) {
            current_repeat++;
            if (current_repeat > max_repeat) {
                max_repeat = current_repeat;
            }
        } else {
            current_repeat = 1;
        }
        
        prev_char = c;
    }
    
    result.gcContent = static_cast<float>(gc_count) / length;
    result.maxRepeat = max_repeat;
    result.isRepetitive = (max_repeat > length / 3) || (max_repeat > 10);
    
    return result;
}

// Decision tree based algorithm selection
vector<int> hybrid_search_cpp(const string& text, const string& pattern) {
    PatternAnalysis analysis = analyze_pattern_fast(pattern);
    int length = analysis.length;
    
    // Decision tree for algorithm selection
    if (length <= 16) {
        // Very short patterns: always bit-parallel
        return bit_parallel_search(text, pattern);
    } else if (length <= 64) {
        // Short patterns: bit-parallel for non-repetitive, KMP for repetitive
        if (!analysis.isRepetitive) {
            return bit_parallel_search(text, pattern);
        } else {
            return kmp_search(text, pattern);
        }
    } else if (length <= 256) {
        // Medium patterns: Boyer-Moore for non-repetitive, KMP for repetitive
        if (!analysis.isRepetitive) {
            return boyer_moore_search(text, pattern);
        } else {
            return kmp_search(text, pattern);
        }
    } else {
        // Long patterns: Boyer-Moore for non-repetitive/low GC, KMP for repetitive/high GC
        if (!analysis.isRepetitive && analysis.gcContent < 0.6) {
            return boyer_moore_search(text, pattern);
        } else {
            return kmp_search(text, pattern);
        }
    }
}

// Extern C functions
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

    int* hybrid_search_cpp_c(const char* text, const char* pattern, int* result_size) {
        string s_text(text);
        string s_pattern(pattern);
        vector<int> matches = hybrid_search_cpp(s_text, s_pattern);
        *result_size = matches.size();
        int* result = new int[matches.size()];
        copy(matches.begin(), matches.end(), result);
        return result;
    }

    void free_result(int* result) {
        delete[] result;
    }
}