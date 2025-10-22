#include "../../include/KMP.hpp"
#include "../../include/FastaReader.hpp"
#include "../../include/BioUtils.hpp"

#include <string>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <fstream>
#include <chrono>
#include <omp.h>
#include <cctype>
using namespace std;

// Minimum characters per thread to justify parallelism (tuneable)
static constexpr size_t MIN_PER_THREAD = 1 << 16; // 64k

std::vector<size_t> KMP::computeLPS(const std::string& pattern) const {
    size_t m = pattern.size();
    std::vector<size_t> lps(m, 0);
    size_t len = 0;
    for (size_t i = 1; i < m; ++i) {
        while (len > 0 && pattern[i] != pattern[len]) len = lps[len - 1];
        if (pattern[i] == pattern[len]) ++len;
        lps[i] = len;
    }
    return lps;
}

size_t KMP::search(const string& pattern, const string& text) const {
    const size_t n = text.size(), m = pattern.size();
    if (m == 0 || n < m) return 0;

    std::vector<size_t> lps = computeLPS(pattern);
    size_t j = 0;
    size_t count = 0;
    for (size_t i = 0; i < n; ++i) {
        while (j > 0 && pattern[j] != text[i]) j = lps[j - 1];
        if (pattern[j] == text[i]) ++j;
        if (j == m) {
            ++count;
            j = lps[j - 1];
        }
    }
    return count;
}

size_t KMP::searchInFasta(const string& pattern, const string& fastaPath) const {
    string dnaSequence = FastaReader::readSequence(fastaPath);
    return search(pattern, dnaSequence);
}

size_t KMP::searchParallel(const std::string& pattern, const std::string& text, int num_threads) const {
     const size_t n = text.size(), m = pattern.size();
    if (m == 0 || n < m) return 0;
    if (num_threads <= 0) num_threads = 1;
    if (n < static_cast<size_t>(num_threads) * MIN_PER_THREAD) num_threads = 1;

    std::vector<size_t> lps = computeLPS(pattern);
    size_t total_count = 0;

    #pragma omp parallel num_threads(num_threads) reduction(+: total_count)
    {
        int tid = omp_get_thread_num();
        size_t chunk = (n + num_threads - 1) / num_threads;
        size_t worker_start = tid * chunk;
        size_t worker_end = std::min(n, (tid + 1) * chunk);
        if (worker_start >= worker_end) {}

        size_t j = 0;
        size_t prefix_from = (worker_start >= (m - 1)) ? (worker_start - (m - 1)) : 0;
        for (size_t k = prefix_from; k < worker_start; ++k) {
            while (j > 0 && pattern[j] != text[k]) j = lps[j - 1];
            if (pattern[j] == text[k]) ++j;
        }

        size_t local_count = 0;
        for (size_t i = worker_start; i < std::min(n, worker_end + (m - 1)); ++i) {
            while (j > 0 && pattern[j] != text[i]) j = lps[j - 1];
            if (pattern[j] == text[i]) ++j;
            if (j == m) {
                size_t pos = i + 1 - m;
                if (pos >= worker_start && pos < worker_end) ++local_count;
                j = lps[j - 1];
            }
        }
        total_count += local_count;
    }
    return total_count;
}

size_t KMP::searchParallelInFasta(const std::string& pattern, const std::string& fastaPath) const {
    int num_threads = 4; // Default to 4 threads
    string dnaSequence = FastaReader::readSequence(fastaPath);
    return searchParallel(pattern, dnaSequence, num_threads);
}

size_t KMP::searchWithReverseComplement(const std::string& pattern, const std::string& text, bool parallel) const {
    std::string rc_pattern = BioUtils::reverseComplement(pattern);
    
    if (parallel) {
        int num_threads = 4;
        // Run searches sequentially but each uses internal parallelism
        size_t count_pattern = searchParallel(pattern, text, num_threads);
        size_t count_rc_pattern = searchParallel(rc_pattern, text, num_threads);
        return count_pattern + count_rc_pattern;
    } else {
        // Sequential version
        size_t count_pattern = search(pattern, text);
        size_t count_rc_pattern = search(rc_pattern, text);
        return count_pattern + count_rc_pattern;
    }
}
