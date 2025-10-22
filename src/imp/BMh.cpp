#include "../../include/BM.hpp"
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

std::vector<size_t> BoyerMooreHorspool::createBadCharTable(const std::string& pattern) const {
    std::vector<size_t> table(256, pattern.size());
    size_t m = pattern.size();
    for (size_t i = 0; i + 1 < m; ++i)
        table[(unsigned char)pattern[i]] = m - 1 - i;
    return table;
}

size_t BoyerMooreHorspool::search(const string& pattern, const string& text) const {
    const size_t n = text.size(), m = pattern.size();
    if (m == 0 || n < m) return 0;

    std::vector<size_t> badChar = createBadCharTable(pattern);
    size_t count = 0;
    size_t s = 0;
    while (s + m <= n) {
        size_t j = m;
        while (j > 0 && pattern[j - 1] == text[s + j - 1]) --j;
        if (j == 0) {
            ++count;
            ++s;
        } else {
            unsigned char mc = static_cast<unsigned char>(text[s + m - 1]);
            size_t shift = badChar[mc];
            if (shift == 0) shift = 1;
            s += shift;
        }
    }
    return count;
}

size_t BoyerMooreHorspool::searchInFasta(const string& pattern, const string& fastaPath) const {
    string dnaSequence = FastaReader::readSequence(fastaPath);
    return search(pattern, dnaSequence);
}

size_t BoyerMooreHorspool::searchParallel(const std::string& pattern, const std::string& text, int num_threads) const {
    const size_t n = text.size(), m = pattern.size();
    if (m == 0 || n < m) return 0;
    if (num_threads <= 0) num_threads = 1;
    if (n < static_cast<size_t>(num_threads) * MIN_PER_THREAD) num_threads = 1;

    std::vector<size_t> badChar = createBadCharTable(pattern);
    size_t total_count = 0;

    #pragma omp parallel num_threads(num_threads) reduction(+: total_count)
    {
        int tid = omp_get_thread_num();
        size_t chunk = (n + num_threads - 1) / num_threads;
        size_t worker_start = tid * chunk;
        size_t worker_end = std::min(n, (tid + 1) * chunk);
        if (worker_start >= worker_end) {}

        size_t local_count = 0;
        size_t s = worker_start;
        while (s + m <= n && s < worker_end) {
            size_t j = m;
            while (j > 0 && pattern[j - 1] == text[s + j - 1]) --j;
            if (j == 0) {
                ++local_count;
                ++s;
            } else {
                unsigned char mc = static_cast<unsigned char>(text[s + m - 1]);
                size_t shift = badChar[mc];
                if (shift == 0) shift = 1;
                s += shift;
            }
        }
        total_count += local_count;
    }
    return total_count;
}
size_t BoyerMooreHorspool::searchParallelInFasta(const std::string& pattern, const std::string& fastaPath) const {
    int num_threads = 4;
    std::string dnaSequence = FastaReader::readSequence(fastaPath);
    return searchParallel(pattern, dnaSequence, num_threads);
}

size_t BoyerMooreHorspool::searchWithReverseComplement(const std::string& pattern, const std::string& text, bool parallel) const {
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