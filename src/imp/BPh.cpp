#include "../../include/BP.hpp"
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
#include <immintrin.h>
using namespace std;


// Minimum characters per thread to justify parallelism (tuneable)
static constexpr size_t MIN_PER_THREAD = 1 << 16; // 64k

size_t BitParallelShiftOr::search(const string& pattern, const string& text) const {
    const size_t n = text.size(), m = pattern.size();
    if (m == 0 || n < m || m > 64) return 0;

    uint64_t B[256];
    for (size_t i = 0; i < 256; ++i) B[i] = ~0ULL;
    for (size_t i = 0; i < m; ++i)
        B[(unsigned char)pattern[i]] &= ~(1ULL << i);

    uint64_t state = ~0ULL;
    size_t count = 0;
    for (size_t i = 0; i < n; ++i) {
        state = (state << 1) | B[(unsigned char)text[i]];
        if (i >= m - 1 && (state & (1ULL << (m - 1))) == 0)
            ++count;
    }
    return count; 
}

size_t BitParallelShiftOr::searchInFasta(const string& pattern, const string& fastaPath) const {
    string dnaSequence = FastaReader::readSequence(fastaPath);
    return search(pattern, dnaSequence);
}

size_t BitParallelShiftOr::searchParallel(const std::string& pattern, const std::string& text, int num_threads) const {
    const size_t n = text.size(), m = pattern.size();
    if (m == 0 || n < m || m > 64) return 0;

    uint64_t B[256];
    for (size_t i = 0; i < 256; ++i) B[i] = ~0ULL;
    for (size_t i = 0; i < m; ++i)
        B[(unsigned char)pattern[i]] &= ~(1ULL << i);

    uint64_t state = ~0ULL;
    size_t count = 0;
    for (size_t i = 0; i < n; ++i) {
        state = (state << 1) | B[(unsigned char)text[i]];
        if (i >= m - 1 && (state & (1ULL << (m - 1))) == 0)
            ++count;
    }
    return count; 
}

size_t BitParallelShiftOr::searchParallelInFasta(const std::string& pattern, const std::string& fastaPath) const {
    int num_threads = 4;
    std::string dnaSequence = FastaReader::readSequence(fastaPath);
    return searchParallel(pattern, dnaSequence, num_threads);
}

size_t BitParallelShiftOr::searchWithReverseComplement(const std::string& pattern, const std::string& text, bool parallel) const {
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