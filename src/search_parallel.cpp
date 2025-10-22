// file: search_parallel.cpp
// Compile with: g++ -O3 -std=c++17 -fopenmp search_parallel.cpp -o search_parallel

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

static constexpr size_t MIN_PER_THREAD = 1 << 16; // 64k

// ---------- FASTA Reader ----------
std::string load_clean_fasta(const std::string &path) {
    std::ifstream f(path);
    if (!f) {
        std::cerr << "Error: cannot open " << path << "\n";
        std::exit(1);
    }
    std::string buf, seq;
    seq.reserve(10'000'000);
    while (std::getline(f, buf)) {
        if (!buf.empty() && buf[0] == '>') continue;
        for (char c : buf) {
            if (c != '\r' && c != '\n') {
                c = std::toupper(static_cast<unsigned char>(c));
                if (c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N')
                    seq.push_back(c);
            }
        }
    }
    return seq;
}

// ---------- Helpers ----------
static std::vector<size_t> createBadCharTable(const std::string &pat) {
    std::vector<size_t> table(256, pat.size());
    size_t m = pat.size();
    for (size_t i = 0; i + 1 < m; ++i)
        table[(unsigned char)pat[i]] = m - 1 - i;
    return table;
}

static std::vector<size_t> computeLPS(const std::string &pat) {
    size_t m = pat.size();
    std::vector<size_t> lps(m, 0);
    size_t len = 0;
    for (size_t i = 1; i < m; ++i) {
        while (len > 0 && pat[i] != pat[len]) len = lps[len - 1];
        if (pat[i] == pat[len]) ++len;
        lps[i] = len;
    }
    return lps;
}

// ---------- Sequential Implementations ----------

size_t BoyerMooreHorspool_searchSequential(const std::string &pattern, const std::string &text) {
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

size_t BitParallelShiftOr_searchSequential(const std::string &pattern, const std::string &text) {
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

size_t KMP_searchSequential(const std::string &pattern, const std::string &text) {
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

// ---------- Parallel Implementations (Existing) ----------

size_t BoyerMooreHorspool_searchParallel(const std::string &pattern, const std::string &text, int num_threads) {
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

size_t BitParallelShiftOr_searchParallel(const std::string &pattern, const std::string &text, int num_threads) {
    const size_t n = text.size(), m = pattern.size();
    if (m == 0 || n < m || m > 64) return 0;
    if (num_threads <= 0) num_threads = 1;
    if (n < static_cast<size_t>(num_threads) * MIN_PER_THREAD) num_threads = 1;

    uint64_t B_global[256];
    for (size_t i = 0; i < 256; ++i) B_global[i] = ~0ULL;
    for (size_t i = 0; i < m; ++i)
        B_global[(unsigned char)pattern[i]] &= ~(1ULL << i);

    size_t total_count = 0;
    #pragma omp parallel num_threads(num_threads) reduction(+: total_count)
    {
        uint64_t B[256];
        std::memcpy(B, B_global, sizeof(B_global));

        int tid = omp_get_thread_num();
        size_t chunk = (n + num_threads - 1) / num_threads;
        size_t worker_start = tid * chunk;
        size_t worker_end = std::min(n, (tid + 1) * chunk);
        if (worker_start >= worker_end) {}

        uint64_t state = ~0ULL;
        size_t prefix_from = (worker_start >= (m - 1)) ? (worker_start - (m - 1)) : 0;
        for (size_t k = prefix_from; k < worker_start; ++k)
            state = (state << 1) | B[(unsigned char)text[k]];

        size_t local_count = 0;
        for (size_t i = worker_start; i < std::min(n, worker_end + (m - 1)); ++i) {
            state = (state << 1) | B[(unsigned char)text[i]];
            if (i >= m - 1) {
                size_t pos = i - (m - 1);
                if ((state & (1ULL << (m - 1))) == 0)
                    if (pos >= worker_start && pos < worker_end) ++local_count;
            }
        }
        total_count += local_count;
    }
    return total_count;
}

size_t KMP_searchParallel(const std::string &pattern, const std::string &text, int num_threads) {
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

// ---------- Main ----------
int main(int argc, char **argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <fasta_file> <pattern> <threads>\n";
        return 1;
    }

    std::string fasta_path = argv[1];
    std::string pattern = argv[2];
    int threads = std::stoi(argv[3]);

    std::cout << "Loading FASTA file: " << fasta_path << " ...\n";
    auto text = load_clean_fasta(fasta_path);
    std::cout << "Sequence length: " << text.size() << " bases\n";

    auto run = [&](auto func, const std::string &name) {
        auto t1 = std::chrono::high_resolution_clock::now();
        size_t count = func(pattern, text);
        auto t2 = std::chrono::high_resolution_clock::now();
        double secs = std::chrono::duration<double>(t2 - t1).count();
        std::cout << name << ": " << count << " matches in " << secs << " s\n";
    };

    auto runParallel = [&](auto func, const std::string &name) {
        auto t1 = std::chrono::high_resolution_clock::now();
        size_t count = func(pattern, text, threads);
        auto t2 = std::chrono::high_resolution_clock::now();
        double secs = std::chrono::duration<double>(t2 - t1).count();
        std::cout << name << ": " << count << " matches in " << secs << " s\n";
    };

    std::cout << "\n=== Sequential Algorithms ===\n";
    run(BoyerMooreHorspool_searchSequential, "Boyer–Moore–Horspool (seq)");
    run(BitParallelShiftOr_searchSequential,  "Bit-Parallel Shift-Or (seq)");
    run(KMP_searchSequential,                 "KMP (seq)");

    std::cout << "\n=== Parallel Algorithms ===\n";
    runParallel(BoyerMooreHorspool_searchParallel, "Boyer–Moore–Horspool (parallel)");
    runParallel(BitParallelShiftOr_searchParallel,  "Bit-Parallel Shift-Or (parallel)");
    runParallel(KMP_searchParallel,                 "KMP (parallel)");

    return 0;
}
