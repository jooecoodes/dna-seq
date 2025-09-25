#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <omp.h>
#include <chrono>
#include <algorithm>
#include <cstdint>
#include <windows.h>
#include <psapi.h>

// Get memory usage in KB
size_t getMemoryUsageKB() {
    PROCESS_MEMORY_COUNTERS_EX pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(),
                             (PROCESS_MEMORY_COUNTERS*)&pmc,
                             sizeof(pmc))) {
        return pmc.WorkingSetSize / 1024; // Current memory in KB
    }
    return 0;
}

// --- Shift-Or Serial ---
int shiftOrSearchSerial(const std::string& text, const std::string& pattern) {
    int n = text.length();
    int m = pattern.length();
    if (m == 0 || n < m) return 0;
    if (m > 64) return 0; // Shift-Or limit

    uint64_t B[256];
    for (int i = 0; i < 256; ++i) B[i] = ~0ULL;
    for (int i = 0; i < m; ++i) {
        B[static_cast<unsigned char>(pattern[i])] &= ~(1ULL << i);
    }

    uint64_t state = ~0ULL;
    int count = 0;
    
    for (int i = 0; i < n; ++i) {
        state = (state << 1) | B[static_cast<unsigned char>(text[i])];
        if (i >= m - 1 && (state & (1ULL << (m - 1))) == 0) count++;
    }
    return count;
}

// --- Shift-Or Parallel ---
int shiftOrSearchParallel(const std::string& text, const std::string& pattern, int num_threads) {
    int n = text.length();
    int m = pattern.length();
    if (m == 0 || n < m) return 0;
    if (m > 64) return 0;

    uint64_t B[256];
    for (int i = 0; i < 256; ++i) B[i] = ~0ULL;
    for (int i = 0; i < m; ++i) {
        B[static_cast<unsigned char>(pattern[i])] &= ~(1ULL << i);
    }

    int total_count = 0;

    #pragma omp parallel num_threads(num_threads) reduction(+:total_count)
    {
        int tid = omp_get_thread_num();
        int chunk_size = (n + num_threads - 1) / num_threads;
        int start = std::max(0, tid * chunk_size - (m - 1));
        int end = std::min(n, (tid + 1) * chunk_size);
        
        uint64_t state = ~0ULL;
        int local_count = 0;
        
        for (int i = start; i < end; ++i) {
            state = (state << 1) | B[static_cast<unsigned char>(text[i])];
            if (i >= m - 1 && i - (m - 1) >= tid * chunk_size && (state & (1ULL << (m - 1))) == 0) {
                local_count++;
            }
        }
        total_count += local_count;
    }
    return total_count;
}

// --- Read genome ---
std::string readGenome(const std::string& filename) {
    std::ifstream file(filename);
    std::string line, genome;
    if (!file.is_open()) return "";
    std::getline(file, line); // skip header
    while (std::getline(file, line)) genome += line;
    return genome;
}

int main() {
    std::string genome = readGenome("../dna/ecoli.fasta");
    if (genome.empty()) {
        std::cerr << "Failed to read genome file." << std::endl;
        return 1;
    }
    std::cout << "Genome length: " << genome.length() << " bp" << std::endl;

    // Benchmark metadata arrays
    std::vector<int> lengths = {64, 128, 256, 512, 1000, 2000};
    std::vector<double> gc_contents = {0.2, 0.5, 0.8};
    std::vector<double> entropies = {0.3, 1.1, 1.9};

    std::ifstream pfile("patterns.txt");
    if (!pfile.is_open()) {
        std::cerr << "Error: could not open patterns.txt" << std::endl;
        return 1;
    }

    std::ofstream csv("../benchmarks/algo_results/bp_results.csv");
    csv << "length,gc_content,entropy,matches,"
        << "serial_count,serial_time,serial_mem,"
        << "parallel_count,parallel_time,parallel_mem,"
        << "speedup,efficiency,overhead\n";

    std::string pattern;
    int idx = 0;

    while (std::getline(pfile, pattern)) {
        if (pattern.empty()) continue;

        // Map index -> metadata
        int len_idx = (idx / (gc_contents.size() * entropies.size())) % lengths.size();
        int gc_idx  = (idx / entropies.size()) % gc_contents.size();
        int h_idx   = idx % entropies.size();

        int len = lengths[len_idx];
        double gc = gc_contents[gc_idx];
        double h  = entropies[h_idx];

        idx++;

        if (pattern.size() > 64) continue; // Shift-Or limit

        std::cout << "\n=== Pattern " << idx 
                  << " | len=" << len 
                  << " GC=" << gc 
                  << " H=" << h << " ===" << std::endl;

        // Serial
        auto start = std::chrono::high_resolution_clock::now();
        int serial_count = shiftOrSearchSerial(genome, pattern);
        auto end = std::chrono::high_resolution_clock::now();
        auto serial_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        size_t serial_mem = getMemoryUsageKB();

        // Parallel (4 threads fixed)
        int num_threads = 4;
        start = std::chrono::high_resolution_clock::now();
        int parallel_count = shiftOrSearchParallel(genome, pattern, num_threads);
        end = std::chrono::high_resolution_clock::now();
        auto parallel_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        size_t parallel_mem = getMemoryUsageKB();

        double speedup = static_cast<double>(serial_time) / parallel_time;
        double efficiency = speedup / num_threads * 100;
        double overhead = parallel_time - (serial_time / num_threads);

        // "matches" = reference count (serial result)
        int matches = serial_count;

        // Write CSV row
        csv << len << "," << gc << "," << h << "," << matches << ","
            << serial_count << "," << serial_time << "," << serial_mem << ","
            << parallel_count << "," << parallel_time << "," << parallel_mem << ","
            << speedup << "," << efficiency << "," << overhead << "\n";
    }

    csv.close();
    return 0;
}
