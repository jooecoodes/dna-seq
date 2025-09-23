#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <omp.h>
#include <chrono>
#include <algorithm>
#include <cstdint>
#define NOMINMAX // to not include the old min max func from windows.h that collides with C++’s proper std::min function template.
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

// Bit-parallel Shift-Or algorithm implementation
int shiftOrSearchSerial(const std::string& text, const std::string& pattern) {
    int n = text.length();
    int m = pattern.length();
    if (m == 0 || n < m) return 0;
    if (m > 64) {
        std::cerr << "Pattern too long for Shift-Or algorithm." << std::endl;
        return 0;
    }

    uint64_t B[256];
    for (int i = 0; i < 256; ++i) {
        B[i] = ~0ULL;
    }
    for (int i = 0; i < m; ++i) {
        B[static_cast<unsigned char>(pattern[i])] &= ~(1ULL << i);
    }

    uint64_t state = ~0ULL;
    int count = 0;
    
    for (int i = 0; i < n; ++i) {
        state = (state << 1) | B[static_cast<unsigned char>(text[i])];
        if (i >= m - 1 && (state & (1ULL << (m - 1))) == 0) {
            count++;
        }
    }
    
    return count;
}

// Parallel Shift-Or implementation
int shiftOrSearchParallel(const std::string& text, const std::string& pattern, int num_threads) {
    int n = text.length();
    int m = pattern.length();
    if (m == 0 || n < m) return 0;
    if (m > 64) {
        std::cerr << "Pattern too long for Shift-Or algorithm." << std::endl;
        return 0;
    }

    uint64_t B[256];
    for (int i = 0; i < 256; ++i) {
        B[i] = ~0ULL;
    }
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

// Function to read genome from FASTA file
std::string readGenome(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;
    std::string genome;
    
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return "";
    }
    
    std::getline(file, line);
    while (std::getline(file, line)) {
        genome += line;
    }
    
    file.close();
    return genome;
}

int main() {
    std::string genome = readGenome("ecoli.fasta");
    if (genome.empty()) {
        std::cerr << "Failed to read genome file." << std::endl;
        return 1;
    }
    
    std::cout << "Genome length: " << genome.length() << " bp" << std::endl;

    std::ifstream pfile("patterns.txt");
    if (!pfile.is_open()) {
        std::cerr << "Error: could not open patterns.txt" << std::endl;
        return 1;
    }

    std::string pattern;
    while (std::getline(pfile, pattern)) {
        if (pattern.empty()) continue;

        std::cout << "\n=== Testing pattern: " << pattern << " ===" << std::endl;

        auto start = std::chrono::high_resolution_clock::now();
        int serial_count = shiftOrSearchSerial(genome, pattern);
        auto end = std::chrono::high_resolution_clock::now();
        auto serial_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        size_t serial_mem = getMemoryUsageKB();
        
        std::cout << "Serial execution:" << std::endl;
        std::cout << "  Matches found: " << serial_count << std::endl;
        std::cout << "  Time: " << serial_time << " ms" << std::endl;
        std::cout << "  Memory: " << serial_mem << " KB" << std::endl;
        
        for (int num_threads = 2; num_threads <= 8; num_threads *= 2) {
            start = std::chrono::high_resolution_clock::now();
            int parallel_count = shiftOrSearchParallel(genome, pattern, num_threads);
            end = std::chrono::high_resolution_clock::now();
            auto parallel_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            size_t parallel_mem = getMemoryUsageKB();
            
            double speedup = static_cast<double>(serial_time) / parallel_time;
            double efficiency = speedup / num_threads * 100;
            double overhead = parallel_time - (serial_time / num_threads);
            
            std::cout << "\nParallel execution with " << num_threads << " threads:" << std::endl;
            std::cout << "  Matches found: " << parallel_count << std::endl;
            std::cout << "  Time: " << parallel_time << " ms" << std::endl;
            std::cout << "  Memory: " << parallel_mem << " KB" << std::endl;
            std::cout << "  Speedup: " << speedup << "x" << std::endl;
            std::cout << "  Efficiency: " << efficiency << "%" << std::endl;
            std::cout << "  Overhead: " << overhead << " ms" << std::endl;
        }
    }

    return 0;
}
