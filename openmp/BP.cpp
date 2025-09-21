#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <omp.h>
#include <chrono>
#include <algorithm>
#include <cstdint>

// Bit-parallel Shift-Or algorithm implementation
int shiftOrSearchSerial(const std::string& text, const std::string& pattern) {
    int n = text.length();
    int m = pattern.length();
    if (m == 0 || n < m) return 0;
    if (m > 64) {
        std::cerr << "Pattern too long for Shift-Or algorithm." << std::endl;
        return 0;
    }

    // Preprocessing: create character masks
    uint64_t B[256];
    for (int i = 0; i < 256; ++i) {
        B[i] = ~0ULL; // Initialize all bits to 1
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

    // Preprocessing: create character masks (shared across threads)
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
        
        if (start < 0) start = 0;
        if (end > n) end = n;
        
        uint64_t state = ~0ULL;
        int local_count = 0;
        
        // Process the chunk
        for (int i = start; i < end; ++i) {
            state = (state << 1) | B[static_cast<unsigned char>(text[i])];
            
            // Only count matches that start within this thread's responsible area
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
    
    // Skip header line
    std::getline(file, line);
    
    // Read genome data
    while (std::getline(file, line)) {
        genome += line;
    }
    
    file.close();
    return genome;
}

int main() {
    // Read E. Coli genome
    std::string genome = readGenome("ecoli.fasta");
    if (genome.empty()) {
        std::cerr << "Failed to read genome file." << std::endl;
        return 1;
    }
    
    std::string pattern = "ATG"; // Start codon
    
    std::cout << "Genome length: " << genome.length() << " bp" << std::endl;
    std::cout << "Pattern: " << pattern << std::endl;
    
    // Serial bit-parallel execution
    auto start = std::chrono::high_resolution_clock::now();
    int serial_count = shiftOrSearchSerial(genome, pattern);
    auto end = std::chrono::high_resolution_clock::now();
    auto serial_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    
    std::cout << "Serial bit-parallel execution:" << std::endl;
    std::cout << "  Matches found: " << serial_count << std::endl;
    std::cout << "  Time: " << serial_time << " ms" << std::endl;
    
    // Parallel execution with different thread counts
    for (int num_threads = 2; num_threads <= 8; num_threads *= 2) {
        start = std::chrono::high_resolution_clock::now();
        int parallel_count = shiftOrSearchParallel(genome, pattern, num_threads);
        end = std::chrono::high_resolution_clock::now();
        auto parallel_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
        double speedup = static_cast<double>(serial_time) / parallel_time;
        double efficiency = speedup / num_threads * 100;
        double overhead = parallel_time - (serial_time / num_threads);
        
        std::cout << "\nParallel execution with " << num_threads << " threads:" << std::endl;
        std::cout << "  Matches found: " << parallel_count << std::endl;
        std::cout << "  Time: " << parallel_time << " ms" << std::endl;
        std::cout << "  Speedup: " << speedup << "x" << std::endl;
        std::cout << "  Efficiency: " << efficiency << "%" << std::endl;
        std::cout << "  Overhead: " << overhead << " ms" << std::endl;
    }
    
    return 0;
}