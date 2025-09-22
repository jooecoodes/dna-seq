#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <omp.h>
#include <chrono>
#include <algorithm>

// Preprocess pattern for bad character rule (Boyer-Moore-Horspool)
std::vector<int> computeBadChar(const std::string& pattern) {
    std::vector<int> badChar(256, pattern.length());
    for (int i = 0; i < pattern.length() - 1; ++i) {
        badChar[static_cast<unsigned char>(pattern[i])] = pattern.length() - 1 - i;
    }
    return badChar;
}

// Serial Boyer-Moore-Horspool implementation
int bmhSearchSerial(const std::string& text, const std::string& pattern) {
    int n = text.length();
    int m = pattern.length();
    if (m == 0 || n < m) return 0;
    
    std::vector<int> badChar = computeBadChar(pattern);
    int count = 0;
    int s = 0;

    while (s <= n - m) {
        int j = m - 1;
        while (j >= 0 && pattern[j] == text[s + j]) {
            j--;
        }
        if (j < 0) {
            count++;
        }
        s += badChar[static_cast<unsigned char>(text[s + m - 1])];
    }
    return count;
}

// Parallel Boyer-Moore-Horspool implementation
int bmhSearchParallel(const std::string& text, const std::string& pattern, int num_threads) {
    int n = text.length();
    int m = pattern.length();
    if (m == 0 || n < m) return 0;
    
    std::vector<int> badChar = computeBadChar(pattern);
    int total_count = 0;

    #pragma omp parallel num_threads(num_threads) reduction(+:total_count)
    {
        int tid = omp_get_thread_num();
        int chunk_size = (n + num_threads - 1) / num_threads;
        int start = tid * chunk_size;
        int end = std::min(start + chunk_size + m - 1, n);
        
        int local_count = 0;
        int s = start;
        while (s <= end - m) {
            int j = m - 1;
            while (j >= 0 && pattern[j] == text[s + j]) {
                j--;
            }
            if (j < 0) {
                local_count++;
            }
            s += badChar[static_cast<unsigned char>(text[s + m - 1])];
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
    
    std::cout << "Genome length: " << genome.length() << " bp" << std::endl;

    // Open the patterns file
    std::ifstream pfile("patterns.txt");
    if (!pfile.is_open()) {
        std::cerr << "Error: could not open patterns.txt" << std::endl;
        return 1;
    }

    std::string pattern;
    while (std::getline(pfile, pattern)) {
        if (pattern.empty()) continue;

        std::cout << "\n=== Testing pattern: " << pattern << " ===" << std::endl;

        // Serial execution
        auto start = std::chrono::high_resolution_clock::now();
        int serial_count = bmhSearchSerial(genome, pattern);
        auto end = std::chrono::high_resolution_clock::now();
        auto serial_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
        std::cout << "Serial execution:" << std::endl;
        std::cout << "  Matches found: " << serial_count << std::endl;
        std::cout << "  Time: " << serial_time << " ms" << std::endl;
        
        // Parallel execution with different thread counts
        for (int num_threads = 2; num_threads <= 8; num_threads *= 2) {
            start = std::chrono::high_resolution_clock::now();
            int parallel_count = bmhSearchParallel(genome, pattern, num_threads);
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
    }

    return 0;
}
