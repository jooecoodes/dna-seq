#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <omp.h>
#include <chrono>

// Preprocess pattern to create LPS array
std::vector<int> computeLPS(const std::string& pattern) {
    int m = pattern.length();
    std::vector<int> lps(m, 0);
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
    return lps;
}

// Serial KMP implementation
int kmpSearchSerial(const std::string& text, const std::string& pattern) {
    int n = text.length();
    int m = pattern.length();
    std::vector<int> lps = computeLPS(pattern);
    int count = 0;
    int i = 0, j = 0;
    
    while (i < n) {
        if (pattern[j] == text[i]) {
            i++;
            j++;
        }
        
        if (j == m) {
            count++;
            j = lps[j - 1];
        } else if (i < n && pattern[j] != text[i]) {
            if (j != 0) {
                j = lps[j - 1];
            } else {
                i++;
            }
        }
    }
    return count;
}

// Parallel KMP implementation
int kmpSearchParallel(const std::string& text, const std::string& pattern, int num_threads) {
    int n = text.length();
    int m = pattern.length();
    std::vector<int> lps = computeLPS(pattern);
    int total_count = 0;
    
    #pragma omp parallel num_threads(num_threads) reduction(+:total_count)
    {
        int local_count = 0;
        int tid = omp_get_thread_num();
        int chunk_size = n / num_threads;
        int start = tid * chunk_size;
        int end = (tid == num_threads - 1) ? n : start + chunk_size;
        
        if (tid > 0) {
            start -= (m - 1);
            if (start < 0) start = 0;
        }
        
        int i = start;
        int j = 0;
        
        while (i < end) {
            if (pattern[j] == text[i]) {
                i++;
                j++;
            }
            
            if (j == m) {
                local_count++;
                j = lps[j - 1];
            } else if (i < end && pattern[j] != text[i]) {
                if (j != 0) {
                    j = lps[j - 1];
                } else {
                    i++;
                }
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
        int serial_count = kmpSearchSerial(genome, pattern);
        auto end = std::chrono::high_resolution_clock::now();
        auto serial_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
        std::cout << "Serial execution:" << std::endl;
        std::cout << "  Matches found: " << serial_count << std::endl;
        std::cout << "  Time: " << serial_time << " ms" << std::endl;
        
        for (int num_threads = 2; num_threads <= 8; num_threads *= 2) {
            start = std::chrono::high_resolution_clock::now();
            int parallel_count = kmpSearchParallel(genome, pattern, num_threads);
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
