#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <omp.h>
#include <chrono>
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
    if (m == 0 || n < m) return 0;

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
    if (m == 0 || n < m) return 0;

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

// Read genome
std::string readGenome(const std::string& filename) {
    std::ifstream file(filename);
    std::string line, genome;
    if (!file.is_open()) return "";
    std::getline(file, line); // skip header
    while (std::getline(file, line)) genome += line;
    return genome;
}

int main() {
    std::string genome = readGenome("ecoli.fasta");
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

    std::ofstream csv("../benchmarks/algo_results/kmp_results.csv");
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

        std::cout << "\n=== Pattern " << idx 
                  << " | len=" << len 
                  << " GC=" << gc 
                  << " H=" << h << " ===" << std::endl;

        // Serial
        auto start = std::chrono::high_resolution_clock::now();
        int serial_count = kmpSearchSerial(genome, pattern);
        auto end = std::chrono::high_resolution_clock::now();
        auto serial_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        size_t serial_mem = getMemoryUsageKB();

        // Parallel (4 threads fixed)
        int num_threads = 4;
        start = std::chrono::high_resolution_clock::now();
        int parallel_count = kmpSearchParallel(genome, pattern, num_threads);
        end = std::chrono::high_resolution_clock::now();
        auto parallel_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        size_t parallel_mem = getMemoryUsageKB();

        double speedup = static_cast<double>(serial_time) / parallel_time;
        double efficiency = speedup / num_threads * 100;
        double overhead = parallel_time - (serial_time / num_threads);

        int matches = serial_count; // reference

        // Write CSV row
        csv << len << "," << gc << "," << h << "," << matches << ","
            << serial_count << "," << serial_time << "," << serial_mem << ","
            << parallel_count << "," << parallel_time << "," << parallel_mem << ","
            << speedup << "," << efficiency << "," << overhead << "\n";
    }

    csv.close();
    return 0;
}
