#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <omp.h>
#include <chrono>
#include <algorithm>
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

// Preprocess pattern for bad character rule (BMH)
std::vector<int> computeBadChar(const std::string& pattern) {
    std::vector<int> badChar(256, pattern.length());
    for (int i = 0; i < pattern.length() - 1; ++i) {
        badChar[static_cast<unsigned char>(pattern[i])] = pattern.length() - 1 - i;
    }
    return badChar;
}

// Serial Boyer-Moore-Horspool
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

// Parallel Boyer-Moore-Horspool
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

// Read genome
std::string readGenome(const std::string& filename) {
    std::ifstream file(filename);
    std::string line, genome;
    if (!file.is_open()) return "";
    std::getline(file, line); // skip header
    while (std::getline(file, line)) genome += line;
    return genome;
}

int main(int argc, char* argv[]) {
    std::string genome_file = (argc > 1) ? argv[1] : "../../dna/ecoli.fasta";
    std::string pattern_file = (argc > 2) ? argv[2] : "../patterns.txt";
    std::string output_csv = (argc > 3) ? argv[3] : "../../benchmarks/algo_results/bmh_results.csv";
    std::string genome = readGenome(genome_file);
    
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

    std::ofstream csv(output_csv);
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
        int serial_count = bmhSearchSerial(genome, pattern);
        auto end = std::chrono::high_resolution_clock::now();
        auto serial_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        size_t serial_mem = getMemoryUsageKB();

        // Parallel (4 threads fixed)
        int num_threads = 4;
        start = std::chrono::high_resolution_clock::now();
        int parallel_count = bmhSearchParallel(genome, pattern, num_threads);
        end = std::chrono::high_resolution_clock::now();
        auto parallel_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        size_t parallel_mem = getMemoryUsageKB();

        double speedup = 0.0;
        double efficiency = 0.0;
        double overhead = 0.0;

        if (parallel_time > 0) {
            speedup = static_cast<double>(serial_time) / parallel_time;
            efficiency = speedup / num_threads * 100.0;
            overhead = parallel_time - (serial_time / num_threads);
        }

        int matches = serial_count; // ground truth matches

        // Write CSV row
        csv << len << "," << gc << "," << h << "," << matches << ","
            << serial_count << "," << serial_time << "," << serial_mem << ","
            << parallel_count << "," << parallel_time << "," << parallel_mem << ","
            << speedup << "," << efficiency << "," << overhead << "\n";
    }

    csv.close();
    return 0;
}
