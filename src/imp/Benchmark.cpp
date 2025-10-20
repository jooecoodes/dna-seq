#include "../../include/Benchmark.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;
using namespace std::chrono;

long long Benchmark::getMemoryUsage() {
    // Linux specific memory usage check
    ifstream statm("/proc/self/statm");
    if (statm.is_open()) {
        long long size, resident, share, text, lib, data, dt;
        statm >> size >> resident >> share >> text >> lib >> data >> dt;
        statm.close();
        return resident * sysconf(_SC_PAGESIZE) / 1024; // Convert to KB
    }
    return 0;
}

Benchmark::BenchmarkResult Benchmark::runAlgorithm(const string& algorithmName, 
                                                  const string& pattern, 
                                                  const string& fastaPath) {
    BenchmarkResult result;
    result.algorithm_name = algorithmName;
    
    // Measure memory before
    size_t memory_before = getMemoryUsage();
    
    auto start = high_resolution_clock::now();
    
    // Run the algorithm
    HybridPicker picker;
    result.matches = picker.pickAndSearch(algorithmName, pattern, fastaPath);
    
    auto end = high_resolution_clock::now();
    
    // Measure memory after
    size_t memory_after = getMemoryUsage();
    
    result.time_ms = duration_cast<milliseconds>(end - start).count();
    result.memory_usage_kb = memory_after - memory_before;
    
    return result;
}

void Benchmark::printResults(const vector<BenchmarkResult>& results) {
    cout << "\n" << string(80, '=') << "\n";
    cout << "BENCHMARK RESULTS\n";
    cout << string(80, '=') << "\n";
    
    cout << left << setw(20) << "ALGORITHM" 
         << setw(15) << "TIME (ms)" 
         << setw(15) << "MEMORY (KB)" 
         << setw(15) << "MATCHES" 
         << "\n";
    cout << string(80, '-') << "\n";
    
    for (const auto& result : results) {
        cout << left << setw(20) << result.algorithm_name
             << setw(15) << result.time_ms
             << setw(15) << result.memory_usage_kb
             << setw(15) << result.matches.size()
             << "\n";
    }
    cout << string(80, '=') << "\n";
}

void Benchmark::runComparison(const string& pattern, const string& fastaPath) {
    cout << "Running algorithm comparison benchmark...\n";
    cout << "Pattern: " << pattern << " (" << pattern.length() << " chars)\n";
    
    vector<BenchmarkResult> results;
    vector<string> algorithms = {"bmh", "kmp", "bithiftor"};
    
    for (const auto& algo : algorithms) {
        if (algo == "bithiftor" && pattern.length() > 64) {
            cout << "Skipping BitParallel (pattern too long: " << pattern.length() << " > 64)\n";
            continue;
        }
        cout << "Running " << algo << "... ";
        results.push_back(runAlgorithm(algo, pattern, fastaPath));
        cout << "Done\n";
    }
    
    printResults(results);
}

void Benchmark::runHybridEvaluation(const string& pattern, const string& fastaPath) {
    cout << "\nEvaluating Hybrid Picker vs Individual Algorithms...\n";
    
    vector<BenchmarkResult> results;
    
    // Test individual algorithms
    vector<string> algorithms = {"bmh", "kmp", "bithiftor"};
    for (const auto& algo : algorithms) {
        if (algo == "bithiftor" && pattern.length() > 64) continue;
        results.push_back(runAlgorithm(algo, pattern, fastaPath));
    }
    
    // Test hybrid auto-picker
    cout << "Running Hybrid Picker (auto-selection)... ";
    auto start = high_resolution_clock::now();
    size_t memory_before = getMemoryUsage();
    
    HybridPicker picker;
    auto hybrid_matches = picker.autoPickAndSearch(pattern, fastaPath);
    
    auto end = high_resolution_clock::now();
    size_t memory_after = getMemoryUsage();
    
    BenchmarkResult hybrid_result;
    hybrid_result.algorithm_name = "HYBRID_AUTO";
    hybrid_result.matches = hybrid_matches;
    hybrid_result.time_ms = duration_cast<milliseconds>(end - start).count();
    hybrid_result.memory_usage_kb = memory_after - memory_before;
    results.push_back(hybrid_result);
    cout << "Done\n";
    
    printResults(results);
    
    // Analysis
    cout << "\nANALYSIS:\n";
    auto fastest = min_element(results.begin(), results.end(), 
        [](const auto& a, const auto& b) { return a.time_ms < b.time_ms; });
    
    cout << "Fastest algorithm: " << fastest->algorithm_name 
         << " (" << fastest->time_ms << " ms)\n";
    
    auto most_efficient = min_element(results.begin(), results.end(),
        [](const auto& a, const auto& b) { return a.memory_usage_kb < b.memory_usage_kb; });
    
    cout << "Most memory efficient: " << most_efficient->algorithm_name
         << " (" << most_efficient->memory_usage_kb << " KB)\n";
}