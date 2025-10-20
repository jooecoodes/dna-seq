#pragma once

#include "HybridPicker.hpp"
#include <chrono>
#include <vector>
#include <string>
#include <map>
#include <iostream>

class Benchmark {
private:
    struct BenchmarkResult {
        std::vector<size_t> matches;
        long long time_ms;
        size_t memory_usage_kb;
        std::string algorithm_name;
    };

public:
    static void runComparison(const std::string& pattern, const std::string& fastaPath);
    static void runHybridEvaluation(const std::string& pattern, const std::string& fastaPath);
    
private:
    static long long getMemoryUsage();
    static BenchmarkResult runAlgorithm(const std::string& algorithmName, 
                                       const std::string& pattern, 
                                       const std::string& fastaPath);
    static void printResults(const std::vector<BenchmarkResult>& results);
};