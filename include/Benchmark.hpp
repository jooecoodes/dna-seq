#pragma once

#include <iostream>
#include <string>
#include <vector>

class Benchmark {
    private:
        struct BenchmarkResult {
            std::vector<size_t> matches;
            long long time_ms;
            size_t memory_usage_kb;
            std::string algorithm_name;
        };
    public:
        static void run(const std::string& pattern, const std::string& fastaPath);

};