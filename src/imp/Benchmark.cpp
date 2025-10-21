#include "../../include/Benchmark.hpp"
#include "../../include/HybridPicker.hpp"

#include <string>
#include <iostream>
#include <chrono>



void Benchmark::run(const std::string& pattern, const std::string& fastaPath) {

    HybridPicker picker; // Just use one picker
    auto benchmarkAlgorithm = [&picker](const std::string& algName,
                                        const std::string& pattern,
                                        const std::string& fastaPath,
                                        bool parallel = false) {
        auto start = std::chrono::steady_clock::now();

        size_t matches = parallel
            ? picker.pickAndSearchParallel(algName, pattern, fastaPath)
            : picker.pickAndSearch(algName, pattern, fastaPath);

        auto end = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        std::cout << "Algorithm: " << algName
                << (parallel ? " (Parallel)" : " (Serial)")
                << ", Matches: " << matches
                << ", Time: " << duration << " µs\n";
    };

    std::cout << "Sequential: " << std::endl;
    benchmarkAlgorithm("bmh", pattern, fastaPath, false);
    benchmarkAlgorithm("kmp", pattern, fastaPath, false);
    benchmarkAlgorithm("bithiftor", pattern, fastaPath, false);

    std::cout << "Parallel: " << std::endl;
    benchmarkAlgorithm("bmh", pattern, fastaPath, true);
    benchmarkAlgorithm("kmp", pattern, fastaPath, true);
    benchmarkAlgorithm("bithiftor", pattern, fastaPath, true);

    // HybridPicker picker;

    // size_t bmhMatches = picker.pickAndSearch("bmh", pattern, fastaPath);
    // size_t kmpMatches = picker.pickAndSearch("kmp", pattern, fastaPath);
    // size_t bpMatches = picker.pickAndSearch("bithiftor", pattern, fastaPath);

    // std::cout << "Benchmark Results:\n";
    // std::cout << "Boyer-Moore-Horspool Matches: " << bmhMatches << "\n";
    // std::cout << "Knuth-Morris-Pratt Matches: " << kmpMatches << "\n";
    // std::cout << "Bit-Parallel Shift-Or Matches: " << bpMatches << "\n";

}
