#include "../include/KMP.hpp"
#include "../include/BM.hpp"
#include "../include/BP.hpp"
#include "../include/HybridPicker.hpp"
#include "../include/Benchmark.hpp"
#include <iostream>

int main() {
    // std::string text = "ATGCTAGCTAGCTAGCTAGC";
    // std::string pattern = "TAGCT";
    
    // // Test KMP
    // KMP kmp;
    // auto kmpResults = kmp.search(pattern, text);
    // std::cout << "KMP found " << kmpResults.size() << " matches" << std::endl;
    
    // // Test Boyer-Moore-Horspool
    // BoyerMooreHorspool bmh;
    // auto bmhResults = bmh.search(pattern, text);
    // std::cout << "Boyer-Moore-Horspool found " << bmhResults.size() << " matches" << std::endl;
    
    // // Test Bit-Parallel (only for short patterns)
    // if (pattern.length() <= 64) {
    //     BitParallelShiftOr bitOr;
    //     auto bitResults = bitOr.search(pattern, text);
    //     std::cout << "Bit-Parallel found " << bitResults.size() << " matches" << std::endl;
    // }

    // std::string pattern, fastaPath;
    
    // std::cout << "DNA Pattern Matcher Benchmark\n";
    // std::cout << "Enter pattern to search: ";
    // std::cin >> pattern;
    
    // std::cout << "Enter FASTA file path: ";
    // std::cin >> fastaPath;
    
    // // Run individual algorithm comparison
    // Benchmark::runComparison(pattern, fastaPath);
    
    // // Run hybrid vs individual evaluation
    // Benchmark::runHybridEvaluation(pattern, fastaPath);
    
    std::string pattern, fastaPath;

    std::cout << "DNA Pattern Matcher Benchmark\n";
    std::cout << "Enter pattern to search: ";
    std::cin >> pattern;

    std::cout << "Enter FASTA file path: ";
    std::cin >> fastaPath;
    
    Benchmark::run(pattern, fastaPath);
    
    return 0;
}