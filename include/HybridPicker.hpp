#pragma once

#include "PatternMatcher.hpp"
#include "BM.hpp"
#include "KMP.hpp"
#include "BP.hpp"
#include <memory>
#include <vector>
#include <string>
#include <stdexcept>

class HybridPicker {
private:
    std::unique_ptr<PatternMatcher> createMatcher(const std::string& algorithmName);
    double calculateShannonEntropy(const std::string& pattern);
    double calculateGCContent(const std::string& pattern);
    
public:
    /**
     * @brief Selects and executes the appropriate pattern matching algorithm
     * @param algorithmName Name of the algorithm: "bmh", "kmp", or "bithiftor"
     * @param pattern The DNA pattern to search for
     * @param fastaPath Path to the FASTA file containing the DNA sequence
     * @return Vector of positions where the pattern was found
     * @throws std::invalid_argument if algorithmName is unknown
     */
    std::vector<size_t> pickAndSearch(const std::string& algorithmName, 
                                     const std::string& pattern, 
                                     const std::string& fastaPath);
    
    /**
     * @brief Automatically selects the best algorithm based on pattern characteristics
     * @param pattern The DNA pattern to search for
     * @param fastaPath Path to the FASTA file containing the DNA sequence
     * @return Vector of positions where the pattern was found
     */
    std::vector<size_t> autoPickAndSearch(const std::string& pattern, 
                                         const std::string& fastaPath);
    
                                         /**
     * @brief Recommends appropriate algorithm based on conditions
     * @return String of the name of the algorithm
     */
    std::string recommendAlgorithm(const std::string& pattern);
    
    /**
     * @brief Gets list of available algorithms
     * @return Vector of algorithm names
     */
    std::vector<std::string> getAvailableAlgorithms() const;
};