#include "../../include/HybridPicker.hpp"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

unique_ptr<PatternMatcher> HybridPicker::createMatcher(const string& algorithmName) {
    if (algorithmName == "bmh") return make_unique<BoyerMooreHorspool>();
    if (algorithmName == "kmp") return make_unique<KMP>();
    if (algorithmName == "bithiftor") return make_unique<BitParallelShiftOr>();
    return nullptr;
}

double HybridPicker::calculateShannonEntropy(const std::string& pattern) {
    if (pattern.empty()) return 0.0;
    
    // Count each nucleotide (case insensitive)
    int counts[256] = {0};
    int total = 0;
    
    for (char c : pattern) {
        char upper = std::toupper(c);
        if (upper == 'A' || upper == 'T' || upper == 'G' || upper == 'C') {
            counts[static_cast<unsigned char>(upper)]++;
            total++;
        }
    }
    
    if (total == 0) return 0.0;
    
    // Calculate entropy
    double entropy = 0.0;
    for (int i = 0; i < 256; i++) {
        if (counts[i] > 0) {
            double p = static_cast<double>(counts[i]) / total;
            entropy -= p * std::log2(p);
        }
    }
    
    // Returns 0.0 (repetitive) to 2.0 (complex)
    return entropy;
}

double HybridPicker::calculateGCContent(const std::string& pattern) {
    if (pattern.empty()) return 0.0;
    
    size_t gc_count = 0;
    size_t total_bases = 0;
    
    for (char c : pattern) {
        char upper = std::toupper(c);
        if (upper == 'G' || upper == 'C') {
            gc_count++;
        }
        if (upper == 'A' || upper == 'T' || upper == 'G' || upper == 'C') {
            total_bases++;
        }
    }
    
    if (total_bases == 0) return 0.0;
    return static_cast<double>(gc_count) / total_bases;
}

vector<size_t> HybridPicker::pickAndSearch(const string& algorithmName, 
                                         const string& pattern, 
                                         const string& fastaPath) {
    auto matcher = createMatcher(algorithmName);
    if (!matcher) {
        throw invalid_argument("Unknown algorithm: " + algorithmName + 
                              ". Available: bmh, kmp, bithiftor");
    }
    return matcher->searchInFasta(pattern, fastaPath);
}

vector<size_t> HybridPicker::autoPickAndSearch(const string& pattern, 
                                             const string& fastaPath) {
    string bestAlgorithm = recommendAlgorithm(pattern);
    cout << "Hybrid Picker selected: " << bestAlgorithm << " algorithm" << endl;
    return pickAndSearch(bestAlgorithm, pattern, fastaPath);
}

string HybridPicker::recommendAlgorithm(const string& pattern) {
    size_t length = pattern.length();
    double entropy = calculateShannonEntropy(pattern);
    double gc_content = calculateGCContent(pattern);
    double repetitiveness = 2.0 - entropy;
    
    // Decision Tree
    if (length <= 64) {
        return "bithiftor";
    }
    if (repetitiveness >= 1.5) {
        return "kmp";
    }
    if (length > 1000 && entropy >= 1.2) {
        return "bmh";
    }
    if (gc_content >= 0.65 && length <= 500) {
        return "bmh";
    }
    if (length > 2000) {
        return "kmp";
    }
    if (repetitiveness >= 1.0) {
        return "kmp";
    }
    return "bmh";
}

vector<string> HybridPicker::getAvailableAlgorithms() const {
    return {"bmh", "kmp", "bithiftor"};
}