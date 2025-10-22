#include "../../include/HybridPicker.hpp"
#include "../../include/BioUtils.hpp"

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

size_t HybridPicker::pickAndSearch(const string& algorithmName, 
                                         const string& pattern, 
                                         const string& fastaPath) {
    auto matcher = createMatcher(algorithmName);
    if (!matcher) {
        throw invalid_argument("Unknown algorithm: " + algorithmName + 
                              ". Available: bmh, kmp, bithiftor");
    }
    return matcher->searchInFasta(pattern, fastaPath);
}

size_t HybridPicker::autoPickAndSearch(const string& pattern, 
                                             const string& fastaPath) {
    string bestAlgorithm = recommendAlgorithm(pattern);
    cout << "Hybrid Picker selected: " << bestAlgorithm << " algorithm" << endl;
    return pickAndSearch(bestAlgorithm, pattern, fastaPath);
}

size_t HybridPicker::pickAndSearchParallel(const string& algorithmName, 
                                         const string& pattern, 
                                         const string& fastaPath) {
    auto matcher = createMatcher(algorithmName);
    if (!matcher) {
        throw invalid_argument("Unknown algorithm: " + algorithmName + 
                              ". Available: bmh, kmp, bithiftor");
    }
    return matcher->searchParallelInFasta(pattern, fastaPath);
}

size_t HybridPicker::autoPickAndSearchParallel(const string& pattern, 
                                             const string& fastaPath) {
    string bestAlgorithm = recommendAlgorithm(pattern);
    cout << "Hybrid Picker selected: " << bestAlgorithm << " algorithm as parallel" << endl;
    return pickAndSearchParallel(bestAlgorithm, pattern, fastaPath);
}


string HybridPicker::recommendAlgorithm(const string& pattern) {
    size_t length = pattern.length();
    double entropy = BioUtils::calculateShannonEntropy(pattern);
    double gc_content = BioUtils::calculateGCContent(pattern);
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

size_t HybridPicker::searchWithReverseComplementHybrid( const string& pattern, 
                                         const string& text, 
                                         const string& algorithmName, 
                                         bool parallel) {
    auto matcher = createMatcher(algorithmName);
    if (!matcher) {
        throw invalid_argument("Unknown algorithm: " + algorithmName + 
                              ". Available: bmh, kmp, bithiftor");
    }
    return matcher->searchWithReverseComplement(pattern, text, parallel);
}