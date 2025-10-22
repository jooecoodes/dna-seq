#pragma once

#include "PatternMatcher.hpp"

#include <string>
#include <vector>

class KMP : public PatternMatcher {
private:
    std::vector<size_t> computeLPS(const std::string& pattern) const;
    
public:
    size_t search(const std::string& pattern, const std::string& text) const override;
    size_t searchInFasta(const std::string& pattern, const std::string& fastaPath) const override;
    size_t searchParallel(const std::string& pattern, const std::string& text, int num_threads) const override;
    size_t searchParallelInFasta(const std::string& pattern, const std::string& fastaPath) const override;
    size_t searchWithReverseComplement(const std::string& pattern, const std::string& text, bool parallel) const override;
};