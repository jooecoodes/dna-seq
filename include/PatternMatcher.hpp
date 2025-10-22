#pragma once

#include <string>
#include <vector>

class PatternMatcher {
public:
    virtual ~PatternMatcher() = default;
    virtual size_t search(const std::string& pattern, const std::string& text) const = 0;
    virtual size_t searchInFasta(const std::string& pattern, const std::string& fastaPath) const = 0;
    virtual size_t searchParallel(const std::string& pattern, const std::string& text, int num_threads) const = 0;
    virtual size_t searchParallelInFasta(const std::string& pattern, const std::string& fastaPath) const = 0;
    virtual size_t searchWithReverseComplement(const std::string& pattern, const std::string& text, bool parallel) const = 0;   
};