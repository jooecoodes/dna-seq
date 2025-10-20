#pragma once

#include <string>
#include <vector>

class PatternMatcher {
public:
    virtual ~PatternMatcher() = default;
    virtual std::vector<size_t> search(const std::string& pattern, const std::string& text) const = 0;
    virtual std::vector<size_t> searchInFasta(const std::string& pattern, const std::string& fastaPath) const = 0;
};