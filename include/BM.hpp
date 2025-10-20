#pragma once

#include "PatternMatcher.hpp"

#include <string>
#include <vector>

class BoyerMooreHorspool : public PatternMatcher {
private:
    std::vector<int> createBadCharTable(const std::string& pattern) const;
    
public:
    std::vector<size_t> search(const std::string& pattern, const std::string& text) const override;
    std::vector<size_t> searchInFasta(const std::string& pattern, const std::string& fastaPath) const override;
};