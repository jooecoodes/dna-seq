#pragma once

#include "PatternMatcher.hpp"

#include <string>
#include <vector>

class KMP : public PatternMatcher {
private:
    std::vector<int> computeLPS(const std::string& pattern) const;
    
public:
    std::vector<size_t> search(const std::string& pattern, const std::string& text) const override;
    std::vector<size_t> searchInFasta(const std::string& pattern, const std::string& fastaPath) const override;
};