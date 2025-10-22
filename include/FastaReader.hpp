#pragma once
#include <string>

class FastaReader {
public:
    static std::string readSequence(const std::string& fastaPath);
};