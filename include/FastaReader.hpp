#pragma once
#include <string>
#include <fstream>
#include <stdexcept>
#include <sstream>

class FastaReader {
public:
    static std::string readSequence(const std::string& fastaPath);
};