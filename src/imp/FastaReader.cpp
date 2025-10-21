#include "../../include/FastaReader.hpp"

#include <iostream>

using namespace std;

string FastaReader::readSequence(const string& fastaPath) {
      std::ifstream f(fastaPath);
    if (!f) {
        std::cerr << "Error: cannot open " << fastaPath << "\n";
        std::exit(1);
    }
    std::string buf, seq;
    seq.reserve(10'000'000);
    while (std::getline(f, buf)) {
        if (!buf.empty() && buf[0] == '>') continue;
        for (char c : buf) {
            if (c != '\r' && c != '\n') {
                c = std::toupper(static_cast<unsigned char>(c));
                if (c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N')
                    seq.push_back(c);
            }
        }
    }
    return seq;
}