#include "../../include/FastaReader.hpp"

using namespace std;

string FastaReader::readSequence(const string& fastaPath) {
    ifstream file(fastaPath);
    if (!file.is_open()) {
        throw runtime_error("Cannot open FASTA file: " + fastaPath);
    }
    
    string line;
    string sequence;
    bool inSequence = false;
    
    while (getline(file, line)) {
        if (line.empty()) continue;
        
        if (line[0] == '>') {
            // Header line - if we were already reading sequence, we're done
            if (inSequence) break;
            inSequence = true;
            continue;
        }
        
        if (inSequence) {
            // Remove any whitespace and add to sequence
            for (char c : line) {
                if (!isspace(c)) {
                    sequence += c;
                }
            }
        }
    }
    
    file.close();
    
    if (sequence.empty()) {
        throw runtime_error("No DNA sequence found in FASTA file: " + fastaPath);
    }
    
    return sequence;
}