#include "../../include/BP.hpp"
#include "../../include/FastaReader.hpp"

#include <vector>
#include <cstdint>
#include <string>

using namespace std;

vector<size_t> BitParallelShiftOr::search(const string& pattern, const string& text) const {
    vector<size_t> matches;
    int n = text.length();
    int m = pattern.length();
    
    // Handle edge cases
    if (m == 0 || n < m) return matches;
    if (m > 64) {
        // For patterns longer than 64, you might want to throw an exception
        // or return empty results since Bit-parallel has this limitation
        return matches;
    }

    // Preprocessing: Build the B table
    uint64_t B[256];
    for (int i = 0; i < 256; ++i) {
        B[i] = ~0ULL; // Initialize all bits to 1
    }
    
    for (int i = 0; i < m; ++i) {
        unsigned char c = pattern[i];
        B[c] &= ~(1ULL << i); // Set the i-th bit to 0 for character c
    }

    // Search phase
    uint64_t state = ~0ULL; // Initialize all bits to 1
    
    for (int i = 0; i < n; ++i) {
        unsigned char text_char = text[i];
        state = (state << 1) | B[text_char];
        
        // Check if we found a match
        if (i >= m - 1 && (state & (1ULL << (m - 1))) == 0) {
            matches.push_back(i - m + 1); // Store starting position
        }
    }
    
    return matches;
}

vector<size_t> BitParallelShiftOr::searchInFasta(const string& pattern, const string& fastaPath) const {
    string dnaSequence = FastaReader::readSequence(fastaPath);
    return search(pattern, dnaSequence);
}

