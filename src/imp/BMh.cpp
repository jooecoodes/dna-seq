#include "../../include/BM.hpp"
#include "../../include/FastaReader.hpp"

#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>

using namespace std;

vector<int> BoyerMooreHorspool::createBadCharTable(const string& pattern) const {
    const int ALPHABET_SIZE = 256;
    vector<int> badChar(ALPHABET_SIZE, pattern.length());
    
    int m = pattern.length();
    for (int i = 0; i < m - 1; i++) {
        badChar[static_cast<unsigned char>(pattern[i])] = m - 1 - i;
    }
    
    return badChar;
}

vector<size_t> BoyerMooreHorspool::search(const string& pattern, const string& text) const {
    vector<size_t> matches;
    int n = text.length();
    int m = pattern.length();
    
    if (m == 0 || n < m) return matches;

    vector<int> badChar = createBadCharTable(pattern);
    
    size_t i = 0;
    while (i <= n - m) {
        int j = m - 1;
        
        // Check from right to left
        while (j >= 0 && pattern[j] == text[i + j]) {
            j--;
        }
        
        if (j < 0) {
            // Pattern found
            matches.push_back(i);
            i++; // Move forward by 1 after match
        } else {
            // Shift based on bad character rule
            unsigned char mismatchChar = text[i + j];
            i += max(1, badChar[mismatchChar] - (m - 1 - j));
        }
    }
    
    return matches;
}

vector<size_t> BoyerMooreHorspool::searchInFasta(const string& pattern, const string& fastaPath) const {
    string dnaSequence = FastaReader::readSequence(fastaPath);
    return search(pattern, dnaSequence);
}