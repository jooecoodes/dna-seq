#include "../../include/KMP.hpp"
#include "../../include/FastaReader.hpp"


#include <vector>
#include <string>
#include <stdexcept>

using namespace std;

vector<int> KMP::computeLPS(const string& pattern) const {
    int m = pattern.length();
    vector<int> lps(m, 0);
    int len = 0; // Length of the previous longest prefix suffix
    int i = 1;

    while (i < m) {
        if (pattern[i] == pattern[len]) {
            len++;
            lps[i] = len;
            i++;
        } else {
            if (len != 0) {
                len = lps[len - 1];
            } else {
                lps[i] = 0;
                i++;
            }
        }
    }
    return lps;
}

vector<size_t> KMP::search(const string& pattern, const string& text) const {
    vector<size_t> matches;
    int n = text.length();
    int m = pattern.length();
    
    if (m == 0 || n < m) return matches;

    vector<int> lps = computeLPS(pattern);
    int i = 0; // Index for text
    int j = 0; // Index for pattern

    while (i < n) {
        if (pattern[j] == text[i]) {
            i++;
            j++;
        }

        if (j == m) {
            matches.push_back(i - j); // Found match at position i-j
            j = lps[j - 1];
        } else if (i < n && pattern[j] != text[i]) {
            if (j != 0) {
                j = lps[j - 1];
            } else {
                i++;
            }
        }
    }
    
    return matches;
}

vector<size_t> KMP::searchInFasta(const string& pattern, const string& fastaPath) const {
    string dnaSequence = FastaReader::readSequence(fastaPath);
    return search(pattern, dnaSequence);
}