#pragma once

#include <string>
#include <algorithm>
#include <cmath>

namespace BioUtils {
    
    inline std::string reverseComplement(const std::string& dna) {
        std::string rc;
        rc.reserve(dna.length());
        
        for (int i = dna.length() - 1; i >= 0; i--) {
            switch (dna[i]) {
                case 'A': rc += 'T'; break;
                case 'T': rc += 'A'; break;
                case 'G': rc += 'C'; break;
                case 'C': rc += 'G'; break;
                case 'a': rc += 't'; break;
                case 't': rc += 'a'; break;
                case 'g': rc += 'c'; break;
                case 'c': rc += 'g'; break;
                case 'N': rc += 'N'; break;
                case 'n': rc += 'n'; break;
                default:  rc += dna[i]; break;  // Handle other IUPAC codes
            }
        }
        return rc;
    }
    
    
    inline double calculateShannonEntropy(const std::string& pattern) {
        if (pattern.empty()) return 0.0;
        
        // Count each nucleotide (case insensitive)
        int counts[256] = {0};
        int total = 0;
        
        for (char c : pattern) {
            char upper = std::toupper(c);
            if (upper == 'A' || upper == 'T' || upper == 'G' || upper == 'C') {
                counts[static_cast<unsigned char>(upper)]++;
                total++;
            }
        }
        
        if (total == 0) return 0.0;
        
        // Calculate entropy
        double entropy = 0.0;
        for (int i = 0; i < 256; i++) {
            if (counts[i] > 0) {
                double p = static_cast<double>(counts[i]) / total;
                entropy -= p * std::log2(p);
            }
        }
        
        // Returns 0.0 (repetitive) to 2.0 (complex)
        return entropy;
    }

    inline double calculateGCContent(const std::string& pattern) {
        if (pattern.empty()) return 0.0;
        
        size_t gc_count = 0;
        size_t total_bases = 0;
        
        for (char c : pattern) {
            char upper = std::toupper(c);
            if (upper == 'G' || upper == 'C') {
                gc_count++;
            }
            if (upper == 'A' || upper == 'T' || upper == 'G' || upper == 'C') {
                total_bases++;
            }
        }
        
        if (total_bases == 0) return 0.0;
        return static_cast<double>(gc_count) / total_bases;
    }

    inline std::string toUpperCaseDNA(const std::string& dna) {
        std::string result = dna;
        std::transform(result.begin(), result.end(), result.begin(), ::toupper);
        return result;
    }
}