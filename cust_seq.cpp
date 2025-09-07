#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <unordered_map>
#include <algorithm>
#include <random>
#include <cstring>
#include <windows.h>
#include <psapi.h>

using namespace std;
using namespace std::chrono;

// Memory tracking variables
size_t total_memory_allocated = 0;
size_t peak_memory_used = 0;
size_t current_memory_used = 0;

// Overload new operator to track memory allocation
void* operator new(size_t size) {
    total_memory_allocated += size;
    current_memory_used += size;
    peak_memory_used = max(peak_memory_used, current_memory_used);
    return malloc(size);
}

// Overload delete operator to track memory deallocation
void operator delete(void* memory, size_t size) {
    current_memory_used -= size;
    free(memory);
}

// Function to get current process memory usage
size_t get_current_memory_usage() {
    PROCESS_MEMORY_COUNTERS_EX pmc;
    GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
    return pmc.PrivateUsage;
}

// Function to reset memory tracking
void reset_memory_tracking() {
    total_memory_allocated = 0;
    peak_memory_used = 0;
    current_memory_used = 0;
}

// Function to generate a random DNA sequence with customizable GC content
string generate_dna_sequence(int length, double gc_content = 0.5) {
    string bases = "ATGC";
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);
    
    string sequence;
    for (int i = 0; i < length; i++) {
        if (dis(gen) < gc_content) {
            sequence += (dis(gen) < 0.5) ? 'G' : 'C';
        } else {
            sequence += (dis(gen) < 0.5) ? 'A' : 'T';
        }
    }
    return sequence;
}

// Function to generate a repetitive DNA sequence
string generate_repetitive_sequence(int length, const string& motif, double mutation_rate = 0.1) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);
    uniform_int_distribution<> base_dis(0, 3);
    string bases = "ATGC";
    
    string sequence;
    int motif_length = motif.length();
    
    for (int i = 0; i < length; i++) {
        int motif_index = i % motif_length;
        char base = motif[motif_index];
        
        // Apply mutation with given probability
        if (dis(gen) < mutation_rate) {
            base = bases[base_dis(gen)];
        }
        
        sequence += base;
    }
    
    return sequence;
}

// Function to create a pattern with specific characteristics
string create_custom_pattern(int length, double gc_content, double repetitiveness, 
                            string& motif_used, double& actual_gc_content) {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> motif_len_dis(2, 10);
    uniform_int_distribution<> base_dis(0, 3);
    string bases = "ATGC";
    
    // Generate a random motif
    int motif_length = motif_len_dis(gen);
    string motif;
    for (int i = 0; i < motif_length; i++) {
        motif += bases[base_dis(gen)];
    }
    motif_used = motif;
    
    // Create pattern based on repetitiveness
    string pattern;
    if (repetitiveness > 0.7) {
        // Highly repetitive pattern
        pattern = generate_repetitive_sequence(length, motif, 0.1);
    } else if (repetitiveness > 0.3) {
        // Moderately repetitive pattern
        pattern = generate_repetitive_sequence(length, motif, 0.3);
    } else {
        // Non-repetitive pattern
        pattern = generate_dna_sequence(length, gc_content);
    }
    
    // Adjust GC content if needed
    if (gc_content > 0) {
        int current_gc = count(pattern.begin(), pattern.end(), 'G') + count(pattern.begin(), pattern.end(), 'C');
        double current_gc_content = static_cast<double>(current_gc) / length;
        int gc_to_change = static_cast<int>((gc_content - current_gc_content) * length);
        
        if (gc_to_change > 0) {
            // Need to add more GC bases
            for (int i = 0; i < gc_to_change && i < length; i++) {
                if (pattern[i] == 'A' || pattern[i] == 'T') {
                    pattern[i] = (rand() % 2 == 0) ? 'G' : 'C';
                }
            }
        } else if (gc_to_change < 0) {
            // Need to add more AT bases
            for (int i = 0; i < -gc_to_change && i < length; i++) {
                if (pattern[i] == 'G' || pattern[i] == 'C') {
                    pattern[i] = (rand() % 2 == 0) ? 'A' : 'T';
                }
            }
        }
    }
    
    // Calculate actual GC content
    int gc_count = count(pattern.begin(), pattern.end(), 'G') + count(pattern.begin(), pattern.end(), 'C');
    actual_gc_content = static_cast<double>(gc_count) / length;
    
    return pattern;
}

// KMP Algorithm
vector<int> kmp_search(const string& text, const string& pattern) {
    int n = text.length();
    int m = pattern.length();
    vector<int> matches;
    if (m == 0) return matches;

    // Precompute LPS array
    vector<int> lps(m, 0);
    int len = 0;
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

    // Search using LPS
    i = 0;
    int j = 0;
    while (i < n) {
        if (pattern[j] == text[i]) {
            j++;
            i++;
        }
        
        if (j == m) {
            matches.push_back(i - j);
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

// Boyer-Moore Algorithm (Bad Character Heuristic)
vector<int> boyer_moore_search(const string& text, const string& pattern) {
    int n = text.length();
    int m = pattern.length();
    vector<int> matches;
    if (m == 0) return matches;

    // Bad character heuristic
    unordered_map<char, int> badChar;
    for (int i = 0; i < m; i++) {
        badChar[pattern[i]] = i;
    }

    int s = 0;
    while (s <= n - m) {
        int j = m - 1;
        
        while (j >= 0 && pattern[j] == text[s + j]) {
            j--;
        }
        
        if (j < 0) {
            matches.push_back(s);
            s += (s + m < n) ? m - badChar[text[s + m]] : 1;
        } else {
            s += max(1, j - badChar[text[s + j]]);
        }
    }
    
    return matches;
}

// Bit-Parallel Algorithm (Shift-Or)
vector<int> bit_parallel_search(const string& text, const string& pattern) {
    int m = pattern.length();
    int n = text.length();
    vector<int> matches;
    if (m == 0 || m > 64) return matches;

    unsigned long long pattern_mask[256];
    for (int i = 0; i < 256; i++) {
        pattern_mask[i] = ~0;
    }
    
    for (int i = 0; i < m; i++) {
        pattern_mask[pattern[i]] &= ~(1ULL << i);
    }

    unsigned long long R = ~1;
    for (int i = 0; i < n; i++) {
        R |= pattern_mask[text[i]];
        R <<= 1;
        
        if ((R & (1ULL << m)) == 0) {
            matches.push_back(i - m + 1);
        }
    }
    
    return matches;
}

// Simple pattern analysis for hybrid selection
bool is_repetitive(const string& pattern) {
    if (pattern.length() <= 10) return false;
    
    int max_repeat = 1;
    int current_repeat = 1;
    for (int i = 1; i < pattern.length(); i++) {
        if (pattern[i] == pattern[i-1]) {
            current_repeat++;
            max_repeat = max(max_repeat, current_repeat);
        } else {
            current_repeat = 1;
        }
    }
    
    return max_repeat > pattern.length() / 3;
}

// Hybrid algorithm selector that returns which algorithm was chosen
pair<vector<int>, string> hybrid_search(const string& text, const string& pattern) {
    int length = pattern.length();
    if (length == 0) return make_pair(vector<int>(), "None");
    
    // Calculate GC content
    int gc_count = 0;
    for (char c : pattern) {
        if (c == 'G' || c == 'C') gc_count++;
    }
    float gc_content = static_cast<float>(gc_count) / length;
    
    // Check repetitiveness
    bool repetitive = is_repetitive(pattern);
    
    // Algorithm selection
    if (length <= 64 && !repetitive) {
        return make_pair(bit_parallel_search(text, pattern), "Bit-Parallel");
    } else if (repetitive || gc_content > 0.6) {
        return make_pair(kmp_search(text, pattern), "KMP");
    } else {
        return make_pair(boyer_moore_search(text, pattern), "Boyer-Moore");
    }
}

// Function to measure execution time and memory usage of an algorithm
void measure_algorithm(const string& name, vector<int> (*algorithm)(const string&, const string&), 
                      const string& text, const string& pattern) {
    // Reset memory tracking
    reset_memory_tracking();
    size_t initial_memory = get_current_memory_usage();
    
    // Measure time and execute algorithm
    auto start = high_resolution_clock::now();
    vector<int> matches = algorithm(text, pattern);
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    
    // Get memory usage
    size_t final_memory = get_current_memory_usage();
    size_t memory_used = final_memory - initial_memory;
    
    // Output results
    cout << "Algorithm: " << name << endl;
    cout << "Matches found: " << matches.size() << endl;
    if (!matches.empty()) {
        cout << "First match at position: " << matches[0] << endl;
    }
    cout << "Execution time: " << duration.count() << " microseconds" << endl;
    cout << "Memory used: " << memory_used << " bytes" << endl;
    cout << "Peak memory during execution: " << peak_memory_used << " bytes" << endl;
    cout << "Total memory allocated: " << total_memory_allocated << " bytes" << endl;
    cout << "----------------------------------------" << endl;
}

// Special function to measure the hybrid algorithm and show which algorithm it chose
void measure_hybrid_algorithm(const string& text, const string& pattern) {
    // Reset memory tracking
    reset_memory_tracking();
    size_t initial_memory = get_current_memory_usage();
    
    // Measure time and execute algorithm
    auto start = high_resolution_clock::now();
    auto result = hybrid_search(text, pattern);
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    
    vector<int> matches = result.first;
    string chosen_algorithm = result.second;
    
    // Get memory usage
    size_t final_memory = get_current_memory_usage();
    size_t memory_used = final_memory - initial_memory;
    
    // Output results
    cout << "Algorithm: Hybrid" << endl;
    cout << "Chosen algorithm: " << chosen_algorithm << endl;
    cout << "Selection criteria:" << endl;
    
    // Calculate GC content for display
    int gc_count = 0;
    for (char c : pattern) {
        if (c == 'G' || c == 'C') gc_count++;
    }
    float gc_content = static_cast<float>(gc_count) / pattern.length();
    bool repetitive = is_repetitive(pattern);
    
    cout << "  - Pattern length: " << pattern.length() << endl;
    cout << "  - GC content: " << gc_content * 100 << "%" << endl;
    cout << "  - Repetitive: " << (repetitive ? "Yes" : "No") << endl;
    
    cout << "Matches found: " << matches.size() << endl;
    if (!matches.empty()) {
        cout << "First match at position: " << matches[0] << endl;
    }
    cout << "Execution time: " << duration.count() << " microseconds" << endl;
    cout << "Memory used: " << memory_used << " bytes" << endl;
    cout << "Peak memory during execution: " << peak_memory_used << " bytes" << endl;
    cout << "Total memory allocated: " << total_memory_allocated << " bytes" << endl;
    cout << "----------------------------------------" << endl;
}

int main() {
    // Get user input for pattern characteristics
    int pattern_length;
    double gc_content, repetitiveness;
    
    cout << "DNA Sequence Matching with Custom Pattern Generation" << endl;
    cout << "====================================================" << endl;
    cout << "Enter pattern length: ";
    cin >> pattern_length;
    cout << "Enter desired GC content for pattern (0.0 to 1.0): ";
    cin >> gc_content;
    cout << "Enter desired repetitiveness for pattern (0.0 to 1.0): ";
    cin >> repetitiveness;
    
    // Validate inputs
    pattern_length = max(1, min(1000, pattern_length));
    gc_content = max(0.0, min(1.0, gc_content));
    repetitiveness = max(0.0, min(1.0, repetitiveness));
    
    // Generate a 5000-character DNA sequence
    string text = generate_dna_sequence(5000, 0.5);
    
    // Generate custom pattern
    string motif_used;
    double actual_gc_content;
    string pattern = create_custom_pattern(pattern_length, gc_content, repetitiveness, 
                                         motif_used, actual_gc_content);
    
    // Ensure the pattern appears at least once in the text
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, text.length() - pattern.length());
    int insert_position = dis(gen);
    text.replace(insert_position, pattern.length(), pattern);
    
    cout << "\nDNA Sequence Matching Analysis" << endl;
    cout << "===============================" << endl;
    cout << "Text length: " << text.length() << " characters" << endl;
    cout << "Pattern length: " << pattern.length() << " characters" << endl;
    cout << "Pattern inserted at position: " << insert_position << endl;
    cout << "Requested GC content: " << gc_content * 100 << "%" << endl;
    cout << "Actual GC content: " << actual_gc_content * 100 << "%" << endl;
    cout << "Requested repetitiveness: " << repetitiveness * 100 << "%" << endl;
    
    if (repetitiveness > 0.3) {
        cout << "Motif used: " << motif_used << endl;
    }
    
    // Calculate pattern characteristics
    bool repetitive = is_repetitive(pattern);
    cout << "Pattern is repetitive: " << (repetitive ? "Yes" : "No") << endl;
    cout << "===============================" << endl;
    cout << endl;
    
    // Test each algorithm
    measure_algorithm("KMP", kmp_search, text, pattern);
    measure_algorithm("Boyer-Moore", boyer_moore_search, text, pattern);
    
    // Only run Bit-Parallel if pattern is 64 characters or less
    if (pattern.length() <= 64) {
        measure_algorithm("Bit-Parallel", bit_parallel_search, text, pattern);
    } else {
        cout << "Bit-Parallel: Skipped (pattern > 64 characters)" << endl;
        cout << "----------------------------------------" << endl;
    }
    
    // Test the hybrid algorithm
    measure_hybrid_algorithm(text, pattern);
    
    return 0;
}