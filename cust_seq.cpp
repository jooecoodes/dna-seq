#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <array>
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
    sequence.reserve(length);  // Pre-allocate memory
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
    sequence.reserve(length);  // Pre-allocate memory
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
    motif.reserve(motif_length);
    for (int i = 0; i < motif_length; i++) {
        motif += bases[base_dis(gen)];
    }
    motif_used = motif;
    
    // Create pattern based on repetitiveness
    string pattern;
    pattern.reserve(length);
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

// KMP Algorithm with optimized memory
vector<int> kmp_search(const string& text, const string& pattern) {
    int n = text.length();
    int m = pattern.length();
    vector<int> matches;
    if (m == 0) return matches;

    // Precompute LPS array with minimal memory
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

// Optimized Boyer-Moore Algorithm with array-based bad character heuristic
vector<int> boyer_moore_search(const string& text, const string& pattern) {
    int n = text.length();
    int m = pattern.length();
    vector<int> matches;
    if (m == 0) return matches;

    // Use array instead of unordered_map for bad character heuristic
    // DNA has only 4 characters, but we'll create a 256 array for all possible chars
    array<int, 256> badChar;
    badChar.fill(-1);
    
    // Fill the actual values
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

// Bit-Parallel Algorithm (Shift-Or) with optimizations
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

// Structure for pattern analysis
struct PatternAnalysis {
    bool is_repetitive;
    bool is_periodic;
    float gc_content;
    int longest_run;
    int distinct_chars;
    int period;
};

// Function to compute the period of a pattern without full LPS (simplified)
int compute_period_simple(const string& pattern) {
    int m = pattern.length();
    if (m == 0) return 0;
    
    // Check for periodicity using a simple method
    for (int p = 1; p <= m / 2; p++) {
        if (m % p != 0) continue;
        bool periodic = true;
        for (int i = p; i < m; i++) {
            if (pattern[i] != pattern[i % p]) {
                periodic = false;
                break;
            }
        }
        if (periodic) return p;
    }
    return m;
}

// Enhanced pattern analysis with memory efficiency
PatternAnalysis analyze_pattern(const string& pattern) {
    PatternAnalysis result;
    int length = pattern.length();
    
    // Calculate GC content
    int gc_count = 0;
    for (char c : pattern) {
        if (c == 'G' || c == 'C') gc_count++;
    }
    result.gc_content = static_cast<float>(gc_count) / length;
    
    // Calculate distinct characters and longest run
    std::array<bool, 256> seen = {false};
    result.distinct_chars = 0;
    result.longest_run = 1;
    int current_run = 1;
    
    for (int i = 0; i < length; i++) {
        if (!seen[pattern[i]]) {
            seen[pattern[i]] = true;
            result.distinct_chars++;
        }
        if (i > 0) {
            if (pattern[i] == pattern[i-1]) {
                current_run++;
                result.longest_run = std::max(result.longest_run, current_run);
            } else {
                current_run = 1;
            }
        }
    }
    
    // Check repetitiveness based on longest run first
    result.is_repetitive = (result.longest_run > length / 3);
    
    // Only compute periodicity if necessary and for shorter patterns
    result.is_periodic = false;
    result.period = length;
    if (result.is_repetitive && length <= 100) {
        // Use simple periodicity check for shorter patterns
        result.period = compute_period_simple(pattern);
        result.is_periodic = (result.period <= length / 2);
    } else if (result.is_repetitive && length > 100) {
        // For long patterns, assume periodic if longest run is significant
        result.is_periodic = (result.longest_run > length / 2);
    }
    
    // Update repetitiveness based on periodicity
    result.is_repetitive = result.is_repetitive || result.is_periodic;
    
    return result;
}

// Enhanced hybrid algorithm selector
pair<vector<int>, string> hybrid_search(const string& text, const string& pattern) {
    int length = pattern.length();
    if (length == 0) return make_pair(vector<int>(), "None");
    
    // Analyze pattern characteristics
    PatternAnalysis analysis = analyze_pattern(pattern);
    
    // Algorithm selection based on realistic thresholds
    if (length <= 64 && !analysis.is_repetitive && analysis.distinct_chars > 2) {
        return make_pair(bit_parallel_search(text, pattern), "Bit-Parallel");
    } else if (analysis.is_repetitive) {
        return make_pair(kmp_search(text, pattern), "KMP");
    } else if (length > 100 && analysis.distinct_chars >= 3) {
        // Boyer-Moore excels for long patterns with good character distribution
        return make_pair(boyer_moore_search(text, pattern), "Boyer-Moore");
    } else {
        // Default to Boyer-Moore for most cases
        return make_pair(boyer_moore_search(text, pattern), "Boyer-Moore");
    }
}

// Function to measure execution time and memory usage of an algorithm
void measure_algorithm(const string& name, vector<int> (*algorithm)(const string&, const string&), 
                      const string& text, const string& pattern) {
    // Reset memory tracking
    reset_memory_tracking();
    size_t initial_memory = get_current_memory_usage();
    
    // Warm-up run to account for CPU caching
    algorithm(text, pattern);
    
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
    
    // Analyze pattern for display
    PatternAnalysis analysis = analyze_pattern(pattern);
    
    // Output results
    cout << "Algorithm: Hybrid" << endl;
    cout << "Chosen algorithm: " << chosen_algorithm << endl;
    cout << "Selection criteria:" << endl;
    cout << "  - Pattern length: " << pattern.length() << endl;
    cout << "  - GC content: " << analysis.gc_content * 100 << "%" << endl;
    cout << "  - Repetitive: " << (analysis.is_repetitive ? "Yes" : "No") << endl;
    cout << "  - Longest run: " << analysis.longest_run << " characters" << endl;
    cout << "  - Distinct characters: " << analysis.distinct_chars << endl;
    
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
    
    // Analyze pattern for display
    PatternAnalysis analysis = analyze_pattern(pattern);
    
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
    
    cout << "Pattern is repetitive: " << (analysis.is_repetitive ? "Yes" : "No") << endl;
    cout << "Longest character run: " << analysis.longest_run << endl;
    cout << "Distinct characters: " << analysis.distinct_chars << endl;
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