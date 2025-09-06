import ctypes
import time
import random
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import platform
from tabulate import tabulate
import seaborn as sns

# Load the shared library
if platform.system() == "Windows":
    lib_name = "dna_matchers.dll"
else:
    lib_name = "dna_matchers.so"

lib_path = os.path.join(os.path.dirname(__file__), lib_name)
try:
    lib = ctypes.CDLL(lib_path)
    print(f"Successfully loaded library: {lib_path}")
except Exception as e:
    print(f"Error loading library: {e}")
    exit(1)

# Define function prototypes
lib.kmp_search_c.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.POINTER(ctypes.c_int)]
lib.kmp_search_c.restype = ctypes.POINTER(ctypes.c_int)

lib.boyer_moore_search_c.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.POINTER(ctypes.c_int)]
lib.boyer_moore_search_c.restype = ctypes.POINTER(ctypes.c_int)

lib.bit_parallel_search_c.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.POINTER(ctypes.c_int)]
lib.bit_parallel_search_c.restype = ctypes.POINTER(ctypes.c_int)

lib.free_result.argtypes = [ctypes.POINTER(ctypes.c_int)]

# Python wrappers for each algorithm
def kmp_search(text, pattern):
    result_size = ctypes.c_int()
    result_ptr = lib.kmp_search_c(text.encode(), pattern.encode(), ctypes.byref(result_size))
    results = [result_ptr[i] for i in range(result_size.value)]
    lib.free_result(result_ptr)
    return results

def boyer_moore_search(text, pattern):
    result_size = ctypes.c_int()
    result_ptr = lib.boyer_moore_search_c(text.encode(), pattern.encode(), ctypes.byref(result_size))
    results = [result_ptr[i] for i in range(result_size.value)]
    lib.free_result(result_ptr)
    return results

def bit_parallel_search(text, pattern):
    result_size = ctypes.c_int()
    result_ptr = lib.bit_parallel_search_c(text.encode(), pattern.encode(), ctypes.byref(result_size))
    results = [result_ptr[i] for i in range(result_size.value)]
    lib.free_result(result_ptr)
    return results

# Enhanced hybrid selector function
def hybrid_search(text, pattern):
    length = len(pattern)
    if length == 0:
        return []
    
    # Calculate GC content
    gc_count = pattern.count('G') + pattern.count('C')
    gc_content = gc_count / length
    
    # Enhanced repetitiveness detection
    repetitive = is_highly_repetitive(pattern)
    
    # Enhanced algorithm selection logic
    if length <= 64 and not repetitive:
        # Bit-parallel for short, non-repetitive patterns
        return bit_parallel_search(text, pattern)
    elif length <= 200 and gc_content < 0.4 and not repetitive:
        # Boyer-Moore for medium-length, low-GC, non-repetitive patterns
        return boyer_moore_search(text, pattern)
    elif repetitive or gc_content > 0.6:
        # KMP for repetitive patterns or high GC content
        return kmp_search(text, pattern)
    else:
        # Default to Boyer-Moore for other cases
        return boyer_moore_search(text, pattern)

# Enhanced repetitiveness detection
def is_highly_repetitive(pattern):
    if len(pattern) <= 10:
        return False
    
    # Check for simple repeats
    max_repeat = 1
    current_repeat = 1
    for i in range(1, len(pattern)):
        if pattern[i] == pattern[i-1]:
            current_repeat += 1
            max_repeat = max(max_repeat, current_repeat)
        else:
            current_repeat = 1
    
    # Consider pattern repetitive if it has long repeats
    if max_repeat > len(pattern) / 3:
        return True
    
    # Check for periodic patterns
    for period in range(2, min(20, len(pattern)//2)):
        is_periodic = True
        for i in range(period, len(pattern)):
            if pattern[i] != pattern[i % period]:
                is_periodic = False
                break
        if is_periodic:
            return True
    
    return False

# Generate DNA sequence for testing
def generate_dna_sequence(length, gc_content=0.5):
    gc = int(gc_content * length)
    at = length - gc
    seq = ['G', 'C'] * (gc // 2) + ['A', 'T'] * (at // 2)
    if gc % 2 == 1:
        seq.append('G')
    if at % 2 == 1:
        seq.append('A')
    random.shuffle(seq)
    return ''.join(seq)

# Enhanced benchmarking function
def benchmark_algorithms(text, patterns, num_runs=7):
    algorithms = {
        'KMP': kmp_search,
        'Boyer-Moore': boyer_moore_search,
        'Bit-Parallel': bit_parallel_search,
        'Hybrid': hybrid_search
    }
    
    # Initialize results dictionary
    results = {}
    for pattern_desc, pattern in patterns.items():
        results[pattern_desc] = {
            'pattern_length': len(pattern),
            'gc_content': (pattern.count('G') + pattern.count('C')) / len(pattern),
            'repetitive': is_highly_repetitive(pattern),
            'algorithms': {algo: {'times': [], 'avg_time': 0, 'std_dev': 0} for algo in algorithms}
        }
        
        # Run each algorithm multiple times
        for algo_name, algo_func in algorithms.items():
            # Skip bit-parallel for patterns longer than 64 bp
            if algo_name == 'Bit-Parallel' and len(pattern) > 64:
                results[pattern_desc]['algorithms'][algo_name] = 'N/A'
                continue
                
            run_times = []
            for _ in range(num_runs):
                start_time = time.perf_counter()
                algo_func(text, pattern)
                end_time = time.perf_counter()
                run_times.append(end_time - start_time)
            
            # Store detailed timing information
            results[pattern_desc]['algorithms'][algo_name]['times'] = run_times
            results[pattern_desc]['algorithms'][algo_name]['avg_time'] = sum(run_times) / num_runs
            results[pattern_desc]['algorithms'][algo_name]['std_dev'] = np.std(run_times) if num_runs > 1 else 0
    
    return results

# Generate comprehensive test patterns
def generate_test_patterns():
    patterns = {}
    
    # Very short patterns (1-20 bp)
    patterns['Very Short Low GC'] = generate_dna_sequence(15, gc_content=0.2)
    patterns['Very Short High GC'] = generate_dna_sequence(15, gc_content=0.8)
    patterns['Very Short Repetitive'] = "AT" * 7 + "A"  # 15 bp
    
    # Short patterns (21-64 bp) - Bit-Parallel optimal range
    patterns['Short Low GC'] = generate_dna_sequence(50, gc_content=0.3)
    patterns['Short High GC'] = generate_dna_sequence(50, gc_content=0.7)
    patterns['Short Repetitive'] = "AT" * 25  # 50 bp
    
    # Medium patterns (65-200 bp)
    patterns['Medium Low GC'] = generate_dna_sequence(150, gc_content=0.3)
    patterns['Medium High GC'] = generate_dna_sequence(150, gc_content=0.7)
    patterns['Medium Repetitive'] = "ATCG" * 37 + "AT"  # 150 bp
    
    # Long patterns (201-500 bp)
    patterns['Long Low GC'] = generate_dna_sequence(350, gc_content=0.3)
    patterns['Long High GC'] = generate_dna_sequence(350, gc_content=0.7)
    patterns['Long Repetitive'] = "ATCG" * 87 + "AT"  # 350 bp
    
    # Very long patterns (501-1000 bp)
    patterns['Very Long Low GC'] = generate_dna_sequence(800, gc_content=0.3)
    patterns['Very Long High GC'] = generate_dna_sequence(800, gc_content=0.7)
    patterns['Very Long Repetitive'] = "ATCG" * 200  # 800 bp
    
    return patterns

# Create comprehensive results table
def create_results_table(results):
    table_data = []
    headers = ["Pattern Type", "Length", "GC Content", "Repetitive", "KMP (ms)", "Boyer-Moore (ms)", "Bit-Parallel (ms)", "Hybrid (ms)", "Best Algorithm", "Speedup"]
    
    for pattern_desc, data in results.items():
        # Get algorithm times
        kmp_data = data['algorithms']['KMP']
        bm_data = data['algorithms']['Boyer-Moore']
        bp_data = data['algorithms']['Bit-Parallel']
        hybrid_data = data['algorithms']['Hybrid']
        
        # Convert to milliseconds and format
        kmp_ms = f"{kmp_data['avg_time'] * 1000:.4f}" if kmp_data != 'N/A' else 'N/A'
        bm_ms = f"{bm_data['avg_time'] * 1000:.4f}" if bm_data != 'N/A' else 'N/A'
        bp_ms = f"{bp_data['avg_time'] * 1000:.4f}" if bp_data != 'N/A' else 'N/A'
        hybrid_ms = f"{hybrid_data['avg_time'] * 1000:.4f}" if hybrid_data != 'N/A' else 'N/A'
        
        # Determine which algorithm was fastest
        valid_times = {}
        if kmp_data != 'N/A':
            valid_times['KMP'] = kmp_data['avg_time']
        if bm_data != 'N/A':
            valid_times['Boyer-Moore'] = bm_data['avg_time']
        if bp_data != 'N/A':
            valid_times['Bit-Parallel'] = bp_data['avg_time']
            
        if valid_times:
            best_algo = min(valid_times, key=valid_times.get)
            best_time = valid_times[best_algo]
            
            # Check if hybrid was the best
            if hybrid_data != 'N/A' and hybrid_data['avg_time'] <= best_time:
                best_algo = 'Hybrid'
                best_time = hybrid_data['avg_time']
            
            # Calculate speedup of hybrid over best individual algorithm
            if hybrid_data != 'N/A' and best_algo != 'Hybrid':
                speedup = valid_times[best_algo] / hybrid_data['avg_time']
                speedup_str = f"{speedup:.2f}x"
            else:
                speedup_str = "N/A"
        else:
            best_algo = 'N/A'
            speedup_str = 'N/A'
        
        table_data.append([
            pattern_desc,
            data['pattern_length'],
            f"{data['gc_content']:.2%}",
            "Yes" if data['repetitive'] else "No",
            kmp_ms,
            bm_ms,
            bp_ms,
            hybrid_ms,
            best_algo,
            speedup_str
        ])
    
    return tabulate(table_data, headers=headers, tablefmt="grid")

# Create enhanced visualization
def create_enhanced_visualization(results):
    # Prepare data for plotting
    pattern_types = list(results.keys())
    algorithms = ['KMP', 'Boyer-Moore', 'Bit-Parallel', 'Hybrid']
    
    # Create a 2x2 grid of subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # Plot 1: Execution time by pattern type
    for algo in algorithms:
        times = []
        for pattern_desc in pattern_types:
            algo_data = results[pattern_desc]['algorithms'][algo]
            if algo_data == 'N/A':
                times.append(0)
            else:
                times.append(algo_data['avg_time'] * 1000)
        
        ax1.plot(pattern_types, times, marker='o', label=algo, linewidth=2)
    
    ax1.set_title('Execution Time by Pattern Type', fontsize=14, fontweight='bold')
    ax1.set_xlabel('Pattern Type')
    ax1.set_ylabel('Time (ms)')
    ax1.tick_params(axis='x', rotation=45)
    ax1.legend()
    ax1.grid(True, linestyle='--', alpha=0.7)
    
    # Plot 2: Speedup of Hybrid over other algorithms
    speedup_data = {algo: [] for algo in algorithms if algo != 'Hybrid'}
    
    for pattern_desc in pattern_types:
        hybrid_data = results[pattern_desc]['algorithms']['Hybrid']
        if hybrid_data == 'N/A':
            continue
            
        hybrid_time = hybrid_data['avg_time']
        for algo in algorithms:
            if algo == 'Hybrid':
                continue
                
            algo_data = results[pattern_desc]['algorithms'][algo]
            if algo_data != 'N/A':
                algo_time = algo_data['avg_time']
                speedup = algo_time / hybrid_time if hybrid_time > 0 else 0
                speedup_data[algo].append(speedup)
            else:
                speedup_data[algo].append(0)
    
    x = np.arange(len(pattern_types))
    width = 0.25
    multiplier = 0
    
    for algo, speedups in speedup_data.items():
        offset = width * multiplier
        rects = ax2.bar(x + offset, speedups, width, label=algo)
        ax2.bar_label(rects, fmt='%.1fx', padding=3, fontsize=8)
        multiplier += 1
    
    ax2.axhline(y=1, color='r', linestyle='--', alpha=0.7)
    ax2.set_title('Speedup of Hybrid over Other Algorithms', fontsize=14, fontweight='bold')
    ax2.set_xlabel('Pattern Type')
    ax2.set_ylabel('Speedup (x times faster)')
    ax2.set_xticks(x + width, pattern_types)
    ax2.tick_params(axis='x', rotation=45)
    ax2.legend(loc='upper left')
    ax2.grid(True, linestyle='--', alpha=0.7)
    
    # Plot 3: Algorithm performance by pattern length
    pattern_lengths = [results[pt]['pattern_length'] for pt in pattern_types]
    algo_performance = {algo: [] for algo in algorithms}
    
    for pattern_desc in pattern_types:
        for algo in algorithms:
            algo_data = results[pattern_desc]['algorithms'][algo]
            if algo_data == 'N/A':
                algo_performance[algo].append(0)
            else:
                algo_performance[algo].append(algo_data['avg_time'] * 1000)
    
    for algo in algorithms:
        ax3.scatter(pattern_lengths, algo_performance[algo], label=algo, s=100, alpha=0.7)
    
    ax3.set_title('Performance vs Pattern Length', fontsize=14, fontweight='bold')
    ax3.set_xlabel('Pattern Length (bp)')
    ax3.set_ylabel('Time (ms)')
    ax3.legend()
    ax3.grid(True, linestyle='--', alpha=0.7)
    
    # Plot 4: Algorithm performance by GC content
    gc_contents = [results[pt]['gc_content'] for pt in pattern_types]
    
    for algo in algorithms:
        ax4.scatter(gc_contents, algo_performance[algo], label=algo, s=100, alpha=0.7)
    
    ax4.set_title('Performance vs GC Content', fontsize=14, fontweight='bold')
    ax4.set_xlabel('GC Content')
    ax4.set_ylabel('Time (ms)')
    ax4.legend()
    ax4.grid(True, linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig('enhanced_benchmark_results.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Create a summary table of when each algorithm wins
    create_algorithm_summary(results)

def create_algorithm_summary(results):
    algo_wins = {'KMP': 0, 'Boyer-Moore': 0, 'Bit-Parallel': 0, 'Hybrid': 0}
    algo_conditions = {}
    
    for pattern_desc, data in results.items():
        valid_times = {}
        for algo in ['KMP', 'Boyer-Moore', 'Bit-Parallel', 'Hybrid']:
            algo_data = data['algorithms'][algo]
            if algo_data != 'N/A':
                valid_times[algo] = algo_data['avg_time']
        
        if valid_times:
            best_algo = min(valid_times, key=valid_times.get)
            algo_wins[best_algo] += 1
            
            # Record the conditions for this win
            if best_algo not in algo_conditions:
                algo_conditions[best_algo] = []
            
            algo_conditions[best_algo].append({
                'pattern_type': pattern_desc,
                'length': data['pattern_length'],
                'gc_content': data['gc_content'],
                'repetitive': data['repetitive'],
                'time': valid_times[best_algo] * 1000
            })
    
    # Print summary
    print("\nAlgorithm Performance Summary:")
    print("=" * 50)
    for algo, wins in algo_wins.items():
        print(f"{algo}: {wins} wins")
    
    print("\nConditions where each algorithm performed best:")
    for algo, conditions in algo_conditions.items():
        print(f"\n{algo}:")
        for cond in conditions:
            print(f"  - {cond['pattern_type']} (Length: {cond['length']}, GC: {cond['gc_content']:.2%}, "
                  f"Repetitive: {cond['repetitive']}, Time: {cond['time']:.4f} ms)")

# Main function
if __name__ == "__main__":
    # Generate a long text (1Mbp bacterial genome)
    print("Generating test text (1Mbp bacterial genome)...")
    text = generate_dna_sequence(1000000, gc_content=0.5)
    
    # Generate test patterns
    print("Generating test patterns...")
    patterns = generate_test_patterns()
    
    # Run benchmarks
    print("Running benchmarks...")
    results = benchmark_algorithms(text, patterns)
    
    # Create and display results table
    print("\nBenchmark Results:")
    table = create_results_table(results)
    print(table)
    
    # Create enhanced visualization
    print("\nGenerating enhanced visualization...")
    create_enhanced_visualization(results)
    
    # Save detailed results to CSV
    print("Saving detailed results to CSV...")
    detailed_data = []
    for pattern_desc, data in results.items():
        row = {
            'Pattern Type': pattern_desc,
            'Length': data['pattern_length'],
            'GC Content': data['gc_content'],
            'Repetitive': data['repetitive'],
        }
        
        for algo in ['KMP', 'Boyer-Moore', 'Bit-Parallel', 'Hybrid']:
            algo_data = data['algorithms'][algo]
            if algo_data == 'N/A':
                row[algo] = 'N/A'
                row[f'{algo}_StdDev'] = 'N/A'
            else:
                row[algo] = algo_data['avg_time'] * 1000
                row[f'{algo}_StdDev'] = algo_data['std_dev'] * 1000
        
        detailed_data.append(row)
    
    df = pd.DataFrame(detailed_data)
    df.to_csv('enhanced_benchmark_results.csv', index=False)
    print("Done! Results saved to 'enhanced_benchmark_results.csv' and 'enhanced_benchmark_results.png'")