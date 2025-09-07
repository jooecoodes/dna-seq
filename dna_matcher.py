import ctypes
import time
import random
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import platform
from tabulate import tabulate

# Load the shared library
if platform.system() == "Windows":
    lib_name = "dna_matchers.dll"
else:
    lib_name = "dna_matchers.so"

lib_path = os.path.join(os.path.dirname(__file__), lib_name)
lib = ctypes.CDLL(lib_path)

# Define function prototypes
lib.kmp_search_c.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.POINTER(ctypes.c_int)]
lib.kmp_search_c.restype = ctypes.POINTER(ctypes.c_int)

lib.boyer_moore_search_c.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.POINTER(ctypes.c_int)]
lib.boyer_moore_search_c.restype = ctypes.POINTER(ctypes.c_int)

lib.bit_parallel_search_c.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.POINTER(ctypes.c_int)]
lib.bit_parallel_search_c.restype = ctypes.POINTER(ctypes.c_int)

lib.hybrid_search_cpp_c.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.POINTER(ctypes.c_int)]
lib.hybrid_search_cpp_c.restype = ctypes.POINTER(ctypes.c_int)

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

def hybrid_search_cpp(text, pattern):
    result_size = ctypes.c_int()
    result_ptr = lib.hybrid_search_cpp_c(text.encode(), pattern.encode(), ctypes.byref(result_size))
    results = [result_ptr[i] for i in range(result_size.value)]
    lib.free_result(result_ptr)
    return results

# Python-based hybrid selector for comparison
def hybrid_search_python(text, pattern):
    length = len(pattern)
    if length == 0:
        return []
    
    # Calculate GC content
    gc_count = pattern.count('G') + pattern.count('C')
    gc_content = gc_count / length
    
    # Check repetitiveness
    repetitive = is_highly_repetitive(pattern)
    
    # Algorithm selection logic
    if length <= 64 and not repetitive:
        return bit_parallel_search(text, pattern)
    elif repetitive or gc_content > 0.6:
        return kmp_search(text, pattern)
    else:
        return boyer_moore_search(text, pattern)

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

# Benchmarking function
def benchmark_algorithms(text, patterns, num_runs=7):
    algorithms = {
        'KMP': kmp_search,
        'Boyer-Moore': boyer_moore_search,
        'Bit-Parallel': bit_parallel_search,
        'Hybrid (Python)': hybrid_search_python,
        'Hybrid (C++)': hybrid_search_cpp
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
                # Warm-up run (not measured)
                algo_func(text, pattern)
                
                # Measured run
                start_time = time.perf_counter()
                algo_func(text, pattern)
                end_time = time.perf_counter()
                run_times.append(end_time - start_time)
            
            # Store median time (more robust than average)
            run_times.sort()
            median_time = run_times[num_runs // 2]
            results[pattern_desc]['algorithms'][algo_name]['avg_time'] = median_time
            results[pattern_desc]['algorithms'][algo_name]['std_dev'] = np.std(run_times) if num_runs > 1 else 0
    
    return results

# Generate comprehensive test patterns
def generate_test_patterns():
    patterns = {}
    
    # Test patterns of various lengths and characteristics
    lengths = [20, 50, 100, 200, 500, 800, 1000]
    gc_levels = [0.3, 0.5, 0.7]
    
    for length in lengths:
        for gc in gc_levels:
            # Non-repetitive patterns
            patterns[f'Len{length}_GC{gc}_NonRep'] = generate_dna_sequence(length, gc_content=gc)
            
            # Repetitive patterns
            if length <= 100:
                patterns[f'Len{length}_GC{gc}_Rep'] = "AT" * (length // 2)
            else:
                patterns[f'Len{length}_GC{gc}_Rep'] = "ATCG" * (length // 4)
    
    return patterns

# Create results table
def create_results_table(results):
    table_data = []
    headers = ["Pattern", "Length", "GC", "Rep", "KMP (μs)", "BM (μs)", "BP (μs)", "Hybrid-Py (μs)", "Hybrid-C++ (μs)", "Best Algo", "C++ Speedup"]
    
    for pattern_desc, data in results.items():
        # Get algorithm times
        kmp_time = data['algorithms']['KMP']['avg_time']
        bm_time = data['algorithms']['Boyer-Moore']['avg_time']
        bp_time = data['algorithms']['Bit-Parallel']['avg_time'] if data['algorithms']['Bit-Parallel'] != 'N/A' else float('inf')
        hybrid_py_time = data['algorithms']['Hybrid (Python)']['avg_time']
        hybrid_cpp_time = data['algorithms']['Hybrid (C++)']['avg_time']
        
        # Convert to microseconds and format
        kmp_us = f"{kmp_time * 1e6:.1f}"
        bm_us = f"{bm_time * 1e6:.1f}"
        bp_us = f"{bp_time * 1e6:.1f}" if bp_time != float('inf') else 'N/A'
        hybrid_py_us = f"{hybrid_py_time * 1e6:.1f}"
        hybrid_cpp_us = f"{hybrid_cpp_time * 1e6:.1f}"
        
        # Determine which algorithm was fastest
        times = {
            'KMP': kmp_time,
            'Boyer-Moore': bm_time,
            'Bit-Parallel': bp_time if bp_time != float('inf') else float('inf'),
            'Hybrid (Python)': hybrid_py_time,
            'Hybrid (C++)': hybrid_cpp_time
        }
        
        # Remove infinite values
        valid_times = {k: v for k, v in times.items() if v != float('inf')}
        best_algo = min(valid_times, key=valid_times.get)
        
        # Calculate speedup of C++ hybrid over Python hybrid
        speedup = hybrid_py_time / hybrid_cpp_time if hybrid_cpp_time > 0 else 0
        speedup_str = f"{speedup:.2f}x" if speedup > 0 else "N/A"
        
        table_data.append([
            pattern_desc,
            data['pattern_length'],
            f"{data['gc_content']:.2f}",
            "Yes" if data['repetitive'] else "No",
            kmp_us,
            bm_us,
            bp_us,
            hybrid_py_us,
            hybrid_cpp_us,
            best_algo,
            speedup_str
        ])
    
    return tabulate(table_data, headers=headers, tablefmt="grid")

# Create visualization focused on hybrid performance
def create_hybrid_visualization(results):
    # Prepare data for plotting
    pattern_types = list(results.keys())
    
    # Create a 2x2 grid of subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # Plot 1: Execution times comparison
    algo_times = {
        'KMP': [],
        'Boyer-Moore': [],
        'Bit-Parallel': [],
        'Hybrid (Python)': [],
        'Hybrid (C++)': []
    }
    
    for pattern_desc in pattern_types:
        data = results[pattern_desc]
        algo_times['KMP'].append(data['algorithms']['KMP']['avg_time'] * 1e6)
        algo_times['Boyer-Moore'].append(data['algorithms']['Boyer-Moore']['avg_time'] * 1e6)
        
        bp_time = data['algorithms']['Bit-Parallel']
        algo_times['Bit-Parallel'].append(bp_time['avg_time'] * 1e6 if bp_time != 'N/A' else 0)
        
        algo_times['Hybrid (Python)'].append(data['algorithms']['Hybrid (Python)']['avg_time'] * 1e6)
        algo_times['Hybrid (C++)'].append(data['algorithms']['Hybrid (C++)']['avg_time'] * 1e6)
    
    x = np.arange(len(pattern_types))
    width = 0.15
    
    for i, (algo, times) in enumerate(algo_times.items()):
        offset = width * (i - 2)  # Center the bars
        ax1.bar(x + offset, times, width, label=algo)
    
    ax1.set_xlabel('Pattern Type')
    ax1.set_ylabel('Time (μs)')
    ax1.set_title('Execution Time by Algorithm')
    ax1.set_xticks(x)
    ax1.set_xticklabels(pattern_types, rotation=45, ha='right')
    ax1.legend()
    ax1.grid(True, linestyle='--', alpha=0.7)
    
    # Plot 2: Speedup of C++ hybrid over Python hybrid
    speedups = []
    for pattern_desc in pattern_types:
        py_time = results[pattern_desc]['algorithms']['Hybrid (Python)']['avg_time']
        cpp_time = results[pattern_desc]['algorithms']['Hybrid (C++)']['avg_time']
        speedups.append(py_time / cpp_time if cpp_time > 0 else 0)
    
    ax2.bar(pattern_types, speedups)
    ax2.axhline(y=1, color='r', linestyle='--')
    ax2.set_xlabel('Pattern Type')
    ax2.set_ylabel('Speedup (C++ / Python)')
    ax2.set_title('Speedup of C++ Hybrid over Python Hybrid')
    ax2.set_xticklabels(pattern_types, rotation=45, ha='right')
    ax2.grid(True, linestyle='--', alpha=0.7)
    
    # Plot 3: Performance by pattern length
    lengths = [data['pattern_length'] for data in results.values()]
    hybrid_py_times = [data['algorithms']['Hybrid (Python)']['avg_time'] * 1e6 for data in results.values()]
    hybrid_cpp_times = [data['algorithms']['Hybrid (C++)']['avg_time'] * 1e6 for data in results.values()]
    
    ax3.scatter(lengths, hybrid_py_times, label='Python Hybrid', alpha=0.7)
    ax3.scatter(lengths, hybrid_cpp_times, label='C++ Hybrid', alpha=0.7)
    ax3.set_xlabel('Pattern Length')
    ax3.set_ylabel('Time (μs)')
    ax3.set_title('Hybrid Performance by Pattern Length')
    ax3.legend()
    ax3.grid(True, linestyle='--', alpha=0.7)
    
    # Plot 4: Performance by GC content
    gc_contents = [data['gc_content'] for data in results.values()]
    ax4.scatter(gc_contents, hybrid_py_times, label='Python Hybrid', alpha=0.7)
    ax4.scatter(gc_contents, hybrid_cpp_times, label='C++ Hybrid', alpha=0.7)
    ax4.set_xlabel('GC Content')
    ax4.set_ylabel('Time (μs)')
    ax4.set_title('Hybrid Performance by GC Content')
    ax4.legend()
    ax4.grid(True, linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig('hybrid_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()

# Main function
if __name__ == "__main__":
    # Generate a long text (1Mbp bacterial genome)
    print("Generating test text (1Mbp bacterial genome)...")
    text = generate_dna_sequence(1000000, gc_content=0.5)
    
    # Generate test patterns
    print("Generating test patterns...")
    patterns = generate_test_patterns()
    
    # Run benchmarks
    print("Running benchmarks (this may take a while)...")
    results = benchmark_algorithms(text, patterns)
    
    # Create and display results table
    print("\nBenchmark Results:")
    table = create_results_table(results)
    print(table)
    
    # Create visualization
    print("\nGenerating visualization...")
    create_hybrid_visualization(results)
    
    # Calculate overall statistics
    total_patterns = len(patterns)
    cpp_faster_count = 0
    speedups = []
    
    for pattern_desc, data in results.items():
        py_time = data['algorithms']['Hybrid (Python)']['avg_time']
        cpp_time = data['algorithms']['Hybrid (C++)']['avg_time']
        
        if cpp_time < py_time:
            cpp_faster_count += 1
            speedups.append(py_time / cpp_time)
    
    print(f"\nOverall Statistics:")
    print(f"C++ hybrid was faster than Python hybrid: {cpp_faster_count}/{total_patterns} ({cpp_faster_count/total_patterns*100:.1f}%)")
    if speedups:
        print(f"Average speedup when C++ was faster: {np.mean(speedups):.3f}x")
        print(f"Median speedup when C++ was faster: {np.median(speedups):.3f}x")
    
    # Save detailed results to CSV
    print("Saving detailed results to CSV...")
    detailed_data = []
    for pattern_desc, data in results.items():
        row = {
            'Pattern': pattern_desc,
            'Length': data['pattern_length'],
            'GC_Content': data['gc_content'],
            'Repetitive': data['repetitive'],
        }
        
        for algo in ['KMP', 'Boyer-Moore', 'Bit-Parallel', 'Hybrid (Python)', 'Hybrid (C++)']:
            algo_data = data['algorithms'][algo]
            if algo_data == 'N/A':
                row[algo] = 'N/A'
                row[f'{algo}_StdDev'] = 'N/A'
            else:
                row[algo] = algo_data['avg_time'] * 1e6
                row[f'{algo}_StdDev'] = algo_data['std_dev'] * 1e6
        
        detailed_data.append(row)
    
    df = pd.DataFrame(detailed_data)
    df.to_csv('hybrid_comparison_results.csv', index=False)
    print("Done! Results saved to 'hybrid_comparison_results.csv' and 'hybrid_comparison.png'")