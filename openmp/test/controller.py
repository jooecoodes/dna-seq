import math
import subprocess
import os
import csv
import sys

# ------------------------------
# Feature Extraction
# ------------------------------
def compute_features(sequence: str):
    """Compute DNA sequence features: length, GC content, entropy."""
    length = len(sequence)
    if length == 0:
        return [0, 0.0, 0.0]

    gc_content = (sequence.count("G") + sequence.count("C")) / length

    # Shannon entropy for repetitiveness
    freqs = [sequence.count(base) / length for base in "ACGT" if sequence.count(base) > 0]
    entropy = -sum(f * math.log2(f) for f in freqs)

    return [length, gc_content, entropy]


# ------------------------------
# Algorithm Selection Heuristic
# ------------------------------
def choose_algorithm(length, gc_content, entropy):
    """Choose the best algorithm based on sequence features."""
    # Rule 1: Very short patterns → Bit-Parallel (BP)
    if length <= 64:
        return "BP"

    # Rule 2: Medium patterns (65–256)
    if length <= 256:
        if entropy <= 1.1:
            return "KMP"   # better for repetitive patterns
        else:
            return "BMH"   # better for diverse patterns

    # Rule 3: Longer patterns (≥ 512)
    if length >= 512:
        if entropy >= 1.1:
            return "BMH"
        else:
            return "KMP"

    # Rule 4: Adjust for extreme GC content
    if gc_content <= 0.2 or gc_content >= 0.8:
        return "KMP"

    # Default
    return "BMH"


# ------------------------------
# Run Chosen Executable
# ------------------------------
def run_algorithm(algo, genome_file, pattern):
    """Run the chosen algorithm executable with genome and pattern."""
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "CPP"))
    exe_map = {
        "BMH": os.path.join(base_dir, "bmh.exe"),
        "KMP": os.path.join(base_dir, "kmp.exe"),
        "BP":  os.path.join(base_dir, "bp.exe"),
    }
    exe_file = exe_map[algo]

    try:
        result = subprocess.run(
            [exe_file, genome_file, pattern],
            capture_output=True,
            text=True,
            check=True
        )
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        return f"Error running {exe_file}: {e.stderr}"


# ------------------------------
# Main Pipeline with CSV output
# ------------------------------
if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python controller.py <genome_file> <patterns.txt> <output.csv>")
        sys.exit(1)

    genome_file = sys.argv[1]
    patterns_file = sys.argv[2]
    output_csv = sys.argv[3]

    with open(patterns_file, "r") as pf, open(output_csv, "w", newline="") as cf:
        writer = csv.writer(cf)
        writer.writerow(["pattern", "length", "gc_content", "entropy", "chosen_algorithm", "result"])

        for line in pf:
            seq = line.strip()
            if not seq:
                continue

            # Step 1: Compute features
            length, gc_content, entropy = compute_features(seq)

            # Step 2: Choose algorithm
            algo = choose_algorithm(length, gc_content, entropy)

            # Step 3: Run algorithm executable
            result = run_algorithm(algo, genome_file, seq)

            # Save row to CSV
            writer.writerow([seq, length, gc_content, entropy, algo, result])

            print(f"Pattern: {seq[:20]}... | Algo: {algo} | Done ✅")
