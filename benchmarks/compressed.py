from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

def plot_results(filename, algo_name=None):
    # Load file
    base_dir = Path("algo_results")
    file = base_dir / filename
    df = pd.read_csv(file)

    # Auto-detect algorithm name from file if not given
    if algo_name is None:
        algo_name = filename.split("_")[0].upper()

    # Compute metrics
    df["Speedup"] = df["SerialTime"] / df["ParallelTime"]
    df["Efficiency"] = df["Speedup"] / df["Threads"]

    # Create subplots
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # --- 1. Speedup vs Threads ---
    for length in df["PatternLength"].unique():
        subset = df[df["PatternLength"] == length]
        axes[0].plot(subset["Threads"], subset["Speedup"], marker="o", label=f"Len {length}")
    axes[0].set_title(f"{algo_name}: Speedup vs Threads")
    axes[0].set_xlabel("Threads")
    axes[0].set_ylabel("Speedup")
    axes[0].legend()
    axes[0].grid(True)

    # --- 2. Execution Time vs Pattern Length ---
    serial = df[df["Threads"] == 1]
    axes[1].plot(serial["PatternLength"], serial["SerialTime"], marker="o", linestyle="--", label="Serial")
    max_threads = df["Threads"].max()
    parallel = df[df["Threads"] == max_threads]
    axes[1].plot(parallel["PatternLength"], parallel["ParallelTime"], marker="s", linestyle="-", label=f"Parallel ({max_threads} threads)")
    axes[1].set_title(f"{algo_name}: Execution Time vs Pattern Length")
    axes[1].set_xlabel("Pattern Length")
    axes[1].set_ylabel("Time (ms)")
    axes[1].legend()
    axes[1].grid(True)

    # --- 3. Efficiency vs Threads ---
    for length in df["PatternLength"].unique():
        subset = df[df["PatternLength"] == length]
        axes[2].plot(subset["Threads"], subset["Efficiency"]*100, marker="o", label=f"Len {length}")
    axes[2].set_title(f"{algo_name}: Efficiency vs Threads")
    axes[2].set_xlabel("Threads")
    axes[2].set_ylabel("Efficiency (%)")
    axes[2].legend()
    axes[2].grid(True)

    plt.tight_layout()
    plt.show()


# Example usage:
# For BM results
# plot_results("bm_results.csv", "BM")

# For BP results
# plot_results("bp_results.csv", "BP")

# For KMP results
plot_results("kmp_results.csv", "KMP")
