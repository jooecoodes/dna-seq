import pandas as pd
import matplotlib.pyplot as plt
import os

# Base folder for your results
base_dir = "algo_results"

# Pick file
bm_file = os.path.join(base_dir, "bm_results.csv")

df = pd.read_csv(bm_file)

plt.figure(figsize=(8,5))

# Serial
serial = df[df["Threads"] == 1]
plt.plot(serial["PatternLength"], serial["SerialTime"], marker="o", linestyle="--", label="Serial")

# Parallel (max threads, e.g., 8)
parallel = df[df["Threads"] == df["Threads"].max()]
plt.plot(parallel["PatternLength"], parallel["ParallelTime"], marker="s", linestyle="-", label=f"Parallel ({df['Threads'].max()} threads)")

plt.xlabel("Pattern Length")
plt.ylabel("Execution Time (ms)")
plt.title("BM: Execution Time vs Pattern Length")
plt.legend()
plt.grid(True)
plt.show()
