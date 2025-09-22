import pandas as pd
import matplotlib.pyplot as plt
import os

# Base folder for your results
base_dir = "algo_results"

# Pick file
bm_file = os.path.join(base_dir, "bm_results.csv")

df = pd.read_csv(bm_file)

# Compute speedup
df["Speedup"] = df["SerialTime"] / df["ParallelTime"]

plt.figure(figsize=(8,5))
for length in df["PatternLength"].unique():
    subset = df[df["PatternLength"] == length]
    plt.plot(subset["Threads"], subset["Speedup"], marker="o", label=f"Length {length}")

plt.xlabel("Number of Threads")
plt.ylabel("Speedup (SerialTime / ParallelTime)")
plt.title("BM: Speedup vs Threads")
plt.legend()
plt.grid(True)
plt.show()
