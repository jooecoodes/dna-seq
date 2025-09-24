import pandas as pd

# Load all result CSVs
files = {
    "BMH": "algo_results/bmh_results.csv",
    "BP": "algo_results/bp_results.csv",
    "KMP": "algo_results/kmp_results.csv"
}

# Read and tag with algorithm name
dfs = []
for algo, path in files.items():
    df = pd.read_csv(path)
    df["algorithm"] = algo
    dfs.append(df)

# Merge into one big dataframe
data = pd.concat(dfs, ignore_index=True)

# Group by experiment setup (ignoring algorithm)
group_cols = ["length", "gc_content", "entropy", "matches"]

def choose_best(group):
    # Best algorithm by lowest parallel_time
    best_row = group.loc[group["parallel_time"].idxmin()]
    
    # Collect all parallel times
    times = group.set_index("algorithm")["parallel_time"].to_dict()
    
    # Sort algorithms by time
    sorted_times = sorted(times.items(), key=lambda x: x[1])
    
    best_time = sorted_times[0][1]
    second_best_time = sorted_times[1][1] if len(sorted_times) > 1 else best_time
    
    margin = second_best_time / best_time if best_time > 0 else None
    
    return pd.Series({
        "length": best_row["length"],
        "gc_content": best_row["gc_content"],
        "entropy": best_row["entropy"],
        "matches": best_row["matches"],
        "best_algorithm": best_row["algorithm"],
        "BMH_parallel_time": times.get("BMH", None),
        "BP_parallel_time": times.get("BP", None),
        "KMP_parallel_time": times.get("KMP", None),
        "margin_of_victory": margin
    })

# Apply grouping
labeled = data.groupby(group_cols).apply(choose_best).reset_index(drop=True)

# Save the labeled dataset
labeled.to_csv("training_data/labeled_best_algorithms.csv", index=False)

print("✅ Labeled dataset saved with margin_of_victory to training_data/labeled_best_algorithms.csv")
