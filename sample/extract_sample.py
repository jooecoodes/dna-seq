#!/usr/bin/env python3
"""
extract_snippets_by_length.py

Extract N snippets per length from a genome FASTA (ecoli.fasta).
Sampling methods: 'systematic' (default), 'stratified', 'random_quota'.

Outputs:
  - snippets.fasta
  - snippets_metadata.csv

Each row in CSV: id, length, start, method, gc, entropy, sequence
"""

import random
import math
import csv
from collections import Counter
from pathlib import Path

# --------------- CONFIG ---------------
FASTA = "../dna/streptomyces_coelicolor.fasta"  # input genome FASTA
OUT_FASTA = "fasta/snippets.fasta"
OUT_CSV = "csv/snippets_metadata.csv"

LENGTHS = [64, 128, 256, 512, 1000, 2000]
PER_LENGTH_K = 9               # 9 per length -> 54 total
SAMPLING_METHOD = "systematic" # "systematic" | "stratified" | "random_quota"
MIN_DISTANCE_FACTOR = 1.0      # used in random_quota: min separation = factor * length
RANDOM_SEED = 42
ALLOW_OVERLAP = False          # if False, try to avoid overlapping snippets for same length
MAX_TRIES = 1000               # attempts when resolving conflicts / ambiguous bases
# --------------------------------------

random.seed(RANDOM_SEED)

# ---------- helpers ----------
def read_fasta_one_seq(path):
    seq_parts = []
    with open(path, "r") as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq_parts.append(line.strip().upper())
    return "".join(seq_parts)

def gc_content(seq):
    if len(seq) == 0:
        return 0.0
    g = seq.count("G")
    c = seq.count("C")
    return (g + c) / len(seq)

def shannon_entropy(seq):
    if len(seq) == 0:
        return 0.0
    counts = Counter(seq)
    total = len(seq)
    probs = [c / total for c in counts.values() if c > 0]
    return -sum(p * math.log2(p) for p in probs if p > 0)

def is_valid_seq(seq):
    # skip ambiguous or N bases
    return all(ch in "ACGT" for ch in seq)

# ---------- sampling functions ----------
def systematic_positions(genome_len, L, k):
    """Return k start positions equally spaced across genome."""
    if k <= 0:
        return []
    step = max(1, (genome_len - L) // k)
    positions = []
    # start with small offset for better coverage
    offset = (step // 2)
    for i in range(k):
        s = offset + i * step
        if s < 0: s = 0
        if s > genome_len - L: s = max(0, genome_len - L)
        positions.append(s)
    # deduplicate and clamp
    positions = sorted(list(dict.fromkeys(positions)))
    return positions[:k]

def stratified_positions(genome_len, L, k):
    """Divide genome into k strata and pick one random start inside each stratum."""
    if k <= 0:
        return []
    stratum_size = max(1, (genome_len - L) // k)
    positions = []
    for i in range(k):
        start_range = i * stratum_size
        end_range = min(start_range + stratum_size + (L-1), genome_len - L)
        if end_range <= start_range:
            s = min(start_range, max(0, genome_len - L))
        else:
            s = random.randint(start_range, end_range)
        positions.append(s)
    return positions[:k]

def random_quota_positions(genome_len, L, k, min_distance):
    """Random sampling enforcing a minimum distance between chosen starts."""
    if k <= 0:
        return []
    positions = []
    tries = 0
    while len(positions) < k and tries < MAX_TRIES:
        s = random.randint(0, max(0, genome_len - L))
        # check min distance
        if all(abs(s - p) >= min_distance for p in positions):
            positions.append(s)
        tries += 1
    # if we couldn't fill because genome small, fallback to systematic
    if len(positions) < k:
        return systematic_positions(genome_len, L, k)
    return positions

# ---------- main ----------
def main():
    genome = read_fasta_one_seq(FASTA)
    if not genome:
        raise SystemExit("Failed to read genome FASTA.")
    n = len(genome)
    print(f"Genome length: {n} bp")

    Path(OUT_FASTA).parent.mkdir(parents=True, exist_ok=True)
    Path(OUT_CSV).parent.mkdir(parents=True, exist_ok=True)

    uid = 0
    rows = []

    for L in LENGTHS:
        k = PER_LENGTH_K
        print(f"Sampling length={L}, k={k}, method={SAMPLING_METHOD}")

        if SAMPLING_METHOD == "systematic":
            starts = systematic_positions(n, L, k)
        elif SAMPLING_METHOD == "stratified":
            starts = stratified_positions(n, L, k)
        elif SAMPLING_METHOD == "random_quota":
            min_dist = int(max(1, MIN_DISTANCE_FACTOR * L))
            starts = random_quota_positions(n, L, k, min_dist)
        else:
            raise ValueError("Unknown sampling method")

        # If we want to avoid overlapping snippets with the same length, enforce min distance
        if not ALLOW_OVERLAP:
            # greedy filter to ensure distance >= L (non-overlap)
            filtered = []
            for s in sorted(starts):
                if all(abs(s - t) >= L for t in filtered):
                    filtered.append(s)
            # if filtering removed too many, fill from systematic fallback
            if len(filtered) < k:
                fallback = systematic_positions(n, L, k * 2)
                for s in fallback:
                    if all(abs(s - t) >= L for t in filtered):
                        filtered.append(s)
                        if len(filtered) >= k: break
            starts = filtered[:k]

        # validate and collect
        collected = 0
        for s in starts:
            if collected >= k:
                break
            seq = genome[s:s+L]
            if not is_valid_seq(seq):
                # try small forward search for a valid window (within MAX_TRIES)
                found = False
                for shift in range(1, min(MAX_TRIES, L)):
                    for sign in (+1, -1):
                        s2 = s + sign * shift
                        if 0 <= s2 <= n - L:
                            seq2 = genome[s2:s2+L]
                            if is_valid_seq(seq2) and all(abs(s2 - ex) >= L for ex in [r[2] for r in rows if r[0].startswith("snip") and r[1]==L]):
                                s = s2
                                seq = seq2
                                found = True
                                break
                    if found:
                        break
                if not found:
                    # skip this candidate
                    continue
            uid += 1
            sid = f"snip{uid:05d}"
            gc = gc_content(seq)
            ent = shannon_entropy(seq)
            rows.append((sid, L, s, SAMPLING_METHOD, gc, ent, seq))
            collected += 1

        # If still short (rare), try random fill
        tries = 0
        while collected < k and tries < MAX_TRIES:
            s = random.randint(0, max(0, n - L))
            seq = genome[s:s+L]
            if is_valid_seq(seq) and all(abs(s - r[2]) >= L for r in rows if r[1] == L):
                uid += 1
                sid = f"snip{uid:05d}"
                rows.append((sid, L, s, SAMPLING_METHOD, gc_content(seq), shannon_entropy(seq), seq))
                collected += 1
            tries += 1
        if collected < k:
            print(f"  [WARN] only collected {collected}/{k} for length {L}")

    # Write outputs
    with open(OUT_FASTA, "w") as fa, open(OUT_CSV, "w", newline='') as csvf:
        writer = csv.writer(csvf)
        writer.writerow(["id", "length", "start", "method", "gc", "entropy", "sequence"])
        for sid, L, start, method, gc, ent, seq in rows:
            fa.write(f">{sid} len={L} start={start} method={method} gc={gc:.4f} ent={ent:.4f}\n")
            fa.write(seq + "\n")
            writer.writerow([sid, L, start, method, f"{gc:.4f}", f"{ent:.4f}", seq])

    print(f"Done. Wrote {len(rows)} snippets to {OUT_FASTA} and {OUT_CSV}")

if __name__ == "__main__":
    main()
