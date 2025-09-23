import math
import random
import itertools

def make_dna(length, gc_content, repetitiveness):
    total_GC = round(length * gc_content)
    total_AT = length - total_GC
    
    best_skew = find_perfect_mix(gc_content, repetitiveness)
    
    num_G = round(best_skew * total_GC)
    num_C = total_GC - num_G
    num_A = round(best_skew * total_AT)
    num_T = total_AT - num_A

    blocks = ['G'] * int(num_G) + ['C'] * int(num_C) + ['A'] * int(num_A) + ['T'] * int(num_T)
    random.shuffle(blocks)
    
    return ''.join(blocks)

def find_perfect_mix(gc_content, target_entropy):
    best_skew = 0.5
    best_diff = abs(calculate_entropy(0.5, gc_content) - target_entropy)
    
    for test_skew in [i/10000.0 for i in range(0, 10001)]:
        current_entropy = calculate_entropy(test_skew, gc_content)
        diff = abs(current_entropy - target_entropy)
        if diff < best_diff:
            best_diff = diff
            best_skew = test_skew
            
    return best_skew

def calculate_entropy(skew, gc_content):
    fG = skew * gc_content
    fC = (1 - skew) * gc_content
    fA = skew * (1 - gc_content)
    fT = (1 - skew) * (1 - gc_content)
    
    entropy = 0.0
    for freq in [fG, fC, fA, fT]:
        if freq > 0:
            entropy -= freq * math.log2(freq)
    return entropy

# benchmark bins
lengths = [64, 128, 256, 512, 1000, 2000]
gc_contents = [0.2, 0.5, 0.8]
entropies = [0.3, 1.1, 1.9]

# === NEW PART: generate patterns automatically ===
patterns = [
    make_dna(length, gc, h)
    for length, gc, h in itertools.product(lengths, gc_contents, entropies)
]

with open("patterns.txt", "w") as f:
    for p in patterns:
        f.write(p + "\n")

print("Generated patterns saved to patterns.txt")
