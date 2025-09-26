import math

def compute_features(sequence):
    length = len(sequence)
    gc_content = (sequence.count("G") + sequence.count("C")) / length
    
    # Shannon entropy for repetitiveness
    freqs = [sequence.count(base)/length for base in "ACGT" if sequence.count(base) > 0]
    entropy = -sum(f * math.log2(f) for f in freqs)
    
    return [length, gc_content, entropy]

seq = "ATGCGCGCTTATCGATCGATC"  # example DNA snippet
features = compute_features(seq)

print("Features:", features)

