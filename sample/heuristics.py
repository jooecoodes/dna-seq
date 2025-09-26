def choose_algorithm(length, gc_content, entropy):
    # Rule 1: Very short patterns → BP
    if length <= 64:
        return "BP"
    
    # Rule 2: Medium patterns
    if length <= 256:
        if entropy <= 1.1:
            return "KMP"
        else:
            return "BMH"
    
    # Rule 3: Longer patterns
    if length >= 512:
        if entropy >= 1.1:
            return "BMH"
        else:
            return "KMP"
    
    # Rule 4: Adjust for extreme GC content
    if gc_content <= 0.2 or gc_content >= 0.8:
        return "KMP"
    
    return "BMH"
