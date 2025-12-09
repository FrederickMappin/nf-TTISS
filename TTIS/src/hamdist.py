def hamdist(str1, str2):
    """
    Returns Hamming distance between two strings
    """
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs += 1
    return diffs
