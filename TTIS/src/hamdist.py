"""
hamdist — Hamming distance between two equal-length strings.
"""


def hamdist(str1: str, str2: str) -> int:
    """Return the Hamming distance between *str1* and *str2*.

    The Hamming distance counts the number of positions at which the
    corresponding characters differ.  Only the overlapping portion
    (determined by :pyfunc:`zip`) is compared, so if the strings have
    different lengths the trailing characters of the longer string are
    silently ignored.

    Parameters
    ----------
    str1, str2 : str
        Two nucleotide (or arbitrary character) strings to compare.

    Returns
    -------
    int
        Number of mismatched positions.
    """
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs += 1
    return diffs
