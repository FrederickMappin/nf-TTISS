"""
gRNAfinder — Resolve a guide index to its ID, orientation and sequence.

Because the Guides list alternates forward / reverse-complement entries
(index 0 = guide 1 fwd, index 1 = guide 1 RC, index 2 = guide 2 fwd, …),
this helper maps an index back to the original guide metadata.
"""

import numpy as np


def gRNAfinder(
    gRNAnumber: int,
    OriginalGuides: np.ndarray,
) -> list:
    """Map a linear guide index to the original guide entry.

    Parameters
    ----------
    gRNAnumber : int
        Index into the expanded (fwd + RC interleaved) guide list.
    OriginalGuides : np.ndarray
        Array of shape ``(N, 2)`` with columns ``[id, sequence]``.

    Returns
    -------
    list
        ``[GuideNumber: int, RC: bool, Sequence: str]``
        *RC* is True when the index corresponds to the reverse-complement
        entry of that guide.
    """
    RC = False
    if gRNAnumber % 2 == 1:
        # Odd index → reverse complement of the preceding forward entry
        GuideNumber = int(OriginalGuides[int((gRNAnumber - 1) / 2), 0])
        Sequence    = OriginalGuides[int((gRNAnumber - 1) / 2), 1]
        RC = True
    else:
        # Even index → forward entry
        GuideNumber = int(OriginalGuides[int(gRNAnumber / 2), 0])
        Sequence    = OriginalGuides[int(gRNAnumber / 2), 1]
    return [GuideNumber, RC, Sequence]
