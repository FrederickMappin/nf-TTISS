import numpy as np

def gRNAfinder(gRNAnumber, OriginalGuides):
    """
    Returns guide number, reverse complement status, and sequence
    """
    RC = False
    if gRNAnumber % 2 == 1:
        GuideNumber = int(OriginalGuides[int((gRNAnumber-1)/2), 0])
        Sequence = OriginalGuides[int((gRNAnumber-1)/2), 1]
        RC = True
    else:
        GuideNumber = int(OriginalGuides[int(gRNAnumber/2), 0])
        Sequence = OriginalGuides[int(gRNAnumber/2), 1]
    return [GuideNumber, RC, Sequence]
