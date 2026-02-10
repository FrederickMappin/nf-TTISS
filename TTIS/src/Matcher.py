"""
Matcher — Best-match guide finder for a single genomic window.

Slides a 20 bp window across the target sequence, computes the Hamming
distance to every guide (forward + reverse complement), and returns the
best hit(s) along with orientation, PAM context, and cut-site distance.
"""

from Bio.Seq import Seq

from src.hamdist import hamdist
from src.gRNAfinder import gRNAfinder
from src.StringChanger import StringChanger


def Matcher(
    TargetString: str,
    CutSite: int,
    Guides: list[str],
    OriginalGuides,
    MismatchThreshold: int,
    CutDist: int,
    DelimiterReadout: bool,
    Delimiter: str,
    pam_length: int = 3,
) -> list:
    """Find the best-matching guide RNA for a genomic window.

    Parameters
    ----------
    TargetString : str
        The 100 bp genomic window sequence to search.
    CutSite : int
        Expected cut-site position for score calculation.
    Guides : list[str]
        Alternating forward / reverse-complement guide sequences.
    OriginalGuides : numpy.ndarray
        Original guide entries ``[[id, sequence], ...]``.
    MismatchThreshold : int
        Maximum Hamming distance to consider a guide a valid match.
    CutDist : int
        Distance (bp) from spacer start to predicted cut site.
    DelimiterReadout : bool
        If True, annotate matching positions with *Delimiter*.
    Delimiter : str
        Character used for matching positions (typically ``-``).
    pam_length : int
        Length of the PAM sequence (default 3 for NGG).

    Returns
    -------
    list
        ``[SingleMatch, MMs, GuideNumberList, CrudeSpacerList,
          RealSpacerList, CutSiteScoreList]``
    """
    BestDist = 25
    GuideNumber = 0
    CutSiteScore = 100
    SingleMatch = True
    MMs = []
    GuideNumberList = []
    CrudeSpacerList = []
    RealSpacerList = []
    CutSiteScoreList = []
    for g in range(len(Guides)):
        SpacerBestDist = 25
        SpacerBestMatch = 0
        for m in range(len(TargetString)-20):
            CurrentString = TargetString[m:m+20]
            SpacerTestDist = hamdist(CurrentString, Guides[g])
            if SpacerTestDist <= SpacerBestDist:
                SpacerBestDist = SpacerTestDist
                SpacerBestMatch = m
        t = SpacerBestMatch
        TestDist = SpacerBestDist
        if TestDist <= MismatchThreshold or TestDist <= BestDist:
            GuideNumber, RC, Sequence = gRNAfinder(g, OriginalGuides)
            if not RC:
                if TargetString[t:t+20+pam_length] != "":
                    AlignedString = TargetString[t:t+20+pam_length]
                    RealString = AlignedString
                else:
                    AlignedString = TargetString[t:t+20] + " PAM INDEX ERROR " + str(t)
                    RealString = TargetString[t:t+20]
                CutSiteScore = t + CutDist - CutSite
            else:
                if TargetString[t-pam_length:t+20] != "":
                    AlignedString = TargetString[t-pam_length:t+20]
                    SpacerSeq = Seq(TargetString[t-pam_length:t+20])
                    RealString = str(SpacerSeq.reverse_complement())
                else:
                    AlignedString = TargetString[t:t+20] + " PAM INDEX ERROR RC " + str(t)
                    SpacerSeq = Seq(TargetString[t:t+20])
                    RealString = str(SpacerSeq.reverse_complement())
                CutSiteScore = t + 20 - CutDist - CutSite
            RealString = StringChanger(RealString, Sequence, DelimiterReadout, Delimiter)
        if (TestDist <= MismatchThreshold and BestDist <= MismatchThreshold) or (TestDist == BestDist):
            SingleMatch = False
            MMs.append(TestDist)
            GuideNumberList.append(GuideNumber)
            CrudeSpacerList.append(AlignedString)
            RealSpacerList.append(RealString)
            CutSiteScoreList.append(CutSiteScore)
        elif (TestDist <= MismatchThreshold) or (TestDist < BestDist):
            BestDist = TestDist
            SingleMatch = True
            MMs = [TestDist]
            GuideNumberList = [GuideNumber]
            CrudeSpacerList = [AlignedString]
            RealSpacerList = [RealString]
            CutSiteScoreList = [CutSiteScore]
    return [SingleMatch, MMs, GuideNumberList, CrudeSpacerList, RealSpacerList, CutSiteScoreList]
