"""
StringChanger — Annotate spacer mismatches relative to the guide.

Two annotation modes:

* **Delimiter mode** (``DelimiterReadout=True``):
  Matching positions are replaced with *Delimiter* (e.g. ``-``),
  so only mismatches remain as uppercase letters.

* **Lowercase mode** (``DelimiterReadout=False``):
  Mismatched positions are lowercased; matches keep their case.
"""


def StringChanger(
    RealSpacer: str,
    guideString: str,
    DelimiterReadout: bool,
    Delimiter: str,
) -> str:
    """Annotate a 20 nt spacer string against its matched guide.

    Parameters
    ----------
    RealSpacer : str
        Extracted spacer sequence (with PAM context possibly appended).
    guideString : str
        The 20 nt reference guide sequence.
    DelimiterReadout : bool
        If True, replace **matching** positions with *Delimiter*.
        If False, **lowercase** mismatched positions.
    Delimiter : str
        Character used for matching positions (typically ``-``).

    Returns
    -------
    str
        Annotated spacer string of at least 20 characters.
    """
    realseq = RealSpacer
    for i in range(20):
        if DelimiterReadout:
            if realseq[i] == guideString[i]:
                realseq = realseq[:i] + Delimiter + realseq[i + 1:]
        else:
            if realseq[i] != guideString[i]:
                realseq = realseq[:i] + realseq[i].lower() + realseq[i + 1:]
    return realseq
