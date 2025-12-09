def StringChanger(RealSpacer, guideString, DelimiterReadout, Delimiter):
    realseq = RealSpacer
    for i in range(20):
        if DelimiterReadout:
            if realseq[i] == guideString[i]:
                newstring = realseq[:i] + Delimiter + realseq[i+1:]
                realseq = newstring
        else:
            if realseq[i] != guideString[i]:
                newstring = realseq[:i] + realseq[i].lower() + realseq[i+1:]
                realseq = newstring
    return realseq
