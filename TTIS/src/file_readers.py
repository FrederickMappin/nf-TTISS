def read_guides(guide_file):
    guides = []
    OriginalGuides = []
    with open(guide_file) as f:
        for line in f:
            parts = line.strip().split(",")
            if len(parts) == 2:
                num, seq = parts
                guides.append(seq)
                OriginalGuides.append([int(num), seq])
    return guides, OriginalGuides

def read_fasta(fasta_file):
    with open(fasta_file) as f:
        seqs = {}
        name = None
        seq = ""
        for line in f:
            if line.startswith(">"):
                if name:
                    seqs[name] = seq
                name = line.strip()[1:]
                seq = ""
            else:
                seq += line.strip()
        if name:
            seqs[name] = seq
    return seqs
