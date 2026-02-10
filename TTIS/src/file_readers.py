"""
file_readers — I/O helpers for guide and FASTA files.

Functions
---------
read_guides(guide_file)
    Parse a CSV guide file (``id,sequence``) and return two lists.
read_fasta(fasta_file)
    Parse a simple FASTA file and return an ordered dict of name→sequence.
"""

from __future__ import annotations


def read_guides(guide_file: str) -> tuple[list[str], list[list]]:
    """Read a guide-RNA CSV file.

    Each line is expected to be ``<numeric_id>,<20nt_sequence>``.

    Parameters
    ----------
    guide_file : str
        Path to the guide CSV file.

    Returns
    -------
    guides : list[str]
        List of 20 nt guide sequences.
    original_guides : list[list]
        List of ``[id: int, sequence: str]`` pairs (used by gRNAfinder).
    """
    guides: list[str] = []
    original_guides: list[list] = []
    with open(guide_file) as f:
        for line in f:
            parts = line.strip().split(",")
            if len(parts) == 2:
                num, seq = parts
                guides.append(seq)
                original_guides.append([int(num), seq])
    return guides, original_guides


def read_fasta(fasta_file: str) -> dict[str, str]:
    """Parse a FASTA file into a name → sequence dictionary.

    Multi-line sequences are concatenated.  Headers are stripped of the
    leading ``>`` character.

    Parameters
    ----------
    fasta_file : str
        Path to the FASTA file.

    Returns
    -------
    dict[str, str]
        Ordered mapping of sequence name to nucleotide string.
    """
    seqs: dict[str, str] = {}
    name: str | None = None
    seq = ""
    with open(fasta_file) as f:
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
