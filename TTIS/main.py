#!/usr/bin/env python3
"""
nf-TTISS Guide Matcher — Entry Point
=====================================

Match CRISPR guide-RNA sequences against genomic windows extracted by the
Nextflow pipeline.  For each window the tool finds the best-matching guide
(forward or reverse-complement), validates the PAM, and writes hits to a
CSV file.

Usage:
    python3 main.py --guidefile <guidefile> --fastafile <fastafile> \
                    --pam <pam_sequence> [--output outputfile.csv]

Arguments:
    --guidefile   Path to the guide CSV file.  Each line: ``<id>,<20nt seq>``
    --fastafile   Path to the FASTA file of 100 bp genomic windows.
    --pam         PAM sequence (e.g. NGG, AGG, TTT). 'N' is a wildcard.
    --output      (Optional) Output CSV filename (default: matcher_results.csv)

Example:
    python3 main.py --guidefile Guides.txt \
                    --fastafile windows.fasta \
                    --pam NGG \
                    --output matcher_results.csv
"""

import sys
import os

# ---------------------------------------------------------------------------
# Ensure the parent directory is on the Python path so that ``src`` can be
# imported both when running from the project root and from within TTIS/.
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import argparse

from Bio.Seq import Seq
import numpy as np

from src.Matcher import Matcher
from src.file_readers import read_guides, read_fasta
from src.write_matcher_output import write_matcher_output_to_csv

# ---------------------------------------------------------------------------
# Matching parameters (hard-coded; adjust here if needed)
# ---------------------------------------------------------------------------
MISMATCH_THRESHOLD = 6       # Max Hamming distance to consider a match
CUT_DIST           = 17      # Expected distance (bp) from spacer start to cut site
DELIMITER_READOUT  = True    # Annotate matches with delimiter characters
DELIMITER          = "-"     # Character for matching positions in RealSpacer


def main() -> None:
    """Parse CLI arguments, run guide matching, and write results."""

    parser = argparse.ArgumentParser(
        description="Match CRISPR guide-RNA sequences to genomic windows and validate PAM.",
    )
    parser.add_argument("--guidefile", required=True, help="Guide file (CSV): id,sequence")
    parser.add_argument("--fastafile", required=True, help="FASTA file of genomic windows")
    parser.add_argument("--pam",       required=True, help="PAM pattern (e.g. NGG). N = any nucleotide.")
    parser.add_argument("--output",    default="matcher_results.csv", help="Output CSV filename")
    args = parser.parse_args()

    pam_length = len(args.pam)

    # ------------------------------------------------------------------
    # 1. Read guide sequences and build forward + reverse-complement list
    # ------------------------------------------------------------------
    guides_fwd, original_guides = read_guides(args.guidefile)
    original_guides = np.array(original_guides)

    guides: list[str] = []
    for g in guides_fwd:
        guides.append(g)                                  # forward
        guides.append(str(Seq(g).reverse_complement()))   # reverse complement

    # ------------------------------------------------------------------
    # 2. Read FASTA windows
    # ------------------------------------------------------------------
    fasta_seqs = read_fasta(args.fastafile)

    # ------------------------------------------------------------------
    # 3. Run Matcher on every window and collect results
    # ------------------------------------------------------------------
    results: dict = {}
    for name, seq in fasta_seqs.items():
        result = Matcher(
            TargetString=seq,
            CutSite=CUT_DIST,
            Guides=guides,
            OriginalGuides=original_guides,
            MismatchThreshold=MISMATCH_THRESHOLD,
            CutDist=CUT_DIST,
            DelimiterReadout=DELIMITER_READOUT,
            Delimiter=DELIMITER,
            pam_length=pam_length,
        )
        print(f"{name}: {result}")
        results[name] = result

    # ------------------------------------------------------------------
    # 4. Write only PAM-validated hits to CSV
    # ------------------------------------------------------------------
    write_matcher_output_to_csv(results, args.pam, args.output)
    print(f"Results written to {args.output}")


if __name__ == "__main__":
    main()
