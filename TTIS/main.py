import sys
import os
# Ensure src is in the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
"""
Usage:
    python3 main.py --guidefile <guidefile> --fastafile <fastafile> --pam <pam_sequence> [--output outputfile.csv]

Arguments:
    --guidefile    Path to the guide CSV file (e.g., PilotGuides.txt)
    --fastafile    Path to the FASTA file (e.g., chr8_reads_windows.fasta)
    --pam          PAM sequence (e.g., AGG, NGG, TTT, etc.)
    --output       (Optional) Output CSV filename (default: matcher_results.csv)

Example:
    python3 main.py --guidefile Guides.txt --fastafile windows.fasta --pam AGG --output matcher_results.csv
"""

import argparse
from Bio.Seq import Seq
import numpy as np
from src.Matcher import Matcher
from src.file_readers import read_guides, read_fasta
from src.write_matcher_output import write_matcher_output_to_csv


def main():
    parser = argparse.ArgumentParser(description="Run matcher and output filtered CSV based on PAM.")
    parser.add_argument("--guidefile", required=True, help="Guide file (CSV)")
    parser.add_argument("--fastafile", required=True, help="FASTA file")
    parser.add_argument("--pam", required=True, help="PAM sequence (e.g. NGG)")
    parser.add_argument("--output", default="matcher_results.csv", help="Output CSV filename")
    args = parser.parse_args()

    guide_file = args.guidefile
    fasta_file = args.fastafile
    pam_sequence = args.pam
    pam_length = len(pam_sequence)
    output_file = args.output

    guides, OriginalGuides = read_guides(guide_file)
    OriginalGuides = np.array(OriginalGuides)
    Guides = []
    for g in guides:
        Guides.append(g)
        Guides.append(str(Seq(g).reverse_complement()))

    fasta_seqs = read_fasta(fasta_file)

    # Parameters
    MismatchThreshold = 6
    CutDist = 17
    DelimiterReadout = True
    Delimiter = "-"
    StringSearchLength = 25

    # Run Matcher on each sequence and collect results
    results = {}
    for name, seq in fasta_seqs.items():
        result = Matcher(seq, CutDist, Guides, OriginalGuides, MismatchThreshold, CutDist, DelimiterReadout, Delimiter, pam_length=pam_length)
        print(f"{name}: {result}")
        results[name] = result

    # Write results to CSV
    write_matcher_output_to_csv(results, pam_sequence, output_file)
    print(f"Results written to {output_file}")

if __name__ == "__main__":
    main()
