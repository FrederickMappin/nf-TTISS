"""
write_matcher_output — Serialize Matcher results to CSV with PAM filtering.

Only windows whose PAM-adjacent bases match the user-supplied pattern are
written to the output.  ``N`` in the PAM pattern acts as a wildcard
matching any nucleotide.
"""

from __future__ import annotations

import csv


def pam_matches(pam_seq: str, pam_pattern: str) -> bool:
    """Check whether *pam_seq* matches *pam_pattern*.

    ``N`` in *pam_pattern* is treated as a wildcard that matches any
    character.  All other characters must match exactly.

    Parameters
    ----------
    pam_seq : str
        Observed PAM nucleotides extracted from the spacer.
    pam_pattern : str
        Expected PAM pattern (e.g. ``NGG``).

    Returns
    -------
    bool
    """
    if len(pam_seq) != len(pam_pattern):
        return False
    for a, b in zip(pam_seq, pam_pattern):
        if b == 'N':
            continue
        if a != b:
            return False
    return True


def write_matcher_output_to_csv(
    results: dict,
    pam_sequence: str,
    output_file: str,
) -> None:
    """Write PAM-validated Matcher results to a CSV file.

    For each window result the PAM bases at the end of the RealSpacer are
    compared to *pam_sequence*.  Only passing rows are emitted.

    Columns
    -------
    Name, SingleMatch, Mismatches, GuideNumber, CrudeSpacer,
    RealSpacer, CutSiteScore, PAM_Match
    """
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Name", "SingleMatch", "Mismatches", "GuideNumber", "CrudeSpacer", "RealSpacer", "CutSiteScore", "PAM_Match"])
        for name, result in results.items():
            # Extract PAM from real_spacer (last len(pam_sequence) bases)
            real_spacer = result[4][0] if isinstance(result[4], list) else result[4]
            pam_in_spacer = real_spacer[-len(pam_sequence):]
            pam_match = pam_matches(pam_in_spacer, pam_sequence)
            if pam_match:
                writer.writerow([
                    name,
                    result[0],
                    result[1][0] if isinstance(result[1], list) else result[1],
                    result[2][0] if isinstance(result[2], list) else result[2],
                    result[3][0] if isinstance(result[3], list) else result[3],
                    real_spacer,
                    result[5][0] if isinstance(result[5], list) else result[5],
                    pam_match
                ])
