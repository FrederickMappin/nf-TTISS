import csv

def pam_matches(pam_seq, pam_pattern):
    """
    Returns True if pam_seq matches pam_pattern, where 'N' in pam_pattern is a wildcard.
    """
    if len(pam_seq) != len(pam_pattern):
        return False
    for a, b in zip(pam_seq, pam_pattern):
        if b == 'N':
            continue
        if a != b:
            return False
    return True

def write_matcher_output_to_csv(results, pam_sequence, output_file):
    """
    Write matcher results to a CSV, including PAM match status.
    Each row: [name, match, mismatches, guide_number, crude_spacer, real_spacer, cut_site_score, pam_match]
    pam_match: True if PAM in real_spacer matches pam_sequence (with 'N' as wildcard)
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
