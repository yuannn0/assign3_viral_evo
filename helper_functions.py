import numpy as np
import math
def jukes_cantor(reference_sequence: str, distant_sequence: str) -> float:
    """The Jukes-Cantor correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Jukes-Cantor corrected genetic distance using Hamming distance.
        For example 1.163.

    """
    if len(reference_sequence) != len(distant_sequence):
        raise ValueError("Sequences must have the same length")
    length = 0
    difference = 0
    for seq1, seq2 in zip(reference_sequence.upper(), distant_sequence.upper()):
        if seq1 not in "ACGT" or seq2 not in "ACGT":
            continue
        if seq1 != seq2:
            difference += 1
        length += 1
    if length == 0:
        raise ValueError ("The length of sequence cannot be 0")
    
    #calculate p
    p = difference / length
    if p >= 0.75:
        raise ValueError ("p value is greater than 0.75")
    
    #calculate d
    d_value = -0.75 * math.log(1 - (4 / 3) * p) * length

    return float(d_value)


def kimura_two_parameter(reference_sequence: str, distant_sequence: str) -> float:
    """The Kimura Two Parameter correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Kimura corrected genetic distance using Hamming distance.
        For example 1.196.

    """
    if len(reference_sequence) != len(distant_sequence):
        raise ValueError ("Sequences must have the same length")
    pairs = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

    length = 0
    transitions = 0
    transversions = 0

    for seq1, seq2 in zip(reference_sequence.upper(), distant_sequence.upper()):
        # # ignore gap "-"
        # if seq1 == "-" or seq2 == "-":
        #     continue

        if seq1 not in "ACGT" or seq2 not in "ACGT":
            continue
            # raise ValueError ("nucleotides error")
        length += 1

        if seq1 != seq2:
            if (seq1, seq2) in pairs:
                transitions += 1
            else:
                transversions += 1

        if length == 0:
            raise ValueError("No valid position to compare")

        # calculate P and Q
        P = transitions / length
        Q = transversions / length

        term1 = 1 - 2 * P - Q
        term2 = 1 - 2 * Q

        if term1 <= 0 or term2 <= 0:
            raise ValueError("K2P is undefined")

        distance_per_nt = (-0.5 * math.log(term1) - 0.25 * math.log(term2))

        return float (distance_per_nt * length)
    
