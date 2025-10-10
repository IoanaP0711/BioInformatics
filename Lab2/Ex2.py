# find in sequence S only the dinucleodites and trinucleodites that exist, without the use of the brute force engine. In order to achive the result one must verify each combination starting from the beggining of sequence S Example: ABAA we have AB, BA, AA, ABA

from typing import List, Tuple

def find_di_tri(seq: str, alphabetical: bool = False) -> Tuple[List[str], List[str]]:
    """
    Return the unique dinucleotides and trinucleotides found in `seq`,
    discovered by a single left-to-right scan (no brute-force generation).

    Parameters
    ----------
    seq : str
        Input sequence (e.g., DNA string like 'ATTGCCCCGAAT' or generic 'ABAA').
    alphabetical : bool
        If True, results are returned sorted alphabetically; otherwise, in order of first appearance.

    Returns
    -------
    (dinucs, trinucs) : (List[str], List[str])
        Lists of unique 2-mers and 3-mers found in the sequence.
    """
    s = seq.strip().upper()

    dinucs, trinucs = [], []
    seen2, seen3 = set(), set()

    # collect dinucleotides
    for i in range(len(s) - 1):
        d = s[i:i+2]
        if d not in seen2:
            seen2.add(d)
            dinucs.append(d)

    # collect trinucleotides
    for i in range(len(s) - 2):
        t = s[i:i+3]
        if t not in seen3:
            seen3.add(t)
            trinucs.append(t)

    if alphabetical:
        dinucs = sorted(dinucs)
        trinucs = sorted(trinucs)

    return dinucs, trinucs


# --- Examples ---
# Example from your prompt:
seq1 = "ABAA"
d2, t3 = find_di_tri(seq1)
print(seq1, "-> dinucs:", d2, "trinucs:", t3)
# ABAA -> dinucs: ['AB', 'BA', 'AA'] trinucs: ['ABA']

# A DNA example:
seq2 = "ATTGCCCCGAAT"
d2, t3 = find_di_tri(seq2)
print(seq2, "-> dinucs:", d2, "trinucs:", t3)

# Same, but alphabetical order if you ever need it:
d2a, t3a = find_di_tri(seq2, alphabetical=True)
print("Alphabetical:", "dinucs:", d2a, "trinucs:", t3a)
