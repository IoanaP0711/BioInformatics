#take an arbitrary DNA sequence from the NCBI, between 1000-3000 nucleotides/letters
#implement a software app that detects repetitions between 3b to 6b in this DNA sequence
#repetitive sequences refer to patters that repeat N times. Min no of repetitions is 2.


from Bio import SeqIO

# ---------------------------------------------------------
# Detect repetitive patterns (motifs) of length 3–6 bp
# that repeat at least 2 times in a DNA sequence
# ---------------------------------------------------------

def find_repetitions(sequence, min_len=3, max_len=6, min_reps=2):
    results = []

    n = len(sequence)
    for L in range(min_len, max_len + 1):              # pattern length (3 → 6)
        seen = {}
        for i in range(n - L + 1):
            pattern = sequence[i:i+L]
            if pattern not in seen:
                seen[pattern] = []
            seen[pattern].append(i)

        # Filter only motifs that repeat min_reps times
        for motif, positions in seen.items():
            if len(positions) >= min_reps:
                results.append((motif, positions))

    return results


# ---------------------------------------------------------
# Read a FASTA file
# ---------------------------------------------------------
def read_fasta(filepath):
    for record in SeqIO.parse(filepath, "fasta"):
        return str(record.seq).upper()
    return ""


# ---------------------------------------------------------
# Main execution
# ---------------------------------------------------------
if __name__ == "__main__":
    fasta_file = "dna.fasta"   # <-- replace with your downloaded NCBI file

    seq = read_fasta(fasta_file)
    print(f"Loaded sequence of length {len(seq)} bp\n")

    repetitions = find_repetitions(seq)

    if not repetitions:
        print("No repetitive patterns found.")
    else:
        print("Detected repetitive patterns:\n")
        for motif, positions in repetitions:
            print(f"Motif: {motif} (length {len(motif)})")
            print(f"Occurrences: {len(positions)}")
            print(f"Positions: {positions}\n")
