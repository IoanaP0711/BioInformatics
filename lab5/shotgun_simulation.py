#!/usr/bin/env python3
"""
A simple shotgun sequencing simulation.

Steps:
1. Read a DNA sequence from a FASTA file (1000–3000 bp recommended).
2. Take 2000 random reads of length 100–150 bp from this sequence.
3. Store the reads in a list and in a text file.
4. Try to rebuild the original sequence using a greedy overlap assembly.
5. Print some stats and discuss the main problem of this approach.
"""

import random
from typing import List, Tuple


# ---------- FASTA READER ----------

def read_fasta(path: str) -> str:
    """
    Read the first sequence from a FASTA file and return it as a string (upper-case).
    Lines starting with '>' are headers and are skipped.
    """
    seq_parts = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_parts.append(line)
    sequence = "".join(seq_parts).upper()
    return sequence


# ---------- READ (FRAGMENT) SAMPLING ----------

def get_random_reads(
    sequence: str,
    num_reads: int = 2000,
    min_len: int = 100,
    max_len: int = 150,
) -> List[str]:
    """
    Take num_reads random fragments from the sequence,
    each with a random length between min_len and max_len.
    """
    print(f"Sampling {num_reads} reads from sequence of length {len(sequence)} ...")
    reads = []
    seq_len = len(sequence)

    for i in range(num_reads):
        read_len = random.randint(min_len, max_len)
        start = random.randint(0, seq_len - read_len)
        frag = sequence[start:start + read_len]
        reads.append(frag)

        # show a preview of the first few reads
        if i < 5:
            print(f"Read {i+1}: start={start}, length={read_len}, seq={frag[:40]}...")

    print(f"Total reads sampled: {len(reads)}")
    return reads


# ---------- GREEDY ASSEMBLY ----------

def overlap(a: str, b: str, min_overlap: int) -> int:
    """
    Return the length of the longest suffix of 'a' that matches
    a prefix of 'b', with at least min_overlap bases.
    """
    max_possible = min(len(a), len(b))
    # check from largest possible overlap down to min_overlap
    for ov in range(max_possible, min_overlap - 1, -1):
        if a[-ov:] == b[:ov]:
            return ov
    return 0


def greedy_assembly(reads: List[str], min_overlap: int = 30) -> List[str]:
    """
    Very naive greedy shortest-superstring style assembly.

    Repeatedly find the pair of reads with the maximum overlap and merge them,
    until we can no longer merge with at least min_overlap bases.
    Returns the list of remaining contigs (ideally one big contig).
    """
    reads = reads[:]  # work on a copy

    print("\nStarting greedy assembly...")
    iteration = 0
    while len(reads) > 1:
        best_i, best_j = None, None
        best_overlap = 0
        best_merged = None

        # find the best pair to merge
        for i in range(len(reads)):
            for j in range(len(reads)):
                if i == j:
                    continue
                ov = overlap(reads[i], reads[j], min_overlap)
                if ov > best_overlap:
                    best_overlap = ov
                    best_i, best_j = i, j
                    best_merged = reads[i] + reads[j][ov:]

        if best_overlap < min_overlap or best_merged is None:
            # no more useful merges
            break

        iteration += 1
        print(
            f"Iteration {iteration}: merging reads {best_i} and {best_j} "
            f"(overlap={best_overlap}), new contig length={len(best_merged)}"
        )

        # build the new list of reads
        new_reads = []
        for k, r in enumerate(reads):
            if k not in (best_i, best_j):
                new_reads.append(r)
        new_reads.append(best_merged)
        reads = new_reads

    print(f"Assembly finished with {len(reads)} contig(s).")
    return reads


# ---------- MAIN SCRIPT ----------

if __name__ == "__main__":
    # 1. Load DNA sequence from FASTA
    fasta_path = "lab5.fasta"   
    dna_sequence = read_fasta(fasta_path)

    print("\n--- Shotgun Simulation ---")
    print(f"Sequence loaded from {fasta_path}")
    print(f"Original length: {len(dna_sequence)} bases")
    print(f"First 60 bases: {dna_sequence[:60]}")
    print(f"Last  60 bases: {dna_sequence[-60:]}\n")

    # (optional) If the sequence is longer than 3000, take a slice
    if len(dna_sequence) > 3000:
        dna_sequence = dna_sequence[:3000]
        print(f"Using only the first 3000 bases for the simulation.")
        print(f"New length: {len(dna_sequence)}\n")

    # 2–3. Generate random reads and save them
    reads = get_random_reads(
        dna_sequence,
        num_reads=2000,
        min_len=100,
        max_len=150,
    )

    with open("random_reads.txt", "w") as out:
        for r in reads:
            out.write(r + "\n")
    print("\nReads written to 'random_reads.txt'.")

    # 4. Rebuild sequence via greedy assembly
    contigs = greedy_assembly(reads, min_overlap=30)
    # choose the longest contig as our assembled genome
    assembled = max(contigs, key=len)

    with open("assembled_sequence.txt", "w") as out:
        out.write(assembled + "\n")
    print(f"\nAssembled sequence saved to 'assembled_sequence.txt'.")
    print(f"Assembled length: {len(assembled)} bases")

    # simple comparison with the original (only for the overlapping part)
    compare_len = min(len(dna_sequence), len(assembled))
    matches = sum(
        1 for a, b in zip(dna_sequence[:compare_len], assembled[:compare_len]) if a == b
    )
    identity = matches / compare_len * 100
    print(f"Identity with original over first {compare_len} bases: "
          f"{identity:.2f}%")

    # 5. Print the conceptual main problem of this algorithm
    print("\nMain problem of this naive shotgun assembly:")
    print(
        "- It becomes confused in repetitive regions and where coverage is low.\n"
        "  Different parts of the genome can look the same, so the greedy\n"
        "  algorithm may join reads in the wrong order or fail to join them at all."
    )
