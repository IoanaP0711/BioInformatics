#download 10 influenza genoms. For each genome plot on a chart the most frequent repetitions.

import os
import requests
import matplotlib.pyplot as plt

# ----------------------------------------
# Download 10 influenza genomes from NCBI
# ----------------------------------------

influenza_urls = [
    "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NC_007366.1&db=nuccore&report=fasta",
    "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NC_007367.1&db=nuccore&report=fasta",
    "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NC_007368.1&db=nuccore&report=fasta",
    "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NC_007369.1&db=nuccore&report=fasta",
    "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NC_007370.1&db=nuccore&report=fasta",
    "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NC_007371.1&db=nuccore&report=fasta",
    "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NC_007372.1&db=nuccore&report=fasta",
    "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NC_007373.1&db=nuccore&report=fasta",
    "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NC_007374.1&db=nuccore&report=fasta",
    "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NC_007375.1&db=nuccore&report=fasta"
]

os.makedirs("genomes", exist_ok=True)

def download_fasta(url, index):
    response = requests.get(url)
    file_path = f"genomes/influenza_{index}.fasta"
    with open(file_path, "w") as f:
        f.write(response.text)
    return file_path

# ----------------------------------------
# Read FASTA (pure Python)
# ----------------------------------------
def read_fasta(file_path):
    seq = ""
    with open(file_path) as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip()
    return seq

# ----------------------------------------
# Find repeated motifs (length 3–6)
# ----------------------------------------
def find_repeats(sequence, min_len=3, max_len=6):
    counts = {}
    n = len(sequence)

    for L in range(min_len, max_len + 1):
        for i in range(n - L + 1):
            motif = sequence[i:i+L]
            counts[motif] = counts.get(motif, 0) + 1

    # keep only motifs repeated >= 2 times
    repeats = {m: c for m, c in counts.items() if c >= 2}

    if not repeats:
        return None, 0

    # return motif with maximum count
    best = max(repeats, key=repeats.get)
    return best, repeats[best]

# ----------------------------------------
# MAIN
# ----------------------------------------
most_frequent_counts = []
labels = []

print("Downloading genomes...")
for i, url in enumerate(influenza_urls):
    path = download_fasta(url, i+1)
    seq = read_fasta(path)

    motif, freq = find_repeats(seq)
    if motif:
        print(f"Genome {i+1}: Most frequent motif = {motif}, count = {freq}")
    else:
        print(f"Genome {i+1}: No repetitions found")

    labels.append(f"Genome {i+1}")
    most_frequent_counts.append(freq)

# ----------------------------------------
# Plot results
# ----------------------------------------
plt.figure(figsize=(12, 6))
plt.bar(labels, most_frequent_counts)
plt.xticks(rotation=45)
plt.xlabel("Influenza Genomes")
plt.ylabel("Most Frequent Motif Count")
plt.title("Most Frequent Repetitions (Motifs 3–6 bp) in 10 Influenza Genomes")
plt.tight_layout()
plt.savefig("influenza_repetitions.png")
plt.show()
