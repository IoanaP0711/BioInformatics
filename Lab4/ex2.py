from collections import Counter
from pathlib import Path
from Bio import SeqIO
from Bio.Data import CodonTable
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme()

def count_codons(fasta_file: str) -> Counter:
    seq = "".join(str(rec.seq).upper() for rec in SeqIO.parse(fasta_file, "fasta"))
    codons = [seq[i:i+3] for i in range(0, len(seq) - 2, 3) if "N" not in seq[i:i+3]]
    codons = [c for c in codons if len(c) == 3]
    return Counter(codons)

def top10_df(codons_counter: Counter) -> pd.DataFrame:
    return pd.DataFrame(codons_counter.most_common(10), columns=["Codon", "Frequency"])

def barplot_save(df: pd.DataFrame, title: str, outpath: str):
    plt.figure(figsize=(8, 5))
    sns.barplot(data=df, x="Codon", y="Frequency", color="C0", edgecolor="black")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()

def top_amino_acids_from_codon_counts(codon_counts: Counter, table_name="Standard", k=3):
    table = CodonTable.unambiguous_dna_by_name[table_name]
    aa_counts = Counter()
    for codon, count in codon_counts.items():
        if codon in table.forward_table:      # skip STOP codons
            aa_counts[table.forward_table[codon]] += count
    return aa_counts.most_common(k)

# --- run (A & B) ---
covid_codons = count_codons("covid.fasta")
flu_codons   = count_codons("influenza.fasta")

df_covid = top10_df(covid_codons)
df_flu   = top10_df(flu_codons)

Path("plots").mkdir(exist_ok=True)
barplot_save(df_covid, "Top 10 Most Frequent Codons – COVID-19",    "plots/covid_codons.png")
barplot_save(df_flu,   "Top 10 Most Frequent Codons – Influenza",  "plots/influenza_codons.png")

print("Saved charts to plots/covid_codons.png and plots/influenza_codons.png")

# --- (C) compare & show the most frequent codons between genomes ---
shared = set(covid_codons) & set(flu_codons)
combined = []
for c in shared:
    combined.append((c, covid_codons[c] + flu_codons[c], covid_codons[c], flu_codons[c]))
combined.sort(key=lambda t: t[1], reverse=True)
top_shared = combined[:10]

print("\nTop shared codons ranked by combined frequency (covid+flu):")
for codon, total, c_cnt, f_cnt in top_shared:
    print(f"{codon}: total={total} (COVID={c_cnt}, Influenza={f_cnt})")

# --- (D) top-3 amino acids for each genome ---
top3_covid = top_amino_acids_from_codon_counts(covid_codons)
top3_flu   = top_amino_acids_from_codon_counts(flu_codons)
print("\nTop 3 amino acids (COVID-19):", top3_covid)
print("Top 3 amino acids (Influenza):", top3_flu)
