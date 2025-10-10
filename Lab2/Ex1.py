# Find the percentage for all the dinucleotide and trinucleotide combinations for the sequence: S="ATTGTCCAATCTGTTG".

# 1. Build a brute force engine to generate all dinucleotide and trinucleotide combinations.
# 2. For each combination, find out the percentage inside the S sequence.
# 3. Show the percentage for each combination in the output of your implementation.


# Ex1_GUI.py
# Simple GUI to calculate dinucleotide and trinucleotide percentages
# in alphabetical order

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
from itertools import product
from collections import Counter

def all_kmers(k, alphabet="ACGT"):
    """Generate all k-length combinations of nucleotides."""
    return ["".join(p) for p in product(alphabet, repeat=k)]

def kmer_stats(seq: str, k: int, alphabet="ACGT"):
    """Return total windows, counts dict, and percentage dict."""
    seq = seq.upper()
    seq = "".join([s for s in seq if s in alphabet])
    total = max(len(seq) - k + 1, 0)
    found = Counter(seq[i:i+k] for i in range(len(seq) - k + 1))
    kmers = all_kmers(k, alphabet)
    counts = {kmer: found.get(kmer, 0) for kmer in kmers}
    perc = {kmer: (counts[kmer] / total * 100 if total else 0.0) for kmer in kmers}
    return total, counts, perc

class KmerApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Dinucleotide & Trinucleotide Percentages")
        self.geometry("700x600")
        self.configure(padx=10, pady=10)
        self.create_widgets()

    def create_widgets(self):
        ttk.Label(self, text="Enter DNA Sequence:", font=("Segoe UI", 11, "bold")).pack(anchor="w")

        self.seq_box = tk.Text(self, height=4, wrap="word")
        self.seq_box.pack(fill="x", pady=5)
        self.seq_box.insert("1.0", "ATTGTCCAATCTGTTG")

        btn_frame = ttk.Frame(self)
        btn_frame.pack(pady=5)
        ttk.Button(btn_frame, text="Compute k=2", command=lambda: self.compute(2)).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Compute k=3", command=lambda: self.compute(3)).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Compute Both", command=self.compute_both).pack(side="left", padx=5)

        ttk.Label(self, text="Results (alphabetical order):", font=("Segoe UI", 10, "bold")).pack(anchor="w", pady=(10, 0))

        # Text areas for results
        self.result2 = self.create_result_box("Dinucleotides (k=2)")
        self.result3 = self.create_result_box("Trinucleotides (k=3)")

    def create_result_box(self, title):
        frame = ttk.LabelFrame(self, text=title)
        frame.pack(fill="both", expand=True, pady=5)

        text_widget = tk.Text(frame, height=10, wrap="none")
        text_widget.pack(fill="both", expand=True, padx=5, pady=5)
        text_widget.config(state="disabled")

        return text_widget

    def compute(self, k):
        seq = self.seq_box.get("1.0", "end").strip()
        if not seq:
            messagebox.showerror("Error", "Please enter a DNA sequence.")
            return
        total, counts, perc = kmer_stats(seq, k)
        result_box = self.result2 if k == 2 else self.result3

        result_box.config(state="normal")
        result_box.delete("1.0", "end")
        result_box.insert("1.0", f"Total windows: {total}\n\n")
        result_box.insert("end", f"{'k-mer':<6}{'count':>8}{'percent':>12}\n")
        result_box.insert("end", "-" * 30 + "\n")

        for kmer in sorted(counts.keys()):
            line = f"{kmer:<6}{counts[kmer]:>8}{perc[kmer]:>11.2f}%\n"
            result_box.insert("end", line)

        result_box.config(state="disabled")

    def compute_both(self):
        self.compute(2)
        self.compute(3)

if __name__ == "__main__":
    app = KmerApp()
    app.mainloop()
