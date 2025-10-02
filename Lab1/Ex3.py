# use AI to design an app with a graphical user interface (GUI) which is able to integrate your alphabet from Ex1 and Ex2. Your app must have a button that allows the user to choose a Fasta file. 
# Fasta files contain a specific biological sequence format. The output should be known on the main window by using a Text box object on something similar.
# # fasta files have the format: 
# 1. the first line = inf line that shows the ID of the sequence and other type of info
# 2. the following lines = rows sequence which can be DNA, RNA or protein sequences(split in 80 characters until the end of the file)
# Use AI to simmulate a Fasta File for the input 

# fasta_gui_app.py
# GUI that integrates:
#  - Ex1: alphabet (unique symbols)
#  - Ex2: relative frequencies
# Supports loading a FASTA file or simulating one (DNA/RNA/Protein).

import os
import random
import tempfile
from collections import Counter
from textwrap import wrap
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

# ---------- Core logic (from Ex1 & Ex2) ----------

def find_alphabet(sequence: str):
    """Return the set of unique symbols in the sequence."""
    return set(sequence)

def calculate_relative_frequencies(sequence: str):
    """Return dict: symbol -> relative frequency (0..1)."""
    counts = Counter(sequence)
    total = len(sequence) if sequence else 1
    return {sym: counts[sym] / total for sym in counts}

# ---------- FASTA utilities ----------

def parse_fasta(path: str):
    """Parse the FIRST record in a FASTA file and return (header, sequence)."""
    header = None
    seq_parts = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is None:
                    header = line[1:]  # drop '>'
                else:
                    # stop at first record only
                    break
            else:
                seq_parts.append(line.upper())
    sequence = "".join(seq_parts).replace(" ", "")
    return header or "(no header found)", sequence

def format_report(header: str, sequence: str):
    """Build a human-readable report for the GUI Text widget."""
    alphabet = sorted(find_alphabet(sequence))
    freqs = calculate_relative_frequencies(sequence)
    # Sort by symbol for stable display
    freq_lines = [f"{sym}: {freqs[sym]:.4f}" for sym in sorted(freqs.keys())]

    lines = [
        f"Header: {header}",
        f"Length: {len(sequence)}",
        f"Alphabet ({len(alphabet)} symbols): {', '.join(alphabet)}",
        "",
        "Relative Frequencies:",
        *freq_lines,
    ]
    return "\n".join(lines)

# ---------- AI-like FASTA simulation ----------

DNA_ALPHABET = "ACGT"
RNA_ALPHABET = "ACGU"
PROTEIN_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"  # 20 canonical amino acids

def simulate_sequence(seq_type: str, length: int) -> str:
    if seq_type == "DNA":
        alphabet = DNA_ALPHABET
    elif seq_type == "RNA":
        alphabet = RNA_ALPHABET
    elif seq_type == "Protein":
        alphabet = PROTEIN_ALPHABET
    else:
        raise ValueError("Unknown seq_type")
    rng = random.Random()
    return "".join(rng.choice(alphabet) for _ in range(length))

def write_simulated_fasta(seq_type: str = "DNA", length: int = 500) -> str:
    """Create a temporary FASTA file with lines wrapped at 80 chars and return its path."""
    header = f"{seq_type}_SIM|id=SIM001|len={length}"
    seq = simulate_sequence(seq_type, length)
    wrapped = "\n".join(wrap(seq, 80))

    fd, path = tempfile.mkstemp(prefix=f"{seq_type}_sim_", suffix=".fasta", text=True)
    os.close(fd)
    with open(path, "w", encoding="utf-8") as f:
        f.write(f">{header}\n{wrapped}\n")
    return path

# ---------- GUI ----------

class FastaApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("FASTA Alphabet & Relative Frequency Analyzer")
        self.geometry("820x580")

        # Top controls frame
        top = ttk.Frame(self, padding=10)
        top.pack(side=tk.TOP, fill=tk.X)

        open_btn = ttk.Button(top, text="Open FASTA…", command=self.open_fasta)
        open_btn.pack(side=tk.LEFT, padx=(0, 10))

        ttk.Label(top, text="Simulate:").pack(side=tk.LEFT)
        self.sim_type = tk.StringVar(value="DNA")
        sim_menu = ttk.OptionMenu(top, self.sim_type, "DNA", "DNA", "RNA", "Protein")
        sim_menu.pack(side=tk.LEFT, padx=5)

        ttk.Label(top, text="Length:").pack(side=tk.LEFT, padx=(10, 0))
        self.sim_len = tk.StringVar(value="500")
        len_entry = ttk.Entry(top, textvariable=self.sim_len, width=8)
        len_entry.pack(side=tk.LEFT, padx=5)

        sim_btn = ttk.Button(top, text="Simulate FASTA", command=self.simulate_and_load)
        sim_btn.pack(side=tk.LEFT, padx=10)

        # Status label
        self.status = tk.StringVar(value="Ready.")
        status_label = ttk.Label(self, textvariable=self.status, anchor="w")
        status_label.pack(side=tk.BOTTOM, fill=tk.X, padx=10, pady=5)

        # Text area with scrollbar
        mid = ttk.Frame(self, padding=10)
        mid.pack(expand=True, fill=tk.BOTH)

        self.text = tk.Text(mid, wrap="word")
        self.text.configure(font=("Consolas", 11))
        yscroll = ttk.Scrollbar(mid, orient="vertical", command=self.text.yview)
        self.text.configure(yscrollcommand=yscroll.set)

        self.text.pack(side=tk.LEFT, expand=True, fill=tk.BOTH)
        yscroll.pack(side=tk.LEFT, fill=tk.Y)

        # Hint message
        self._set_text(
            "Open a FASTA file or click 'Simulate FASTA'.\n\n"
            "FASTA format reminder:\n"
            "  - First line starts with '>' and contains ID/metadata\n"
            "  - Following lines are the sequence (DNA/RNA/Protein),\n"
            "    typically wrapped to 80 characters per line.\n"
        )

    def _set_text(self, content: str):
        self.text.config(state=tk.NORMAL)
        self.text.delete("1.0", tk.END)
        self.text.insert(tk.END, content)
        self.text.config(state=tk.DISABLED)

    def open_fasta(self):
        path = filedialog.askopenfilename(
            title="Select FASTA file",
            filetypes=[("FASTA files", "*.fasta *.fa *.faa *.fna *.ffn *.frn"), ("All files", "*.*")]
        )
        if not path:
            return
        self.load_and_display(path)

    def simulate_and_load(self):
        try:
            length = int(self.sim_len.get())
            if length <= 0:
                raise ValueError
        except ValueError:
            messagebox.showerror("Invalid length", "Please enter a positive integer length.")
            return

        seq_type = self.sim_type.get()
        path = write_simulated_fasta(seq_type, length)
        self.load_and_display(path, simulated=True)

    def load_and_display(self, path: str, simulated: bool = False):
        try:
            header, seq = parse_fasta(path)
            if not seq:
                messagebox.showwarning("Empty sequence", "No sequence data found in the file.")
                return
            report = format_report(header, seq)
            self._set_text(report)
            tag = "Simulated" if simulated else "Loaded"
            self.status.set(f"{tag} FASTA: {os.path.basename(path)} (length {len(seq)})")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to process file:\n{e}")
            self.status.set("Error.")

if __name__ == "__main__":
    # Optional: seed for reproducibility of simulations
    random.seed()
    app = FastaApp()
    app.mainloop()
