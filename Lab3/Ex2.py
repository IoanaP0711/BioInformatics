#implement an app that uses a sliding window of 9 positions to scan a DNA sequence from a FASTA file.
# each sliding window provides a value for the components of a vector "P"
# plot vector P on a chart, where the horizontal axis is the length of the sequence and the vertical axis is the melting temp.  

import math
from collections import Counter
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter.ttk import Button, Label, Frame, Spinbox, Radiobutton

# Matplotlib for plotting inside Tkinter
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

# ------------------- Core Tm utilities -------------------
def clean_dna(seq: str) -> str:
    """Uppercase, remove spaces/newlines, convert U->T, keep only A/C/G/T."""
    s = "".join(seq.upper().split()).replace("U", "T")
    allowed = set("ACGT")
    s = "".join(ch for ch in s if ch in allowed)
    return s

def gc_percent(seq: str) -> float:
    c = Counter(seq)
    return 100.0 * (c["G"] + c["C"]) / len(seq) if seq else 0.0

def tm_wallace(seq: str) -> float:
    """Wallace rule: Tm = 2*(A+T) + 4*(G+C)"""
    c = Counter(seq)
    return 2.0 * (c["A"] + c["T"]) + 4.0 * (c["G"] + c["C"])

def tm_salt_adjusted(seq: str, na_mM: float = 50.0) -> float:
    """
    Tm = 81.5 + 16.6*log10([Na+]) + 0.41*(%GC) - 600/length
    [Na+] in mM (converted to M inside).
    """
    if not seq:
        return float("nan")
    na_M = max(na_mM, 0.0001) / 1000.0  # mM -> M; avoid log10(0)
    L = len(seq)
    return 81.5 + 16.6*math.log10(na_M) + 0.41*gc_percent(seq) - (600.0/L)

def read_fasta_sequence(path: Path) -> str:
    """
    Reads first sequence from a FASTA file (concatenates wrapped lines).
    Ignores header '>' lines and blank lines.
    """
    seq_parts = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_parts.append(line)
    return clean_dna("".join(seq_parts))

def sliding_windows(seq: str, k: int):
    """Yield (start_index, window_seq) for each k-length window."""
    for i in range(len(seq) - k + 1):
        yield i, seq[i:i+k]

# ------------------- GUI App -------------------
class SlidingTmApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Sliding-Window DNA Melting Temperature (Tm)")
        self.geometry("920x600")

        # State
        self.sequence = ""
        self.window_size = tk.IntVar(value=9)
        self.na_mM = tk.DoubleVar(value=50.0)
        self.formula = tk.StringVar(value="salt")  # 'salt' or 'wallace'

        # Top controls
        top = Frame(self)
        top.pack(fill="x", pady=6, padx=8)

        Button(top, text="Open FASTA…", command=self.load_fasta).grid(row=0, column=0, padx=4)
        Label(top, text="Window size:").grid(row=0, column=1, padx=(12,4))
        Spinbox(top, from_=3, to=200, textvariable=self.window_size, width=6).grid(row=0, column=2, padx=4)

        Label(top, text="[Na+] (mM):").grid(row=0, column=3, padx=(12,4))
        Spinbox(top, from_=0, to=1000, increment=5, textvariable=self.na_mM, width=8).grid(row=0, column=4, padx=4)

        Radiobutton(top, text="Salt-adjusted", variable=self.formula, value="salt").grid(row=0, column=5, padx=(16,4))
        Radiobutton(top, text="Wallace", variable=self.formula, value="wallace").grid(row=0, column=6, padx=4)

        Button(top, text="Compute & Plot", command=self.compute_and_plot).grid(row=0, column=7, padx=(16,4))

        # Info labels
        info = Frame(self)
        info.pack(fill="x", padx=8)
        self.lbl_file = Label(info, text="No file loaded")
        self.lbl_file.pack(anchor="w")
        self.lbl_stats = Label(info, text="")
        self.lbl_stats.pack(anchor="w", pady=(2,8))

        # Matplotlib Figure
        self.fig = Figure(figsize=(8.8, 4.6), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_xlabel("Position in sequence (start of window)")
        self.ax.set_ylabel("Melting temperature (°C)")
        self.ax.set_title("Vector P: Tm per sliding window")
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().pack(fill="both", expand=True, padx=8, pady=8)

    # ---------- Actions ----------
    def load_fasta(self):
        path = filedialog.askopenfilename(
            title="Choose FASTA file",
            filetypes=[("FASTA", "*.fasta *.fa *.fna *.ffn *.frn *.fas"), ("All files", "*.*")]
        )
        if not path:
            return
        try:
            seq = read_fasta_sequence(Path(path))
            if not seq:
                raise ValueError("No valid A/C/G/T bases found after cleaning.")
            self.sequence = seq
            self.lbl_file.config(text=f"Loaded: {path}")
            self.lbl_stats.config(text=f"Length: {len(seq)} bases | GC%: {gc_percent(seq):.2f}%")
            # Clear plot
            self.ax.clear()
            self.ax.set_xlabel("Position in sequence (start of window)")
            self.ax.set_ylabel("Melting temperature (°C)")
            self.ax.set_title("Vector P: Tm per sliding window")
            self.canvas.draw_idle()
        except Exception as e:
            messagebox.showerror("Error", f"Failed to read FASTA:\n{e}")

    def compute_and_plot(self):
        if not self.sequence:
            messagebox.showwarning("No sequence", "Please open a FASTA file first.")
            return
        k = int(self.window_size.get())
        if k < 1 or k > len(self.sequence):
            messagebox.showwarning("Window size", f"Window must be between 1 and {len(self.sequence)}.")
            return

        # Build vector P
        P = []
        X = []
        for i, win in sliding_windows(self.sequence, k):
            if self.formula.get() == "wallace":
                tm = tm_wallace(win)
            else:
                tm = tm_salt_adjusted(win, self.na_mM.get())
            P.append(tm)
            X.append(i + 1)  # positions start at 1 for nicer axis

        # Plot
        self.ax.clear()
        self.ax.plot(X, P, marker="", linewidth=1.5)
        self.ax.set_xlabel("Position in sequence (start of 9-nt window)" if k == 9
                           else f"Position in sequence (start of {k}-nt window)")
        self.ax.set_ylabel("Melting temperature (°C)")
        self.ax.set_title("Vector P: Tm per sliding window")
        self.ax.grid(True, linewidth=0.5, alpha=0.4)
        self.canvas.draw_idle()

        # Update stats
        if P:
            self.lbl_stats.config(
                text=f"Length: {len(self.sequence)} | Windows: {len(P)} | "
                     f"P(min/avg/max): {min(P):.2f} / {sum(P)/len(P):.2f} / {max(P):.2f} °C"
            )

if __name__ == "__main__":
    app = SlidingTmApp()
    app.mainloop()
