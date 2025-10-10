#design an app by using AI, which contains a GUI that allows the user to select a FASTA file.
# The content of the FASTA file should be analyzed by using a sliding window of 30 positions.
# The content for each sliding window should be used in order to compute the relative frequencies of the nucleotides in the alphabet of the sequence.
# The output of the app should be a chart containing 4 signals, one for each nucleotide in the alphabet of the sequence.
#


import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import ttk
import os
import math

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure


# ----------------------------
# FASTA utilities
# ----------------------------
def read_fasta_sequence(path: str) -> str:
    """
    Reads a FASTA file and returns the concatenated sequence (uppercase),
    ignoring header lines ('>') and whitespace.
    """
    seq_parts = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line:
                continue
            if line.startswith(">"):
                continue
            seq_parts.append("".join(line.strip().split()).upper())
    return "".join(seq_parts)


def choose_alphabet_and_normalize(seq: str):
    """
    Decide the 4-letter alphabet and normalize the sequence accordingly.
    If U present and T not present -> alphabet A,C,G,U (keep U).
    Else -> alphabet A,C,G,T (map U->T).
    Returns: alphabet (list of 4 chars), normalized_seq (str)
    """
    has_u = "U" in seq
    has_t = "T" in seq
    if has_u and not has_t:
        alphabet = ["A", "C", "G", "U"]
        normalized = seq  # keep U as U
    else:
        alphabet = ["A", "C", "G", "T"]
        normalized = seq.replace("U", "T")  # map U to T if any
    return alphabet, normalized


# ----------------------------
# Sliding-window frequencies using prefix sums
# ----------------------------
def window_relative_frequencies(seq: str, k: int, alphabet):
    """
    Computes relative frequency time-series per symbol in `alphabet`
    using a sliding window of size k.

    Ambiguous letters (not in alphabet) are ignored in a window's denominator.
    If a window contains 0 valid letters, frequency is NaN for that window.

    Returns:
      positions: list[int] - window start indices (1-based)
      freqs: dict[symbol] -> list[float] of length num_windows
    """
    n = len(seq)
    if k <= 0:
        raise ValueError("Window size must be positive.")
    if k > n:
        return [], {sym: [] for sym in alphabet}

    # Build prefix sums for each symbol and for "valid" positions
    # prefix_counts[s][i] = number of s in seq[:i]
    prefix_counts = {s: [0] * (n + 1) for s in alphabet}
    prefix_valid = [0] * (n + 1)  # counts letters that are in alphabet

    for i, ch in enumerate(seq, start=1):
        for s in alphabet:
            prefix_counts[s][i] = prefix_counts[s][i - 1]
        prefix_valid[i] = prefix_valid[i - 1]

        if ch in prefix_counts:
            prefix_counts[ch][i] += 1
            prefix_valid[i] += 1
        # else: ambiguous, do not increment

    num_windows = n - k + 1
    positions = list(range(1, num_windows + 1))  # 1-based start index

    freqs = {s: [math.nan] * num_windows for s in alphabet}

    for idx, start in enumerate(positions):
        end = start + k - 1
        valid = prefix_valid[end] - prefix_valid[start - 1]
        if valid == 0:
            # all ambiguous in this window
            for s in alphabet:
                freqs[s][idx] = math.nan
            continue
        for s in alphabet:
            count_s = prefix_counts[s][end] - prefix_counts[s][start - 1]
            freqs[s][idx] = count_s / valid

    return positions, freqs


# ----------------------------
# GUI Application
# ----------------------------
class FastaWindowApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("FASTA Sliding-Window Nucleotide Frequencies")
        self.geometry("1000x700")

        self.seq = ""
        self.file_path = None
        self.alphabet = ["A", "C", "G", "T"]
        self.positions = []
        self.freqs = {s: [] for s in self.alphabet}

        self._build_ui()

    def _build_ui(self):
        # Top controls frame
        top = ttk.Frame(self, padding=10)
        top.pack(side=tk.TOP, fill=tk.X)

        # File picker
        self.file_label_var = tk.StringVar(value="No file selected")
        ttk.Button(top, text="Open FASTA...", command=self.on_open).pack(side=tk.LEFT)
        ttk.Label(top, textvariable=self.file_label_var).pack(side=tk.LEFT, padx=10)

        # Window size
        ttk.Label(top, text="Window size:").pack(side=tk.LEFT, padx=(20, 5))
        self.win_entry = ttk.Entry(top, width=6)
        self.win_entry.insert(0, "30")
        self.win_entry.pack(side=tk.LEFT)

        ttk.Button(top, text="Analyze & Plot", command=self.on_analyze).pack(side=tk.LEFT, padx=10)

        # Info line
        self.info_var = tk.StringVar(value="")
        info_frame = ttk.Frame(self, padding=(10, 0))
        info_frame.pack(side=tk.TOP, fill=tk.X)
        ttk.Label(info_frame, textvariable=self.info_var).pack(side=tk.LEFT)

        # Matplotlib Figure
        self.fig = Figure(figsize=(9.5, 5.8), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_xlabel("Window start (1-based index)")
        self.ax.set_ylabel("Relative frequency")
        self.ax.set_title("Sliding-window nucleotide frequencies")
        self.ax.grid(True)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Bottom status bar
        self.status_var = tk.StringVar(value="Pick a FASTA file to begin.")
        status = ttk.Label(self, textvariable=self.status_var, anchor="w", padding=(10, 5))
        status.pack(side=tk.BOTTOM, fill=tk.X)

    def on_open(self):
        path = filedialog.askopenfilename(
            title="Select FASTA file",
            filetypes=[("FASTA files", "*.fasta *.fa *.fna *.ffn *.faa *.frn"),
                       ("All files", "*.*")]
        )
        if not path:
            return
        try:
            seq = read_fasta_sequence(path)
            if not seq:
                messagebox.showwarning("Empty", "The FASTA appears to contain no sequence.")
                return
            self.file_path = path
            self.seq = seq
            self.alphabet, self.seq = choose_alphabet_and_normalize(self.seq)

            base = os.path.basename(path)
            self.file_label_var.set(base)
            self.status_var.set(f"Loaded {base}: length={len(self.seq)} bases. Alphabet={''.join(self.alphabet)}")
            self.info_var.set("")
            self.clear_plot()
        except Exception as e:
            messagebox.showerror("Error reading FASTA", str(e))

    def on_analyze(self):
        if not self.seq:
            messagebox.showinfo("No sequence", "Please select a FASTA file first.")
            return
        # get window size
        try:
            k = int(self.win_entry.get())
            if k <= 0:
                raise ValueError
        except Exception:
            messagebox.showerror("Invalid window", "Window size must be a positive integer.")
            return

        positions, freqs = window_relative_frequencies(self.seq, k, self.alphabet)

        if not positions:
            messagebox.showinfo("Too short", f"Sequence length ({len(self.seq)}) is shorter than window size ({k}).")
            return

        self.positions = positions
        self.freqs = freqs

        # draw
        self.redraw_plot()

        # info text
        windows = len(positions)
        self.info_var.set(
            f"Alphabet={''.join(self.alphabet)} | Sequence length={len(self.seq)} | Window={k} | Windows={windows}"
        )
        self.status_var.set("Done.")

    def clear_plot(self):
        self.ax.cla()
        self.ax.set_xlabel("Window start (1-based index)")
        self.ax.set_ylabel("Relative frequency")
        self.ax.set_title("Sliding-window nucleotide frequencies")
        self.ax.grid(True)
        self.canvas.draw_idle()

    def redraw_plot(self):
        self.ax.cla()
        # One chart, one line per symbol; rely on default colors
        for s in self.alphabet:
            self.ax.plot(self.positions, self.freqs[s], label=s)

        self.ax.set_ylim(0.0, 1.0)
        self.ax.set_xlabel("Window start (1-based index)")
        self.ax.set_ylabel("Relative frequency")
        self.ax.set_title("Sliding-window nucleotide frequencies")
        self.ax.grid(True)
        self.ax.legend(title="Symbol")
        self.canvas.draw_idle()


if __name__ == "__main__":
    app = FastaWindowApp()
    app.mainloop()
