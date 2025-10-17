#establish a threshold that acts as a cutoff value, for each signal.
# plot a chart that shows the cutoff values as horizontal lines across the chart.


import argparse
import math
from collections import Counter
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# ---------- Core ----------
def clean_dna(seq: str) -> str:
    s = "".join(seq.upper().split()).replace("U", "T")
    return "".join(ch for ch in s if ch in "ACGT")

def read_fasta_sequence(path: Path) -> str:
    parts = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            parts.append(line)
    return clean_dna("".join(parts))

def gc_percent(seq: str) -> float:
    c = Counter(seq)
    return 100.0 * (c["G"] + c["C"]) / len(seq) if seq else 0.0

def tm_wallace(seq: str) -> float:
    c = Counter(seq)
    return 2.0 * (c["A"] + c["T"]) + 4.0 * (c["G"] + c["C"])

def tm_salt_adjusted(seq: str, na_mM: float = 50.0) -> float:
    if not seq:
        return float("nan")
    L = len(seq)
    na_M = max(na_mM, 0.0001) / 1000.0  # mM -> M; avoid log10(0)
    return 81.5 + 16.6 * math.log10(na_M) + 0.41 * gc_percent(seq) - (600.0 / L)

def sliding_tm_vectors(seq: str, k: int, na_mM: float):
    """Return (X, P_wallace, P_salt) for k-length windows."""
    if k < 1 or len(seq) < k:
        raise ValueError(f"Window must be 1..{len(seq)} (got {k}).")
    X, Pw, Ps = [], [], []
    for i in range(len(seq) - k + 1):
        win = seq[i:i+k]
        X.append(i + 1)                # 1-indexed start position
        Pw.append(tm_wallace(win))
        Ps.append(tm_salt_adjusted(win, na_mM))
    return np.array(X), np.array(Pw, float), np.array(Ps, float)

# --- helper to convert a boolean mask into contiguous [start,len] segments
def mask_to_segments(X: np.ndarray, mask: np.ndarray):
    """
    Given positions X (1..N) and a boolean mask of same length, return
    a list of (start, length) segments where mask is True, measured in x units.
    """
    segs = []
    if len(X) == 0:
        return segs
    # We assume X increases by 1 each step (start-of-window positions)
    in_seg = False
    start = None
    for i, flag in enumerate(mask):
        if flag and not in_seg:
            in_seg = True
            start = X[i]
        elif not flag and in_seg:
            in_seg = False
            stop = X[i-1]
            segs.append((int(start), int(stop - start + 1)))
    if in_seg:
        stop = X[-1]
        segs.append((int(start), int(stop - start + 1)))
    return segs

# ---------- Plot like the board ----------
def plot_board_style(X, Pw, Ps, k):
    # thresholds (cutoffs) — using the mean of each signal
    w_cut = float(np.mean(Pw))
    s_cut = float(np.mean(Ps))

    # Build figure with small “bar rows” under the main plot
    fig = plt.figure(figsize=(12, 6))
    gs = fig.add_gridspec(nrows=3, height_ratios=[4.0, 0.6, 0.6], hspace=0.07)
    ax = fig.add_subplot(gs[0])
    ax_w = fig.add_subplot(gs[1], sharex=ax)
    ax_s = fig.add_subplot(gs[2], sharex=ax)

    # Main curves
    line_w, = ax.plot(X, Pw, linewidth=1.8, label="Wallace Tm")
    line_s, = ax.plot(X, Ps, linewidth=1.8, label="Salt-adjusted Tm")

    # Dashed cutoff lines (board-like)
    ax.axhline(w_cut, linestyle="--", color=line_w.get_color(), alpha=0.9,
               label=f"Wallace cutoff (mean) = {w_cut:.2f}°C")
    ax.axhline(s_cut, linestyle="--", color=line_s.get_color(), alpha=0.9,
               label=f"Salt cutoff (mean) = {s_cut:.2f}°C")

    ax.set_ylabel("Melting temperature (°C)")
    ax.set_title(f"Sliding-Window Tm (k={k}) with Cutoff Lines and Activity Bars")
    ax.grid(alpha=0.35)
    ax.legend(ncols=2, fontsize=9)

    # ---- Bottom horizontal bars (where signal >= cutoff) ----
    # Wallace bars
    mask_w = Pw >= w_cut
    segs_w = mask_to_segments(X, mask_w)
    ax_w.set_yticks([]); ax_w.set_ylim(0, 1); ax_w.set_ylabel("W", rotation=0, labelpad=10, va="center")
    for start, length in segs_w:
        ax_w.broken_barh([(start, length)], (0.25, 0.5), facecolor=line_w.get_color(), alpha=0.85)
    ax_w.grid(alpha=0.25, axis="x")

    # Salt bars
    mask_s = Ps >= s_cut
    segs_s = mask_to_segments(X, mask_s)
    ax_s.set_yticks([]); ax_s.set_ylim(0, 1); ax_s.set_ylabel("S", rotation=0, labelpad=10, va="center")
    for start, length in segs_s:
        ax_s.broken_barh([(start, length)], (0.25, 0.5), facecolor=line_s.get_color(), alpha=0.85)
    ax_s.grid(alpha=0.25, axis="x")

    # shared x label
    ax_s.set_xlabel(f"Position in sequence (start of {k}-nt window)")
    # Make x tick labels on middle axis invisible to avoid clutter
    plt.setp(ax_w.get_xticklabels(), visible=False)

    fig.tight_layout()
    plt.show()

    # Console summary
    print(f"\nWindows: {len(X)} (k={k})")
    print(f"Wallace cutoff (mean): {w_cut:.2f} °C | segments: {segs_w}")
    print(f"Salt    cutoff (mean): {s_cut:.2f} °C | segments: {segs_s}")

# ---------- Main ----------
def main():
    ap = argparse.ArgumentParser(description="Board-style plot with cutoff lines and bottom activity bars.")
    ap.add_argument("fasta", type=str, help="Path to FASTA file")
    ap.add_argument("--window", type=int, default=9, help="Sliding window size (default 9)")
    ap.add_argument("--na-mM", type=float, default=50.0, help="[Na+] in mM for salt-adjusted formula")
    args = ap.parse_args()

    seq = read_fasta_sequence(Path(args.fasta))
    if not seq:
        raise SystemExit("No valid A/C/G/T sequence found in the FASTA file.")
    if len(seq) < args.window:
        raise SystemExit(f"Sequence length {len(seq)} < window size {args.window}.")
    print(f"Loaded len={len(seq)} | GC%={gc_percent(seq):.2f}%")

    X, Pw, Ps = sliding_tm_vectors(seq, args.window, args.na_mM)
    plot_board_style(X, Pw, Ps, args.window)

if __name__ == "__main__":
    main()
