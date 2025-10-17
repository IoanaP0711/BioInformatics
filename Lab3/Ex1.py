import math
from collections import Counter

# ---------- Core helpers ----------
def clean_and_validate(seq: str) -> str:
    """Uppercase, strip spaces/newlines, ensure only A/C/G/T. Automatically converts U → T."""
    s = "".join(seq.upper().split()).replace("U", "T")  # Convert RNA bases to DNA
    allowed = set("ACGT")
    bad = sorted({ch for ch in s if ch not in allowed})
    if not s:
        raise ValueError("Sequence is empty.")
    if bad:
        raise ValueError(f"Invalid characters found: {', '.join(bad)}. Use only A/C/G/T.")
    return s

def gc_percent(seq: str) -> float:
    """Compute GC percentage."""
    c = Counter(seq)
    return 100.0 * (c["G"] + c["C"]) / len(seq)

# ---------- Tm formulas ----------
def tm_wallace(seq: str) -> float:
    """Wallace rule for short primers/oligos: Tm = 2*(A+T) + 4*(G+C)."""
    s = clean_and_validate(seq)
    c = Counter(s)
    return 2.0 * (c["A"] + c["T"]) + 4.0 * (c["G"] + c["C"])

def tm_salt_adjusted(seq: str, na_mM: float = 50.0) -> float:
    """
    Salt-adjusted formula (basic):
      Tm = 81.5 + 16.6*log10([Na+]) + 0.41*(%GC) - 600/length
    [Na+] is in mol/L. We accept millimolar for convenience.
    """
    s = clean_and_validate(seq)
    L = len(s)
    if L == 0:
        raise ValueError("Sequence length must be > 0.")
    na_M = max(na_mM, 0.0001) / 1000.0  # convert mM -> M; avoid log10(0)
    return 81.5 + 16.6 * math.log10(na_M) + 0.41 * gc_percent(s) - (600.0 / L)

# ---------- Simple CLI ----------
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Calculate DNA melting temperature (Tm).")
    parser.add_argument("sequence", nargs="?", help="DNA sequence (A/C/G/T). If omitted, you’ll be prompted.")
    parser.add_argument("--na-mM", type=float, default=50.0, help="Sodium concentration in mM (default: 50).")
    parser.add_argument("--gui", action="store_true", help="Open a simple GUI instead of CLI.")
    args = parser.parse_args()

    if args.gui:
        # ---------- Minimal Tkinter GUI ----------
        import tkinter as tk
        from tkinter import messagebox

        def on_calc():
            try:
                seq = seq_var.get()
                na = float(na_var.get())
                seq_clean = clean_and_validate(seq)
                t1 = tm_wallace(seq_clean)
                t2 = tm_salt_adjusted(seq_clean, na)
                out_var.set(
                    f"Input sequence (U→T converted): {seq_clean}\n"
                    f"Length: {len(seq_clean)}\n"
                    f"GC%: {gc_percent(seq_clean):.2f}%\n"
                    f"Wallace Tm: {t1:.2f} °C\n"
                    f"Salt-adjusted Tm: {t2:.2f} °C (Na+={na} mM)"
                )
            except Exception as e:
                messagebox.showerror("Error", str(e))

        root = tk.Tk()
        root.title("DNA Melting Temperature (Tm)")
        tk.Label(root, text="DNA sequence (A/C/G/T or RNA with U):").grid(row=0, column=0, sticky="w", padx=6, pady=6)
        seq_var = tk.StringVar()
        tk.Entry(root, width=50, textvariable=seq_var).grid(row=0, column=1, padx=6, pady=6)

        tk.Label(root, text="[Na+] (mM):").grid(row=1, column=0, sticky="w", padx=6, pady=6)
        na_var = tk.StringVar(value=str(args.na_mM))
        tk.Entry(root, width=12, textvariable=na_var).grid(row=1, column=1, sticky="w", padx=6, pady=6)

        tk.Button(root, text="Compute Tm", command=on_calc).grid(row=2, column=0, columnspan=2, pady=8)
        out_var = tk.StringVar()
        tk.Label(root, textvariable=out_var, justify="left").grid(row=3, column=0, columnspan=2, padx=6, pady=6)
        root.mainloop()
    else:
        # CLI path
        seq = args.sequence or input("Enter DNA sequence (A/C/G/T or RNA with U): ").strip()
        try:
            seq_clean = clean_and_validate(seq)
            t1 = tm_wallace(seq_clean)
            t2 = tm_salt_adjusted(seq_clean, args.na_mM)
            print(f"Input sequence (U→T converted): {seq_clean}")
            print(f"Length: {len(seq_clean)}")
            print(f"GC%: {gc_percent(seq_clean):.2f}%")
            print(f"Wallace Tm: {t1:.2f} °C")
            print(f"Salt-adjusted Tm: {t2:.2f} °C (Na+={args.na_mM} mM)")
        except Exception as e:
            print("Error:", e)
