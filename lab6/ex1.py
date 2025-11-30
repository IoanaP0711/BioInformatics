import random
import io

from Bio import Entrez, SeqIO
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------
# User parameters
# ------------------------------------------------------------------------------

# Email for NCBI Entrez 
EMAIL = "postelnecuioana@gmail.com"

# NCBI search term: E. coli sequences with length in [1000, 3000]
SEARCH_TERM = "Escherichia coli[Organism] AND 1000:3000[SLEN]"

# Lane label on the gel
LANE_LABEL = "Random fragments"

# Number of fragments to sample from the sequence
N_FRAG = 10

# Minimum and maximum fragment size (bp) to sample
FRAG_MIN = 100
FRAG_MAX = 3000

# Gel drawing parameters
GEL_HEIGHT = 500
GEL_WIDTH = 270
AGAROSE_PCT = 1.0  # used only to tweak the mapping slope

# Random seed for reproducibility
RANDOM_SEED = 42

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------


def sample_fragment(seq_len, length_bp):
    """
    Sample a fragment from a sequence of length `seq_len` with the requested length `length_bp`.
    If the requested length > seq_len, it is clipped to seq_len.
    Returns:
        (start, end, frag_len)
    """
    L = min(int(length_bp), seq_len)
    # If seq_len == L, start will be 0
    start = random.randint(0, max(0, seq_len - L))
    end = start + L
    return (start, end, L)


# ------------------------------------------------------------------------------
# Main script
# ------------------------------------------------------------------------------

def main():
    # Initialize Entrez
    Entrez.email = EMAIL

    # Search for candidate E. coli sequences
    print(f"[INFO] Searching NCBI with term: {SEARCH_TERM!r}")
    search = Entrez.esearch(db="nucleotide", term=SEARCH_TERM, retmax=50)
    result = Entrez.read(search)
    ids = result["IdList"]

    if not ids:
        raise RuntimeError("No records found for the given search term.")

    print(f"[INFO] Found {len(ids)} candidate sequences.")
    random.seed(RANDOM_SEED)
    chosen_id = random.choice(ids)
    print(f"[INFO] Randomly chosen sequence ID: {chosen_id}")

    # Fetch the chosen sequence in FASTA format
    handle = Entrez.efetch(
        db="nucleotide", id=chosen_id, rettype="fasta", retmode="text"
    )
    fasta_text = handle.read()
    handle.close()

    # Parse the FASTA
    seq_record = SeqIO.read(io.StringIO(fasta_text), "fasta")
    seq_len = len(seq_record.seq)
    print(f"[INFO] Chosen sequence length: {seq_len} nt")

    # Generate fragment lengths
    # We build a pool of possible lengths from FRAG_MIN up to min(FRAG_MAX, seq_len).
    length_pool = list(range(FRAG_MIN, min(FRAG_MAX, seq_len) + 1))
    if len(length_pool) < N_FRAG:
        raise RuntimeError(
            f"Not enough possible lengths in [{FRAG_MIN}, {FRAG_MAX}] for a sequence of length {seq_len}."
        )

    chosen_lengths = random.sample(length_pool, k=N_FRAG)
    print(f"[INFO] Chosen fragment lengths (before sampling positions): {chosen_lengths}")

    # Sample fragment positions
    fragments = [sample_fragment(seq_len, L) for L in chosen_lengths]
    # We only need the lengths for plotting on a gel
    sizes_bp = np.array(sorted(int(L) for (_, _, L) in fragments), dtype=float)

    # --------------------------------------------------------------------------
    # Map fragment sizes to gel positions
    # --------------------------------------------------------------------------
    # We use a log10-based mapping reminiscent of DNA migration in agarose gels:
    #   y ≈ a - b * (log10(bp) - log_min) / (log_max - log_min)
    #
    # so that large fragments (e.g. ~3000 bp) remain near the top,
    # and small fragments (e.g. ~100 bp) run closer to the bottom.

    a = GEL_HEIGHT * 0.85
    # A "shape" factor depending on gel concentration (just a simple tweak)
    b = GEL_HEIGHT * (0.55 + 0.15 * AGAROSE_PCT)

    log_min = np.log10(min(sizes_bp.min(), 150.0))
    log_max = np.log10(max(sizes_bp.max(), 3000.0))

    def y_from_bp(bp):
        """Convert fragment size (bp) to a vertical position y in the gel."""
        # Normalize
        y = a - b * (np.log10(bp) - log_min) / (log_max - log_min + 1e-9)
        # Clip to stay nicely within the gel area
        return float(np.clip(y, 40, GEL_HEIGHT - 30))

    # Compute positions
    y_positions = [y_from_bp(size) for size in sizes_bp]

    # --------------------------------------------------------------------------
    # Plot the gel
    # --------------------------------------------------------------------------
    plt.figure(figsize=(3.3, 7.6), dpi=150)
    ax = plt.gca()
    ax.set_xlim(0, GEL_WIDTH)
    ax.set_ylim(GEL_HEIGHT, 0)  # invert y-axis like a real gel

    # Background (gel)
    ax.add_patch(
        plt.Rectangle((0, 0), GEL_WIDTH, GEL_HEIGHT, color="black", zorder=0)
    )

    # Single lane for our random fragments
    lane_x_center = GEL_WIDTH * 0.65
    lane_half_width = 22

    ax.add_patch(
        plt.Rectangle(
            (lane_x_center - lane_half_width, 10),
            2 * lane_half_width,
            GEL_HEIGHT - 20,
            color="#222222",
            zorder=1,
        )
    )

    # Plot the bands
    band_height = 4
    for y in y_positions:
        ax.add_patch(
            plt.Rectangle(
                (lane_x_center - lane_half_width + 2, y - band_height / 2),
                2 * lane_half_width - 4,
                band_height,
                color="#EEEEEE",
                zorder=2,
            )
        )

    # Lane label at the top
    ax.text(
        lane_x_center,
        20,
        LANE_LABEL,
        color="white",
        ha="center",
        va="bottom",
        fontsize=9,
        rotation=0,
    )

    # Ladder on the left side
    ladder_x_center = GEL_WIDTH * 0.20
    ladder_half_width = 18
    ax.add_patch(
        plt.Rectangle(
            (ladder_x_center - ladder_half_width, 10),
            2 * ladder_half_width,
            GEL_HEIGHT - 20,
            color="#222222",
            zorder=1,
        )
    )

    # Example ladder sizes in bp
    ladder_sizes = [3000, 1500, 1000, 700, 500, 400, 300, 200, 100]
    for bp in ladder_sizes:
        y = y_from_bp(bp)
        ax.add_patch(
            plt.Rectangle(
                (ladder_x_center - ladder_half_width + 2, y - band_height / 2),
                2 * ladder_half_width - 4,
                band_height,
                color="#DDDDDD",
                zorder=2,
            )
        )
        ax.text(
            ladder_x_center - ladder_half_width - 4,
            y,
            f"{bp}",
            color="white",
            ha="right",
            va="center",
            fontsize=7,
        )

    # Final plot tweaks
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)
    plt.tight_layout()

    # Show the gel
    plt.show()

    # Print summary
    print("------------------------------------------------------------")
    print(f"NCBI ID: {chosen_id}")
    print(f"Sequence length: {seq_len} nt")
    print("Distinct fragment sizes (bp):", [int(x) for x in sizes_bp])
    print("------------------------------------------------------------")


if __name__ == "__main__":
    main()
