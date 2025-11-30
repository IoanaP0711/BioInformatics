import random
import io
import textwrap

from Bio import Entrez, SeqIO
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------
# User parameters
# ------------------------------------------------------------------------------


EMAIL = "postelnecuioana@gmail.com"

# Number of influenza genomes we want to fetch
N_GENOMES = 10

# Random seed for reproducible choices
RANDOM_SEED = 11

# General gel figure parameters
BASE_FIG_DPI = 150

# For combined gel
COMBINED_GEL_HEIGHT = 1000  # pixels
LANE_WIDTH = 36
LANE_SPACING = 60
LADDER_X = 90

# Ladder sizes (bp) to display on left side
LADDER_SIZES = [12000, 10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 700, 500, 400, 300, 200, 100]

# Font sizes
LABEL_FONT = 8
ANNOT_LABEL_FONT = 7
TICK_FONT = 7

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------


def esearch_ids(queries, need):
    """
    Try each query in `queries` against the NCBI Nucleotide database
    until we collect at least `need` unique IDs (or exhaust all queries).
    """
    seen = []
    for q in queries:
        print(f"[INFO] ESearch with query: {q!r}")
        h = Entrez.esearch(db="nucleotide", term=q, retmax=500)
        r = Entrez.read(h)
        ids = r.get("IdList", [])
        print(f"[INFO]  Found {len(ids)} IDs")
        for i in ids:
            if i not in seen:
                seen.append(i)
            if len(seen) >= need:
                return seen[:need]
    return seen


def fetch_fasta_concat(nid):
    """
    Fetch a nucleotide record in FASTA format from NCBI and:
      - if there are multiple records, concatenate them;
      - otherwise, return the single sequence.

    Returns:
        label, seq_string
    """
    h = Entrez.efetch(db="nucleotide", id=nid, rettype="fasta", retmode="text")
    text = h.read()
    h.close()

    recs = list(SeqIO.parse(io.StringIO(text), "fasta"))
    if not recs:
        return nid, ""
    if len(recs) == 1:
        label = recs[0].id
        seq = str(recs[0].seq).upper()
    else:
        label = recs[0].id
        seq = "".join(str(r.seq).upper() for r in recs)

    return label, seq


def ecoRI_fragments(seq):
    """
    Simulate digestion by the restriction enzyme EcoRI:
      - Recognition motif: GAATTC
      - Cut between G and A (offset 1 from motif start)

    Returns:
        list of fragment lengths in bp
    """
    motif = "GAATTC"
    cut_offset = 1

    s = 0
    sites = []
    while True:
        i = seq.find(motif, s)
        if i == -1:
            break
        sites.append(i + cut_offset)
        s = i + 1

    cuts = [0] + sorted(sites) + [len(seq)]
    fragments = [cuts[i + 1] - cuts[i] for i in range(len(cuts) - 1)]
    return fragments


def y_from_bp(bp, a, b, log_min, log_max, H):
    """
    Convert bp to a vertical coordinate using a log10 mapping.
    The result is clipped to [40, H-40] in the gel.
    """
    y = a - b * (np.log10(bp) - log_min) / (log_max - log_min + 1e-9)
    return float(np.clip(y, 40, H - 40))


def plot_ladder(ax, x_left, a, b, log_min, log_max, H, fs_tick=TICK_FONT):
    """
    Draw the DNA ladder and its labels along the left side.
    """
    ladder_w = 30
    ax.add_patch(
        plt.Rectangle(
            (x_left, 10),
            ladder_w,
            H - 20,
            color="#111111",
            zorder=1,
        )
    )

    band_h = 4
    for bp in LADDER_SIZES:
        y = y_from_bp(bp, a, b, log_min, log_max, H)
        ax.add_patch(
            plt.Rectangle(
                (x_left + 2, y - band_h / 2),
                ladder_w - 4,
                band_h,
                color="#DDDDDD",
                zorder=2,
            )
        )
        # text to the left of ladder
        ax.text(
            x_left - 6,
            y,
            f"{bp}",
            color="white",
            ha="right",
            va="center",
            fontsize=fs_tick,
        )


def plot_lane(ax, x_center, H, sizes, a, b, log_min, log_max, color="white", annotate=False, fsize=ANNOT_LABEL_FONT):
    """
    Draw one gel lane at horizontal position x_center.
    If annotate=True, also label each fragment size (bp) at the right margin.
    """
    lane_half = LANE_WIDTH // 2
    lane_h = H - 20

    # lane rectangle
    ax.add_patch(
        plt.Rectangle(
            (x_center - lane_half, 10),
            2 * lane_half,
            lane_h,
            color="#191919",
            zorder=1,
        )
    )

    band_h = 4
    sizes_sorted = sorted(sizes, reverse=False)

    for bp in sizes_sorted:
        y = y_from_bp(bp, a, b, log_min, log_max, H)
        ax.add_patch(
            plt.Rectangle(
                (x_center - lane_half + 2, y - band_h / 2),
                2 * lane_half - 4,
                band_h,
                color=color,
                zorder=2,
            )
        )

        # Optionally annotate each band with its size (for detailed gels)
        if annotate:
            ax.text(
                x_center + lane_half + 5,
                y,
                f"{bp}",
                color="white",
                ha="left",
                va="center",
                fontsize=fsize,
            )


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------


def main():
    Entrez.email = EMAIL

    # Define several queries to gather influenza genomes
    queries = [
        "Influenza A virus[ORGN] AND 13000:16000[SLEN]",
        "(influenza[ALL]) AND 12000:17000[SLEN]",
        "Influenza B virus[ORGN] AND 14000:16000[SLEN]",
        "Influenza C virus[ORGN] AND 12000:16000[SLEN]",
    ]

    # Get up to N_GENOMES distinct IDs
    ids = esearch_ids(queries, N_GENOMES)
    if len(ids) < N_GENOMES:
        raise RuntimeError(
            f"Not enough influenza genome IDs found; needed {N_GENOMES}, got {len(ids)}."
        )

    ids = ids[:N_GENOMES]
    print(f"[INFO] Using these {N_GENOMES} IDs: {ids}")

    random.seed(RANDOM_SEED)

    labels = []
    genomes = []

    # Fetch sequences
    for nid in ids:
        label, seq = fetch_fasta_concat(nid)
        if not seq:
            print(f"[WARN] Empty sequence for ID {nid}, skipping.")
            continue
        labels.append(label)
        genomes.append(seq)
        print(f"[INFO] Got sequence for {label} (len={len(seq)})")

    # Digest with EcoRI
    digests = []
    gen_lengths = []
    for seq in genomes:
        frags = ecoRI_fragments(seq)
        digests.append(frags)
        gen_lengths.append(len(seq))

    # Collect all fragments (for global min/max)
    all_frags = [f for sizes in digests for f in sizes if f > 0]
    if not all_frags:
        raise RuntimeError("No fragments were produced. Check motif, queries or sequences.")

    # Global log10 min/max for mapping in combined gel
    overall_min = max(50, min(all_frags))
    overall_max = max(all_frags)
    log_min = np.log10(overall_min)
    log_max = np.log10(overall_max)

    # y = a - b * normalized(log10(bp))
    H = COMBINED_GEL_HEIGHT
    a = H * 0.9
    b = H * 0.7

    # --------------------------------------------------------------------------
    # Combined gel with ladder and all 10 lanes
    # --------------------------------------------------------------------------
    W = 240 + (len(digests) + 1) * LANE_SPACING

    fig_w = W / BASE_FIG_DPI
    fig_h = H / BASE_FIG_DPI

    plt.figure(figsize=(fig_w, fig_h), dpi=BASE_FIG_DPI)
    ax = plt.gca()
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)  # inverted like a real gel

    # Background
    ax.add_patch(
        plt.Rectangle(
            (0, 0),
            W,
            H,
            color="black",
            zorder=0,
        )
    )

    # Ladder on the left
    plot_ladder(ax, LADDER_X, a, b, log_min, log_max, H, fs_tick=TICK_FONT)

    # Each genome lane
    x0 = LADDER_X + 110
    for i, sizes in enumerate(digests):
        x_center = x0 + i * LANE_SPACING
        plot_lane(
            ax,
            x_center=x_center,
            H=H,
            sizes=sizes,
            a=a,
            b=b,
            log_min=log_min,
            log_max=log_max,
            color="white",
            annotate=False,
        )

        # Label at the top (wrapped)
        label = labels[i]
        short_label = "\n".join(textwrap.wrap(label, width=18))
        ax.text(
            x_center,
            18,
            short_label,
            color="white",
            ha="center",
            va="bottom",
            fontsize=LABEL_FONT,
            rotation=0,
        )

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)
    plt.tight_layout()

    combined_png = "gel_combined.png"
    plt.savefig(combined_png, bbox_inches="tight", dpi=BASE_FIG_DPI)
    print(f"[INFO] Saved combined gel: {combined_png}")

    # --------------------------------------------------------------------------
    # Individual gels for each genome (with annotated fragment sizes)
    # --------------------------------------------------------------------------
    for i, (lab, sizes, total_len) in enumerate(zip(labels, digests, gen_lengths), start=1):
        frag_list = [s for s in sizes if s > 0]
        local_min = max(50, min(frag_list))
        local_max = max(frag_list)
        log_min_loc = np.log10(local_min)
        log_max_loc = np.log10(local_max)

        H_i = 900
        W_i = 500
        a_i = H_i * 0.88
        b_i = H_i * 0.75

        fig_w_i = W_i / BASE_FIG_DPI
        fig_h_i = H_i / BASE_FIG_DPI

        plt.figure(figsize=(fig_w_i, fig_h_i), dpi=BASE_FIG_DPI)
        ax_i = plt.gca()
        ax_i.set_xlim(0, W_i)
        ax_i.set_ylim(H_i, 0)

        # Background
        ax_i.add_patch(
            plt.Rectangle(
                (0, 0),
                W_i,
                H_i,
                color="black",
                zorder=0,
            )
        )

        x_center_i = W_i * 0.4
        plot_lane(
            ax_i,
            x_center=x_center_i,
            H=H_i,
            sizes=frag_list,
            a=a_i,
            b=b_i,
            log_min=log_min_loc,
            log_max=log_max_loc,
            color="white",
            annotate=True,
            fsize=ANNOT_LABEL_FONT,
        )

        title = f"{lab}\nEcoRI digest ({len(frag_list)} fragments, {total_len} nt)"
        wrapped_title = "\n".join(textwrap.wrap(title, width=50))
        ax_i.text(
            W_i * 0.5,
            25,
            wrapped_title,
            color="white",
            ha="center",
            va="bottom",
            fontsize=LABEL_FONT,
        )

        ax_i.set_xticks([])
        ax_i.set_yticks([])
        ax_i.set_frame_on(False)
        plt.tight_layout()

        out_name = f"gel_{i}_{lab}.png"
        # Clean filename a bit (remove spaces or problematic chars)
        out_name = out_name.replace(" ", "_").replace("/", "_")

        plt.savefig(out_name, bbox_inches="tight", dpi=BASE_FIG_DPI)
        print(f"[INFO] Saved lane gel: {out_name}")

    # --------------------------------------------------------------------------
    # Print summary about longest genome
    # --------------------------------------------------------------------------
    longest_idx = int(np.argmax(gen_lengths))
    print("------------------------------------------------------------")
    print("Genome_labels:", labels)
    print("Genome_lengths_nt:", gen_lengths)
    print("Longest_genome_index_1based:", longest_idx + 1)
    print("Longest_genome_label:", labels[longest_idx])
    print("Longest_genome_length_nt:", gen_lengths[longest_idx])
    print("------------------------------------------------------------")


if __name__ == "__main__":
    main()
