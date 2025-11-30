import matplotlib.pyplot as plt
import matplotlib.patches as patches
import re
from Bio import Entrez, SeqIO
from collections import defaultdict

# ==========================================
# 1. SETUP & INPUTS
# ==========================================
Entrez.email = "your.email@example.com"  # REQUIRED: Change this to your email

# Your enzyme list
restriction_enzymes = {
    "EcoRI":   {"pattern": "GAATTC", "cut_offset": 1},
    "BamHI":   {"pattern": "GGATCC", "cut_offset": 1},
    "HindIII": {"pattern": "AAGCTT", "cut_offset": 1},
    "TaqI":    {"pattern": "TCGA",   "cut_offset": 1},
    "HaeIII":  {"pattern": "GGCC",   "cut_offset": 2} 
}

# ==========================================
# 2. FETCH 10 VARIANTS (NCBI)
# ==========================================
def fetch_influenza_variants(n=10):
    print(f"Fetching {n} Influenza A (H3N2) HA gene variants from NCBI...")
    
    try:
        # 1. Search for IDs
        search_handle = Entrez.esearch(db="nucleotide", term=term, retmax=n, sort="relevance")
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        id_list = search_results["IdList"]
        
        if not id_list:
            raise ValueError("No sequences found.")

        # 2. Fetch actual sequences
        fetch_handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="fasta", retmode="text")
        records = list(SeqIO.parse(fetch_handle, "fasta"))
        fetch_handle.close()
        
        # Clean up names for the plot
        data = []
        for record in records:
            # Extract a short name (e.g., "A/New York/392/2004") from description
            match = re.search(r'Influenza A virus \((.*?)\)', record.description)
            short_name = match.group(1) if match else record.id[:10]
            data.append({"name": short_name, "seq": str(record.seq).upper()})
            
        print(f"Successfully downloaded {len(data)} variants.")
        return data

    except Exception as e:
        print(f"Error fetching data: {e}")
        return []

# ==========================================
# 3. DIGESTION LOGIC
# ==========================================
def digest_sequence(dna_seq, enzymes):
    """Digests a single sequence with ALL enzymes and pools fragments."""
    all_fragments = []
    
    for name, data in enzymes.items():
        pattern = data["pattern"]
        offset = data["cut_offset"]
        
        cut_sites = [m.start() + offset for m in re.finditer(pattern, dna_seq)]
        all_sites = [0] + sorted(cut_sites) + [len(dna_seq)]
        fragments = [all_sites[i+1] - all_sites[i] for i in range(len(all_sites)-1)]
        
        # We only care about fragment sizes for the gel
        all_fragments.extend(fragments)
        
    return sorted(all_fragments, reverse=True)

# ==========================================
# 4. ANALYSIS: REMOVE COMMON BANDS
# ==========================================
def filter_common_bands(variants_data):
    """
    1. Identifies bands present in ALL variants (Common).
    2. Removes them from each variant's list.
    """
    
    # Step A: Digest everyone first
    for variant in variants_data:
        variant['fragments'] = digest_sequence(variant['seq'], restriction_enzymes)

    common_candidates = set(variants_data[0]['fragments'])
    
    for variant in variants_data[1:]:
        current_fragments = set(variant['fragments'])
        # Intersection: Keep only bands present in BOTH sets
        common_candidates = common_candidates.intersection(current_fragments)
        
    print(f"\nidentified {len(common_candidates)} common bands shared by all variants (to be hidden).")
    print(f"Common sizes: {sorted(list(common_candidates))}")

    # Step C: Remove these common bands from everyone
    for variant in variants_data:
        original = variant['fragments']
        # Filter out if it's in the common set
        filtered = [f for f in original if f not in common_candidates]
        variant['unique_fragments'] = filtered
        
    return variants_data

# ==========================================
# 5. VISUALIZATION: MERGED GEL
# ==========================================
def plot_merged_gel(variants_data):
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_facecolor('black')
    fig.patch.set_facecolor('#1a1a1a')
    
    lanes = [v['name'] for v in variants_data]
    x_positions = range(len(lanes))
    
    # Plotting
    for i, variant in enumerate(variants_data):
        fragments = variant['unique_fragments']
        lane_x = i
        
        for frag_len in fragments:
            # Draw band
            # Red tint for unique bands to make them look "warning" / distinct
            ax.hlines(y=frag_len, xmin=lane_x - 0.35, xmax=lane_x + 0.35, 
                      colors='#ff5555', alpha=0.9, linewidth=2)
            
    # Styling
    ax.set_xticks(x_positions)
    ax.set_xticklabels(lanes, color='white', fontsize=10, rotation=45, ha='right')
    ax.set_ylabel("Fragment Size (bp)", color='white')
    ax.set_title("Differential Fingerprinting: Influenza A (H3N2) HA Variants\n(Conserved bands removed; only differences shown)", color='white', fontsize=14, pad=20)
    
    ax.set_yscale('log')
    ax.set_ylim(50, 4000) # HA gene is ~1700bp, fragments will be smaller
    
    # Spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_color('white')
    ax.spines['left'].set_color('white')
    ax.tick_params(colors='white')

    # Add wells
    for x in x_positions:
        rect = patches.Rectangle((x-0.35, 3800), 0.7, 200, linewidth=1, edgecolor='gray', facecolor='#222')
        ax.add_patch(rect)

    plt.tight_layout()
    plt.show()

# ==========================================
# MAIN
# ==========================================
if __name__ == "__main__":
    # 1. Get 10 Variants
    data = fetch_influenza_variants(10)
    
    if data:
        # 2. & 3. Digest and Filter
        processed_data = filter_common_bands(data)
        
        # 4. Visualize
        plot_merged_gel(processed_data)
