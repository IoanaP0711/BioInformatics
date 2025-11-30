import matplotlib.pyplot as plt
import matplotlib.patches as patches
import re
import random
from Bio import Entrez, SeqIO

restriction_enzymes = {
    "EcoRI":   {"pattern": "GAATTC", "cut_offset": 1},
    "BamHI":   {"pattern": "GGATCC", "cut_offset": 1},
    "HindIII": {"pattern": "AAGCTT", "cut_offset": 1},
    "TaqI":    {"pattern": "TCGA",   "cut_offset": 1},
    "HaeIII":  {"pattern": "GGCC",   "cut_offset": 2} 
}

# ==========================================
# 2. DATA ACQUISITION (NCBI)
# ==========================================

def get_dna_sequence():
    """
    Fetches a sequence from NCBI. 
    If it fails (no internet/bad key), generates a random one.
    """
    try:
        # Fetching a fragment of the Human Beta-Globin gene (HBB)
        # We ask for a range to keep it between 1000-3000bp as requested
        print("Attempting to fetch DNA from NCBI...")
        handle = Entrez.efetch(db="nucleotide", id="NM_000518", rettype="fasta", retmode="text", seq_start=1, seq_stop=2500)
        record = SeqIO.read(handle, "fasta")
        handle.close()
        
        sequence = str(record.seq).upper()
        print(f"Successfully fetched {record.id}: {len(sequence)} bp")
        return sequence, record.description

    except Exception as e:
        print(f"NCBI fetch failed or skipped ({e}). Generating random DNA sequence instead.")
        # Fallback: Generate random DNA between 1000 and 3000 bp
        length = random.randint(1000, 3000)
        bases = ["A", "T", "C", "G"]
        return "".join(random.choice(bases) for _ in range(length)), "Randomly Generated DNA"

# ==========================================
# 3. DIGESTION LOGIC
# ==========================================

def digest_dna(dna_seq, enzymes):
    results = {}
    
    for name, data in enzymes.items():
        pattern = data["pattern"]
        offset = data["cut_offset"]
        
        # Find all occurrences of the pattern
        # re.finditer finds all matches. We add the offset to the start index to get the cut site.
        cut_sites = [m.start() + offset for m in re.finditer(pattern, dna_seq)]
        
        # Add start (0) and end (len) of DNA to calculate fragment lengths
        # We sort them to ensure calculations are linear
        all_sites = [0] + sorted(cut_sites) + [len(dna_seq)]
        
        # Calculate lengths: Distance between current site and next site
        fragment_lengths = [all_sites[i+1] - all_sites[i] for i in range(len(all_sites)-1)]
        
        # Sort fragments descending (standard for gel visualization data lists)
        fragment_lengths.sort(reverse=True)
        
        results[name] = {
            "num_cuts": len(cut_sites),
            "cut_positions": cut_sites,
            "fragments": fragment_lengths
        }
        
    return results

# ==========================================
# 4. GEL SIMULATION (VISUALIZATION)
# ==========================================

def simulate_gel(digestion_results, sequence_name):
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Visual Styling to match the "Dark" look of a real UV gel
    ax.set_facecolor('black')
    fig.patch.set_facecolor('#1a1a1a') # Dark grey border
    
    # Create Lanes
    enzymes = list(digestion_results.keys())
    lanes = ["Ladder"] + enzymes
    x_positions = range(len(lanes))
    
    # 1. Draw the Ladder (Lane 0) - A standard 100bp ladder simulation
    ladder_bands = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1500, 2000, 3000]
    for band in ladder_bands:
        # Draw a glowing white line
        ax.hlines(y=band, xmin=-0.3, xmax=0.3, colors='white', alpha=0.8, linewidth=2)
        
    # 2. Draw Enzyme Digests (Lanes 1-5)
    print("\n--- DETAILED OUTPUT ---")
    print(f"Sequence: {sequence_name}")
    
    for i, enzyme in enumerate(enzymes):
        lane_idx = i + 1
        data = digestion_results[enzyme]
        fragments = data['fragments']
        
        # Print text stats to console
        print(f"\nEnzyme: {enzyme}")
        print(f"  - Cuts: {data['num_cuts']}")
        print(f"  - Fragments (bp): {fragments}")
        
        for frag_len in fragments:
            # We add a little random "wobble" to x to simulate real gel imperfections
            wobble = random.uniform(-0.02, 0.02)
            
            # Draw the band
            # Thicker lines for larger DNA (more intensity usually)
            intensity = 0.9 if frag_len > 500 else 0.6
            ax.hlines(y=frag_len, xmin=lane_idx - 0.3 + wobble, xmax=lane_idx + 0.3 + wobble, 
                      colors='white', alpha=intensity, linewidth=2.5)

    # Chart Configuration
    ax.set_xticks(x_positions)
    ax.set_xticklabels(lanes, color='white', fontsize=12, rotation=45)
    ax.set_ylabel("Base Pairs (bp)", color='white', fontsize=12)
    ax.set_title(f"In Silico Restriction Digest\n{sequence_name}", color='white', fontsize=14, pad=20)
    
    # Handle Y-Axis
    # Real gels are logarithmic-ish (small stuff runs fast/far, big stuff runs slow/short).
    # Matplotlib plots 0 at bottom by default. 
    # To look like a gel, we want High Numbers at Top, Low Numbers at Bottom.
    # Matplotlib standard Y-axis already does this (0 -> Max), which fits the "Ladder" visualization visually 
    # if we treat the Y-axis as "Size".
    ax.set_yscale('log') # Log scale is essential for DNA gels
    ax.set_ylim(50, 3500) # Set limits to frame our data well
    
    # Remove top and right spines for clean look
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('white')
    ax.spines['bottom'].set_color('white')
    ax.tick_params(axis='x', colors='white')
    ax.tick_params(axis='y', colors='white')

    # Add "Wells" at the top
    for x in x_positions:
        rect = patches.Rectangle((x-0.3, 3300), 0.6, 200, linewidth=1, edgecolor='gray', facecolor='#222')
        ax.add_patch(rect)
    
    plt.tight_layout()
    plt.show()

# ==========================================
# MAIN EXECUTION
# ==========================================

if __name__ == "__main__":
    # 1. Get DNA
    dna_sequence, description = get_dna_sequence()
    
    # 2. Perform Digestion
    results = digest_dna(dna_sequence, restriction_enzymes)
    
    # 3. Output Text and Simulate Gel
    simulate_gel(results, description)