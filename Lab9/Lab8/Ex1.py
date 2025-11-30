#make an artificial DNA sequence of 200-400 bases in length, in which to simulate 3-4 transposable elements.
#implement a software app to detect the positions of these transposable elements(start/end) within the created DNA sequence 
# NOTE! intersect 2-3 transpons to see how the algorithm behaves.


import random
import re

def get_random_bases(length):
    """Generates a random DNA sequence of length N."""
    return ''.join(random.choice('ACGT') for _ in range(length))

def reverse_complement(seq):
    """Returns the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

class TransposonSimulator:
    def __init__(self, length=300):
        self.dna_sequence = get_random_bases(length)
        # Define a signature Inverted Repeat (IR) for our transposons
        # Whiteboard example suggests looking for specific patterns
        self.ir_seq = "CCGG" 
        self.ir_seq_rc = reverse_complement(self.ir_seq) # CCGG -> CCGG (Palindrome) or similar
        
        # Let's use a non-palindromic IR to be clear:
        self.ir_start = "TACT"
        self.ir_end = reverse_complement(self.ir_start) # AGTA
        
        self.insertions_log = []

    def insert_transposon(self, position, length=30):
        """
        Simulates a transposon insertion:
        1. Duplicates target site (Direct Repeats)
        2. Inserts IR_START + CONTENT + IR_END
        """
        # Ensure position is within bounds
        if position >= len(self.dna_sequence): return
        
        # 1. Simulate Target Site Duplication (Direct Repeats - e.g., 3 bases)
        target_site_len = 3
        target_site = self.dna_sequence[position:position+target_site_len]
        
        # 2. Create Transposon Body
        body_length = length - (len(self.ir_start) + len(self.ir_end))
        content = get_random_bases(body_length)
        transposon = self.ir_start + content + self.ir_end
        
        # 3. Construct the new sequence: 
        # Old_Pre + [Target_Copy1 (DR)] + Transposon + [Target_Copy2 (DR)] + Old_Post
        # Note: usually the target site exists, splits, and fills gaps. 
        # Simplification: We just insert the transposon + one copy of target site to simulate DRs.
        
        insertion_block = transposon + target_site 
        
        self.dna_sequence = (
            self.dna_sequence[:position+target_site_len] + 
            insertion_block + 
            self.dna_sequence[position+target_site_len:]
        )
        
        print(f"Inserted TE at original index {position}. Sequence length is now {len(self.dna_sequence)}")

    def detect_transposons(self):
        """
        Scans the DNA sequence for the defined Inverted Repeats (Start and End).
        """
        print("\n--- Running Detection Algorithm ---")
        sequence = self.dna_sequence
        start_pattern = self.ir_start
        end_pattern = self.ir_end
        
        # Find all start indices
        starts = [m.start() for m in re.finditer(start_pattern, sequence)]
        # Find all end indices (point to the end of the tag)
        ends = [m.end() for m in re.finditer(end_pattern, sequence)]
        
        detected = []
        
        # Logic to match Starts with Ends
        # NOTE: This handles the "Intersection" behavior.
        # If we have Start1... Start2... End2... End1
        # A simple greedy algo might match Start2-End2, but Start1-End1 contains it.
        
        for s in starts:
            for e in ends:
                # Minimal length check (e.g., transposon must be at least 10bp)
                if e > s + 10:
                    fragment = sequence[s:e]
                    # We look for the shortest valid pair first, or the longest?
                    # Let's simply report all valid pairings found
                    detected.append((s, e, fragment))

        if not detected:
            print("No transposons detected.")
        else:
            print(f"Found {len(detected)} potential candidates (including nested/overlaps):")
            for i, (s, e, frag) in enumerate(detected):
                print(f"Candidate {i+1}: Start {s}, End {e}, Length {len(frag)}")
                # Print preview
                print(f"   Seq: {frag[:10]}...{frag[-10:]}")

# --- MAIN EXECUTION ---

# 1. Make artificial DNA sequence (200-400 bases)
sim = TransposonSimulator(length=300)
print(f"Original Host DNA Length: {len(sim.dna_sequence)}")

# 2. Simulate 3-4 transposable elements
# We insert from back to front so indices don't mess up as easily during manual placement,
# OR we just insert dynamically.

# Insertion 1 (Normal)
sim.insert_transposon(position=50, length=40)

# Insertion 2 & 3 (INTERSECTING / NESTED)
# To test "how the algorithm behaves" (as per prompt), we put one INSIDE another.
# We insert a big one...
sim.insert_transposon(position=150, length=80)
# ...and then we insert another one right in the middle of the one we just made (approx pos 170)
sim.insert_transposon(position=170, length=30)

print(f"Final Sequence Length: {len(sim.dna_sequence)}")

# 3. Implement detection software
sim.detect_transposons()
