# download from NCBI a total of 3 bacterial genomes. Use these genomes as an input for your previous app.
# also modify your app in order to be able to handle the amount of information from these genomes. 
# your app must detect transposable elements and the results from the output must show their position and length. 
# in this case, the inverted repeats are unknown, unlike the previous assignments, in which they were known.
# !! the following thesis must be taken into consideration : 1. transposoms involving, 2. transposoms overlapping 
#!! the  minimum size of the inverted repeats must be of 4 bases and the max size must be of 6 bases

import sys
import csv
from Bio import Entrez, SeqIO
from Bio.Seq import Seq

# --- CONFIGURATION ---
Entrez.email = "student@university.edu" 

class GenomeAnalyzer:
    def __init__(self):
        self.genomes = {} 

    def download_genomes(self, accession_ids):
        print(f"--- Downloading {len(accession_ids)} genomes from NCBI ---")
        for acc in accession_ids:
            try:
                print(f"Fetching {acc}...")
                handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
                record = SeqIO.read(handle, "fasta")
                handle.close()
                self.genomes[acc] = str(record.seq).upper()
                print(f" > Downloaded {acc}: {len(self.genomes[acc])} bp")
            except Exception as e:
                print(f" ! Error downloading {acc}: {e}")

    def find_transposons(self, sequence, genome_name, min_ir=4, max_ir=6, min_dist=50, max_dist=1000):
        hits = [] 
        seq_len = len(sequence)
        scan_limit = min(seq_len, 20000) 

        print(f"Scanning {genome_name} (first {scan_limit} bases)...")

        id_counter = 1
        
        for k in range(min_ir, max_ir + 1):
            for i in range(0, scan_limit - k):
                start_seq = sequence[i : i+k]
                end_seq_target = str(Seq(start_seq).reverse_complement())
                
                search_start = i + min_dist
                search_end = min(i + max_dist, seq_len)
                
                if search_start >= seq_len: break

                window = sequence[search_start : search_end]
                rel_pos = window.find(end_seq_target)
                
                if rel_pos != -1:
                    match_start = i
                    match_end = search_start + rel_pos + k
                    
                    hits.append({
                        'genome': genome_name,
                        'id': f"{genome_name}_{id_counter}",
                        'start': match_start,
                        'end': match_end,
                        'len': match_end - match_start,
                        'ir_seq': start_seq,
                        'ir_len': k,
                        'notes': "Normal" # Default value
                    })
                    id_counter += 1
        
        unique_hits = self._filter_overlaps(hits)
        return unique_hits

    def _filter_overlaps(self, hits):
        # Sort by length of transposon
        hits.sort(key=lambda x: x['len'], reverse=True)
        kept = []
        for h in hits:
            is_duplicate = False
            for k in kept:
                if h['start'] == k['start'] and h['end'] == k['end']:
                    is_duplicate = True
                    break
            if not is_duplicate:
                kept.append(h)
        return kept

    def analyze_relationships(self, hits):
        """
        Updates the 'notes' field in the hits list with relationship info.
        """
        hits.sort(key=lambda x: x['start'])
        
        for i in range(len(hits)):
            for j in range(i + 1, len(hits)):
                t1 = hits[i]
                t2 = hits[j]
                
                s1, e1 = t1['start'], t1['end']
                s2, e2 = t2['start'], t2['end']

                # Case 1: NESTED
                if s1 < s2 and e2 < e1:
                    msg = f"NESTED inside {t1['id']}"
                    t2['notes'] = msg # Mark the inner one
                    print(f"[NESTED] {t2['id']} is inside {t1['id']}")

                # Case 2: OVERLAPPING
                elif s1 < s2 and s2 < e1 and e1 < e2:
                    msg = f"OVERLAPS with {t1['id']}"
                    t2['notes'] = msg
                    print(f"[OVERLAPPING] {t2['id']} overlaps {t1['id']}")

    def export_to_csv(self, all_hits, filename="transposon_results.csv"):
        """
        Writes all collected data to a CSV file.
        """
        if not all_hits:
            print("No data to write to CSV.")
            return

        print(f"\n--- Saving results to {filename} ---")
        
        # Define columns
        fieldnames = ['genome', 'id', 'start', 'end', 'len', 'ir_seq', 'notes']
        
        with open(filename, mode='w', newline='') as csv_file:
            writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
            
            writer.writeheader()
            for hit in all_hits:
                # We only write the fields we defined
                row = {field: hit[field] for field in fieldnames}
                writer.writerow(row)
        
        print("Save complete.")

# --- MAIN APP EXECUTION ---

if __name__ == "__main__":
    app = GenomeAnalyzer()
    
    # 3 small genomes/plasmids for testing
    accessions = ["NC_001416", "NC_002516", "NC_001318"] 
    app.download_genomes(accessions)
    
    all_results = []

    for acc_id, sequence in app.genomes.items():
        # 1. Find
        found_tes = app.find_transposons(sequence, acc_id, min_ir=4, max_ir=6, min_dist=20, max_dist=600)
        
        # 2. Analyze Relationships
        if len(found_tes) > 1:
            app.analyze_relationships(found_tes)
            
        # 3. Add to master list
        all_results.extend(found_tes)

    # 4. Export everything to CSV
    app.export_to_csv(all_results)