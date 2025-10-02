# design an app which is able to calculate the relative frequencies of each symbol from the alphabet of a given sequence. Use the sequence S ="ATTGCCCCGAAT"

from collections import Counter

def find_alphabet(sequence):
    # Using a set to store unique symbols
    alphabet = set(sequence)
    return alphabet

def calculate_relative_frequencies(sequence):
    # Count occurrences of each symbol
    counts = Counter(sequence)
    total = len(sequence)
    # Calculate relative frequency for each symbol
    frequencies = {symbol: counts[symbol] / total for symbol in counts}
    return frequencies

if __name__ == "__main__":
    # Example sequence
    S = "ATTGCCCCGAAT"
    
    # Find alphabet
    alphabet = find_alphabet(S)
    print("Sequence:", S)
    print("Alphabet:", alphabet)
    
    # Calculate relative frequencies
    rel_freqs = calculate_relative_frequencies(S)
    print("Relative Frequencies:")
    for symbol, freq in rel_freqs.items():
        print(f"  {symbol}: {freq:.3f}")
