# design an app which is able to find the alphabet of a given sequence. The alphabet means the unique symbols from which the sequence is made of.
# eg: sequence S ="ATTGCCCCGAAT"
# find the alphabet of the sequence S

def find_alphabet(sequence):
    # Using a set to store unique symbols
    alphabet = set(sequence)
    return alphabet

if __name__ == "__main__":
    # Example sequence
    S = "ATTGGCCCCGAAT"
    result = find_alphabet(S)
    
    print("Sequence:", S)
    print("Alphabet:", result)