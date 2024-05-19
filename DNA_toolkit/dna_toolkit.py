from collections import Counter
nucleotides = ['A', 'C', 'G', 'T']

# creating functions for  dna string validation and freq calculation 

def validateSeq(dna_seq: str):
    """Takes dna string, and check if it is valid, menaing in nucleotide sequience"""
    tempseq = dna_seq.upper()
    for nuc in tempseq:
        if nuc not in nucleotides:
            return False 
    return tempseq

def countNucFrequency(seq: str) -> dict:
    """calculating the frequency of char in str"""
    tempFreqDict = Counter(seq)
    return dict(tempFreqDict)

