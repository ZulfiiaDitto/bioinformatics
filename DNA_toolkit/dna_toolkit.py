from collections import Counter
from bio_structure import *

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

def transcription(seq: str) -> str:
    """ takes string of dna and return the rna string of it. 
    Refere to transcription process in biology. 
    Note: T doe snot excist in rna, it is replaced by the U  """
    return seq.replace('T', 'U')

def reverse_complement(seq: str) -> str:
    """ Takes the string DNA and return the reversed dna string 
     by replacing the A to T and C to G and reversing the string after all """
   # return ''.join([dna_reverse_complement[nuc] for nuc in seq])[::-1]
    mapping = str.maketrans('ATCG', 'TAGC')
    return seq.translate(mapping)[::-1]

def gc_content(seq: str) -> float:
    """GC Content in dna/rna sequence"""
    return round((seq.count('C') + seq.count('G'))/len(seq) * 100)

def gc_content_subset(seq: str, k = 20 ) -> list:
    """GC content in dna/rna sub-sequence lenkgh k, k = 20 by default"""
    res = [gc_content(seq[i:i+k]) for i in range(0, len(seq) - k+1, k)]
    # for i in range(0, len(seq) - k+1, k):
    #     subseq = seq[i:i+k]
    #     res.append(gc_content(subseq))
    return res

def translate_seq(seq: str, init_pos = 0) -> list:
    """Translate dna seq into aminoacid seq"""
    return [DNA_Codons[seq[pos:pos+3]] for pos in range(init_pos, len(seq)-2, 3) ]

def codon_usage(seq : str, aminoacid : str) ->dict:
    """provide the frequency of each codon encoding a given aminoacid in a dna seq """
    #tmplist = [seq[i:i+3] for i in range(0,len(seq)-2, 3)  if DNA_Codons[seq[i:i+3] == aminoacid]]
    tmplist = []
    for i in range(0, len(seq) - 2, 3):
        if DNA_Codons[seq[i:i + 3]] == aminoacid:
            tmplist.append(seq[i:i + 3])

    freqDct = dict(Counter(tmplist))
    totalWight = sum(freqDct.values())
    for seq in freqDct:
        freqDct[seq] = round(freqDct[seq]/totalWight, 2)
    return freqDct


