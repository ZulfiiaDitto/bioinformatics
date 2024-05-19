def transcription(seq: str) -> str:
    """ takes string of dna and return the rna string of it. 
    Refere to transcription process in biology. 
    Note: T doe snot excist in rna, it is replaced by the U  """
    return seq.replace('T', 'U')


dna = 'TATCGAGATTTATGGTCGTCGAAATACACGGACCTGATTACCAGTCCTGGGATCGTCGGTCTGGGAACCTCATCAGGCTCCACAGAAGGGTACATGGGCTGGAACAACTATCGTTCTTTCCCTAGAGGATATGTATAAGATCTGTTCGCAGATGTTTATCATATGCCGCATGCTAGCAATGCCACGTATTATACACTGATCGGGATACCTCAGAGTTTAGCGGCCGACTGAGACGCTCCGAGCATATTGGTCACGTCGCACTACGTATACAGCGTTTCGACTGCCGCGCCTTGTTCGATCCAGAGATAAAACCCAAGACGTGCCTCCCTCATTAGGACGGGATCCGCCACCCTGCAAAGGATACAATCCACTCTTATTCCTTGAATTCTAGTGGGGGGGACGCACCCGCGCTAGGCACACTCGTAAGACAAGCTCGGTTAACGCGTTATCTTTGCTTTAGCTTGTTTTGACCCTCTATACTACTCAACGGTAGGGTAGGCCCCATAGTTGCCGACGGACCTGATGACTGTCAGTGTAGGAGACTACAAGTCCTTCCTCTTATGGACAATACCTCACCGGGGGAGCAGTGTGATCCATGTTAGTGCCTGCTATATTAGTCGGCACGTGTTGGAACTCTACTACCGGGAGGGGTCAGTTCGTCCAAGGAGCCTTGGGTACCCCTGCCGGAGTACCCCGTCAAACCAAACGGTGTTCACCCCTTTTCAGCGCCAATACACCAGGAACGAGAGAGCGTAGCTTTCACTGAGGAGTCCATAGACGTCCAGAAGCCTGATCTAACTAGTTCCCTGCGCTTTCCAAAAGCGTCCTGCCTATTAGTCACGAGAGTGGTTAAAGCCTAGAACTAAAGCCACCCGACCCCGAGTTCTTCATGTAAAACTCGTGGTATCTGGTTAGGCCCAGTGTGAAAGCGTTTCT'
print(transcription(dna))