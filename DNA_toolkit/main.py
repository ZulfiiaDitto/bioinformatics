from dna_toolkit import *
import random
# testing enviroment 

randDNA = ''.join([random.choice(nucleotides) for nuc in range(20)])
dna = validateSeq(randDNA)

print(countNucFrequency(dna))   