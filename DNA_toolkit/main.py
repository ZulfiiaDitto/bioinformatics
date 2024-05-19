from dna_toolkit import *
import random
# testing enviroment to make sure fucntions are working  

randDNA = ''.join([random.choice(nucleotides) for nuc in range(20)])
dna = validateSeq(randDNA)

print(countNucFrequency(dna))   