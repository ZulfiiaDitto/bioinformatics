from dna_toolkit import *
import random
from  structures import * 
# testing enviroment to make sure fucntions are working  

randDNA = ''.join([random.choice(nucleotides) for nuc in range(20)])
dna = validateSeq(randDNA)
print(f'Sequence: {dna}\n')
print(f'[1]. Sequence lenght {len(dna)}\n')
print(f'[2]. Nucleotide frequeny : {countNucFrequency(dna)}\n')
print(f'[3]. DNA/RNA transcriptions: {transcription(dna)}\n')
print(f'[4]. DNA string + reverse Complement\n5` {reverse_complement(dna)} 3`')
print(f'   {"".join(["|" for _ in range(len(dna))])}')
print(f'3` {reverse_complement(dna)} 5`\n')