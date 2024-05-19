from dna_toolkit import *
import random
from bio_structure import * 
from utilities import colored 

# testing enviroment to make sure fucntions are working  

randDNA = ''.join([random.choice(nucleotides) for nuc in range(20)])
dna = validateSeq(randDNA)
print(f'Sequence: {colored(dna)}\n')
print(f'[1]. Sequence lenght {len(dna)}\n')
print(colored(f'[2]. Nucleotide frequeny : {countNucFrequency(dna)}\n'))
print(f'[3]. DNA/RNA transcriptions: {colored(transcription(dna))}\n')
print(f'[4]. DNA string + reverse Complement\n5` {colored(reverse_complement(dna))} 3`')
print(f'   {"".join(["|" for _ in range(len(dna))])}')
print(f'3` {colored(reverse_complement(dna))} 5`\n')

print(colored(dna))