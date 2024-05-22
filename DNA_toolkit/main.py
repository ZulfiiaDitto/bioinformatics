from dna_toolkit import *
import random
from bio_structure import * 
from utilities import colored 

# testing enviroment to make sure fucntions are working  

randDNA = ''.join([random.choice(nucleotides) for nuc in range(100)])
dna = validateSeq(randDNA)

# print(f'Sequence: {colored(dna)}\n')
# print(f'[1]. Sequence lenght {len(dna)}\n')
# print(colored(f'[2]. Nucleotide frequeny : {countNucFrequency(dna)}\n'))
# print(f'[3]. DNA/RNA transcriptions: {colored(transcription(dna))}\n')
# print(f'[4]. DNA string + reverse Complement\n5` {colored((dna))} 3`')
# print(f'   {"".join(["|" for _ in range(len(dna))])}')
# print(f'3` {colored(reverse_complement(dna)[::-1])} 5` [Complement] \n ')
# print(f"5` {colored(reverse_complement(dna))} 3` [Reversed Complement] \n")
# print(f'[5]. GC content : {gc_content(dna)}')
# print(f'[6]. GC content in subsequence k = 5: {gc_content_subset(dna, 5)} \n ')
# print(f'[7].  Aminoacid seq from DNA : {translate_seq(dna, 0)} \n')

# print(f'[8]. Codon frequency (L): {codon_usage(dna, "L")} \n')
# print(f'[9]. reading frames :')
# for i in gen_reading_frames(dna):
#     print(i)

# original protein from the database 
result = 'MYPESTTGSPARLSLRQTGSPGMIYRNLGKSGLRVSCLGLGTWVTFGGQITDEMAEQLMTLAYDNGINLFDTAEVYAAGKAEVVLGNIIKKKGWRRSSLVITTKIFWGGKAETERGLSRKHIIEGLKASLERLQLEYVDVVFANRPDPNTPMEETVRAMTHVINQGMAMYWGTSRWSSMEIMEAYSVARQFNLTPPICEQAEYHMFQREKVEVQLPELFHKIGVGAMTWSPLACGIVSGKYDSGIPPYSRASLKGYQWLKDKILSEEGRRQQAKLKELQAIAERLGCTLPQLAIAWCLRNEGVSSVLLGASNADQLMENIGAIQVLPKLSSSIIHEIDSILGNKPYSKKDYRS'
print(f'[10]. All prots in 6 open reading frames:')
for prot in all_proteins_from_all_pr_frames(NM_172130, 0,0, True):
    if prot == result :
        print(f'{prot}')

