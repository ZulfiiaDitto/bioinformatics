from Bio import SeqIO

from collections import defaultdict

dna_reverse_complement = {'A':'T', 
                          'T': 'A', 
                          'G': 'C', 
                          'C' : 'G'}
# question 1
filename = "dna2.fasta"
def Identify_len_file(filename):
    """Parse the fasta files and return dct"""
    dct = {}
    with open(filename, 'r') as file:
        for record in SeqIO.parse(file, "fasta"):
            dct[record.id] = record.seq
    result = len(dct)
    print(f"There {result} in a file")
    return dct

# question 2
def max_min_seq(dct):
    """input dct of parsed fasta file"""
    dctLen = dct.copy()
    for k,v in dctLen.items():
        dctLen[k]= len(v)

    sortedDCT = {k: v for k,v in sorted(dctLen.items(),
                                    key=lambda item: item[1])}

    print(f"the min sequence lengh {list(sortedDCT.values())[0]},\
          id of the sequence {list(sortedDCT.keys())[0]}") 
    print(f"the max sequence lengh {list(sortedDCT.values())[-1]},\
          the id of the sequence {list(sortedDCT.keys())[-1]}")
    return sortedDCT

# answers 
dct = Identify_len_file(filename)
print(max_min_seq(dct))


# question 3 
def reverse_complement(seq: str) -> str:
    """ Takes the string DNA and return the reversed dna string 
     by replacing the A to T and C to G and reversing the string after all """
    return  ''.join([dna_reverse_complement[nuc] for nuc in seq])[::-1]


def open_frame_coordinates(seq, keys, position):
    start_codon = ["ATG"]
    stop_codon = ['TAA', 'TAG', 'TGA']
    frame = []
    seq_pos = seq[position:]
    start, end = [], []
    lookup = []
    flag  = 0 
    for i in range(0, len(seq_pos)-3, 3):
        if seq_pos[i:i+3] == "ATG":
            start.append(i)

    for i in range(0, len(seq_pos)-3, 3):
            if seq_pos[i:i+3] in stop_codon:
                end.append(i+3)
               # print(end)
    for i in start:
        for j in end:
            if i < j and i > flag:
                seq_sub = seq_pos[i:j]
                flag = j
                lookup.append((i, j, position, len(seq_sub))) # collects coordinates, position and sequence lenght
                frame.append(seq_sub) # collects actual sequence 
    return (frame, lookup)

### answers 
seq = dct['gi|142022655|gb|EQ086233.1|16']
result_seq = open_frame_coordinates(seq, 'gi|142022655|gb|EQ086233.1|16', 2)[1]
print("\nresults for the sequence gi|142022655|gb|EQ086233.1|16 ",result_seq)



def combine_seq_orf(dct):
    reverse_dct = {}
    for k,v in dct.items():
        reverse_dct[k] = reverse_complement(v)
    orf = defaultdict(list)
    index_len = defaultdict(list)
    for k,v in dct.items():
        orf[k].extend(open_frame_coordinates(v, k, position=0)[0])
        orf[k].extend(open_frame_coordinates(v,k, position=1)[0])
        orf[k].extend(open_frame_coordinates(v,k, position=2)[0])
        #orf[k].extend(open_frame_coordinates(reverse_dct[k],k, position=0)[0])
        # orf[k].extend(open_frame_coordinates(reverse_dct[k],k, position=1)[0])
        # orf[k].extend(open_frame_coordinates(reverse_dct[k],k, position=2)[0])

        index_len[k].extend(open_frame_coordinates(v, k, position=0)[1])
        index_len[k].extend(open_frame_coordinates(v, k, position=1)[1])
        index_len[k].extend(open_frame_coordinates(v, k, position=2)[1])
        # index_len[k].extend(open_frame_coordinates(reverse_dct[k],k, position=0)[1])
        # index_len[k].extend(open_frame_coordinates(reverse_dct[k],k, position=1)[1])
        # index_len[k].extend(open_frame_coordinates(reverse_dct[k],k, position=2)[1])
    return orf, index_len

orf_dct = combine_seq_orf(dct)[0]
orf_len = combine_seq_orf(dct)[1]


def max_seq_in_orf(dct_len, position = None):
    if position:
        max_len_in_position = max([i[j][3] for i in list(dct_len.values())\
                                    for j in range(len(i)) \
                                    if i[j][2] == position])
        return (f"max orf length in fasta file in position {position} is {max_len_in_position}")
    else:
        max_in_dct = max([i[j][3] for i in list(dct_len.values()) for j in range(len(i))])
        return (f"max orf len in fasta file is {max_in_dct}")

#position = 1 -> is python coordinates: starts with 0 = 1, 1=2, 2=3 
print("Answer for question 4 is ", max_seq_in_orf(orf_len, 1), "\n") # answer - 1458

print("Answer for question 5")
for i in list(orf_len.values()):
    for j in range(len(i)):
        if i[j][2]==2:
            print(i[j])
    break

result = [i for i in list(orf_len.values()) for j in range(len(i)) if i[j][2] == 2]
#print("here",result)


position = 0
maxima = max_seq_in_orf(orf_len)#, position)
print("\nThe answer for question 6 is",maxima) # 2394

# question 4
# step1 reading whole fasta file again and append sequences to get the one huge lst of sequence
# step 2 read all repeats in big list 
# step 3 count the repeats 

def read_fast_manual(filename):
    sequences = []
    seq = ""
    with open(filename, "r") as file:
        line = file.readlines()
        for f in line:
            #print(f)
            if not f.startswith('>'):
                f = f.replace("\n", "")
                seq += f
                #print(seq)
            else:
                sequences.append(seq)
                seq = ""
    return sequences 

fasta_manual = read_fast_manual(filename=filename)

def repeats_in_seq(seq, n):
    """making all substring from string"""
    n = n*2
    repeats = list()
    length = len(seq)
    for k in range(0, length-n, n):
        substr = seq[k:k+n]
        repeats.append(substr)
    return repeats


from collections import Counter
def find_repeats_in_fasta(lst_seq, n , most_common_number = None):
    n = n*2
    repeats = []
    for seq in lst_seq:
        repeats.extend(repeats_in_seq(seq, n))
    result = Counter(repeats)
    if most_common_number:
        return result.most_common(most_common_number)
    return result

#### answers 
six_seq = find_repeats_in_fasta(fasta_manual, 6)
print(f"\nthe most common seq of len 6 occurs - {six_seq.most_common(1)} times\n")
#print(f"{six_seq}")


#print(six_seq)


seven_seq = find_repeats_in_fasta(fasta_manual, 7, 5)
print(f"the most common seq of len 7 , {seven_seq}")
print("\n")

tw = find_repeats_in_fasta(fasta_manual, 12)
max_value = tw.most_common(1)[0][1] 
num_seq = {k:v for k,v in tw.items() if v == max_value}
print(f"max of ccurence of lenght 12 is {max_value}, it occure in {len(num_seq)} seq-s")

from collections import defaultdict
from Bio import SeqIO

def find_repeats_in_fasta(file_path, repeat_length):
    repeat_counts = defaultdict(int)

    for record in SeqIO.parse(file_path, "fasta"):
        sequence = str(record.seq)
        for i in range(len(sequence) - repeat_length + 1):
            repeat = sequence[i:i + repeat_length]
            repeat_counts[repeat] += 1

    most_frequent_repeat = max(repeat_counts, key=repeat_counts.get)
    most_frequent_count = repeat_counts[most_frequent_repeat]

    return repeat_counts, most_frequent_repeat, most_frequent_count

# Example usage
file_path = filename
repeat_length = 6
repeat_counts, most_frequent_repeat, most_frequent_count = find_repeats_in_fasta(file_path, repeat_length)
#print(f"Repeat counts: {repeat_counts}")
print(f"Answer for question 8 is Most frequent repeat: {most_frequent_repeat} occurs {most_frequent_count} times")

repeat_length = 12
repeat_counts, most_frequent_repeat, most_frequent_count = find_repeats_in_fasta(file_path, repeat_length)
most_frequent_count = max(repeat_counts.values())
most_frequent_repeats = [repeat for repeat, count in repeat_counts.items() if count == most_frequent_count]

print(f"Answer for question 9 is Most frequent repeat: {most_frequent_repeat} occurs {most_frequent_count} times")
#print(most_frequent_repeats)
print(f"Number of different 12-base sequences that occur {most_frequent_count} times: {len(most_frequent_repeats)}")
print(f"Most frequent 12-base sequences: {most_frequent_repeats}")






