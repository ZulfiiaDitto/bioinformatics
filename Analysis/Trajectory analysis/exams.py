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
    for i in range(position, len(seq)-3, 3):
        if seq[i:i+3] == "ATG":
            start.append(i)

    for i in range(position, len(seq)-3, 3):
            if seq[i:i+3] in stop_codon:
                end.append(i+3)
               # print(end)
    for i in start:
        for j in end:
            if i < j and i > flag:
                seq_pos = seq[i:j]
                flag = j
                lookup.append((i, j, position, len(seq_pos))) # collects coordinates, position and sequence lenght
                frame.append(seq_pos) # collects actual sequence 
    return (frame, lookup)

seq = dct['gi|142022655|gb|EQ086233.1|16']
result_seq = open_frame_coordinates(seq, 'gi|142022655|gb|EQ086233.1|16', 2)[1]
print("results for the sequence gi|142022655|gb|EQ086233.1|16 ",result_seq)



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
        orf[k].extend(open_frame_coordinates(reverse_dct[k],k, position=0)[0])
        orf[k].extend(open_frame_coordinates(reverse_dct[k],k, position=1)[0])
        orf[k].extend(open_frame_coordinates(reverse_dct[k],k, position=2)[0])

        index_len[k].extend(open_frame_coordinates(v, k, position=0)[1])
        index_len[k].extend(open_frame_coordinates(v, k, position=1)[1])
        index_len[k].extend(open_frame_coordinates(v, k, position=2)[1])
        index_len[k].extend(open_frame_coordinates(reverse_dct[k],k, position=0)[1])
        index_len[k].extend(open_frame_coordinates(reverse_dct[k],k, position=1)[1])
        index_len[k].extend(open_frame_coordinates(reverse_dct[k],k, position=2)[1])
    return orf, index_len

orf_dct = combine_seq_orf(dct)[0]
orf_len = combine_seq_orf(dct)[1]
#print(orf_len)
#print('lookup',orf_len['gi|142022655|gb|EQ086233.1|16'])
#print([i for i in list(orf_len.values()) if i and  i[2] == 1 ])

def max_seq_in_orf(dct_len, position = None):
    if position:
        max_len_in_position = max([j[3] for i in list(dct_len.values()) for j in i \
                                   for k in range(len(j)) if j[2] == position])
        print (f"max orf length in fasta file in position {position} is{max_len_in_position}")
    else:
        max_in_dct = max([j[3] for i in list(dct_len.values()) for j in i \
                          for k in range(len(j)) ])
    print (f"max orf len in fasta file is {max_in_dct}")

position = 1
maxima = max_seq_in_orf(orf_len, position)
print(maxima)

# question 4

def repeats_in_seq(seq, n):
    repeats = defaultdict(int)
    length = len(seq)
    for k in range(0, length-n, n):
        first_str = seq[k:k+n]
    #print("first",first_str)
        for i in range(k, length-n, n):
            repeat = seq[i:i+n]
        #print("second",repeat)
            if repeat == first_str:
                repeats[repeat] += 1
    result = {repeat: count for repeat, count in repeats.items() if count > 1 and repeat != ''}
    result_sorted = {k: v for k,v in sorted(result.items(),
                                    key=lambda item: item[1])}
    return result_sorted



def find_repeats_in_fasta(dct, n ):
    repeats = {}
    for k,v in dct.items():
        repeats[k] = repeats_in_seq(v, n)
    #print('repeats are here')
    return repeats

#print(repeats_in_seq(string,6))
#rr = find_repeats_in_fasta(dct, 12)
#print(rr['gi|142022655|gb|EQ086233.1|346'])


