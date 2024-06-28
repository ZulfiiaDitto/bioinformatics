from Bio import SeqIO

from collections import defaultdict

dna_reverse_complement = {'A':'T', 
                          'T': 'A', 
                          'G': 'C', 
                          'C' : 'G'}
# question 1
filename = "dna.example.fasta"
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
#print(max_min_seq(dct))

seq = dct['gi|142022655|gb|EQ086233.1|43']

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
    start, end = None, None
    lookup = None
    for i in range(position, len(seq_pos), 3):
        #print(seq_pos[i:i+3])
        if seq_pos[i:i+3] == "ATG":
            start = i
        if seq_pos[i:i+3] in stop_codon:
            end = i+3
    if start and end:
        seq_pos = seq_pos[start:end]
        if len(seq_pos):
            lookup = (start, end, position, len(seq_pos))
            frame.append(seq_pos)
    #print(lookup)
    return frame, lookup


def combine_seq_orf(dct):
    reverse_dct = {}
    for k,v in dct.items():
        reverse_dct[k] = reverse_complement(v)
    orf = defaultdict(list)
    index_len = defaultdict()
    for k,v in dct.items():
        orf[k].extend(open_frame_coordinates(v, k, position=0)[0])
        orf[k].extend(open_frame_coordinates(v,k, position=1)[0])
        orf[k].extend(open_frame_coordinates(v,k, position=2)[0])
        orf[k].extend(open_frame_coordinates(reverse_dct[k],k, position=0)[0])
        orf[k].extend(open_frame_coordinates(reverse_dct[k],k, position=1)[0])
        orf[k].extend(open_frame_coordinates(reverse_dct[k],k, position=2)[0])

        index_len[k] =(open_frame_coordinates(v, k, position=0)[1])
        index_len[k] =(open_frame_coordinates(v, k, position=1)[1])
        index_len[k] =(open_frame_coordinates(v, k, position=2)[1])
        index_len[k] =(open_frame_coordinates(reverse_dct[k],k, position=0)[1])
        index_len[k] =(open_frame_coordinates(reverse_dct[k],k, position=1)[1])
        index_len[k] =(open_frame_coordinates(reverse_dct[k],k, position=2)[1])


    return orf, index_len

orf_dct = combine_seq_orf(dct)[0]
orf_len = combine_seq_orf(dct)[1]

def max_seq_in_orf(dct):
    max_one = {}
    for k,v in dct.items():
        max_one[k] = max(set([len(i) for i in v]))
    dct_max = {k: v for k,v in sorted(max_one.items(),
                                    key=lambda item: item[1])}
    values = list(dct_max.values())[-1]
    keys = list(dct_max.keys())[-1]
    print(f"the max orf in {keys} and len of it {values}")
    return max_one



#result = combine_seq_orf(dct)
#print(result)

maxima = max_seq_in_orf(orf_dct)
print(maxima)

print(orf_len['gi|142022655|gb|EQ086233.1|229'])

