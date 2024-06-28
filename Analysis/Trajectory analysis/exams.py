from Bio import SeqIO
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

