
# read data form FASTA file

def readfile(filepath):
    """Reading file and returning list of lines"""
    with open(filepath, 'r') as f:
        return [l.strip() for l in f.readlines()]
    
def gc_content(seq):
    """GC Content in dna/rna sequence"""
    return round(((seq.count('C') + seq.count('G'))/len(seq) * 100), 6)

#1. sttoring file content in a list 
fastaFile = readfile('rosalind_problem/test_data/gc_content.txt')
#2 dct for data and label 
fastadct = {}
#3 string for holding currentl label 
fastalabel = ''

for line in fastaFile:
    if '>' in line:
        fastalabel = line
        fastadct[fastalabel] = ''
    else: fastadct[fastalabel] += line

fastaGC_content = { key: gc_content(value) for (key,value )in fastadct.items()}

# used sorted dct for to return the first key and values 
maxCG_content = sorted(fastaGC_content.items(), key=lambda item: item[1], reverse=True)
print(maxCG_content)

maxCG_contentValue = max(fastaGC_content, key = fastaGC_content.get)
result = f'{maxCG_contentValue[1:]}\n{fastaGC_content[maxCG_contentValue]}'

print(result)
