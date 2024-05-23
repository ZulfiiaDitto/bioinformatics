def colored(seq):
    bcolors ={
        'A' : '\033[92m', 
        'C' : '\033[94m', 
        'G' : '\033[93m',
        'T' : '\033[91m',
        'U' : '\033[91m',
        'reset' : '\033[0;0m', 
        }
    tmpseq = ''
    for nuc in seq: 
        if nuc in bcolors:
            tmpseq +=bcolors[nuc] + nuc
        else: tmpseq += bcolors['reset'] + nuc

    return tmpseq + '\033[0;0m'

def readTextFile(filePath):
    """Reads text file """
    with open(filePath, 'r') as f:
        return ''.join([l.strip() for l in f.readline()])

def writeTextFile(filePath, seq, mode = 'w'):
    """write into text file """
    with open(filePath, 'w') as f:
        f.write(seq+'\n')

def readFasta(filePath):
    """read fasta file"""
    with open(filePath, 'r') as f:
        fastafile = [l.strip() for l in f.readlines()]
    fastadct ={}
    fastalabel = ''
    for line in fastafile:
        if '>' in line:
            fastalabel = 'line'
            fastadct[fastalabel] = ''
        else:
            fastadct[fastalabel] = line
    return fastadct
