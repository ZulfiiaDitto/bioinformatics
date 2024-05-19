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