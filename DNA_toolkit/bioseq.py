from bio_structure import DNA_Codons, NUCLEOTIDE_BASE, RNA_Codons
import random
from collections import Counter

class bio_seq:
    """DNA sequence class, Default value: ATCG, DNA, No label"""
    def __init__(self, seq = "ATCG", seq_type = "DNA", label = 'No label'):
       self.seq = seq
       self.label = label
       self.seq_type = seq_type
       self.is_valid = self.__validate()
       assert self.is_valid, f"Provided data does not seem to be correct {self.seq_type} sequence"

    def __validate(self): # note: two underscores before the name make the method private 
        return set(NUCLEOTIDE_BASE[self.seq_type]).issuperset(self.seq)
    
    def show_seq_info(self):
        """Return 4 strings. Full sequence information"""
        return f"[Label]: {self.label}\n[Sequence]: {self.seq}\n[Biotype]: {self.seq_type}\n[Lenght]: {len(self.seq)}"
    
    def generate_rand_seq(self, length = 10, seq_type = 'DNA'):
        """Generate a random DNA sequence, provided the length"""
        seq = ''.join([random.choice(NUCLEOTIDE_BASE[seq_type])
                       for x in range(length)])
        self.__init__(seq, seq_type, "Randomly generated sequence")


    def countNucFrequency(self) -> dict:
        """calculating the frequency of char in str"""
        tempFreqDict = Counter(self.seq)
        return dict(tempFreqDict)
    
    def transcription(self) -> str:
        """ takes string of dna and return the rna string of it. 
        Refere to transcription process in biology. 
        Note: T doe snot excist in rna, it is replaced by the U  """
        if self.seq_type == "DNA":
            return self.seq.replace('T', 'U')
        return "Not valid sequence type"
    
    def reverse_complement(self) -> str:
        """ Takes the string DNA and return the reversed dna string 
        by replacing the A to T and C to G and reversing the string after all """
        if self.seq_type == "DNA":
            mapping = str.maketrans('ATCG', 'TAGC')
            return self.seq.translate(mapping)[::-1]
        else:
            mapping = str.maketrans('AUCG', 'UAGC')
        return self.seq.translate(mapping)[::-1]

    
    def gc_content(self) -> float:
        """GC Content in dna/rna sequence"""
        return round((self.seq.count('C') + self.seq.count('G'))/len(self.seq) * 100)
    
    def gc_content_subset(self, k = 20) -> list:
        """GC content in dna/rna sub-sequence lenkgh k, k = 20 by default"""
        res = []
        for i in range(0, len(self.seq) - k+1, k):
            subseq = self.seq[i:i+k]
            res.append(round((subseq.count('C') + subseq.count('G'))/len(subseq) * 100))
        return res
    
    def translate_seq(self, init_pos = 0) -> list:
        """Translate dna seq into aminoacid seq"""
        if self.seq_type == "DNA": return [DNA_Codons[self.seq[pos:pos+3]] for pos in range(init_pos, len(self.seq)-2, 3) ]
        else: return [RNA_Codons[self.seq[pos:pos+3]] for pos in range(init_pos, len(self.seq)-2, 3)]
       
    def codon_usage(self, aminoacid : str) ->dict:
        """provide the frequency of each codon encoding a given aminoacid in a dna seq """
        tmplist = []
        if self.seq_type == "DNA":
            for i in range(0, len(self.seq) - 2, 3):
                if DNA_Codons[self.seq[i:i + 3]] == aminoacid:
                    tmplist.append(self.seq[i:i + 3])
        elif self.seq_type == "RNA":
            for i in range(0, len(self.seq) - 2, 3):
                if RNA_Codons[self.seq[i:i + 3]] == aminoacid:
                    tmplist.append(self.seq[i:i + 3])
        freqDct = dict(Counter(tmplist))
        totalWight = sum(freqDct.values())
        for seq in freqDct:
            freqDct[seq] = round(freqDct[seq]/totalWight, 2)
        return freqDct

    def gen_reading_frames(self):
        """Generating the six open reading frames of a dna seq, including reverse compliment"""
        frames = []
        frames.append(self.translate_seq(0))
        frames.append(self.translate_seq(1))
        frames.append(self.translate_seq(2))
        tmp_seq = bio_seq(self.reverse_complement(), self.seq_type)
        frames.append(tmp_seq.translate_seq(0))
        frames.append(tmp_seq.translate_seq(1))
        frames.append(tmp_seq.translate_seq(2))
        del tmp_seq
        return frames


    def proteins_from_rf(self, aa_seq):
        """ compute all posible proteins in seq nad return list of them """
        current_pr = []
        protein = []
        for aa in aa_seq:
            if aa == '_':
                # STOP accumulating if Stop codon was found
                if current_pr:
                    for p in current_pr:
                        protein.append(p)
                    current_pr = []
            else:
                # START accumulating AA if M ( it is a strat codon)
                if aa == 'M':
                    current_pr.append("")
                for i in range(len(current_pr)):
                    current_pr[i] += aa
        return protein

    def all_proteins_from_all_pr_frames(self, start =0, end = 0, ordered = False):
        if end > start:
            tmp = bio_seq(self.seq[start:end], self.seq_type)
            rfs = tmp.gen_reading_frames()
        else: rfs = self.gen_reading_frames()
        res = []
        for rf in rfs:
            prots = self.proteins_from_rf(rf)
            for p in prots:
                res.append(p)
        if ordered:
            return sorted(res, key = len, reverse=True)
        return res