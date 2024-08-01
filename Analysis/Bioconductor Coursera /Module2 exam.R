library(AnnotationHub)
library(Biostrings)
library(BSgenome)
library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# question1 
freq = alphabetFrequency(Hsapiens$chr22) # initiate object 
total_len = sum(freq[c('A','C', 'G', 'T')])
gc_freq = sum(freq[c('C', 'G')])

gc_freq/total_len
# 0.4798807 

# question2 
ah = AnnotationHub()
result = query(ah, c("H3K27me3", "E003", "narrowPeak"))
result_record = result[["AH29892"]]
chr22 = subset(result_record, seqnames == 'chr22')
chr22_view = Views(Hsapiens, chr22)
chr22_gc = letterFrequency(chr22_view, 'GC', as.prob = TRUE)
mean(chr22_gc)
# 0.528866 


# question3 
signal = mcols(chr22_view)$signalValue
cor(signal, chr22_gc)
# 0.004467924 

# quesiton4 




