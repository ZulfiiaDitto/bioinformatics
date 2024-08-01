# question 1 
library(AnnotationHub)
ah = AnnotationHub()
# pull only the humans 
homo  = query(ah, 'Homo sapiens')
# pull the Islands 
result = query(homo, "CpG Islands")
unique(result$genome)
# autosomes, pull the first one 
pull_autosome = result[['AH5086']] #GRanges object with 28691 ranges and 1 metadata column
filter = c(paste('chr', 1:22, sep = ''))

split_result = split(pull_autosome, seqnames(pull_autosome))
autosomes <- split_result[filter]
unlist(autosomes)

# question 2 
library(GenomicRanges)
autosomes[4] # cuz it has levels so level 4 -> will be chr 4 
# GRanges object with 1031 ranges and 1 metadata column:
unique(autosomes@unlistData@seqnames@values)
autosomes[12]
 
# question3 

epiFiles <- query(ah, "EpigenomeRoadMap")
ah_h3 = query(epiFiles, "H3K4me3")

ah_h3_record <- ah_h3[["AH29884"]]
split_record = split(ah_h3_record, seqnames(ah_h3_record))
record <- split_record[filter]
values_ah_h3 = unlist(record)
sum(width(unlist(record))) # 41135164

# question 4 

ah_h3k4 = query(epiFiles, "H3K27me3")
ah_h3k4
ah_h3k4_record <- ah_h3k4[["AH29892"]]
split_record_ah_h3k4 = split(ah_h3k4_record, seqnames(ah_h3k4_record))
record_ah_h3k4 <- split_record_ah_h3k4[filter]
record_ah_h3k4_sub = record_ah_h3k4[seqnames(record_ah_h3k4) %in% filter]
#record_ah_h3k4_auto = subset(record_ah_h3k4 , seqnames %in% filter)
values_h3k4 = unlist(record_ah_h3k4_sub)

mean(values_h3k4$signalValue) # 4.770728

# question5 

bival = intersect(values_ah_h3, values_h3k4)
bival

sum(width(bival))  # 10289096

# question 6 

bival_autosome = findOverlaps(bival, unlist(autosomes))

length(unique(queryHits(bival_autosome)))/length(bival) # 0.5383644

# question7 

bival_autosome_intersect = intersect(bival, unlist(autosomes))

rate_bi_auto = sum(width(bival_autosome_intersect))/sum(width(unlist(autosomes)))

rate_bi_auto
# 0.241688

# question8 

cpg_10kb <- resize(unlist(autosomes), width = 20000 + width(unlist(autosomes)), fix = "center")
cpg_10kb_bivalent <- intersect(cpg_10kb, bival)
sum(width(cpg_10kb_bivalent))

#9782086

# question9 

genome = keepSeqlevels(pull_autosome, filter, pruning.mode = 'coarse')
gene_size = sum(as.numeric(seqlengths(genome))) # 2881033286 

sum(width(unlist(autosomes)))/gene_size # 0.007047481

# question 10 

# build matrix 
overlapMat <- matrix(0,, ncol = 2, nrow = 2)
colnames(overlapMat) <- c("in", "out")
rownames(overlapMat) <- c("in", "out")

overlapMat[1,1] <- sum(width(intersect(bival, unlist(autosomes))))
overlapMat[1,2] <- sum(width(setdiff(bival, unlist(autosomes))))
overlapMat[2,1] <- sum(width(setdiff(unlist(autosomes), bival)))
overlapMat[2,2] <- gene_size - sum(overlapMat)

# odd ration 
oddsRatio <- overlapMat[1,1] * overlapMat[2,2] / (overlapMat[2,1] * overlapMat[1,2])
oddsRatio





