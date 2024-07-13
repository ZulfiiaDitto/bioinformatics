# quesiton 1 how many sequences in genome 
# answer 7 
grep -c ">" wu_0.v7.fas  

# question 2
#What was the name of the third sequence in the genome file? Give the name only, without the “>” sign.
more wu_0.v7.fas|  grep ">" | head -n 3 

# answer Chr3 

# question 3 the last sequence in genome file? 
more wu_0.v7.fas|  grep ">" | tail -n 1 
# answer - mitochondria

# question4 How many index files did the operation create?
#1 mkdir index 
#2 bowtie2-build wu_0.v7.fas index/wu_0_idx
#3 cd index 
#4 ls - 6 files 

# quesiton5 What is the 3-character extension for the index files created?
# answer - bt2

# question6 How many reads were in the original fastq file?
grep -c "/1" wu_0_A_wgs.fastq
# answer 147487

# question7 7. How many matches (alignments) were reported for the original (full-match) setting? Exclude lines in the file containing unmapped reads.
bowtie2 -p 4 -x index/wu_0_idx -U wu_0_A_wgs.fastq -S out.full.sam
# answer: 147354 (all reads) - 9636(aligned zero time) = 137718 (137719)

# question 8 
#How many matches (alignments) were reported with the local-match setting? Exclude lines in the file containing unmapped reads.
bowtie2 -p 4 -x index/wu_0_idx -U wu_0_A_wgs.fastq -S out.full.sam
# answer: 147354 - 6823 = 140531

# question9 How many reads were mapped in the scenario in Question 7?
# 93780 + 43938 = 137718

# quesiton 10 How many reads were mapped in the scenario in Question 8?
# 88935 + 51596 = 140531

# question 11
# How many reads had multiple matches in the scenario in Question 7?
# answer 43938 (29.82%) aligned >1 times

# question 12 
#How many reads had multiple matches in the scenario in Question 8? Use the format above. You can တ†nd this in the bowtie2 summary; note that by default bowtie2 only reports the best match for each read.
# answer 51596 (35.01%) aligned >1 times

# question 13 
#How many alignments contained insertions and/or deletions, in the scenario in Question 7?
cat out.full.sam | grep -v "^#" | cut -f6 | grep -c 'D' 
#1395
cat out.full.sam | grep -v '^#' | cut -f6 | grep -c  'I'    
#1429
cat out.full.sam | grep -v '^#' | cut -f6 | grep 'D' | grep 'I' | wc -l
#42
# answer: 1395 + 1429 + 42 = 2866

# question 14 14. How many alignments contained insertions and/or deletions, in the scenario in Question 8?
cat out.local.sam | grep -v '^#' | cut -f6 |grep -c 'D' 
#1453
cat out.local.sam | grep -v '^#' | cut -f6 |grep -c 'I' 
#1197
cat out.local.sam | grep -v '^#' | cut -f6 |grep 'D' | grep 'I' | wc -l 
#81
# answer : 1453 + 1197 + 81 = 2731


# question 15 How many entries were reported for Chr3?














