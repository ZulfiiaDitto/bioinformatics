# installing SRA tool kit from the ncbi on mac, 
# here is link on how to do it: https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit#1-fetch-the-tar-file-from-the-canonical-location-at-ncbi

# install the samtools and bedtools on mac utilizing home brewer 

# unzip the file gunzip gencommand_proj2_data.tar.gz
# open the archive tar -xvf gencommand_proj2_data.tar


# question 1
#How many alignments does the set contain? 
samtools flagstat athal_wu_0_A.bam # first line sshould have it 
# or 
samtools view athal_wu_0_A.bam | wc -l 
# answer 221372

#Question 2
#How many alignments show the read’s mate unmapped? 
# answer 65521
samtools view athal_wu_0_A.bam | cut –f7 | grep –c "*"

# question 3 How many deletions?
# answer 2451
samtools view athal_wu_0_A.bam | cut –f6 | grep  "D" 

# Question 4
#How many alignments show the read’s mate mapped to the same chromosome?
# answer 150913
samtools view athal_wu_0_A.bam | cut –f7 | grep –c "=" # 


# question 5 How many alignments are spliced?
# answer 0
samtools view athal_wu_0_A.bam | cut –f6 | grep –c ‘N’ # 0 

# question 6 How many alignments does the set contain? 
# answer  7081
samtools sort athal_wu_0_A.bam athal_wu_0_A.sorted
samtools index athal_wu_0_A.sorted.bam
samtools view athal_wu_0_A.sorted.bam "Chr3:11777000-11794000" | wc -l 

#or 
samtools view –b athal_wu_0_A.sorted.bam “Chr3:11777000-11794000” > athal_wu_0_A.region.bam
samtools flagstat athal_wu_0_A.region.bam

# 7 How many alignments show the read’s mate unmapped? 
# answer 1983
 samtools view athal_wu_0_A.sorted.bam "Chr3:11777000-11794000" | cut -f7 | grep -c "*"


# Question 8
# How many alignments contain a deletion (D)?
# answer 31
samtools view athal_wu_0_A.region.bam | cut –f6 | grep –c ‘D’

# question 9 How many alignments show the read’s mate mapped to the same chromosome?
# answer 4670
samtools view athal_wu_0_A.bam | cut –f7 | grep –c ‘=’ 


# question 10 
# How many alignments are spliced?
# answer 0 
samtools view athal_wu_0_A.bam | cut –f6 | grep –c ‘N’ # 0


# Question 11
#How many sequences are in the genome file? 
# answer 7 
samtools view -H athal_wu_0_A.bam | grep "SN:" | wc -l

#Question 12
#What is the length of the first sequence in the genome file?
# answer 29923332 
samtools view –H athal_wu_0_A.bam | grep “SN:” | more # 29923332 - look at the LN in first line 


# quesiton 13 What alignment tool was used?
# answer stampy
samtools view –H athal_wu_0_A.bam | grep “^@PG”   

#Question 14
#What is the read identifier (name) for the first alignment?
# answer GAII05_0002:1:113:7822:3886#0
% samtools view athal_wu_0_A.bam | head -1 | cut –f1  

# Question 15
# What is the start position of this read’s mate on the genome? Give this as ‘chrom:pos’ if the read was mapped, or ‘*” if unmapped.
# answer should be Chr3:11699950 - but I got erro, so I dont know what they mean 



# quesiton 16 
#How many overlaps (each overlap is reported on one line) are reported?
# answer  3101
bedtools intersect –abam athal_wu_0_A.region.bam –b athal_wu_0_A_annot.gtf –bed -wo > overlaps.bed
wc –l overlaps.bed


# question 17 How many of these are 10 bases or longer?
# answer 2899
cut –f22 overlaps.bed | sort –nrk1 > lengths 
#or 
cut –f22 overlaps.bed | sort –nrk1 | grep –n “^9” | head -1

#Question 18
# How many alignments overlap the annotations?
# answer 3101
cut –f1-12 overlaps.bed | sort –u | wc -l . 


# Question 19
# Conversely, how many exons have reads mapped to them?
# answer 21
cut –f13-21 overlaps.bed | sort –u | wc -l 

# Question 20
# If you were to convert the transcript annotations in the file “athal_wu_0_A_annot.gtf” into BED format, how many BED records would be generated?
# answer 4 
cut –f9 athal_wu_0_A.annot.gtf | cut –d ‘ ‘ –f1,2 | sort –u | wc -l # 4 


