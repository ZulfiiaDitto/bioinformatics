# Exam # 1 Module #1 Unix commands basic 
# first create the new directory mkdir <name>
# cd into new directory cd <name>
# unzip the file gunzip gencommand_proj1_data.tar.gz
# open the archive tar -xvf gencommand_proj1_data.tar

# how many chr in a file ?
grep -c '>' apple.genome # 3 

# How many unique genes?
cut -f1 apple.genes | sort -u | wc -l  #5453

# how many transcript variants?
cut -f2 apple.genes | wc -l  # 5456

# how many genes has single variant?
cut -f1  apple.genes | uniq -c | grep -c " 1 " # 5450

# how many genes has more than one variant?
cut -f1  apple.genes | uniq -c | grep -c -v " 1 " # 3 -> which make sense 5450 + 3 -> number of genes

# how many genes on "+" strand? 
grep "+" apple.genes | cut -f1 | sort -u | wc -l  #2662

# how many genes on "-" strand?
cut -f1,4 apple.genes | sort -u |grep "-"| wc -l #2791

# how many transcripts on each chromosome 
cut -f3 apple.genes | uniq -c 
# 1625 chr1
# 2059 chr2
# 1772 chr3

# how many genes in each chormosome?
# well... not that each -> no group by have to calculate individually 
cut -f1,3 apple.genes | sort -u | grep -c "chr1"    
#1624 chr1
cut -f1,3 apple.genes | sort -u | grep -c "chr2" 
#2058 chr2
cut -f1,3 apple.genes | sort -u | grep -c "chr3" 
#1771 chr3

# How many genes are in common between condition A and condition B?
#1 prep the files and pull uniq genes from both into separate files 
cut -f1 apple.conditionA | sort -u > condA.genes 
cut -f1 apple.conditionB | sort -u > condB.genes
#2 apply comm command to define the number common genes 
comm -1 -2 condA.genes condB.genes | wc -l #2410

# How many genes are specific to condition A?
comm -2 -3 condA.genes condB.genes | wc -l #1205

# How many specific to condition B?
comm -1 -3 condA.genes condB.genes | wc -l # 1243

# how many genes are in common into 3 conditions?
# first creat eth uniq genes for condition C 
cut -f1 apple.conditionC | sort -u > condC.genes
# second export the common genes b/w condition A and B 
comm -1 -2 condA.genes condB.genes > condab.genes
#third let define all common genes b/n all conditions 
comm -1 -2 condab.genes condC.genes | wc -l #1608