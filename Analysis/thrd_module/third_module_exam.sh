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




