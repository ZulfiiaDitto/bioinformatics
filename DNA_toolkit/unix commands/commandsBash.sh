"""command need to be run in a tb1.fasta directory, eighter give full file name in comand below.
The command return the colored red Nucteotid, whcih doe snot belong to DNA nucleotides

grep -v "^>" tb1.fasta | grep --color -i "[^ATCG]"

The following command execute the commandas sequentially (;) 
and utilize && (and command) to print success 

 cat tb1.fasta; grep --color -i "[^ATCG]" tb1.fasta && echo "suscess" 

head {path to the document}. I am in working directory of the file, this command will print 5 first rows 
The bed file extension is tab delimeter file (\t)

head -n 5 Mus_musculus.GRCm38.75_chr1.bed



"""

