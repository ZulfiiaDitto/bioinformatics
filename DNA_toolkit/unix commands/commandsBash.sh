"""command need to be run in a tb1.fasta directory, eighter give full file name in comand below.
The command return the colored red Nucteotid, whcih doe snot belong to DNA nucleotides

.sh
grep -v "^>" tb1.fasta | grep --color -i "[^ATCG]"

The following command execute the commandas sequentially (;) 
and utilize && (and command) to print success 

.sh
 cat tb1.fasta; grep --color -i "[^ATCG]" tb1.fasta && echo "suscess" 

head {path to the document}. I am in working directory of the file, this command will print 5 first rows 
The bed file extension is tab delimeter file (\t)

.sh
head -n 5 Mus_musculus.GRCm38.75_chr1.bed

.sh
(head -n 2; tail -n 2) Mus_musculus.GRCm38.75_chr1.bed

count of the rows -> returns # of lines, words, characters

.sh
wc Mus_musculus.GRCm38.75_chr1.bed

to output only needed output use flag -l, -w, -c

.sh
wc -w  Mus_musculus.GRCm38.75_chr1.bed

to identify the size of the file 

.sh
ls -lh  Mus_musculus.GRCm38.75_chr1.bed

To print number of columns use awk 
.sh
tail -n +6 Mus_musculus.GRCm38.75_chr1.bed | awk -F "\t" '{print NF; exit}'

regular expression with grep, printing only 5 first rows
.sh
grep -E -o 'gene_id "\w+"' Mus_musculus.GRCm38.75_chr1.gtf | head -n 5

cleaning of the above and outpur in new file
.sh
grep -E -o 'gene_id "(\w+)"' Mus_musculus.GRCm38.75_chr1.gtf | cut -f2 -d " " | sed 's/"//g' | sort | uniq > mm_gene_id.txt

Use of bioawk library
printing the gen length 
bioawk -c gff '$3 ~ /gene/ && $2 ~ /protein_coding/ {print $seqname, $end-$start}' Mus_musculus.GRCm38.75_chr1.gtf | head -n 5

Process of fasta files using bioawk 
bioawk -c fastx '{print ">"$name"\n"$seq}' contam.fastq | head -n 4

using the revcomp() function -> printing reverse complement

bioawk -c fastx '{print ">"$name"\n"revcomp($seq)}' contam.fastq | head -n 4

Getting the summary reference form ncbi 
curl https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt > assembly_ref.txt


lets review bacilus cereus strain

cat assembly_ref.txt | grep cereus | grep 'ATCC 10876' --color=auto

To get link for in last column -> it is putput 3 links -> take first one and naviagte to the url 
cat assembly_ref.txt | grep cereus | grep 'ATCC 10876'| cut -f 20

After you navigate to url, lets get the link for gff file

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/160/895/GCF_000160895.1_ASM16089v1/GCF_000160895.1_ASM16089v1_genomic.gff.gz

ok, if you dont have wget, install homebrew. and install wget thought the folllowing comand 
brew install wget 

tring to make the change 
"" 




