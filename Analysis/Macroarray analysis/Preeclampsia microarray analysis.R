# analysis used publically available data from EMBL-EBI
# R analysis happend on Jun 11 2024

# getting data from the the database
# data pubmed ID 16860862  - pre-eclampsia data 
# Includes 47,000 genes screening, two chanel microarray 
# expression of the genes from 
# placenta of 13 normal and 14 pathological pregnancy, 

# 1. read the table with ref data 
# downloaded files have isue and need to be pre-processed: 
# run the following bash command, whitch replaced "#" to the tab delimeter and output cleanRefFile.txt 
#  sed 's/#/\t/g' "E-GEOD-4707.sdrf.txt" > 'CleanRefFile.txt' and then 
#  sed 's/ftp:/https:/g' CleanRefFile.txt > Reference.txt


df <- read.table('Reference.txt', sep="\t", header = TRUE)

url <- df$Comment..Derived.ArrayExpress.FTP.file.[1]

print(url)
# 2. get data from the ftp.file

temp <- tempfile()
download.file(url,temp)
con <- unz(temp, "E-GEOD-4707-processed-data-1618157355.txt")

data <- read.table(con, header = TRUE,  sep="\t")
unlink(temp)

df[c('Characteristics..DiseaseState.' , 'Scan.Name')]



