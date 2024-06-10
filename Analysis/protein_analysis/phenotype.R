library(seqinr)
df <- read.table('phenotypes.txt')
df[1:10, 1:5]
# correlations
correlations <- cor(df)
correlations[1:10, 1:5]
# re calculation of correlations 

correlations <- cor(df, use ='pair')
correlations[1:10, 1:5]

# pulling out the high score correlation 

highcor <- NULL
for(r in rownames(correlations)) {
  for(c in colnames(correlations)) {
    if(correlations[r,c] > 0.7) {
      connection <- c(r , "to", c)
      highcor <- rbind(highcor, connection)
    }
  }
}

highcor[1:10]
