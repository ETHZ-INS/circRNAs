cat("* Loading packages ...\n")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
})

setwd("~/Documents/R_stuff/")
# import fuchs datasets
D_fuchs <- read.table("1746_D.exon_counts.slim.bed", header = FALSE, sep = "\t")
colnames(D_fuchs) <- c("Chr","CircStart","CircEnd","Name","Counts","Strand",
                     "Start","End","Color","ExonNr","ExonLen","ExonStart")
H_fuchs <- read.table("1746_H.exon_counts.slim.bed", header = FALSE, sep = "\t")
colnames(H_fuchs) <- c("Chr","CircStart","CircEnd","Name","Counts","Strand",
                       "Start","End","Color","ExonNr","ExonLen","ExonStart")
# import salmon datasets
D_salmon <- read.table("1746_D_quant_slim.sf", header = TRUE, sep = "\t")
H_salmon <- read.table("1746_H_quant_slim.sf", header = TRUE, sep = "\t")

# replace simple transcript names with concatenated names for FUCHS datasets
D_fuchs$Name <- paste(D_fuchs$Name, D_fuchs$ExonLen, D_fuchs$ExonStart, sep = ":")
D_fuchs$Name <- paste(D_fuchs$Name, "(", D_fuchs$Strand, ")", sep = "")
H_fuchs$Name <- paste(H_fuchs$Name, H_fuchs$ExonLen, H_fuchs$ExonStart, sep = ":")
H_fuchs$Name <- paste(H_fuchs$Name, "(", H_fuchs$Strand, ")", sep = "")
# change indices of data.frames to name values
row.names(D_fuchs) <- D_fuchs$Name
row.names(D_salmon) <- D_salmon$Name
row.names(H_fuchs) <- H_fuchs$Name
row.names(H_salmon) <- H_salmon$Name

# get Salmon v FUCHS intersect for both samples
int_D <- intersect(D_salmon$Name, D_fuchs$Name)
int_H <- intersect(H_salmon$Name, H_fuchs$Name)

# get counts
counts_D_fuchs <- data.frame( "Counts" = D_fuchs[int_D, "Counts"] )
counts_D_salmon <- data.frame( "TPM" = D_salmon[int_D, "TPM"] )
counts_H_fuchs <- data.frame( "Counts" = H_fuchs[int_H, "Counts"] )
counts_H_salmon <- data.frame( "TPM" = H_salmon[int_H, "TPM"] )

rownames(counts_D_fuchs) <- int_D
rownames(counts_D_salmon) <- int_D
rownames(counts_H_fuchs) <- int_H
rownames(counts_H_salmon) <- int_H

# data test for H samples
cor(log1p(counts_H_fuchs), log1p(counts_H_salmon))

# transcripts common to all datasets
common <- intersect(int_D, int_H)

# store all unique Names; this determines the length of the resulting table
common_id_H <- !is.na(match(int_H,common))
int_uni <- c(int_D, int_H[!common_id_H])

# create empty dataframe
df_1746 <- data.frame("D_fuchs" = rep_len(NA, length.out = length(int_uni)), 
                      "D_salmon" = rep_len(NA, length.out = length(int_uni)), 
                      "H_fuchs" = rep_len(NA, length.out = length(int_uni)),
                      "H_salmon" = rep_len(NA, length.out = length(int_uni))) 
rownames(df_1746) <- int_uni

# function to assign counts to empty dataframe into correct row
assignment <- function(counts){
  # get indices of counts transcript names inside dataframe
  id <- match(rownames(counts), rownames(df_1746))
  # generate empty array spanning length of dataframe
  n <- array(rep_len(NA, length.out = nrow(df_1746)))
  # fill up n with counts at correct positions
  for (i in 1:length(id)){
    n[id[i]] <- counts[i,]
  }
  return(n)
}

# assignment
df_1746$D_fuchs <- assignment(counts_D_fuchs)
df_1746$D_salmon <- assignment(counts_D_salmon)
df_1746$H_fuchs <- assignment(counts_H_fuchs)
df_1746$H_salmon <- assignment(counts_H_salmon)

# export
write.table(df_1746, "1746_fuchs_salmon.txt", sep="\t", quote = FALSE) 
