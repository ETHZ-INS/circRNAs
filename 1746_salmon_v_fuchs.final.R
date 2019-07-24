cat("* Loading packages ...\n")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
})

cat("* Importing files ...\n")
setwd("/mnt/ins/mouse_RNAseR_hippocampus/circtools/03_reconstruct/")
# import fuchs datasets
D_fuchs <- read.table("./1746_D/1746_D.exon_counts.bed", header = FALSE, sep = "\t")
colnames(D_fuchs) <- c("Chr","CircStart","CircEnd","Name","Counts","Strand",
                     "Start","End","Color","ExonNr","ExonLen","ExonStart")
H_fuchs <- read.table("./1746_H/1746_H.exon_counts.bed", header = FALSE, sep = "\t")
colnames(H_fuchs) <- c("Chr","CircStart","CircEnd","Name","Counts","Strand",
                       "Start","End","Color","ExonNr","ExonLen","ExonStart")

setwd("/mnt/schratt/tgermade_test/salmon/mouse_RNAseR_hippocampus/")
# import salmon datasets
D_salmon <- read.table("./transcripts_quant_1746_D/quant.sf", header = TRUE, sep = "\t")
H_salmon <- read.table("./transcripts_quant_1746_H/quant.sf", header = TRUE, sep = "\t")

#cat("* Constructing table ...\n")
# replace simple transcript names with concatenated names for FUCHS datasets
D_fuchs$Name <- paste(D_fuchs$Name, D_fuchs$ExonLen, D_fuchs$ExonStart, sep = ":")
D_fuchs$Name <- paste(D_fuchs$Name, "(", D_fuchs$Strand, ")", sep = "")
H_fuchs$Name <- paste(H_fuchs$Name, H_fuchs$ExonLen, H_fuchs$ExonStart, sep = ":")
H_fuchs$Name <- paste(H_fuchs$Name, "(", H_fuchs$Strand, ")", sep = "")

cat("** Testing for & removing duplicates ...\n")
# There are duplicates in Name feature! How much do their counts differ?
D_fuchs$Counts[duplicated(D_fuchs$Name) | duplicated(D_fuchs$Name, fromLast = TRUE)]
D_salmon$TPM[duplicated(D_salmon$Name) | duplicated(D_salmon$Name, fromLast = TRUE)]
H_fuchs$Counts[duplicated(H_fuchs$Name) | duplicated(H_fuchs$Name, fromLast = TRUE)]
H_salmon$TPM[duplicated(H_salmon$Name) | duplicated(H_salmon$Name, fromLast = TRUE)]

# 3 out of 4 samples have duplicates; for all 3 its the very same transcript
dup_D_fuchs <- D_fuchs$Name[duplicated(D_fuchs$Name)]
dup_D_salmon <- D_salmon$Name[duplicated(D_salmon$Name)]
dup_H_fuchs <- H_fuchs$Name[duplicated(H_fuchs$Name)]
dup_H_salmon <- H_salmon$Name[duplicated(H_salmon$Name)]
all(sapply(list(dup_D_fuchs, dup_D_salmon, dup_H_salmon), 
           function(x) x == dup_H_salmon))

# Drop duplicates
D_fuchs <- D_fuchs[!duplicated(D_fuchs$Name),]
D_salmon <- D_salmon[!duplicated(D_salmon$Name),]
H_salmon <- H_salmon[!duplicated(H_salmon$Name),]

cat("* Constructing table ...\n")
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
cor(log1p(counts_D_fuchs), log1p(counts_D_salmon))

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

cat("* Exporting ...\n")
# export
write.table(df_1746, "/mnt/schratt/tgermade_test/1746_fuchs_salmon.txt", sep="\t", quote = FALSE) 

cat("* Done\n")
