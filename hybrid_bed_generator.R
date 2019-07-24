cat("* Loading packages ...\n")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
})

setwd("~/Documents/R_stuff/")

# import datasets
dcc <- read.table("/mnt/ins/mouse_ribozero_hippocampus/circtools/01_detect/mouse_ribozero_hippocampus/output/CircRNACount", header = TRUE, sep = "\t")
bed <- read.table("/mnt/schratt/tgermade_test/salmon/mouse_RNAseR_hippocampus/1746_combined.exon_counts.bed", header = FALSE, sep = "\t")
colnames(bed) <- c("Chr","Start","End","Name","Counts","Strand",
                       "ExonStart","ExonEnd","Color","ExonNr","ExonLen","RelExonStart")
salmonA <- read.table("/mnt/schratt/tgermade_test/salmon/mouse_ribozero_hippocampus/transcripts_quant_SRX1641435_1746/quant.sf", header = TRUE, sep = "\t")
salmonB <- read.table("/mnt/schratt/tgermade_test/salmon/mouse_ribozero_hippocampus/transcripts_quant_SRX1641436_1746/quant.sf", header = TRUE, sep = "\t")
# modify bed names to be compatible with the salmon quant.sf names
bed$Name <- paste(bed$Name, bed$ExonLen, bed$RelExonStart, sep = ":")
bed$Name <- paste(bed$Name, "(", bed$Strand, ")", sep = "")

dcc_rnase <- read.table("/mnt/ins/mouse_RNAseR_hippocampus/circtools/01_detect/mouse_RNAseR_hippocampus/output/CircRNACount", header = TRUE, sep = "\t")

# extract transcript coordinates
dcc2 <- paste(dcc$Chr, dcc$Start, dcc$End)
bed2 <- paste(bed$Chr, bed$Start, bed$End)

dcc_rnase2 <- paste(dcc_rnase$Chr, dcc_rnase$Start, dcc_rnase$End)


# find matching transcripts between dcc and bed
m <- lapply( dcc2, FUN=function(x) which(bed2==x) )

# find matching transcripts between ribo-dcc and rnaser-dcc
m2 <- lapply( dcc2, FUN=function(x) which(dcc_rnase2==x) )
m2 <- lapply( dcc_rnase2, FUN=function(x) which(dcc2==x) )


# for all transcripts with > 1 matches, extract salmon TPM counts
m[unlist( lapply( m, FUN=function(x) length(x)>1 ) )]

for( i in m ){
  if( length(i)>1 ){
    tpmA <- lapply(i, function(x) salmonA$TPM[salmonA$Name == bed$Name[x]] )
    tpmB <- lapply(i, function(x) salmonB$TPM[salmonB$Name == bed$Name[x]] )
    print( lapply( unlist(tpmA), FUN=function(x) x/sum(unlist(tpmA)) ) )
    print( lapply( unlist(tpmB), FUN=function(x) x/sum(unlist(tpmB)) ) )
    print(lapply(i, function(x) as.character(salmonA$Name[salmonA$Name == bed$Name[x]]) ))
  } #else {print("nope")}
}

# check: are we catching the highest count RiboZero transcripts in RNAseR?
ll <- lapply(split(dcc,sapply(m2,length)),FUN=function(x) x[,4])
ll <- lapply(split(dcc,sapply(m2,length)),FUN=function(x) x[,5])
boxplot(lapply(ll[1:2],log1p))
table(sapply(m2, length))
ll2 <- lapply(split(dcc,sapply(m2,length)),FUN=function(x) rowSums(x[,-1:-3]))
boxplot(lapply(ll2[1:2],log1p))


# check: are we catching the highest count RNAseR transcripts in RiboZero?
ll <- lapply(split(dcc_rnase,sapply(m2,length)),FUN=function(x) x[,4])
ll <- lapply(split(dcc_rnase,sapply(m2,length)),FUN=function(x) x[,5])
boxplot(lapply(ll[1:2],log1p))
table(sapply(m2, length))
ll2 <- lapply(split(dcc_rnase,sapply(m2,length)),FUN=function(x) rowSums(x[,-1:-3]))
boxplot(lapply(ll2[1:2],log1p))

