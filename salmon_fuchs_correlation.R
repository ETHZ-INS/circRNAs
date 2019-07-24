cat("* Loading packages ...\n")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
})

setwd("~/Documents/R_stuff/")
#fuchs <- import("1746_D.exon_counts.slim.bed")
fuchs <- read.table("1746_D.exon_counts.slim.bed", header = FALSE, sep = "\t")
colnames(fuchs) <- c("Chr","CircStart","CircEnd","Name","Counts","Strand",
                     "Start","End","Color","ExonNr","ExonLen","ExonStart")
salmon <- read.table("1746_D_quant_slim.sf", header = TRUE, sep = "\t")


# 1)

# replace simple transcript names with concatenated names (like salmon)
fuchs$Name <- paste(fuchs$Name, fuchs$ExonLen, fuchs$ExonStart, sep = ":")
fuchs$Name <- paste(fuchs$Name, "(", fuchs$Strand, ")", sep = "")
# change indices of data.frames to name values
row.names(fuchs) <- fuchs$Name
row.names(salmon) <- salmon$Name

# get intersect of both sets
int <- intersect(salmon$Name, fuchs$Name)

# get counts
s_fuchs <- fuchs[int, "Counts"]
s_salmon <- salmon[int, "TPM"]
## data check
length(s_fuchs) == length(s_salmon)
all( fuchs[int, "Names"] == salmon[int, "Names"] )

# data tests 
plot(s_salmon, s_fuchs, abline(lm(s_fuchs ~ s_salmon)), log="xy")
lm(s_fuchs ~ s_salmon)
cor1 <- cor(log1p(s_fuchs), log1p(s_salmon))


# 2)

# change Name feature, excluding transcript names
fuchs$Name <- paste(fuchs$ExonLen, fuchs$ExonStart, sep = ":")
fuchs$Name <- paste(fuchs$Name, "(", fuchs$Strand, ")", sep = "")
salmon$Name <- substring(salmon$Name, first=20)
# change indices of data.frames to name values
row.names(fuchs) <- fuchs$Name
row.names(salmon) <- salmon$Name

# get intersect of both sets
int <- intersect(salmon$Name, fuchs$Name)

# get counts
s_fuchs <- fuchs[int, "Counts"]
s_salmon <- salmon[int, "TPM"]
## data check
length(s_fuchs) == length(s_salmon)
all( fuchs[int, "Names"] == salmon[int, "Names"] )

# data tests 
plot(s_salmon, s_fuchs, abline(lm(s_fuchs ~ s_salmon)), log="xy")
lm(s_fuchs ~ s_salmon)
cor2 <- cor(log1p(s_fuchs), log1p(s_salmon))

# are the two tests identical?
cor1==cor2

#############################################################################

# function to count occurrences of transcript names
count <- function(x, n){ length((which(x == n))) }

# select single occurrences in salmon & sort
num_salmon <- apply( as.array(1:length(salmon$Name)), 1, 
                     function(x) count(salmon$Name, salmon$Name[x]) )
uni_salmon <- salmon[num_salmon==1,]
uni_sort_salmon <- uni_salmon[order(uni_salmon$Name),]

# select single occurrences in fuchs & sort
num_fuchs <- apply( as.array(1:length(fuchs)), 1, 
                    function(x) count(fuchs$name, fuchs$name[x]) )
uni_fuchs <- fuchs[num_fuchs==1]
uni_sort_fuchs <- uni_fuchs[order(uni_fuchs$name)]

# get intersect of both
uni_int <- intersect(uni_sort_salmon$Name, uni_sort_fuchs$name)

# counts for all shared and unique transcripts
uni_s_salmon <- uni_sort_salmon$TPM[!is.na(match(uni_sort_salmon$Name, uni_int))]
uni_s_fuchs <- uni_sort_fuchs$score[!is.na(match(uni_sort_fuchs$name, uni_int))]

# test the order and occurrence of the transcripts in both modified datasets
test_salmon <- uni_sort_salmon$Name[!is.na(match(uni_sort_salmon$Name, uni_int))]
test_fuchs <- uni_sort_fuchs$name[!is.na(match(uni_sort_fuchs$name, uni_int))]
all( test_salmon == test_fuchs )

# plot
plot(uni_s_salmon, uni_s_fuchs, abline(lm(uni_s_fuchs ~ uni_s_salmon)), log="xy")
lm(uni_s_fuchs ~ uni_s_salmon)
cor(log1p(uni_s_fuchs), log1p(uni_s_salmon))
# extremely high correlation but only very few comparisons
length(uni_s_salmon)

