# Script to merge exon count bed files (made with FUCHS) of same batch to create 
# a combined index that can be used for salmon quantification.

# arg[1] <- first bed file
# arg[2] <- second bed file
# arg[3] <- output file path & name
# arg[4] <- number of cores

# input files and output file info from command line
arg = commandArgs(trailingOnly=TRUE)
if (length(arg) < 3) {
  stop("At least 2 bed files and 1 output name must be supplied.\n\n
* Example: Rscript bed_index_merger.R [1st bed file] [2nd bed file] [output name] [nr of cores (optional)]\n", call.=FALSE)
}

cat("* Loading packages ...\n")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(BiocParallel)
})

# default nr of cores if not given in command line
cores <- as.numeric(arg[4])
if (length(arg)==3) {
  cores = detectCores() - 2
}
cat( sprintf("Using %i cores.\n", cores) )

cat("* Importing ...\n")
#label <- paste("b", 1:length(arg), sep = "")
#apply(as.array(arg), 2, FUN=function(x) assign(label, x) )
b1 <- import(arg[1])
b2 <- import(arg[2])

#setwd("~/Documents/R_stuff/")
#b1 <- import("1746_D.exon_counts.slim.bed")
#b2 <- import("1746_H.exon_counts.slim.bed")

cat("* Finding mismatches btw ends of transcripts & ends of outermost exons (blocks) ...\n")
error_sel <- function(x) {
  err <- which(
    start(x) != sapply(blocks(x), FUN=function(x) min(start(x))) |
      end(x) != sapply(blocks(x), FUN=function(x) max(end(x)))
  )
  return( x[-err] )
}
b1_sel <- error_sel(b1)
b2_sel <- error_sel(b2)

# we concatenate them:
b <- c(blocks(b1_sel), blocks(b2_sel))

# filter out transcripts with negative exon widths
#length( b[any(width(ranges(b)) < 0)] ) # see how many instances were found
#b <- b[!any(width(ranges(b)) < 0)]

cat("* Finding overlapping transcripts ...\n")
o <- as.matrix(findOverlaps(b,b))
o <- o[o[,1]!=o[,2],]
o_slim <- head(o) # replace o with o_slim for tests
#b[unique(as.vector(o))] # these are all transcripts that in some way overlap
bed_no <- b[-as.vector(o)] # these are all transcripts that don't overlap


thr <- 2
# function to select transcripts for combined BED file
transcript_sel <- function(x) {
  # get differences in block sequences for each overlapping transcript pair
  d <- setdiff( b[x[1]], b[x[2]] )
  # get block range widths of the dissimilar sequences
  w <- width( ranges(d) )
  # sum up the differences for each block
  s <- sum( w )
  # filter by similarity threshold and scores
  if ( s <= thr ) {
    id <- which.max( score(b[x]) ) 
  } else { 
    id <- c(1, 2) 
  }
  return( x[id] )
}
# call function
cat("* Selecting transcripts for combined bed file ...\n")
#sel <- unlist( apply(o, 1, transcript_sel) )
sel <- unlist( bplapply( 1:nrow(o), FUN = function(i){ 
          x <- o[i,] 
          transcript_sel(x) 
          }, BPPARAM = MulticoreParam(cores, progressbar = TRUE) ) 
        )
# generate GRangesList containing selected overlapping transcripts
bed_o <- GRangesList( b[unique(sel)] )

# generate a combined list with non-overlapping transcripts
bed_comb <- c(bed_o, bed_no)
# get rid of duplicate entries
bout <- asBED(bed_comb)
bed_comb_uni <- bout[!duplicated(bout)]

cat( "* Number of transcripts found to be identical (based on threshold):",
      length(b) - length(bed_comb_uni), "\n" )

# export BED file
export( bed_comb_uni, arg[3] )
cat("* Finished.\n")
