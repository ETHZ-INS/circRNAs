#############################################################################
#' salmon dataset renaming
#' 
#' function to rename Salmon datasets
#'
#' @param x List of Salmon datasets
#'
#' @return List of names, one for each dataset in list
#' 
dataset_naming_s <- function (x){
  name <-  names(x)
  sample <- lapply( names(x), FUN = function(x){ 
    if( grepl("_D_", x) ){
      gsub(".*", "young", x)
    } else if( grepl("_H_", x) ) {
      gsub(".*", "old", x)
    } else if( grepl("SRX1641435", x) ){
      gsub(".*", "s1", x)
    } else if( grepl("SRX1641436", x) ){
      gsub(".*", "s2", x)
    } else if ( grepl("SRR5715023", x) ){
      gsub(".*", "s1", x)
    } else if ( grepl("SRR5715024", x) ){
      gsub(".*", "s2", x)
    } else if ( grepl("SRR5715025", x) ){
      gsub(".*", "s3", x)
    } else if ( grepl("SRR5715026", x) ){
      gsub(".*", "s4", x)
    } else if ( grepl("SRR5715027", x) ){
      gsub(".*", "s5", x)
    }
  } )
  gentrome <- lapply( names(x), FUN = function(x){ 
    if(grepl("no_gentrome", x)){
      gsub(".*", "ng", x)
    } else {
      gsub(".*", "g", x)
    }
  } )
  valmap <- lapply( names(x), FUN = function(x){ 
    if(grepl("no_valMap", x)){
      gsub(".*", "nv", x)
    } else {
      gsub(".*", "v", x)
    }
  } )
  return( paste(sample, gentrome, valmap, sep = "_") )
}


#############################################################################
#' salmon dataset cleanup
#'
#' function to aggregate duplicates & assign transcript names to rows
#'
#' @param List of Salmon datasets 
#'
#' @return List of Salmon datasets
#' 
dataset_cleaning <- function(x){
  x <- aggregate(x[,4:5], by = list(Name = x$Name), FUN = sum)
  row.names(x) <- x$Name
  return(x)
}


#############################################################################
#' function to modify FUCHS transcript names (make them fit Salmon names)
#'
#' @param List of FUCHS datasets 
#'
#' @return List of FUCHS datasets
#' 
transcript_naming_f <- function(x){
  x$Name <- paste(x$Name, "(", x$Strand, ")", sep = "")
  return(x)
}


#############################################################################
#' salmon polya false positive analysis
#'
#' function to assess non-cero circRNA counts in PolyA data
#'
#' @param salmon polya dataset 
#' @param fuchs rnaser dataset, merged
#'
#' @return list of false positive statistics
#' 
valMap_analysis <- function(polya, fuchs){
  int <- merge( polya, fuchs, by="Name" )
  return( 
    list( "g0_nr" = sum( int$TPM!=0 ), 
          "g0_median" = median( int$TPM[int$TPM!=0] ), 
          "g0_max" = max( int$TPM[int$TPM!=0]),
          "g0_circ.ratio" = sum(int$TPM!=0)/sum(polya$TPM!=0),
          "g1_nr" = sum( int$TPM>1 ), 
          "g1_median" = median( int$TPM[int$TPM>1] ), 
          "g1_max" = max( int$TPM[int$TPM>1]),
          "g1_circ.ratio" = sum(int$TPM>1)/sum(polya$TPM>1) )
  )
}


#############################################################################
#' salmon dataset stats
#'
#' function to analyse circ vs. lin RNA count stats for dataset
#'
#' @param salmon A salmon dataset
#' @param fuchs A fuchs dataset
#'
#' @return Dataframe of count stats with no cutoffs, cutoffs at TPM=0 and TPM=1
#' 
lin_v_circ <- function(salmon, fuchs){
  ## allocation
  comb <- merge(salmon, fuchs, by="Name")
  circ <- salmon$TPM[salmon$Name %in% comb$Name]
  lin <- salmon$TPM[!salmon$Name %in% comb$Name]
  salmon.g0 <- salmon[salmon$TPM!=0,]$TPM
  circ.g0 <- salmon[salmon$TPM!=0,][salmon[salmon$TPM!=0,]$Name %in% comb$Name,]$TPM
  lin.g0 <- salmon[salmon$TPM!=0,][!salmon[salmon$TPM!=0,]$Name %in% comb$Name,]$TPM
  salmon.g1 <- salmon[salmon$TPM>1,]$TPM
  circ.g1 <- salmon[salmon$TPM>1,][salmon[salmon$TPM>1,]$Name %in% comb$Name,]$TPM
  lin.g1 <- salmon[salmon$TPM>1,][!salmon[salmon$TPM>1,]$Name %in% comb$Name,]$TPM
  ## return stats: all identified tx, circ tx, lin tx (tx nr., tx nr. ratio, median counts, max counts)
  df <- data.frame( 
    "all" = c(nrow(salmon), 1, 
              median(salmon$TPM), max(salmon$TPM)), 
    "circ" = c(length(circ), length(circ)/nrow(salmon), 
               median(circ), max(circ)), 
    "lin" = c(length(lin), length(lin)/nrow(salmon), 
              median(lin), max(lin)),
    "all.g0" = c(length(salmon.g0), 1, 
              median(salmon.g0), max(salmon.g0)), 
    "circ.g0" = c(length(circ.g0), length(circ.g0)/length(salmon.g0), 
               median(circ.g0), max(circ.g0)), 
    "lin.g0" = c(length(lin.g0), length(lin.g0)/length(salmon.g0), 
              median(lin.g0), max(lin.g0)),
    "all.g1" = c(length(salmon.g1), 1, 
                 median(salmon.g1), max(salmon.g1)), 
    "circ.g1" = c(length(circ.g1), length(circ.g1)/length(salmon.g1), 
                  median(circ.g1), max(circ.g1)), 
    "lin.g1" = c(length(lin.g1), length(lin.g1)/length(salmon.g1), 
                 median(lin.g1), max(lin.g1))
  )
  rownames(df) <- c("nc_tx", "nc_tx_ratio", "med_tpm", "max_tpm")
  return(df)
}


#############################################################################
#' selected circ stats
#'
#' function to generate selected circRNA ratios for a dataset
#'
#' @param stats List of dataset stats, generated by `salmon dataset stats` function
#' @param n_threshold Vector of threshold names
#' @param n_sample Vector of sample names
#'
#' @return Vector of circRNA fractions per sample of dataset
#' 
circ_ratio_analysis <- function(stats, n_threshold, n_sample){
  ratios <- lapply(n_threshold, FUN=function(g){
    lapply(n_sample, FUN=function(s){
      stats["nc_tx",grepl(paste0(s,".*circ.*",g), colnames(stats))] /
        stats["nc_tx",grepl(paste0(s,".*all.*",g), colnames(stats))]
    })
  })
  return( unlist(ratios) )
}
