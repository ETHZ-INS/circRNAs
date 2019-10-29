## ---- packages
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(LSD)
  library(ggplot2)
  library(cowplot)
  library(ggpubr)
})


###################################### ---- functions
source("/home/tomas/Documents/Schratt_Lab/schratt_rotation_functions.R")


###################################### ---- import circbase
# CircBase import
circbase <- read.table("/home/tomas/Documents/Schratt_Lab/circbase_mouse_mm10.bed", 
                       header = TRUE, sep = "\t")
colnames(circbase) <- c("Chr","Start","End","Name","Counts","Strand",
                        "ExonStart","ExonEnd","Color","ExonNr","ExonLen","RelExonStart")

setwd("/home/tomas/Documents/Schratt_Lab")
save(circbase, file="circbase.RData")

###################################### ---- import salmon rnaser
# Salmon RNAseR import & tweaking

## names of all directories containing Salmon quantifications
setwd("/mnt/schratt/tgermade_test/salmon/")
s_dirs <- list.files(pattern="quant.sf", recursive=TRUE, full.names=TRUE)

## create lists of datasets, separated by Our/Dieterich
list_s_rnaser_our <- lapply(s_dirs[grepl("our/",s_dirs)], header=TRUE, 
                            sep="\t", FUN=read.table)
list_s_rnaser_diet <- lapply(s_dirs[grepl("dieterich/",s_dirs)], header=TRUE, 
                             sep="\t", FUN=read.table)

## rename datasets
names(list_s_rnaser_our) <- dirname(s_dirs[grepl("our/",s_dirs)])
names(list_s_rnaser_our) <- dataset_naming_s(list_s_rnaser_our)

names(list_s_rnaser_diet) <- dirname(s_dirs[grepl("dieterich/",s_dirs)])
names(list_s_rnaser_diet) <- dataset_naming_s(list_s_rnaser_diet)

## clean datasets
list_s_rnaser_our <- lapply( list_s_rnaser_our, dataset_cleaning )
list_s_rnaser_diet <- lapply( list_s_rnaser_diet, dataset_cleaning )

setwd("/home/tomas/Documents/Schratt_Lab")
save(list_s_rnaser_diet, list_s_rnaser_our, file="list_s_rnaser.RData")


###################################### ---- import salmon polya
# Salmon PolyA import & tweaking

## names of all directories containing Salmon quantifications
setwd("/mnt/schratt/tgermade_test/salmon/mouse_hippocampus_polya_bohacek/")
dirs_polya <- list.files(pattern="quant.sf", recursive=TRUE, full.names=TRUE)
list_s_polya <- lapply(dirs_polya, header=TRUE, sep="\t", FUN=read.table)

## rename dataset
names(list_s_polya) <- dirname(dirs_polya)
names(list_s_polya) <- dataset_naming_s(list_s_polya)

setwd("/home/tomas/Documents/Schratt_Lab")
save(list_s_polya, file="list_s_polya.RData")


###################################### ---- import salmon ribozero
# Salmon RiboZero import & tweaking

setwd("/mnt/schratt/tgermade_test/salmon/mouse_hippocampus_ribozero/")
dirs_riboz <- list.files(pattern="quant.sf", recursive=TRUE, full.names=TRUE)

## create lists of datasets, separated by Our/Dieterich
list_s_riboz_our <- lapply(dirs_riboz[grepl("our_",dirs_riboz)], header=TRUE, 
                           sep="\t", FUN=read.table)
list_s_riboz_diet <- lapply(dirs_riboz[grepl("dieterich_",dirs_riboz)], header=TRUE, 
                            sep="\t", FUN=read.table)

## rename datasets
names(list_s_riboz_our) <- dirname(dirs_riboz[grepl("our_",dirs_riboz)])
names(list_s_riboz_our) <- dataset_naming_s(list_s_riboz_our)

names(list_s_riboz_diet) <- dirname(dirs_riboz[grepl("dieterich_",dirs_riboz)])
names(list_s_riboz_diet) <- dataset_naming_s(list_s_riboz_diet)

## clean datasets
list_s_riboz_our <- lapply( list_s_riboz_our, dataset_cleaning )
list_s_riboz_diet <- lapply( list_s_riboz_diet, dataset_cleaning )

setwd("/home/tomas/Documents/Schratt_Lab")
save(list_s_riboz_diet, list_s_riboz_our, file="list_s_riboz.RData")


###################################### ---- import fuchs rnaser
# FUCHS RNAseR import & tweaking

## find all directories containing Our FUCHS quantifications & Import
setwd("/mnt/schratt/tgermade_test/mouse_rnaser_hippocampus/1746_D_H_corrected/")
f_our_dirs <- list.files(pattern="corrected.bed", recursive=TRUE, full.names=TRUE)
list_f_rnaser_our <- lapply(f_our_dirs, header=TRUE, sep="\t", FUN=read.table)

## find all directories containing Dieterich FUCHS quantifications & Import
setwd("/mnt/schratt/tgermade_test/mouse_rnaser_hippocampus_Dieterich/")
f_diet_dirs <- list.files(pattern="corrected.bed", recursive=TRUE, full.names=TRUE)
list_f_rnaser_diet <- lapply(f_diet_dirs, header=TRUE, sep="\t", FUN=read.table)

## name datasets & their columns
list_f_rnaser_our <- lapply(list_f_rnaser_our, setNames, 
                            nm = c("Chr","CircStart","CircEnd","Name","Counts","Strand",
                                   "Start","End","Color","ExonNr","ExonLen","ExonStart"))
list_f_rnaser_diet <- lapply(list_f_rnaser_diet, setNames, 
                             nm = c("Chr","CircStart","CircEnd","Name","Counts","Strand",
                                    "Start","End","Color","ExonNr","ExonLen","ExonStart"))
## rename FUCHS datasets
names(list_f_rnaser_our) <- c("young", "old")
names(list_f_rnaser_diet) <- c("young", "old")

## rename transcripts
list_f_rnaser_diet <- lapply( list_f_rnaser_diet, transcript_naming_f )
list_f_rnaser_our <- lapply( list_f_rnaser_our, transcript_naming_f )

setwd("/home/tomas/Documents/Schratt_Lab")
save(list_f_rnaser_diet, list_f_rnaser_our, file="list_f_rnaser.RData")


###################################### ---- salmon fuchs
# Salmon/FUCHS RNAseR
merge_s_f_young <- merge(list_f_rnaser_our$young, list_s_rnaser_our$young_g_v, by = "Name")
merge_s_f_old <- merge(list_f_rnaser_our$old, list_s_rnaser_our$old_g_v, by = "Name")
cor_s_f_young <- cor(log1p(merge_s_f_young$TPM), log1p(merge_s_f_young$Counts))
cor_s_f_old <- cor(log1p(merge_s_f_old$TPM), log1p(merge_s_f_old$Counts))


setwd("/home/tomas/Documents/Schratt_Lab")
save(cor_s_f_young, cor_s_f_old, file="list_f_s_cors.RData")


###################################### ---- salmon parameters lists
# lists of merged datasets
list_merged_s <- list(
  "young_v" = merge(list_s_rnaser_our$young_ng_v, list_s_rnaser_our$young_g_v, by="Name", 
                      suffix = c("_ng", "_g")),
  "young_nv" = merge(list_s_rnaser_our$young_ng_nv, list_s_rnaser_our$young_g_nv, by="Name",
                       suffix = c("_ng", "_g")),
  "old_v" = merge(list_s_rnaser_our$old_ng_v, list_s_rnaser_our$old_g_v, by="Name",
                    suffix = c("_ng", "_g")),
  "old_nv" = merge(list_s_rnaser_our$old_ng_nv, list_s_rnaser_our$old_g_nv, by="Name",
                     suffix = c("_ng", "_g")),
  "young_g" = merge(list_s_rnaser_our$young_g_nv, list_s_rnaser_our$young_g_v, by="Name", 
                      suffix = c("_nv", "_v")),
  "old_g" = merge(list_s_rnaser_our$old_g_nv, list_s_rnaser_our$old_g_v, by="Name", 
                    suffix = c("_nv", "_v"))
)
list_merged_f <- list(
  "young_g_v" = merge(list_f_rnaser_our$young, list_s_rnaser_our$young_g_v, by="Name"),
  "young_g_nv" = merge(list_f_rnaser_our$young, list_s_rnaser_our$young_g_nv, by="Name"),
  "old_g_v" = merge(list_f_rnaser_our$old, list_s_rnaser_our$old_g_v, by="Name"),
  "old_g_nv" =  merge(list_f_rnaser_our$old, list_s_rnaser_our$old_g_nv, by="Name"),
  "young_ng_v" = merge(list_f_rnaser_our$young, list_s_rnaser_our$young_ng_v, by="Name"),
  "young_ng_nv" = merge(list_f_rnaser_our$young, list_s_rnaser_our$young_ng_nv, by="Name"),
  "old_ng_v" = merge(list_f_rnaser_our$old, list_s_rnaser_our$old_ng_v, by="Name"),
  "old_ng_nv" =  merge(list_f_rnaser_our$old, list_s_rnaser_our$old_ng_nv, by="Name")
)

setwd("/home/tomas/Documents/Schratt_Lab")
save(list_merged_s, list_merged_f, file="list_merged.RData")


######## this can be added to text ########
round( mean(c(cor(log1p(list_merged_s$young_v$TPM_ng), log1p(list_merged_s$young_v$TPM_g)), 
              cor(log1p(list_merged_s$old_v$TPM_ng), log1p(list_merged_s$old_v$TPM_g)))),2 )
round( sd(c(cor(log1p(list_merged_s$young_v$TPM_ng), log1p(list_merged_s$young_v$TPM_g)), 
            cor(log1p(list_merged_s$old_v$TPM_ng), log1p(list_merged_s$old_v$TPM_g)))),5 )
round( mean(c(cor(log1p(list_merged_s$young_nv$TPM_ng), log1p(list_merged_s$young_nv$TPM_g)), 
              cor(log1p(list_merged_s$old_nv$TPM_ng), log1p(list_merged_s$old_nv$TPM_g)))),2 )
round( sd(c(cor(log1p(list_merged_s$young_nv$TPM_ng), log1p(list_merged_s$young_nv$TPM_g)), 
            cor(log1p(list_merged_s$old_nv$TPM_ng), log1p(list_merged_s$old_nv$TPM_g)))),5 )


###################################### ---- salmon polya valmap
# Salmon PolyA: valMap analysis

## union of FUCHS datasets (young + old): all identified circRNAs
comb_f <- merge(list_f_rnaser_our$young, list_f_rnaser_our$old, by="Name", all=TRUE)

## get circRNA stats for PolyA data
nonc_polya <- lapply( list_s_polya, FUN=function(x) valMap_analysis(x, comb_f) )

## create dataframe of stats
polya_circ_nr <- unlist(nonc_polya)[grepl("nr",names(unlist(nonc_polya)))]
polya_circ_med <- unlist(nonc_polya)[grepl("median",names(unlist(nonc_polya)))]
polya_circ_max <- unlist(nonc_polya)[grepl("max",names(unlist(nonc_polya)))]
polya_circ_ratio <- unlist(nonc_polya)[grepl("circ.ratio",names(unlist(nonc_polya)))]

polya_circ <- data.frame("nr" = polya_circ_nr, 
                         "median" = polya_circ_med, 
                         "max" = polya_circ_max,
                         "circ.ratio" = polya_circ_ratio)

## get mean stat values over all 5 samples, separated by young/old & v/nv & g0/g1
polya_circ_stats <- data.frame(
  "valMap.g0.mean" = colMeans( polya_circ[grepl(".*_v.g0", rownames(polya_circ)),] ),
  "valMap.g0.sd" = sapply( polya_circ[grepl(".*_v.g0", rownames(polya_circ)),], sd ),
  "no_valMap.g0.mean" = colMeans(polya_circ[grepl(".*_nv.g0", rownames(polya_circ)),]),
  "no_valMap.g0.sd" = sapply( polya_circ[grepl(".*_nv.g0", rownames(polya_circ)),], sd ),
  "valMap.g1.mean" = colMeans( polya_circ[grepl(".*_v.g1", rownames(polya_circ)),] ),
  "valMap.g1.sd" = sapply( polya_circ[grepl(".*_v.g1", rownames(polya_circ)),], sd ),
  "no_valMap.g1.mean" = colMeans(polya_circ[grepl(".*_nv.g1", rownames(polya_circ)),]),
  "no_valMap.g1.sd" = sapply( polya_circ[grepl(".*_nv.g1", rownames(polya_circ)),], sd )
)

### create dataframe for plot
polya_circ_stats_df <- data.frame(
  "opt" = c("S w vM", "S w/o vM", "S w vM", 
            "S w/o vM"),
  "thr" = c("> 0","> 0","> 1","> 1"),
  "val" = as.numeric(polya_circ_stats["circ.ratio", grepl(".mean",colnames(polya_circ_stats))]),
  "sd" = as.numeric(polya_circ_stats["circ.ratio", grepl(".sd",colnames(polya_circ_stats))])
)

## get mean stat ratios: no_valMap / valMap
polya_stats_ratio <- as.table( 
  polya_circ_stats[,grepl("no_valMap.mean",colnames(polya_circ_stats))] / 
    polya_circ_stats[,grepl("^valMap.mean",colnames(polya_circ_stats))] )
names(polya_stats_ratio) <- c("nr","median","max", "circ.ratio")

## ratio: how many FUCHS transcripts do we find in either v/nv Salmon?
sel_mean <- grepl("mean", colnames(polya_circ_stats))
polya_ratio_fuchs <- polya_circ_stats["nr",] / nrow(comb_f)
### create dataframe to plot
polya_ratio_fuchs <- data.frame(
  "opt" = c("S w vM", "S w/o vM", "S w vM", 
            "S w/o vM","FUCHS"),
  "thr" = c("0","0","1","1","0"),
  "val" = c(as.numeric(polya_ratio_fuchs[grepl(".mean",colnames(polya_ratio_fuchs))]),1),
  "sd" = c(as.numeric(polya_ratio_fuchs[grepl(".sd",colnames(polya_ratio_fuchs))]), NA)
)

## p value for difference btw number of detected circRNA in PolyA by Salmon
polya_v.g0_ttest <- t.test(polya_circ_nr[grepl("_v.g0", names(polya_circ_nr))], 
       polya_circ_nr[grepl("_nv.g0", names(polya_circ_nr))])$p.value
polya_v.g1_ttest <- t.test(polya_circ_nr[grepl("_v.g1", names(polya_circ_nr))], 
                           polya_circ_nr[grepl("_nv.g1", names(polya_circ_nr))])$p.value


setwd("/home/tomas/Documents/Schratt_Lab")
save(comb_f, polya_circ_stats, polya_stats_ratio, polya_ratio_fuchs, polya_circ_stats_df,
     file="salmon_valmap_stats.RData")


###################################### ---- salmon lin vs. circ ratios
# Salmon RNAseR/PolyA: lin vs. circ tx analysis

## create lin vs circ TPMs stats
s_rnaser_tx_stats <- as.data.frame(
  lapply( list_s_rnaser_our, FUN = function(x) lin_v_circ(x, comb_f) )
)
s_polya_tx_stats <- as.data.frame(
  lapply( list_s_polya, FUN = function(x) lin_v_circ(x, comb_f) )
)
sel_s_rnaser_tx_stats <- s_rnaser_tx_stats[grepl("_g_nv", colnames(s_rnaser_tx_stats))]
sel_s_polya_tx_stats <- s_polya_tx_stats[grepl("_g_nv", colnames(s_polya_tx_stats))]

# generate circRNA ratios for RNAseR & PolyA datasets
circ_r_rnaser <- circ_ratio_analysis(sel_s_rnaser_tx_stats, 
                                     c("g0","g1"), c("young","old"))
circ_r_polya <- circ_ratio_analysis(sel_s_polya_tx_stats, 
                                    c("g0","g1"), c("s1","s2","s3","s4","s5"))

# generate dataframe of circRNA ratio values
circ_ratio_s <- data.frame(
  "val" = c(circ_r_rnaser, circ_r_polya),
  "thr" = c(rep("0",2), rep("1",2), rep("0",5), rep("1",5)),
  "set" = factor( c(rep("RNAseR",4), rep("PolyA",10)), levels=c("RNAseR","PolyA") )
)
circ_ratio_s$set_f = factor(circ_ratio_s$set, levels=c("RNAseR","PolyA"))

setwd("/home/tomas/Documents/Schratt_Lab")
save(s_rnaser_tx_stats, s_polya_tx_stats, circ_ratio_s, file="salmon_tx_stats.RData")



ggplot(circ_ratio_s, aes(group=thr, y=val, fill=factor(thr))) + 
  geom_boxplot() + facet_wrap(~set, scale="free_y") + theme_cowplot() + 
  labs(y = "circRNA fraction") +
  guides(fill=guide_legend(title="TPM threshold")) 
+
  geom_point(position=position_dodge(width=0.75),aes(group=thr))


bp_rnaser <- ggplot(circ_ratio_s[1:4,], aes(x=thr, y=val, fill=thr)) + 
  geom_boxplot() + theme_cowplot() + 
  labs(x="RNAseR", y = "circRNA fraction") +
  guides(fill=guide_legend(title="TPM threshold")) +
  geom_point(position=position_dodge(width=0.75),aes(group=thr))