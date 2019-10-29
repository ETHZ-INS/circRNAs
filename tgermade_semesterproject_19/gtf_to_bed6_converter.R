##############################################################################################
### Circtools reconstruct: Convert gene annotation GTF file into exon annotation BED6 file ###
##############################################################################################

library(GenomicRanges)
library(rtracklayer)

# import gff file (make sure you have server mounted)
## mouse
#gff <- import.gff("/reference/Mus_musculus/GENCODE/GRCm38.p5/Annotation/Genes/genes.gtf")
## rat
#gff <- import.gff("/mnt/reference/Rattus_norvegicus/Ensembl/Rnor_6.0/Annotation/Genes/genes.gtf")
## rat (using assembled gene annotation)
gff <- import.gff("/mnt/reference/Rattus_norvegicus/Ensembl/Rnor_6.0/Annotation/Genes/assembled.gtf")

# filter out exons
exon <- gff[gff$type == 'exon']
# create unique names
exon_name <- paste(exon$transcript_id, exon$exon_number, 
                   exon$gene_name, start(exon), sep = '_')
## test if the names are truly unique
## length(unique(exon_name)) == length(exon_name)

# create dataframe
exon_df <- data.frame(seqnames(exon), start(exon), end(exon), 
                 exon_name, score(exon), strand(exon))
## uncomment this depending on the naming convention used in your aligned STAR BAM files
## it gets rid of the chr prefix:
#exon_df$seqnames.exon. <- gsub("chr", "", exon_df$seqnames.exon.)
## check the seqnames:
#unique(exon_df$seqnames.exon.)
# get rid of potential NAs in score
exon_df$score.exon.[is.na(exon_df$score.exon.)] <- '.'

# Output BED6 file
## mouse
#write.table(exon_df, file='/mnt/schratt/tgermade_test/mouse_GRCm38_p5.GENCODE.exons.bed', 
#            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
## rat
#write.table(exon_df, file='/mnt/schratt/tgermade_test/rat_Rnor_6.0.Ensembl.exons.bed',
#            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
## rat (using assembled gene annotation)
write.table(exon_df, file='/mnt/schratt/tgermade_test/rat_Rnor_6.0.Ensembl.exons_assembled.bed',
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

