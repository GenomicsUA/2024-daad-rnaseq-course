#!/usr/bin/Rscript

# adapted from https://github.com/dpryan79/Answers/tree/master/SEQanswers_42420

library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)
library(tidyverse)

#GTFfile = "../../data/Mus_musculus_c57bl6nj.C57BL_6NJ_v1.112.gtf"
#FASTAfile = "../../data/Mus_musculus_c57bl6nj.C57BL_6NJ_v1.dna.toplevel.fa"

GTFfile = "../../data/Mus_musculus.GRCm39.112.gtf"
FASTAfile = "../../data/Mus_musculus.GRCm39.dna.toplevel.fa"

# Load the annotation and reduce it
# asRangedData = F
GTF <- import.gff(GTFfile, format="gtf", genome="GRCm38", feature.type="exon")
grl <- reduce(split(GTF, elementMetadata(GTF)$gene_id))
reducedGTF <- unlist(grl, use.names=T)
elementMetadata(reducedGTF)$gene_id <- rep(names(grl), elementNROWS(grl))

#Open the fasta file
FASTA <- FaFile(FASTAfile)
open(FASTA)

#Add the GC numbers
elementMetadata(reducedGTF)$nGCs <- letterFrequency(getSeq(FASTA, reducedGTF), "GC")[,1]
elementMetadata(reducedGTF)$widths <- width(reducedGTF)

#Create a list of the ensembl_id/GC/length
calc_GC_length <- function(x) {
    nGCs = sum(elementMetadata(x)$nGCs)
    width = sum(elementMetadata(x)$widths)
    c(width, nGCs/width)
}
# takes time ~ 10 min
output <- t(sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), calc_GC_length))
df_output <- output %>% as.data.frame() %>% rownames_to_column("ensembl_gene_id")

colnames(df_output) <- c("ensembl_gene_id", "Length", "GC")

write_tsv(df_output, "../../99_technical/GC_lengths.tsv")
