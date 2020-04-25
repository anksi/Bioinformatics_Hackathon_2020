{
library(R.utils)
library(tidyverse)
library(DESeq2)
library(readxl)
# library(cmapR)
}

# setwd("D:/GitHub/Bioinformatics_Hackathon_2020/taoyu_mei")


# download preprocessed RNA-Seq data from GTEx ----------------------------

download.file("https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz", 
              destfile = "./GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz")

R.utils::gunzip("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz")
geneCounts <- read_tsv("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct",
                       skip = 2)

# geneCounts <- read.gct("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
# geneCounts <- cmapR::parse_gctx("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")


# subset the count matrix by patients included ---------------------------

pheno1 <- read_excel("../simons_eda/GTEx_pancreas_liver_images_liverfat_pancreasfat.xlsx")
pheno2 <- read_excel("../simons_eda/GTEx_pancreas_liver_images_liverfat_pancreasfat_seq.xlsx")
# the first file include all patients with images thus fat percentage
# the second file only include patients with RNA-Seq data, so it's a subset of the first

# pheno1$Tissue.Sample.ID_pancreas[1] %in% colnames(geneCounts)

# SAMID <- c(pheno2$SAMPID_pancreas, pheno2$SAMPID_liver)
# SAMID <- intersect(SAMID, colnames(geneCounts))
# geneCounts_subsetted <- geneCounts[c('Name', 'Description', SAMID)]


# download.file("https://www.dropbox.com/s/0ezaxy5e56nbp34/count_matrix_target_subset.tsv.gz",
#               destfile = "count_matrix_target_subset.tsv.gz")
# R.utils::gunzip("count_matrix_target_subset.tsv.gz")
# geneCounts_subsetted <- read_tsv("count_matrix_target_subset.tsv")

geneCounts_pancreas <- geneCounts[c('Name', 'Description', 
                                    intersect(pheno2$SAMPID_pancreas, colnames(geneCounts)))]

geneCounts_liver <- geneCounts[c('Name', 'Description', 
                                 intersect(pheno2$SAMPID_liver, colnames(geneCounts)))]

# make ready for differential expression analysis -------------------------




