library(R.utils)
library(tidyverse)
library(DESeq2)
library(readxl)
# library(CePa)
library(cmapR)

setwd("D:/GitHub/Bioinformatics_Hackathon_2020/taoyu_mei")


# download preprocessed RNA-Seq data from GTEx ----------------------------

download.file("https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz", 
              destfile = "./GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz")

R.utils::gunzip("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz", remove = FALSE)
# geneCounts <- read_tsv("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
# geneCounts <- read.gct("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
geneCounts <- parse_gctx("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")

# subset the count matrix by patients included ---------------------------

pheno1 <- read_excel("../simons_eda/GTEx_pancreas_liver_images_liverfat_pancreasfat.xlsx")
pheno2 <- read_excel("../simons_eda/GTEx_pancreas_liver_images_liverfat_pancreasfat_seq.xlsx")
