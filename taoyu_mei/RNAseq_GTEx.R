library(R.utils)
library(tidyverse)
library(DESeq2)


# download preprocessed RNA-Seq data from GTEx

download.file("https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz", 
              destfile = "./GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz")

R.utils::gunzip("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz", remove = FALSE)
# geneCounts <- read_tsv(gzfile("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"))

