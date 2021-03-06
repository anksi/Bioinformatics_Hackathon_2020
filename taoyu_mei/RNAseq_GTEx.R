{
library(R.utils)
library(tidyverse)
library(DESeq2)
library(readxl)
# library(cmapR)
library(stringr)
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

# make 2 colData for the 2 organs; add a group variable according to age group

# pancreas
colDataPancreas <- pheno2[c('SAMPID_pancreas', 'Sex', 'Age.Bracket', 'Hardy.Scale', 
                            'Fat,Percentage_liver', 'Fat,Percentage_pancreas')]
colDataPancreas <- filter(colDataPancreas, 
                          SAMPID_pancreas %in% colnames(geneCounts_pancreas))
colDataPancreas <- unique.data.frame(colDataPancreas)  # remove duplicate

colDataPancreas$Age.Bracket <- gsub("^[60|70].*", "old", colDataPancreas$Age.Bracket)
colDataPancreas$Age.Bracket <- gsub("^[40|50].*", "mid_age", colDataPancreas$Age.Bracket)
colDataPancreas$Age.Bracket <- gsub("^[20|30].*", "young", colDataPancreas$Age.Bracket)

row_names <- colDataPancreas$SAMPID_pancreas
colDataPancreas <-as.matrix(colDataPancreas[-1])
row.names(colDataPancreas) <- row_names


# liver
colDataLiver <- pheno2[c('SAMPID_liver', 'Sex', 'Age.Bracket', 'Hardy.Scale', 
                            'Fat,Percentage_liver', 'Fat,Percentage_pancreas')]
colDataLiver <- filter(colDataLiver, 
                          SAMPID_liver %in% colnames(geneCounts_liver))
colDataLiver <- unique.data.frame(colDataLiver)

colDataLiver$Age.Bracket <- gsub("^[60|70].*", "old", colDataLiver$Age.Bracket)
colDataLiver$Age.Bracket <- gsub("^[40|50].*", "mid_age", colDataLiver$Age.Bracket)
colDataLiver$Age.Bracket <- gsub("^[20|30].*", "young", colDataLiver$Age.Bracket)

row_names <- colDataLiver$SAMPID_liver
colDataLiver <-as.matrix(colDataLiver[-1])
row.names(colDataLiver) <- row_names


# convert the count matrices' row names to fit DESeq2
countMatrixPancreas <- as.matrix(geneCounts_pancreas[c(-1, -2)])
row.names(countMatrixPancreas) <- geneCounts_pancreas$Name

countMatrixLiver <- as.matrix(geneCounts_liver[c(-1, -2)])
row.names(countMatrixLiver) <- geneCounts_liver$Name


# check the sample existence and order

# all(rownames(colDataLiver) %in% colnames(countMatrixLiver))
all(rownames(colDataLiver) == colnames(countMatrixLiver))
all(rownames(colDataPancreas) == colnames(countMatrixPancreas))

saveRDS(colDataLiver, "colDataLiver.rds")
saveRDS(countMatrixLiver, "countMatrixLiver.rds")
saveRDS(colDataPancreas, "colDataPancreas.rds")
saveRDS(countMatrixPancreas, "countMatrixPancreas.rds")


# create a dds object



# overview of the gene expression before DE analysis ----------------------
vsd <- vst(dds, blind=FALSE)

# PCA


# heatmap 


# hierarchical clustering tree 



# DE analysis -------------------------------------------------------------







