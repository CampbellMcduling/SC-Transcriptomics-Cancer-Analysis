### DATA PREPROCESSING ###

# read in data on expression
bc_dat <- read.table(gzfile("GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt.gz"), header = TRUE)
# read in data on sample info
bc_info <- read.table(gzfile("GSE75688_final_sample_information.txt.gz"),header = TRUE)
str(bc_info)

# remove pooled samples and only keep single cell data from sample info
sc_info <- bc_info[which(bc_info$type == "SC"),]
# remove pooled samples and only keep single cell data from expression data
sc_dat <- bc_dat[-(57913:57915),-(4:17)]
# order sample info to match order of expression data (alphabetical)
sc_info <- sc_info[order(sc_info$sample),]
rownames(sc_info) <- 1:515

# the sample info only contains high quality samples where the expression matrix
# contains all samples so match the expression data samples with sample info
# by removing poor quality samples.
sc_match <- sc_dat[,colnames(sc_dat) %in% sc_info$sample]
sc_match <- cbind(sc_dat[,1:3], sc_match)
rownames(sc_match) <- sc_match[,1]
sc_match <- sc_match[,-(1:3)]

## Pre-process data using Seurat
library(Seurat)
# create seurat object
seurat_dat <- CreateSeuratObject(counts = sc_match)
# reduce dataset to only 10000 most differentially expressed genes
var_gene <- FindVariableFeatures(norm_sc, selection.method = 'vst', nfeatures = 10000)
# scale data
scale_sc_red <- ScaleData(object = var_gene)
# extract scaled and reduced expression matric from Seurat object
sc_pro_red <- GetAssayData(object = scale_sc_red, slot = "scale.data")
sc_pro_red <- as.data.frame(sc_pro_red)




