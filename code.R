# Load required libraries
library(Seurat)
library(Matrix)

# Set working directory
setwd("D:/Projects/sc_atlas_human_aging_tissues/GSE157827")

# Load Aged Sample (AD1)
aged_data <- Read10X(data.dir = ".", gene.column = 2,
                     files = c("GSM4775561_AD1_matrix.mtx.gz",
                               "GSM4775561_AD1_features.tsv.gz",
                               "GSM4775561_AD1_barcodes.tsv.gz"))

aged <- CreateSeuratObject(counts = aged_data, project = "Aged_AD1", min.cells = 3, min.features = 200)
aged$condition <- "aged"

# Load Young/Control Sample (NC3)
young_data <- Read10X(data.dir = ".", gene.column = 2,
                      files = c("GSM4775573_NC3_matrix.mtx.gz",
                                "GSM4775573_NC3_features.tsv.gz",
                                "GSM4775573_NC3_barcodes.tsv.gz"))

young <- CreateSeuratObject(counts = young_data, project = "Young_NC3", min.cells = 3, min.features = 200)
young$condition <- "young"

# Merge into one object for comparison
combined <- merge(young, y = aged, add.cell.ids = c("Young", "Aged"), project = "PBMC_Aging")

# Preview
combined
