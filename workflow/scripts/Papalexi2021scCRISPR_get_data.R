
### libraries
library(Seurat)
library(SeuratData)
library(Matrix)

### config
matrix_path <- file.path(snakemake@output[["matrix"]])
barcodes_path <- file.path(snakemake@output[["barcodes"]])
features_path <- file.path(snakemake@output[["features"]])
metadata_path <- file.path(snakemake@output[["metadata"]])
metadata_selected_path <- file.path(snakemake@output[["metadata_selected"]])

# load data from SeuratData
SeuratData::InstallData('thp1.eccite')
data("thp1.eccite")

# Get counts for each modality
rna_counts <- GetAssayData(thp1.eccite, assay = "RNA", slot = "counts")
adt_counts <- GetAssayData(thp1.eccite, assay = "ADT", slot = "counts")
gdo_counts <- GetAssayData(thp1.eccite, assay = "GDO", slot = "counts")

# Combine count matrices
combined_counts <- rbind(rna_counts, adt_counts, gdo_counts)

# Create the features.tsv file with modality information
rna_features <- data.frame(
  id = rownames(rna_counts),
  name = rownames(rna_counts),
  type = "Gene Expression"
)
adt_features <- data.frame(
  id = rownames(adt_counts),
  name = rownames(adt_counts),
  type = "Antibody Capture"
)
gdo_features <- data.frame(
  id = rownames(gdo_counts),
  name = rownames(gdo_counts),
  type = "CRISPR Guide Capture"
)

combined_features <- rbind(rna_features, adt_features, gdo_features)

### Save the combined data in 10x format

# counts
writeMM(combined_counts, file = matrix_path)

# barcodes
write.table(
  colnames(combined_counts),
  file = barcodes_path,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# features
write.table(
  combined_features,
  file = features_path,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# save metadata
metadata <- thp1.eccite@meta.data
write.csv(metadata, file = metadata_path)


# save selected metadata
metadata <- metadata[,c("guide_ID","gene","replicate")]
write.csv(metadata, file = metadata_selected_path)


print("Data successfully saved in the MTX format.")