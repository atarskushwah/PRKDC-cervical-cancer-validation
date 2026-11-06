library(Seurat)
library(GEOquery)

extracted_path <- "data/processed/GSE224327_raw"
if (!dir.exists(extracted_path)) stop("Missing extracted 10X directory: ", extracted_path)

cat("Loading 10X multi-sample single-cell data...\n")
# Read10X automatically scans the directory for barcodes, features, and matrix files
raw_data <- Read10X(data.dir = extracted_path)

# If it's a list (multiple samples), create a combined Seurat object or take the main one
if(is.list(raw_data)) {
  sc_obj <- CreateSeuratObject(counts = raw_data[[1]], project = "GSE224327")
} else {
  sc_obj <- CreateSeuratObject(counts = raw_data, project = "GSE224327")
}

saveRDS(sc_obj, "data/processed/gse224327_seurat_object.rds")
cat("SUCCESS: Seurat object saved successfully!\n")
