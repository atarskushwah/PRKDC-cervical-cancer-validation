# ==============================================================================
# scripts/06a_gse224327_download.R
# ==============================================================================
library(httr)
library(utils)

cat("=============================================================\n")
cat("--- 1. Initializing Chunked Stream Download ---\n")
cat("=============================================================\n")

# Define the target directory and file paths
sc_dir <- "data/raw/GSE224327_scRNA"
if (!dir.exists(sc_dir)) dir.create(sc_dir, recursive = TRUE)

tar_file <- file.path(sc_dir, "GSE224327_RAW.tar")
geo_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE224327&format=file"

# We use httr::GET with write_disk() to stream the data in small packets directly 
# to the SSD, completely bypassing RAM to prevent memory overflow crashes.
cat(sprintf("Streaming single-cell data from NCBI GEO to %s...\n", tar_file))
cat("This may take several minutes depending on your connection. Please wait.\n")

res <- GET(
  url = geo_url,
  write_disk(tar_file, overwrite = TRUE),
  progress() # This will show you a live progress bar in the console!
)

if (status_code(res) == 200) {
  cat("\nSUCCESS: Data stream completed flawlessly.\n")
  
  cat("\n--- 2. Unpacking Single-Cell Archive ---\n")
  # Extract the downloaded .tar file into the single-cell directory
  untar(tar_file, exdir = sc_dir)
  cat(sprintf("Extracted raw matrix files to: '%s'\n", sc_dir))
  
  # Optional: Clean up the massive .tar file to save hard drive space
  file.remove(tar_file)
  cat("Cleaned up temporary .tar archive.\n")
  
} else {
  cat(sprintf("\nERROR: Download failed with HTTP status %s\n", status_code(res)))
}
cat("=============================================================\n")