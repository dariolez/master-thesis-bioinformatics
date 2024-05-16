# AUTHOR: Darío González
# DESCRIPTION: Merging all samples into one Seurat object

# IMPORTS
library(BUSpaRse) # github (devtools)
library(Seurat)
library(SeuratWrappers) # github (devtools)
library(zeallot) # For %<-% that unpacks lists in the Python manner
library(DropletUtils)
library(plyr)
library(tidyverse)
library(Matrix)

# MAIN

# Get directories
directories <- dir(path = "~/TFM/2019.12.31_Meyer/data",
                   pattern = "raw",
                   full.names = TRUE,
                   recursive = TRUE,
                   include.dirs = TRUE)

# Load run uuid to file name
bundle2filename <- read.table("~/TFM/2019.12.31_Meyer/bundle_filename_10x.tsv", header = TRUE, sep = "\t", quote = "")

# Read in the matrices
read_velocity <- function(dir=getwd(),
                          spliced_file="spliced.mtx",
                          unspliced_file="unspliced.mtx",
                          barcodes="barcodes.txt",
                          genes="genes.txt",
                          ...) {
    spliced <- ReadMtx(mtx = paste0(dir, "/", spliced_file),
                       cells = paste0(dir, "/", barcodes),
                       features = paste0(dir, "/", genes),
                       ...)
    unspliced <- ReadMtx(mtx = paste0(dir, "/", unspliced_file),
                         cells = paste0(dir, "/", barcodes),
                         features = paste0(dir, "/", genes),
                         ...)
    return(list(spliced, unspliced))
}

# Make Seurat object from all individual spliced/unspliced files:
counter <- 0

for (run_dir in directories) {
    # Update counter
    counter <- counter + 1
    
    # State the file being processed
    cat("\n", counter, "Processing", run_dir, "\n")
    
    # Load the count matrices for spliced and unspliced reads
    c(spliced, unspliced) %<-%  read_velocity(dir = run_dir,
                                              spliced_file = "spliced.mtx",
                                              unspliced_file = "unspliced.mtx",
                                              barcodes = "barcodes.tsv",
                                              genes = "features.tsv",
                                              cell.column = 1,
                                              feature.column = 2,
                                              mtx.transpose = FALSE)

    # Filter matrices so that barcodes that appear in both remain
    # sf = spliced filtered, uf = unspliced filtered
    bcs_use <- intersect(colnames(spliced), colnames(unspliced))
    sf <- spliced[, bcs_use]
    uf <- unspliced[, bcs_use]
    
    # Create Seurat Object
    seu <- CreateSeuratObject(sf, assay = "spliced")
    seu[["unspliced"]] <- CreateAssayObject(uf)
    
    # Rename cells to avoid barcode collision when merging
    bundle_row <- grep(str_extract(run_dir, "([0-9a-z]+-){4}[0-9a-z]+"), bundle2filename$bundle_uuid)
    newbarcodes <- paste0(colnames(seu), "-1-", rep(bundle2filename[bundle_row, "file_name"], length(colnames(seu))))
    seu <- RenameCells(seu, new.names = newbarcodes)
    
    # Merge Seurat objects
    if (counter == 1) {
        myesoph <- seu
    } else {
        myesoph <- merge(myesoph, seu)
    }
    
    rm(seu)
}

myesoph

# Add feature metadata
# Load an Ensembl gene ID to HGNC symbol `.genes.map.txt` file. 
# All files for spliced and unspliced samples are identical so it doesn't matter which one we load
ensembl2symbol <- read.table(paste0(directories[1], "/", "features.tsv"), header = FALSE, sep = "\t")

# Rows are left ordered after merges (I checked), so we can add Ensembl IDs directly without more checks
#myesoph@assays$spliced <- AddMetaData(myesoph@assays$spliced, metadata = ensembl2symbol[[1]], col.name = "ensembl_id")
#myesoph@assays$unspliced <- AddMetaData(myesoph@assays$unspliced, metadata = ensembl2symbol[[1]], col.name = "ensembl_id")

# If we want to check before assigning Ensembl IDs to symbols
mysymbols_match <- match(rownames(myesoph), make.unique(ensembl2symbol[[2]]))
myesoph@assays$spliced <- AddMetaData(myesoph@assays$spliced, metadata = ensembl2symbol[[1]][mysymbols_match], col.name = "ensembl_id")
myesoph@assays$unspliced <- AddMetaData(myesoph@assays$unspliced, metadata = ensembl2symbol[[1]][mysymbols_match], col.name = "ensembl_id")

# Save merged Seurat Object
saveRDS(myesoph, file = "~/TFM/2019.12.31_Meyer/results/myesoph_star.rds")
