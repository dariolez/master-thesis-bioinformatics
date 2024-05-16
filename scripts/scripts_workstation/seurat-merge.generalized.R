# AUTHOR: Darío González
# DESCRIPTION: Merging all samples into one Seurat object

############################################
#
# This script isn't completely functional. Use the specific version for each tool results
#
############################################

# PARSE COMMAND LINE OPTIONS
library(optparse)

option_list <- list(
    make_option(c("-d", "--dir"), action="store", default=getwd(), help="Base directory with sample folders. [Default current wd]"),
    make_option(c("-p", "--pattern"), action="store", help="Folder where the count matrices reside."),
    make_option(c("-s", "--spliced"), action="store", help=".mtx file with spliced counts"),
    make_option(c("-u", "--unspliced"), action="store", help=".mtx file with unspliced counts"),
    make_option(c("-b", "--barcodes"), action="store", help="File with barcodes"),
    make_option(c("-g", "--genes"), action="store", help="File with gene names"), 
    make_option(c("-t", "--tool"), action="store", help="Tool that generated the count matrices. One of ['bustools', 'kb-count', 'star']"),
    make_option(c("-y", "--symbols"), action="store_true", help="Use gene symbols for row names"),
    make_option(c("-m", "--map-id-symbol", dest="id2gene", action="store", help="File with Ensembl ID to gene symbol mapping"),
    make_option(c("-w", "--write"), action="store", help="Name of the .rds file to write")
)

args <- parse_args(OptionParser(option_list=option_list))

# Check arguments
if (! args$tool %in% c('bustools', 'kb-count', 'star')) {
    stop("Invalid tool. You must use one of ['bustools', 'kb-count', 'star']")
}

if (is.null(args$tool) || is.null(args$write) || is.null(args$dir) || is.null(args$pattern) ||
    is.null(args$spliced) || is.null(args$unspliced)) {
    stop("Missing arguments.")
}

if (args$tool %in% c('kb-count', 'star') && (is.null(args$barcodes) || is.null(args$genes)) ) {
    stop("For tools 'kb-count' and 'star', arguments 'barcodes' and 'genes' are required.")    
}

if (args$symbol && is.null(args$id2gene)) {
    stop("A file mapping gene symbols to Ensemble IDs has to be supplied through '--map-id-symbol' when using '--symbol'.")
}

if (args$tool == "bustools") { transpose_mtx <- TRUE }
if (args$tool == "kb-count") { transpose_mtx <- TRUE }
if (args$tool == "star") { 
    transpose_mtx <- FALSE 
    if (args$symbols) { 
        feature.column <- 2 
    } else {
        feature.column <- 1
    }
}

# IMPORTS

library(BUSpaRse) # github (devtools)
library(Seurat)
library(SeuratWrappers) # github (devtools)
library(zeallot) # For %<-% that unpacks lists in the Python manner
library(DropletUtils)
library(plyr)
library(tidyverse)
library(Matrix)
library(GGally)
library(scales)
library(plotly)

# FUNCTIONS

# Finding duplicate gene names
find_dup_genes <- function(x, mode = "count") {
    
    dup_genes_rows <- which(duplicated(rownames(x)) | duplicated(rownames(x), fromLast = TRUE))
    dup_genes <- rownames(x)[dup_genes_rows]
    
    if (mode == "count") {
        dup_genes <- as.data.frame(table(dup_genes))
    } else if (mode == "unique") {
        dup_genes <- unique(dup_genes)
    } else if (mode == "all") {
        return(dup_genes)
    } else {
        stop("Unrecognised mode. Only 'count', 'unique' and 'all' are valid.\n")
    }
    
    return(dup_genes)
}

# MAIN

# Get directories
directories <- dir(path = args$dir,
                   pattern = args$pattern,
                   full.names = TRUE,
                   recursive = TRUE,
                   include.dirs = TRUE)

# Directories for kb-python 0.27.3
#directories <- dir(path = "~/TFM/2019.12.31_Meyer/data", pattern = "bus_kb-devel_k0.50.0", full.names = TRUE, recursive = TRUE, include.dirs = TRUE)

# Directories For kb-python 0.28.0
#directories <- dir(path = "~/TFM/2019.12.31_Meyer/data", pattern = "counts_unfiltered", full.names = TRUE, recursive = TRUE, include.dirs = TRUE)

# Directories for STAR
#directories <- dir(path = "~/TFM/2019.12.31_Meyer/data", pattern = "raw", full.names = TRUE, recursive = TRUE, include.dirs = TRUE)

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

# Test merges



# Make Seurat object from all individual spliced/unspliced files:
counter <- 0

for (run_dir in directories) {
    # Update counter
    counter <- counter + 1
    
    # State the file being processed
    cat("\n", counter, "Processing", run_dir, "\n")
    
    # Load the count matrices for spliced and unspliced reads
    
    if (args$tool == "bustools") {
        c(spliced, unspliced) %<-% read_velocity_output(spliced_dir = run_dir,
                                                        spliced_name = "spliced",
                                                        unspliced_dir = run_dir,
                                                        unspliced_name = "unspliced")
    } else {
        c(spliced, unspliced) %<-%  read_velocity(dir = run_dir,
                                                  spliced_file = args$spliced,
                                                  unspliced_file = args$unspliced,
                                                  barcodes = args$barcodes,
                                                  genes = args$genes,
                                                  cell.column = 1,
                                                  feature.column = feature.column,
                                                  mtx.transpose = transpose_mtx)
    }

    # Filter matrices so that barcodes that appear in both remain
    # sf = spliced filtered, uf = unspliced filtered
    bcs_use <- intersect(colnames(spliced), colnames(unspliced))
    sf <- spliced[, bcs_use]
    uf <- unspliced[, bcs_use]
    
    # Create Seurat Object
    seu <- CreateSeuratObject(sf, assay = "spliced")
    seu[["unspliced"]] <- CreateAssayObject(uf)
    
    # Rename cells to avoid barcode collision when merging
    bundle_row <- grep(basename(dirname(run_dir)), bundle2filename$bundle_uuid)
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

if (args$symbol) {
    # Add feature metadata
    # Load an Ensembl gene ID to HGNC symbol `.genes.map.txt` file. 
    # All files for spliced and unspliced samples are identical so it doesn't matter which one we load
    ensembl2symbol <- read.table(paste0(directories[1], args$id2gene), header = FALSE, sep = "\t")

    # Rows are left ordered after merges (I checked with `is.unsorted()`), so we can add Ensembl IDs directly without more checks
    myesoph@assays$spliced <- AddMetaData(myesoph@assays$spliced, metadata = ensembl2symbol[[1]], col.name = "ensembl_id")
    myesoph@assays$unspliced <- AddMetaData(myesoph@assays$unspliced, metadata = ensembl2symbo[[1]], col.name = "ensembl_id")

    # If we want to check before assigning Ensembl IDs to symbols
    #mysymbols_match <- match(rownames(myesoph), make.unique(ensembl2symbol[[2]]))
    #myesoph@assays$spliced <- AddMetaData(myesoph@assays$spliced, metadata = ensembl2symbol[[1]][mysymbols_match], col.name = "ensembl_id")
    #myesoph@assays$unspliced <- AddMetaData(myesoph@assays$unspliced, metadata = ensembl2symbol[[1]][mysymbols_match], col.name = "ensembl_id")
}

# Save merged Seurat Object
saveRDS(myesoph, file = args$args)
