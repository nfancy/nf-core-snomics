#!/usr/bin/env Rscript
# Check the samplesheet and file paths for cellranger
# 

##  ............................................................................
##  Load packages                                                           ####
suppressMessages(library(argparse))
suppressMessages(library(dplyr))

##  ............................................................................
##  Parse command-line arguments                                            ####

# create parser object
parser <- ArgumentParser()

required <- parser$add_argument_group("Required", "required arguments")

required$add_argument(
  "--input_file",
  help = "full path cellbender output dir",
  metavar = "dir",
  required = TRUE
)

required$add_argument(
  "--ensembl_mapping",
  help = "full path to ensembl_mapping file",
  metavar = "ensembl_mapping.tsv",
  required = TRUE
)

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args                                                        ####
args <- parser$parse_args()


#Read cellbender output_filtered.h5 file
  
input_file <-  args$input_file
ensemble_mapping <- read.delim(args$ensembl_mapping)


####function####
Read_CellBender_h5_Mat <- function(
    file_name,
    use.names = TRUE,
    unique.features = TRUE,
    h5_group_name = NULL,
    feature_slot_name = "features"
) {
  # Check hdf5r installed
  hdf5r_check <- SeuratObject::PackageCheck("hdf5r", error = FALSE)
  if (!hdf5r_check[1]) {
    cli_abort(message = c(
      "Please install the {.val hdf5r} package to use {.code Read_CellBender_h5_Mat} and read HDF5 files.",
      "i" = "This can be accomplished with the following commands: ",
      "----------------------------------------",
      "{.field `install.packages({symbol$dquote_left}hdf5r{symbol$dquote_right})`}",
      "----------------------------------------"
    ))
  }
  
  # Check file
  if (!file.exists(file_name)) {
    cli_abort(message = "File: {.val {file_name}} not found.")
  }
  
  # Check feature_slot_name is acceptable
  if (!feature_slot_name %in% c("features", "genes")) {
    cli_abort(message = c("{.code feature_slot_name} must be one of {.val features} or {.val genes}.",
                          "i" = "If unsure, check contents of H5 file {.code rhdf5::h5ls('{file_name}')}."))
  }
  
  # Read file
  infile <- hdf5r::H5File$new(filename = file_name, mode = "r")
  
  # Get list of H5 contents
  h5_dataset_list <- hdf5r::list.datasets(infile)
  
  # Check feature_slot_name is correct
  if (!length(x = grep(pattern = feature_slot_name, x = h5_dataset_list, value = TRUE)) > 0) {
    cli::cli_abort(message = c("{.code feature_slot_name}: {.val {feature_slot_name}} not found in H5 file.",
                               "i" = "Check contents of H5 file {.code rhdf5::h5ls('{file_name}')} to confirm correct {.code feature_slot_name}."))
  }
  
  # Assign feature slot name
  if (feature_slot_name == "features") {
    if (use.names) {
      feature_slot <- 'features/name'
    }
    else {
      feature_slot <- 'features/id'
    }
  }
  
  if (feature_slot_name == "genes") {
    if (use.names) {
      feature_slot <- 'gene_names'
    }
    else {
      feature_slot <- 'genes'
    }
  }
  
  # add name check
  group_names <- names(x = infile)
  
  if (!is.null(x = h5_group_name) && !h5_group_name %in% group_names) {
    cli::cli_abort(message = c("{.code h5_group_name} {.val {h5_group_name}} not found.",
                               "i" = "Check H5 file group names {.code rhdf5::h5ls('{file_name}')}."))
  }
  
  # Read in data
  if ("matrix" %in% group_names) {
    counts <- infile[["matrix/data"]]
    indices <- infile[["matrix/indices"]]
    indptr <- infile[["matrix/indptr"]]
    shp <- infile[["matrix/shape"]]
    features <- infile[[paste0("matrix/", feature_slot)]][]
    barcodes <- infile[["matrix/barcodes"]]
  } else {
    if (length(x = group_names) == 1) {
      counts <- infile[[paste0(group_names, '/data')]]
      indices <- infile[[paste0(group_names, '/indices')]]
      indptr <- infile[[paste0(group_names, '/indptr')]]
      shp <- infile[[paste0(group_names, '/shape')]]
      features <- infile[[paste0(group_names, '/', feature_slot)]][]
      barcodes <- infile[[paste0(group_names, '/barcodes')]]
    } else {
      # check subgroups
      if (is.null(x = h5_group_name)) {
        cli::cli_abort(message = c("H5 file contains multiple sub-groups.",
                                   "i" = "Please provide {.code h5_group_name} specifying which subgroup contains count data."))
      } else {
        counts <- infile[[paste0(h5_group_name, '/data')]]
        indices <- infile[[paste0(h5_group_name, '/indices')]]
        indptr <- infile[[paste0(h5_group_name, '/indptr')]]
        shp <- infile[[paste0(h5_group_name, '/shape')]]
        features <- infile[[paste0(h5_group_name, '/', feature_slot)]][]
        barcodes <- infile[[paste0(h5_group_name, '/barcodes')]]
      }
    }
  }
  
  # Create sparse matrix
  sparse.mat <- Matrix::sparseMatrix(
    i = indices[] + 1,
    p = indptr[],
    x = as.numeric(x = counts[]),
    dims = shp[],
    repr = "T"
  )
  
  if (unique.features) {
    features <- make.unique(names = features)
  }
  
  rownames(x = sparse.mat) <- features
  colnames(x = sparse.mat) <- barcodes[]
  sparse.mat <- Seurat::as.sparse(x = sparse.mat)
  
  infile$close_all()
  
  return(sparse.mat)
}

####function end####

cli::cli_text("Reading: {.path {input_file}}")
  
mat <- Read_CellBender_h5_Mat(file_name = input_file)
  
  
# Convert gene name to Ensembl gene id
ensembl_id <- ensemble_mapping[match(rownames(mat), 
                                     ensemble_mapping$external_gene_name), 
                                    "ensembl_gene_id"]
    idx <- which(is.na(ensembl_id))
    mat <- mat[-idx, ]
    ensembl_id <- ensembl_id[-idx]
    rownames(mat) <- ensembl_id

dest_path <- "cellbender_feature_bc_matrix"
dir.create(dest_path, recursive = TRUE)

cli::cli_text("Creating: {.path {dest_path}}")
  
DropletUtils::write10xCounts(path = dest_path,
                               x = mat,
                               barcodes = colnames(mat),
                               gene.id = rownames(mat),
                               gene.symbol = rownames(mat),
                               gene.type = "Gene Expression",
                               overwrite = TRUE,
                               type = "auto",
                               genome = "unknown",
                               version = "3")

cli::cli_text("Sparse matrix generation is successful")


