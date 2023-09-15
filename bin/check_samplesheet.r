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
  "--input",
  help = "full path to the sample sheet csv file",
  metavar = "samplesheet.csv",
  required = TRUE
)

required$add_argument(
  "--aligner",
  help = "name of the cellranger software to be used as aligner",
  metavar = "cellranger",
  required = TRUE
)

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args                                                        ####
args <- parser$parse_args()

if (!file.exists(args$input)) {
  stop("The input samplesheet was not found.")
}

input <- read.csv(args$input)

if( args$aligner == "cellranger" ){
  expected_cols <- c("sample", "fastq_1", "fastq_2")
} else if (args$aligner == "cellranger_arc") {
  expected_cols <- c("fastqs", "sample", "library_type")
}

assertthat::assert_that(all(expected_cols %in% colnames(input)),
                        msg = sprintf("%s columns are expected", 
                                      paste0(expected_cols, collapse = ","))
                        )

# check sample paths exist

check_exists <- function(filepath) {
  RCurl::url.exists(filepath) |
    file.exists(filepath) |
    any(startsWith(filepath, c("gs://", "s3://")))
}

file_list <- c(input$fastq_1, input$fastq_2)

dir_exists <- purrr::map_lgl(file_list, ~ check_exists(.))

if (!all(dir_exists)) {
  cat("The following paths were not found: -\n")
  print(file_list[!dir_exists])
  stop("Fastq paths specified in the input samplesheet were not found.")
} else {
  cat("✓ All paths specified in the input samplesheet were found.\n")
}

cat("Checks passed!\n")

# # write the same manifest back out
write.csv(x = input, 
          file = "checked_samplesheet.csv", 
          quote = FALSE,
          row.names = FALSE)

