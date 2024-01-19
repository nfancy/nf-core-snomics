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
  "--metrics_summary",
  help = "full path to the metrics_summary.csv file",
  metavar = "metrics_summary.csv",
  required = TRUE
)

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args                                                        ####
args <- parser$parse_args()

metrics_summary <- read.csv(args$metrics_summary)
colnames(metrics_summary) <- tolower(colnames(metrics_summary))

n_cells <- gsub(",", "", metrics_summary[["estimated.number.of.cells"]]) %>%
           as.numeric()

cat(n_cells)