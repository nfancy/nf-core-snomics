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
  expected_cols <- c("sample", "fastqs")
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
    dir.exists(filepath) |
    any(startsWith(filepath, c("gs://", "s3://")))
}

dir_exists <- purrr::map_lgl(input[["fastqs"]], ~ check_exists(.))

if (!all(dir_exists)) {
  cat("The following paths were not found: -\n")
  print(input[!dir_exists, ])
  stop("Directory paths specified in the input samplesheet were not found.")
} else {
  cat("✓ All paths specified in the input samplesheet were found.\n")
}


# check both R1 and R2 fastq files are present in the dir
check_fastq_exists <- function(dirpath) {
  
  files <- list.files(path = dirpath, pattern = "R1|R2", full.names = TRUE)
  
  if(length(files) == 0){
  cat("The following direcotry is empty: -\n")
  print(dirpath)
  return(FALSE)
  } else {
  
  unique_files <- unique(gsub("_R[1-2]_.*","",files))
  file_list <- vector('list',length(unique_files))
  for(f in seq_along(unique_files)){
    files_sub <- stringr::str_detect(files, unique_files[f])
    file_pairs <- files[files_sub]
    l <- length(file_pairs)
    if (l == 2) { file_list[[f]] <- TRUE} 
    else if (l == 1) { 
      
      cat("The following file pairs were not found: -\n")
      print(file_pairs)
      file_list[[f]] <- FALSE 
    } 
  }
  file_list <- unlist(file_list)
  if (!all(file_list)){ return(FALSE) } 
  else { return(TRUE) }
  }
}
  

fastq_pair_exists <- purrr::map_lgl(input[["fastqs"]],
                                    ~ check_fastq_exists(dirpath = .))

if (!all(fastq_pair_exists)) {
  cat("The following file pairs were not found: -\n")
  print(input[!fastq_pair_exists, ])
  stop("Fastq files were not found.")
} else {
  cat("✓ All directories specified in the input samplesheet contain fastq files.\n")
}

cat("Checks passed!\n")


file_l <- lapply(input[["fastqs"]], function(x){
  list.files(path = x, pattern = "R1|R2", full.names = TRUE)
})
names(file_l) <- basename(input[["fastqs"]])
dt <- stack(file_l)
dt <- dt %>%
  mutate(fastq = ifelse(grepl("R1", values), "fastq_1", "fastq_2"))

fastq_input <- dt %>% 
  tidyr::pivot_wider(names_from = "fastq",
              values_from = "values") %>%
  tidyr::unnest(cols = c(fastq_1, fastq_2)) %>%
  dplyr::rename(sample = ind)

# write input for fastqc

write.csv(x = fastq_input, 
          file = "checked_samplesheet.csv", 
          quote = FALSE,
          row.names = FALSE)



# # write the same manifest back out
# write.csv(input,
#           "checked_samplesheet.csv",
#           quote = FALSE,
#           row.names = FALSE)
