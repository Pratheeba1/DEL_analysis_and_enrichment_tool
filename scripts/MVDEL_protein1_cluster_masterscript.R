rm(list = ls())

# Prefer SLURM array ID, fallback to command line argument
task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")

if (nzchar(task_id)) {
  which_protein <- as.integer(task_id)
} else {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) stop("Usage: Rscript MVDEL_protein1_cluster.R <sample_index>")
  which_protein <- as.integer(args[1])
}

if (is.na(which_protein)) stop("sample_index must be an integer")

library(readxl)
library(ShortRead)
library(Biostrings)
library(seqinr)
library(magrittr)
library(dplyr)
library(gtools)
library(stringr)

# helper functions
analyze_dataset <- source("ClusterCode/analyze_datasets.R")$value
get_sequence_start <- source("ClusterCode/get_sequence_start.R")$value
get_indices_dataset <- source("ClusterCode/get_indices_dataset.R")$value
check_single_code <- source("ClusterCode/check_single_code.R")$value
check_single_code_for_insertion <- source("ClusterCode/check_single_code_for_insertion.R")$value
check_single_code_for_miss_element <- source("ClusterCode/check_single_code_for_miss_element.R")$value
check_for_combined_errors <- source("ClusterCode/check_for_combined_errors.R")$value
analyze_sequence_onesheet <- source("ClusterCode/analyze_sequence_onesheet.R")$value
sheet_info <- source("ClusterCode/sheet_info.R")$value
comparing_all_patterns_onesheet <- source("ClusterCode/comparing_all_patterns_onesheet.R")$value

# structure sheet
structure_sheet <- read_excel("../data/structure_sheet.xlsx")

# trim text columns to avoid hidden space problems
structure_sheet <- structure_sheet %>%
  mutate(across(where(is.character), trimws))

# sample folders
sample_dirs <- list.dirs("../data", full.names = TRUE, recursive = FALSE)
sample_dirs <- sample_dirs[dir.exists(sample_dirs)]
sample_dirs <- sort(sample_dirs)

if (which_protein < 1 || which_protein > length(sample_dirs)) {
  stop("sample_index out of range")
}

sample_dir <- sample_dirs[which_protein]
sample_name <- basename(sample_dir)
out_file <- paste0("../results/result_dataset_", sample_name, ".RData")

# skip if already done
if (file.exists(out_file)) {
  cat("Skipping", sample_name, "- already finished\n")
  quit(save = "no", status = 0)
}

# match structure sheet row by sample/protein name, NOT by order
protein_row <- which(structure_sheet$ProteinCode1_name == sample_name)

if (length(protein_row) != 1) {
  stop(paste(
    "Could not uniquely match sample",
    sample_name,
    "in structure_sheet. Matches found:",
    length(protein_row)
  ))
}

cat("=================================\n")
cat("Processing sample:", sample_name, "\n")
cat("Sample folder index:", which_protein, "\n")
cat("Structure sheet row:", protein_row, "\n")
cat("=================================\n")

# read all possible FASTQ files: R1, R2, or single-end
fastq_files <- list.files(sample_dir, pattern = "\\.fq\\.gz$", full.names = TRUE)

if (length(fastq_files) == 0) {
  stop(paste("No .fq.gz file found in", sample_dir))
}

fastq_files <- sort(fastq_files)

sheet_inf <- sheet_info("MVDEL_1", structure_sheet, protein_row, prime_length = 8)

best_result <- NULL
best_file <- NULL
best_score <- -1

for (fastq_file in fastq_files) {
  
  cat("=================================\n")
  cat("Testing FASTQ:", fastq_file, "\n")
  cat("=================================\n")
  
  data_fastq <- readFastq(fastq_file)
  matching_indices <- get_indices_dataset(length(data_fastq), step = 100000)
  
  cat("Number of reads:", length(data_fastq), "\n")
  cat("Number of chunks:", nrow(matching_indices), "\n")
  
  all_results <- vector("list", nrow(matching_indices))
  
  for (p in seq_len(nrow(matching_indices))) {
    cat("Processing", basename(fastq_file), "chunk", p, "of", nrow(matching_indices), "\n")
    
    sequences_dataset <- data_fastq@sread[
      matching_indices[p, 1]:matching_indices[p, 2]
    ]
    
    chunk_result <- analyze_dataset(
      sequences_dataset = sequences_dataset,
      sheet_inf = sheet_inf,
      thresholds_two_errors = 20,
      error_threshold = 2
    )
    
    all_results[[p]] <- chunk_result
    
    rm(sequences_dataset, chunk_result)
    gc()
  }
  
  result_tmp <- bind_rows(all_results)
  
  score <- sum(result_tmp$Region1Code_smile != "No Smile", na.rm = TRUE) +
    sum(result_tmp$Region2Code_smile != "No Smile", na.rm = TRUE) +
    sum(result_tmp$Region3Code_smile != "No Smile", na.rm = TRUE)
  
  cat("FASTQ:", fastq_file, "\n")
  cat("Good decoding score:", score, "\n")
  
  if (score > best_score) {
    best_score <- score
    best_result <- result_tmp
    best_file <- fastq_file
  }
  
  rm(data_fastq, result_tmp, all_results)
  gc()
}

result_dataset <- best_result

attr(result_dataset, "sample_name") <- sample_name
attr(result_dataset, "structure_sheet_row") <- protein_row
attr(result_dataset, "selected_fastq") <- best_file
attr(result_dataset, "decoding_score") <- best_score

save(result_dataset, file = out_file)

cat("=================================\n")
cat("Saved:", out_file, "\n")
cat("Sample:", sample_name, "\n")
cat("Structure sheet row:", protein_row, "\n")
cat("Selected FASTQ:", best_file, "\n")
cat("Best decoding score:", best_score, "\n")
cat("Rows:", nrow(result_dataset), " Cols:", ncol(result_dataset), "\n")
cat("=================================\n")