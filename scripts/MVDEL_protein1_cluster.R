### Masterscript zum Durchlauf eines Datensatzes auf dem Computecluster
rm(list = ls())
p <- as.integer(Sys.getenv("PBS_ARRAYID"))

### Laden der Pakete
library(readxl)
library(ShortRead)
library(Biostrings)
library(seqinr)
library(magrittr)
library(dplyr)
library(gtools)
library(stringr)

which_protein <- 1

data_fastq <- readFastq(paste0("/SHK/DEL Daten/MVDEL/01_fastq",which_protein,".1.fastq.gz"))

analyze_datasets <- source("/SHK/DEL Daten/DELCode/ClusterCode/analyze_datasets.R")$value
get_sequence_start <-source("/SHK/DEL Daten/DELCode/ClusterCode/get_sequence_start.R")$value
get_indices_dataset <- source("/SHK/DEL Daten/DELCode/ClusterCode/get_indices_dataset.R")$value
check_single_code <- source("/SHK/DEL Daten/DELCode/ClusterCode/check_single_code.R")$value
check_single_code_for_insertion <- source("/SHK/DEL Daten/DELCode/ClusterCode/check_single_code_for_insertion.R")$value
check_single_code_for_miss_element <- source("/SHK/DEL Daten/DELCode/ClusterCode/check_single_code_for_miss_element.R")$value
check_for_combined_errors <- source("/SHK/DEL Daten/DELCode/ClusterCode/check_for_combined_errors.R")$value
analyze_sequence_onesheet <- source("/SHK/DEL Daten/DELCode/ClusterCode/analyze_sequence_onesheet.R")$value
sheet_info <- source("/SHK/DEL Daten/DELCode/ClusterCode/sheet_info.R")$value
comparing_all_patterns_onesheet <- source("/SHK/DEL Daten/DELCode/ClusterCode/comparing_all_patterns_onesheet.R")$value

structure_sheet <- read_excel("/SHK/DEL Daten/DEL Daten/MVDEL/CEGAT_MVDEL_Pool3_StructureSheet.xlsx")

sheet_inf <- sheet_info("MVDEL_1",structure_sheet,which_protein,prime_length = 8)
matching_indices <- get_indices_dataset(length(data_fastq),step = 100000)
seq_data <- data_fastq

sequences_dataset <- data_fastq@sread[matching_indices[p, 1]:matching_indices[p, 2]]
result_dataset<- analyze_dataset(sequences_dataset = sequences_dataset,
                                 sheet_inf = sheet_inf,
                                 thresholds_two_errors = 20,
                                 error_threshold = 2)

save(result_dataset,
     file = paste("result_dataset_protein",which_protein,"_", p,".RData", sep = ""))

rm(list = ls())
gc()

