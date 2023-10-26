################################################################################
##########                            Init                            ##########
################################################################################
# Packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(Biobase)
})

# Shared functions
source('code/shared.R')
source('code/functions.R')

################################################################################
##########               Iterate over FASTQ folders                   ##########
################################################################################

# It is possible to construct a config file with BCL folders without included 
# reactions. Since this function is fast, and we may be interested in running 
# the pipeline on BCL folders without considering the reactions, I am iterating
# over all BLC folders in the config file.
sequencing_sheet <- dplyr::as_tibble(snakemake@config[["sequencing_sheet"]])

# Do we want all folders to have the same timepoint?
timestamp <- lubridate::now("UTC")

print_header <- function(this_run, timestamp)
{
  bcl <- this_run[["bcl_folder"]]
  sequencing_id <- this_run[["sequencing_id"]]
  override_cycles <- this_run[["override_cycles"]]
  
  cat(rep('#', 80), '\n',
      '#####                                                                      #####\n',
      '#####                        Collecting FASTQ stats                        #####\n',
      '#####                                                                      #####\n',
      '#####                             3/6                                      #####\n',
      '#####           ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▒▒▒▒▒▒▒▒░░░░░░░░░░░░░░░░░░░░░░░░            #####\n',
      '#####                                                                      #####\n',
      rep('#',80), '\n',
      '#####\n',
      '##\n',
      '##\tSCOP ID:      \t\t', snakemake@config[['scop_id']], '\n',
      '##\tCOMUNEQAID ID:\t\t', snakemake@config[['com_id']], '\n',
      '##\n',
      '##\tsequencing_id:\t\t- ', sequencing_id, '\n',
      '##\tbcl_folder:   \t\t- ', bcl,'\n',
      '##\toverride_cycles:\t- ', override_cycles, '\n',
      '##\n',
      '##\n',
      '##\t', format_time(timestamp, "pretty"), '\n',
      '##\n',
      '##\n',
      '###\n',
      rep('#', 80), '\n',
      '#\n',
      sep = '')
}

fastq_summaries <- function(bcl, config = snakemake@config) {
  
}

for (bcl in sequencing_sheet[["bcl_folder"]]) {
  this_run <- filter(sequencing_sheet, bcl_folder == bcl)
  
  out_folder <- file.path(fastq_path(), bcl, "metadata")
  
  
  read.tib <-
    tibble(
      index = character(),
      reads = numeric(),
      bcl_folder = character(),
      class = character()
    )
  
  unknown.barcodes.tib <-
    tibble(
      dual_index = character(),
      reads = numeric(),
      bcl_folder = character()
    )
  
  fastq_stats_path <- file.path(
    project_path(),
    snakemake@config[['fastq_path']],
    bcl,
    'metadata')
  
  log_file <- file.path(
    fastq_stats_path,
    paste0(date_and_time(),'.log'))
  
  dir.create(fastq.stats.path, recursive = T, showWarnings = F)
  
  # Sink
  sink.file <- file(log.file, 'w')
  sink(file = sink.file, append = T,
       type = 'output', split = T)
  
  # Startup message
  print_header()
  
  cat(paste0('#\t\tFASTQ set:\t',bcl,'\n'),
      sep = '')
  
  cat(paste0('#\t..\timporting demultiplexing stats..\n'),
      sep = '')
  
  reads <- readr::read_csv(demult.stats.path, 
                  col_select = c("index" = "SampleID", "reads" = "# Reads"), 
                  col_types = c("index" = readr::col_character(), 
                                "reads" = readr::col_integer())) |>
    dplyr::group_by(index) |>
    dplyr::summarise(reads = sum(reads))
  
  types <- merge_sheets("reaction2library", "library_sheet", config = config) |>
    select(-c("library_id", "reaction_id"))
  types <- rbind(types, c("Undetermined", "Undetermined"))
  merge.data.frame(types, reads)
  
  
  stats_folder <- file.path(project_path(),'scRNAseq','dry-lab','FASTQ',bcl,bcl.convert.version,'Reports')
  
  demult.stats.path <- file.path(project_path(),'scRNAseq','dry-lab','FASTQ',bcl,bcl.convert.version,'Reports','Demultiplex_Stats.csv')
  read.tib <- read_csv(demult.stats.path,
                           col_select = c(index = 'SampleID', reads = '# Reads'),
                           col_types = c('c','d')) %>%
    group_by(index) %>%
    summarise(reads = sum(reads)) %>% 
    mutate(bcl_folder = bcl,
           class = case_when(str_detect(index, '^SI-TT-') ~ '10x',
                             str_detect(index, '^D7') ~ 'HTO',
                             TRUE ~ 'Undetermined'))
  
  write.table(x = read.tib,
              file = file.path(fastq.stats.path,'read-demultiplexing.csv'),
              sep = ',',
              row.names = F)
  
  cat(paste0('#\t..\timporting unknown barcode stats..\n',
             '#\t..\n'),
      sep = '')
  unknown.barcodes.path <- file.path(project_path(),'scRNAseq','dry-lab','FASTQ',bcl,bcl.convert.version,'Reports','Top_Unknown_Barcodes.csv')
  unknown.barcodes.tib <- read_csv(unknown.barcodes.path,
                                   col_select = c(lane = 'Lane', index = index, index2 = index2, reads = '# Reads'),
                                   col_types = c('c','c','c','d')) %>%
    mutate(dual_index = paste(index, index2, sep = '_')) %>%
    group_by(dual_index) %>% 
    summarise(reads = sum(reads)) %>% 
    arrange(desc(reads)) %>% 
    mutate(bcl_folder = bcl)
  
  write.table(x = unknown.barcodes.tib,
              file = file.path(fastq.stats.path,'unknown-barcodes.csv'),
              sep = ',',
              row.names = F)
  cat('#\t..done\n')
  
  cat('#\n',
      '#\n',
      '###\n',
      rep('#',80),'\n',
      sep = '')
  
  sink()
}


file.create(snakemake@output[[1]])