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

# Paths
project.path <-
  file.path(snakemake@config[['project_path']],
            snakemake@config[['scop_id']])

# Metadata
suppressMessages({
  rnx2lib.libsheet.lib2seq.seqsheet <-
    full_join(as_tibble(snakemake@config[['reaction2library']]),
              as_tibble(snakemake@config[['library_sheet']])) %>% 
    full_join(as_tibble(snakemake@config[['library2sequencing']])) %>% 
    full_join(as_tibble(snakemake@config[['sequencing_sheet']]))
})


################################################################################
##########               Iterate over FASTQ folders                   ##########
################################################################################

for (bcl in unique(rnx2lib.libsheet.lib2seq.seqsheet[['bcl_folder']])) {
  
  sequencing_id <- as_tibble(snakemake@config[['sequencing_sheet']]) %>% 
    filter(bcl_folder == bcl) %>% 
    select(sequencing_id) %>% 
    unlist()
  override_cycles <- as_tibble(snakemake@config[['sequencing_sheet']]) %>% 
    filter(bcl_folder == bcl) %>% 
    select(override_cycles) %>% 
    unlist()
  
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
  
  fastq.stats.path <- file.path(
    snakemake@config[['project_path']],
    snakemake@config[['scop_id']],
    snakemake@config[['fastq_path']],
    bcl,
    'metadata')
  
  log.file <- file.path(
    fastq.stats.path,
    paste0(date.and.time,'.log'))
  
  dir.create(fastq.stats.path, recursive = T, showWarnings = F)
  
  # Sink
  sink.file <- file(log.file, 'w')
  sink(file = sink.file, append = T,
       type = 'output', split = T)
  
  # Startup message
  cat(rep('#',80),'\n',
      '#####                                                                      #####\n',
      '#####                        Collecting FASTQ stats                        #####\n',
      '#####                                                                      #####\n',
      '#####                             3/6                                      #####\n',
      '#####           ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▒▒▒▒▒▒▒▒░░░░░░░░░░░░░░░░░░░░░░░░            #####\n',
      '#####                                                                      #####\n',
      rep('#',80),'\n',
      '#####\n',
      '##\n',
      '##\tSCOP ID:      \t\t',snakemake@config[['scop_id']],'\n',
      '##\tCOMUNEQAID ID:\t\t',snakemake@config[['com_id']],'\n',
      '##\n',
      '##\tsequencing_id:\t\t- ',sequencing_id,'\n',
      '##\tbcl_folder:   \t\t- ',bcl,'\n',
      '##\toverride_cycles:\t- ',override_cycles,'\n',
      '##\n',
      '##\n',
      '##\t',date.and.time.pretty,'\n',
      '##\n',
      '##\n',
      '###\n',
      rep('#',80),'\n',
      '#\n',
      sep = '')
  
  cat(paste0('#\t\tFASTQ set:\t',bcl,'\n'),
      sep = '')
  
  cat(paste0('#\t..\timporting demultiplexing stats..\n'),
      sep = '')
  
  demult.stats.path <- file.path(project.path,'scRNAseq','dry-lab','FASTQ',bcl,bcl.convert.version,'Reports','Demultiplex_Stats.csv')
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
  unknown.barcodes.path <- file.path(project.path,'scRNAseq','dry-lab','FASTQ',bcl,bcl.convert.version,'Reports','Top_Unknown_Barcodes.csv')
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