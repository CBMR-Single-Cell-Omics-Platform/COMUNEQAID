################################################################################
##########                            Init                            ##########
################################################################################
# Packages
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(sctransform)
  library(DropletUtils)
  library(Matrix)
  library(mclust)
  library(kableExtra)
  library(intrinsicDimension)
})

# Shared functions
source('code/shared.R')

# Paths
project.path <- file.path(
  snakemake@config[['project_path']],
  snakemake@config[['scop_id']])

com.path <- file.path(
  project.path,
  'output',
  snakemake@config[['com_id']],
  'COMUNEQAID_v2.0')

aggr.path <- file.path(
  com.path,
  'aggregated')

aggr.stats.path <- file.path(
  project.path,
  snakemake@config[['out_path']],
  snakemake@config[['com_id']],
  'metadata',
  'aggregated')

# Files
log.file <- file.path(
  aggr.path,
  paste0(date.and.time,'.log'))

unfiltered.file <- file.path(
  aggr.path,
  'seurat_unfiltered.rds')

filtered.file <- file.path(
  aggr.path,
  'seurat_filtered.rds')

unfiltered.velo.file <- file.path(
  aggr.path,
  'seurat_velocity_unfiltered.rds')

filtered.velo.file <- file.path(
  aggr.path,
  'seurat_velocity_filtered.rds')

dir.create(aggr.stats.path, recursive = T, showWarnings = F)

# Metadata
suppressMessages({
  rnx.sheet <-
    as_tibble(snakemake@config[['reaction_sheet']])
})

if (length(rnx.sheet[['reaction_id']]) == 1) {
  message <- paste0('################################################################################\n',
                    'Attempted to aggregate single reaction - skipping process\n',
                    '################################################################################\n')
  
  # Write the message to a text file
  writeLines(message,
             file.path(
               aggr.path,
               paste0('_skipped_pipestance_',date.and.time,'.log')))
  cat(message)
}

if (length(rnx.sheet[['reaction_id']]) > 1) {
  
  if (file.exists(unfiltered.file) &
      file.exists(filtered.file)){
    message <- paste0('################################################################################\n',
                      'Seurat objects already exists - skipping pipestance\n',
                      '################################################################################\n')
    
    # Write the message to a text file
    writeLines(message,
               file.path(
                 aggr.path,
                 paste0('_skipped_pipestance_',date.and.time,'.log')))
    cat(message)
  }
  
  if (!file.exists(unfiltered.file) &
      !file.exists(filtered.file)){
    
    dir.create(aggr.path, recursive = T, showWarnings = F)
    
    ################################################################################
    ##########                 Preparing aggregated output                ##########
    ################################################################################
    
    sink.file <- file(log.file, 'w')
    sink(file = sink.file,
         append = T,
         type = 'output',
         split = T)
    
    # Startup message
    cat(rep('#',80),'\n',
        '#####                                                                      #####\n',
        '#####                     Aggregating Seurat objects (aggr)                #####\n',
        '#####                                                                      #####\n',
        '#####                                             5/6                      #####\n',
        '#####           ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▒▒▒▒▒▒▒▒░░░░░░░░            #####\n',
        '#####                                                                      #####\n',
        rep('#',80),'\n',
        '#####\n',
        '##\n',
        '##\tSCOP ID:      \t\t',snakemake@config[['scop_id']],'\n',
        '##\tCOMUNEQAID ID:\t\t',snakemake@config[['com_id']],'\n',
        '##\n',
        '##\n',
        '##\t',date.and.time.pretty,'\n',
        '##\n',
        '##\n',
        '###\n',
        rep('#',80),'\n',
        '#\n',
        sep = '')
    
    cat('#\treading seurat objects..\n',
        sep = '')
    
    seur.unfiltered.list <- lapply(rnx.sheet[['reaction_id']], function(rnx) {
      q.protocol <- rnx.sheet[rnx.sheet[['reaction_id']] == rnx,][['seq_type']]
      tmp.path <- file.path(com.path, rnx)
      if (q.protocol == 'cells') {
        tmp.file <- file.path(
          tmp.path,
          'seurat_unfiltered.rds')
      }
      if (q.protocol == 'nuclei') {
        tmp.file <- file.path(
          tmp.path,
          'seurat_velocity_unfiltered.rds')
      }
      cat('#\t..\t\t',rnx,' (unfiltered)\n')
      readRDS(tmp.file)
      })
    
    seur.filtered.list <- lapply(rnx.sheet[['reaction_id']], function(rnx) {
      q.protocol <- rnx.sheet[rnx.sheet[['reaction_id']] == rnx,][['seq_type']]
      tmp.path <- file.path(com.path, rnx)
      if (q.protocol == 'cells') {
        tmp.file <- file.path(
          tmp.path,
          'seurat_filtered.rds')
      }
      if (q.protocol == 'nuclei') {
        tmp.file <- file.path(
          tmp.path,
          'seurat_velocity_filtered.rds')
      }
      cat('#\t..\t\t',rnx,' (filtered)\n')
      readRDS(tmp.file)
    })
    cat('#\t..\n',
        '#\t..done\n',
        sep = '')
    
    cat('#\n',
        '#\tremoving hto assays..\n',
        sep = '')
    for (i in seq(seur.filtered.list)) {
      seur.unfiltered.list[[i]][['HTO']] <- NULL
      seur.filtered.list[[i]][['HTO']] <- NULL
    }
    cat('#\t..\n',
        '#\t..done\n',
        sep = '')
    
    cat('#\n',
        '#\tmerging objects..\n',
        sep = '')
    # Merge objects
    suppressWarnings({
      seur.unfiltered.comb <- merge(seur.unfiltered.list[[1]],
                                    seur.unfiltered.list[2:length(seur.unfiltered.list)])
      seur.filtered.comb <- merge(seur.filtered.list[[1]],
                                  seur.filtered.list[2:length(seur.filtered.list)])
    })
    cat('#\t..\n',
        '#\t..done\n',
        sep = '')
    
    cat('#\n',
        '#\tprocessing unfiltered object..\n',
        sep = '')
    seur.unfiltered.comb <- process_seur(seur.unfiltered.comb)
    seur.unfiltered.comb[['SCT_snn_res.0.8']] <- NULL
    
    seur.unfiltered.comb[['SCT_tsne1']] <- seur.unfiltered.comb@reductions$tsne@cell.embeddings[,1]
    seur.unfiltered.comb[['SCT_tsne2']] <- seur.unfiltered.comb@reductions$tsne@cell.embeddings[,2]
    seur.unfiltered.comb[['SCT_umap1']] <- seur.unfiltered.comb@reductions$umap@cell.embeddings[,1]
    seur.unfiltered.comb[['SCT_umap2']] <- seur.unfiltered.comb@reductions$umap@cell.embeddings[,2]
    
    cat('#\t..\n',
        '#\t..done\n',
        sep = '')
    
    cat('#\n',
        '#\tprocessing filtered object..\n',
        sep = '')
    seur.filtered.comb <- process_seur(seur.filtered.comb)
    seur.filtered.comb[['SCT_snn_res.0.8']] <- NULL
    
    seur.filtered.comb[['SCT_tsne1']] <- seur.filtered.comb@reductions$tsne@cell.embeddings[,1]
    seur.filtered.comb[['SCT_tsne2']] <- seur.filtered.comb@reductions$tsne@cell.embeddings[,2]
    seur.filtered.comb[['SCT_umap1']] <- seur.filtered.comb@reductions$umap@cell.embeddings[,1]
    seur.filtered.comb[['SCT_umap2']] <- seur.filtered.comb@reductions$umap@cell.embeddings[,2]
    
    cat('#\t..\n',
        '#\t..done\n',
        sep = '')
    
    cat('#\n',
        '#\twriting metadata..\n',
        sep = '')
    cat('#\t..\t\"seurat_unfiltered_metadata.csv\"\n')
    write.table(x = seur.unfiltered.comb[[]],
                file = file.path(aggr.stats.path,'seurat_unfiltered_metadata.csv'),
                sep = ',',
                row.names = T)
    
    seur.unfiltered.comb[['SCT_tsne1']] <- NULL
    seur.unfiltered.comb[['SCT_tsne2']] <- NULL
    seur.unfiltered.comb[['SCT_umap1']] <- NULL
    seur.unfiltered.comb[['SCT_umap2']] <- NULL
    
    cat('#\t..\t\"seurat_filtered_metadata.csv\"\n')
    write.table(x = seur.filtered.comb[[]],
                file = file.path(aggr.stats.path,'seurat_filtered_metadata.csv'),
                sep = ',',
                row.names = T)
    
    seur.filtered.comb[['SCT_tsne1']] <- NULL
    seur.filtered.comb[['SCT_tsne2']] <- NULL
    seur.filtered.comb[['SCT_umap1']] <- NULL
    seur.filtered.comb[['SCT_umap2']] <- NULL
    
    cat('#\t..\n',
        '#\t..done\n',
        sep = '')
    
    cat('#\n',
        '#\twriting seurat objects..\n',
        sep = '')
    

    if (!is.null(seur.unfiltered.comb[['spliced']])) {
      
      cat('#\t..\t\"seurat_velocity_unfiltered.rds\"\n')
      saveRDS(seur.unfiltered.comb, unfiltered.velo.file)
      
      seur.unfiltered.comb[['spliced']] <- NULL
      seur.unfiltered.comb[['unspliced']] <- NULL
      
      cat('#\t..\t\"seurat_velocity_filtered.rds\"\n')
      saveRDS(seur.filtered.comb, filtered.velo.file)
      
      seur.filtered.comb[['spliced']] <- NULL
      seur.filtered.comb[['unspliced']] <- NULL
    }
    
    cat('#\t..\t\"seurat_unfiltered.rds\"\n')
    saveRDS(seur.unfiltered.comb, unfiltered.file)
    
    cat('#\t..\t\"seurat_filtered.rds\"\n')
    saveRDS(seur.filtered.comb, filtered.file)
    
    cat('#\t..\n',
        '#\t..done\n',
        '#\n',
        '#\n',
        '###\n',
        rep('#',80),'\n',
        sep = '')
    
  }
}
sink()
file.create(snakemake@output[[1]])