################################################################################
##########                            Init                            ##########
################################################################################
# Packages
suppressPackageStartupMessages({
  library(foreach)
  library(parallel)
  library(doParallel)
  library(Seurat)
  library(tximport)       
  library(SingleCellExperiment)
  library(future)
  library(tidyverse)
  library(scater)
  library(scran)
  library(sctransform)
  library(DropletUtils)
  library(Matrix)
  library(mclust)
  library(data.table)
  library(ggrepel)
  library(rjson)
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
  snakemake@config[['pipeline_version']])

# Metadata
suppressMessages({
  rnx2lib.libsheet.lib2seq.seqsheet <-
    full_join(as_tibble(snakemake@config[['reaction2library']]),
              as_tibble(snakemake@config[['library_sheet']])) %>%
    full_join(as_tibble(snakemake@config[['library2sequencing']])) %>%
    full_join(as_tibble(snakemake@config[['sequencing_sheet']]))
  rnx.sheet <-
    as_tibble(snakemake@config[['reaction_sheet']])
})


################################################################################
##########       Processing single nucleus RNA sequencing data        ##########
################################################################################

####################  Iterate over reactions  ####################
cl <- makeCluster(length(row_number(rnx.sheet)))
registerDoParallel(cl)

foreach(q.rnx = rnx.sheet[["reaction_id"]],
        .combine = 'c',
        .packages = c('Seurat','tximport','SingleCellExperiment','future',
                      'tidyverse','scater','scran','sctransform',
                      'DropletUtils','Matrix','mclust','intrinsicDimension')
) %dopar% {

#for (q.rnx in rnx.sheet[["reaction_id"]]) {
  
  ####################  Prep processing of rnx i  ####################
  # Reaction vars
  rnx.tib <- filter(rnx.sheet, reaction_id == q.rnx)
  
  q.protocol <- rnx.tib[['seq_type']]
  q.cutoff <- rnx.tib[['cutoff']]
  q.organism <- rnx.tib[['align']]
  
  q.low <- filter(as_tibble(snakemake@config[['reaction_sheet']]), reaction_id == q.rnx)[['quantile_low']][[1]]
  q.high <- filter(as_tibble(snakemake@config[['reaction_sheet']]), reaction_id == q.rnx)[['quantile_high']][[1]]
  
  if (is.character(q.low) & is.character(q.high)) {
    q.low <- as.numeric(unlist(strsplit(q.low, " ")))
    q.high <- as.numeric(unlist(strsplit(q.high, " ")))
  }
  
  if (is.null(q.low) & is.null(q.high)) {
    q.low <- 1
    q.high <- 0.01
  }
  
  rnx.path <- file.path(
    com.path,
    q.rnx)
  
  mat.stats.path <- file.path(
    project.path,
    snakemake@config[['out_path']],
    snakemake@config[['com_id']],
    'metadata',
    q.rnx)
  
  log.file <- file.path(
    rnx.path,
    paste0(date.and.time,'.log'))
  
  unfiltered.file <- file.path(
    rnx.path,
    'seurat_unfiltered.rds')
  
  filtered.file <- file.path(
    rnx.path,
    'seurat_filtered.rds')
  
  unfiltered.velo.file <- file.path(
    rnx.path,
    'seurat_velocity_unfiltered.rds')
  
  filtered.velo.file <- file.path(
    rnx.path,
    'seurat_velocity_filtered.rds')
  
  if (file.exists(unfiltered.file) &
      file.exists(filtered.file)) {
    
    message <- paste0('################################################################################\n',
                      'Seurat objects already exists - skipping reaction_id: ',q.rnx,'\n',
                      '################################################################################\n')
    
    # Write the message to a text file
    writeLines(message,
               file.path(
                 rnx.path,
                 paste0('_skipped_pipestance_',date.and.time,'.log')))
    cat(message)
  }
  
  if (!file.exists(unfiltered.file) &
      !file.exists(filtered.file)) {
    dir.create(rnx.path, recursive = T, showWarnings = F)
    dir.create(mat.stats.path, recursive = T, showWarnings = F)
    
    # Library vars
    lib.tib <- filter(rnx2lib.libsheet.lib2seq.seqsheet,
                      reaction_id == q.rnx)
    
    q.library.10x <- lib.tib %>% 
      filter(library_type == '10x') %>% 
      select(library_id) %>% 
      unlist()
    q.bcl.10x <- lib.tib %>% 
      filter(library_type == '10x') %>% 
      select(bcl_folder) %>% 
      unlist()
    
    mat.files.10x <- file.path(project.path,snakemake@config[['out_path']],snakemake@config[['com_id']],salmon.version.alevin.fry.version,q.rnx,q.library.10x,'res')
    if (!file.exists(mat.files.10x)) { stop('Not able to locate quants_mat.gz (gene expression)')}
    
    if ('hto' %in% lib.tib[['library_type']]) {
      q.library.hto <- lib.tib %>% 
        filter(library_type == 'hto') %>% 
        select(library_id) %>% 
        unlist()
      q.bcl.hto <- lib.tib %>% 
        filter(library_type == 'hto') %>% 
        select(bcl_folder) %>% 
        unlist()
      
      mat.files.hto <- file.path(project.path,snakemake@config[['out_path']],snakemake@config[['com_id']],salmon.version.alevin.fry.version,q.rnx,q.library.hto,'res')
      if (!file.exists(mat.files.hto)) { stop('Not able to locate quants_mat.gz (HTO)')}
    }
    
    sink.file <- file(log.file, 'w')
    sink(file = sink.file,
         append = T,
         type = 'output',
         split = T)
    
    # Startup message
    cat(rep('#',80),'\n',
        '#####                                                                      #####\n',
        '#####                        Making Seurat object (rnx)                    #####\n',
        '#####                                                                      #####\n',
        '#####                                     4/6                              #####\n',
        '#####           ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▒▒▒▒▒▒▒▒░░░░░░░░░░░░░░░░            #####\n',
        '#####                                                                      #####\n',
        rep('#',80),'\n',
        '#####\n',
        '##\n',
        '##\tSCOP ID:      \t\t',snakemake@config[['scop_id']],'\n',
        '##\tCOMUNEQAID ID:\t\t',snakemake@config[['com_id']],'\n',
        '##\n',
        '##\treaction ID:  \t\t- ',q.rnx,'\n',
        '##\torganism:     \t\t- ',q.organism,'\n',
        '##\tprotocol:     \t\t- ',q.protocol,'\n',
        '##\tcutoff:       \t\t- ',q.cutoff,'\n',
        '##\n',
        '##\n',
        '##\t',date.and.time.pretty,'\n',
        '##\n',
        '##\n',
        '###\n',
        rep('#',80),'\n',
        sep = '')
    
    
    ####################  Read RNA matrices ####################
    cat('#\n',
        '#\tReading count matrices..\n',
        sep = '')
    
    cat('#\t..\t',q.library.10x,'..\n',
        sep = '')
    if (q.protocol == 'nuclei') {
      counts.all <- load_fry(frydir = mat.files.10x, output_list = T, verbose = T)
      counts.rna <- counts.all[['RNA']]
      counts.spl <- counts.all[['spliced']]
      counts.uns <- counts.all[['unspliced']]
    }
    if (q.protocol == 'cells') {
      counts.rna <- counts.all <- load_fry(frydir = mat.files.10x)
    }
    
    ####################  Read HTO matrices ####################
    if ('hto' %in% lib.tib[['library_type']]) {
      cat('#\t..\n',
          '#\t..\t',q.library.hto,'..\n',
          sep = '')
      counts.hto <- load_fry(frydir = mat.files.hto)
      rownames(counts.hto) <- gsub("_", "-", rownames(counts.hto))
      
      # Saving col sums as metadata
      cat('#\t..\twriting metadata..\n')
      cat('#\t..\t\t\"hto_cell-barcode_stats.csv\"\n',
          sep = '')
      col.sums.hto <- Matrix::colSums(counts.hto)
      hto.tib <- tibble(
        Barcode = names(col.sums.hto),
        nUMI_HTO = col.sums.hto)
      
      write.table(x = hto.tib,
                  file = file.path(mat.stats.path,'hto_cell-barcode_stats.csv'),
                  sep = ',',
                  row.names = F)
    }
    
    cat('#\t..\n',
        '#\t..done\n',
        '#\n',
        sep = '')
    
    
    ####################  Cell calling ####################
    cat('#\tcalling cell barcodes based on RNA library..\n',
        '#\t..\n',
        sep = '')
    
    call.droplets.list <- call_droplets(counts = counts.all,
                                        protocol = q.protocol,
                                        cutoff = q.cutoff)
    
    valid.bcs <- call.droplets.list[[1]]
    bc.df <- call.droplets.list[[2]]
    bc.annot <- call.droplets.list[[3]]
    
    if ('hto' %in% lib.tib[['library_type']]) {
      valid.bcs <- intersect(valid.bcs,colnames(counts.hto))
    }
    
    ####################  Filter matrices  ####################
    cat('#\t..\n',
        '#\t..\tsubsetting matricies..\n',
        sep = '')
    
    cat('#\t..\t\t',q.library.10x,'..\n',
        sep = '')
    cat('#\t..\t\t\tRNA counts (USA)..\n')
    counts.rna <- counts.rna[,valid.bcs]
    if (q.protocol == 'nuclei') {
      cat('#\t..\t\t\tSpliced (& Ambiguous) counts (SA)..\n')
      counts.spl <- counts.spl[,valid.bcs]
      cat('#\t..\t\t\tUnspliced counts (U)..\n')
      counts.uns <- counts.uns[,valid.bcs]
    }
    if ('hto' %in% lib.tib[['library_type']]) {
      cat('#\t..\n',
          '#\t..\t\t',q.library.hto,'..\n',
          sep = '')
      cat('#\t..\t\t\tHTO counts..\n')
      counts.hto <- counts.hto[,valid.bcs]
    }
    
    
    ####################  Write barcode rank info as metadata  ####################
    cat('#\t..\n',
        '#\t..\twriting metadata..\n',
        sep = '')
    cat('#\t..\t\t\"rna_cell-barcode_stats.csv\"\n',
        sep = '')
    write.table(x = bc.df,
                file = file.path(mat.stats.path,'rna_cell-barcode_stats.csv'),
                sep = ',',
                row.names = F)
    cat('#\t..\t\t\"rna_cell-barcode_calling_metadata.csv\"\n',
        sep = '')
    write.table(x = bc.annot,
                file = file.path(mat.stats.path,'rna_cell-barcode_calling_metadata.csv'),
                sep = ',',
                row.names = F)
    cat('#\t..\n',
        '#\t..done\n',
        '#\n',
        sep = '')
    
    
    #################### Create Seurat  ####################
    cat('#\tcreate seurat object..\n',
        '#\t..\tfrom RNA assay..\n',
        sep = '')
    seur.full <- CreateSeuratObject(counts.rna, project = q.rnx)
    cat('#\t..\t\tname:\t\"RNA\"\n',
        sep = '')
    seur.full[['reactionID']] <- q.rnx
    
    if (q.protocol == 'nuclei') {
      cat('#\t..\tadding spliced assay..\n',
          sep = '')
      seur.full[['spliced']] <- CreateAssayObject(counts.spl)
      cat('#\t..\t\tname:\t\"spliced\"\n',
          sep = '')
      cat('#\t..\tadding unspliced assay..\n',
          sep = '')
      seur.full[['unspliced']] <- CreateAssayObject(counts.uns)
      cat('#\t..\t\tname:\t\"unspliced\"\n',
          sep = '')
    }
    
    #################### Demultiplex HTO library  ####################
    if ('hto' %in% lib.tib[['library_type']]) {
      cat('#\t..\tadding HTO assay..\n',
          sep = '')
      seur.full[['HTO']] <- CreateAssayObject(counts.hto)
      cat('#\t..\t\tname:\t\"HTO\"\n',
          sep = '')
      
      cat('#\t..\n',
          '#\t..done\n',
          '#\n',
          sep = '')
      
      cat('#\tdemultiplexing hto library..\n')
      
      cat('#\t..\tclr-normalizing hto assay..\n')
      seur.full <- NormalizeData(seur.full, assay = 'HTO', normalization.method = 'CLR', margin = 2, verbose = F)
      
      cat('#\t..\tdemultiplexing HTOs..\n')
      
      # Setting the threshold based on the quantiles of the negative and the positive clusters.
      cat('#\t..\t\tq.low = ',q.low,'\n',
          '#\t..\t\tq.high = ',q.high,'\n',
          sep = '')
      seur.full <- HTODemux.mcl(seur.full, q_l = q.low, q_h = q.high)
      
      cat('#\t..\t\tdoublets:\t',table(seur.full[['HTO_globalClass']])[names(table(seur.full[['HTO_globalClass']])) == 'Doublet'],'\n',
          '#\t..\t\tnegatives:\t',table(seur.full[['HTO_globalClass']])[names(table(seur.full[['HTO_globalClass']])) == 'Negative'],'\n',
          '#\t..\t\tsinglets:\t',table(seur.full[['HTO_globalClass']])[names(table(seur.full[['HTO_globalClass']])) == 'Singlet'],'\n',
          sep = '')
      
      tmp.data <- GetAssayData(object = seur.full, assay = 'HTO', slot = 'data')
      hto.data.wide = data.frame(t(tmp.data))
      colnames(hto.data.wide) <- gsub("[.]", "-", colnames(hto.data.wide), perl = T)
      hto.data.long = data.table::melt(data.table::setDT(hto.data.wide,keep.rownames = T), id.vars = 'rn', value.name = 'expression',variable.name = 'sample_id')
      hto.data.long <- rename(hto.data.long, 'CB' = 'rn')
      hto.data.long[, sample_id := as.character(sample_id)]
      
      if (length(q.low) == 1 && length(q.high) == 1) {
        hto_mcl.cutoff = data.table::data.table(cut_off = future.apply::future_apply(tmp.data,1,function(x) select_hash_cutoff_mcl(x, q_l = q.low, q_h = q.high), future.seed = T), hto_name = colnames(hto.data.wide)[-1], key = "hto_name")
      }
      
      if (length(q.low) > 1 && length(q.high) > 1) {
        feature_indices <- seq_len(nrow(tmp.data))
        
        hto.cutoff.metadata <- data.frame('cut_off' = numeric())
        
        for (idx in feature_indices) {
          hto.cutoff.metadata[idx,] <- select_hash_cutoff_mcl(tmp.data[idx,], q.low[idx], q.high[idx])
        }
        
        hto_mcl.cutoff = data.table::data.table(hto.cutoff.metadata, hto_name = colnames(hto.data.wide)[-1], key = "hto_name")
      }
      
      hto.data.long[, above_cutoff := expression > hto_mcl.cutoff[sample_id, cut_off]]
      
      
      cat('#\t..\n')
      seur.full <- process_seur(seur.full,assay = 'HTO')
      
      cat('#\t..\n',
          '#\t..\twriting metadata..\n',
          sep = '')
      cat('#\t..\t\t\"hto_sample_calling_metadata.csv\"\n')
      write.table(x = hto.data.long,
                  file = file.path(mat.stats.path,'hto_sample_calling_metadata.csv'),
                  sep = ',',
                  row.names = F)
      cat('#\t..\t\t\"hto_cutoff_metadata.csv\"\n')
      write.table(x = hto_mcl.cutoff,
                  file = file.path(mat.stats.path,'hto_cutoff_metadata.csv'),
                  sep = ',',
                  row.names = F)
      cat('#\t..\n',
          '#\t..done\n',
          '#\n',
          sep = '')
      
      # identify intra-hash doublets
      seur.full[['HTO_doubletBool']] <- if_else(seur.full[['sampleID']] == 'Doublet', true = T, false = F)
      
      if (TRUE %in% seur.full[[]][['HTO_doubletBool']]) {
        cat('#\tidentifying intra-HTO doublets..\n',
            '#\n',
            sep = '')
        
        cat('#\t..\tconverting to SingleCellExperiment..\n')
        sce.full <- as.SingleCellExperiment(seur.full)
        
        cat('#\t..\tlog-normalizing rna assay..\n')
        sce.full <- logNormCounts(sce.full)
        
        cat('#\t..\tmodeling per-gene variance..\n')
        dec.hash <- modelGeneVar(sce.full)
        
        cat('#\t..\trunning PCA on top 1000 variable genes..\n')
        top.hash <- getTopHVGs(dec.hash, n = 1000)
        sce.full <- runPCA(sce.full, subset_row = top.hash, ncomponents = 20)
        
        # Recovering the intra-sample doublets:
        cat('#\t..\trecovering intra-sample doublets..\n')
        hashed.doublets <- scDblFinder::recoverDoublets(sce.full,
                                                        use.dimred = 'PCA',
                                                        doublets = sce.full[['HTO_doubletBool']],
                                                        samples = table(sce.full[['sampleID']]))
        
        sce.full[['RNA_doubletNeighborProportion']] <- hashed.doublets[['proportion']]
        sce.full[['RNA_recoveredDoubletBool']] <- hashed.doublets[['predicted']]
        
        cat('#\t..\trecovering doublet neighbors..\n')
        sce.full[['RNA_doubletNeighborBool']] <- dub_cutoff(sce.full)
        
        seur.full[['HTO_doubletBool']] <- sce.full[['HTO_doubletBool']]
        seur.full[['RNA_recoveredDoubletBool']] <- sce.full[['RNA_recoveredDoubletBool']]
        seur.full[['RNA_doubletNeighborBool']] <- sce.full[['RNA_doubletNeighborBool']]
        seur.full[['RNA_doubletNeighborProportion']] <- sce.full[['RNA_doubletNeighborProportion']]
        
        cat('#\t..\n',
            '#\t..done\n',
            '#\n',
            sep = '')
      }
        
      if (TRUE %!in% seur.full[[]][['HTO_doubletBool']]) {
        seur.full[['RNA_recoveredDoubletBool']] <- FALSE
        seur.full[['RNA_doubletNeighborBool']] <- FALSE
        seur.full[['RNA_doubletNeighborProportion']] <- 0
      }
    }
    
    if ('hto' %!in% lib.tib[['library_type']]) {
      cat('#\t..\n',
          '#\t..done\n',
          '#\n',
          sep = '')
      
      sce.full <- as.SingleCellExperiment(seur.full)
      
      cat('#\t..\tidentifying intra-hash doublets\n')
      sce.full <- logNormCounts(sce.full)
      top.rna <- getTopHVGs(modelGeneVar(sce.full))
      sce.full <- runPCA(sce.full, subset_row = top.rna, n = 1000, ncomponents = 20)
      
      dbl.dens <- scDblFinder::computeDoubletDensity(sce.full, subset.row = top.rna,
                                                     d = ncol(reducedDim(sce.full)))
      
      seur.full[['DoubletScore']] <- dbl.dens
    }
    
    cat('#\toutput..\n')
    
    seur.full <- process_seur(seur.full, assay = 'RNA')
    
    seur.full[['SCT_tsne1']] <- seur.full@reductions$tsne@cell.embeddings[,1]
    seur.full[['SCT_tsne2']] <- seur.full@reductions$tsne@cell.embeddings[,2]
    seur.full[['SCT_umap1']] <- seur.full@reductions$umap@cell.embeddings[,1]
    seur.full[['SCT_umap2']] <- seur.full@reductions$umap@cell.embeddings[,2]
    
    if ('hto' %in% lib.tib[['library_type']]) {
      seur.full[['HTO_tsne1']] <- seur.full@reductions$hto.tsne@cell.embeddings[,1]
      seur.full[['HTO_tsne2']] <- seur.full@reductions$hto.tsne@cell.embeddings[,2]
      seur.full[['HTO_umap1']] <- seur.full@reductions$hto.umap@cell.embeddings[,1]
      seur.full[['HTO_umap2']] <- seur.full@reductions$hto.umap@cell.embeddings[,2]
      
      cat('#\t..\n',
          '#\t..\tsubsetting HTO/intra-HTO singlets..\n',
          sep = '')
      seur.sing <- subset(
        seur.full,
        HTO_globalClass == 'Singlet' &
          RNA_recoveredDoubletBool == F &
          RNA_doubletNeighborBool == F)
      
      
      
      cat('#\t..\t\tremoved HTO negatives:\t\t',table(seur.full[['HTO_globalClass']])['Negative'],'\n',
          '#\t..\t\tremoved HTO doublets:\t\t',table(seur.full[['HTO_globalClass']])['Doublet'],'\n',
          '#\t..\t\tremoved intra-HTO doublets:\t',table(seur.full[['RNA_doubletNeighborBool']])['TRUE'],'\n',
          '#\t..\t\tretaining ',length(colnames(seur.sing)),' ',q.protocol,' \n',
          '#\t..\n',
          sep = '')
      
      seur.sing <- process_seur(seur.sing, assay = 'RNA')
      
      seur.sing[['SCT_tsne1']] <- seur.sing@reductions$tsne@cell.embeddings[,1]
      seur.sing[['SCT_tsne2']] <- seur.sing@reductions$tsne@cell.embeddings[,2]
      seur.sing[['SCT_umap1']] <- seur.sing@reductions$umap@cell.embeddings[,1]
      seur.sing[['SCT_umap2']] <- seur.sing@reductions$umap@cell.embeddings[,2]
      
      seur.sing@meta.data <- seur.sing[[c('reactionID','sampleID','nCount_RNA','nFeature_RNA','nCount_HTO','nFeature_HTO','nCount_spliced','nFeature_spliced','nCount_unspliced','nFeature_unspliced','nCount_SCT','nFeature_SCT','seurat_clusters','SCT_tsne1','SCT_tsne2','SCT_umap1','SCT_umap2','HTO_tsne1','HTO_tsne2','HTO_umap1','HTO_umap2')]]
    }
    
    seur.full@meta.data <- seur.full[[c('reactionID','sampleID','nCount_RNA','nFeature_RNA','nCount_HTO','nFeature_HTO','nCount_spliced','nFeature_spliced','nCount_unspliced','nFeature_unspliced','nCount_SCT','nFeature_SCT','HTO_primaryID','HTO_secondaryID','HTO_margin','HTO_calledFeatures','HTO_globalClass','HTO_doubletBool','RNA_recoveredDoubletBool','RNA_doubletNeighborBool','RNA_doubletNeighborProportion','seurat_clusters','SCT_tsne1','SCT_tsne2','SCT_umap1','SCT_umap2','HTO_tsne1','HTO_tsne2','HTO_umap1','HTO_umap2')]]
    
    cat('#\t..\n',
        '#\t..\twriting metadata..\n',
        sep = '')
    cat('#\t..\t\t\"seurat_unfiltered_metadata.csv\"\n')
    write.table(x = seur.full[[]],
                file = file.path(mat.stats.path,'seurat_unfiltered_metadata.csv'),
                sep = ',',
                row.names = T)
    cat('#\t..\t\t\"seurat_filtered_metadata.csv\"\n')
    write.table(x = seur.sing[[]],
                file = file.path(mat.stats.path,'seurat_filtered_metadata.csv'),
                sep = ',',
                row.names = T)
    
    seur.full@meta.data <- seur.full[[c('reactionID','sampleID','nCount_RNA','nFeature_RNA','nCount_HTO','nFeature_HTO','nCount_spliced','nFeature_spliced','nCount_unspliced','nFeature_unspliced','nCount_SCT','nFeature_SCT','HTO_primaryID','HTO_secondaryID','HTO_margin','HTO_calledFeatures','HTO_globalClass','HTO_doubletBool','RNA_recoveredDoubletBool','RNA_doubletNeighborBool','RNA_doubletNeighborProportion','seurat_clusters')]]
    seur.sing@meta.data <- seur.sing[[c('reactionID','sampleID','nCount_RNA','nFeature_RNA','nCount_HTO','nFeature_HTO','nCount_spliced','nFeature_spliced','nCount_unspliced','nFeature_unspliced','nCount_SCT','nFeature_SCT','seurat_clusters')]]
    
    cat('#\t..\n',
        '#\t..\twriting seurat objects..\n',
        sep = '')
    if (q.protocol == 'nuclei') {
      cat('#\t..\t\t\"seurat_velocity_unfiltered.rds\"\n')
      saveRDS(seur.full, unfiltered.velo.file)
      cat('#\t..\t\t\"seurat_velocity_filtered.rds\"\n')
      saveRDS(seur.sing, filtered.velo.file)
      
      seur.full[['spliced']] <- NULL
      seur.full[['unspliced']] <- NULL
      seur.full@meta.data <- seur.full[[c('reactionID','sampleID','nCount_RNA','nFeature_RNA','nCount_HTO','nFeature_HTO','nCount_SCT','nFeature_SCT','HTO_primaryID','HTO_secondaryID','HTO_margin','HTO_calledFeatures','HTO_globalClass','HTO_doubletBool','RNA_recoveredDoubletBool','RNA_doubletNeighborBool','RNA_doubletNeighborProportion','seurat_clusters')]]
      
      seur.sing[['spliced']] <- NULL
      seur.sing[['unspliced']] <- NULL
      seur.sing@meta.data <- seur.sing[[c('reactionID','sampleID','nCount_RNA','nFeature_RNA','nCount_HTO','nFeature_HTO','nCount_SCT','nFeature_SCT','seurat_clusters')]]
      
    }
    
    cat('#\t..\t\t\"seurat_unfiltered.rds\"\n')
    saveRDS(seur.full, unfiltered.file)
    cat('#\t..\t\t\"seurat_filtered.rds\"\n')
    saveRDS(seur.sing, filtered.file)
    
    cat('#\t..\n',
        '#\t..done\n',
        sep = '')
    
    cat('#\n',
        '#\n',
        '###\n',
        rep('#',80),'\n',
        sep = '')
    
    sink()
  }
  NULL
}
file.create(snakemake@output[[1]])
