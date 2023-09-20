set.seed(102)

# Date and time
current.time <- Sys.time()
date.and.time <- format(current.time, "%Y-%m-%d_%H-%M-%S")
date.and.time.pretty <- format(current.time, "%Y-%m-%d %H:%M:%S")

# # Software versions
get_software_version <- function(software) {
  commands <- list(
    salmon = "salmon --version",
    `alevin-fry` = "alevin-fry --version",
    `bcl-convert` = "bcl-convert --version"
  )
  
  if (!software %in% names(commands)) {
    stop(paste("Unsupported software:", software))
  }
  
  command_with_redirect <- paste(commands[[software]], "2>&1")
  output <- system(command_with_redirect, intern = T)
  
  # For bcl-convert, assuming it's a two-line output and version is on the first line
  if (software == "bcl-convert") {
    words <- unlist(strsplit(output[1], " "))
  } else {
    words <- unlist(strsplit(output, " "))
  }
  
  software_name <- tolower(words[1])
  
  if (software == "bcl-convert") {
    version_full <- words[3]
    version_segments <- unlist(strsplit(version_full, "\\."))
    
    version <- paste0("v", paste(tail(version_segments, 3), collapse="."))
  } else {
    version <- paste0("v", words[2])
  }
  
  return(paste0(software_name, "_", version))
}

bcl.convert.version <- get_software_version('bcl-convert')
salmon.version = get_software_version('salmon')
alevin.fry.version = get_software_version('alevin-fry')
salmon.version.alevin.fry.version <- paste0(salmon.version, '_', alevin.fry.version)


# Data processing
load_fry <- function(frydir, which_counts = c('U','S','A'), verbose = FALSE, output_list = F) {
  # read in metadata
  cat('#\t..\t\tread in metadata..\n')
  meta_info = rjson::fromJSON(file = file.path(frydir, 'quant.json'))
  ng = meta_info[['num_genes']]
  usa_mode = meta_info[['usa_mode']]
  
  if (usa_mode) {
    if (length(which_counts) == 0) {
      stop('Please at least provide one status in \'U\' \'S\' \'A\' ')
    }
    if (verbose) {
      #message('processing input in USA mode, will return ', paste(which_counts, collapse = '+'))
      cat('#\t..\t\t','processing input in USA mode, will return ', paste(which_counts,collapse = '+'),'..\n',sep = '')
      
    }
  } else if (verbose) {
    #message('processing input in standard mode, will return spliced count')
    cat('#\t..\t\tprocessing input in standard mode, will return spliced count..\n')
  }
  
  # read in count matrix
  cat('#\t..\t\tread in count matrix..\n')
  af_raw = readMM(file = file.path(frydir, 'alevin', 'quants_mat.mtx'))
  
  # if usa mode, each gene gets 3 rows, so ng/3
  if (usa_mode) {
    ng = as.integer(ng/3)
  }
  
  # read in gene name file and cell barcode file
  cat('#\t..\t\tread in gene name file..\n')
  afg = read.csv(file.path(frydir, 'alevin', 'quants_mat_cols.txt'), strip.white = TRUE, header = FALSE, nrows = ng, col.names = c('gene_ids'))
  cat('#\t..\t\tread in cell barcode file..\n')
  afc = read.csv(file.path(frydir, 'alevin', 'quants_mat_rows.txt'), strip.white = TRUE, header = FALSE, col.names = c('barcodes'))
  
  # if in usa_mode, sum up counts in different status according to which_counts
  if (output_list) {
    which_counts = c('U','S','A')
    cat('#\t..\t\tsum up counts according to status affiliation..\n')
    if (usa_mode) {
      rd = list('S' = seq(1, ng), 'U' =  seq(ng + 1, 2*ng), 'A' =  seq(2*ng + 1, 3*ng))
      o = af_raw[, rd[[which_counts[1]]]]
      for (wc in which_counts[-1]) {
        o = o + af_raw[, rd[[wc]]]
      }
    } else {
      o = af_raw
    }
    
    res.rna <- t(o)
    colnames(res.rna) <- afc[['barcodes']]
    rownames(res.rna) <- afg[['gene_ids']]
    
    which_counts = c('S','A')
    if (usa_mode) {
      rd = list('S' = seq(1, ng), 'U' =  seq(ng + 1, 2*ng), 'A' =  seq(2*ng + 1, 3*ng))
      o = af_raw[, rd[[which_counts[1]]]]
      for (wc in which_counts[-1]) {
        o = o + af_raw[, rd[[wc]]]
      }
    } else {
      o = af_raw
    }
    
    res.spl <- t(o)
    colnames(res.spl) <- afc[['barcodes']]
    rownames(res.spl) <- afg[['gene_ids']]
    
    which_counts = c('U')
    if (usa_mode) {
      rd = list('S' = seq(1, ng), 'U' =  seq(ng + 1, 2*ng), 'A' =  seq(2*ng + 1, 3*ng))
      o = af_raw[, rd[[which_counts[1]]]]
      for (wc in which_counts[-1]) {
        o = o + af_raw[, rd[[wc]]]
      }
    } else {
      o = af_raw
    }
    
    res.uns <- t(o)
    colnames(res.uns) <- afc[['barcodes']]
    rownames(res.uns) <- afg[['gene_ids']]
    
    return(list('RNA' = res.rna,
                'spliced' = res.spl,
                'unspliced' = res.uns))
  }else {
    if (usa_mode) {
      rd = list('S' = seq(1, ng), 'U' =  seq(ng + 1, 2*ng), 'A' =  seq(2*ng + 1, 3*ng))
      o = af_raw[, rd[[which_counts[1]]]]
      for (wc in which_counts[-1]) {
        o = o + af_raw[, rd[[wc]]]
      }
    } else {
      o = af_raw
    }
    res.counts <- t(o)
    colnames(res.counts) <- afc[['barcodes']]
    rownames(res.counts) <- afg[['gene_ids']]
    
    cat('#\t..\n')
    return(res.counts)
  }
}

call_droplets <- function(counts, protocol, cutoff) {
  
  if(cutoff == 'auto') {
    skip.mod <- F
    force.mod <- F
  }
  
  if(cutoff == 'correct') {
    skip.mod <- F
    force.mod <- T
  }
  
  if(cutoff == 'skip') {
    skip.mod <- T
    force.mod <- F
  }
  
  cat('#\t..\tcalling ',protocol,'..\n',
      sep = '')
  
  if (protocol == 'cells') {
    bc.calls <- barcodeRanks(counts)
  }
  if (protocol == 'nuclei') {
    bc.calls <- barcodeRanks(counts[['RNA']])
  }
  
  bc.df <- get_knee_df(counts = counts,
                       inflection = metadata(bc.calls)[['inflection']],
                       protocol = protocol)
  
  bc.annot <- tibble(Knee = metadata(bc.calls)[['knee']],
                     Inflection = metadata(bc.calls)[['inflection']],
                     Above_knee = max(bc.df[['Rank']][bc.df[['nUMI_RNA']] > Knee]),
                     Rank_cutoff = max(bc.df[['Rank']][bc.df[['nUMI_RNA']] > Inflection]),
                     Multimodal = F,
                     Knee.upd = NA,
                     Inflection.upd = NA,
                     Manual = F,
                     Manual.thresh = NA,
                     Knee.manual = NA,
                     Inflection.manual = NA)
  
  bcs.knee <- bc.df[['Barcode']][bc.df[['nUMI_RNA']] > bc.annot[['Knee']]]
  valid.bcs <- bc.df[['Barcode']][bc.df[['nUMI_RNA']] > bc.annot[['Inflection']]]
  cat('#\t..\t\t ',protocol,':\t',length(valid.bcs),'\n',
      sep = '')
  
  cat('#\t..\n',
      '#\t..\ttesting for bimodality..\n',
      sep = '')
  dip <- diptest::dip.test(bc.df[['nUMI_RNA']][bcs.knee])
  cat('#\t..\t\tp-value = ',dip[['p.value']],'\n',
      '#\t..\t\tforce.mod = ',force.mod,'\n',
      sep = '')
  
  if (skip.mod == F) {
    if (dip[['p.value']] < 0.05 || force.mod == T) {
      cat('#\t..\n',
          '#\t..\tfitting bimodal distribution..\n',
          sep = '')
      tmp.fit <- mclust::Mclust(data = bc.df[['nUMI_RNA']][valid.bcs], G = 1:2, verbose = F)
      cat('#\t..\t\tlog-likelihood:\t',tmp.fit$loglik,'\n',
          '#\t..\t\tobservations:\t',tmp.fit$n,'\n',
          '#\t..\t\test.params:\t',tmp.fit$df,'\n',
          '#\t..\t\tBIC:\t\t',tmp.fit$bic,'\n',
          '#\t..\t\tICL:\t\t',tmp.fit$icl,'\n',
          sep = '')
      valid.bcs <- valid.bcs[tmp.fit[['classification']] == 2]
      cat('#\t..\n',
          '#\t..\tcalling ',protocol,' again (only mod2-BCs)..\n',
          sep = '')
      
      if (protocol == 'cells') {
        bc.calls.modal <- barcodeRanks(counts[,valid.bcs])
      }
      if (protocol == 'nuclei') {
        bc.calls.modal <- barcodeRanks(counts[['RNA']][,valid.bcs])
      }
      bc.df <- get_knee_df(counts = counts,
                           inflection = metadata(bc.calls.modal)[['inflection']],
                           protocol = protocol)
      
      bc.annot[['Multimodal']] <- T
      bc.annot[['Knee.upd']] <- metadata(bc.calls.modal)[['knee']]
      bc.annot[['Inflection.upd']] <- metadata(bc.calls.modal)[['inflection']]
      bc.annot[['Rank_cutoff']] <- max(bc.df[['Rank']][bc.df[['nUMI_RNA']] > bc.annot[['Inflection.upd']]])
      
      valid.bcs <- bc.df[['Barcode']][bc.df[['nUMI_RNA']] > bc.annot[['Inflection.upd']]]
      cat('#\t..\t\t',protocol,':\t\t',length(valid.bcs),'\n',
          sep = '')
    }
  }
  
  if (is.numeric(cutoff)) {
    cat('#\t..\n',
        '#\t..\tcalling ',protocol,' again with lower bound = ',cutoff,'..\n',
        sep = '')
    
    if (protocol == 'cells') {
      bc.calls.manual <- barcodeRanks(counts, lower = cutoff)
    }
    if (protocol == 'nuclei') {
      bc.calls.manual <- barcodeRanks(counts[['RNA']], lower = cutoff)
    }
    bc.df <- get_knee_df(counts = counts,
                         inflection = metadata(bc.calls.manual)[['inflection']],
                         protocol = protocol)
    
    bc.annot[['Manual']] <- T
    bc.annot[['Manual.thresh']] <- cutoff
    bc.annot[['Knee.manual']] <- metadata(bc.calls.manual)[['knee']]
    bc.annot[['Inflection.manual']] <- metadata(bc.calls.manual)[['inflection']]
    bc.annot[['Rank_cutoff']] <- max(bc.df[['Rank']][bc.df[['nUMI_RNA']] > bc.annot[['Inflection.manual']]])
    
    valid.bcs <- bc.df[['Barcode']][bc.df[['nUMI_RNA']] > bc.annot[['Inflection.manual']]]
    cat('#\t..\t\t',protocol,':\t\t',length(valid.bcs),'\n',
        sep = '')
  }
  return(list(valid.bcs, bc.df, bc.annot))
}

HTODemux.mcl <- function(object, assay = "HTO", q_l = 1, q_h = 0.005, seed = 42){
  # A function to find the threshold for each hastag, the P, and singlet, doublets and negative.
  # The input is the HTO data matrix (normalized across cells).
  
  assay <- assay %||% DefaultAssay(object = object)
  data <- GetAssayData(object = object, assay = assay, slot = 'data')
  hto.cutoff.metadata = data.frame(cut_off = future.apply::future_apply(data,1,function(x) select_hash_cutoff_mcl(x, q_l = q_l, q_h = q_h), future.seed = T))
  hto.cutoff.metadata$Multi_modal = apply(data,1,function(x) is_multimodal(x))
  hto_mcl.p = t(apply(data,1,function(x) hash_mcl_p(x, seed = seed, q_l = q_l, q_h = q_h)))
  discrete <- data
  discrete[discrete > 0] <- 0
  for (iter in rownames(x = data)) {
    values <- data[iter, ]
    cutoff <- hto.cutoff.metadata[iter,'cut_off']
    discrete[iter, names(x = which(x = values > cutoff))] <- 1
  }
  classification.metadata <- HTO_classifcation(discrete, hto_mcl.p, assay)
  object <- AddMetaData(object = object, metadata = classification.metadata)
  Idents(object) <- paste(assay, "calledFeatures", sep = '_')
  doublets <- rownames(x = object[[]])[which(object[[paste(assay, 'globalClass', sep = "_")]] == "Doublet")]
  Idents(object = object, cells = doublets) <- "Doublet"
  Idents(object) = factor(Idents(object), levels = c('Doublet', 'Negative', rownames(object@assays$HTO)))
  object$sampleID <- Idents(object = object)
  return(object)
}

process_seur <- function(seur, n.pca.dims = 50, assay = 'RNA', verbose = F) {
  cat('#\t..\tprocessing ',assay,'-assay..\n', sep = '')
  suppressWarnings({
    if (assay == 'RNA') {
      cat('#\t..\t\trunning sctransform..\n')
      seur <- SCTransform(seur, assay = 'RNA', method = 'qpoisson', verbose = verbose)
      cat('#\t..\t\trunning PCA..\n')
      seur <- RunPCA(seur, assay = 'SCT', verbose = verbose)
      cat('#\t..\t\testimating global dims..\n')
      dims <- round(
        as.numeric(
          maxLikGlobalDimEst(
            data = seur@reductions[['pca']][, 1:n.pca.dims],
            k = 20)))
      cat('#\t..\t\t\tdims:\t',dims,'\n')
      cat('#\t..\t\trunning tSNE..\n')
      seur <- RunTSNE(seur, assay = 'SCT', dims = seq(dims), verbose = verbose)
      cat('#\t..\t\trunning UMAP..\n')
      seur <- RunUMAP(seur, assay = 'SCT', dims = seq(dims), verbose = verbose)
      cat('#\t..\t\trunning FindNeighbors..\n')
      seur <- FindNeighbors(seur, assay = 'SCT', dims = seq(dims), verbose = verbose)
      cat('#\t..\t\trunning FindClusters..\n')
      seur <- FindClusters(seur, assay = 'SCT', verbose = verbose)
    }
    if (assay == 'HTO') {
      cat('#\t..\t\trunning PCA..\n')
      seur <- ScaleData(seur,
                        assay = 'HTO',
                        features = rownames(seur[['HTO']]@counts),
                        verbose = verbose)
      seur <- RunPCA(seur,
                     assay = 'HTO',
                     features = rownames(seur[['HTO']]@counts),
                     reduction.name = 'hto.pca',
                     reduction.key = 'htoPC_',
                     approx=FALSE,
                     verbose = verbose)
      
      cat('#\t..\t\trunning tSNE..\n')
      seur <- RunTSNE(seur,
                      distance.matrix = as.matrix(dist(t(GetAssayData(object = seur, assay = 'HTO')))),
                      reduction.name = 'hto.tsne',
                      reduction.key = 'htotSNE_',
                      verbose = verbose)
      
      cat('#\t..\t\trunning UMAP..\n')
      seur <- RunUMAP(seur,
                      reduction = "hto.pca",
                      dims = seq(rownames(seur@assays$HTO)),
                      reduction.name = 'hto.umap',
                      reduction.key = 'htoUMAP_',
                      metric='correlation',
                      verbose = verbose)
    }
  })
  return(seur)
}

dub_cutoff <- function(x) {
  
  # calculate number of neighbors at each proportion that are doublets
  data.frame("proportion" = x$RNA_doubletNeighborProportion) %>% 
    group_by(proportion) %>% 
    summarize(n_cells = n()) %>% 
    mutate(pct_cells = n_cells / sum(n_cells)) -> data 
  
  
  # find point at which we gain very few doublets as proportion increases 
  cut <- data$proportion[PCAtools::findElbowPoint(variance = sort(data$n_cells, decreasing = T)) + 1]
  vec <- if_else(x$RNA_doubletNeighborProportion <= cut, F, T)
  
  data[['cut']] <- cut
  
  cat('#\t..\n',
      '#\t..\twriting metadata..\n',
      sep = '')
  cat('#\t..\t\t\"rna_doublet_neighbor_proportion_metadata.csv\"\n')
  write.table(x = data,
              file = file.path(mat.stats.path,'rna_doublet_neighbor_proportion_metadata.csv'),
              sep = ',',
              row.names = F)
  
  return(vec)
  
}

# Helper functions
'%!in%' <- function(x,y)!('%in%'(x,y))

get_knee_df <- function(counts, inflection, protocol) {
  total <- rank <- NULL
  if (protocol == 'cells') {
    mat.rna <- counts
    tibble(Barcode = colnames(mat.rna),
           nUMI_RNA = Matrix::colSums(mat.rna),
           Rank = row_number(desc(nUMI_RNA))) %>%
      mutate(State = ifelse(nUMI_RNA > inflection, 'Called', 'Uncalled')) %>%
      arrange(Rank) -> bc.df
  }
  if (protocol == 'nuclei') {
    mat.rna <- counts[['RNA']]
    mat.spl <- counts[['spliced']]
    mat.uns <- counts[['unspliced']]
    
    tibble(Barcode = colnames(mat.rna),
           nUMI_RNA = Matrix::colSums(mat.rna),
           nUMI_spliced = Matrix::colSums(mat.spl),
           nUMI_unspliced = Matrix::colSums(mat.uns),
           Rank = row_number(desc(nUMI_RNA))) %>%
      mutate(State = ifelse(nUMI_RNA > inflection, 'Called', 'Uncalled')) %>%
      arrange(Rank) -> bc.df
  }
  return(bc.df)
}

is_multimodal <- function(x, p.cutoff = 1e-2) {
  # Test if the expression distribution violates unimodal distribution.
  p = diptest::dip.test(x)$p.value
  return(p < p.cutoff)
}

select_hash_cutoff_mcl <- function(x, q_l = 1, q_h = 0.01, seed = 42) {
  # Model HTO data as a mixture of two Gaussian distributions (for normalized [across cells] data)
  # And select HTO cutoff based on mclust (Model based clustering).
  assertthat::assert_that(class(x) == "numeric")
  assertthat::is.number(seed)
  assertthat::assert_that(length(seed) == 1)
  set.seed(seed)
  km <- mclust::Mclust(data = x, G = 2, verbose = F)
  cl <- km$classification
  cl_center = km$parameters$mean
  high_cl <- which(cl_center == max(cl_center))
  low_cl <- which(cl_center != max(cl_center))
  # q_l and q_h are the quantiles for negative and postive cluster, respectively. 
  cutoff <- max(quantile(x[cl == low_cl], q_l), quantile(x[cl == high_cl], q_h))
  # The higher the cut off, the less false positive (the more false negative).
  return(cutoff)
}

hash_mcl_p <- function(x, seed = 3030, q_l = 1, q_h = 0.001) {
  assertthat::assert_that(class(x) == "numeric")
  assertthat::is.number(seed)
  assertthat::assert_that(length(seed) == 1)
  set.seed(seed)
  km <- mclust::Mclust(data = x, G = 2, verbose = F)
  cl <- km$classification
  cl_center = km$parameters$mean
  high_cl <- which(cl_center == max(cl_center))
  low_cl <- which(cl_center != max(cl_center))
  p.high_cl <- km$z[,high_cl]
  # Correct assignment error from Mclust
  p.high_cl[which(x < max(quantile(cl_center[low_cl], q_l), quantile(cl_center[low_cl], q_h)))] = 0
  names(p.high_cl) = names(x)
  return(p.high_cl)
}

HTO_classifcation = function(discrete, hto_mcl.p, assay){
  # Based on HTODemux (Seurat)
  npositive <- colSums(x = discrete)
  classification.global <- npositive
  classification.global[npositive == 0] <- "Negative"
  classification.global[npositive == 1] <- "Singlet"
  classification.global[npositive > 1] <- "Doublet"
  donor.id = rownames(x = discrete)
  hash.max <- apply(X = hto_mcl.p, MARGIN = 2, FUN = max) # This returns the probability of the most likely HashID (based on the Hashtag distribution among cells)
  hash.maxID <- as.character(donor.id[apply(X = hto_mcl.p, MARGIN = 2, FUN = which.max)])
  hash.second <- apply(X = hto_mcl.p, MARGIN = 2, FUN = function(x) sort(x,decreasing = T)[2])
  hash.secondID <- as.character(donor.id[apply(X = hto_mcl.p, MARGIN = 2, FUN = function(x) order(x,decreasing = T)[2])])
  hash.margin <- hash.max - hash.second
  doublet_id <- sapply(X = 1:length(x = hash.maxID), FUN = function(x) {
    return(paste(sort(x = c(hash.maxID[x], hash.secondID[x])), 
                 collapse = "_"))
  })
  classification <- classification.global
  classification[classification.global == "Negative"] <- "Negative"
  classification[classification.global == "Singlet"] <- hash.maxID[which(x = classification.global == "Singlet")]
  classification[classification.global == "Doublet"] <- doublet_id[which(x = classification.global == "Doublet")]
  classification.metadata <- data.frame(hash.maxID, hash.secondID, hash.margin, classification, classification.global)
  colnames(x = classification.metadata) <- paste(assay, c("primaryID", "secondaryID", "margin", "calledFeatures", "globalClass"), sep = "_")
  return(classification.metadata)
}

# Visualization 
plotting.font <- 'Helvetica'
base.size = 10

my.cols <- c(
  'Doublet' = '#E56786', 'Negative' = '#00ACB9', 'Singlet' = '#7E9D00',
  'Called'    = '#C1E24D', 'Uncalled'  = '#FFB2FF',
  'Highlight' = '#de3163', 'Lowlight'  = '#f4b8c9',
  'Greyneric' = '#878686')

softPallet <- function(n){
  colorspace::sequential_hcl(n = n,
                             h1 = 50,
                             h2 = 300,
                             c1 = 70,
                             c2 = 70,
                             l1 = 85,
                             l2 = 85)
}

exponent_format <- function() {
  function(x) {
    parse(text = gsub("e\\+?", " %*% 10^", scales::scientific_format()(x)))
  }
}

make_summary_table <- function(){
  
  stat.tib.10x <- tibble(
    # Raw - reads from fastqs
    'Cells (Loaded)' = double(),
    'Reads (Raw)' = double(),
    
    # Called cells - (Barcode ranks)
    'Cells (Called)' = double(),  'Cells % Loaded (Called)' = character(),
    'Reads (Called)' = double(),  'Reads % Raw (Called)' = character(), 'Reads/Cell (Called)' = double(),
    'UMIs (Called)' = double(), 'UMIs/Cell (Called)' = double(),
    
    # Called singlets (inter-HTO)
    'Cells (inter-HTO)' = double(),  'Cells % Called (inter-HTO)' = character(),
    'Reads (inter-HTO)' = double(),  'Reads % Called (inter-HTO)' = character(), 'Reads/Cell (inter-HTO)' = double(),
    'UMIs (inter-HTO)' = double(), 'UMIs/Cell (inter-HTO)' = double(),
    
    # Called singlets (intra-HTO)
    'Cells (intra-HTO)' = double(),  'Cells % inter-HTO (intra-HTO)' = character(),
    'Reads (intra-HTO)' = double(),  'Reads % inter-HTO (intra-HTO)' = character(), 'Reads/Cell (intra-HTO)' = double(),
    'UMIs (intra-HTO)' = double(), 'UMIs/Cell (intra-HTO)' = double()
  )
  
  for (rnx.i in row_number(rnx.sheet)) {
    
    rnx.tib <- rnx.sheet[rnx.i,]
    
    q.rnx <- rnx.tib[['reaction_id']]
    q.protocol <- rnx.tib[['seq_type']]
    q.cutoff <- rnx.tib[['cutoff']]
    q.organism <- rnx.tib[['align']]
    
    lib.tib <- filter(rnxsheet.rnx2lib,
                      reaction_id == q.rnx)
    rnx.omnisheet <- filter(omnisheet, reaction_id == q.rnx)
    
    q.index.10x <- rnx.omnisheet %>%
      filter(library_type == '10x') %>% 
      select(index) %>% 
      unique() %>% 
      unlist()
    
    q.bcl <- rnx.omnisheet %>%
      select(bcl_folder) %>% 
      unique() %>% 
      unlist()
    
    q.loaded <- rnx.omnisheet %>% 
      filter(library_type == '10x') %>% 
      group_by(reaction_id) %>% 
      summarise(loaded_cells_rnx = sum(loaded_cells)) %>% 
      filter(reaction_id == q.rnx) %>% 
      select(loaded_cells_rnx) %>% 
      unlist()
    
    q.lib.10x <- lib.tib %>% 
      filter(library_type == '10x') %>% 
      select(library_id) %>% 
      unlist()
    
    demult.stats.path <- file.path(
      project.path,
      'scRNAseq',
      'dry-lab',
      'FASTQ',
      q.bcl,
      bcl.convert.version,
      'Reports')
    
    demultiplex.stats <- read_csv(file.path(demult.stats.path,
                                            'Demultiplex_Stats.csv'))
    
    q.reads <- demultiplex.stats %>%
      group_by(SampleID) %>% 
      summarise(total_reads = sum(`# Reads`)) %>%
      filter(SampleID == q.index.10x) %>% 
      select(total_reads) %>% 
      unlist()
    
    mat.stats.path <- file.path(
      project.path,
      snakemake@config[['out_path']],
      snakemake@config[['com_id']],
      'metadata',
      q.rnx)
    
    # Init
    stats.10x <- list()
    
    # Collect barcode ranks data
    bc.df <- read.csv(file.path(mat.stats.path,'rna_cell-barcode_stats.csv'))
    valid.cbs <- bc.df[['Barcode']][bc.df[['State']] == 'Called']
    
    # Read featureDump.txt
    #mat.files.10x <- file.path(dir.proj,'scRNAseq','03_PipelineOut',com.ID,'10x',pool.10x,'res')
    mat.files.10x <- file.path(project.path,snakemake@config[['out_path']],snakemake@config[['com_id']],salmon.version.alevin.fry.version,q.rnx,q.lib.10x,'res')
    feature.dump.10x <- suppressMessages(read_delim(file.path(mat.files.10x,'featureDump.txt'), delim = '\t'))
    
    stats.10x[['Reads (Raw)']] <- q.reads
    stats.10x[['Cells (Loaded)']] <- q.loaded
    
    # CALLED CELLS - BARCODE RANKS
    # 10x
    stats.10x[['Cells (Called)']] <- length(unique(feature.dump.10x[feature.dump.10x[['CB']] %in% valid.cbs,][['CB']]))
    stats.10x[['Cells % Loaded (Called)']] <- paste(as.character(round(stats.10x[['Cells (Called)']] / stats.10x[['Cells (Loaded)']] * 100, 1)),'%')
    
    stats.10x[['Reads (Called)']] <- sum(feature.dump.10x[feature.dump.10x[['CB']] %in% valid.cbs,][['MappedReads']])
    stats.10x[['Reads % Raw (Called)']] <- paste(as.character(round(stats.10x[['Reads (Called)']] / stats.10x[['Reads (Raw)']] * 100, 1)),'%')
    stats.10x[['Reads/Cell (Called)']] <- round(median(feature.dump.10x[feature.dump.10x[['CB']] %in% valid.cbs,][['MappedReads']]))
    
    stats.10x[['UMIs (Called)']] <- round(sum(feature.dump.10x[feature.dump.10x[['CB']] %in% valid.cbs,][['DeduplicatedReads']]))
    stats.10x[['UMIs/Cell (Called)']] <- round(median(feature.dump.10x[feature.dump.10x[['CB']] %in% valid.cbs,][['DeduplicatedReads']]))
    
    if ('hto' %in% lib.tib[['library_type']]) {
      
      # CALLED SINGLETS (INTER-HTO)
      meta.data <- read.csv(file.path(mat.stats.path,'seurat_unfiltered_metadata.csv'))
      cells.negat.inter <- rownames(meta.data[meta.data[['HTO_globalClass']] == 'Negative',])
      cells.doubl.inter <- rownames(meta.data[meta.data[['HTO_globalClass']] == 'Doublet',])
      cells.singl.inter <- rownames(meta.data[meta.data[['HTO_globalClass']] == 'Singlet',])
      
      # 10x
      stats.10x[['Cells (inter-HTO)']] <- length(unique(feature.dump.10x[feature.dump.10x[['CB']] %in% cells.singl.inter,][['CB']]))
      stats.10x[['Cells % Called (inter-HTO)']] <- paste(as.character(round(stats.10x[['Cells (inter-HTO)']] / stats.10x[['Cells (Called)']] * 100, 1)),'%')
      
      stats.10x[['Reads (inter-HTO)']] <- sum(feature.dump.10x[feature.dump.10x[['CB']] %in% cells.singl.inter,][['MappedReads']])
      stats.10x[['Reads % Called (inter-HTO)']] <- paste(as.character(round(stats.10x[['Reads (inter-HTO)']] / stats.10x[['Reads (Called)']] * 100, 1)),'%')
      stats.10x[['Reads/Cell (inter-HTO)']] <- round(median(feature.dump.10x[feature.dump.10x[['CB']] %in% cells.singl.inter,][['MappedReads']]))
      
      stats.10x[['UMIs (inter-HTO)']] <- round(sum(feature.dump.10x[feature.dump.10x[['CB']] %in% cells.singl.inter,][['DeduplicatedReads']]))
      stats.10x[['UMIs/Cell (inter-HTO)']] <- round(median(feature.dump.10x[feature.dump.10x[['CB']] %in% cells.singl.inter,][['DeduplicatedReads']]))
      
      # CALLED SINGLETS (INTRA-HTO)
      cells.singl.intra <- rownames(meta.data[meta.data[['HTO_globalClass']] == 'Singlet' & meta.data[['RNA_doubletNeighborBool']] == F,])
      
      stats.10x[['Cells (intra-HTO)']] <- length(unique(feature.dump.10x[feature.dump.10x[['CB']] %in% cells.singl.intra,][['CB']]))
      stats.10x[['Cells % inter-HTO (intra-HTO)']] <- paste(as.character(round(stats.10x[['Cells (intra-HTO)']] / stats.10x[['Cells (inter-HTO)']] * 100, 1)),'%')
      
      stats.10x[['Reads (intra-HTO)']] <- sum(feature.dump.10x[feature.dump.10x[['CB']] %in% cells.singl.intra,][['MappedReads']])
      stats.10x[['Reads % inter-HTO (intra-HTO)']] <- paste(as.character(round(stats.10x[['Reads (intra-HTO)']] / stats.10x[['Reads (inter-HTO)']] * 100, 1)),'%')
      stats.10x[['Reads/Cell (intra-HTO)']] <- round(median(feature.dump.10x[feature.dump.10x[['CB']] %in% cells.singl.intra,][['MappedReads']]))
      
      stats.10x[['UMIs (intra-HTO)']] <- round(sum(feature.dump.10x[feature.dump.10x[['CB']] %in% cells.singl.intra,][['DeduplicatedReads']]))
      stats.10x[['UMIs/Cell (intra-HTO)']] <- round(median(feature.dump.10x[feature.dump.10x[['CB']] %in% cells.singl.intra,][['DeduplicatedReads']]))
    }
    # Bind to collective df
    stat.tib.10x <- bind_rows(stat.tib.10x,stats.10x)
  }
  
  stat.df.10x <- as.data.frame(stat.tib.10x)
  rownames(stat.df.10x) <- rnx.sheet[['reaction_id']]
  
  stat.df.10x %>%
    kbl(escape = F, col.names = c('','',
                                  'remaining','% loaded','remaining','% sequenced','per cell','remaining','per cell',
                                  'remaining','% called','remaining','% called','per cell','remaining','per cell',
                                  'remaining','% singlets','remaining','% singlets','per cell','remaining','per cell'),
        format.args = list(big.mark = ','),align = rep('c',25),
        booktabs = T) %>% 
    row_spec(0, color = '#dedede', angle = 25, align = 'center', font_size = '12') %>%
    kable_material_dark(full_width = F, html_font = 'helvetica', lightable_options = 'basic') %>%
    column_spec(1, bold = T, color = '#dedede') %>% 
    column_spec(column = c(4,5,6,7,8,9,10,18,19,20,21,22,23,24), color = '#ffbd69') %>% 
    column_spec(column = c(11,12,13,14,15,16,17), color = '#adffff') %>% 
    add_header_above(c(
      ' ' = 1,
      'Cells' = 1,
      'Reads' = 1,
      'Cells' = 2, 'Reads' = 3, 'UMIs' = 2,
      'Cells' = 2, 'Reads' = 3, 'UMIs' = 2,
      'Cells' = 2, 'Reads' = 3, 'UMIs' = 2),
      color = '#dedede') %>% 
    add_header_above(c(' ' = 1,
                       'Loaded' = 1,
                       'Sequenced' = 1,
                       'Called Cells' = 7,
                       'Between-HTO doublets and negatives removed' = 7,
                       'Within-HTO doublets removed' = 7),
                     align = 'center',
                     bold = T, color = c('#dedede','#dedede','#dedede',
                                                  '#ffbd69','#adffff','#ffbd69','#adffff','#ffbd69','#adffff')) %>% 
                                                    add_header_above(c('Summary stats (RNA expression)' = 3,
                                                                       ' ' = 21),
                                                                     align = 'center',
                                                                     color = '#dedede',
                                                                     font_size = 20,
                                                                     bold = T) %>% 
    save_kable(file = file.path(summary.path,'rna_summary.html'), self_contained = T)
  
  if ('hto' %in% lib.tib[['library_type']]) {
    stat.tib.hto <- tibble(
      # Raw - reads from fastqs
      'Cells (Loaded)' = double(),
      'Reads (Raw)' = double(),
      
      # Called cells - (Barcode ranks)
      'Cells (Called)' = double(),  'Cells % Loaded (Called)' = character(),
      'Reads (Called)' = double(),  'Reads % Raw (Called)' = character(), 'Reads/Cell (Called)' = double(),
      'UMIs (Called)' = double(), 'UMIs/Cell (Called)' = double(),
      
      # Called singlets (inter-HTO)
      'Cells (inter-HTO)' = double(),  'Cells % Called (inter-HTO)' = character(), 
      'Reads (inter-HTO)' = double(),  'Reads % Called (inter-HTO)' = character(), 'Reads/Cell (inter-HTO)' = double(),
      'UMIs (inter-HTO)' = double(), 'UMIs/Cell (inter-HTO)' = double(),
      'Cells (Negative %)' = character(),    
      'Cells (Doublet %)' = character(),    
      'Cells (Singlet %)' = character(),
      'Reads - HTO (Negative %)' = character(),    
      'Reads - HTO (Doublet %)' = character(),    
      'Reads - HTO (Singlet %)' = character(),
      'Reads - RNA (Negative %)' = character(),    
      'Reads - RNA (Doublet %)' = character(),    
      'Reads - RNA (Singlet %)' = character()
    )
    
    for (rnx.i in row_number(rnx.sheet)) {
      
      rnx.tib <- rnx.sheet[rnx.i,]
      
      q.rnx <- rnx.tib[['reaction_id']]
      q.protocol <- rnx.tib[['seq_type']]
      q.cutoff <- rnx.tib[['cutoff']]
      q.organism <- rnx.tib[['align']]
      
      lib.tib <- filter(rnxsheet.rnx2lib,
                        reaction_id == q.rnx)
      rnx.omnisheet <- filter(omnisheet, reaction_id == q.rnx)
      
      q.index.hto <- rnx.omnisheet %>%
        filter(library_type == 'hto') %>% 
        select(index) %>% 
        unique() %>% 
        unlist()
      
      q.bcl <- rnx.omnisheet %>%
        select(bcl_folder) %>% 
        unique() %>% 
        unlist()
      
      q.loaded <- rnx.omnisheet %>% 
        filter(library_type == 'hto') %>% 
        group_by(reaction_id) %>% 
        summarise(loaded_cells_rnx = sum(loaded_cells)) %>% 
        filter(reaction_id == q.rnx) %>% 
        select(loaded_cells_rnx) %>% 
        unlist()
      
      q.lib.hto <- lib.tib %>% 
        filter(library_type == 'hto') %>% 
        select(library_id) %>% 
        unlist()
      
      demult.stats.path <- file.path(
        project.path,
        'scRNAseq',
        'dry-lab',
        'FASTQ',
        q.bcl,
        bcl.convert.version,
        'Reports')
      
      demultiplex.stats <- read_csv(file.path(demult.stats.path,
                                              'Demultiplex_Stats.csv'))
      
      q.reads <- demultiplex.stats %>%
        group_by(SampleID) %>% 
        summarise(total_reads = sum(`# Reads`)) %>%
        filter(SampleID == q.index.hto) %>% 
        select(total_reads) %>% 
        unlist()
      
      mat.stats.path <- file.path(
        project.path,
        snakemake@config[['out_path']],
        snakemake@config[['com_id']],
        'metadata',
        q.rnx)
      
      # Init
      stats.hto <- list()
      
      # Collect barcode ranks data
      bc.df <- read.csv(file.path(mat.stats.path,'rna_cell-barcode_stats.csv'))
      valid.cbs <- bc.df[['Barcode']][bc.df[['State']] == 'Called']
      
      # Read featureDump.txt
      mat.files.hto <- file.path(project.path,snakemake@config[['out_path']],snakemake@config[['com_id']],salmon.version.alevin.fry.version,q.rnx,q.lib.hto,'res')
      featDump.hto <- suppressMessages(read_delim(file.path(mat.files.hto,'featureDump.txt'), delim = '\t'))
      
      stats.hto[['Reads (Raw)']] <- q.reads
      stats.hto[['Cells (Loaded)']] <- q.loaded
      
      # CALLED CELLS - BARCODE RANKS
      # HTO
      stats.hto[['Cells (Called)']] <- length(unique(featDump.hto[featDump.hto[['CB']] %in% valid.cbs,][['CB']]))
      stats.hto[['Cells % Loaded (Called)']] <- paste(as.character(round(stats.hto[['Cells (Called)']] / stats.hto[['Cells (Loaded)']] * 100, 1)),'%')
      
      stats.hto[['Reads (Called)']] <- sum(featDump.hto[featDump.hto[['CB']] %in% valid.cbs,][['MappedReads']])
      stats.hto[['Reads % Raw (Called)']] <- paste(as.character(round(stats.hto[['Reads (Called)']] / stats.hto[['Reads (Raw)']] * 100, 1)),'%')
      stats.hto[['Reads/Cell (Called)']] <- round(median(featDump.hto[featDump.hto[['CB']] %in% valid.cbs,][['MappedReads']]))
      
      stats.hto[['UMIs (Called)']] <- round(sum(featDump.hto[featDump.hto[['CB']] %in% valid.cbs,][['DeduplicatedReads']]))
      stats.hto[['UMIs/Cell (Called)']] <- round(median(featDump.hto[featDump.hto[['CB']] %in% valid.cbs,][['DeduplicatedReads']]))
      
      # CALLED SINGLETS (INTER-HTO)
      meta.data <- read.csv(file.path(mat.stats.path,'seurat_unfiltered_metadata.csv'))
      cells.negat.inter <- rownames(meta.data[meta.data[['HTO_globalClass']] == 'Negative',])
      cells.doubl.inter <- rownames(meta.data[meta.data[['HTO_globalClass']] == 'Doublet',])
      cells.singl.inter <- rownames(meta.data[meta.data[['HTO_globalClass']] == 'Singlet',])
      
      # HTO
      stats.hto[['Cells (Negative %)']] <- paste(as.character(round(length(cells.negat.inter) / length(valid.cbs) * 100, 1)),'%')
      stats.hto[['Cells (Doublet %)']] <- paste(as.character(round(length(cells.doubl.inter) / length(valid.cbs) * 100, 1)),'%')
      stats.hto[['Cells (Singlet %)']] <- paste(as.character(round(length(cells.singl.inter) / length(valid.cbs) * 100, 1)),'%')
      
      stats.hto[['Reads - HTO (Negative %)']] <- paste(as.character(round(sum(featDump.hto[featDump.hto[['CB']] %in% cells.negat.inter,][['MappedReads']]) / stats.hto[['Reads (Called)']] * 100, 1)),'%')
      stats.hto[['Reads - HTO (Doublet %)']] <- paste(as.character(round(sum(featDump.hto[featDump.hto[['CB']] %in% cells.doubl.inter,][['MappedReads']]) / stats.hto[['Reads (Called)']] * 100, 1)),'%')
      stats.hto[['Reads - HTO (Singlet %)']] <- paste(as.character(round(sum(featDump.hto[featDump.hto[['CB']] %in% cells.singl.inter,][['MappedReads']]) / stats.hto[['Reads (Called)']] * 100, 1)),'%')
      
      stats.hto[['Reads - RNA (Negative %)']] <- paste(as.character(round(sum(feature.dump.10x[feature.dump.10x[['CB']] %in% cells.negat.inter,][['MappedReads']]) / stats.10x[['Reads (Called)']] * 100, 1)),'%')
      stats.hto[['Reads - RNA (Doublet %)']] <- paste(as.character(round(sum(feature.dump.10x[feature.dump.10x[['CB']] %in% cells.doubl.inter,][['MappedReads']]) / stats.10x[['Reads (Called)']] * 100, 1)),'%')
      stats.hto[['Reads - RNA (Singlet %)']] <- paste(as.character(round(sum(feature.dump.10x[feature.dump.10x[['CB']] %in% cells.singl.inter,][['MappedReads']]) / stats.10x[['Reads (Called)']] * 100, 1)),'%')
      
      stats.hto[['Cells (inter-HTO)']] <- length(unique(featDump.hto[featDump.hto[['CB']] %in% cells.singl.inter,][['CB']]))
      stats.hto[['Cells % Called (inter-HTO)']] <- paste(as.character(round(stats.hto[['Cells (inter-HTO)']] / stats.hto[['Cells (Called)']] * 100, 1)),'%')
      
      stats.hto[['Reads (inter-HTO)']] <- sum(featDump.hto[featDump.hto[['CB']] %in% cells.singl.inter,][['MappedReads']])
      stats.hto[['Reads % Called (inter-HTO)']] <- paste(as.character(round(stats.hto[['Reads (inter-HTO)']] / stats.hto[['Reads (Called)']] * 100, 1)),'%')
      stats.hto[['Reads/Cell (inter-HTO)']] <- round(median(featDump.hto[featDump.hto[['CB']] %in% cells.singl.inter,][['MappedReads']]))
      
      stats.hto[['UMIs (inter-HTO)']] <- round(sum(featDump.hto[featDump.hto[['CB']] %in% cells.singl.inter,][['DeduplicatedReads']]))
      stats.hto[['UMIs/Cell (inter-HTO)']] <- round(median(featDump.hto[featDump.hto[['CB']] %in% cells.singl.inter,][['DeduplicatedReads']]))
      
      # Bind to collective df
      stat.tib.hto <- bind_rows(stat.tib.hto,stats.hto)
    }
    
    stat.df.hto <- as.data.frame(stat.tib.hto)
    rownames(stat.df.hto) <- rnx.sheet[['reaction_id']]
    
    stat.df.hto %>%
      kbl(escape = F, col.names = c('','',
                                    'remaining','% loaded','remaining','% sequenced','per cell','remaining','per cell',
                                    'remaining','% called','remaining','% called','per cell','remaining','per cell',
                                    'Negative %','Doublet %','Singlet %',
                                    'Negative %','Doublet %','Singlet %',
                                    'Negative %','Doublet %','Singlet %'),
          format.args = list(big.mark = ','),
          booktabs = T) %>% 
      row_spec(0, color = '#dedede', angle = 25,align = 'center', font_size = '12') %>%
      kable_material_dark(full_width = F, html_font = 'helvetica', lightable_options = 'basic') %>%
      column_spec(1, bold = T, color = '#dedede') %>%
      column_spec(column = c(4,5,6,7,8,9,10,18,19,20,21,22,23,24,25,26), color = '#ffbd69') %>% 
      column_spec(column = c(11,12,13,14,15,16,17), color = '#adffff') %>%
      add_header_above(c(' ' = 1,
                         'Cells' = 1,
                         'Reads' = 1,
                         'Cells' = 2, 'Reads' = 3, 'UMIs' = 2,
                         'Cells' = 2, 'Reads' = 3, 'UMIs' = 2,
                         'Cells' = 3, 'Reads (HTO)' = 3, 'Reads (RNA)' = 3),
                       color = '#dedede') %>% 
      add_header_above(c(' ' = 1,
                         'Loaded' = 1,
                         'Sequenced' = 1,
                         'Called Cells' = 7,
                         'Between-HTO doublets and negatives removed' = 7,
                         'Classification proportions' = 9),
                       align = 'center',
                       bold = T, color = c('#dedede','#dedede','#dedede',
                                                    '#ffbd69','#adffff','#ffbd69')) %>% 
                                                      add_header_above(c('Summary stats (HTO expression)' = 4,
                                                                         ' ' = 22),
                                                                       align = 'center',
                                                                       color = '#dedede',
                                                                       font_size = 20,
                                                                       bold = T) %>% 
      save_kable(file = file.path(summary.path,'hto_summary.html'), self_contained = T)
  }
  
}

make_plot_top10bcl <- function(df) {
  
  df <- df[1:10,]
  
  p <- ggplot(df) +
    geom_bar(mapping = aes(x = reorder(dual_index, reads), y = reads),
             stat = 'identity',
             fill = softPallet(1)) +
    labs(x = 'Index pair', y = 'Reads') +
    scale_y_continuous(labels = exponent_format()) +
    coord_flip() +
    theme_minimal(base_size = base.size) +
    theme(text = element_text(family = plotting.font))
  
  p <- patchwork::wrap_plots(p)
  return(p)
}

make_plot_bclDistribution <- function(df) {
  
  p <- ggplot(df, aes(x = index, y = reads, fill = class)) +
    geom_bar(stat = 'identity') +
    coord_flip() +
    labs(x = 'Index', y = 'Reads') +
    scale_y_continuous(labels = exponent_format()) +
    scale_fill_manual(name = '',
                      values = softPallet(3)) +
    theme_minimal(base_size = base.size) +
    theme(text = element_text(family = plotting.font), legend.position = 'bottom')
  
  p <- patchwork::wrap_plots(p)
  return(p)
}

make_plot_barcodeRanks <- function(df, annot, ref, rnx.name, rnx.type) {
  
  knee <- annot[['Knee']]
  inflection <- annot[['Inflection']]
  
  breaks <- pretty(range(df[['nUMI_RNA']]), n = nclass.scott(df[['nUMI_RNA']]), min.n = 1)
  bwidth <- breaks[2] - breaks[1]
  
  # rank plot
  p1 <- ggplot(df) +
    geom_vline(xintercept = 80000, linetype = 5, color = my.cols[['Greyneric']], size = .5) +
    geom_line(data = ref, mapping = aes(x = Rank, y = nUMI_RNA), alpha = 0.75, linewidth = 1, linetype = 3, color = my.cols[['Greyneric']]) +
    geom_line(mapping = aes(x = Rank, y = nUMI_RNA, col = State), linewidth = 2)
  
  # frequency plot
  p2 <- ggplot(df) +
    geom_histogram(aes(x = nUMI_RNA, fill = State), bins = length(breaks))
  
  if (annot[['Multimodal']] == F) {
    
    p1 <- p1 +
      geom_hline(aes(yintercept = knee, linetype = 'Knee'), color = my.cols[['Highlight']], linewidth = .75) +
      scale_linetype_manual(name = '',
                            limits = c('Reference','Knee','Optimal plateau drop'),
                            values = c(3, 2, 5),
                            guide = guide_legend(nrow = 2, override.aes = list(color = c(my.cols[['Greyneric']], my.cols[['Highlight']], my.cols[['Greyneric']]))))
    
    p2 <- p2 +
      geom_vline(aes(xintercept = knee, linetype = 'Knee'), size = .75, color = my.cols[['Highlight']]) +
      scale_linetype_manual(name = '',
                            limits = c('Knee'),
                            values = c(2),
                            guide = guide_legend(override.aes = list(color = c(my.cols[['Highlight']])))) +
      ggtitle(rnx.name, subtitle = paste(annot$Rank_cutoff, rnx.type))
    
  }
  if (annot[['Multimodal']] == T) {
    
    knee.new <- annot[['Knee.upd']]
    inflection.new <- annot[['Inflection.upd']]
    
    p1 <- p1 +
      geom_hline(aes(yintercept = knee.new, linetype = 'Knee (new)'), color = my.cols[['Highlight']], size = .75) +
      geom_hline(aes(yintercept = knee, linetype = 'Knee (old)'), color = my.cols[['Lowlight']], size = .75) +
      scale_linetype_manual(name = '',
                            limits = c('Reference', 'Knee (new)', 'Knee (old)', 'Optimal plateau drop'),
                            values = c(3, 2, 2, 5),
                            guide = guide_legend(nrow = 2, override.aes = list(color = c(my.cols[['Greyneric']],my.cols[['Highlight']], my.cols[['Lowlight']], my.cols[['Greyneric']]))))
    
    p2 <- p2 +
      geom_vline(aes(xintercept = knee.new, linetype = 'Knee (new)'), size = .75, color = my.cols[['Highlight']]) +
      geom_vline(aes(xintercept = knee, linetype = 'Knee (old)'), size = .75, color = my.cols[['Lowlight']]) +
      scale_linetype_manual(name = '',
                            limits = c('Knee (new)','Knee (old)'),
                            values = c(2, 2),
                            guide = guide_legend(override.aes = list(color = c(my.cols[['Highlight']], my.cols[['Lowlight']])))) +
      ggtitle(rnx.name, subtitle = paste(annot$Rank_cutoff, rnx.type))
    
  }
  if (annot[['Manual']] == T) {
    
    knee.new <- annot[['Knee.manual']]
    inflection.new <- annot[['Inflection.manual']]
    
    p1 <- p1 +
      geom_hline(aes(yintercept = knee.new, linetype = 'Knee (new)'), color = my.cols[['Highlight']], size = 1) +
      geom_hline(aes(yintercept = knee, linetype = 'Knee (old)'), color = my.cols[['Lowlight']], size = 1) +
      scale_linetype_manual(name = '',
                            limits = c('Reference', 'Knee (new)', 'Knee (old)', 'Optimal plateau drop'),
                            values = c(3, 2, 2, 5),
                            guide = guide_legend(nrow = 1, override.aes = list(color = c(my.cols[['Greyneric']],my.cols[['Highlight']], my.cols[['Lowlight']], my.cols[['Greyneric']]))))
    
    p2 <- p2 +
      geom_vline(aes(xintercept = knee.new, linetype = 'Knee (new)'), size = 1, color = my.cols[['Highlight']]) +
      geom_vline(aes(xintercept = knee, linetype = 'Knee (old)'), size = 1, color = my.cols[['Lowlight']]) +
      scale_linetype_manual(name = '',
                            limits = c('Knee (new)','Knee (old)'),
                            values = c(2, 2),
                            guide = guide_legend(override.aes = list(color = c(my.cols[['Highlight']], my.cols[['Lowlight']])))) +
      ggtitle(rnx.name, subtitle = paste(annot$Rank_cutoff, rnx.type, '\t|\tlower threshold:',annot[['Manual.thresh']]))
    
  }
  
  p1 <- p1 +
    scale_color_manual(name = '',
                       labels = c('Called', 'Uncalled'),
                       limits = c('Called', 'Uncalled'),
                       values = c(my.cols[['Called']], my.cols[['Uncalled']]),
                       guide = guide_legend(nrow = 1)) +
    scale_x_log10(labels = exponent_format(), limits = c(1,max(df$Rank))) +
    scale_y_log10(limits = c(1,max(df$nUMI_RNA))) +
    labs(x = 'Rank', y = '') +
    annotation_logticks(sides = 'b') +
    theme_minimal(base_size = base.size) +
    theme(text = element_text(family = plotting.font), legend.position = 'bottom', axis.line.y = element_line(size = .5, color = my.cols[['Greyneric']], linetype = 'solid'), axis.text.y = element_blank(), plot.margin = margin(l = 0))
  
  p2 <- p2 +
    scale_fill_manual(name = '',
                      values = c(my.cols[['Called']],my.cols[['Uncalled']]),
                      guide = "none") +
    labs(x = 'UMI counts', y = 'Frequency') +
    scale_x_log10(labels = exponent_format(), limits = c(1,max(df$nUMI_RNA))) +
    scale_y_log10(labels = exponent_format(), limits = c(1,max(df$Rank))) +
    coord_flip() +
    annotation_logticks() +
    theme_minimal(base_size = base.size) +
    theme(text = element_text(family = plotting.font), legend.position = 'none', plot.margin = margin(r = 0))
  
  p.empty <- ggplot() +
    theme_minimal()
  combined <- patchwork::wrap_plots(p2, p1) + theme(legend.justification = c(3.5,0))
  
  return(combined)
}

make_plot_htoThresh <- function(hto.sample.calling.metadata, hto.cutoff.metadata, rnx.name, q_l = 1, q_h = 0.001) {
  
  hto.cutoff.metadata[['hto']] <- hto.cutoff.metadata[['hto_name']]
  
  p <- ggplot(hto.sample.calling.metadata, aes(x = expression, fill = above_cutoff)) +
    geom_histogram(bins = 100) +
    facet_wrap(~sample_id, scales = 'free',ncol = 3) +
    xlab('Expression (CLR-norm)') +
    ylab('UMI counts (HTO)') +
    scale_y_log10(labels = exponent_format()) +
    scale_fill_manual(name = '',
                      labels = c('HTO negative','HTO positive'),
                      values = c(my.cols[['Uncalled']],my.cols[['Called']]),
                      guide = guide_legend(override.aes = list(fill = c(my.cols[['Uncalled']],my.cols[['Called']])))) +
    ggtitle(rnx.name, subtitle = paste0('Q_negative = ', q_l, '; Q_positive = ', q_h)) +
    theme_minimal(base_size = base.size) +
    theme(text = element_text(family = plotting.font), legend.position = 'bottom')
  
  p <- patchwork::wrap_plots(p)
  return(p)
}

make_plot_htoVlnGlobal <- function(meta.data, rnx.name) {
  
  p <- ggplot(meta.data,
              aes(x = HTO_globalClass,
                  y = nCount_HTO,
                  fill = HTO_globalClass)) +
    geom_violin() +
    geom_jitter(alpha = .25) +
    scale_fill_manual(values = c(my.cols[['Doublet']],
                                 my.cols[['Negative']],
                                 my.cols[['Singlet']])) +
    scale_y_log10(labels = exponent_format()) +
    annotation_logticks(sides = 'l') +
    labs(x = 'HTO classification (Global)',
         y = 'UMI counts (HTO)') +
    theme_minimal(base_size = base.size) +
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = 22.5)) +
    ggtitle(rnx.name)

  p <- patchwork::wrap_plots(p)  
  return(p)
}

make_plot_htoVlnIndividual <- function(meta.data, rnx.name) {
  
  meta.data[['sampleID']] <- factor(meta.data[['sampleID']],
                                    levels = c(c('Doublet','Negative'),setdiff(meta.data[['sampleID']],c('Doublet','Negative'))))
  
  tmp.color.palette <- c(my.cols[['Doublet']],
                         my.cols[['Negative']],
                         softPallet(length(levels(meta.data[['sampleID']])) - 2))
  
  p <- ggplot(meta.data,
              aes(x = sampleID,
                  y = nCount_HTO,
                  fill = sampleID)) +
    geom_violin() +
    geom_jitter(alpha = .25) +
    scale_fill_manual(values = tmp.color.palette,
                      limits = levels(meta.data[['sampleID']])) +
    scale_y_log10(labels = exponent_format()) +
    annotation_logticks(sides = 'l') +
    labs(x = 'HTO classification (Individual)',
         y = 'UMI counts (HTO)') +
    theme_minimal(base_size = base.size) +
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = 22.5)) +
    ggtitle(rnx.name)
  
  p <- patchwork::wrap_plots(p)
  return(p)
}

make_plot_htotsne <- function(meta.data, choosen.class = 'Individual', rnx.name) {
  
  if (choosen.class == 'Individual') {
    meta.data[['sampleID']] <- factor(meta.data[['sampleID']],
                                      levels = c(c('Doublet','Negative'),
                                                 setdiff(meta.data[['sampleID']],
                                                         c('Doublet','Negative'))))
    tmp.color.palette <- c(my.cols[['Doublet']],
                           my.cols[['Negative']],
                           softPallet(length(levels(meta.data[['sampleID']])) - 2))
    
    p <- ggplot(meta.data) +
      geom_point(aes(x = HTO_tsne1, y = HTO_tsne2, col = sampleID, alpha = .25)) +
      scale_color_manual(values = tmp.color.palette) +
      labs(x = 'tSNE 1', y = 'tSNE 2', color = 'HTO classification') +
      theme_minimal(base_size = base.size) +
      guides(alpha = 'none') +
      theme(text = element_text(family = plotting.font), legend.position = 'bottom') +
      ggtitle(rnx.name)
  }
  
  if (choosen.class == 'Global') {
    
    tmp.color.palette <- c(my.cols[['Doublet']],
                           my.cols[['Negative']],
                           my.cols[['Singlet']])
    
    p <- ggplot(meta.data) +
      geom_point(aes(x = HTO_tsne1, y = HTO_tsne2, col = HTO_globalClass, alpha = .25)) +
      scale_color_manual(values = tmp.color.palette) +
      labs(x = 'tSNE 1', y = 'tSNE 2', color = 'HTO classification') +
      theme_minimal(base_size = base.size) +
      guides(alpha = 'none') +
      theme(text = element_text(family = plotting.font), legend.position = 'bottom') +
      #NoGrid() +
      ggtitle(rnx.name)
  }
  
  p <- patchwork::wrap_plots(p)
  return(p)
}

make_plot_InfDoubl_std <- function(meta.data, reduction = 'tSNE', rnx.name) {
  
  opacity.vec <- ifelse(meta.data[['RNA_recoveredDoubletBool']] == T, 1, .25)
  
  p <- ggplot(meta.data)
  
  if (reduction == 'tSNE') {
    p <- p + geom_point(aes(x = SCT_tsne1, y = SCT_tsne2, col = RNA_recoveredDoubletBool), alpha = opacity.vec) +
      labs(x = 'tSNE 1', y = 'tSNE 2', color = 'Intra-HTO doublet')
  }
  if (reduction == 'UMAP') {
    p <- p + geom_point(aes(x = SCT_umap1, y = SCT_umap2, col = RNA_recoveredDoubletBool), alpha = opacity.vec) +
      labs(x = 'UMAP 1', y = 'UMAP 2', color = 'Intra-HTO doublet')
  }
  
  p <- p +
    scale_color_manual(values = c(my.cols[['Negative']],
                                  my.cols[['Doublet']])) +
    theme_minimal(base_size = base.size) +
    guides(alpha = 'none') +
    theme(text = element_text(family = plotting.font), legend.position = 'bottom') +
    ggtitle(rnx.name)
  
  p <- patchwork::wrap_plots(p)
  return(p)
}

make_plot_InfDoubl_prop <- function(meta.data, reduction = 'tSNE', rnx.name) {
  
  p <- ggplot(meta.data)
  
  if (reduction == 'tSNE') {
    p <- p + geom_point(aes(x = SCT_tsne1, y = SCT_tsne2, col = RNA_doubletNeighborProportion)) +
      labs(x = 'tSNE 1', y = 'tSNE 2', color = 'Doublet neighbours proportion')
  }
  if (reduction == 'UMAP') {
    p <- p + geom_point(aes(x = SCT_umap1, y = SCT_umap2, col = RNA_doubletNeighborProportion)) +
      labs(x = 'UMAP 1', y = 'UMAP 2', color = 'Doublet neighbour proportion')
  }
  
  p <- p +
    scale_color_viridis() +
    theme_minimal(base_size = base.size) +
    guides(alpha = 'none') +
    theme(text = element_text(family = plotting.font), legend.position = 'bottom') +
    ggtitle(rnx.name)
  
  p <- patchwork::wrap_plots(p)
  return(p)
}

make_plot_InfDoubl_knee <- function(doub.data, rnx.name) {
  
  doub.data <- mutate(doub.data,
                      above_cut = proportion > cut)
  
  p <- ggplot(doub.data) +
    geom_point(aes(x = proportion, y = n_cells, color = above_cut)) +
    labs(x = 'Proportion - doublet neighbors', y = 'Cells', color = 'Inferential doublet') +
    theme_minimal(base_size = base.size) +
    scale_color_manual(name = '',
                       values = c(my.cols[['Negative']],my.cols[['Doublet']]),
                       guide = guide_legend(override.aes = list(color = c(my.cols[['Negative']],my.cols[['Doublet']])))) +
    theme(text = element_text(family = plotting.font), legend.position = 'bottom') +
    ggtitle(rnx.name)
  
  p <- patchwork::wrap_plots(p)
  return(p)
  
}

make_plot_InfDoubl_cut <- function(meta.data, reduction = 'tSNE', rnx.name) {
  
  opacity.vec <- ifelse(meta.data$RNA_doubletNeighborBool == T, 1, .25)
  
  p <- ggplot(meta.data)
  
  if (reduction == 'tSNE') {
    p <- p +
      geom_point(aes(x = SCT_tsne1, y = SCT_tsne2, col = RNA_doubletNeighborBool), alpha = opacity.vec) +
      labs(x = 'tSNE 1', y = 'tSNE 2', color = 'Inferential doublet')
  }
  if (reduction == 'UMAP') {
    p <- p +
      geom_point(aes(x = SCT_umap1, y = SCT_umap2, col = RNA_doubletNeighborBool), alpha = opacity.vec) +
      labs(x = 'UMAP 1', y = 'UMAP 2', color = 'Inferential doublet')
  }
  
  p <- p +
    scale_color_manual(values = c(my.cols[['Negative']],
                                  my.cols[['Doublet']])) +
    theme_minimal(base_size = base.size) +
    guides(alpha = 'none') +
    theme(text = element_text(family = plotting.font), legend.position = 'bottom') +
    ggtitle(rnx.name)
  
  p <- patchwork::wrap_plots(p)
  return(p)
}

make_plot_htoVlnAll <- function(meta.data, rnx.name) {
  
  mutate(meta.data,
         doublet_type = factor(ifelse(HTO_doubletBool == T, 'Inter-HTO doublet',
                                      ifelse(RNA_recoveredDoubletBool == T, 'Intra-HTO doublet (profile)',
                                             ifelse(RNA_doubletNeighborBool == T, 'Intra-HTO doublet (neighbours)',
                                                    ifelse(HTO_globalClass == 'Negative', 'Negative',
                                                           'Singlet')))),
                               levels = c('Negative','Singlet','Inter-HTO doublet','Intra-HTO doublet (profile)','Intra-HTO doublet (neighbours)'))) -> meta.data
  
  p <- ggplot(meta.data[meta.data$HTO_globalClass != 'Negative',],
              aes(x = doublet_type,
                  y = nCount_RNA,
                  fill = doublet_type)) +
    geom_violin() +
    geom_jitter(alpha = .25) +
    scale_fill_manual(values = c(my.cols[['Singlet']],my.cols[['Doublet']],my.cols[['Doublet']],my.cols[['Doublet']]),
                      limits = c('Singlet','Inter-HTO doublet','Intra-HTO doublet (profile)','Intra-HTO doublet (neighbours)')) +
    scale_y_log10(labels = exponent_format()) +
    annotation_logticks(sides = 'l') +
    labs(x = '', y = 'UMI counts') +
    theme_minimal(base_size = base.size) +
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = 22.5)) +
    ggtitle(rnx.name)
  
  p <- patchwork::wrap_plots(p)
  return(p)
}

make_plot_rnaVln <- function(meta.data, rnx.name) {
  
  tmp.color.palette <- softPallet(length(unique(meta.data[['sampleID']])))
  
  p <- ggplot(meta.data,
              aes(x = sampleID,
                  y = nCount_RNA,
                  fill = sampleID)) +
    geom_violin() +
    geom_jitter(alpha = .25) +
    scale_fill_manual(values = tmp.color.palette,
                      limits = unique(meta.data[['sampleID']])) +
    scale_y_log10(labels = exponent_format()) +
    annotation_logticks(sides = 'l') +
    labs(x = 'Sample',
         y = 'UMI counts') +
    theme_minimal(base_size = base.size) +
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = 22.5)) +
    ggtitle(rnx.name)
  
  p <- patchwork::wrap_plots(p)
  return(p)
}

make_plot_umap_class <- function(meta.data, rnx.name) {
  
  tmp.color.palette <- softPallet(length(unique(meta.data[['sampleID']])))
  
  p <- ggplot(meta.data) +
    geom_point(aes(x = SCT_umap1, y = SCT_umap2, col = sampleID, alpha = .25)) +
    scale_color_manual(values = tmp.color.palette,guide = guide_legend(nrow = 3)) +
    labs(x = 'UMAP 1', y = 'UMAP 2', color = 'HTO classification') +
    theme_minimal(base_size = base.size) +
    guides(alpha = 'none') +
    theme(text = element_text(family = plotting.font), legend.position = 'bottom') +
    ggtitle(rnx.name, subtitle = 'SCT assay')
  
  p <- patchwork::wrap_plots(p)
  return(p)
}

# Data
bc.df.ref <- read.csv('data/bc-df-ref.csv')