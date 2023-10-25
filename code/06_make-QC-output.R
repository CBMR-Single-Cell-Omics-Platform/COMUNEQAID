################################################################################
##########                            Init                            ##########
################################################################################
# Packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(viridis)
  library(kableExtra)
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

aggr.path <- file.path(
  com.path,
  'aggregated')

plot.path <- file.path(
  com.path,
  'plots')

summary.path <- file.path(
  com.path,
  'summary')

# Files
log.file <- file.path(
  plot.path,
  paste0(date.and.time,'.log'))


# Metadata
suppressMessages({
  omnisheet <- 
    as_tibble(snakemake@config[['sample_sheet']]) %>% 
    full_join(as_tibble(snakemake@config[['sample2reaction']])) %>% 
    full_join(as_tibble(snakemake@config[['reaction_sheet']])) %>% 
    full_join(as_tibble(snakemake@config[['reaction2library']])) %>% 
    full_join(as_tibble(snakemake@config[['library_sheet']])) %>% 
    full_join(as_tibble(snakemake@config[['library2sequencing']])) %>% 
    full_join(as_tibble(snakemake@config[['sequencing_sheet']]))

  rnxsheet.rnx2lib <-
    full_join(as_tibble(snakemake@config[['reaction_sheet']]),
              as_tibble(snakemake@config[['reaction2library']]))
  
  seq.sheet <-
    as_tibble(snakemake@config[['sequencing_sheet']])
  
  rnx.sheet <-
    as_tibble(snakemake@config[['reaction_sheet']])
  
  samplesheet.sample2rnx <- full_join(as_tibble(snakemake@config[['sample_sheet']]),
            as_tibble(snakemake@config[['sample2reaction']]))
})

dir.create(plot.path, recursive = T, showWarnings = F)
dir.create(summary.path, recursive = T, showWarnings = F)

sink.file <- file(log.file, 'w')
sink(file = sink.file,
     append = T,
     type = 'output',
     split = T)

# Startup message
cat(rep('#',80),'\n',
    '#####                                                                      #####\n',
    '#####                              Make QC output                          #####\n',
    '#####                                                                      #####\n',
    '#####                                                     6/6              #####\n',
    '#####           ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▒▒▒▒▒▒▒▒            #####\n',
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

# Table generation
suppressMessages({
  make_summary_table()
})

cat('#\tmaking FASTQ plots..\n')
# FASTQ generation
#for (bcl in snakemake@config[['sequencing_sheet']][['bcl_folder']]) {
for (seq.id in seq.sheet[["sequencing_id"]]) {
  
  bcl.folder <- filter(seq.sheet, sequencing_id == seq.id)[["bcl_folder"]]
  
  cat('#\t..\t',bcl.folder,':\n',
      sep = '')
  
  fastq.stats.path <- file.path(
    snakemake@config[['project_path']],
    snakemake@config[['scop_id']],
    snakemake@config[['fastq_path']],
    bcl.folder,
    'metadata')
  
  df.unknown.barcodes <- read.csv(file.path(fastq.stats.path,'unknown-barcodes.csv'))
  p.unknown.barcodes <- make_plot_top10bcl(df.unknown.barcodes)
  
  cat('#\t..\t\t- \"00_top-undetermined-barcodes.png\"\n')
  ggsave(filename = paste0('00_top-undetermined-barcodes.png'),
         plot = p.unknown.barcodes,
         path = file.path(plot.path,bcl.folder),
         width = 6,
         height = 4,
         dpi = 'retina')
  
  df.read.demultiplexing <- read_csv(file.path(fastq.stats.path,'read-demultiplexing.csv'))
  p.read.demultiplexing <- plot_demux_stats(df.read.demultiplexing)
  
  cat('#\t..\t\t- \"00_read-index-distribution.png\"\n')
  ggsave(filename = paste0('00_read-index-distribution.png'),
         plot = p.read.demultiplexing,
         path = file.path(plot.path,bcl.folder),
         width = 6,
         height = 4,
         dpi = 'retina')
}




df <- data.frame(Type = character(),
                 CountSum = double(),
                 CountPer = double(),
                 CountTot = double(),
                 Reaction = character())

df.hto <- data.frame(Library = character(),
                     Class = character(),
                     CountPer = double(),
                     Labelpos = double(),
                     LibraryID = character())

cat('#\tmaking reaction plots..\n')
for (rnx in snakemake@config[['reaction_sheet']][['reaction_id']]) {
  
  dir.create(file.path(plot.path,rnx), recursive = T, showWarnings = F)
  
  cat('#\t..\t',rnx,':\n',
      sep = '')
  
  lib.tib <- filter(rnxsheet.rnx2lib,
                    reaction_id == rnx)
  
  mat.stats.path <- file.path(
    project.path,
    snakemake@config[['out_path']],
    snakemake@config[['com_id']],
    'metadata',
    rnx)
  
  # Cell Calling (RNA)
  ## Barcode rank plot
  df.cell.barcode.stats <- read.csv(file.path(mat.stats.path,'rna_cell-barcode_stats.csv'))
  df.rna.cell.barcode.calling.metadata <- read.csv(file.path(mat.stats.path,'rna_cell-barcode_calling_metadata.csv'))
  
  p.celular.barcode.calling <- make_plot_barcodeRanks(df.cell.barcode.stats,
                                            df.rna.cell.barcode.calling.metadata,
                                            bc.df.ref,
                                            rnx,
                                            snakemake@config[['reaction_sheet']][['seq_type']][snakemake@config[['reaction_sheet']][['reaction_id']] == rnx])
  
  cat('#\t..\t\t- \"01_barcode-calling.png\"\n')
  ggsave(filename = '01_barcode-calling.png',
         plot = p.celular.barcode.calling,
         path = file.path(plot.path,rnx),
         width = 8,
         height = 5,
         dpi = 'retina')
  
  # Cell calling (HTO)
  ## Inter-HTO calling
  ### Individual HTO distributions
  hto.sample.calling.metadata <- read.csv(file.path(mat.stats.path,'hto_sample_calling_metadata.csv'))
  hto.cutoff.metadata <- read.csv(file.path(mat.stats.path,'hto_cutoff_metadata.csv'))
  
  p.hto.barcode.calling <- make_plot_htoThresh(hto.sample.calling.metadata,hto.cutoff.metadata,rnx)
  
  cat('#\t..\t\t- \"02_antibody-calling.png\"\n')
  ggsave(filename = '02_antibody-calling.png',
         plot = p.hto.barcode.calling,
         path = file.path(plot.path,rnx),
         width = 6,
         height = 5,
         dpi = 'retina')
  
  ### HTO expression
  #### Global classification
  seurat.unfiltered.metadata <- read.csv(file.path(mat.stats.path,'seurat_unfiltered_metadata.csv'), colClasses=c("sampleID"="character"))
  
  p.violin.global <- make_plot_htoVlnGlobal(seurat.unfiltered.metadata, rnx)
  
  cat('#\t..\t\t- \"02_HTO-counts_violin_global.png\"\n')
  ggsave(filename = '02_HTO-counts_violin_global.png',
         plot = p.violin.global,
         path = file.path(plot.path,rnx),
         width = 6,
         height = 6,
         dpi = 'retina')
  
  #### Sample classification {.tabset}
  p.violin.sample <- make_plot_htoVlnIndividual(seurat.unfiltered.metadata, rnx)
  
  cat('#\t..\t\t- \"02_HTO-counts_violin_sample.png\"\n')
  ggsave(filename = '02_HTO-counts_violin_sample.png',
         plot = p.violin.sample,
         path = file.path(plot.path,rnx),
         width = 6,
         height = 6,
         dpi = 'retina')

  
  ### Clustering
  #### Global classification
  p.hto.tsne.global <- make_plot_htotsne(seurat.unfiltered.metadata, 'Global', rnx)
  
  cat('#\t..\t\t- \"02_HTO-counts_tSNE_global.png\"\n')
  ggsave(filename = '02_HTO-counts_tSNE_global.png',
         plot = p.hto.tsne.global,
         path = file.path(plot.path,rnx),
         width = 6,
         height = 6,
         dpi = 'retina')
  
  #### Sample classification
  p.hto.tsne.sample <- make_plot_htotsne(seurat.unfiltered.metadata, 'Individual', rnx)
  
  cat('#\t..\t\t- \"02_HTO-counts_tSNE_sample.png\"\n')
  ggsave(filename = '02_HTO-counts_tSNE_sample.png',
         plot = p.hto.tsne.sample,
         path = file.path(plot.path,rnx),
         width = 6,
         height = 6,
         dpi = 'retina')
  
  ## Intra-HTO doublet calling
  ### Intra-HTO doublets from doublet RNA profile
  p.recovered.doublets.tsne <- make_plot_InfDoubl_std(seurat.unfiltered.metadata, reduction = 'tSNE', rnx)
  p.recovered.doublets.umap <- make_plot_InfDoubl_std(seurat.unfiltered.metadata, reduction = 'UMAP', rnx)
  
  cat('#\t..\t\t- \"03_RNA-counts_tSNE_recoverd-doublets.png\"\n')
  ggsave(filename = '03_RNA-counts_tSNE_recoverd-doublets.png',
         plot = p.recovered.doublets.tsne,
         path = file.path(plot.path,rnx),
         width = 6,
         height = 6,
         dpi = 'retina')
  
  cat('#\t..\t\t- \"03_RNA-counts_UMAP_recoverd-doublets.png\"\n')
  ggsave(filename = '03_RNA-counts_UMAP_recoverd-doublets.png',
         plot = p.recovered.doublets.umap,
         path = file.path(plot.path,rnx),
         width = 6,
         height = 6,
         dpi = 'retina')
  
  ### Proportion of neighboring doublets
  #p.doublet.neighbors.tsne <- make_plot_InfDoubl_prop(seurat.unfiltered.metadata, reduction = 'tSNE', rnx)
  p.doublet.neighbors.umap <- make_plot_InfDoubl_prop(seurat.unfiltered.metadata, reduction = 'UMAP', rnx)
  
  #cat('#\t..\t\t- \"04_RNA-counts_tSNE_doublet_neighbors.png\"\n')
  #ggsave(filename = '04_RNA-counts_tSNE_doublet_neighbors.png',
  #       plot = p.doublet.neighbors.tsne,
  #       path = file.path(plot.path,rnx),
  #       width = 6,
  #       height = 6,
  #       dpi = 'retina')
  
  cat('#\t..\t\t- \"03_RNA-counts_UMAP_doublet_neighbors.png\"\n')
  ggsave(filename = '03_RNA-counts_UMAP_doublet_neighbors.png',
         plot = p.doublet.neighbors.umap,
         path = file.path(plot.path,rnx),
         width = 6,
         height = 6,
         dpi = 'retina')

  ### Doublets from neighbor proportion threshold
  doub.data <- read.csv(file.path(mat.stats.path,'rna_doublet_neighbor_proportion_metadata.csv'))
  
  p.doublet.neighbour.proportion <- make_plot_InfDoubl_knee(doub.data, rnx)
  
  cat('#\t..\t\t- \"03_RNA_neighbor-doublet-calling.png\"\n')
  ggsave(filename = '03_RNA_neighbor-doublet-calling.png',
         plot = p.doublet.neighbour.proportion,
         path = file.path(plot.path,rnx),
         width = 6,
         height = 4,
         dpi = 'retina')
  
  ### Intra-HTO doublets profile + proportion
  #p.neigbour.doublets.tsne <- make_plot_InfDoubl_cut(seurat.unfiltered.metadata, reduction = 'tSNE', rnx)
  p.neigbour.doublets.umap <- make_plot_InfDoubl_cut(seurat.unfiltered.metadata, reduction = 'UMAP', rnx)
  
  #cat('#\t..\t\t- \"04_RNA-counts_tSNE_neighbor-doublets.png\"\n')
  #ggsave(filename = '04_RNA-counts_tSNE_neighbor-doublets.png',
  #       plot = p.neigbour.doublets.tsne,
  #       path = file.path(plot.path,rnx),
  #       width = 6,
  #       height = 4,
  #       dpi = 'retina')
  
  cat('#\t..\t\t- \"03_RNA-counts_UMAP_neighbor-doublets.png\"\n')
  ggsave(filename = '03_RNA-counts_UMAP_neighbor-doublets.png',
         plot = p.neigbour.doublets.umap,
         path = file.path(plot.path,rnx),
         width = 6,
         height = 6,
         dpi = 'retina')
  
  ### RNA Violin plot - HTO
  p.violin.doublet.distribution <- make_plot_htoVlnAll(seurat.unfiltered.metadata, rnx)
  
  cat('#\t..\t\t- \"03_RNA-counts_violin_doublet_distribution.png\"\n')
  ggsave(filename = '03_RNA-counts_violin_doublet_distribution.png',
         plot = p.violin.doublet.distribution,
         path = file.path(plot.path,rnx),
         width = 6,
         height = 6,
         dpi = 'retina')
  
  # Summarize all
  ## Violin
  seurat.filtered.metadata <- read.csv(file.path(mat.stats.path,'seurat_filtered_metadata.csv'), colClasses=c("sampleID"="character"))
  
  
  p.rna.violin.sample <- make_plot_rnaVln(seurat.filtered.metadata, rnx)
  
  cat('#\t..\t\t- \"04_RNA-counts_violin_sample.png\"\n')
  ggsave(filename = '04_RNA-counts_violin_sample.png',
         plot = p.rna.violin.sample,
         path = file.path(plot.path,rnx),
         width = 6,
         height = 6,
         dpi = 'retina')
  
  ## UMAP {.tabset}
  p.sct.umap.sample <- make_plot_umap_class(seurat.filtered.metadata, rnx)
  
  cat('#\t..\t\t- \"04_SCT_UMAP_sample.png\"\n')
  ggsave(filename = '04_SCT_UMAP_sample.png',
         plot = p.sct.umap.sample,
         path = file.path(plot.path,rnx),
         width = 6,
         height = 6,
         dpi = 'retina')
  
  
  
  cbs.called <- df.cell.barcode.stats[df.cell.barcode.stats[['State']] == 'Called',][['Barcode']]
  cbs.uncalled <- df.cell.barcode.stats[df.cell.barcode.stats[['State']] == 'Uncalled',][['Barcode']]
  
  
  if ('nuclei' %in% lib.tib[['seq_type']]) {
    tmp.counts.spl.called <- sum(df.cell.barcode.stats[df.cell.barcode.stats[['State']] == 'Called',][['nUMI_spliced']])
    tmp.counts.spl.uncalled <- sum(df.cell.barcode.stats[df.cell.barcode.stats[['State']] == 'Uncalled',][['nUMI_spliced']])
    tmp.counts.uns.called <- sum(df.cell.barcode.stats[df.cell.barcode.stats[['State']] == 'Called',][['nUMI_unspliced']])
    tmp.counts.uns.uncalled <- sum(df.cell.barcode.stats[df.cell.barcode.stats[['State']] == 'Uncalled',][['nUMI_unspliced']])
    tmp.counts.tot.called <- tmp.counts.spl.called + tmp.counts.uns.called
    tmp.counts.tot.uncalled <- tmp.counts.spl.uncalled + tmp.counts.uns.uncalled
    
    tmp.percen.spl.called <- round(tmp.counts.spl.called/tmp.counts.tot.called*100)
    tmp.percen.spl.uncalled <- round(tmp.counts.spl.uncalled/tmp.counts.tot.uncalled*100)
    tmp.percen.uns.called <- round(tmp.counts.uns.called/tmp.counts.tot.called*100)
    tmp.percen.uns.uncalled <- round(tmp.counts.uns.uncalled/tmp.counts.tot.uncalled*100)
    
    tmp.df <- data.frame(ReadType = c('Spliced','Unspliced','Spliced','Unspliced'),
                         CellType = c('Called','Called','Uncalled','Uncalled'),
                         CountSum = c(tmp.counts.spl.called,tmp.counts.uns.called,tmp.counts.spl.uncalled,tmp.counts.uns.uncalled),
                         CountPer = c(tmp.percen.spl.called,tmp.percen.uns.called,tmp.percen.spl.uncalled,tmp.percen.uns.uncalled),
                         CountTot = c(tmp.counts.tot.called,tmp.counts.tot.called,tmp.counts.tot.uncalled,tmp.counts.tot.uncalled),
                         Reaction = c(rnx,rnx,rnx,rnx))
    df <- rbind(df,tmp.df)
  }
  
  if ('hto' %in% lib.tib[['library_type']]) {
    
    q.lib.10x <- lib.tib %>% 
      filter(library_type == '10x') %>% 
      select(library_id) %>% 
      unlist()
    
    q.lib.hto <- lib.tib %>% 
      filter(library_type == 'hto') %>% 
      select(library_id) %>% 
      unlist()
    
    hto.cell.barcode.stats <- read.csv(file.path(mat.stats.path,'hto_cell-barcode_stats.csv'))
    #meta.data <- read.csv(file.path(dir.rnx.outs,pool.10x,'seurat-full-metadata.csv'))
    
    bcs.negative <- rownames(seurat.unfiltered.metadata[seurat.unfiltered.metadata[['HTO_globalClass']] == 'Negative',])
    bcs.singlet <- rownames(seurat.unfiltered.metadata[seurat.unfiltered.metadata[['HTO_globalClass']] == 'Singlet',])
    bcs.doublet <- rownames(seurat.unfiltered.metadata[seurat.unfiltered.metadata[['HTO_globalClass']] == 'Doublet',])
    
    tmp.counts.rna.called <- sum(df.cell.barcode.stats[df.cell.barcode.stats[['State']] == 'Called',][['nUMI_RNA']])
    tmp.counts.rna.uncalled <- sum(df.cell.barcode.stats[df.cell.barcode.stats[['State']] == 'Uncalled',][['nUMI_RNA']])
    tmp.counts.rna.negative <- sum(df.cell.barcode.stats[df.cell.barcode.stats[['Barcode']] %in% bcs.negative,][['nUMI_RNA']])
    tmp.counts.rna.singlet <- sum(df.cell.barcode.stats[df.cell.barcode.stats[['Barcode']] %in% bcs.singlet,][['nUMI_RNA']])
    tmp.counts.rna.doublet <- sum(df.cell.barcode.stats[df.cell.barcode.stats[['Barcode']] %in% bcs.doublet,][['nUMI_RNA']])
    tmp.counts.rna.total <- tmp.counts.rna.called + tmp.counts.rna.uncalled
    
    #bc.data
    tmp.counts.hto.called <- sum(hto.cell.barcode.stats[hto.cell.barcode.stats[['Barcode']] %in% cbs.called,][['nUMI_HTO']])
    tmp.counts.hto.uncalled <- sum(hto.cell.barcode.stats[hto.cell.barcode.stats[['Barcode']] %in% cbs.uncalled,][['nUMI_HTO']])
    tmp.counts.hto.negative <- sum(hto.cell.barcode.stats[hto.cell.barcode.stats[['Barcode']] %in% bcs.negative,][['nUMI_HTO']])
    tmp.counts.hto.singlet <- sum(hto.cell.barcode.stats[hto.cell.barcode.stats[['Barcode']] %in% bcs.singlet,][['nUMI_HTO']])
    tmp.counts.hto.doublet <- sum(hto.cell.barcode.stats[hto.cell.barcode.stats[['Barcode']] %in% bcs.doublet,][['nUMI_HTO']])
    tmp.counts.hto.total <- tmp.counts.hto.called + tmp.counts.hto.uncalled
    
    tmp.percen.rna.uncalled <- round(tmp.counts.rna.uncalled/tmp.counts.rna.total*100)
    tmp.percen.rna.negative <- round(tmp.counts.rna.negative/tmp.counts.rna.total*100)
    tmp.percen.rna.singlet <- round(tmp.counts.rna.singlet/tmp.counts.rna.total*100)
    tmp.percen.rna.doublet <- round(tmp.counts.rna.doublet/tmp.counts.rna.total*100)
    
    tmp.percen.hto.uncalled <- round(tmp.counts.hto.uncalled/tmp.counts.hto.total*100)
    tmp.percen.hto.negative <- round(tmp.counts.hto.negative/tmp.counts.hto.total*100)
    tmp.percen.hto.singlet <- round(tmp.counts.hto.singlet/tmp.counts.hto.total*100)
    tmp.percen.hto.doublet <- round(tmp.counts.hto.doublet/tmp.counts.hto.total*100)
    
    tmp.labpos.rna.uncalled <- tmp.counts.rna.uncalled/2
    tmp.labpos.rna.negative <- tmp.counts.rna.uncalled + tmp.counts.rna.negative/2
    tmp.labpos.rna.singlet <- tmp.counts.rna.uncalled + tmp.counts.rna.negative + tmp.counts.rna.singlet/2
    tmp.labpos.rna.doublet <- tmp.counts.rna.uncalled + tmp.counts.rna.negative + tmp.counts.rna.singlet + tmp.counts.rna.doublet/2
    
    tmp.labpos.hto.uncalled <- tmp.counts.hto.uncalled/2
    tmp.labpos.hto.negative <- tmp.counts.hto.uncalled + tmp.counts.hto.negative/2
    tmp.labpos.hto.singlet <- tmp.counts.hto.uncalled + tmp.counts.hto.negative + tmp.counts.hto.singlet/2
    tmp.labpos.hto.doublet <- tmp.counts.hto.uncalled + tmp.counts.hto.negative + tmp.counts.hto.singlet + tmp.counts.hto.doublet/2
    
    tmp.df.hto <- data.frame(Library = c('RNA','RNA','RNA','RNA',
                                         'HTO','HTO','HTO','HTO'),
                             Class = c('Uncalled','Negative','Singlet','Doublet',
                                       'Uncalled','Negative','Singlet','Doublet'),
                             CountSum = c(tmp.counts.rna.uncalled,tmp.counts.rna.negative,tmp.counts.rna.singlet,tmp.counts.rna.doublet,
                                          tmp.counts.hto.uncalled,tmp.counts.hto.negative,tmp.counts.hto.singlet,tmp.counts.hto.doublet),
                             CountPer = c(tmp.percen.rna.uncalled,tmp.percen.rna.negative,tmp.percen.rna.singlet,tmp.percen.rna.doublet,
                                          tmp.percen.hto.uncalled,tmp.percen.hto.negative,tmp.percen.hto.singlet,tmp.percen.hto.doublet),
                             Labelpos = c(tmp.labpos.rna.uncalled,tmp.labpos.rna.negative,tmp.labpos.rna.singlet,tmp.labpos.rna.doublet,
                                          tmp.labpos.hto.uncalled,tmp.labpos.hto.negative,tmp.labpos.hto.singlet,tmp.labpos.hto.doublet),
                             LibraryID = c(q.lib.10x ,q.lib.10x,q.lib.10x,q.lib.10x,
                                         q.lib.hto,q.lib.hto,q.lib.hto,q.lib.hto))
    
    df.hto <- rbind(df.hto,tmp.df.hto)
    
    mat.files.10x <- file.path(project.path,snakemake@config[['out_path']],snakemake@config[['com_id']],salmon.version.alevin.fry.version,rnx,q.lib.10x,'res')
    mat.files.hto <- file.path(project.path,snakemake@config[['out_path']],snakemake@config[['com_id']],salmon.version.alevin.fry.version,rnx,q.lib.hto,'res')
    
    feature.dump.10x <- suppressMessages(read_delim(file.path(mat.files.10x,'featureDump.txt'), delim = '\t'))
    featDump.hto <- suppressMessages(read_delim(file.path(mat.files.hto,'featureDump.txt'), delim = '\t'))
    
    cells.negat.inter <- rownames(seurat.unfiltered.metadata[seurat.unfiltered.metadata[['HTO_globalClass']] == 'Negative',])
    cells.doubl.inter <- rownames(seurat.unfiltered.metadata[seurat.unfiltered.metadata[['HTO_globalClass']] == 'Doublet',])
    cells.singl.inter <- rownames(seurat.unfiltered.metadata[seurat.unfiltered.metadata[['HTO_globalClass']] == 'Singlet',])
    
    common.CBs <- intersect(feature.dump.10x$CB,featDump.hto$CB)
    
    feature.dump.10x %>% 
      filter(CB %in% common.CBs) %>% 
      arrange(CB) -> feature.dump.10x
    
    featDump.hto %>% 
      filter(CB %in% common.CBs) %>% 
      arrange(CB) -> featDump.hto
    
    common.CBs.sort <- sort(common.CBs)
    
    tmpdf <- ifelse(common.CBs.sort %in% cells.negat.inter, 'Negative',
                    ifelse(common.CBs.sort %in% cells.doubl.inter, 'Doublet',
                           ifelse(common.CBs.sort %in% cells.singl.inter, 'Singlet', 'Uncalled')))
    
    p <- ggplot(mapping = aes(x = feature.dump.10x$DeduplicatedReads,
                              y = featDump.hto$DeduplicatedReads,
                              col = tmpdf)) +
      geom_point() + scale_x_log10() + scale_y_log10() +
      scale_color_manual(values = c(my.cols[['Doublet']],my.cols[['Negative']],my.cols[['Singlet']],my.cols[['Uncalled']])) +
      annotation_logticks() +
      labs(x = 'UMI counts (RNA)',
           y = 'UMI counts (HTO)',
           col = '') +
      theme_minimal(base_size = base.size) +
      theme(legend.position = 'bottom') +
      ggtitle(rnx)
    
    p <- patchwork::wrap_plots(p)
    cat('#\t..\t\t- \"04_UMI-UMI_scatter_global.png\"\n')
    ggsave(filename = '04_UMI-UMI_scatter_global.png',
           plot = p,
           path = file.path(plot.path,rnx),
           width = 8,
           height = 6,
           dpi = 'retina')
  }
}

if (length(rnx.sheet[['reaction_id']]) > 1) {
  
  rnx <- 'aggregated'
  
  mat.stats.path <- file.path(
    project.path,
    snakemake@config[['out_path']],
    snakemake@config[['com_id']],
    'metadata',
    rnx)
  
  # Summarize all
  ## Violin
  
  seurat.filtered.metadata <- read.csv(file.path(mat.stats.path,'seurat_filtered_metadata.csv'))
  
  
  p.rna.violin.sample <- make_plot_rnaVln(seurat.filtered.metadata, rnx)
  
  cat('#\t..\t\t- \"04_RNA-counts_violin_sample.png\"\n')
  ggsave(filename = '04_RNA-counts_violin_sample.png',
         plot = p.rna.violin.sample,
         path = file.path(plot.path,rnx),
         width = 6,
         height = 6,
         dpi = 'retina')
  
  ## UMAP {.tabset}
  p.sct.umap.sample <- make_plot_umap_class(seurat.filtered.metadata, rnx)
  
  cat('#\t..\t\t- \"04_SCT_UMAP_sample.png\"\n')
  ggsave(filename = '04_SCT_UMAP_sample.png',
         plot = p.sct.umap.sample,
         path = file.path(plot.path,rnx),
         width = 6,
         height = 6,
         dpi = 'retina')
}
if (length(rnx.sheet[['reaction_id']]) == 1) {
  rnx <- rnx.sheet[['reaction_id']]
}

if ('nuclei' %in% lib.tib[['seq_type']]) {
  
  df[['Labelpos']] <- ifelse(df[['ReadType']] == 'Unspliced',
                             df[['CountSum']]/2, df[['CountTot']] - df[['CountSum']]/2)
  
  p.distribution.splice <- ggplot(data = df, aes(x = Reaction, y = CountSum, fill = ReadType)) +
    geom_bar(stat = 'identity') + 
    geom_text(aes(label = paste0(CountPer,'%'),y = Labelpos),size = 3) +
    coord_flip() + ggtitle('Spliced/Unspliced') +
    xlab('Reaction') + ylab('UMI counts') + scale_fill_discrete(name = 'Read type') +
    facet_wrap(~CellType, nrow = 2) +
    scale_fill_manual(values = softPallet(2)) +
    theme_minimal(base_size = base.size) +
    theme(text = element_text(family = plotting.font))
  
  p.distribution.splice <- patchwork::wrap_plots(p.distribution.splice)
  cat('#\t..\t\t- \"05_UMI-distribution_splice-class.png\"\n')
  ggsave(filename = '05_UMI-distribution_splice-class.png',
         plot = p.distribution.splice,
         path = file.path(plot.path,rnx),
         width = 12,
         height = 6,
         dpi = 'retina')
}

if ('hto' %in% lib.tib[['library_type']]) {
  p.distribution.hto <- ggplot(data = df.hto,
                  aes(x = LibraryID,
                      y = CountSum,
                      fill = factor(Class, levels = c('Doublet','Singlet','Negative','Uncalled')))) +
    geom_bar(stat = 'identity') + 
    geom_text(aes(label = paste0(CountPer, '%'),
                  y = Labelpos),
              size = 3) +
    facet_wrap(~Library) +
    coord_flip() + ggtitle('UMI distribution by HTO class') +
    xlab('Library ID') + ylab('UMI counts') + scale_fill_discrete(name = 'Classification') +
    scale_fill_manual(values = c(my.cols[['Doublet']],my.cols[['Singlet']],my.cols[['Negative']],my.cols[['Uncalled']])) +
    labs(fill = 'Classification') + 
    theme_minimal(base_size = base.size) +
    theme(text = element_text(family = plotting.font))
  
  p.distribution.hto <- patchwork::wrap_plots(p.distribution.hto)
  cat('#\t..\t\t- \"05_UMI-distribution_HTO-class.png\"\n')
  ggsave(filename = '05_UMI-distribution_HTO-class.png',
         plot = p.distribution.hto,
         path = file.path(plot.path,rnx),
         width = 12,
         height = 6,
         dpi = 'retina')
  
}
cat('#\t..\n',
    '#\t..done\n',
    '#\n',
    '#\n',
    '###\n',
    rep('#',80),'\n',
    sep = '')

sink()
file.create(snakemake@output[[1]])
