# COMUNEQAID 

**COmbine MUltiple single Nuclei sEQuencing runs AID**  
_single cell/nucleus RNA sequencing workflow_  
_Version 2.0_

**Authors:** Oliver Knights Møller & Lars Roed Ingerslev

**Date:** September 2023

## Updates

### Version 2.0
- bcl2fastq has been replaced with BCL-convert v4.0.3.
- The salmon & alevin-fry has been upgraded to v1.10.2 and v0.8.2, respectivly.
- Updated output structure. (filtered, unfiltered, velocity).
- Updated QC plots
- New folder structure (see section: "Folder structure")

### Version 1.2
- The salmon tool has been upgraded from v1.5.0 to v1.9.0 and alevin-fry has been upgraded from v0.3.0 to v0.7.0. These upgrades fixed a weird behavior where select genes would not be detected in the same abundance as with other mappers.

### Version 1.1
- Back-end Workflow has been changed substantially. It is now easier to execute more jobs at one time by comma-separating the comIDs when executing. It is also easier and very quick to combine multiple seqs of same samples.
- HTO demultiplexing is done with a new and optimised strategy. Instead of executing the Seurat function; HTODemux() across multiple quantiles and picking option that has maximum singlet-yield, instead we now model each distinct HTO as a mixture of two normal distributions (a negative cluster from background noise and a positive cluster from true HTO signal).


## Introduction
The COMUNEQAID pipeline takes single-cell RNA sequencing (Illumina NovaSeq 6000) output and provides:
- ready-to-analyse Seurat objects
- comprehensive QC-output
- intermediary files, metadata & logs

## Output
### Seurat objects
The following set of Seurat objects will be available for each reaction, as well as for an aggregated object (merged assays from all reactions), if aplicaple.
- 'seurat_unfiltered.rds'
- 'seurat_filtered.rds'
- 'seurat_velocity_unfiltered.rds'
- 'seurat_velocity_filtered.rds'

All Seurat objects contains an RNA assay (unspliced, spliced & ambiguous UMI counts). If HTO protocol has been applied, The HTO assay will be available in all objects (except for the aggregated objects).
Velocity objects additionally contains a spliced assay (spliced & ambigous UMI counts) and an unspliced assay (unspliced UMI counts). These can be applied in RNA velocity analysis, be used for diagnistics, or to just play around with.

Filtered objects aims to only contain singlets. If HTO library was used to keep track of sample identity, the filtered object will only contain singlets, as classified accordimng to the HTO library, and subsequently according to calling of inferential doublets.
Unfiltered object contains all droplets that has been classified as a cell-containing droplet by the barcode ranks-calling.

If HTO workflow hasnt been used, filtered and uflitered contains the same CBs (all droplets that has been classified as a cell-containing droplet by the barcode ranks-calling.)

#### Seurat Objects Metadata Columns

##### **seurat_filtered.rds**

These are the per cell metadata columns available for the **seurat_filtered.rds** object:

- `reactionID`: [10x reaction ID]
- `sampleID`: [Sample ID, classified by HTO library *if HTO workflow]
- `nCount_RNA`: [RNA assay UMI count]
- `nFeature_RNA`: [RNA assay feature count]
- `nCount_HTO`: [HTO assay UMI count*if HTO workflow]
- `nFeature_HTO`: [HTO assay feature count *if HTO workflow]
- `nCount_SCT`: [SCT assay UMI count]
- `nFeature_SCT`: [SCT assay feature count]
- `seurat_clusters`: [UMAP cluster identity, determined by Seurat::FindClusters]

##### **seurat_unfiltered.rds**

If HTO workflow, the **seurat_unfiltered.rds** object contains all the metadata columns present in **seurat_filtered.rds** and additionally includes:
- `HTO_primaryID`: [Primary HTO]
- `HTO_secondaryID`: [Secondary HTO]
- `HTO_margin`: [Delta expression (primary - secondary HTO)]
- `HTO_calledFeatures`: [Called HTO features ("_"-separated)]
- `HTO_globalClass`: [Negative/Singlet/Doublet]
- `HTO_doubletBool`: [HTO library doublet *bool*]
- `RNA_recoveredDoubletBool`: [Recovered doublet *bool*]
- `RNA_doubletNeighborBool`: [Neighbour doublet *bool*]
- `RNA_doubletNeighborProportion`: [Proportion of neighbours that are recovered doublets]

##### **seurat_velocity_unfiltered.rds** & **seurat_velocity_filtered.rds**
Like **seurat_unfiltered.rds** & **seurat_filtered.rds**, but with the addition of the following columns:
- `nCount_spliced`: [Spliced assay UMI count]
- `nFeature_spliced`: [Spliced assay feature count]
- `nCount_unspliced`: [Unspliced assay UMI count]
- `nFeature_unspliced`: [Unspliced assay feature count]

### Metadata

### Plots
#### FASTQ
...pending...
#### Reaction
...pending...
#### Aggregated
...pending...
### Summary
#### ***rna_summary.html*** summary
Contains, per reaction, cell, read and UMI (RNA assay) stats across the following stages:
- `Loaded` [The number of cells loaded]
- `Sequenced` [The number of sequenced RNA-reads]
- `Called Cells` [Cell/RNA-read-stats from remaining cells after cell calling]
- `Between-HTO doublets and negatives removed` [Cell/RNA-read-stats from remaining cells after removing doublets and negatives based on HTO-library]
- `Within-HTO doublets removed` [Cell/RNA-read-stats from remaining cells after  removing inferential doublets based on profile from known doublets]

#### ***hto_summary.html*** summary
Contains, per reaction, cell, read and UMI (HTO assay) stats across the following stages:
- `Loaded` [The number of cells loaded]
- `Sequenced` [The number of sequenced HTO-reads]
- `Called Cells` [Cell/HTO-read-stats from remaining cells after cell calling]
- `Between-HTO doublets and negatives removed` [Cell/HTO-read-stats from remaining cells after removing doublets and negatives based on HTO-library]
- `Classification proportions` [Negative/Doublet/Singlet proportions after previous step]

## Folder structure
- {scop id}
  - output
    - {comuneqaid id}
      - COMUNEQAID_v2.0
        - {reaction id}
          - `{yyyy-mm-dd_hh-mm-ss}.log`
          - `seurat_filtered.rds`
          - `seurat_unfiltered.rds`
          - `seurat_velocity_filtered.rds`
          - `seurat_velocity_unfiltered.rds`
        - aggregated  *(if multiple reaction_id's)*
          - `{yyyy-mm-dd_hh-mm-ss}.log`
          - `seurat_filtered.rds`
          - `seurat_unfiltered.rds`
          - `seurat_velocity_filtered.rds`
          - `seurat_velocity_unfiltered.rds`
        - plots
          - {bcl folder}
            - `00_read-index-distribution.png`
            - `00_top-undetermined-barcodes.png`
          - {reaction id}
            - `{yyyy-mm-dd_hh-mm-ss}.log`
            - `01_barcode-calling.png`
            - `02_antibody-calling.png`
            - `02_HTO-counts_tSNE_global.png`
            - `02_HTO-counts_tSNE_sample.png`
            - `02_HTO-counts_violin_global.png`
            - `02_HTO-counts_violin_sample.png`
            - `03_RNA-counts_tSNE_recoverd-doublets.png`
            - `03_RNA-counts_UMAP_doublet_neighbors.png`
            - `03_RNA-counts_UMAP_neighbor-doublets.png`
            - `03_RNA-counts_UMAP_recoverd-doublets.png`
            - `03_RNA-counts_violin_doublet_distribution.png`
            - `03_RNA_neighbor-doublet-calling.png`
            - `04_RNA-counts_violin_sample.png`
            - `04_SCT_UMAP_sample.png`
            - `04_UMI-UMI_scatter_global.png`
          - aggregated *(if multiple reaction_id's)*
            - `04_RNA-counts_violin_sample.png`
            - `04_SCT_UMAP_sample.png`
            - `05_UMI-distribution_HTO-class.png`
            - `05_UMI-distribution_splice-class.png`
        - summary
          - hto_summary.html
          - rna_summary.html
        - COMUNEQAID-pipeline
          - `{comuneqaid id}.yml`
          - code
          - `COMUNEQAID_python.yml`
          - `COMUNEQAID_R.yml`
          - data
          - `modules.txt`
          - `README.md`
          - `Snakefile`
  - scRNAseq
    - dry-lab
      - FASTQ
        - `{bcl folder}`
          - `bcl-convert_{version}`
            - `{yyyy-mm-dd_hh-mm-ss}.log`
            - `{FASTQ files}`
            - Logs
            - Reports
            - `sample-sheet_{bcl folder}.csv`
          - metadata
            - `{yyyy-mm-dd_hh-mm-ss}.log`
            - `read-demultiplexing.csv`
            - `unknown-barcodes.csv`
      - pipeline-out
        - `{comuneqaid id}`
          - metadata
            - `{reaction id}`
              - `hto_cell-barcode_stats.csv`
              - `hto_cutoff_metadata.csv`
              - `hto_sample_calling_metadata.csv`
              - `rna_cell-barcode_calling_metadata.csv`
              - `rna_cell-barcode_stats.csv`
              - `rna_doublet_neighbor_proportion_metadata.csv`
              - `seurat_filtered_metadata.csv`
              - `seurat_unfiltered_metadata.csv`
            - aggregated *(if multiple reaction_id's)*
              - `seurat_filtered_metadata.csv`
              - `seurat_unfiltered_metadata.csv`
          - `salmon_{version}_alevin-fry_{version}`
            - `{reaction id}`
              - `{library id (10x)}`
                - `alevin-fry.log`
                - map
                - quant
                - res
                  - alevin
                    - `quants_mat_cols.txt`
                    - `quants_mat.mtx`
                    - `quants_mat_rows.txt`
                  - `featureDump.txt`
                  - `quant.json`
              - `{library id (HTO)}`
                - `alevin-fry.log`
                - `features.tsv`
                - `hto_index`
                - map
                - quant
                - res
                  - alevin
                    - `quants_mat_cols.txt`
                    - `quants_mat.mtx`
                    - `quants_mat_rows.txt`
                  - `featureDump.txt`
                  - `quant.json`
                - `t2g_hto.tsv`


## Workflow
### BCL to Fastq conversion
Binary base call (BCL) files are converted to fastq files using BCL-convert (v4.0.3).

### Building references
#### RNA
Reference is a build with Cellranger[ref] and modified into discriminating between mature mRNA and pre-mRNA using the ***make_splici_txome*** function from the alevin-fry[ref], and has been built for *homo sapiens* (GRCh38), *mus musculus* (GRCm38), *rattus-norvegicus* (Rnor_6.0) and *macaca-mulatta* (Mmul 10) using the salmon index tool (version 1.10.2)[ref]

#### HTO
Antibody reference is created with salmon index (v1.10.2)[ref]. The features flag is invoked to tell salmon to index based on a tap-separated file (containing sample name and HTO-sequence).

### Mapping & quantification
#### RNA
salmon alevin (v1.10.2)[ref] is run in sketch mode, which does pseudoalignment with structural constraints. Library type & protocol is defined as ISR & 10x chromium v3, respectively.

#### HTO
salmon alevin (v1.20.2)[ref] is run in sketch mode, which does pseudoalignment with structural constraints. Library type is defined as ISR & protocol is defined by defining the read geometry as read 2, position 1 – 15, barcode geometry as read 1, position 1 – 16 and UMI geometry as read 1, position 17 – 26.

#### RNA \& HTO
The subsequent mapping steps are caried out within the alevin-fry framework (v0.8.2)[ref] (an extension of salmon alevin).

First, a permit list of valid cell barcodes is generated with ***alevin-fry generate-permit-list***. Cell calling strategy is initially set as unfiltered mapping and quantification based on the February 2018 10x whitelist [ref]. Additionally, the expected orientation of alignments is set to both directions, in order to capture both sense and anti-sense reads.

The valid barcodes from the previous step are then corrected, using ***alevin-fry collate***, according to the information given in the permit list step.

Finally, the reads are quantified and output as a count matrix using ***alevin-fry quant***. The resolution strategy, by which molecules will be counted, is set to *cr-like*. *Cr-like* is currently the only strategy within the alevin-fry framework that supports splice aware mapping.

### Post-quantification filtering
#### Loading data
Counts are read into R with the ***load_fry*** function from the alevin-fry tutorial [ref].
- `RNA counts` [unspliced + spliced + ambiguous]
- `Spliced counts` [spliced + ambiguous]
- `Unspliced counts` [unspliced]
- `HTO counts`

#### Cell calling
Cell calling is done on the unfiltered count matrix of the RNA counts (U+S+A). First, barcode ranks are found with the ***barcodeRanks*** function from the DropletUtils package ({version})[ref]. Hartigan’s dip test is performed on counts above the knee point calculated above. If the resulting P-value >= 0.05, the pipeline moves forward with cells above the calculated inflection. If the resulting P-value is < 0.05, the hypothesis that the data is unimodal is rejected and a bimodal model is then fitted to the counts from cells above the original inflection with the ***Mclust*** function from the mclust package ({version})[ref]. Barcode ranks are then calculated again on the subset of cells assigned to mod 2 from before and the pipeline moves forward with cells above the re-calculated inflection. Optional a manual cutoff can be used to set a lower bound (in the ***barcodeRanks*** function) on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets.

#### Cell classification
##### HTO library
If HTOs are **not** a part of the experiment:
Seurat objects is converted to a singleCellExperiment object. RNA counts are log-normalized and PCA is performed.
An individual doublet score for each cell is then computed with the computeDoubletDensity function from the scDblFinder package.

If HTOs are a part of the experiment:
Cells overlapping between RNA and HTO library are imported into Seurat ({version})[ref]. The HTO read count assay is first CLR-normalized across cells (margin = 2) to account for the sequencing depth variation. For each HTO, the normalized HTO data is then modelled as a mixture of two normal distributions (a negative cluster from background noise and a positive cluster from true HTO signal) using the mclust R package[ref]. The normalized HTO read counts is plotted independently to visually examine the distribution and cells with higher number of normalized HTO read count than that of the negative HTO cluster are assigned with the corresponding HTO label. By modelling each HTO independently, our approach takes the cell hashing antibody affinity variation into consideration. Cells are finally demultiplexed into sample identity, where cells without HTO assignment, with only one HTO assignment and with more than one HTO assignment are labelled as negative, singlets and doublets, respectively.

##### Inferential doublets
Seurat objects is converted to a singleCellExperiment ({version})[ref] object. RNA counts are log-normalized and PCA is performed.
Intra-sample doublets are inferred from inter-sample doublets with the ***recoverDoublets*** function from the scDblFinder package ({version})[ref]. Afterwards, elbowcalling on the proportion vector (proportion of neighbours that are doublets) is used to label doublets from high doublet-density microclusters.

#### Output
Object is output into two unfiltered & two filtered Seurat objects, a standard object and a velocity object (containing spliced and unspliced assay). Unfiltered objects contain all called droplets (RNA calling) with all metadata generated throughout proccessing. Filtered objects contain only singlets and has a more concise set of metadata columns.
All objects are proccessed by running ***SCTransform*** on the RNA assay (*method = 'qpoisson'*). This creates the single cell transform assay (SCT), which is used downstream for visualization. The following commands are then run on the SCt assay.
RunPCA()
***maxLikGlobalDimEst*** is run on 50 PCA dimentions to estimate optimal range of dimentions for globaltSNE and UMAP.
RunTSNE()
RunUMAP()
FindNeighbors()
FindClusters()

#### Aggregation
If the experiment contains more than one reaction, all reactions are merged together (filtered and unfiltered, seperatly) and procced and saved in the same maner as the individual reaction objects
After processing individual sequence reactions, all results are aggregated into a single object (except if the experiment only has one reaction).
