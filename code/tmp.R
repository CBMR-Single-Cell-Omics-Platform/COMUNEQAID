top_barcodes <- unknown_barcodes(stats_folder, n = 10)

rna <- readr::read_csv(file = here::here("data/Dual_Index_Kit_TT_Set_A.csv"),
                       skip = 3,
                       col_select = c("name"   = "index_name",
                                      "index"  = "index(i7)",
                                      "index2" = "index2_workflow_b(i5)"),
                       col_types = c("name"   = readr::col_character(),
                                     "index"  = readr::col_character(),
                                     "index2" = readr::col_character()))

hto <- readr::read_csv(file = here::here("data/TruSeq_I7_indexes.csv"),
                       col_select = c("name"   = "index_name",
                                      "index"  = "index(i7)"),
                       col_types = c("name"   = readr::col_character(),
                                     "index"  = readr::col_character())
                       )
hto[["index2"]] <- NA

no_index <- tibble(name = c("no index", "no signal"),
                   "index" = c(NA, "GGGGGGGGGG")
                   "index2" = c("AGATCTCGGT", "GGGGGGGGGG")
)
reference <- rbind(rna, hto, no_index)

guess_index_id <- function(idx, reference_seq, reference_name) {
  quick_haming <- function(x, ref) {
    x_split <- strsplit(x, split = "")[[1]]
    ref_split <- do.call("cbind", strsplit(ref, split = ""))
    
    dist <- x_split != ref_split
    dist[is.na(dist)] <- TRUE
    dist <- colSums(dist)
    which((dist < 2) & (dist == min(dist)))
  }
  quick_haming <- Vectorize(quick_haming, vectorize.args = "x")
  
  lapply(quick_haming(idx, reference_seq), \(x) reference_name[x])
}

plot_unknown_barcodes <- function(top_barcodes, reference = NULL) {
  # If we have both index and index2 combine them
  plot_data <- top_barcodes |>
    tidyr::unite(col = "idx",
          starts_with("index"),
          sep = "+",
          remove = FALSE)
  
  # We want the most common index across the lanes at the top
  total_reads <- plot_data |>
    group_by(idx) |>
    summarise(reads = sum(reads))
  idx_order <- fct_reorder(total_reads[["idx"]],
                           total_reads[["reads"]]) |>
    levels()
  
  plot_data <- plot_data |>
    mutate(idx = factor(x = idx, levels = idx_order, ordered = TRUE))
  
  p <- ggplot(plot_data, ggplot2::aes(y = idx, x = reads)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_x_continuous(name = "Reads", labels = exponent_format()) +
    ggplot2::ylab("Index") +
    ggplot2::facet_wrap(~ Lane) +
    theme_comuneqaid() +
    ggplot2::theme(legend.position = 'bottom')
  
  if (!is.null(reference)) {
    # We try to guess what the indexes might be
    plot_data |>
      mutate(across())
    
  }
}

