top_barcodes <- unknown_barcodes(stats_folder, n = 10)

prepare_reference <- function(data_folder = "data/") {
  
  rna <- readr::read_csv(
    file = file.path(data_folder, "Dual_Index_Kit_TT_Set_A.csv"),
    skip = 3,
    col_select = c("name"   = "index_name",
                   "index"  = "index(i7)",
                   "index2" = "index2_workflow_b(i5)"),
    col_types = c("name"   = readr::col_character(),
                  "index"  = readr::col_character(),
                  "index2" = readr::col_character()))
  
  hto <- readr::read_csv(
    file = file.path(data_folder, "TruSeq_I7_indexes.csv"),
    col_select = c("name"   = "index_name",
                   "index"  = "index(i7)"),
    col_types = c("name"   = readr::col_character(),
                  "index"  = readr::col_character())
  )
  hto[["index2"]] <- NA
  
  no_index <- tibble::tibble(
    name = c("no index", "no signal"),
    "index" = c(NA, "GGGGGGGGGG"),
    "index2" = c("AGATCTCGGT", "GGGGGGGGGG")
  )
  reference <- rbind(rna, hto, no_index)
  reference
}

guess_index_id <- function(idx, reference_seq, 
                           reference_name, max_distance = 2) {
  quick_haming <- function(x, ref) {
    x_split <- strsplit(x, split = "")[[1]]
    ref_split <- do.call("cbind", strsplit(ref, split = ""))
    
    dist <- x_split != ref_split
    dist[is.na(dist)] <- TRUE
    dist <- colSums(dist)
    which((dist < max_distance) & (dist == min(dist)))
  }
  quick_haming <- Vectorize(quick_haming, vectorize.args = "x")
  
  lapply(quick_haming(idx, reference_seq), \(x) reference_name[x])
}

plot_unknown_barcodes <- function(top_barcodes, reference = NULL) {
  # If we have both index and index2 combine them
  plot_data <- top_barcodes |>
    tidyr::unite(col = "idx",
          dplyr::starts_with("index"),
          sep = "+",
          remove = FALSE, 
          na.rm = TRUE)
  
  # We want the most common index across the lanes at the top
  total_reads <- plot_data |>
    dplyr::group_by(idx) |>
    dplyr::summarise(reads = sum(reads))
  idx_order <- forcats::fct_reorder(total_reads[["idx"]],
                                    total_reads[["reads"]]) |>
    levels()
  
  plot_data <- plot_data |>
    dplyr::mutate(idx = factor(x = idx, levels = idx_order, ordered = TRUE))
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(y = idx, x = reads)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_x_continuous(name = "Reads", labels = exponent_format()) +
    ggplot2::ylab("Index") +
    ggplot2::facet_wrap(~ Lane) +
    theme_comuneqaid() +
    ggplot2::theme(legend.position = 'bottom',
                   axis.text.x = ggplot2::element_text(angle = 45, 
                                                       vjust = 1, 
                                                       hjust=1))
    
  if (!is.null(reference)) {
    # We try to guess what the indexes might be
    if (!("name" %in% colnames(reference))) {
      warning("reference lacks index names, skipping index guessing")
    } else {
      # Remove NA columns to handle missing indexes
      no_na <- plot_data[colSums(!is.na(plot_data)) > 0]
      
      index_cols <- grep(pattern = "index", colnames(no_na), value = TRUE)
      guesses <- vector(mode = "list", length = length(index_cols))
      names(guesses) <- index_cols
      for (i in index_cols) {
        if (!(i %in% colnames(reference))) {
          warning("reference lacks ", i, 
                  ", skipping index guessing for this index")
        } else {
          ids <- guess_index_id(no_na[[i]], reference[[i]], reference[["name"]])
          # Multiple indexes can be close, we paste them together
          ids <- lapply(ids, paste0, collapse = ", ") |> unlist()
          ids[ids == ""] <- "no match"
          names(ids) <- no_na[["idx"]]
          guesses[[i]] <- ids
        }
      }
      
    }
    guess_df <- data.frame(
      idx = names(guesses[[1]]),
      guess = Reduce(\(x, y){paste(x, y, sep = " + ")}, guesses)
    )
    guess_df <- guess_df |>
      dplyr::mutate(idx = factor(x = idx, levels = idx_order, ordered = TRUE))
    
      
    p_guess <- ggplot2::ggplot(
      data = guess_df, 
      mapping = ggplot2::aes(x = 1, y = idx, label = guess)) +
      ggplot2::geom_text() +
      ggplot2::theme_void()
    
    p <- patchwork::wrap_plots(p, p_guess, nrow = TRUE, widths = c(0.7, 0.3))
  }
  p
}

