#' Get bcl-convert stats
#'
#' @param stats_folder string, path to bcl-convert stats folder
#'
#' @return
#' @export
#'
#' @examples
demux_stats <- function(stats_folder) {
  stats_file <- file.path(stats_folder, "Demultiplex_Stats.csv")
  if (!file.exists(stats_file)) {
    stop(stats_file, " does not exist.")
  }
  reads <- readr::read_csv(stats_file, 
                           col_select = c(
                             "Lane"          = "Lane",
                             "index"         = "SampleID",
                             "Reads"         = "# Reads",
                             "Perfect Index" = "# Perfect Index Reads",
                             "One Mismatch"  = "# One Mismatch Index Reads",
                             "Two Mismatch"  = "# Two Mismatch Index Reads"
                           ), 
                           col_types = c(
                             "Lane"          = readr::col_character(),
                             "SampleID"      = readr::col_character(),
                             "Reads"         = readr::col_integer(),
                             "Perfect Index" = readr::col_integer(),
                             "One Mismatch"  = readr::col_integer(),
                             "Two Mismatch"  = readr::col_integer())
  )
  
  types <- merge_sheets("reaction2library", "library_sheet", config = config) |>
    dplyr::select(-c("library_id"))
  undetermined <- data.frame("reaction_id"  = "Undetermined",
                             "library_type" = "Unknown",
                             "index"        = "Undetermined")
  types <- rbind(types, undetermined)
  demux_stats <- merge.data.frame(types, reads, by = "index")
  demux_stats
}


#' Plot demultiplexing stats
#'
#' @param demux_table output from demux_stats()
#'
#' @return
#' @export
#'
#' @examples
plot_demux_stats <- function(demux_table) {
  plot_data <- demux_table |>
    dplyr::select(-c("Reads")) |>
    tidyr::pivot_longer(cols = c("Perfect Index", 
                                 "One Mismatch", 
                                 "Two Mismatch")) |>
    dplyr::mutate(name = factor(name, levels = c("Two Mismatch", 
                                                 "One Mismatch", 
                                                 "Perfect Index"))) |>
    dplyr::mutate(index = ifelse(index == "Undetermined", "", index)) |>
    dplyr::filter(value > 0)
  
  # We want 10x/RNA and HTO in the beginning and Unknown in the end.
  # Other types may come in the future.
  order <- factor(unique(plot_data[["library_type"]]))
  order <- relevel(order, "Unknown")
  levels(order) <- rev(levels(order))
  for (i in c("hto", "HTO", "10x", "RNA")) {
    if (i %in% levels(order)) {
      order <- relevel(order, i)
    }
  }
  
  plot_data[["library_type"]] <- factor(plot_data[["library_type"]],
                                        levels = levels(order),
                                        ordered = TRUE)
  
  cols <- softPallet(3)
  
  fill_scale <- ggplot2::scale_fill_manual(
    name = ggplot2::element_blank(),
    values = c(
      "Perfect Index" = cols[1], 
      "One Mismatch"  = cols[2], 
      "Two Mismatch"  = cols[3]
    ),
    breaks = c(
      "Perfect Index", "One Mismatch", "Two Mismatch"
    ),
    drop = TRUE)
  
  ggplot2::ggplot(plot_data, 
                  ggplot2::aes(x = value, y = reaction_id, fill = name)) +
    ggplot2::geom_bar(position = ggplot2::position_stack(), stat = "identity") +
    ggplot2::geom_text(x = 0, 
                       mapping = ggplot2::aes(label = index), 
                       colour = "black", 
                       hjust = 0) +
    ggplot2::scale_x_continuous(name = "Reads", labels = exponent_format()) +
    ggplot2::ylab("Reaction ID") +
    ggplot2::facet_grid(library_type ~ Lane, 
                        scales = "free_y", 
                        space = "free_y") +
    fill_scale +
    theme_comuneqaid() +
    ggplot2::theme(legend.position = 'bottom',
                   axis.text.x = ggplot2::element_text(angle = 45, 
                                                       vjust = 1, 
                                                       hjust=1))
}

#' Get unknown barcodes
#'
#' @param stats_folder string, path to bcl-convert stats folder
#'
#' @return
#' @export
#'
#' @examples
unknown_barcodes <- function(stats_folder, n = 10) {
  barcodes_file <- file.path(stats_folder, "Top_Unknown_Barcodes.csv")
  if (!file.exists(barcodes_file)) {
    stop(barcodes_file, " does not exist.")
  }
  
  # Use starts_with to handle cases where only one index is used.
  readr::read_csv(barcodes_file,
                  col_select = c("Lane" = "Lane",
                                 tidyselect::starts_with("index"), 
                                 "reads"  = '# Reads'),
                  col_types = c("index"  = readr::col_character(), 
                                "index2" = readr::col_character(), 
                                "reads"  = readr::col_integer())) |>
    dplyr::group_by(Lane) |>
    dplyr::slice_max(reads, n = n)
}

#' Read and format index files
#'
#' @param data_folder string, path to folder where known index sequences are
#' stored.
#'
#' @return
#' @export
#'
#' @examples
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

#' Guess index ID
#'
#' @param idx string, indexes sequenced
#' @param reference_seq known indexes
#' @param reference_name names of known indexes
#' @param max_distance maximum hamming distance, default = 2
#'
#' @return
#' @export
#'
#' @examples
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

#' Plot unknown barcodes
#'
#' @param top_barcodes data.frame prepared by `unknown_barcodes`
#' @param reference optional data.frame prepared by `prepare_reference`
#'
#' @return
#' @export
#'
#' @examples
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
