#' Print a POSIXct in a machine or pretty format
#'
#' @param x A POSIXct object
#' @param type string "pretty" or "machine"
#'
#' @return a string
#' @export
#'
#' @examples
#' format_time(lubridate::now(), "pretty")
#' format_time(lubridate::now(), "machine")
format_time <- function(x, type) {
  if (type == "pretty") {
    format(x, "%Y-%m-%d %H:%M:%S %Z")
  } else if (type == "machine") {
    format(x, "%Y-%m-%d_%H-%M-%S")
  } else {
    stop("Type must be 'pretty' or 'machine'")
  }
}

#' Path to project
#'
#' @param config config part of a snakemake object
#'
#' @return string, path to project home dir
#' @export
#'
#' @examples
#' \dontrun{
#' project_path()
#' }
project_path <- function(config = snakemake@config) {
  file.path(config[['project_path']],
            config[['scop_id']])
}

#' Path to project
#'
#' @param config config part of a snakemake object 
#'
#' @return string, path to FASTQ folder
#' @export
#'
#' @examples
fastq_path <- function(config = snakemake@config) {
  file.path(project_path(),
            snakemake@config[['fastq_path']])
}

#' Merge experimental tables
#'
#' @param ... names of sheets to merge
#' @param config config part of a snakemake object
#'
#' @return data.frame with merged tables
#' @details
#' Each table to be merged should have one column in common with the 
#' previous table. This is not checked for by the function, so be careful when
#' merging tables.
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' merge_sheets("reaction2library", "library_sheet", 
#' "library2sequencing", "sequencing_sheet")
#' }
merge_sheets <- function(..., config = snakemake@config) {
  table_names <- unlist(list(...))
  
  # Tables need to be present in config
  present <- table_names %in% names(config)
  if (!all(present)) {
    stop(paste0(table_names[!present], collapse = ", "), " not in config")
  }
  
  reduce(.x = config[table_names],
         .f = merge.data.frame)
}

#' Get bcl-convert stats
#'
#' @param stats_folder string, path to bcl-convert stats folder
#'
#' @return
#' @export
#'
#' @examples
demux_counts <- function(stats_folder) {
  demult_stats <- file.path(stats_folder, "Demultiplex_Stats.csv")
  if (!file.exists(demult_stats)) {
    stop(demult_stats, " does not exist.")
  }
  reads <- readr::read_csv(demult_stats, 
                           col_select = c("index" = "SampleID", 
                                          "reads" = "# Reads"), 
                           col_types = c("index" = readr::col_character(), 
                                         "reads" = readr::col_integer())) |>
    dplyr::group_by(index) |>
    dplyr::summarise(reads = sum(reads))
  
  types <- merge_sheets("reaction2library", "library_sheet", config = config) |>
    select(-c("library_id", "reaction_id"))
  types <- rbind(types, c("Undetermined", "Undetermined"))
  demux_stats <- merge.data.frame(types, reads)
  demux_stats
}


demux_counts <- function(stats_folder) {
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
    select(-c("library_id"))
  undetermined <- data.frame("reaction_id"  = "Undetermined",
                             "library_type" = "Unknown",
                             "index"        = "Undetermined")
  types <- rbind(types, undetermined)
  demux_stats <- merge.data.frame(types, reads, by = "index")
  demux_stats
}

demux_stats <- demux_counts(stats_folder)
plot_demux_stats <- function(demux_stats) {
  
  plot_data <- demux_stats |>
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
  
  fill_scale <- scale_fill_manual(
    name = element_blank(),
    values = c(
      "Perfect Index" = cols[1], 
      "One Mismatch"  = cols[2], 
      "Two Mismatch"  = cols[3]
    ),
    breaks = c(
      "Perfect Index", "One Mismatch", "Two Mismatch"
    ),
    drop = TRUE)
  
  ggplot(plot_data, aes(x = reaction_id, y = value, fill = name)) +
    geom_bar(position = position_stack(), stat = "identity") +
    coord_flip() +
    geom_text(y = 0, mapping = aes(label = index), colour = "black", hjust = 0) +
    scale_y_continuous(name = "Reads", labels = exponent_format()) +
    xlab("Reaction ID") +
    facet_grid(library_type ~ Lane, scales = "free_y", space = "free_y") +
    fill_scale +
    theme_comuneqaid() +
    theme(legend.position = 'bottom')
}


#' Get bcl-convert stats
#'
#' @param stats_folder string, path to bcl-convert stats folder
#'
#' @return
#' @export
#'
#' @examples
unknown_barcodes <- function(stats_folder) {
  barcodes_file <- file.path(stats_folder, "Top_Unknown_Barcodes.csv")
  if (!file.exists(barcodes_file)) {
    stop(barcodes_file, " does not exist.")
  }
  
  read_csv(barcodes_file,
           col_select = c("index"  = "index", 
                          "index2" = "index2", 
                          "reads"  = '# Reads'),
           col_types = c("index"  = readr::col_character(), 
                         "index2" = readr::col_character(), 
                         "reads"  = readr::col_integer())) |>
    group_by(index, index2) |> 
    summarise("reads" = sum(reads), .groups = "keep") |> 
    arrange(desc(reads))
}

#' Exponent format
#'
#' @return
#' @export
#'
#' @examples
exponent_format <- function() {
  function(x) {
    parse(text = gsub("e\\+?", " %*% 10^", scales::scientific_format()(x)))
  }
}


theme_comuneqaid <- function(base_size = 10, base_family = "Helvetica") {
  theme_bw(base_size = base_size, base_family = base_family) 
}
