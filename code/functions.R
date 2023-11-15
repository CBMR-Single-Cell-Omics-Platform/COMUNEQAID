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
  
  purrr::reduce(.x = config[table_names],
         .f = merge.data.frame)
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

#' Common theme for COMUNEQAID
#'
#' @param base_size numeric, base size
#' @param base_family string, text family
#'
#' @return
#' @export
#'
#' @examples
theme_comuneqaid <- function(base_size = 10, base_family = "Helvetica") {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) 
}
