# WARNING - Generated by {fusen} from dev/flat_first.Rmd: do not edit by hand

#' helper function for preprocess_bed to quickly load bed files
#'
#' @param path bed file path
#' @return
#' tibble
#' @importFrom vroom vroom
#' @noRd
read_bed <- function(path) {
  vroom::vroom(
    path,
    col_names = c("chr", "start", "stop"),
    col_types = list(chr = "c", start = "i", stop = "i"),
    col_select = c("chr", "start", "stop"),
  )
}
