# WARNING - Generated by {fusen} from dev/flat_first.Rmd: do not edit by hand

#' hw convolution of the h interaction function with the kernel of the 2
#' processus
#'
#' @param data results of compute_hawkes_histogram
#' @param width a vector of beds interval width (in the same orde as the data columns)
#' @param K (default: 10) size of the histogram bins
#' @param delta (default: 1e4) max range to consider
#' 
#' @return
#' tibble
#' @export
#'
#' @importFrom tidyr unnest
#' @importFrom dplyr mutate filter
#' @importFrom ggplot2 geom_vline geom_hline geom_ribbon aes facet_wrap
#' @importFrom ggplot2 theme_bw element_text
plot_convolution <- function(data, width, K = 10, delta = 1e4){
  data <- convolve_with_kernel(
    data = data, width = width, K = K, delta = delta
    ) |> 
    tidyr::unnest(convolution) |>
    dplyr::mutate(
      title = paste(
        name, "->", vs_name, "|", paste(
          setdiff(beds$names, c(name, vs_name)), collapse = ", "
        )
      )
    )
  ggplot2::ggplot() +
  ggplot2::geom_vline(
    xintercept = 0, color = "gray", size = 0.5, linetype = "dashed") +
  ggplot2::geom_hline(yintercept = 0, color = "gray", size = 0.5) +
  ggplot2::geom_ribbon(
    data = data |> dplyr::filter(convolution >= 0),
    ggplot2::aes(x = position, ymin = 0, ymax = convolution),
    fill = "#EABDBA", color = "black") +
  ggplot2::geom_ribbon(
    data = data |> dplyr::filter(convolution <= 0),
    ggplot2::aes(x = position, ymin = convolution, ymax = 0),
    fill = "#8EBBF5", color = "black") +
  ggplot2::facet_wrap(~title, scales = "free_y") +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(
    angle = 45, vjust = 1, hjust = 1
  ))
}
