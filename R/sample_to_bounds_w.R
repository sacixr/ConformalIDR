#' Produce a weighted CDF.
#'
#' @param w A list of normalized (between 0 and 1, summing to 1) weights, of length n+1.
#'
#' @returns A list containing the lower, upper and crisp CDF.
sample_to_bounds_w <- function(w) {
  n <- length(w)
  cum_w <- cumsum(w_norm)

  lower <- c(0, cum_w[-n])
  upper <- cum_w
  crisp <- c(0, cum_w[-n])

  return(list(cdf_lower = lower, cdf_upper = upper, cdf_crisp = crisp))
}
