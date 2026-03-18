#' Convert a weight vector into cumulative distribution bounds
#'
#' Function takes a numeric vector of weights and constructs lower and upper
#' cumulative distribution function (CDF) bounds, along with a crisp variant.
#' The weights are normalized to sum to 1 and then cumulatively summed to form
#' interval bounds.
#'
#' @param w A numeric vector of non-negative non-normalized weights. Must have a positive sum.
#'
#' @return A list with three numeric vectors:
#' \describe{
#'   \item{cdf_lower}{Lower bounds of the cumulative distribution.}
#'   \item{cdf_upper}{Upper bounds of the cumulative distribution.}
#'   \item{cdf_crisp}{Crisp version of the CDF (identical to \code{cdf_lower}).}
#' }
#'
#' @export
sample_to_bounds_w <- function(w) {
  w_norm <- w / sum(w)
  cum_w <- cumsum(w_norm)

  lower <- c(0, cum_w[-length(w_norm)], cum_w[length(w_norm)-1])
  upper <- c(cum_w[1], cum_w)
  crisp <- c(0, cum_w)

  return(list(cdf_lower = lower, cdf_upper = upper, cdf_crisp = crisp))
}
