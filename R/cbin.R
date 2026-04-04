#' Conformal binning
#'
#' Derive conformal predictive systems based on binning. Either k-means clustering or
#' a regression tree is used to bin the covariates, and separate predictive systems are
#' then constructed for each bin.
#'
#'
#' @inheritParams cops
#' @param binning binning method used to construct the bins; either \code{"k-means"} (default), \code{"tree"}, \code{"hclust"} or \code{"dbscan"}
#' @param k number of bins (for k-means clustering)
#' @param cp complexity parameter (for regression tree binning)
#' @param weights non-negative weights assigned to each training covariate and the
#' new covariate. If omitted, equal weights are used (the default).
#' @param cluster_on binning on covariates `x` or observations `y`
#' @param eps epsilon parameter (for DBSCAN binning)
#' @param minPts minPts parameter (for DBSCAN binning)
#'
#'
#' @return
#' \code{conformal_bin()} returns a list of objects of class \code{"cops"} containing
#' the conformal binning fit. See \code{\link{cops}} for details.
#'
#'
#' @details
#'
#' Conformal binning generates conformal predictive systems based on binning algorithms.
#' Binning algorithms partition the covariate space into a finite number of bins,
#' find the bin in which the new covariate lies, and then use the empirical distribution
#' of the labels in this bin as a prediction for the new label. Conformal binning
#' provides a means to generate predictive system bands that are guaranteed to contain
#' an auto-calibrated forecast distribution. The information in the forecast depends
#' on what binning is performed.
#'
#' The bins can either be constructed using k-means clustering (\code{binning = "k-means"}),
#' or using a regression tree (\code{binning = "tree"}). k-means finds the bins that minimise
#' the variance of the covariates assigned to the same bin, while the regression tree
#' partitions the covariates such that the variance of the corresponding labels
#' assigned to the same bin is minimised.
#'
#' Since the regression tree depends on the labels corresponding to the covariates,
#' it can only be applied in the split conformal setting, and not in the full conformal setting.
#' The regression tree is fit to the estimation data \code{x_set} and \code{y_set}.
#' However, k-means clustering can be applied in both the split and full conformal settings.
#' If \code{x_set} is provided, then the bins are estimated using this estimation data.
#' Otherwise, the bins are estimated from the training data \code{x} (and the new covariate \code{x_out}).
#'
#' The k-means clustering approach requires a choice for the number of bins, which can be specified
#' using the argument \code{k}. While the regression tree cannot return a specific number of bins,
#' the number of bins can indirectly be influenced by the choice of a complexity parameter \code{cp}.
#' The complexity parameter specifies the minimum improvement in fit required for a split in the
#' regression tree to be performed; a smaller complexity parameter therefore results in a more complex
#' tree with more bins.
#'
#' The k-means algorithm is implemented using \code{\link[stats]{kmeans}} and \code{\link[clue]{cl_predict}},
#' while the regression tree leverages \code{\link[rpart]{rpart}}.
#' Further details can be found in the respective help pages.
#'
#'
#' @references
#'
#' Allen, S., Gavrilopolous, G., Henzi, A. and J. Ziegel (2025):
#' `In-sample calibration yields conformal calibration guarantees',
#' \emph{arXiv pre-print} arXiv:2503.03841
#' \doi{10.48550/arXiv.2503.03841}
#'
#'
#' @seealso \code{\link{cops}} \code{\link{cidr}} \code{\link{lspm}}
#'
#' @author Sam Allen
#'
#' @examples
#'
#' n <- 1000
#' x <- rnorm(n)
#' y <- rnorm(n, x, exp(x))
#'
#' N <- 100
#' x_out <- rnorm(N)
#' y_out <- rnorm(N, x_out, exp(x_out))
#' plot(x_out, y_out)
#'
#' ## full conformal
#'
#' fit <- conformal_bin(x, y, x_out, y_out, k = 10)
#' fit1 <- conformal_bin(x, y, x_out, y_out)
#' fit5 <- conformal_bin(x, y, x_out, y_out, k = 5)
#' fit100 <- conformal_bin(x, y, x_out, y_out, k = 100)
#'
#' # compare for different k
#' crps_vec <- c(
#'   "k = 1" = sapply(fit1, function(x) x$crps) |> mean(),
#'   "k = 5" = sapply(fit5, function(x) x$crps) |> mean(),
#'   "k = 10" = sapply(fit, function(x) x$crps) |> mean(),
#'   "k = 100" = sapply(fit100, function(x) x$crps) |> mean()
#' )
#' print(crps_vec)
#'
#'
#' ## split conformal (with k means)
#'
#' N0 <- 500
#' x_est <- rnorm(N0)
#' y_est <- rnorm(N0, x_est, exp(x_est))
#'
#' fit <- conformal_bin(x, y, x_out, y_out, x_est, y_est, k = 10)
#' fit1 <- conformal_bin(x, y, x_out, y_out, x_est, y_est)
#' fit5 <- conformal_bin(x, y, x_out, y_out, x_est, y_est, k = 5)
#' fit100 <- conformal_bin(x, y, x_out, y_out, x_est, y_est, k = 100)
#'
#' # compare for different k
#' crps_vec_sp <- c(
#'   "k = 1" = sapply(fit1, function(x) x$crps) |> mean(),
#'   "k = 5" = sapply(fit5, function(x) x$crps) |> mean(),
#'   "k = 10" = sapply(fit, function(x) x$crps) |> mean(),
#'   "k = 100" = sapply(fit100, function(x) x$crps) |> mean()
#' )
#' print(crps_vec_sp)
#'
#'
#' ## split conformal (with regression tree)
#'
#' fit <- conformal_bin(x, y, x_out, y_out, x_est, y_est, binning = "tree")
#' fit001 <- conformal_bin(x, y, x_out, y_out, x_est, y_est, binning = "tree", cp = 0.0001)
#' fit01 <- conformal_bin(x, y, x_out, y_out, x_est, y_est, binning = "tree", cp = 0.001)
#' fit05 <- conformal_bin(x, y, x_out, y_out, x_est, y_est, binning = "tree", cp = 0.005)
#' fit1 <- conformal_bin(x, y, x_out, y_out, x_est, y_est, binning = "tree", cp = 0.1)
#'
#' # compare for different cp
#' crps_vec_sp_tr <- c(
#'   "cp = 0.0001" = sapply(fit001, function(x) x$crps) |> mean(),
#'   "cp = 0.001" = sapply(fit01, function(x) x$crps) |> mean(),
#'   "cp = 0.005" = sapply(fit05, function(x) x$crps) |> mean(),
#'   "cp = 0.01" = sapply(fit, function(x) x$crps) |> mean(),
#'   "cp = 0.1" = sapply(fit1, function(x) x$crps) |> mean()
#' )
#' print(crps_vec_sp_tr)
#'
#' print(c(crps_vec_sp, crps_vec_sp_tr))
#'
#' @name cbin
#' @export
conformal_bin <- function(x, y, x_out, y_out = NULL, x_est = NULL, y_est = NULL, online = FALSE, weights = NULL, binning = c("kmeans", "tree", "hclust", "dbscan"), cut, k = 1, cp = 0.01, eps = 0.5, minPts = 1, cluster_on = c("x", "y")) {
  binning <- match.arg(binning)
  cluster_on <- match.arg(cluster_on)
  check_cops_args(x, y, x_out, y_out, x_est, y_est, "cbin", online, weights, binning = binning, k = k, cp = cp)

  eval <- !is.null(y_out)
  split <- !is.null(x_est)

  if (split) {
    # split conformal - estimate clusters on the estimation sample

    if (binning == "kmeans") {
      if (cluster_on == "x") {
        if (k == 0) { # if k is not provided, estimate an optimal k
          k <- elbow_method(as.matrix(x_est), NULL, "kmeans")
        }
        out <- stats::kmeans(as.matrix(x_est), k)
        tr_cl <- clue::cl_predict(out, as.matrix(x)) |> as.vector()
        ts_cl <- clue::cl_predict(out, as.matrix(x_out)) |> as.vector()
      } else {
        train_control <- trainControl(method = "cv", number = 5, savePredictions = TRUE)
        model_est <- train(y_est ~., data = data.frame(x_est, y_est), method = "knn", tuneGrid = data.frame(k = 3), trControl = train_control)$pred$pred
        out <- stats::kmeans(as.matrix(model_est), k)

        model_y <- train(y ~., data = data.frame(x, y), method = "knn", tuneGrid = data.frame(k = 3), trControl = train_control)$pred$pred
        tr_cl <- clue::cl_predict(out, as.matrix(model_y)) |> as.vector()

        model_out <- train(y_out ~., data = data.frame(x, y), method = "knn", tuneGrid = data.frame(k = 3), trControl = train_control)$pred$pred
        ts_cl <- clue::cl_predict(out, as.matrix(model_out)) |> as.vector()
      }
    } else if (binning == "tree") {
      dat <- data.frame(y = y_est, x = x_est)
      tree <- rpart::rpart(y ~ x, data = dat, method = "anova", control = rpart::rpart.control(cp = cp))
      leaf_ind <- tree$frame$var == "<leaf>"
      k <- sum(leaf_ind)
      tree$frame$yval[!leaf_ind] <- NA
      tree$frame$yval[leaf_ind] <- 1:k
      tr_cl <- stats::predict(tree, newdata = data.frame(x = x)) |> unname()
      ts_cl <- stats::predict(tree, newdata = data.frame(x = x_out)) |> unname()
    } else if (binning == "hclust") {
      hclust_cl <- opt_hclust(x_est, "ward.D2", cut, k)
      x_est <- as.matrix(x_est)
      k <- length(unique(hclust_cl))
      tr_cl <- class::knn(x_est, as.matrix(x), cl = factor(hclust_cl), k = k) # since hclust doesn't naturally classify, we use knn
      ts_cl <- class::knn(x_est, as.matrix(x_out), cl = factor(hclust_cl), k = k)
    } else if (binning == "dbscan") {
      x_est <- as.matrix(x_est)
      dbscan_cl <- fpc::dbscan(x_est, eps = eps, MinPts = minPts)
      k <- length(unique(dbscan_cl$cluster))
      tr_cl <- fpc::predict.dbscan(dbscan_cl, x_est, as.matrix(x))
      y <- y[tr_cl != 0] # remove outliers
      tr_cl <- tr_cl[tr_cl != 0]
      ts_cl <- fpc::predict.dbscan(dbscan_cl, x_est, as.matrix(x_out))
      x_out <- x_out[ts_cl != 0]
      ts_cl <- ts_cl[ts_cl != 0]
    }

    if (is.null(weights)) {
      points <- sapply(1:k, function(i) y[tr_cl == i] |> sort())
      points <- lapply(points, function(x) c(-Inf, x, Inf))
      points[[k + 1]] <- c(-Inf, y, Inf) |> sort()
      cdfs <- lapply(points, function(x) sample_to_bounds(length(x) - 1))
    } else {
      points <- vector("list", k + 1)
      weights_p <- vector("list", k + 1)

      for (i in 1:k) {
        ord <- y[tr_cl == i] |> order()
        # sort the weights in line with the points
        # append max weight for padding
        weights_p[[i]] <- c((weights[tr_cl == i])[ord], max(weights)+1)
        points[[i]] <- y[tr_cl == i] |> sort()
      }

      points <- lapply(points, function(x) c(-Inf, x, Inf))
      ord <- y |> order()
      points[[k+1]] <- c(-Inf, y[ord], Inf)
      weights_p[[k+1]] <- weights[ord] # sort for unconditional bin
      cdfs <- lapply(weights_p, function(i) sample_to_bounds_w(i))
    }

    eval <- !is.null(y_out)
    out_list <- lapply(seq_along(x_out), function(i) {
      bin <- ts_cl[i]
      # if a bin is not observed in the calibration set, return the unconditional distribution
      if (sum(y[tr_cl == bin]) == 0) bin <- k + 1
      out <- c(points = points[bin], cdfs[[bin]], list(x = x, y = y, x_out = x_out[i]))
      out <- structure(out, class = "cops")

      if (eval) {
        pcal <- pit(out, y_out[i])
        score <- crps(out, y_out[i])
        thick <- thickness(out)
        out <- c(out, list(y_out = y_out[i], pit = pcal, crps = score, thick = thick))
      }
      out <- c(out, bins = list(x = tr_cl, x_out = ts_cl[i]))
      structure(out, class = "cops")
    })


  } else {
    # full conformal - estimate clusters using the new training covariate x_n+1

    out_list <- lapply(seq_along(x_out), function(i) {
      xn1 <- c(x, x_out[i])
      n1 <- length(xn1)

      if (binning == "kmeans") {
        if (k == 0) {
          k <- elbow_method(as.matrix(xn1), NULL, "kmeans")
        }
        if (cluster_on == "x") {
          out <- stats::kmeans(as.matrix(xn1), k)
          tr_cl <- out$cluster[-n1]
          ts_cl <- out$cluster[n1]
        } else {
          yn1 <- knn.reg(train=xn1, y=y, k=3)$pred |> as.matrix()
          yn1[is.na(yn1)] <- 0
          out <- stats::kmeans(yn1, k)
          tr_cl <- out$cluster[-n]
          ts_cl <- out$cluster[n]
        }
      } else if (binning == "hclust") {
        hclust_cl <- opt_hclust(xn1, "ward.D2", cut, k)
        k <- length(unique(hclust_cl))
        tr_cl <- hclust_cl[-n1]
        ts_cl <- hclust_cl[n1]
      } else if (binning == "dbscan") {
        dbscan_cl <- fpc::dbscan(as.matrix(xn1), eps = eps, MinPts = minPts)
        k <- length(unique(dbscan_cl$cluster))
        tr_cl <- dbscan_cl$cluster[-n1]
        tr_cl <- tr_cl[tr_cl != 0]
        ts_cl <- dbscan_cl$cluster[n1]
      }

      if (is.null(weights)) {
        points <- sapply(1:k, function(cl) {
          if (sum(tr_cl == cl) == 0) {
            y |> sort()
          } else {
            y[tr_cl == cl] |> sort()
          }
        })
        points <- lapply(points, function(z) c(-Inf, z, Inf))
        cdfs <- lapply(points, function(z) sample_to_bounds(length(z) - 1))
      } else {
        points <- vector("list", k + 1)
        weights_p <- vector("list", k + 1)

        for (i in 1:k) {
          ord <- y[tr_cl == i] |> order()
          weights_p[[i]] <- c((weights[tr_cl == i])[ord], max(weights)+1)
          points[[i]] <- y[tr_cl == i] |> sort()
        }

        points <- lapply(points, function(x) c(-Inf, x, Inf))
        ord <- y |> order()
        points[[k+1]] <- c(-Inf, y[ord], Inf)
        weights_p[[k+1]] <- weights[ord]
        cdfs <- lapply(weights_p, function(i) sample_to_bounds_w(i))
      }

      out <- c(points = points[ts_cl], cdfs[[ts_cl]], list(x = x, y = y, x_out = x_out[i]))
      out <- structure(out, class = "cops")

      if (eval) {
        pcal <- pit(out, y_out[i])
        score <- crps(out, y_out[i])
        thick <- thickness(out)
        out <- c(out, list(y_out = y_out[i], pit = pcal, crps = score, thick = thick))
      }
      out <- c(out, bins = list(x = tr_cl, x_out = ts_cl))
      structure(out, class = "cops")
    })
  }

  return(out_list)
}
