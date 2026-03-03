# select the optimal numbers of clusters using a dendogram
dendo <- function(x, x_cl) {
  diffs <- diff(x_cl$height) # get differences between merges
  which_max_jump <- which.max(diffs) # index of the largest jump
  opt_height <- (x_cl$height[which_max_jump] + x_cl$height[which_max_jump + 1]) / 2 # optimal height is halfway of the largest max jump
}

elbow_method <- function(x, x_cl, method) {
  plot_val <- matrix(0, nrow=1, ncol=2) # matrix for x (k), and y (ratio)
  for (i in 2:10) { # for now, just test up to 10 clusters
    if (method == "hclust") {
      cut_val <- cutree(x_cl, i) # cut at i clusters
    } else {
      cut_val <- stats::kmeans(x, i)$cluster # perform kmeans for those clusters
    }
    cl_stats <- fpc::cluster.stats(x, cut_val) # fetch the between and within values
    # need some way to handle if avg.within is 0, and alt
    ratio <- cl_stats$average.between/cl_stats$average.within # plot the ratio
    plot_val <- rbind(plot_val, c(i, ratio))
  }
  spline <- stats::smooth.spline(plot_val[,1], plot_val[,2], spar=0.5) # form a smooth spline over values
  sec_deriv <- predict(spline$fit, plot_val[,2], deriv=2) # fetch the second derivative
  opt_k <- which.max(sec_deriv[["y"]]) # the max 2nd deriv is the "elbow" so our opt k
}

# simple hierarchical clustering function
# returns assigned clusters for x
opt_hclust <- function(x, method, cut, k) {
  if (!inherits(x, "dist")) {
    x <- dist(x)
    x[is.na(x)] <- Inf # replace NAs with inf
  }

  x_cl <- hclust(x, method=method) # apply hierarchical clustering
  if (k == 0) {
    if (cut == "dendo") {
      h <- dendo(x, x_cl)
      return(cutree(tree=x_cl, h=h))
    } else {
      k <- elbow_method(x, x_cl, "hclust")
    }
  }
  return(cutree(tree=x_cl, k=k)) # if k already set, cut at k clusters
}
