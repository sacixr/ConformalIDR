# fix all these tests:
# instead of checking if just the first object has class "cops", check that they all do
test_that("full and split conformal kmeans works", {
  n <- 1000
  x <- rnorm(n)
  y <- rnorm(n, x, exp(x))

  N <- 100
  x_out <- rnorm(N)
  y_out <- rnorm(N, x_out, exp(x_out))

  N0 <- 100
  x_est <- rnorm(N0)
  y_est <- rnorm(N0, x_est, exp(x_est))

  full <- conformal_bin(x, y, x_out, y_out, binning = "kmeans", k = 10)
  split <- conformal_bin(x, y, x_out, y_out, x_est, y_est, binning = "kmeans", k = 10)

  # check that the output for both is as expected
  expect_true(all(sapply(full, inherits, "cops")))
  expect_gte(length(full), 100)
  expect_true(all(
    sapply(full, function(x) is.list(x) && length(x) <= 13)
  ))

  expect_true(all(sapply(split, inherits, "cops")))
  expect_gte(length(split), 100)
  expect_true(all(
    sapply(split, function(x) is.list(x) && length(x) <= 13)
  ))
})

test_that("full and split conformal hclust works", {
  n <- 1000
  x <- rnorm(n)
  y <- rnorm(n, x, exp(x))

  N <- 100
  x_out <- rnorm(N)
  y_out <- rnorm(N, x_out, exp(x_out))

  N0 <- 100
  x_est <- rnorm(N0)
  y_est <- rnorm(N0, x_est, exp(x_est))

  full <- conformal_bin(x, y, x_out, y_out, binning = "hclust", k = 10)
  split <- conformal_bin(x, y, x_out, y_out, x_est, y_est, binning = "hclust", k = 10)

  expect_true(all(sapply(full, inherits, "cops")))
  expect_gte(length(full), 100)
  expect_true(all(
    sapply(full, function(x) is.list(x) && length(x) <= 13)
  ))

  expect_true(all(sapply(split, inherits, "cops")))
  expect_gte(length(split), 100)
  expect_true(all(
    sapply(split, function(x) is.list(x) && length(x) <= 13)
  ))
})

test_that("full and split conformal dbscan clustering works", {
  n <- 1000
  x <- rnorm(n)
  y <- rnorm(n, x, exp(x))

  N <- 100
  x_out <- rnorm(N)
  y_out <- rnorm(N, x_out, exp(x_out))

  N0 <- 100
  x_est <- rnorm(N0)
  y_est <- rnorm(N0, x_est, exp(x_est))

  full <- conformal_bin(x, y, x_out, y_out, binning = "dbscan", k = 10)
  split <- conformal_bin(x, y, x_out, y_out, x_est, y_est, binning = "dbscan", k = 10)

  expect_true(all(sapply(full, inherits, "cops")))
  expect_lte(length(full), 100)
  expect_true(all(
    sapply(full, function(x) is.list(x) && length(x) <= 13)
  ))

  expect_true(all(sapply(split, inherits, "cops")))
  expect_lte(length(split), 100)
  expect_true(all(
    sapply(split, function(x) is.list(x) && length(x) <= 13)
  ))
})

test_that("full and split conformal y-binning using kmeans works", {
  n <- 1000
  x <- rnorm(n)
  y <- rnorm(n, x, exp(x))

  N <- 100
  x_out <- rnorm(N)
  y_out <- rnorm(N, x_out, exp(x_out))

  N0 <- 100
  x_est <- rnorm(N0)
  y_est <- rnorm(N0, x_est, exp(x_est))

  full <- conformal_bin(x, y, x_out, y_out, binning = "dbscan", k = 10, cluster_on = "y")
  split <- conformal_bin(x, y, x_out, y_out, x_est, y_est, binning = "dbscan", k = 10, cluster_on = "y")

  expect_true(all(sapply(full, inherits, "cops")))
  expect_lte(length(full), 100)
  expect_true(all(
    sapply(full, function(x) is.list(x) && length(x) <= 13)
  ))

  expect_true(all(sapply(split, inherits, "cops")))
  expect_lte(length(split), 100)
  expect_true(all(
    sapply(split, function(x) is.list(x) && length(x) <= 13)
  ))
})

test_that("weighted split conformal bin works", {
  n <- 1000
  x <- rnorm(n)
  y <- rnorm(n, x, exp(x))

  N <- 100
  x_out <- rnorm(N)
  y_out <- rnorm(N, x_out, exp(x_out))

  N0 <- 100
  x_est <- rnorm(N0)
  y_est <- rnorm(N0, x_est, exp(x_est))

  n_weights <- c(seq(1, n), seq(1,N))
  split_w <- conformal_bin(x, y, x_out, y_out, x_est, y_est, binning = "kmeans", k = 10, weights = n_weights)

  expect_true(all(sapply(split_w, inherits, "cops")))
  expect_lte(length(split_w), 100)
  expect_true(all(
    sapply(split_w, function(x) is.list(x) && length(x) <= 13)
  ))
})

