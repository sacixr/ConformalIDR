test_that("sample_to_bounds_w works", {
  n <- 1000
  w <- seq(1, n)
  cdfs <- sample_to_bounds_w(w)

  expect_equal(length(cdfs), 3)

  names <- names(cdfs)
  expect_true(all(
    names[1] == "cdf_lower",
    names[2] == "cdf_upper",
    names[3] == "cdf_crisp"
  ))

  expect_true(all(
    sapply(cdfs, function(x) x >= 0)
  ))

  expect_true(all(
    sapply(cdfs, function(x) length(x) == n+1)
  ))
})

test_that("sample_to_bounds_w produces the same when weights are equal", {
  n <- 5
  w <- rep(1, n) # equal weights
  cdfs <- sample_to_bounds(n)
  cdfs_w <- sample_to_bounds_w(w)
  expect_equal(cdfs, cdfs_w)
})
