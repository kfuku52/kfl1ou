test_that("adjust_data repairs short-edge non-ultrametric trees", {
  tree <- ape::read.tree(text = "((A:0,B:0):1,(C:0.5,D:0.5001):0.5);")
  Y <- matrix(c(-1, 0, 1, 2), ncol = 1, dimnames = list(tree$tip.label, "trait1"))

  dat <- capture_silently(
    adjust_data(
      tree,
      Y,
      normalize = FALSE,
      quietly = FALSE,
      min.edge.length = 1e-3
    )
  )

  expect_true(isTRUE(ape::is.ultrametric(dat$tree)))
  expect_true(all(dat$tree$edge.length >= 1e-3))

  fit <- capture_silently(
    estimate_shift_configuration(
      dat$tree,
      dat$Y,
      max.nShifts = 1,
      criterion = "AICc",
      quietly = TRUE
    )
  )
  expect_s3_class(fit, "l1ou")
})

test_that("adjust_data drops all-missing tips and invariant traits", {
  env <- new.env(parent = emptyenv())
  data("lizard.tree", package = "kfl1ou", envir = env)
  data("lizard.traits", package = "kfl1ou", envir = env)

  keep <- env$lizard.tree$tip.label[seq_len(8)]
  tree <- ape::drop.tip(env$lizard.tree, setdiff(env$lizard.tree$tip.label, keep))
  Y <- as.matrix(cbind(env$lizard.traits[keep, 1], const = 1))
  rownames(Y) <- keep
  Y[1, ] <- NA_real_

  dat <- capture_silently(adjust_data(tree, Y, normalize = FALSE, quietly = FALSE))

  expect_equal(dat$removed.tips, keep[[1]])
  expect_equal(dat$removed.traits, "const")
  expect_equal(dim(dat$Y), c(7, 1))
  expect_false(keep[[1]] %in% rownames(dat$Y))
})

test_that("adjust_data lets univariate inputs with missing tips proceed", {
  env <- new.env(parent = emptyenv())
  data("lizard.tree", package = "kfl1ou", envir = env)
  data("lizard.traits", package = "kfl1ou", envir = env)

  keep <- env$lizard.tree$tip.label[seq_len(8)]
  tree <- ape::drop.tip(env$lizard.tree, setdiff(env$lizard.tree$tip.label, keep))
  Y <- as.matrix(env$lizard.traits[keep, 1, drop = FALSE])
  Y[1, 1] <- NA_real_

  dat <- capture_silently(adjust_data(tree, Y, normalize = FALSE, quietly = FALSE))
  fit <- suppressWarnings(
    capture_silently(
      estimate_shift_configuration(
        dat$tree,
        dat$Y,
        max.nShifts = 1,
        criterion = "AICc",
        quietly = TRUE
      )
    )
  )

  expect_equal(dat$removed.tips, keep[[1]])
  expect_s3_class(fit, "l1ou")
})

test_that("alpha_upper_bound ignores isolated tiny external branches", {
  tree <- ape::read.tree(
    text = "((A:1e-12,B:1e-12):0.999999999999,((C:0.5,D:0.5):0.25,(E:0.5,F:0.5):0.25):0.25);"
  )

  expect_true(isTRUE(ape::is.ultrametric(tree)))
  expect_equal(kfl1ou:::alpha_upper_bound(tree), log(2) / 0.5, tolerance = 1e-12)
})
