test_that("rescale_matrix_and_error rescales input_error with the trait matrix", {
  Y <- matrix(c(1, 2, 3, 4), ncol = 2)
  rownames(Y) <- c("a", "b")
  colnames(Y) <- c("trait1", "trait2")
  input_error <- matrix(c(1, 2, 3, 4), ncol = 2, dimnames = dimnames(Y))

  rescaled <- kfl1ou:::rescale_matrix_and_error(Y, input_error)

  scale_values <- sqrt(colSums(Y^2))
  scale_values[scale_values == 0] <- 1
  mult <- 0.1 * nrow(Y)
  expected_Y <- mult * scale(Y, center = TRUE, scale = scale_values)
  expected_input_error <- sweep(input_error, 2, (mult / scale_values)^2, FUN = "*")

  expect_equal(rescaled$Y, expected_Y, tolerance = 1e-12)
  expect_equal(rescaled$input_error, expected_input_error, tolerance = 1e-12)
})

test_that("normalize_input_error reorders rows and columns and expands shared errors", {
  dat <- small_lizard_data(n_tips = 6, traits = 1:2, normalize = FALSE)
  input_error <- matrix(
    seq_len(nrow(dat$Y)),
    ncol = 1,
    dimnames = list(rev(rownames(dat$Y)), NULL)
  )

  normalized <- kfl1ou:::normalize_input_error(dat$tree, dat$Y, input_error)

  shared_error <- stats::setNames(seq_len(nrow(dat$Y)), rev(rownames(dat$Y)))
  expected <- cbind(shared_error[rownames(dat$Y)], shared_error[rownames(dat$Y)])
  rownames(expected) <- rownames(dat$Y)
  colnames(expected) <- colnames(dat$Y)

  expect_equal(normalized, expected)
})

test_that("normalize_input_error only allows missing errors for missing observations", {
  dat <- small_lizard_data(n_tips = 6, traits = 1:2, normalize = FALSE)
  Y <- dat$Y
  Y[1, 1] <- NA
  Y[2, 2] <- NA

  input_error <- matrix(
    0.01,
    nrow = nrow(Y),
    ncol = ncol(Y),
    dimnames = dimnames(Y)
  )
  input_error[1, 1] <- NA
  input_error[2, 2] <- NA

  normalized <- kfl1ou:::normalize_input_error(dat$tree, Y, input_error)
  expect_true(is.na(normalized[1, 1]))
  expect_true(is.na(normalized[2, 2]))

  bad_input_error <- input_error
  bad_input_error[3, 1] <- NA
  expect_error(
    kfl1ou:::normalize_input_error(dat$tree, Y, bad_input_error),
    "input_error has missing values for observed traits"
  )
})

test_that("get_trait_input_error matches a reduced trait tree after missing-data filtering", {
  dat <- small_lizard_data(n_tips = 6, traits = 1:2, normalize = FALSE)
  input_error <- matrix(
    seq_len(length(dat$Y)),
    nrow = nrow(dat$Y),
    ncol = ncol(dat$Y),
    dimnames = list(rev(rownames(dat$Y)), rev(colnames(dat$Y)))
  )
  normalized <- kfl1ou:::normalize_input_error(dat$tree, dat$Y, input_error)

  available <- rep(TRUE, nrow(dat$Y))
  available[c(2, 5)] <- FALSE
  trait_tree <- ape::drop.tip(dat$tree, rownames(dat$Y)[!available])
  trait_tree <- ape::reorder.phylo(trait_tree, "postorder")

  trait_error <- kfl1ou:::get_trait_input_error(
    list(input_error = normalized),
    idx = 2,
    tree = trait_tree,
    available = available
  )

  expected <- normalized[available, 2]
  names(expected) <- rownames(normalized)[available]
  expected <- expected[trait_tree$tip.label]

  expect_equal(trait_error, expected)
  expect_equal(names(trait_error), trait_tree$tip.label)
})

test_that("prepare_trait_phylolm_data drops missing tips and aligns named input_error", {
  dat <- small_lizard_data(n_tips = 6, normalize = FALSE)
  Y <- dat$Y
  Y[c(2, 4), 1] <- NA
  input_error <- stats::setNames(seq(0.01, 0.06, length.out = nrow(Y)), rownames(Y))

  prepared <- kfl1ou:::prepare_trait_phylolm_data(dat$tree, Y, input_error = input_error)

  expect_equal(rownames(prepared$Y), prepared$tree$tip.label)
  expect_false(any(is.na(prepared$Y[, 1])))
  expect_equal(names(prepared$input_error), prepared$tree$tip.label)
  expect_equal(prepared$input_error, input_error[prepared$tree$tip.label])

  expect_error(
    kfl1ou:::prepare_trait_phylolm_data(dat$tree, Y, input_error = unname(input_error)),
    "input_error must have names"
  )
})

test_that("sqrt_OU_covariance validates observation-error arguments", {
  dat <- small_lizard_data(n_tips = 6, normalize = FALSE)

  expect_error(
    sqrt_OU_covariance(
      dat$tree,
      sigma2_error = -0.1,
      check.order = FALSE,
      check.ultrametric = FALSE
    ),
    "sigma2_error must be non-negative"
  )
  expect_error(
    sqrt_OU_covariance(
      dat$tree,
      sigma2 = 0,
      sigma2_error = 0.1,
      check.order = FALSE,
      check.ultrametric = FALSE
    ),
    "sigma2 must be strictly positive"
  )
  expect_error(
    sqrt_OU_covariance(
      dat$tree,
      sigma2 = 1,
      input_error = stats::setNames(rep(0.01, nrow(dat$Y)), paste0("x", seq_len(nrow(dat$Y)))),
      check.order = FALSE,
      check.ultrametric = FALSE
    ),
    "input_error names do not match the tree tip labels"
  )
  expect_error(
    sqrt_OU_covariance(
      dat$tree,
      sigma2 = 1,
      input_error = stats::setNames(c(rep(0.01, nrow(dat$Y) - 1), NA_real_), dat$tree$tip.label),
      check.order = FALSE,
      check.ultrametric = FALSE
    ),
    "input_error cannot contain missing values"
  )
})
