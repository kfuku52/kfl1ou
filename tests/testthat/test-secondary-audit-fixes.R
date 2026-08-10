test_that("tree and shift validators reject malformed topology and indices", {
  bad <- structure(list(
    edge = matrix(c(3, 1, 2, 4, 4, 2), ncol = 2, byrow = TRUE),
    tip.label = c("a", "b"),
    Nnode = 2L,
    edge.length = rep(1, 3)
  ), class = "phylo")

  expect_error(
    kfl1ou:::validate_l1ou_tree(bad),
    "tips must be leaves|cycle|disconnected"
  )
  internal.leaf <- structure(list(
    edge = matrix(c(3, 1, 3, 2, 3, 4), ncol = 2, byrow = TRUE),
    tip.label = c("a", "b"),
    Nnode = 2L,
    edge.length = rep(1, 3)
  ), class = "phylo")
  expect_error(
    kfl1ou:::validate_l1ou_tree(internal.leaf),
    "internal tree node"
  )

  tree <- reorder(ape::read.tree(text = "((a:1,b:1):1,c:2);"), "postorder")
  oversized <- tree
  oversized$Nnode <- .Machine$integer.max
  expect_error(
    kfl1ou:::validate_l1ou_tree(oversized),
    "valid positive number"
  )
  oversized.edge <- tree
  oversized.edge$edge[1, 1] <- 1e100
  expect_error(
    kfl1ou:::validate_l1ou_tree(oversized.edge),
    "invalid or disconnected"
  )
  expect_equal(kfl1ou:::validate_shift_configuration(tree, integer()), integer())
  expect_error(
    kfl1ou:::validate_shift_configuration(tree, 1.5),
    "integer edge indices"
  )
  expect_error(
    kfl1ou:::validate_shift_configuration(tree, c(1, 1)),
    "unique"
  )
  expect_error(
    shift_tip_partition(tree, 1.5),
    "unique, finite integer edge indices"
  )
  expect_error(
    convert_shifts2regions(tree, c(1, 1), c(1, 2)),
    "unique, finite integer edge indices"
  )
  expect_error(
    equivalent_shift_configurations(
      tree, 1L, candid.edges = 1.9, max.configurations = 20L
    ),
    "candid.edges"
  )
  expect_error(
    equivalent_shift_configurations(
      tree, 1L, candid.edges = c(1L, NA), max.configurations = 20L
    ),
    "candid.edges"
  )
  expect_error(
    equivalent_shift_configurations(
      tree, 1L, candid.edges = NaN, max.configurations = 20L
    ),
    "candid.edges"
  )
  expect_error(
    equivalent_shift_configurations(
      tree, 1L, candid.edges = list(NA), max.configurations = 20L
    ),
    "candid.edges"
  )
})

test_that("missing-data simulation can generate complete latent responses", {
  dat <- small_lizard_data(n_tips = 10, traits = 1:2)
  dat$Y[1, 1] <- NA_real_
  dat$Y[2, 2] <- NA_real_
  input.error <- matrix(
    0.01, nrow(dat$Y), ncol(dat$Y), dimnames = dimnames(dat$Y)
  )
  input.error[is.na(dat$Y)] <- NA_real_

  fit <- fit_OU(
    dat$tree, dat$Y, integer(),
    input_error = input.error,
    alpha.lower = 0.5, alpha.upper = 0.5
  )

  for (engine in c("tree", "dense")) {
    complete <- simulate(
      fit, nsim = 1, seed = 1, preserve.missing = FALSE, engine = engine
    )[[1]]
    masked <- simulate(
      fit, nsim = 1, seed = 1, preserve.missing = TRUE, engine = engine
    )[[1]]

    expect_false(anyNA(complete))
    expect_identical(is.na(masked), is.na(dat$Y))
  }
})

test_that("fixed-configuration APIs reject unsupported missingness patterns", {
  dat <- small_lizard_data(n_tips = 8, traits = 1:2)
  univariate <- dat$Y[, 1, drop = FALSE]
  univariate[1, 1] <- NA_real_
  expect_error(
    fit_OU(dat$tree, univariate, integer()),
    "univariate response"
  )
  expect_error(
    configuration_ic(dat$tree, univariate, integer()),
    "univariate response"
  )

  multivariate <- dat$Y
  multivariate[1, ] <- NA_real_
  expect_error(
    fit_OU(dat$tree, multivariate, integer()),
    "no observed traits"
  )
  expect_error(
    configuration_ic(dat$tree, multivariate, integer()),
    "no observed traits"
  )
})

test_that("masked simulation does not require unidentifiable latent shift means", {
  dat <- small_lizard_data(n_tips = 12, traits = 1:2)
  design <- kfl1ou:::generate_design_matrix(dat$tree, "simpX")
  sizes <- colSums(design > 0)
  shift <- which(sizes >= 2L & sizes <= nrow(dat$Y) - 2L)[[1L]]
  affected <- design[, shift] > 0
  dat$Y[affected, 2L] <- NA_real_
  fit <- fit_OU(
    dat$tree, dat$Y, shift,
    alpha.lower = 0.5, alpha.upper = 0.5
  )

  expect_true(is.na(fit$shift.means[1L, 2L]))
  for (engine in c("tree", "dense")) {
    masked <- simulate(
      fit, nsim = 1, seed = 2, preserve.missing = TRUE, engine = engine
    )[[1L]]
    expect_identical(is.na(masked), is.na(dat$Y))
  }
  expect_error(
    simulate(fit, nsim = 1, seed = 2, preserve.missing = FALSE),
    "inconsistent|finite expected response"
  )
})

test_that("replicate pooling uses within-group degrees of freedom", {
  tree <- reorder(
    ape::read.tree(text = "((a:1,b:1):1,c:2);"), "postorder"
  )
  Y <- matrix(
    c(0, 2, 0, 2, 1), ncol = 1,
    dimnames = list(NULL, "trait")
  )
  species <- c("a", "a", "b", "b", "c")

  fit <- fit_l1ou_replicates(
    tree, Y, species,
    shift.configuration = integer(),
    alpha.lower = 0.5, alpha.upper = 0.5
  )

  expect_equal(
    unname(fit$replicate.summary$pooled.within.species.variance),
    2
  )
  expect_equal(fit$replicate.summary$input.error["c", "trait"], 2)

  Y[2, 1] <- Inf
  expect_error(
    fit_l1ou_replicates(
      tree, Y, species,
      shift.configuration = integer(),
      alpha.lower = 0.5, alpha.upper = 0.5
    ),
    "infinite or NaN"
  )
})

test_that("tree ensembles exclude zero-weight trees before fitting", {
  dat <- small_lizard_data(n_tips = 8)
  invalid <- dat$tree
  invalid$edge.length[[1]] <- -1

  expect_error(
    fit_l1ou_tree_ensemble(
      list(dat$tree, dat$tree), dat$Y,
      tree.weights = rep(.Machine$double.xmax, 2)
    ),
    "finite, non-negative"
  )

  expect_error(
    fit_l1ou_tree_ensemble(
      list(invalid, dat$tree), dat$Y,
      tree.weights = c(1, 0),
      max.nShifts = 0, search.strategy = "exhaustive"
    ),
    "positive-weight tree"
  )

  fit <- fit_l1ou_tree_ensemble(
    list(invalid, dat$tree), dat$Y,
    tree.weights = c(0, 1),
    max.nShifts = 0, search.strategy = "exhaustive"
  )
  expect_equal(fit$weights, 1)
  expect_equal(fit$excluded.zero.weight, 1L)
  expect_true(all(is.finite(fit$tip.coassignment)))
})

test_that("rate-shift candidates require valid non-reference clades", {
  dat <- small_lizard_data(n_tips = 8)
  fit <- fit_OU(
    dat$tree, dat$Y, integer(),
    alpha.lower = 0.5, alpha.upper = 0.5
  )

  for (edge in list(0, -1, 1.5, ape::Nedge(dat$tree))) {
    expect_error(
      compare_diffusion_rate_shift(
        fit, candidate.edges = edge, min.clade.size = 2, nsim = 0
      ),
      "candidate.edges"
    )
  }
})

test_that("shared-alpha covariance comparisons are labelled non-nested", {
  dat <- small_lizard_data(n_tips = 12, traits = 1:3)
  comparison <- compare_trait_covariance(
    dat$tree, dat$Y,
    nboot = 0, alpha.structure = "shared", optimizer.starts = 1
  )

  expected.df <- attr(logLik(comparison$alternative), "df") -
    attr(logLik(comparison$null), "df")
  expect_false(comparison$nested)
  expect_equal(comparison$df, expected.df)
  expect_false(comparison$calibrated)
  expect_match(comparison$calibration.note, "not nested")
})

test_that("regularized covariance comparisons print sensitivity status", {
  dat <- small_lizard_data(n_tips = 10, traits = 1:3)
  comparison <- compare_trait_covariance(
    dat$tree, dat$Y,
    nboot = 0,
    covariance.regularization = "shrinkage",
    regularization.lambda = 0.2,
    optimizer.starts = 1
  )
  printed <- capture.output(print(comparison))

  expect_identical(comparison$comparison, "regularized sensitivity")
  expect_false(comparison$calibrated)
  expect_match(paste(printed, collapse = "\n"), "regularized sensitivity")
  expect_match(paste(printed, collapse = "\n"), "covariance penalty")
})

test_that("diagnostics distinguish raw and trait-whitened correlations", {
  residual <- matrix(
    c(1, 1.1, 2, 2.2, 3, 2.8, 4, 4.3), ncol = 2,
    dimnames = list(NULL, c("a", "b"))
  )
  covariance <- matrix(c(2, 1.2, 1.2, 1), 2, 2)
  standardized <- kfl1ou:::l1ou_standardized_trait_correlation(
    residual, covariance
  )

  expect_equal(dim(standardized), c(2, 2))
  expect_true(all(is.finite(standardized)))
  expect_false(isTRUE(all.equal(
    standardized, stats::cor(residual), tolerance = 1e-10
  )))
})
