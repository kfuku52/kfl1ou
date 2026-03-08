test_that("configuration_ic and fit_OU sanitize inverted alpha bounds", {
  dat <- small_lizard_data(n_tips = 10)

  expect_warning(
    score <- configuration_ic(
      dat$tree,
      dat$Y,
      shift.configuration = integer(),
      alpha.lower = 2,
      alpha.upper = 1
    ),
    "alpha.upper must be equal or greater than alpha.lower"
  )
  expect_true(is.finite(score))

  expect_warning(
    fit <- fit_OU(
      dat$tree,
      dat$Y,
      shift.configuration = integer(),
      alpha.lower = 2,
      alpha.upper = 1
    ),
    "alpha.upper must be equal or greater than alpha.lower"
  )
  expect_s3_class(fit, "l1ou")
  expect_equal(fit$alpha, 2)
})

test_that("backward convergent-regime search handles zero-shift models", {
  dat <- small_lizard_data(n_tips = 10)
  fit <- fit_OU(dat$tree, dat$Y, shift.configuration = integer())

  conv <- estimate_convergent_regimes(fit, criterion = "AICc", method = "backward")

  expect_s3_class(conv, "l1ou")
  expect_length(conv$shift.configuration, 0)
  expect_true(is.finite(conv$score))
})

test_that("rr convergent-regime search returns early for single-shift models", {
  dat <- small_lizard_data(n_tips = 12)
  fit <- fit_OU(dat$tree, dat$Y, shift.configuration = 1)

  conv <- estimate_convergent_regimes(fit, criterion = "AICc", method = "rr")

  expect_s3_class(conv, "l1ou")
  expect_equal(conv$shift.configuration, fit$shift.configuration)
  expect_true(is.finite(conv$score))
})

test_that("find_convergent_regimes tolerates rank-deficient design matrices", {
  dat <- small_lizard_data(n_tips = 16)

  expect_warning(
    out <- kfl1ou:::find_convergent_regimes(
      dat$tree,
      dat$Y,
      alpha = 1e-7,
      criterion = "AICc",
      regimes = as.list(c(4, 25))
    ),
    "ridge penalty"
  )

  expect_type(out$lambda, "double")
  expect_true(ncol(out$beta) >= 1)
  expect_equal(dim(out$M), c(1, 2))
})

test_that("estimate_shift_configuration honors quietly for LASSO progress", {
  dat <- small_lizard_data(n_tips = 12)

  out <- capture.output(
    est <- estimate_shift_configuration(
      dat$tree,
      dat$Y,
      max.nShifts = 1,
      quietly = TRUE
    )
  )

  expect_s3_class(est, "l1ou")
  expect_false(any(grepl("Starting first LASSO|Starting second LASSO", out)))
})

test_that("parallel bootstrap keeps zero-shift replicates", {
  dat <- small_lizard_data(n_tips = 10)
  fit <- fit_OU(dat$tree, dat$Y, shift.configuration = integer())

  local_mocked_bindings(
    estimate_shift_configuration = function(...) {
      list(shift.configuration = integer())
    },
    l1ou_mclapply = function(X, FUN, ..., mc.cores = 1L) {
      lapply(X, FUN, ...)
    },
    .package = "kfl1ou"
  )

  out <- l1ou_bootstrap_support(
    fit,
    nItrs = 3,
    multicore = TRUE,
    nCores = 2,
    quietly = TRUE
  )

  expect_length(out$all.shifts, 3)
  expect_true(all(vapply(out$all.shifts, length, integer(1)) == 0L))
  expect_true(all(is.finite(out$detection.rate)))
  expect_true(all(out$detection.rate == 0))
})

test_that("sequential bootstrap skips failed replicates and errors when all fail", {
  dat <- small_lizard_data(n_tips = 10)
  fit <- fit_OU(dat$tree, dat$Y, shift.configuration = integer())

  calls <- 0L
  local_mocked_bindings(
    estimate_shift_configuration = function(...) {
      calls <<- calls + 1L
      if(calls == 1L){
        stop("boom")
      }
      list(shift.configuration = integer())
    },
    .package = "kfl1ou"
  )

  out <- l1ou_bootstrap_support(
    fit,
    nItrs = 2,
    multicore = FALSE,
    quietly = TRUE
  )

  expect_length(out$all.shifts, 1)
  expect_length(out$all.shifts[[1]], 0)
  expect_true(all(out$detection.rate == 0))

  local_mocked_bindings(
    estimate_shift_configuration = function(...) {
      stop("boom")
    },
    .package = "kfl1ou"
  )

  expect_error(
    l1ou_bootstrap_support(
      fit,
      nItrs = 2,
      multicore = FALSE,
      quietly = TRUE
    ),
    "all bootstrap replicates failed"
  )
})

test_that("profile.l1ou keeps scores aligned with configurations", {
  model <- structure(
    list(
      profile = list(
        scores = c(5, 1, 2, 3),
        configurations = list(c(7, 8), integer(), 4, 9)
      )
    ),
    class = "l1ou"
  )

  out <- profile(model)

  expect_equal(out$nShifts, c(0L, 1L, 2L))
  expect_equal(out$scores, c(1, 2, 5))
  expect_equal(out$shift.configurations[[1]], integer())
  expect_equal(out$shift.configurations[[2]], 4)
  expect_equal(out$shift.configurations[[3]], c(7, 8))
})

test_that("profile.l1ou falls back to fitted model when profile is absent", {
  dat <- small_lizard_data(n_tips = 10)
  fit <- fit_OU(dat$tree, dat$Y, shift.configuration = integer())

  out <- profile(fit)

  expect_equal(out$nShifts, 0L)
  expect_equal(out$scores, fit$score)
  expect_length(out$shift.configurations, 1)
  expect_length(out$shift.configurations[[1]], 0)
})

test_that("fit_OU keeps trait-specific shift means in multivariate fits", {
  data(lizard.tree, lizard.traits)
  dat <- adjust_data(lizard.tree, lizard.traits[, 1:2], quietly = TRUE)
  fit <- fit_OU(dat$tree, dat$Y, shift.configuration = 1L)

  expected <- vapply(seq_along(fit$alpha), function(i) {
    edge.max <- apply(
      kfl1ou:::generate_design_matrix(dat$tree, "orgX", alpha = fit$alpha[[i]])[, 1, drop = FALSE],
      2,
      max
    )
    fit$shift.values[1, i] * edge.max
  }, numeric(1))

  expect_equal(drop(fit$shift.means), expected, tolerance = 1e-8)
})

test_that("multivariate second search only uses two grplasso passes", {
  dat <- small_lizard_data(n_tips = 12, traits = 1:2)
  old_run_grplasso <- getFromNamespace("run_grplasso", "kfl1ou")
  calls <- 0L

  local_mocked_bindings(
    run_grplasso = function(...) {
      calls <<- calls + 1L
      old_run_grplasso(...)
    },
    .package = "kfl1ou"
  )

  fit <- suppressWarnings(
    estimate_shift_configuration(
      dat$tree,
      dat$Y,
      max.nShifts = 2,
      quietly = TRUE
    )
  )

  expect_s3_class(fit, "l1ou")
  expect_equal(calls, 2L)
})

test_that("multivariate fit can parallelize over traits", {
  dat <- small_lizard_data(n_tips = 12, traits = 1:2)
  calls <- 0L

  local_mocked_bindings(
    l1ou_supports_multicore = function() TRUE,
    l1ou_mclapply = function(X, FUN, ..., mc.cores = 1L) {
      calls <<- calls + 1L
      lapply(X, FUN, ...)
    },
    .package = "kfl1ou"
  )

  opt <- list(
    criterion = "AICc",
    root.model = "OUfixedRoot",
    quietly = TRUE,
    alpha.starting.value = NA_real_,
    alpha.upper.bound = kfl1ou:::alpha_upper_bound(dat$tree),
    alpha.lower.bound = NA_real_,
    nCores = 2L,
    parallel.computing = FALSE,
    use.saved.scores = FALSE,
    multivariate.missing = FALSE,
    measurement_error = FALSE,
    input_error = NULL
  )

  fit <- fit_OU(dat$tree, dat$Y, shift.configuration = integer(), l1ou.options = opt)

  expect_s3_class(fit, "l1ou")
  expect_gte(calls, 1L)
})

test_that("parallel candidate search initializes a worker-local score cache once", {
  erase_calls <- 0L
  saved_score_flags <- logical()

  local_mocked_bindings(
    do_backward_correction = function(tree, Y, shift.configuration, opt) {
      saved_score_flags <<- c(saved_score_flags, isTRUE(opt$use.saved.scores))
      list(score = length(shift.configuration), shift.configuration = shift.configuration)
    },
    erase_configuration_score_db = function() {
      erase_calls <<- erase_calls + 1L
      invisible(NULL)
    },
    l1ou_mclapply = function(X, FUN, ..., mc.cores = 1L) {
      lapply(X, FUN, ...)
    },
    .package = "kfl1ou"
  )

  opt <- list(parallel.computing = TRUE, nCores = 2L, use.saved.scores = FALSE)
  out <- kfl1ou:::evaluate_candidate_configurations(
    tree = NULL,
    Y = NULL,
    candidate.configurations = list(1L, c(1L, 2L)),
    opt = opt
  )

  expect_equal(erase_calls, 1L)
  expect_true(all(saved_score_flags))
  expect_equal(out$shift.configuration, 1L)
})

test_that("candidate collection orders repeated shifts by numeric frequency", {
  dat <- small_lizard_data(n_tips = 12, traits = 1:2)
  n_edges <- Nedge(dat$tree)
  n_traits <- ncol(dat$Y)
  n_solutions <- 12L
  coeffs <- matrix(0, nrow = n_edges * n_traits, ncol = n_solutions)

  mark_edge <- function(col, edge) {
    coeffs[edge, col] <<- 1
    coeffs[edge + n_edges, col] <<- 1
  }

  mark_edges <- function(col, edges) {
    for (edge in edges) {
      mark_edge(col, edge)
    }
  }

  for (idx in 1:10) {
    mark_edges(idx, c(2L, idx + 2L))
  }
  mark_edge(11L, 1L)
  mark_edges(12L, c(1L, 2L))

  local_mocked_bindings(
    correct_unidentifiability = function(tree, shift.configuration, opt) shift.configuration,
    .package = "kfl1ou"
  )

  candidates <- kfl1ou:::collect_candidate_configurations(
    dat$tree,
    dat$Y,
    list(call = "grplasso", coefficients = coeffs),
    list(max.nShifts = 2L)
  )

  expect_equal(candidates[[length(candidates)]], c(2L, 1L))
})
