make_bootstrap_test_fit <- function(n.tips = 9L, n.traits = 1L,
                                    full = FALSE, seed = 601L) {
  set.seed(seed)
  tree <- reorder(ape::rcoal(n.tips), "postorder")
  Y <- matrix(
    rnorm(n.tips * n.traits),
    nrow = n.tips,
    ncol = n.traits,
    dimnames = list(tree$tip.label, paste0("t", seq_len(n.traits)))
  )
  fit <- fit_OU(
    tree,
    Y,
    integer(),
    criterion = "BIC",
    trait.covariance = if (full) "full" else "diagonal",
    alpha.lower = 0.6,
    alpha.upper = 0.6,
    optimizer.starts = 1L,
    compute.hessian = FALSE
  )
  list(tree = tree, fit = fit)
}

test_that("bootstrap result finalization records support and failures", {
  result <- kfl1ou:::finalize_bootstrap_results(
    list(1L, c(1L, 2L), NA),
    numeric(4L),
    attempted = 3L,
    failure.messages = c("", "", "optimizer failed")
  )

  expect_equal(result$detection.rate, c(1, 0.5, 0, 0))
  expect_equal(result$attempted, 3L)
  expect_equal(result$successful, 2L)
  expect_equal(result$failed, 1L)
  expect_equal(unname(result$failure.messages[["optimizer failed"]]), 1L)
  expect_error(
    kfl1ou:::finalize_bootstrap_results(list(NA), numeric(2L)),
    "all bootstrap replicates failed"
  )
})

test_that("univariate residual bootstrap covers sequential and parallel paths", {
  dat <- make_bootstrap_test_fit()
  set.seed(602L)
  before <- .Random.seed
  calls <- 0L

  local_mocked_bindings(
    estimate_shift_configuration = function(...) {
      calls <<- calls + 1L
      list(shift.configuration = if (calls == 1L) 1L else integer())
    },
    .package = "kfl1ou"
  )
  sequential <- l1ou_bootstrap_support(
    dat$fit,
    nItrs = 2L,
    multicore = FALSE,
    quietly = TRUE,
    type = "residual",
    seed = 77L
  )

  expect_identical(.Random.seed, before)
  expect_equal(sequential$detection.rate[[1L]], 0.5)
  expect_equal(sequential$attempted, 2L)
  expect_equal(sequential$successful, 2L)
  expect_identical(sequential$type, "residual")
  expect_identical(sequential$seed, 77L)

  local_mocked_bindings(
    estimate_shift_configuration = function(...) {
      list(shift.configuration = 2L)
    },
    l1ou_supports_multicore = function() TRUE,
    l1ou_mclapply = function(X, FUN, ..., mc.cores = 1L) {
      lapply(X, FUN, ...)
    },
    with_l1ou_thread_limit = function(n, expr) force(expr),
    .package = "kfl1ou"
  )
  parallel <- l1ou_bootstrap_support(
    dat$fit,
    nItrs = 2L,
    multicore = TRUE,
    nCores = 2L,
    quietly = TRUE,
    type = "residual",
    seed = 78L
  )

  expect_equal(parallel$detection.rate[[2L]], 1)
  expect_equal(parallel$successful, 2L)
})

test_that("diagonal multivariate residual bootstrap preserves missingness", {
  dat <- make_bootstrap_test_fit(n.traits = 2L, seed = 603L)
  dat$fit$Y[1L, 1L] <- NA_real_
  dat$fit$Y[2L, 2L] <- NA_real_
  generated <- list()

  local_mocked_bindings(
    estimate_shift_configuration = function(tree, Y, ...) {
      generated[[length(generated) + 1L]] <<- Y
      list(shift.configuration = integer())
    },
    .package = "kfl1ou"
  )
  sequential <- l1ou_bootstrap_support(
    dat$fit,
    nItrs = 2L,
    multicore = FALSE,
    quietly = TRUE,
    type = "residual",
    seed = 79L
  )

  expect_equal(sequential$successful, 2L)
  expect_length(generated, 2L)
  expect_true(all(vapply(generated, function(Y) {
    identical(is.na(Y), is.na(dat$fit$Y))
  }, logical(1))))

  complete <- make_bootstrap_test_fit(n.traits = 2L, seed = 604L)
  local_mocked_bindings(
    estimate_shift_configuration = function(...) {
      list(shift.configuration = 1L)
    },
    l1ou_supports_multicore = function() TRUE,
    l1ou_mclapply = function(X, FUN, ..., mc.cores = 1L) {
      lapply(X, FUN, ...)
    },
    with_l1ou_thread_limit = function(n, expr) force(expr),
    .package = "kfl1ou"
  )
  parallel <- l1ou_bootstrap_support(
    complete$fit,
    nItrs = 2L,
    multicore = TRUE,
    nCores = 2L,
    quietly = TRUE,
    type = "residual",
    seed = 80L
  )

  expect_equal(parallel$detection.rate[[1L]], 1)
  expect_equal(parallel$successful, 2L)
})

test_that("full covariance residual bootstrap covers both execution modes", {
  dat <- make_bootstrap_test_fit(
    n.tips = 10L,
    n.traits = 2L,
    full = TRUE,
    seed = 605L
  )

  local_mocked_bindings(
    estimate_shift_configuration = function(...) {
      list(shift.configuration = integer())
    },
    .package = "kfl1ou"
  )
  sequential <- l1ou_bootstrap_support(
    dat$fit,
    nItrs = 2L,
    multicore = FALSE,
    quietly = TRUE,
    type = "residual",
    seed = 81L
  )
  expect_equal(sequential$successful, 2L)
  expect_true(all(sequential$detection.rate == 0))

  local_mocked_bindings(
    estimate_shift_configuration = function(...) {
      list(shift.configuration = 1L)
    },
    l1ou_supports_multicore = function() TRUE,
    l1ou_mclapply = function(X, FUN, ..., mc.cores = 1L) {
      lapply(X, FUN, ...)
    },
    with_l1ou_thread_limit = function(n, expr) force(expr),
    .package = "kfl1ou"
  )
  parallel <- l1ou_bootstrap_support(
    dat$fit,
    nItrs = 2L,
    multicore = TRUE,
    nCores = 2L,
    quietly = TRUE,
    type = "residual",
    seed = 82L
  )

  expect_equal(parallel$detection.rate[[1L]], 1)
  expect_equal(parallel$successful, 2L)
})

test_that("bootstrap validates arguments and falls back from unsupported forks", {
  dat <- make_bootstrap_test_fit(seed = 606L)

  expect_error(l1ou_bootstrap_support(list()), "not of class")
  expect_error(l1ou_bootstrap_support(dat$fit, nItrs = 0L), "nItrs")
  expect_error(l1ou_bootstrap_support(dat$fit, nCores = 0L), "nCores")
  expect_error(l1ou_bootstrap_support(dat$fit, multicore = NA), "multicore")

  local_mocked_bindings(
    l1ou_supports_multicore = function() FALSE,
    estimate_shift_configuration = function(...) {
      list(shift.configuration = integer())
    },
    .package = "kfl1ou"
  )
  expect_warning(
    fallback <- l1ou_bootstrap_support(
      dat$fit,
      nItrs = 1L,
      multicore = TRUE,
      quietly = TRUE,
      seed = 83L
    ),
    "unavailable"
  )
  expect_equal(fallback$successful, 1L)
})
