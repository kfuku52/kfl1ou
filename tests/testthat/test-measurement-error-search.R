test_that("estimate_shift_configuration supports univariate measurement error", {
  dat <- small_lizard_data(n_tips = 12)

  est <- suppressWarnings(
    capture_silently(
      estimate_shift_configuration(
        dat$tree,
        dat$Y,
        max.nShifts = 1,
        measurement_error = TRUE,
        quietly = TRUE
      )
    )
  )

  expect_s3_class(est, "l1ou")
  expect_length(est$sigma2_error, 1)
  expect_true(is.finite(est$sigma2_error))
  expect_true(est$sigma2_error >= 0)
  expect_true(est$nShifts <= 1)
})

test_that("estimate_shift_configuration supports multivariate measurement error", {
  dat <- small_lizard_data(n_tips = 12, traits = 1:2)

  est <- suppressWarnings(
    capture_silently(
      estimate_shift_configuration(
        dat$tree,
        dat$Y,
        max.nShifts = 1,
        measurement_error = TRUE,
        quietly = TRUE
      )
    )
  )

  expect_s3_class(est, "l1ou")
  expect_equal(ncol(est$Y), 2)
  expect_length(est$sigma2_error, 2)
  expect_true(all(is.finite(est$sigma2_error)))
  expect_true(all(est$sigma2_error >= 0))
})

test_that("estimate_shift_configuration supports multivariate missing data with measurement error", {
  dat <- small_lizard_data(n_tips = 12, traits = 1:2)
  dat$Y[1, 1] <- NA
  dat$Y[2, 2] <- NA

  est <- suppressWarnings(
    capture_silently(
      estimate_shift_configuration(
        dat$tree,
        dat$Y,
        max.nShifts = 1,
        measurement_error = TRUE,
        quietly = TRUE
      )
    )
  )

  expect_s3_class(est, "l1ou")
  expect_true(isTRUE(est$l1ou.options$multivariate.missing))
  expect_length(est$l1ou.options$tree.list, 2)
  expect_equal(length(est$l1ou.options$tree.list[[1]]$tip.label), nrow(dat$Y) - 1)
  expect_equal(length(est$l1ou.options$tree.list[[2]]$tip.label), nrow(dat$Y) - 1)
  expect_length(est$sigma2_error, 2)
  expect_true(all(is.finite(est$sigma2_error)))
})

test_that("bootstrap uses the measurement-error covariance path", {
  dat <- small_lizard_data(n_tips = 10)
  est <- capture_silently(
    estimate_shift_configuration(
      dat$tree,
      dat$Y,
      max.nShifts = 1,
      measurement_error = TRUE,
      quietly = TRUE
    )
  )

  bs <- capture_silently(
    l1ou_bootstrap_support(
      est,
      nItrs = 1,
      multicore = FALSE,
      quietly = TRUE
    )
  )

  expect_equal(length(bs$detection.rate), nrow(dat$tree$edge))
  expect_equal(length(bs$all.shifts), 1)
  expect_true(all(is.finite(bs$detection.rate)))
  expect_true(all(bs$detection.rate >= 0 & bs$detection.rate <= 1))
})

test_that("convergent-regime fitting rejects measurement error explicitly", {
  dat <- small_lizard_data(n_tips = 10)

  expect_error(
    fit_OU(
      dat$tree,
      dat$Y,
      shift.configuration = c(1),
      cr.regimes = list(c(0), c(1)),
      measurement_error = TRUE
    ),
    "does not yet support measurement_error or input_error"
  )
})

test_that("fit_OU supports fixed-alpha input_error fits", {
  dat <- small_lizard_data(n_tips = 10)
  input_error <- named_input_error(dat$Y, value = 0.01)
  alpha <- 0.4

  fit <- fit_OU(
    dat$tree,
    dat$Y,
    shift.configuration = c(),
    input_error = input_error,
    alpha.lower = alpha,
    alpha.upper = alpha,
    alpha.starting.value = alpha
  )
  manual <- kfl1ou:::dense_known_input_error_gls_fit(
    dat$tree,
    dat$Y,
    preds = matrix(1, nrow = nrow(dat$Y), ncol = 1),
    model = "OUfixedRoot",
    lower.bound = alpha,
    upper.bound = alpha,
    starting.value = alpha,
    input_error = input_error,
    coefficient_names = "(Intercept)"
  )

  expect_s3_class(fit, "l1ou")
  expect_equal(unname(fit$alpha), alpha, tolerance = 1e-12)
  expect_equal(unname(fit$intercept), unname(manual$coefficients[[1]]), tolerance = 1e-10)
  expect_equal(unname(fit$sigma2), unname(manual$sigma2), tolerance = 1e-10)
  expect_equal(unname(fit$logLik), unname(manual$logLik), tolerance = 1e-10)
})

test_that("estimate_shift_configuration supports input_error", {
  dat <- small_lizard_data(n_tips = 10)
  input_error <- named_input_error(dat$Y, value = 0.01)

  est <- capture_silently(
    estimate_shift_configuration(
      dat$tree,
      dat$Y,
      max.nShifts = 0,
      input_error = input_error,
      quietly = TRUE
    )
  )

  expect_s3_class(est, "l1ou")
  expect_equal(est$nShifts, 0)
  expect_true(all(is.finite(est$sigma2)))
})

test_that("fit_OU supports input_error with measurement_error at fixed alpha", {
  dat <- small_lizard_data(n_tips = 10)
  input_error <- named_input_error(dat$Y, value = 0.01)
  alpha <- 0.4

  fit <- fit_OU(
    dat$tree,
    dat$Y,
    shift.configuration = c(),
    input_error = input_error,
    measurement_error = TRUE,
    alpha.lower = alpha,
    alpha.upper = alpha,
    alpha.starting.value = alpha
  )
  direct <- direct_joint_error_intercept_fit(
    dat$tree,
    dat$Y,
    alpha = alpha,
    input_error = input_error
  )

  expect_s3_class(fit, "l1ou")
  expect_equal(unname(fit$alpha), alpha, tolerance = 1e-12)
  expect_equal(unname(fit$intercept), unname(direct$coefficients), tolerance = 1e-8)
  expect_equal(unname(fit$sigma2), unname(direct$sigma2), tolerance = 1e-6)
  expect_equal(unname(fit$sigma2_error), unname(direct$sigma2_error), tolerance = 1e-8)
  expect_equal(unname(fit$logLik), unname(direct$logLik), tolerance = 1e-8)
})

test_that("estimate_shift_configuration supports input_error with measurement_error", {
  dat <- small_lizard_data(n_tips = 12, traits = 1:2)
  dat$Y[1, 1] <- NA
  dat$Y[2, 2] <- NA
  input_error <- matrix(
    0.01,
    nrow = nrow(dat$Y),
    ncol = ncol(dat$Y),
    dimnames = dimnames(dat$Y)
  )

  est <- capture_silently(
    estimate_shift_configuration(
      dat$tree,
      dat$Y,
      max.nShifts = 1,
      input_error = input_error,
      measurement_error = TRUE,
      quietly = TRUE
    )
  )

  expect_s3_class(est, "l1ou")
  expect_true(isTRUE(est$l1ou.options$multivariate.missing))
  expect_length(est$sigma2_error, 2)
  expect_true(all(is.finite(est$sigma2_error)))
  expect_true(all(est$sigma2_error >= 0))
})
