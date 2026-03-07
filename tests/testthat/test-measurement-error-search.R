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

test_that("input_error handling is version-gated against the installed phylolm", {
  dat <- small_lizard_data(n_tips = 10)
  input_error <- named_input_error(dat$Y, value = 0.01)

  if (kfl1ou:::phylolm_supports_input_error()) {
    fit <- fit_OU(dat$tree, dat$Y, shift.configuration = c(), input_error = input_error)
    est <- capture_silently(
      estimate_shift_configuration(
        dat$tree,
        dat$Y,
        max.nShifts = 0,
        input_error = input_error,
        quietly = TRUE
      )
    )

    expect_s3_class(fit, "l1ou")
    expect_s3_class(est, "l1ou")
    expect_true(all(is.finite(fit$sigma2)))
  } else {
    expect_error(
      fit_OU(dat$tree, dat$Y, shift.configuration = c(), input_error = input_error),
      "does not support input_error"
    )
    expect_error(
      estimate_shift_configuration(
        dat$tree,
        dat$Y,
        max.nShifts = 0,
        input_error = input_error,
        quietly = TRUE
      ),
      "does not support input_error"
    )
  }
})
