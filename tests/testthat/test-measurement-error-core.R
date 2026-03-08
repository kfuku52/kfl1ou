test_that("sqrt_OU_covariance adds scalar and tip-specific observation error under BM", {
  dat <- small_lizard_data(n_tips = 8, normalize = FALSE)
  input_error <- stats::setNames(seq(0.01, 0.08, length.out = nrow(dat$Y)), rownames(dat$Y))

  res <- sqrt_OU_covariance(
    dat$tree,
    alpha = 0,
    sigma2 = 2,
    sigma2_error = 0.4,
    input_error = input_error,
    check.order = FALSE,
    check.ultrametric = FALSE
  )

  expected <- ape::vcv(dat$tree)
  diag(expected) <- diag(expected) + 0.4 / 2 + input_error[dat$tree$tip.label] / 2

  expect_equal(
    unname(res$sqrtSigma %*% t(res$sqrtSigma)),
    unname(expected),
    tolerance = 1e-8
  )
  expect_equal(
    unname(res$sqrtInvSigma %*% t(res$sqrtInvSigma)),
    unname(solve(expected)),
    tolerance = 1e-8
  )
})

test_that("fit_OU with measurement error matches direct phylolm for a no-shift model", {
  testthat::skip_if_not_installed("phylolm")
  dat <- small_lizard_data(n_tips = 20)
  fit <- fit_OU(dat$tree, dat$Y, shift.configuration = c(), measurement_error = TRUE)

  direct_data <- data.frame(
    trait = dat$Y[, 1],
    int = 1,
    row.names = rownames(dat$Y)
  )
  direct <- suppressWarnings(
    phylolm::phylolm(
      trait ~ int - 1,
      data = direct_data,
      phy = dat$tree,
      model = "OUfixedRoot",
      measurement_error = TRUE,
      upper.bound = kfl1ou:::alpha_upper_bound(dat$tree)
    )
  )

  expect_equal(unname(fit$alpha), unname(direct$optpar), tolerance = 1e-10)
  expect_equal(unname(fit$sigma2), unname(direct$sigma2), tolerance = 1e-10)
  expect_equal(unname(fit$sigma2_error), unname(direct$sigma2_error), tolerance = 1e-10)
  expect_equal(unname(fit$logLik), unname(direct$logLik), tolerance = 1e-10)
  expect_equal(unname(fit$intercept), unname(direct$coefficients[[1]]), tolerance = 1e-10)
})

test_that("configuration_ic reuses the measurement-error fit consistently", {
  dat <- small_lizard_data(n_tips = 16)

  fit <- fit_OU(dat$tree, dat$Y, shift.configuration = c(), measurement_error = TRUE)
  numeric_score <- configuration_ic(
    dat$tree,
    dat$Y,
    shift.configuration = c(),
    measurement_error = TRUE
  )
  fit_from_ic <- configuration_ic(
    dat$tree,
    dat$Y,
    shift.configuration = c(),
    measurement_error = TRUE,
    fit.OU.model = TRUE
  )

  expect_equal(numeric_score, fit$score, tolerance = 1e-10)
  expect_equal(fit_from_ic$score, fit$score, tolerance = 1e-10)
  expect_equal(fit_from_ic$sigma2_error, fit$sigma2_error, tolerance = 1e-10)
})

test_that("print output reports sigma2_error when measurement error is estimated", {
  dat <- small_lizard_data(n_tips = 12)
  fit <- fit_OU(dat$tree, dat$Y, shift.configuration = c(), measurement_error = TRUE)

  printed <- capture.output(print(fit))
  summarized <- capture.output(summary(fit))

  expect_true(any(grepl("sigma2_error", printed, fixed = TRUE)))
  expect_true(any(grepl("sigma2_error", summarized, fixed = TRUE)))
})
