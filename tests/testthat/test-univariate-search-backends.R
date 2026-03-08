test_that("univariate estimate_shift_configuration works with the internal lasso path", {
  dat <- small_lizard_data(n_tips = 12)

  fit <- suppressWarnings(
    estimate_shift_configuration(
      dat$tree,
      dat$Y,
      max.nShifts = 2,
      quietly = TRUE
    )
  )

  expect_s3_class(fit, "l1ou")
  expect_true(is.finite(fit$score))
  expect_true(fit$nShifts <= 2)
})

test_that("univariate estimate_shift_configuration works with the internal stepwise path", {
  dat <- small_lizard_data(n_tips = 12)

  fit <- suppressWarnings(
    estimate_shift_configuration(
      dat$tree,
      dat$Y,
      max.nShifts = 2,
      lars.alg = "stepwise",
      quietly = TRUE
    )
  )

  expect_s3_class(fit, "l1ou")
  expect_true(is.finite(fit$score))
  expect_true(fit$nShifts <= 2)
})
