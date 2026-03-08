test_that("backward convergent-regime search supports joint observation-error fits", {
  dat <- small_lizard_data(n_tips = 10)
  input_error <- named_input_error(dat$Y, value = 0.01)

  fit <- fit_OU(
    dat$tree,
    dat$Y,
    shift.configuration = integer(),
    measurement_error = TRUE,
    input_error = input_error
  )

  conv <- estimate_convergent_regimes(
    fit,
    criterion = "AICc",
    method = "backward"
  )

  expect_s3_class(conv, "l1ou")
  expect_length(conv$shift.configuration, 0)
  expect_true(is.finite(conv$score))
})

test_that("fit_OU with convergent regimes supports joint observation-error fits", {
  dat <- small_lizard_data(n_tips = 10)
  input_error <- named_input_error(dat$Y, value = 0.01)

  fit <- fit_OU(
    dat$tree,
    dat$Y,
    shift.configuration = integer(),
    cr.regimes = list(0L),
    measurement_error = TRUE,
    input_error = input_error
  )

  expect_s3_class(fit, "l1ou")
  expect_true(is.finite(fit$cr.score))
})

test_that("rr convergent-regime search still rejects observation-error fits", {
  dat <- small_lizard_data(n_tips = 10)
  input_error <- named_input_error(dat$Y, value = 0.01)

  fit <- fit_OU(
    dat$tree,
    dat$Y,
    shift.configuration = integer(),
    measurement_error = TRUE,
    input_error = input_error
  )

  expect_error(
    estimate_convergent_regimes(fit, criterion = "AICc", method = "rr"),
    "does not yet support measurement_error or input_error"
  )
})
