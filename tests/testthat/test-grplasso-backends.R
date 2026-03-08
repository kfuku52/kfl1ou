make_group_lasso_fixture <- function(seed = 1L) {
  set.seed(seed)
  x <- matrix(rnorm(24 * 12), nrow = 24, ncol = 12)
  y <- rnorm(24)
  group <- rep(1:6, each = 2)
  lambda_max <- kfl1ou:::run_package_grplasso_lambdamax(x, y, group)
  lambda <- lambda_max * c(1, 0.75, 0.5)

  list(x = x, y = y, group = group, lambda = lambda)
}

test_that("package grplasso backend uses standardize FALSE semantics", {
  fx <- make_group_lasso_fixture()

  direct_lambda <- grplasso::lambdamax(
    fx$x,
    fx$y,
    model = grplasso::LinReg(),
    index = fx$group,
    standardize = FALSE,
    center = FALSE
  )
  direct_path <- grplasso::grplasso(
    fx$x,
    y = fx$y,
    standardize = FALSE,
    center = FALSE,
    lambda = fx$lambda,
    model = grplasso::LinReg(),
    index = fx$group,
    control = grplasso::grpl.control(tol = 1e-8, trace = 0, save.y = FALSE)
  )

  expect_equal(kfl1ou:::run_package_grplasso_lambdamax(fx$x, fx$y, fx$group), direct_lambda, tolerance = 1e-12)

  wrapped_path <- kfl1ou:::run_package_grplasso_path(
    fx$x,
    fx$y,
    fx$group,
    fx$lambda,
    tol = 1e-8
  )

  expect_equal(wrapped_path$coefficients, as.matrix(direct_path$coefficients), tolerance = 1e-10)
})

test_that("vendored grplasso backend matches package backend", {
  fx <- make_group_lasso_fixture()

  package_path <- kfl1ou:::run_package_grplasso_path(
    fx$x,
    fx$y,
    fx$group,
    fx$lambda,
    tol = 1e-8
  )
  vendored_path <- kfl1ou:::run_vendored_grplasso_path(
    fx$x,
    fx$y,
    fx$group,
    fx$lambda,
    tol = 1e-8
  )

  expect_equal(unname(vendored_path$coefficients), unname(package_path$coefficients), tolerance = 1e-7)
})

test_that("cpp grplasso backend matches vendored backend", {
  fx <- make_group_lasso_fixture()

  vendored_path <- kfl1ou:::run_vendored_grplasso_path(
    fx$x,
    fx$y,
    fx$group,
    fx$lambda,
    tol = 1e-8
  )
  cpp_path <- kfl1ou:::run_cpp_grplasso_path(
    fx$x,
    fx$y,
    fx$group,
    fx$lambda,
    tol = 1e-8
  )

  expect_lt(
    max(abs(unname(cpp_path$coefficients) - unname(vendored_path$coefficients))),
    1e-5
  )
})

test_that("multivariate estimate_shift_configuration agrees across grplasso backends", {
  dat <- small_lizard_data(n_tips = 12, traits = 1:2)

  fit_with_backend <- function(backend) {
    local_mocked_bindings(
      resolve_grplasso_backend = function(opt) backend,
      .package = "kfl1ou"
    )

    suppressWarnings(
      estimate_shift_configuration(
        dat$tree,
        dat$Y,
        max.nShifts = 2,
        quietly = TRUE
      )
    )
  }

  fit_package <- fit_with_backend("package")
  fit_vendor <- fit_with_backend("vendor-r")
  fit_cpp <- fit_with_backend("cpp")

  expect_equal(fit_vendor$shift.configuration, fit_package$shift.configuration)
  expect_equal(fit_cpp$shift.configuration, fit_package$shift.configuration)
  expect_equal(fit_vendor$score, fit_package$score, tolerance = 1e-8)
  expect_equal(fit_cpp$score, fit_package$score, tolerance = 1e-8)
})
