test_that("OU half-life conversion helpers cover limits and invalid inputs", {
  expect_equal(alpha_to_half_life(c(0, log(2))), c(Inf, 1))
  expect_equal(half_life_to_alpha(c(1, Inf)), c(log(2), 0))

  expect_error(alpha_to_half_life(-1), "non-negative")
  expect_error(alpha_to_half_life(Inf), "finite")
  expect_error(half_life_to_alpha(0), "positive")
  expect_error(half_life_to_alpha(NA_real_), "positive")
})

test_that("configuration lookup and prediction methods cover public contracts", {
  profiled <- structure(
    list(
      profile = list(
        scores = c(1, 2),
        configurations = list(integer(), 3L)
      )
    ),
    class = "l1ou"
  )
  expect_identical(get_shift_configuration(profiled, 0L), integer())
  expect_identical(get_shift_configuration(profiled, 1L), 3L)
  expect_error(get_shift_configuration(profiled, 2L), "no configuration")
  expect_error(get_shift_configuration(list(), 0L), "inherit")

  response <- matrix(1:4, ncol = 1)
  optimum <- matrix(5:8, ncol = 1)
  fitted <- structure(list(mu = response, optima = optimum), class = "l1ou")
  expect_identical(predict(fitted), response)
  expect_identical(predict(fitted, type = "optimum"), optimum)
  expect_error(predict(fitted, type = "unknown"), "arg")
})

test_that("result print methods return their inputs invisibly", {
  cases <- list(
    covariance = structure(
      list(statistic = 2.5, p.value = 0.2, successful = 9L),
      class = "l1ou_covariance_comparison"
    ),
    diagnostics = structure(
      list(
        convergence = TRUE,
        boundary = FALSE,
        likelihood.engine = "dense",
        trait.covariance.condition = 1,
        measurement.error.identifiability.warning = TRUE,
        predictive.simulations = 2L,
        predictive.p.value = c(maximum = 0.5),
        rate.heterogeneity.warning = TRUE,
        whitening.note = "Whitening used observed contrasts."
      ),
      class = "l1ou_diagnostics"
    ),
    model_average = structure(
      list(
        models = list(list(), list()),
        failed = 0L,
        scores = c(1, 2),
        weights = c(0.7, 0.3),
        configurations = list(integer(), 2L)
      ),
      class = "l1ou_model_average"
    ),
    rate_heterogeneity = structure(
      list(statistic = 1.5, p.value = 0.4, nsim = 10L),
      class = "l1ou_rate_heterogeneity"
    ),
    rate_shift = structure(
      list(best.edge = 2L, multiplier = 1.4, statistic = 3, p.value = 0.1),
      class = "l1ou_rate_shift_comparison"
    ),
    shift_uncertainty = structure(
      list(
        successful = 3L,
        attempted = 4L,
        failed = 1L,
        configurations = data.frame(
          configuration = c("none", "2"),
          count = c(2L, 1L),
          probability = c(2 / 3, 1 / 3)
        ),
        partitions = NULL
      ),
      class = "l1ou_shift_uncertainty"
    ),
    tree_ensemble = structure(
      list(
        successful = 2L,
        failed = 0L,
        weights = c(0.6, 0.4),
        shift.counts = c(0L, 1L),
        partition.weights = c(root = 0.6, shifted = 0.4)
      ),
      class = "l1ou_tree_ensemble"
    )
  )

  for (result in cases) {
    output <- capture.output(returned <- print(result))
    expect_gt(length(output), 0L)
    expect_identical(returned, result)
  }
})

test_that("restored model print and plot methods delegate to l1ou methods", {
  dat <- small_lizard_data(n_tips = 8)
  fit <- suppressWarnings(
    capture_silently(
      fit_OU(
        dat$tree,
        dat$Y,
        shift.configuration = integer(),
        criterion = "BIC"
      )
    )
  )
  fit$restoration <- list(
    removed.tips = "omitted-tip",
    ambiguous.shifts = TRUE,
    representative = "tipward",
    shift.edge.paths = list()
  )
  class(fit) <- c("restored_l1ou", "l1ou")

  output <- capture.output(returned <- print(fit))
  expect_match(paste(output, collapse = "\n"), "restoration note")
  expect_identical(returned, fit)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(
    plot(
      fit,
      plot.bar = FALSE,
      edge.shift.ann = FALSE,
      edge.label.ann = FALSE
    )
  )
})
