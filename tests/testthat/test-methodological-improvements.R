method_small_data <- function(n=8L, p=1L, seed=701L){
  set.seed(seed)
  tree <- reorder(ape::rcoal(n), "postorder")
  Y <- matrix(
    rnorm(n * p), n, p,
    dimnames=list(tree$tip.label, paste0("t", seq_len(p)))
  )
  list(tree=tree, Y=Y)
}

test_that("automatic small-tree search certifies the exhaustive optimum", {
  dat <- method_small_data()
  fit <- estimate_shift_configuration(
    dat$tree, dat$Y, max.nShifts=2L, criterion="BIC",
    alpha.lower=1, alpha.upper=1, quietly=TRUE
  )

  expect_identical(fit$search.diagnostics$strategy, "exhaustive")
  expect_true(fit$search.diagnostics$globally.optimal)
  expect_equal(fit$search.diagnostics$configuration.space.size, 92)
  expect_equal(fit$search.diagnostics$coverage, 1)
  expect_equal(length(fit$profile$configurations), 92L)
})

test_that("forced ensemble search preserves RNG and exposes partial coverage", {
  dat <- method_small_data(n=10L, seed=702L)
  set.seed(99L)
  before <- .Random.seed
  fit <- estimate_shift_configuration(
    dat$tree, dat$Y, max.nShifts=2L, criterion="BIC",
    search.strategy="ensemble", ensemble.replicates=3L,
    ensemble.seed=12L, quietly=TRUE
  )

  expect_identical(.Random.seed, before)
  expect_identical(fit$search.diagnostics$strategy, "ensemble")
  expect_false(fit$search.diagnostics$globally.optimal)
  expect_true(fit$search.diagnostics$coverage < 1)
})

test_that("full-selection intervals work from fit_OU and include nonselection", {
  dat <- method_small_data(n=10L, seed=703L)
  Z <- kfl1ou:::generate_design_matrix(dat$tree, "simpX")
  edge <- which(colSums(Z) >= 2L & colSums(Z) <= 6L)[[1L]]
  fit <- fit_OU(
    dat$tree, dat$Y, edge, criterion="BIC",
    alpha.lower=1, alpha.upper=1, search.max.nShifts=1L
  )
  parameter <- grep("^shift\\[", names(coef(fit)), value=TRUE)[[1L]]
  interval <- confint(
    fit, parm=parameter, selection="full", nsim=4L, seed=704L
  )

  expect_equal(attr(interval, "successful"), 4L)
  expect_equal(unname(attr(interval, "effective")[[parameter]]), 4L)
  expect_equal(attr(interval, "nonselection.value"), 0)
  expect_true(attr(interval, "reselection.probability")[[parameter]] <= 1)
  expect_true(is.finite(attr(interval, "partition.reselection.probability")))
})

test_that("root-complement shifts collapse to the same tip partition", {
  dat <- method_small_data(seed=705L)
  root <- length(dat$tree$tip.label) + 1L
  root.edges <- which(dat$tree$edge[, 1L] == root)
  expect_length(root.edges, 2L)

  expect_identical(
    shift_partition_key(dat$tree, root.edges[[1L]]),
    shift_partition_key(dat$tree, root.edges[[2L]])
  )
  equivalent <- equivalent_shift_configurations(
    dat$tree, root.edges[[1L]], max.nShifts=1L,
    candid.edges=root.edges,
    max.configurations=100L
  )
  expect_true(all(root.edges %in% unlist(equivalent)))

  summary <- summarize_shift_uncertainty(
    list(all.shifts=as.list(root.edges)), dat$tree
  )
  expect_equal(nrow(summary$configurations), 2L)
  expect_equal(nrow(summary$partitions), 1L)
  expect_equal(summary$partitions$probability, 1)
  expect_equal(unname(diag(summary$tip.coassignment)), rep(1, 8L))
})

test_that("multivariate support can retain trait-specific branch effects", {
  coefficients <- matrix(
    c(1, 0, 0, 0, 0, 0), nrow=6L, ncol=1L
  )
  groups <- rep(1:2, each=3L)
  any.support <- kfl1ou:::grplasso_support_summary(
    coefficients, groups, nVariables=3L, support.threshold=1L
  )
  majority.support <- kfl1ou:::grplasso_support_summary(
    coefficients, groups, nVariables=3L, support.threshold="majority"
  )
  expect_equal(any.support$df.vec, 1L)
  expect_equal(majority.support$df.vec, 0L)
})

test_that("full covariance pBIC and analytic shrinkage are finite", {
  dat <- method_small_data(n=10L, p=2L, seed=706L)
  fit <- fit_OU(
    dat$tree, dat$Y, integer(), criterion="pBIC",
    trait.covariance="full", covariance.regularization="shrinkage",
    alpha.lower=.7, alpha.upper=.7, optimizer.starts=1L,
    compute.hessian=FALSE
  )
  expect_true(is.finite(fit$score))
  expect_true(fit$regularization.lambda > 0)
  expect_lte(fit$regularization.lambda, .95)
  expect_false(fit$information.criterion.calibrated)
})

test_that("predictive diagnostics use projected residual contrasts", {
  dat <- method_small_data(n=10L, seed=707L)
  fit <- fit_OU(
    dat$tree, dat$Y, integer(), criterion="BIC",
    alpha.lower=.8, alpha.upper=.8
  )
  diagnostics <- diagnose_l1ou(
    fit, nsim=4L, seed=708L, min.clade.size=2L
  )
  expect_length(diagnostics$standardized.residuals, 9L)
  expect_equal(diagnostics$predictive.simulations, 4L)
  expect_true(all(is.finite(diagnostics$predictive.p.value)))
  rate <- check_rate_heterogeneity(
    fit, nsim=3L, seed=709L, min.clade.size=2L
  )
  expect_s3_class(rate, "l1ou_rate_heterogeneity")
  expect_true(is.finite(rate$p.value))
})

test_that("branchwise covariance and rate-shift comparison are coherent", {
  dat <- method_small_data(n=10L, seed=710L)
  alpha <- .9
  expect_equal(
    kfl1ou:::branchwise_ou_covariance(dat$tree, alpha),
    kfl1ou:::multivariate_ou_base_covariance(
      dat$tree, alpha, "OUfixedRoot"
    ),
    tolerance=1e-10
  )
  fit <- fit_OU(
    dat$tree, dat$Y, integer(), criterion="BIC",
    alpha.lower=alpha, alpha.upper=alpha
  )
  comparison <- compare_diffusion_rate_shift(
    fit, min.clade.size=2L, nsim=2L, seed=711L
  )
  expect_s3_class(comparison, "l1ou_rate_shift_comparison")
  expect_true(is.finite(comparison$statistic))
  expect_equal(comparison$successful, 2L)
})

test_that("replicates and tree ensembles propagate structural uncertainty", {
  dat <- method_small_data(n=8L, p=2L, seed=712L)
  replicate.Y <- rbind(
    dat$Y,
    dat$Y + matrix(rnorm(length(dat$Y), sd=.1), nrow=8L)
  )
  replicated <- fit_l1ou_replicates(
    dat$tree, replicate.Y, rep(dat$tree$tip.label, 2L),
    shift.configuration=integer(), criterion="BIC",
    alpha.lower=.7, alpha.upper=.7
  )
  expect_equal(unname(replicated$replicate.summary$counts), matrix(2, 8, 2))
  expect_true(all(replicated$replicate.summary$input.error >= 0))

  ensemble <- fit_l1ou_tree_ensemble(
    list(dat$tree, dat$tree), dat$Y,
    max.nShifts=0L, criterion="BIC"
  )
  expect_s3_class(ensemble, "l1ou_tree_ensemble")
  expect_equal(ensemble$successful, 2L)
  expect_equal(ensemble$pairwise.ari, matrix(1, 2, 2))
})

test_that("measurement error and model assumptions have explicit sensitivity checks", {
  dat <- method_small_data(n=8L, seed=713L)
  fitted.error <- fit_OU(
    dat$tree, dat$Y, integer(), criterion="BIC",
    measurement_error=TRUE, alpha.lower=.8, alpha.upper=.8,
    search.max.nShifts=0L
  )
  error.profile <- profile_measurement_error_l1ou(
    fitted.error, multipliers=c(0, 1)
  )
  expect_equal(error.profile$multiplier, c(0, 1))
  expect_true(all(is.finite(error.profile$logLik)))

  sensitivity <- sensitivity_l1ou(
    fitted.error, root.models="OUfixedRoot", alpha.upper=c(.8, 1.2),
    criteria="BIC", selection="full", keep.fits=TRUE
  )
  expect_s3_class(sensitivity, "l1ou_sensitivity")
  expect_equal(nrow(sensitivity), 2L)
  expect_true(all(is.finite(sensitivity$score)))
  expect_length(attr(sensitivity, "fits"), 2L)
})

test_that("covariance comparison can repeat shift discovery in each bootstrap", {
  dat <- method_small_data(n=8L, p=2L, seed=714L)
  comparison <- compare_trait_covariance(
    dat$tree, dat$Y, nboot=1L, seed=715L,
    optimizer.starts=1L, selection="full", search.max.nShifts=0L
  )
  expect_s3_class(comparison, "l1ou_covariance_comparison")
  expect_identical(comparison$selection, "full")
  expect_equal(comparison$attempted, 1L)
  expect_equal(comparison$successful, 1L)
})
