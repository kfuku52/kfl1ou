simulate_improvement_traits <- function(n.tips=20L, alpha=.7, seed=1L){
  set.seed(seed)
  tree <- reorder(ape::rcoal(n.tips), "postorder")
  Omega <- matrix(c(1, .5, .5, .8), 2)
  covariance <- kfl1ou:::multivariate_ou_dense_covariance(
    tree, alpha, Omega, "OUfixedRoot"
  )
  Y <- matrix(drop(t(chol(covariance)) %*% rnorm(2*n.tips)), n.tips, 2,
              dimnames=list(tree$tip.label, c("x", "y")))
  list(tree=tree, Y=Y)
}

test_that("pruning likelihood equals the dense observed Gaussian likelihood", {
  set.seed(201)
  tree <- reorder(ape::rcoal(13), "postorder")
  Omega <- matrix(c(1, .3, .3, .7), 2)
  error <- matrix(runif(26, 0, .02), 13, 2)
  residuals <- matrix(rnorm(26), 13, 2,
                      dimnames=list(tree$tip.label, c("x", "y")))
  residuals[c(2, 8), 1] <- NA
  residuals[4, 2] <- NA

  for(root in c("OUfixedRoot", "OUrandomRoot")){
    alpha <- c(.35, .9)
    opt <- list(root.model=root, input_error=error)
    covariance <- kfl1ou:::multivariate_ou_observed_covariance(
      tree, residuals, alpha, Omega, opt
    )
    y <- as.vector(residuals)
    y <- y[!is.na(y)]
    factor <- chol(covariance)
    dense <- -.5 * (length(y) * log(2*pi) +
      2 * sum(log(diag(factor))) +
      sum(forwardsolve(t(factor), y)^2))
    pruning <- kfl1ou:::pruning_multivariate_ou_loglik(
      tree, residuals, alpha, Omega, root, error
    )
    expect_equal(pruning, dense, tolerance=1e-8)
  }
})

test_that("full covariance supports trait-specific adaptation rates", {
  set.seed(202)
  tree <- reorder(ape::rcoal(24), "postorder")
  alpha <- c(.3, 1.1)
  Omega <- matrix(c(1, .45, .45, .8), 2)
  covariance <- kfl1ou:::multivariate_ou_dense_covariance(
    tree, alpha, Omega, "OUfixedRoot"
  )
  Y <- matrix(drop(t(chol(covariance)) %*% rnorm(48)), 24, 2,
              dimnames=list(tree$tip.label, c("slow", "fast")))
  fit <- fit_OU(
    tree, Y, integer(), criterion="BIC", trait.covariance="full",
    alpha.structure="diagonal", alpha.lower=alpha, alpha.upper=alpha,
    likelihood.engine="pruning", optimizer.starts=1,
    compute.hessian=FALSE
  )

  expect_equal(unname(fit$alpha), alpha, tolerance=1e-12)
  expect_equal(fit$likelihood.engine, "pruning")
  expect_true(all(eigen(fit$trait.covariance, symmetric=TRUE)$values > 0))
  expect_true(all(eigen(fit$tip.trait.covariance, symmetric=TRUE)$values > 0))
  expect_false(isTRUE(all.equal(
    fit$trait.correlation, fit$tip.trait.correlation, tolerance=1e-6
  )))
})

test_that("shrinkage permits high-dimensional trait covariance estimation", {
  set.seed(203)
  tree <- reorder(ape::rcoal(8), "postorder")
  p <- 8L
  Omega <- diag(p) + .15
  covariance <- kfl1ou:::multivariate_ou_dense_covariance(
    tree, .6, Omega, "OUfixedRoot"
  )
  Y <- matrix(drop(t(chol(covariance)) %*% rnorm(8*p)), 8, p,
              dimnames=list(tree$tip.label, paste0("t", seq_len(p))))

  expect_error(
    fit_OU(tree, Y, integer(), criterion="BIC", trait.covariance="full",
           alpha.lower=.6, alpha.upper=.6),
    "not estimable"
  )
  fit <- fit_OU(
    tree, Y, integer(), criterion="BIC", trait.covariance="full",
    covariance.regularization="shrinkage", alpha.lower=.6, alpha.upper=.6,
    compute.hessian=FALSE
  )
  expect_true(fit$regularization.lambda > 0)
  expect_true(all(eigen(fit$trait.covariance, symmetric=TRUE)$values > 0))
  fitted.covariance <- kfl1ou:::multivariate_ou_dense_covariance(
    tree, fit$alpha, fit$trait.covariance, "OUfixedRoot"
  )
  factor <- chol(fitted.covariance)
  residual <- as.vector(fit$residuals)
  direct.logLik <- -.5 * (length(residual) * log(2*pi) +
    2 * sum(log(diag(factor))) +
    sum(forwardsolve(t(factor), residual)^2))
  # The p = n shrinkage case is intentionally close to singular.  Different
  # BLAS/LAPACK implementations accumulate the equivalent matrix-normal and
  # dense Cholesky likelihoods in a different order, so compare at a tolerance
  # appropriate for this ill-conditioned validation case.
  expect_equal(fit$joint.logLik, direct.logLik, tolerance=1e-5)
})

test_that("standard model API, diagnostics and simulation are coherent", {
  dat <- simulate_improvement_traits(n.tips=18L, alpha=.7, seed=204L)
  fit <- fit_OU(
    dat$tree, dat$Y, integer(), criterion="BIC",
    trait.covariance="full", alpha.lower=.7, alpha.upper=.7
  )
  expect_equal(fitted(fit), fit$mu)
  expect_equal(residuals(fit), fit$residuals)
  expect_equal(nobs(fit), nrow(fit$Y))
  expect_equal(fit$observed.entries, length(fit$Y))
  expect_true(is.finite(as.numeric(logLik(fit))))
  expect_equal(evolutionary_vcov(fit), fit$trait.covariance)
  expect_warning(coefficient.vcov <- vcov(fit), "unavailable")
  expect_true(all(is.na(coefficient.vcov)))
  expect_length(coef(fit), 2L)
  simulated <- simulate(fit, nsim=2L, seed=1)
  expect_length(simulated, 2L)
  expect_equal(dim(simulated[[1L]]), dim(fit$Y))
  diagnostics <- diagnose_l1ou(fit)
  expect_length(diagnostics$standardized.residuals, length(fit$Y))

  interval <- confint(
    fit, parm="alpha:shared", method="parametric", nsim=2L, seed=2
  )
  expect_true(all(is.finite(interval)))
})

test_that("covariance comparison and shift uncertainty return calibrated structures", {
  dat <- simulate_improvement_traits(n.tips=16L, seed=205L)
  comparison <- compare_trait_covariance(
    dat$tree, dat$Y, nboot=0L, optimizer.starts=1L
  )
  expect_s3_class(comparison, "l1ou_covariance_comparison")
  expect_true(comparison$statistic >= 0)
  expect_true(is.na(comparison$p.value))

  bootstrap <- list(
    detection.rate=numeric(nrow(dat$tree$edge)),
    all.shifts=list(integer(), 1L, c(1L, 2L), 1L)
  )
  uncertainty <- summarize_shift_uncertainty(bootstrap, dat$tree)
  expect_equal(unname(uncertainty$edge.inclusion[[1L]]), .75)
  expect_equal(uncertainty$co.selection[1L, 2L], .25)
  expect_equal(sum(uncertainty$configurations$probability), 1)
})

test_that("model averaging propagates evaluated configuration uncertainty", {
  set.seed(206)
  tree <- reorder(ape::rcoal(12), "postorder")
  Y <- matrix(rnorm(12), 12, 1,
              dimnames=list(tree$tip.label, "trait"))
  model <- fit_OU(tree, Y, integer(), criterion="BIC")
  Z <- kfl1ou:::generate_design_matrix(tree, "simpX")
  candidate <- which(colSums(Z) >= 3 & colSums(Z) <= 8)[[1L]]
  model$profile <- list(
    scores=c(model$score, model$score + 2),
    configurations=list(integer(), candidate)
  )
  averaged <- model_average_l1ou(model, delta.max=5)
  expect_s3_class(averaged, "l1ou_model_average")
  expect_equal(sum(averaged$weights), 1)
  expect_equal(dim(averaged$mu), dim(Y))
  expect_true(averaged$edge.inclusion[[candidate]] > 0)
})
