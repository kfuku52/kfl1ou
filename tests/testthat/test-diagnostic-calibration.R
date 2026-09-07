test_that("named single-column observation errors are shared across traits", {
  dat <- small_lizard_data(n_tips=8, traits=1:2)
  shared <- matrix(seq(.01, .08, length.out=8), ncol=1,
                   dimnames=list(rev(rownames(dat$Y)), "sampling_variance"))
  expected <- shared[rownames(dat$Y), rep(1L, 2L), drop=FALSE]
  colnames(expected) <- colnames(dat$Y)
  for(input in list(shared, as.data.frame(shared))){
    actual <- fit_OU(dat$tree, dat$Y, integer(), criterion="BIC",
                     alpha.lower=.8, alpha.upper=.8, input_error=input)
    reference <- fit_OU(dat$tree, dat$Y, integer(), criterion="BIC",
                        alpha.lower=.8, alpha.upper=.8, input_error=expected)
    expect_equal(actual$l1ou.options$input_error, expected)
    expect_equal(as.numeric(logLik(actual)), as.numeric(logLik(reference)))
  }
})

test_that("clade diagnostic simulations remove the fitted mean space", {
  dat <- small_lizard_data(n_tips=12, traits=1:2)
  dat$Y[c(2, 5), 2] <- NA_real_
  Z <- kfl1ou:::generate_design_matrix(dat$tree, "simpX")
  edge <- which(colSums(Z) >= 3 & colSums(Z) <= 8)[1]
  fit <- fit_OU(dat$tree, dat$Y, edge, criterion="BIC",
                alpha.lower=.8, alpha.upper=.8)
  observed <- is.finite(as.vector(fit$residuals))
  covariance <- kfl1ou:::l1ou_simulation_covariance(fit)[observed, observed]
  factor <- t(chol(covariance))
  design <- kfl1ou:::l1ou_diagnostic_mean_design(fit)[observed, , drop=FALSE]
  whitened.design <- forwardsolve(factor, design)
  precision <- solve(covariance)
  expected <- kfl1ou:::with_l1ou_seed(851L, replicate(4L, {
    z <- rnorm(sum(observed))
    # GLS residuals after fitting the same mean model to a simulated response.
    residual <- drop(factor %*% lm.fit(whitened.design, z)$residuals)
    gradient <- drop(precision %*% residual)
    observed.indices <- which(observed)
    tip.scores <- vapply(seq_len(nrow(dat$Y)), function(i){
      indices <- match(i + (0:1) * nrow(dat$Y), observed.indices, nomatch=0L)
      indices <- indices[indices > 0L]
      g <- gradient[indices]
      drop(crossprod(g, solve(precision[indices, indices, drop=FALSE], g))) /
        length(indices)
    }, numeric(1))
    kfl1ou:::l1ou_clade_variance_statistic(dat$tree, tip.scores, 2L)
  }))
  actual <- diagnose_l1ou(fit, nsim=4L, seed=851L, min.clade.size=2L)
  expect_equal(unname(actual$predictive.statistics[, "clade.variance.ratio"]),
               unname(expected), tolerance=1e-9)
})

test_that("measurement-error fractions use tip variance for marginal models", {
  dat <- small_lizard_data(n_tips=8, traits=1:2)
  for(root in c("OUfixedRoot", "OUrandomRoot")){
    fit <- fit_OU(dat$tree, dat$Y, integer(), criterion="BIC",
                  alpha.lower=2, alpha.upper=2, root.model=root)
    # Set a known positive shared error to exercise the diagnostic denominator.
    fit$sigma2_error <- c(.1, .2)
    actual <- diagnose_l1ou(fit)
    expected <- vapply(1:2, function(j){
      process <- kfl1ou:::multivariate_ou_base_covariance(
        dat$tree, fit$alpha[j], root
      ) * fit$sigma2[j]
      fit$sigma2_error[j] / (fit$sigma2_error[j] + process[1, 1])
    }, numeric(1))
    expect_equal(unname(actual$measurement.error.fraction), unname(expected))
  }
})
