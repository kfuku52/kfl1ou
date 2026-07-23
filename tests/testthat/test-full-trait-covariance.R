simulate_correlated_ou_traits <- function(n.tips=40L, alpha=0.7,
                                          Omega=matrix(c(1, 0.6, 0.6, 0.8), 2L),
                                          means=NULL, seed=1L,
                                          root.model="OUfixedRoot"){
  set.seed(seed)
  tree <- reorder(ape::rcoal(n.tips), "postorder")
  covariance <- kfl1ou:::multivariate_ou_base_covariance(
    tree, alpha, root.model
  )
  if(is.null(means)){
    means <- matrix(0, nrow=n.tips, ncol=nrow(Omega))
  }
  noise <- drop(t(chol(kronecker(Omega, covariance))) %*%
    stats::rnorm(n.tips * nrow(Omega)))
  Y <- means + matrix(noise, nrow=n.tips, ncol=nrow(Omega))
  dimnames(Y) <- list(tree$tip.label, paste0("trait", seq_len(ncol(Y))))
  list(tree=tree, Y=Y, alpha=alpha, Omega=Omega)
}

test_that("full trait covariance recovers correlated OU innovations", {
  dat <- simulate_correlated_ou_traits(n.tips=80L, seed=42L)
  fit <- fit_OU(
    dat$tree,
    dat$Y,
    integer(),
    criterion="BIC",
    trait.covariance="full",
    alpha.lower=dat$alpha,
    alpha.upper=dat$alpha
  )

  expect_equal(unname(fit$alpha), rep(dat$alpha, 2L), tolerance=1e-12)
  expect_equal(
    fit$trait.correlation[1L, 2L],
    stats::cov2cor(dat$Omega)[1L, 2L],
    tolerance=0.15
  )
  expect_true(all(eigen(fit$trait.covariance, symmetric=TRUE)$values > 0))
  expect_equal(fit$residuals, fit$Y - fit$mu, tolerance=1e-12)
  expect_length(fit$joint.logLik, 1L)
})

test_that("full trait covariance uses the observed joint marginal with missing data", {
  dat <- simulate_correlated_ou_traits(n.tips=25L, seed=44L)
  input.error <- matrix(
    0.01,
    nrow=nrow(dat$Y), ncol=ncol(dat$Y),
    dimnames=dimnames(dat$Y)
  )
  dat$Y[c(2L, 5L), 1L] <- NA
  dat$Y[c(7L, 9L, 11L), 2L] <- NA
  fit <- fit_OU(
    dat$tree,
    dat$Y,
    integer(),
    criterion="BIC",
    trait.covariance="full",
    input_error=input.error,
    alpha.lower=dat$alpha,
    alpha.upper=dat$alpha
  )

  expect_true(is.finite(fit$joint.logLik))
  expect_equal(is.na(fit$residuals), is.na(dat$Y))
  expect_true(all(is.finite(fit$trait.covariance)))
  expect_true(abs(fit$trait.correlation[1L, 2L]) > 0.1)
})

test_that("full trait covariance supports estimated observation error", {
  dat <- simulate_correlated_ou_traits(n.tips=32L, seed=49L)
  set.seed(50L)
  error.variance <- c(0.03, 0.06)
  dat$Y <- dat$Y + matrix(
    stats::rnorm(length(dat$Y), sd=rep(sqrt(error.variance), each=nrow(dat$Y))),
    nrow=nrow(dat$Y), ncol=ncol(dat$Y)
  )
  fit <- fit_OU(
    dat$tree,
    dat$Y,
    integer(),
    criterion="BIC",
    trait.covariance="full",
    measurement_error=TRUE,
    alpha.lower=dat$alpha,
    alpha.upper=dat$alpha
  )

  expect_length(fit$sigma2_error, 2L)
  expect_true(all(is.finite(fit$sigma2_error)))
  expect_true(all(fit$sigma2_error > 0))
  expect_true(all(eigen(fit$trait.covariance, symmetric=TRUE)$values > 0))
})

test_that("full trait covariance optimizes alpha under a random-root model", {
  dat <- simulate_correlated_ou_traits(
    n.tips=50L, alpha=0.8, seed=51L, root.model="OUrandomRoot"
  )
  fit <- fit_OU(
    dat$tree,
    dat$Y,
    integer(),
    criterion="BIC",
    root.model="OUrandomRoot",
    trait.covariance="full",
    alpha.upper=5
  )

  expect_equal(fit$l1ou.options$root.model, "OUrandomRoot")
  expect_true(all(is.finite(fit$alpha)))
  expect_true(fit$alpha[[1L]] > 0 && fit$alpha[[1L]] <= 5)
  expect_true(fit$trait.correlation[1L, 2L] > 0)
})

test_that("full covariance is used during multivariate shift selection", {
  set.seed(45L)
  tree <- reorder(ape::rcoal(18L), "postorder")
  alpha <- 0.7
  Omega <- matrix(c(1, 0.45, 0.45, 0.8), 2L)
  covariance <- kfl1ou:::multivariate_ou_base_covariance(
    tree, alpha, "OUfixedRoot"
  )
  Z <- kfl1ou:::generate_design_matrix(tree, "simpX")
  shift <- which.min(abs(colSums(Z) - 6L))
  means <- cbind(2 * Z[, shift], -1.5 * Z[, shift])
  noise <- drop(t(chol(kronecker(Omega, covariance))) %*% stats::rnorm(36L))
  Y <- means + matrix(noise, nrow=18L, ncol=2L)
  dimnames(Y) <- list(tree$tip.label, c("x", "z"))

  fit <- estimate_shift_configuration(
    tree,
    Y,
    max.nShifts=1L,
    criterion="BIC",
    trait.covariance="full",
    alpha.upper=3,
    quietly=TRUE
  )

  expect_equal(unname(fit$shift.configuration), shift)
  expect_equal(fit$l1ou.options$trait.covariance, "full")
  expect_true(is.finite(fit$trait.correlation[1L, 2L]))
})

test_that("convergent models retain the full trait covariance likelihood", {
  dat <- simulate_correlated_ou_traits(n.tips=24L, seed=46L)
  Z <- kfl1ou:::generate_design_matrix(dat$tree, "simpX")
  candidates <- which(colSums(Z) >= 3L & colSums(Z) <= 8L)
  shifts <- candidates[c(1L, length(candidates))]
  regimes <- list(0L, shifts)
  fit <- fit_OU(
    dat$tree,
    dat$Y,
    shifts,
    cr.regimes=regimes,
    criterion="BIC",
    trait.covariance="full",
    alpha.lower=dat$alpha,
    alpha.upper=dat$alpha
  )
  states <- kfl1ou:::convergent_regime_states(dat$tree, shifts, regimes)

  expect_true(isTRUE(fit$convergent))
  expect_length(unique(fit$optima[states$tip.regime == 1L, 1L]), 1L)
  expect_true(is.finite(fit$joint.logLik))
  expect_equal(fit$score, fit$cr.score)
})

test_that("full covariance bootstrap preserves the joint mode", {
  dat <- simulate_correlated_ou_traits(n.tips=14L, seed=47L)
  fit <- estimate_shift_configuration(
    dat$tree,
    dat$Y,
    max.nShifts=0L,
    criterion="BIC",
    trait.covariance="full",
    quietly=TRUE
  )
  bootstrap <- l1ou_bootstrap_support(
    fit, nItrs=2L, multicore=FALSE, quietly=TRUE
  )

  expect_length(bootstrap$all.shifts, 2L)
  expect_equal(bootstrap$detection.rate, rep(0, nrow(dat$tree$edge)))
})

test_that("localization-aware pBIC is supported for full covariance", {
  dat <- simulate_correlated_ou_traits(n.tips=12L, seed=48L)
  fit <- fit_OU(
    dat$tree,
    dat$Y,
    integer(),
    trait.covariance="full",
    criterion="pBIC",
    alpha.lower=dat$alpha,
    alpha.upper=dat$alpha
  )
  expect_true(is.finite(fit$score))
  expect_match(fit$score.note, "shift-location")
})
