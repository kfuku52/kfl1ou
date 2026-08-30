review_data <- function(n=12L, p=1L, seed=3L){
  set.seed(seed)
  tree <- ape::reorder.phylo(ape::rcoal(n), "postorder")
  Y <- matrix(rnorm(n * p), n, p,
              dimnames=list(tree$tip.label, paste0("t", seq_len(p))))
  list(tree=tree, Y=Y)
}

review_convergent_data <- function(){
  tree <- ape::reorder.phylo(ape::read.tree(text=paste0(
    "(((a:1,b:1):1,(c:0.2,d:0.2):1.8):1,",
    "((e:0.1,f:0.1):1,(g:0.5,h:0.5):0.6):1.9);"
  )), "postorder")
  shifts <- which(tree$edge[, 2L] %in% match(c("a", "e"), tree$tip.label))
  regimes <- list(0L, shifts)
  W <- kfl1ou:::generate_prediction_vec(tree, shifts, regimes, .7)
  set.seed(1)
  Y <- matrix(W %*% c(0, 3) + rnorm(8, sd=.1), 8, 1,
              dimnames=list(tree$tip.label, "t1"))
  list(tree=tree, Y=Y, shifts=shifts, regimes=regimes, W=W)
}

test_that("ensemble paths actually run and retain their candidates", {
  dat <- review_data(30L)
  fit <- estimate_shift_configuration(
    dat$tree, dat$Y, max.nShifts=2L, criterion="BIC",
    search.strategy="ensemble", ensemble.replicates=4L
  )
  plain <- estimate_shift_configuration(
    dat$tree, dat$Y, max.nShifts=2L, criterion="BIC", search.strategy="lasso"
  )
  expect_equal(fit$search.diagnostics$ensemble$attempted, 8L)
  expect_equal(fit$search.diagnostics$ensemble$successful, 8L)
  expect_equal(fit$search.diagnostics$ensemble$failed, 0L)
  keys <- function(fit) vapply(fit$profile$configurations,
                              kfl1ou:::configuration_score_key, character(1))
  expect_true(all(keys(plain) %in% keys(fit)))
  expect_gt(length(keys(fit)), length(keys(plain)))
  expect_lte(fit$score, plain$score + 1e-8)
})

test_that("one-column paths retain dimensions and path failures are reported", {
  dat <- review_data(8L)
  run <- kfl1ou:::run_univariate_sparse_path
  local_mocked_bindings(run_univariate_sparse_path=function(XX, YY, opt){
    if(nrow(XX) < 7L) stop("deliberate subsample failure")
    run(XX, YY, opt)
  }, .package="kfl1ou")
  fit <- estimate_shift_configuration(
    dat$tree, dat$Y, candid.edges=1L, max.nShifts=1L,
    criterion="BIC", search.strategy="ensemble", ensemble.replicates=2L
  )
  expect_equal(fit$search.diagnostics$ensemble$attempted, 4L)
  expect_equal(fit$search.diagnostics$ensemble$failed, 4L)
  expect_equal(unname(fit$search.diagnostics$ensemble$failure.messages), 4L)
  expect_match(names(fit$search.diagnostics$ensemble$failure.messages),
               "deliberate subsample failure")
})

review_parallel_profiles <- function(){
  dat <- review_data(10L)
  args <- list(tree=dat$tree, Y=dat$Y, max.nShifts=2L, criterion="BIC",
               search.strategy="exhaustive", alpha.lower=.7, alpha.upper=.7)
  serial <- do.call(estimate_shift_configuration, c(args, list(nCores=1L)))
  parallel <- do.call(estimate_shift_configuration, c(args, list(nCores=2L)))
  expect_equal(parallel$profile, serial$profile, tolerance=1e-9)
  expect_gt(length(parallel$profile$configurations), 1L)
  serial.average <- model_average_l1ou(serial, delta.max=10)
  parallel.average <- model_average_l1ou(parallel, delta.max=10)
  expect_equal(parallel.average$configurations, serial.average$configurations)
  expect_equal(parallel.average$weights, serial.average$weights, tolerance=1e-9)
  expect_equal(parallel.average$mu, serial.average$mu, tolerance=1e-9)
}

test_that("fork search preserves complete score profiles and model weights", {
  skip_if_not(kfl1ou:::l1ou_supports_multicore())
  review_parallel_profiles()
})

test_that("socket search preserves complete score profiles and model weights", {
  skip_if_not(kfl1ou:::l1ou_supports_socket_cluster())
  # A source-loaded namespace cannot be imported by fresh R workers. This
  # contract runs against the built package in R CMD check on every platform.
  skip_if_not(file.exists(system.file("Meta", "package.rds", package="kfl1ou")))
  local_mocked_bindings(l1ou_supports_multicore=function() FALSE,
                         .package="kfl1ou")
  review_parallel_profiles()
})

test_that("default full-covariance search excludes inestimable configurations", {
  dat <- review_data(6L, 3L)
  fit <- estimate_shift_configuration(
    dat$tree, dat$Y, criterion="BIC", trait.covariance="full",
    alpha.lower=.7, alpha.upper=.7, compute.hessian=FALSE
  )
  expect_true(is.finite(fit$score))
  expect_equal(fit$search.diagnostics$estimability$requested.max.nShifts, 5L)
  expect_equal(fit$search.diagnostics$estimability$effective.max.nShifts, 2L)
  expect_match(fit$search.diagnostics$estimability$reason, "rank")
  expect_true(all(lengths(fit$profile$configurations) <= 2L))
  expect_true(fit$search.diagnostics$globally.optimal)
  invalid <- which(dat$tree$edge[, 2L] <= 6L)[1:3]
  expect_error(fit_OU(
    dat$tree, dat$Y, invalid, criterion="BIC", trait.covariance="full",
    alpha.lower=.7, alpha.upper=.7
  ), "not estimable")
})

test_that("simulation preserves positive covariance under unit changes", {
  dat <- review_data()
  fit <- fit_OU(dat$tree, dat$Y, integer(), criterion="BIC",
                 alpha.lower=.7, alpha.upper=.7)
  tiny <- fit
  scale <- 1e-9
  for(field in c("Y", "mu", "residuals", "optima", "intercept")){
    tiny[[field]] <- tiny[[field]] * scale
  }
  tiny$sigma2 <- tiny$sigma2 * scale^2
  for(engine in c("tree", "dense")){
    ordinary <- simulate(fit, nsim=2L, seed=11L, engine=engine)
    small <- simulate(tiny, nsim=2L, seed=11L, engine=engine)
    expect_equal(small[[1L]] / scale, ordinary[[1L]], tolerance=1e-12)
    expect_gt(max(abs(small[[1L]] - tiny$mu)), 0)
    expect_false(identical(small[[1L]], small[[2L]]))
  }
  singular <- matrix(1e-30, 2, 2)
  factor <- kfl1ou:::l1ou_gaussian_factor(singular)
  expect_equal(tcrossprod(factor) / 1e-30, matrix(1, 2, 2), tolerance=1e-12)
  expect_equal(kfl1ou:::draw_zero_mean_gaussian(matrix(0, 2, 2)), c(0, 0))
})

test_that("known errors are identical with or without a trivial convergence constraint", {
  for(p in 1:2){
    dat <- review_data(10L, p)
    error <- if(p == 1L) rep(.5, 10L) else matrix(.5, 10, p)
    for(estimated in c(FALSE, TRUE)){
      args <- list(tree=dat$tree, Y=dat$Y, shift.configuration=integer(),
                   criterion="BIC", alpha.lower=.7, alpha.upper=.7,
                   input_error=error, measurement_error=estimated)
      ordinary <- do.call(fit_OU, args)
      constrained <- do.call(fit_OU, c(args, list(cr.regimes=list(0L))))
      expect_equal(constrained$logLik, ordinary$logLik, tolerance=1e-6)
      expect_equal(constrained$mu, ordinary$mu, tolerance=1e-5)
      expect_equal(constrained$score, -2 * sum(constrained$logLik) +
                     constrained$information.parameter.count * log(10),
                   tolerance=1e-6)
    }
  }
})

test_that("conditional refits and rate comparisons keep convergent mean constraints", {
  dat <- review_convergent_data()
  fit <- fit_OU(dat$tree, dat$Y, dat$shifts, cr.regimes=dat$regimes,
                 criterion="BIC", alpha.lower=.7, alpha.upper=.7)
  sensitivity <- sensitivity_l1ou(
    fit, root.models="OUfixedRoot", alpha.upper=.7,
    criteria="BIC", keep.fits=TRUE
  )
  retained <- attr(sensitivity, "fits")[[1L]]
  expect_true(retained$convergent)
  expect_equal(retained$convergent.regimes, fit$convergent.regimes)
  expect_equal(as.numeric(logLik(retained)), as.numeric(logLik(fit)), tolerance=1e-8)
  design <- kfl1ou:::l1ou_diagnostic_mean_design(fit)
  expect_lt(norm(qr.resid(qr(design), dat$W), "F"), 1e-10)
  expect_lt(norm(qr.resid(qr(design), fit$mu), "F"), 1e-10)
  rate <- compare_diffusion_rate_shift(fit, min.clade.size=2L)
  expect_equal(rate$null.logLik, as.numeric(logLik(fit)), tolerance=1e-7)

  error.fit <- fit_OU(dat$tree, dat$Y, dat$shifts, cr.regimes=dat$regimes,
                       criterion="BIC", alpha.lower=.7, alpha.upper=.7,
                       measurement_error=TRUE)
  profile <- profile_measurement_error_l1ou(error.fit, multipliers=1)
  fixed <- fit_OU(dat$tree, dat$Y, dat$shifts, cr.regimes=dat$regimes,
                   criterion="BIC", alpha.lower=.7, alpha.upper=.7,
                   input_error=rep(error.fit$sigma2_error, nrow(dat$Y)))
  expect_equal(profile$logLik, as.numeric(logLik(fixed)), tolerance=1e-7)
})

test_that("full sensitivity inherits candidate-generation settings", {
  dat <- review_data(8L)
  fit <- estimate_shift_configuration(
    dat$tree, dat$Y, max.nShifts=1L, candid.edges=1:2, criterion="BIC",
    alpha.lower=.7, alpha.upper=.7, search.strategy="ensemble",
    ensemble.replicates=2L, ensemble.fraction=.6, ensemble.seed=42L,
    rescale=FALSE
  )
  refit <- attr(sensitivity_l1ou(
    fit, root.models="OUfixedRoot", alpha.upper=.7, criteria="BIC",
    selection="full", keep.fits=TRUE
  ), "fits")[[1L]]
  fields <- c("candid.edges", "search.strategy", "ensemble.replicates",
               "ensemble.fraction", "ensemble.seed", "rescale", "alpha.lower.bound")
  expect_equal(refit$l1ou.options[fields], fit$l1ou.options[fields])
  expect_equal(refit$profile, fit$profile, tolerance=1e-9)
})

test_that("convergent designs and known errors survive trait-specific missingness", {
  dat <- review_convergent_data()
  Y <- cbind(dat$Y, t2=dat$Y[, 1L] * .8 + seq_len(8L) / 10)
  Y[2L, 1L] <- NA
  Y[4L, 2L] <- NA
  fit <- fit_OU(dat$tree, Y, dat$shifts, cr.regimes=dat$regimes,
                 criterion="BIC", alpha.lower=.7, alpha.upper=.7,
                 input_error=matrix(.05, 8, 2))
  design <- kfl1ou:::l1ou_diagnostic_mean_design(fit)
  observed <- is.finite(as.vector(fit$mu))
  expect_lt(norm(qr.resid(qr(design[observed, , drop=FALSE]),
                          matrix(as.vector(fit$mu)[observed], ncol=1L)), "F"), 1e-9)
  retained <- kfl1ou:::refit_l1ou_model(fit)
  expect_equal(retained$mu, fit$mu, tolerance=1e-7)
  expect_equal(retained$score, fit$score, tolerance=1e-7)
  unconstrained.bounds <- sensitivity_l1ou(
    fit, root.models="OUfixedRoot", alpha.upper=.9, criteria="BIC", keep.fits=TRUE
  )
  expect_true(all(as.numeric(attr(unconstrained.bounds, "fits")[[1L]]$alpha) >= .7))
})

test_that("full selection repeats convergence and conditional refits keep a fixed alpha", {
  dat <- review_data(8L)
  fit <- estimate_shift_configuration(
    dat$tree, dat$Y, max.nShifts=0L, criterion="BIC", alpha.upper=1.5
  )
  conv <- estimate_convergent_regimes(fit, criterion="BIC", fixed.alpha=TRUE)
  fixed <- kfl1ou:::refit_l1ou_model(conv)
  expect_equal(fixed$alpha, conv$alpha, tolerance=1e-8)
  selected <- kfl1ou:::refit_l1ou_model(conv, selection="full")
  expect_true(selected$convergent)
  expect_equal(selected$convergent.regimes, list(0L))
  expect_equal(selected$l1ou.options$max.nShifts, 0L)
})

test_that("failed alpha grid rows retain their position and missing values", {
  dat <- review_data()
  fit <- fit_OU(dat$tree, dat$Y, integer(), root.model="OUrandomRoot",
                 alpha.lower=.7, alpha.upper=.7)
  for(grid in list(c(.7, 0), c(0, .7), c(.7, 0, 1), c(.7, 1, 0), c(0, 0))){
    profile <- profile_alpha_l1ou(fit, grid, keep.fits=TRUE)
    expect_equal(nrow(profile), length(grid))
    expect_length(attr(profile, "fits"), length(grid))
    expect_identical(is.na(profile$logLik), grid == 0)
    expect_identical(is.na(profile$score), grid == 0)
    expect_identical(nzchar(profile$error), grid == 0)
    expect_true(all(is.na(profile$delta.deviance[grid == 0])))
  }
})

test_that("reusable options never carry stale tree geometry", {
  dat <- review_data()
  fit <- fit_OU(dat$tree, dat$Y, 1L, criterion="BIC", alpha.lower=.7, alpha.upper=.7)
  expect_null(fit$l1ou.options$Z)
  expect_null(fit$l1ou.options$edge.age)
  legacy <- fit$l1ou.options
  legacy$Z <- kfl1ou:::generate_design_matrix(dat$tree, "simpX")
  legacy$edge.age <- kfl1ou:::tree_edge_ages(dat$tree)
  changed <- dat$tree
  changed$edge.length <- changed$edge.length * 10
  changed$tip.label <- rev(changed$tip.label)
  Y <- dat$Y[changed$tip.label, , drop=FALSE]
  reused <- fit_OU(changed, Y, 1L, l1ou.options=legacy)
  fresh <- fit_OU(changed, Y, 1L, criterion="BIC", alpha.lower=.7, alpha.upper=.7)
  expect_equal(reused$score, fresh$score, tolerance=1e-10)
  expect_equal(reused$shift.values, fresh$shift.values, tolerance=1e-10)
  expect_equal(reused$mu, fresh$mu, tolerance=1e-10)
})

test_that("parametric intervals name and align selected parameters", {
  dat <- review_data()
  fit <- fit_OU(dat$tree, dat$Y, integer(), criterion="BIC", alpha.lower=.7, alpha.upper=.7)
  all <- confint(fit, nsim=3L, seed=5L)
  parameters <- rev(names(kfl1ou:::l1ou_inference_parameters(fit)))[1:2]
  subset <- confint(fit, parm=parameters, nsim=3L, seed=5L)
  expect_identical(rownames(subset), parameters)
  expect_equal(as.numeric(subset), as.numeric(all[parameters, , drop=FALSE]))
})

test_that("dense diagnostic limits account for missing cells before allocation", {
  dat <- review_data(10L, 2L)
  dat$Y[cbind(1:10, rep(1:2, 5))] <- NA
  fit <- fit_OU(dat$tree, dat$Y, integer(), criterion="BIC", alpha.lower=.7, alpha.upper=.7)
  local_mocked_bindings(l1ou_simulation_covariance=function(...){
    stop("dense covariance must not be allocated")
  }, .package="kfl1ou")
  diagnostic <- diagnose_l1ou(fit, max.dense.dimension=12L, nsim=3L)
  expect_null(diagnostic$standardized.residuals)
  expect_match(diagnostic$whitening.note, "full response dimension 20")
})

test_that("fixed configurations avoid an all-edge design and search retains final Hessian", {
  dat <- review_data(8L, 2L)
  calls <- integer()
  original <- kfl1ou:::generate_design_matrix
  local_mocked_bindings(generate_design_matrix=function(...){
    value <- original(...)
    calls <<- c(calls, ncol(value))
    value
  }, .package="kfl1ou")
  fixed <- fit_OU(dat$tree, dat$Y, integer(), criterion="BIC", trait.covariance="full",
                   alpha.structure="diagonal", alpha.lower=.7, alpha.upper=.7,
                   likelihood.engine="pruning", optimizer.starts=1L,
                   compute.hessian=FALSE)
  expect_true(all(calls == 0L))
  with.hessian <- estimate_shift_configuration(
    dat$tree, dat$Y, max.nShifts=1L, criterion="BIC", trait.covariance="full",
    alpha.structure="diagonal", alpha.lower=.7, alpha.upper=.7,
    likelihood.engine="dense", optimizer.starts=1L
  )
  without <- estimate_shift_configuration(
    dat$tree, dat$Y, max.nShifts=1L, criterion="BIC", trait.covariance="full",
    alpha.structure="diagonal", alpha.lower=.7, alpha.upper=.7,
    likelihood.engine="dense", optimizer.starts=1L, compute.hessian=FALSE
  )
  expect_equal(with.hessian$profile, without$profile, tolerance=1e-9)
  expect_true(is.matrix(with.hessian$optimization$hessian))
  expect_null(without$optimization$hessian)
})
