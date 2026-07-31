cr_audit_data <- function(n=14L, p=1L, seed=4101L){
  set.seed(seed)
  tree <- ape::reorder.phylo(ape::rcoal(n), "postorder")
  Z <- kfl1ou:::generate_design_matrix(tree, "simpX")
  sizes <- colSums(Z > 0)
  shift <- which(sizes >= 2L & sizes <= n - 2L)[[1L]]
  Y <- matrix(
    stats::rnorm(n * p), nrow=n, ncol=p,
    dimnames=list(tree$tip.label, paste0("t", seq_len(p)))
  )
  list(tree=tree, Y=Y, Z=Z, shift=shift)
}


test_that("nested convergent designs transfer weight only from the parent regime", {
  tree <- ape::reorder.phylo(
    ape::read.tree(text="((a:1,b:1):1,c:2);"),
    "postorder"
  )
  shifts <- c(3L, 1L)
  regimes <- list(0L, 3L, 1L)

  design <- kfl1ou:::generate_prediction_vec(
    tree, shifts, regimes, alpha=.7
  )
  reordered <- kfl1ou:::generate_prediction_vec(
    tree, rev(shifts), regimes, alpha=.7
  )

  expect_equal(rowSums(design), rep(1, 3), tolerance=1e-12)
  expect_equal(drop(design %*% rep(2, 3)), rep(2, 3), tolerance=1e-12)
  expect_identical(dimnames(reordered), dimnames(design))
  expect_equal(reordered, design, tolerance=1e-12)
})


test_that("whitened convergent designs whiten the background regime", {
  tree <- ape::reorder.phylo(
    ape::read.tree(text="(((a:1,b:1):1,c:2):1,d:3);"),
    "postorder"
  )
  shifts <- c(3L, 1L)
  alpha <- 0.7
  design <- kfl1ou:::generate_prediction_vec(
    tree,
    shifts,
    list(c(0L, shifts)),
    alpha=alpha,
    designMatrix=TRUE
  )
  inverse.factor <- t(sqrt_OU_covariance(
    tree, alpha=alpha, root.model="OUfixedRoot"
  )$sqrtInvSigma)

  expect_equal(
    rowSums(design),
    drop(inverse.factor %*% rep(1, length(tree$tip.label))),
    tolerance=1e-12
  )
})


test_that("phylolm_CR returns the design evaluated at its fitted alpha", {
  dat <- cr_audit_data(n=16L)
  terminal <- which(dat$tree$edge[, 2L] <= length(dat$tree$tip.label))
  shifts <- terminal[1:2]
  regimes <- list(0L, shifts)
  starting.alpha <- .2
  preds <- kfl1ou:::generate_prediction_vec(
    dat$tree, shifts, regimes, alpha=starting.alpha
  )
  response <- drop(
    kfl1ou:::generate_prediction_vec(
      dat$tree, shifts, regimes, alpha=.8
    ) %*% c(.2, 1.5) +
      stats::rnorm(nrow(dat$Y), sd=.05)
  )
  names(response) <- dat$tree$tip.label

  fit <- suppressWarnings(kfl1ou:::phylolm_CR(
    response ~ preds - 1,
    phy=dat$tree,
    model="OUfixedRoot",
    sc=shifts,
    cr=regimes,
    starting.value=starting.alpha,
    lower.bound=.01,
    upper.bound=2
  ))
  final.design <- kfl1ou:::generate_prediction_vec(
    dat$tree, shifts, regimes, alpha=fit$optpar
  )
  reference <- kfl1ou:::phylolm(
    response ~ final.design - 1,
    phy=dat$tree,
    model="OUfixedRoot",
    starting.value=list(alpha=fit$optpar),
    lower.bound=fit$optpar,
    upper.bound=fit$optpar
  )

  expect_equal(fit$X, final.design, tolerance=1e-12)
  expect_equal(fit$logLik, reference$logLik, tolerance=1e-8)
  expect_equal(
    as.numeric(fit$fitted.values), as.numeric(reference$fitted.values),
    tolerance=1e-8
  )
})


test_that("RR convergent search whitens its intercept predictor", {
  dat <- cr_audit_data(n=12L)
  candidates <- which(
    colSums(dat$Z > 0) >= 2L &
      colSums(dat$Z > 0) <= nrow(dat$Z) - 2L
  )[1:2]
  captured.design <- NULL

  local_mocked_bindings(
    l1ou_require_genlasso=function() NULL,
    getFromNamespace=function(x, ns){
      function(y, X, D, ...){
        captured.design <<- X
        list(
          beta=matrix(0, nrow=ncol(X), ncol=1L),
          lambda=1
        )
      }
    },
    coef=function(object, ...){
      list(beta=rep(0, nrow(object$beta)))
    },
    .package="kfl1ou"
  )

  kfl1ou:::find_convergent_regimes(
    dat$tree, dat$Y, alpha=.7, criterion="AICc",
    regimes=as.list(candidates)
  )
  inverse.factor <- t(sqrt_OU_covariance(
    dat$tree, alpha=.7, root.model="OUfixedRoot"
  )$sqrtInvSigma)
  expected <- drop(inverse.factor %*% rep(1, nrow(dat$Y)))

  expect_equal(
    captured.design[, ncol(captured.design)],
    expected,
    tolerance=1e-12
  )
  expect_false(isTRUE(all.equal(expected, rep(1, length(expected)))))
})


test_that("fixed-alpha convergent fits retain OU weights and fixed bounds", {
  dat <- cr_audit_data(n=16L)
  terminal <- which(dat$tree$edge[, 2L] <= length(dat$tree$tip.label))
  shifts <- terminal[1:2]
  regimes <- list(0L, shifts)
  alpha <- .5
  weighted <- kfl1ou:::generate_prediction_vec(
    dat$tree, shifts, regimes, alpha=alpha
  )
  opt <- list(
    root.model="OUfixedRoot",
    fixed.alpha=TRUE,
    measurement_error=FALSE,
    input_error=NULL,
    alpha.upper.bound=2
  )

  marginal <- kfl1ou:::phylolm_interface_CR(
    dat$tree, dat$Y[, 1, drop=FALSE], regimes,
    alpha=alpha, fixed.alpha=TRUE, opt=opt,
    shift.configuration=shifts
  )
  public <- suppressWarnings(fit_OU(
    dat$tree, dat$Y, shifts,
    cr.regimes=regimes,
    criterion="AICc",
    alpha.lower=alpha,
    alpha.upper=alpha,
    alpha.starting.value=alpha,
    optimizer.starts=1L,
    compute.hessian=FALSE
  ))

  expect_equal(
    as.numeric(marginal$X), as.numeric(weighted),
    tolerance=1e-12
  )
  expect_true(any(weighted > 0 & weighted < 1))
  expect_equal(unname(public$alpha), alpha, tolerance=1e-12)
  expect_true(isTRUE(public$l1ou.options$fixed.alpha))
})


test_that("convergent BIC and AICc use fitted trait dimensions", {
  tree <- ape::reorder.phylo(ape::stree(10L, type="star"), "postorder")
  tree$edge.length <- rep(1, nrow(tree$edge))
  Y <- matrix(
    0, nrow=10, ncol=2,
    dimnames=list(tree$tip.label, c("t1", "t2"))
  )
  opt <- list(
    shift.configuration=c(1L, 2L),
    fixed.alpha=FALSE,
    input_error=NULL
  )

  local_mocked_bindings(
    get_data=function(tree, Y, shift.configuration, opt, idx){
      n <- c(7L, 8L)[[idx]]
      list(
        tr=tree,
        y=matrix(seq_len(n), ncol=1),
        s.c=shift.configuration,
        augmented.s.c=shift.configuration,
        input_error=NULL
      )
    },
    phylolm_interface_CR=function(tr, Y, ...){
      if(nrow(Y) == 7L){
        list(logLik=-10, p=4, n=7)
      } else{
        list(logLik=-20, p=3, n=8)
      }
    },
    .package="kfl1ou"
  )

  bic <- kfl1ou:::cmp_BIC_CR(
    tree, Y, list(0L, c(1L, 2L)), alpha=c(.3, .4), opt=opt
  )
  aicc <- kfl1ou:::cmp_AICc_CR(
    tree, Y, list(0L, c(1L, 2L)), alpha=c(.3, .4), opt=opt
  )
  expected.bic <- 2 * log(10) +
    20 + 4 * log(7) +
    40 + 3 * log(8)
  total.p <- 2 + 4 + 3
  expected.aicc <- 60 + 2 * total.p +
    2 * total.p * (total.p + 1) / (15 - total.p - 1)

  expect_equal(bic, expected.bic, tolerance=1e-12)
  expect_equal(aicc, expected.aicc, tolerance=1e-12)
})


test_that("all explicit candidate edges reach the univariate sparse path", {
  dat <- cr_audit_data(n=10L)
  seen <- integer()

  local_mocked_bindings(
    run_univariate_sparse_path=function(XX, YY, opt){
      seen <<- c(seen, ncol(XX))
      list(
        beta=matrix(0, nrow=1L, ncol=ncol(XX)),
        call="l1ou_audit_path"
      )
    },
    .package="kfl1ou"
  )

  suppressWarnings(estimate_shift_configuration(
    dat$tree, dat$Y,
    max.nShifts=1L,
    candid.edges=seq_len(Nedge(dat$tree)),
    search.strategy="lasso",
    ensemble.replicates=0L,
    optimizer.starts=1L,
    compute.hessian=FALSE,
    quietly=TRUE
  ))

  expect_gt(length(seen), 0L)
  expect_true(all(seen == Nedge(dat$tree)))
})


test_that("fixed-configuration APIs reject non-integer shift indices", {
  dat <- cr_audit_data(n=10L)

  expect_error(
    kfl1ou:::generate_prediction_vec(
      dat$tree, 1.9, list(0L, 1L), alpha=.5
    ),
    "unique, finite integer edge indices"
  )
  expect_error(
    fit_OU(dat$tree, dat$Y, 1.9),
    "unique, finite integer edge indices"
  )
  expect_error(
    configuration_ic(dat$tree, dat$Y, 1.9),
    "unique, finite integer edge indices"
  )
})


test_that("convergent regimes reject lossy or ambiguous edge assignments", {
  dat <- cr_audit_data(n=10L)
  shift <- dat$shift

  expect_error(
    fit_OU(dat$tree, dat$Y, shift, cr.regimes=list(0L, shift + 0.9)),
    "finite integer edge indices"
  )
  expect_error(
    fit_OU(dat$tree, dat$Y, shift, cr.regimes=list(c(0L, shift), shift)),
    "exactly one convergent regime"
  )
  expect_error(
    fit_OU(dat$tree, dat$Y, shift, cr.regimes=list(shift)),
    "background/intercept"
  )
  expect_error(
    fit_OU(dat$tree, dat$Y, shift, cr.regimes=list(0L)),
    "do not match"
  )
  expect_error(
    fit_OU(dat$tree, dat$Y, shift, cr.regimes=c(0L, shift)),
    "supplied as a list"
  )

  named.shifts <- c(2L, 1L)
  names(named.shifts) <- c("shared", "shared")
  expect_identical(
    kfl1ou:::normalize_convergent_regimes(NULL, named.shifts),
    list(0L, c(1L, 2L))
  )
})


test_that("exact multivariate alpha bounds propagate through shift search", {
  dat <- cr_audit_data(n=10L, p=2L)
  captured.opt <- NULL

  local_mocked_bindings(
    resolve_shift_search_strategy=function(tree, opt){
      list(
        requested="exhaustive", strategy="exhaustive",
        configuration.space.size=1
      )
    },
    enumerate_shift_configurations=function(tree, opt){
      list(integer())
    },
    estimate_shift_configuration_from_candidates=function(
        tree, Y, candidate.configurations, opt){
      captured.opt <<- opt
      list(
        score=0,
        profile=list(scores=0, configurations=list(integer()))
      )
    },
    .package="kfl1ou"
  )

  suppressWarnings(estimate_shift_configuration(
    dat$tree, dat$Y,
    criterion="BIC",
    max.nShifts=1L,
    search.strategy="exhaustive",
    trait.covariance="full",
    alpha.structure="diagonal",
    alpha.lower=c(0, 0),
    alpha.upper=c(0, 0),
    alpha.starting.value=c(0, 0),
    optimizer.starts=1L,
    compute.hessian=FALSE,
    quietly=TRUE
  ))

  expect_true(isTRUE(captured.opt$fixed.alpha))
  expect_identical(captured.opt$alpha.lower.bound, c(0, 0))
  expect_identical(captured.opt$alpha.upper.bound, c(0, 0))
})


test_that("alpha-zero shift fits keep mean shifts finite and mark optima unknown", {
  diagonal <- cr_audit_data(n=16L, p=1L)
  diagonal$Y[, 1] <- 2 * diagonal$Z[, diagonal$shift] +
    stats::rnorm(nrow(diagonal$Y), sd=.1)
  fit.diagonal <- fit_OU(
    diagonal$tree, diagonal$Y, diagonal$shift,
    criterion="BIC",
    alpha.lower=0,
    alpha.upper=0,
    alpha.starting.value=0,
    optimizer.starts=1L,
    compute.hessian=FALSE
  )

  full <- cr_audit_data(n=18L, p=2L, seed=4102L)
  full$Y[, 1] <- 1.5 * full$Z[, full$shift] +
    stats::rnorm(nrow(full$Y), sd=.2)
  full$Y[, 2] <- -.8 * full$Z[, full$shift] +
    stats::rnorm(nrow(full$Y), sd=.2)
  fit.full <- fit_OU(
    full$tree, full$Y, full$shift,
    criterion="BIC",
    trait.covariance="full",
    alpha.lower=0,
    alpha.upper=0,
    alpha.starting.value=0,
    likelihood.engine="dense",
    optimizer.starts=1L,
    compute.hessian=FALSE
  )

  for(fit in list(fit.diagonal, fit.full)){
    expect_true(all(is.finite(fit$shift.means)))
    expect_true(all(is.na(fit$shift.values)))
    expect_false(any(is.infinite(fit$optima), na.rm=TRUE))
    expect_false(any(is.nan(fit$optima)))
    expect_true(all(is.finite(fit$mu)))
  }
})


test_that("missing input errors are zero only in latent unobserved cells", {
  Y <- matrix(
    c(1, NA, 2, 3, 4, NA), nrow=3, ncol=2,
    dimnames=list(paste0("s", 1:3), c("t1", "t2"))
  )
  input.error <- matrix(.1, nrow=3, ncol=2, dimnames=dimnames(Y))
  input.error[is.na(Y)] <- NA_real_
  result <- kfl1ou:::multivariate_observation_error_vector(
    Y, list(input_error=input.error), sigma2.error=c(.2, .3)
  )
  result <- matrix(result, nrow=3, ncol=2)

  expect_true(all(is.finite(result)))
  expect_equal(result[is.na(Y)], c(.2, .3))
  expect_equal(result[!is.na(Y)], input.error[!is.na(Y)] +
    matrix(rep(c(.2, .3), each=3), 3, 2)[!is.na(Y)])
})
