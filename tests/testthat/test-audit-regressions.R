audit_small_traits <- function(n=14L, p=2L, seed=901L){
  set.seed(seed)
  tree <- ape::reorder.phylo(ape::rcoal(n), "postorder")
  Y <- matrix(rnorm(n * p), n, p,
              dimnames=list(tree$tip.label, paste0("t", seq_len(p))))
  list(tree=tree, Y=Y)
}

test_that("fixed alpha never borrows covariance curvature for Wald intervals", {
  dat <- audit_small_traits()
  fit <- fit_OU(
    dat$tree, dat$Y, integer(), criterion="BIC",
    trait.covariance="full", alpha.structure="diagonal",
    alpha.lower=c(.35, 1.1), alpha.upper=c(.35, 1.1),
    likelihood.engine="dense", optimizer.starts=1L
  )
  interval <- confint(
    fit, parm=c("alpha:t1", "alpha:t2"), method="wald"
  )
  expect_true(all(is.na(interval)))
  expect_length(fit$optimization$alpha.parameter.index, 0L)
})

test_that("warning options are restored after likelihood errors", {
  dat <- audit_small_traits(n=8L, p=1L)
  opt <- list(
    measurement_error=FALSE, input_error=NULL,
    shift.configuration=integer(), fixed.alpha=TRUE,
    root.model="invalid-root", alpha.upper.bound=2
  )
  old <- getOption("warn")
  on.exit(options(warn=old), add=TRUE)
  options(warn=1)
  expect_error(
    kfl1ou:::phylolm_interface_CR(
      dat$tree, dat$Y, conv.regimes=list(0L), alpha=.5,
      fixed.alpha=TRUE, opt=opt
    )
  )
  expect_equal(getOption("warn"), 1)
})

test_that("root polytomies are resolved safely", {
  tree <- ape::stree(3L, type="star")
  tree$edge.length <- rep(1, nrow(tree$edge))
  result <- sqrt_OU_covariance(tree)
  expected <- ape::vcv.phylo(tree)
  dimnames(expected) <- NULL
  expect_equal(
    result$sqrtSigma %*% t(result$sqrtSigma), expected,
    tolerance=1e-10
  )
})

test_that("covariance native boundary rejects malformed edge lists", {
  valid <- cbind(
    ancestor = c(2, 2),
    descendant = c(0, 1),
    length = c(1, 1)
  )
  result <- kfl1ou:::cmp_sqrt_OU_covariance(valid, 2L, 0)
  expect_equal(dim(result$sqrtSigma), c(2L, 2L))
  expect_true(all(is.finite(result$sqrtSigma)))

  bad_node <- valid
  bad_node[1L, "descendant"] <- 3
  expect_error(
    kfl1ou:::cmp_sqrt_OU_covariance(bad_node, 2L, 0),
    "invalid node identifiers or lengths"
  )

  duplicate_tip <- valid
  duplicate_tip[, "descendant"] <- 0
  expect_error(
    kfl1ou:::cmp_sqrt_OU_covariance(duplicate_tip, 2L, 0),
    "invalid tree topology"
  )

  invalid_length <- valid
  invalid_length[1L, "length"] <- Inf
  expect_error(
    kfl1ou:::cmp_sqrt_OU_covariance(invalid_length, 2L, 0),
    "invalid node identifiers or lengths"
  )

  expect_error(
    kfl1ou:::cmp_sqrt_OU_covariance(valid, .Machine$integer.max, 0),
    "too large or inconsistent"
  )
})

test_that("multivariate model plots exercise all annotation paths", {
  dat <- audit_small_traits(n=8L)
  design <- kfl1ou:::generate_design_matrix(dat$tree, "simpX")
  candidate <- which(
    colSums(design) >= 2L & colSums(design) <= nrow(design) - 2L
  )[[1L]]
  fit <- fit_OU(
    dat$tree,
    dat$Y,
    candidate,
    criterion="BIC",
    alpha.lower=.5,
    alpha.upper=.5
  )
  fit$tree$edge.label <- paste0("edge-", seq_len(nrow(fit$tree$edge)))
  output <- tempfile(fileext=".pdf")
  grDevices::pdf(output)
  on.exit(grDevices::dev.off(), add=TRUE)

  expect_error(
    plot(
      fit,
      edge.label.ann=TRUE,
      edge.label.pos=.5,
      plot.bar=TRUE,
      asterisk=TRUE
    ),
    NA
  )
  expect_true(file.exists(output))
})

test_that("explicit engines and shrinkage semantics are stable", {
  dat <- audit_small_traits(n=12L)
  dense <- fit_OU(
    dat$tree, dat$Y, integer(), criterion="BIC",
    trait.covariance="full", alpha.lower=.6, alpha.upper=.6,
    likelihood.engine="dense", optimizer.starts=1L,
    covariance.regularization="shrinkage"
  )
  missing <- dat$Y
  missing[1, 1] <- NA_real_
  dense.missing <- fit_OU(
    dat$tree, missing, integer(), criterion="BIC",
    trait.covariance="full", alpha.lower=.6, alpha.upper=.6,
    likelihood.engine="dense", optimizer.starts=1L,
    covariance.regularization="shrinkage"
  )
  expect_equal(dense$likelihood.engine, "dense")
  expect_equal(dense.missing$likelihood.engine, "dense")
  expect_equal(dense$regularization.lambda,
               dense.missing$regularization.lambda)
  expect_false(dense$information.criterion.calibrated)
})

test_that("seeded simulation restores RNG and tree simulation matches shape", {
  dat <- audit_small_traits(n=10L)
  fit <- fit_OU(
    dat$tree, dat$Y, integer(), criterion="BIC",
    trait.covariance="full", alpha.lower=.7, alpha.upper=.7,
    optimizer.starts=1L
  )
  set.seed(123)
  before <- .Random.seed
  simulations <- simulate(fit, nsim=3L, seed=99)
  expect_identical(.Random.seed, before)
  expect_length(simulations, 3L)
  expect_equal(dim(simulations[[1L]]), dim(dat$Y))
})

test_that("normalization scale and original-time parameters are retained", {
  dat <- audit_small_traits(n=8L, p=1L)
  original.height <- max(ape::node.depth.edgelength(dat$tree)[seq_len(8L)])
  adjusted <- adjust_data(dat$tree, dat$Y, quietly=TRUE)
  expect_equal(adjusted$tree.scale, original.height)
  fit <- fit_OU(adjusted$tree, adjusted$Y, integer(), criterion="BIC")
  converted <- original_time_parameters(fit)
  expect_equal(converted$tree.scale, original.height)
  expect_equal(converted$alpha, as.numeric(fit$alpha) / original.height)
})

test_that("invalid data and duplicate labels fail at the public boundary", {
  dat <- audit_small_traits(n=8L)
  bad <- dat$Y
  bad[1, 1] <- Inf
  expect_error(adjust_data(dat$tree, bad), "infinite")
  bad[1, 1] <- NaN
  expect_error(adjust_data(dat$tree, bad), "NaN")
  duplicate <- dat$tree
  duplicate$tip.label[2] <- duplicate$tip.label[1]
  expect_error(adjust_data(duplicate, dat$Y), "unique")

  input.error <- matrix(.01, nrow(dat$Y), ncol(dat$Y),
                        dimnames=dimnames(dat$Y))
  input.error[1, 1] <- Inf
  expect_error(
    fit_OU(dat$tree, dat$Y, integer(), criterion="BIC",
           input_error=input.error),
    "finite variances"
  )
})

test_that("malformed trees fail before native traversal", {
  dat <- audit_small_traits(n=8L, p=1L)
  malformed <- dat$tree
  malformed$edge[2, 2] <- malformed$edge[1, 2]
  expect_error(adjust_data(malformed, dat$Y), "edge matrix")
})

test_that("alpha profiles include the Brownian boundary", {
  dat <- audit_small_traits(n=10L, p=1L)
  fit <- fit_OU(dat$tree, dat$Y, integer(), criterion="BIC")
  fixed <- fit_OU(
    dat$tree, dat$Y, integer(), criterion="BIC",
    alpha.lower=.5, alpha.upper=.5
  )
  expect_equal(fixed$parameter.count, fit$parameter.count - 1L)
  expect_equal(kfl1ou:::l1ou_alpha_parameter_names(fit), "alpha:t1")
  profile <- profile_alpha_l1ou(fit, c(0, .5))
  expect_equal(profile$shared, c(0, .5))
  expect_true(is.finite(profile$logLik[[1L]]))
})

test_that("effective sample size rejects unsafe edge indices", {
  tree <- reorder(ape::stree(4, type="balanced"), "pruningwise")
  tree$edge.length <- rep(1, nrow(tree$edge))

  expect_error(
    kfl1ou:::effective.sample.size(tree, edges=c(1, 1)),
    "unique"
  )
  expect_error(
    kfl1ou:::effective.sample.size(tree, edges=nrow(tree$edge) + 1L),
    "non-root"
  )
  expect_true(all(is.finite(
    kfl1ou:::effective.sample.size(tree, edges=c(1L, 2L))
  )))
})

test_that("effective sample size native boundary validates dimensions and topology", {
  valid <- list(
    N = 2L,
    n = 2L,
    pN = 1L,
    root = 3L,
    transa = 1,
    transb = c(1, 1),
    des = c(1L, 2L),
    anc = c(3L, 3L),
    edge = 3L
  )

  result <- do.call(kfl1ou:::effective_sample_size_c, valid)
  expect_length(result, 1L)
  expect_true(is.finite(result))

  bad_dimensions <- valid
  bad_dimensions$N <- 3L
  expect_error(
    do.call(kfl1ou:::effective_sample_size_c, bad_dimensions),
    "dimensions or root"
  )

  overflowing_dimensions <- valid
  overflowing_dimensions$n <- .Machine$integer.max
  expect_error(
    do.call(kfl1ou:::effective_sample_size_c, overflowing_dimensions),
    "dimensions or root"
  )

  duplicate_tip <- valid
  duplicate_tip$des <- c(1L, 1L)
  expect_error(
    do.call(kfl1ou:::effective_sample_size_c, duplicate_tip),
    "topology or edge lengths"
  )

  invalid_length <- valid
  invalid_length$transb[[1L]] <- Inf
  expect_error(
    do.call(kfl1ou:::effective_sample_size_c, invalid_length),
    "topology or edge lengths"
  )

  missing_sentinel <- valid
  missing_sentinel$edge <- 1L
  expect_error(
    do.call(kfl1ou:::effective_sample_size_c, missing_sentinel),
    "root-edge sentinel"
  )
})

test_that("threepoint native boundary validates dimensions and topology", {
  valid <- list(
    N = 2L,
    n = 2L,
    pN = 1L,
    dY = 1L,
    dX = 1L,
    root = 3L,
    transa = 1,
    transb = c(1, 1),
    des = c(1L, 2L),
    anc = c(3L, 3L),
    y = c(1, 2),
    X = c(1, 1)
  )

  result <- do.call(kfl1ou:::threepoint_l1ou_c, valid)
  expect_length(result, 7L)
  expect_true(all(is.finite(result)))

  bad_dimensions <- valid
  bad_dimensions$N <- 3L
  expect_error(
    do.call(kfl1ou:::threepoint_l1ou_c, bad_dimensions),
    "dimensions or root"
  )

  bad_vectors <- valid
  bad_vectors$X <- 1
  expect_error(
    do.call(kfl1ou:::threepoint_l1ou_c, bad_vectors),
    "response or design dimensions"
  )

  oversized <- valid
  oversized$dY <- .Machine$integer.max
  expect_error(
    do.call(kfl1ou:::threepoint_l1ou_c, oversized),
    "working memory is too large"
  )

  duplicate_tip <- valid
  duplicate_tip$des <- c(1L, 1L)
  expect_error(
    do.call(kfl1ou:::threepoint_l1ou_c, duplicate_tip),
    "topology is invalid"
  )

  non_finite <- valid
  non_finite$y[[1L]] <- Inf
  expect_error(
    do.call(kfl1ou:::threepoint_l1ou_c, non_finite),
    "response contains non-finite"
  )
})
