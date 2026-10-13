test_that("AIC uses the declared parameter count without the AICc correction", {
  tree <- ape::read.tree(text="((a:1,b:1):1,(c:1,d:1):1);")
  trait.results <- list(list(fit=list(p=3, n=4, logLik=-7, optpar=.5, sigma2=2)))
  score <- kfl1ou:::summarize_trait_fit_results(tree, 1L, list(), trait.results, "AIC")
  expect_equal(score$score, 22)
  expect_true(all(is.na(kfl1ou:::summarize_trait_fit_results(tree, 1L, list(), trait.results, "AICc"))))
  expect_equal(kfl1ou:::multivariate_full_information_score(list(p=8, logLik=-20, sample.size=10), 2, "AIC"), 60)
})

test_that("public fixed-layout AIC differs from AICc only by its correction", {
  tree <- ape::reorder.phylo(ape::read.tree(text="(((a:1,b:1):1,(c:1,d:1):1):1,((e:1,f:1):1,(g:1,h:1):1):1);"), "postorder")
  Y <- matrix(c(.2, 1, -.3, .5, 2, 1.4, .7, -.1), ncol=1, dimnames=list(tree$tip.label,"x"))
  a <- fit_OU(tree, Y, integer(), criterion="AIC", alpha.lower=.4, alpha.upper=.4, compute.hessian=FALSE)
  b <- fit_OU(tree, Y, integer(), criterion="AICc", alpha.lower=.4, alpha.upper=.4, compute.hessian=FALSE)
  expect_equal(a$score, -2 * sum(a$logLik) + 4)
  expect_equal(b$score-a$score, 12/5)
  selected <- estimate_shift_configuration(tree, Y, criterion="AIC", max.nShifts=1, alpha.lower=.4, alpha.upper=.4, quietly=TRUE)
  expect_true(is.finite(selected$score))
})


test_that("full covariance AIC can be refitted by the sensitivity interface", {
  set.seed(11)
  tree <- ape::reorder.phylo(ape::rcoal(12), "postorder")
  Y <- matrix(rnorm(24),12,2,dimnames=list(tree$tip.label,c("x","y")))
  fit <- fit_OU(tree,Y,integer(),criterion="AIC",trait.covariance="full",
                alpha.lower=.7,alpha.upper=.7,optimizer.starts=1L,compute.hessian=FALSE)
  expect_true(is.finite(fit$score))
  sensitivity <- sensitivity_l1ou(fit,root.models="OUfixedRoot",alpha.upper=.7,criteria="AIC")
  expect_true(all(is.finite(sensitivity$score)))
  conv <- estimate_convergent_regimes(fit,criterion="AIC")
  expect_true(is.finite(conv$score))
})

test_that("AIC counts mean, covariance and location parameters actually fitted", {
  set.seed(37)
  tree <- ape::reorder.phylo(ape::rcoal(16), "postorder")
  Y <- matrix(rnorm(32), 16, 2,
              dimnames=list(tree$tip.label, c("x", "y")))
  shift <- match(1L, tree$edge[, 2L])
  known.error <- setNames(seq(.01, .04, length.out=16), tree$tip.label)
  missing <- Y
  missing[1L, "y"] <- NA_real_
  # Two means per complete trait, one residual variance per trait, and one
  # shared shift location. An absent shifted tip removes that trait's mean.
  cases <- list(
    list(data=Y, error=NULL, estimated.error=FALSE, count=7L),
    list(data=Y, error=known.error, estimated.error=FALSE, count=7L),
    list(data=Y, error=known.error, estimated.error=TRUE, count=9L),
    list(data=missing, error=NULL, estimated.error=FALSE, count=6L)
  )
  for(case in cases){
    fit <- fit_OU(tree, case$data, shift, criterion="AIC",
                  alpha.lower=.7, alpha.upper=.7, input_error=case$error,
                  measurement_error=case$estimated.error, compute.hessian=FALSE)
    expect_equal(fit$information.parameter.count, case$count)
    expect_equal(fit$score, -2 * sum(fit$logLik) + 2 * case$count)
  }
  # Four mean coefficients, three entries in a 2 x 2 symmetric covariance,
  # and one shift location. A shared estimated alpha adds one parameter.
  for(fixed in c(TRUE, FALSE)){
    fit <- fit_OU(tree, Y, shift, criterion="AIC", trait.covariance="full",
                  alpha.lower=if(fixed) .7 else .2, alpha.upper=.7,
                  optimizer.starts=1L, compute.hessian=FALSE)
    count <- 8L + as.integer(!fixed)
    expect_equal(fit$information.parameter.count, count)
    expect_equal(fit$score, -2 * sum(fit$logLik) + 2 * count)
  }
})

test_that("convergent AIC counts shared optima and fixed alpha correctly", {
  set.seed(38)
  tree <- ape::reorder.phylo(ape::rcoal(16), "postorder")
  Y <- matrix(rnorm(16), 16, 1, dimnames=list(tree$tip.label, "x"))
  shifts <- match(c(1L, 2L), tree$edge[, 2L])
  # Two optimum coefficients and a residual variance, plus two locations.
  for(fixed in c(TRUE, FALSE)){
    fit <- fit_OU(tree, Y, shifts, cr.regimes=list(0L, shifts), criterion="AIC",
                  alpha.lower=if(fixed) .7 else .2, alpha.upper=.7,
                  optimizer.starts=1L, compute.hessian=FALSE)
    count <- 5L + as.integer(!fixed)
    expect_equal(fit$information.parameter.count, count)
    expect_equal(fit$score, -2 * sum(fit$logLik) + 2 * count)
  }
})
