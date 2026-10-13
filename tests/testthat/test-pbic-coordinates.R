pbic_coordinate_fixture <- function(){
  tree <- ape::reorder.phylo(ape::read.tree(text=paste0(
    "(((t0:1,t1:1):1,(t2:1,t3:1):1):1,",
    "((t4:1,t5:1):1,(t6:1,t7:1):1):1);"
  )), "postorder")
  values <- c(1.5549029785536523, 1.3137883823887941,
              -.21225628554971665, -.24002282307176648,
              1.7035058708492061, 1.8521601652123756,
              .7166573577624096, .2547434554367455)
  Y <- matrix(values, ncol=1L,
              dimnames=list(paste0("t", 0:7), "value"))
  shifts <- match(match(c("t0", "t4"), tree$tip.label), tree$edge[, 2L])
  list(tree=tree, Y=Y[tree$tip.label, , drop=FALSE], shifts=shifts)
}

# Independent branchwise expectation recursion, not generate_prediction_vec().
pbic_reference_design <- function(tree, shifts, regimes, alpha){
  n <- length(tree$tip.label)
  depth <- ape::node.depth.edgelength(tree)
  root <- setdiff(tree$edge[, 1L], tree$edge[, 2L])
  weights <- matrix(0, n + tree$Nnode, length(regimes))
  group <- integer(n + tree$Nnode)
  group[root] <- which(vapply(regimes, function(x) 0L %in% x, logical(1)))
  weights[root, group[root]] <- 1
  for(edge in order(depth[tree$edge[, 1L]])){
    parent <- tree$edge[edge, 1L]
    child <- tree$edge[edge, 2L]
    group[child] <- if(edge %in% shifts){
      which(vapply(regimes, function(x) edge %in% x, logical(1)))
    } else group[parent]
    decay <- exp(-alpha * tree$edge.length[edge])
    weights[child, ] <- decay * weights[parent, ]
    weights[child, group[child]] <- weights[child, group[child]] - expm1(-alpha * tree$edge.length[edge])
  }
  weights[seq_len(n), , drop=FALSE]
}

pbic_dense_reference <- function(fit, regimes, fixed=FALSE, error=NULL,
                                 estimated.error=FALSE){
  tree <- fit$tree
  n <- length(tree$tip.label)
  shared <- ape::vcv.phylo(tree)[tree$tip.label, tree$tip.label]
  distance <- outer(diag(shared), diag(shared), "+") - 2 * shared
  alpha <- unname(fit$alpha[[1L]])
  Sigma <- fit$sigma2[[1L]]/(2 * alpha) * exp(-alpha * distance)
  if(fit$l1ou.options$root.model == "OUfixedRoot"){
    Sigma <- Sigma * -expm1(-2 * alpha * shared)
  }
  if(!is.null(error)) Sigma <- Sigma + diag(error[tree$tip.label])
  if(estimated.error) Sigma <- Sigma + diag(fit$sigma2_error[[1L]], n)
  X <- pbic_reference_design(tree, unname(fit$shift.configuration), regimes, alpha)
  whitened <- forwardsolve(t(chol(Sigma)), X)
  # abs(R diagonal) avoids dependence on QR sign conventions.
  log.information <- 2 * sum(log(abs(diag(qr.R(qr(whitened))))))
  log.information <- log.information + ncol(X) * log(var(as.numeric(fit$Y)))
  -2 * sum(fit$logLik) + 2 * length(fit$shift.configuration) *
    log(nrow(tree$edge) - 1L) +
    (2L - as.integer(fixed) + as.integer(estimated.error)) * log(n) + log.information
}

pbic_public_fit <- function(dat, shifts=dat$shifts, regimes=NULL,
                            root="OUfixedRoot", fixed=FALSE, error=NULL, upper=log(2)/3){
  suppressWarnings(fit_OU(
    dat$tree, dat$Y, shifts, cr.regimes=regimes, criterion="pBIC",
    root.model=root, alpha.lower=if(fixed) .4 else 1e-7/3,
    alpha.upper=if(fixed) .4 else upper,
    compute.hessian=FALSE, optimizer.starts=1L, input_error=error
  ))
}


test_that("pBIC equals the optimum-coordinate dense information formula", {
  dat <- pbic_coordinate_fixture()
  for(root in c("OUfixedRoot", "OUrandomRoot")){
    for(fixed in c(FALSE, TRUE)){
      fit <- pbic_public_fit(dat, root=root, fixed=fixed)
      expect_equal(fit$score, pbic_dense_reference(fit, as.list(c(0L, dat$shifts)), fixed),
                   tolerance=2e-6, info=paste(root, fixed))
    }
  }
})


test_that("free and singleton pBIC agree for zero, disconnected and nested shifts", {
  dat <- pbic_coordinate_fixture()
  ancestor <- ape::getMRCA(dat$tree, c("t0", "t1"))
  nested <- c(match(ancestor, dat$tree$edge[, 2L]), dat$shifts[[1L]])
  for(shifts in list(integer(), dat$shifts, nested)){
    for(root in c("OUfixedRoot", "OUrandomRoot")){
      for(error in list(NULL, setNames(seq(.01, .08, length.out=8L), dat$tree$tip.label))){
        free <- pbic_public_fit(dat, shifts, root=root, fixed=TRUE, error=error)
        singleton <- pbic_public_fit(dat, shifts, as.list(c(0L, shifts)),
                                     root=root, fixed=TRUE, error=error)
        expect_equal(singleton$logLik, free$logLik, tolerance=1e-6)
        expect_equal(singleton$score, free$score, tolerance=1e-6)
      }
    }
  }
})


test_that("merged pBIC evaluates information at the constrained fitted alpha", {
  dat <- pbic_coordinate_fixture()
  dat$Y[, 1L] <- c(.65104646774656, -.904417507797141, -.0885200215018778,
                   .75323985426968, -.889038820862497, .446587397080797,
                   -.731843206541243, -.256717829917734)
  regimes <- list(0L, dat$shifts)
  for(root in c("OUfixedRoot", "OUrandomRoot")){
    fit <- pbic_public_fit(dat, regimes=regimes, root=root, upper=3)
    expect_equal(fit$score, pbic_dense_reference(fit, regimes), tolerance=1e-6)
    permuted <- pbic_public_fit(dat, shifts=rev(dat$shifts), regimes=rev(regimes), root=root, upper=3)
    expect_equal(permuted$score, fit$score, tolerance=1e-6)
  }
})


test_that("estimated-alpha singleton scores and missing-trait mappings agree", {
  dat <- pbic_coordinate_fixture()
  for(root in c("OUfixedRoot", "OUrandomRoot")){
    free <- pbic_public_fit(dat, root=root)
    singleton <- pbic_public_fit(dat, regimes=as.list(c(0L, dat$shifts)), root=root)
    expect_equal(singleton$score, free$score, tolerance=2e-6)
  }
  dat$Y <- cbind(dat$Y, second=rev(dat$Y[, 1L]) + seq_len(8L)/20)
  dat$Y["t0", "second"] <- NA_real_
  free <- pbic_public_fit(dat, fixed=TRUE)
  singleton <- pbic_public_fit(dat, regimes=as.list(c(0L, dat$shifts)), fixed=TRUE)
  expect_equal(singleton$score, free$score, tolerance=1e-6)
  expect_equal(singleton$logLik, free$logLik, tolerance=1e-6)
})


test_that("shifted exact Brownian fits do not acquire negative infinite pBIC", {
  dat <- pbic_coordinate_fixture()
  fit <- fit_OU(dat$tree, dat$Y, dat$shifts, criterion="pBIC",
                alpha.lower=0, alpha.upper=0, compute.hessian=FALSE)
  expect_identical(fit$score, Inf)
  expect_true(all(is.finite(fit$mu)))
  null <- fit_OU(dat$tree, dat$Y, integer(), criterion="pBIC",
                 alpha.lower=0, alpha.upper=0, compute.hessian=FALSE)
  expect_true(is.finite(null$score))
})

test_that("pBIC dense information agrees for unequal branch ages and observation error", {
  set.seed(39)
  dat <- list(tree=ape::reorder.phylo(ape::rcoal(12), "postorder"))
  dat$Y <- matrix(rnorm(12), 12, 1, dimnames=list(dat$tree$tip.label, "x"))
  dat$shifts <- match(c(1L, 3L), dat$tree$edge[, 2L])
  error <- setNames(seq(.01, .04, length.out=12), dat$tree$tip.label)
  for(root in c("OUfixedRoot", "OUrandomRoot")){
    for(regimes in list(NULL, list(0L, dat$shifts))){
      for(estimated.error in c(FALSE, TRUE)){
        fit <- fit_OU(dat$tree, dat$Y, dat$shifts, cr.regimes=regimes,
                      criterion="pBIC", root.model=root,
                      alpha.lower=.4, alpha.upper=.4, input_error=error,
                      measurement_error=estimated.error, compute.hessian=FALSE)
        design.regimes <- if(is.null(regimes)) as.list(c(0L, dat$shifts)) else regimes
        expect_equal(fit$score,
                     pbic_dense_reference(fit, design.regimes, fixed=TRUE,
                                          error=error, estimated.error=estimated.error),
                     tolerance=2e-6)
      }
    }
  }
})
