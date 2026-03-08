small_lizard_data <- function(n_tips = 12, traits = 1, normalize = TRUE) {
  env <- new.env(parent = emptyenv())
  data("lizard.tree", package = "kfl1ou", envir = env)
  data("lizard.traits", package = "kfl1ou", envir = env)

  keep <- env$lizard.tree$tip.label[seq_len(n_tips)]
  tree <- ape::drop.tip(env$lizard.tree, setdiff(env$lizard.tree$tip.label, keep))
  Y <- as.matrix(env$lizard.traits[keep, traits, drop = FALSE])

  adjust_data(tree, Y, normalize = normalize, quietly = TRUE)
}

capture_silently <- function(expr) {
  invisible(capture.output(result <- eval.parent(substitute(expr))))
  result
}

named_input_error <- function(Y, value = 0.01) {
  stats::setNames(rep(value, nrow(Y)), rownames(Y))
}

direct_joint_error_intercept_fit <- function(tree, Y, alpha, input_error,
                                             root.model = "OUfixedRoot") {
  y <- as.numeric(as.matrix(Y)[, 1])
  X <- matrix(1, nrow = length(y), ncol = 1)
  base_cov <- {
    re <- kfl1ou:::sqrt_OU_covariance(
      tree,
      alpha = alpha,
      root.model = root.model,
      sigma2 = 1,
      check.order = FALSE,
      check.ultrametric = FALSE
    )
    tcrossprod(re$sqrtSigma)
  }
  tip_error <- input_error[tree$tip.label]

  objective <- function(par) {
    sigma2 <- exp(par[[1]])
    sigma2_error <- par[[2]]
    Sigma <- sigma2 * base_cov
    diag(Sigma) <- diag(Sigma) + tip_error + sigma2_error
    chol_sigma <- tryCatch(chol(Sigma), error = function(e) NULL)
    if (is.null(chol_sigma)) {
      return(.Machine$double.xmax / 1000)
    }

    Xt <- forwardsolve(t(chol_sigma), X)
    yt <- drop(forwardsolve(t(chol_sigma), matrix(y, ncol = 1)))
    XX <- crossprod(Xt)
    Xy <- crossprod(Xt, yt)
    beta <- drop(solve(XX, Xy))
    rss <- sum((yt - drop(Xt %*% beta))^2)
    logdet <- 2 * sum(log(diag(chol_sigma)))
    as.numeric(length(y) * log(2 * pi) + logdet + rss)
  }

  opt <- stats::optim(
    c(log(max(stats::var(y), 1e-8)), 0),
    objective,
    method = "L-BFGS-B",
    lower = c(log(.Machine$double.eps), 0),
    upper = c(Inf, Inf)
  )

  sigma2 <- exp(opt$par[[1]])
  sigma2_error <- opt$par[[2]]
  Sigma <- sigma2 * base_cov
  diag(Sigma) <- diag(Sigma) + tip_error + sigma2_error
  chol_sigma <- chol(Sigma)
  Xt <- forwardsolve(t(chol_sigma), X)
  yt <- drop(forwardsolve(t(chol_sigma), matrix(y, ncol = 1)))
  beta <- drop(solve(crossprod(Xt), crossprod(Xt, yt)))

  list(
    coefficients = beta,
    sigma2 = sigma2,
    sigma2_error = sigma2_error,
    logLik = -objective(opt$par) / 2
  )
}
