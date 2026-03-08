univariate_beta_path_result <- function(beta, backend) {

    beta <- as.matrix(beta)
    if (nrow(beta) == 0) {
        beta <- matrix(0, nrow = 1, ncol = ncol(beta))
    }

    list(
        beta = beta,
        call = paste0("l1ou_", backend, "_path")
    )
}

run_univariate_lasso_path <- function(XX, YY, opt) {

    if (ncol(XX) == 0L || opt$max.nShifts <= 0L) {
        return(univariate_beta_path_result(matrix(0, nrow = 1, ncol = ncol(XX)), "lasso"))
    }

    sol <- run_grplasso(
        grpX = XX,
        grpY = as.numeric(YY),
        nVariables = 1L,
        grpIdx = seq_len(ncol(XX)),
        opt = opt
    )

    beta.path <- t(as.matrix(sol$coefficients))
    if (nrow(beta.path) == 0 || any(beta.path[1, ] != 0)) {
        beta.path <- rbind(0, beta.path)
    }

    univariate_beta_path_result(beta.path, "lasso")
}

run_univariate_stepwise_path <- function(XX, YY, opt) {

    XX <- as.matrix(XX)
    YY <- as.numeric(YY)
    nPred <- ncol(XX)
    max.steps <- min(as.integer(opt$max.nShifts), nPred)

    beta.path <- matrix(0, nrow = max.steps + 1L, ncol = nPred)
    active <- integer()
    ignored <- integer()
    residual <- YY

    for (step in seq_len(max.steps)) {
        inactive <- setdiff(seq_len(nPred), c(active, ignored))
        if (length(inactive) == 0L) {
            beta.path <- beta.path[seq_len(step), , drop = FALSE]
            break
        }

        ranked <- inactive[order(abs(drop(crossprod(XX[, inactive, drop = FALSE], residual))),
                                 decreasing = TRUE)]
        chosen <- NA_integer_
        for (candidate in ranked) {
            proposal <- c(active, candidate)
            proposal.x <- XX[, proposal, drop = FALSE]
            if (qr(proposal.x)$rank < ncol(proposal.x)) {
                ignored <- c(ignored, candidate)
                next
            }
            chosen <- candidate
            break
        }

        if (is.na(chosen)) {
            beta.path <- beta.path[seq_len(step), , drop = FALSE]
            break
        }

        active <- c(active, chosen)
        fit <- stats::lm.fit(x = XX[, active, drop = FALSE], y = YY)
        beta <- numeric(nPred)
        beta[active] <- fit$coefficients
        beta[is.na(beta)] <- 0
        beta.path[step + 1L, ] <- beta
        residual <- YY - drop(XX %*% beta)
    }

    univariate_beta_path_result(beta.path, "stepwise")
}

run_univariate_sparse_path <- function(XX, YY, opt) {

    if (identical(opt$lars.alg, "stepwise")) {
        return(run_univariate_stepwise_path(XX, YY, opt))
    }

    run_univariate_lasso_path(XX, YY, opt)
}
