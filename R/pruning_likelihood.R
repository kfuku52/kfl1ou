ou_branch_innovation_covariance <- function(alpha, trait.covariance,
                                             branch.length){
    alpha <- as.numeric(alpha)
    p <- nrow(trait.covariance)
    alpha <- rep(alpha, length.out=p)
    rates <- outer(alpha, alpha, "+")
    multiplier <- matrix(branch.length, p, p)
    positive <- rates > .Machine$double.eps
    multiplier[positive] <- -expm1(-rates[positive] * branch.length) /
        rates[positive]
    covariance <- trait.covariance * multiplier
    0.5 * (covariance + t(covariance))
}

ou_stationary_trait_covariance <- function(alpha, trait.covariance){
    alpha <- rep(as.numeric(alpha), length.out=nrow(trait.covariance))
    rates <- outer(alpha, alpha, "+")
    if(any(rates <= .Machine$double.eps)){
        stop("OUrandomRoot requires strictly positive alpha values.")
    }
    covariance <- trait.covariance / rates
    0.5 * (covariance + t(covariance))
}

pruning_spd_components <- function(covariance){
    covariance <- 0.5 * (covariance + t(covariance))
    factor <- tryCatch(chol(covariance), error=function(e) NULL)
    if(is.null(factor)){
        covariance <- regularize_positive_definite(covariance, 1e-12)
        factor <- chol(covariance)
    }
    list(
        inverse=chol2inv(factor),
        log.determinant=2 * sum(log(diag(factor)))
    )
}

transform_pruning_message <- function(message, transition, innovation){
    q <- pruning_spd_components(innovation)
    precision <- 0.5 * (message$precision + t(message$precision))
    combined <- pruning_spd_components(precision + q$inverse)
    middle <- q$inverse - q$inverse %*% combined$inverse %*% q$inverse
    list(
        precision=0.5 * (t(transition) %*% middle %*% transition +
                          t(t(transition) %*% middle %*% transition)),
        linear=drop(t(transition) %*% q$inverse %*%
                    combined$inverse %*% message$linear),
        constant=message$constant - 0.5 * q$log.determinant -
            0.5 * combined$log.determinant +
            0.5 * drop(crossprod(
                message$linear,
                combined$inverse %*% message$linear
            ))
    )
}

tip_pruning_message <- function(value, transition, innovation,
                                observation.error){
    observed <- !is.na(value)
    p <- length(value)
    if(!any(observed)){
        return(list(
            precision=matrix(0, p, p), linear=numeric(p), constant=0
        ))
    }
    selection <- diag(p)[observed, , drop=FALSE]
    covariance <- selection %*% innovation %*% t(selection)
    diag(covariance) <- diag(covariance) + observation.error[observed]
    components <- pruning_spd_components(covariance)
    loading <- selection %*% transition
    y <- value[observed]
    list(
        precision=t(loading) %*% components$inverse %*% loading,
        linear=drop(t(loading) %*% components$inverse %*% y),
        constant=-0.5 * length(y) * log(2 * pi) -
            0.5 * components$log.determinant -
            0.5 * drop(crossprod(y, components$inverse %*% y))
    )
}

# Exact Gaussian pruning likelihood for diagonal-drift multivariate OU models.
pruning_multivariate_ou_loglik <- function(tree, residuals, alpha,
                                            trait.covariance, root.model,
                                            observation.error=NULL){
    residuals <- as.matrix(residuals)
    n <- length(tree$tip.label)
    p <- ncol(residuals)
    if(nrow(residuals) != n){
        stop("residuals must have one row per tree tip.")
    }
    if(!is.null(rownames(residuals))){
        residuals <- residuals[tree$tip.label, , drop=FALSE]
    }
    alpha <- rep(as.numeric(alpha), length.out=p)
    if(length(alpha) != p || any(!is.finite(alpha)) || any(alpha < 0)){
        stop("alpha must contain one non-negative value per trait.")
    }
    if(!all(dim(trait.covariance) == c(p, p))){
        stop("trait.covariance must be a square matrix with one row per trait.")
    }
    if(is.null(observation.error)){
        observation.error <- matrix(0, n, p)
    } else{
        observation.error <- matrix(
            as.numeric(observation.error), nrow=n, ncol=p
        )
    }
    if(any(!is.finite(observation.error)) || any(observation.error < 0)){
        stop("observation.error must contain finite non-negative variances.")
    }

    ordered <- ape::reorder.phylo(tree, "postorder")
    messages <- vector("list", n + ordered$Nnode)
    zero.message <- function(){
        list(precision=matrix(0, p, p), linear=numeric(p), constant=0)
    }
    transition.for <- function(length){
        diag(exp(-alpha * length), p)
    }
    for(edge.index in seq_len(nrow(ordered$edge))){
        parent <- ordered$edge[edge.index, 1L]
        child <- ordered$edge[edge.index, 2L]
        branch.length <- ordered$edge.length[[edge.index]]
        transition <- transition.for(branch.length)
        innovation <- ou_branch_innovation_covariance(
            alpha, trait.covariance, branch.length
        )
        incoming <- if(child <= n){
            tip_pruning_message(
                residuals[child, ], transition, innovation,
                observation.error[child, ]
            )
        } else{
            transform_pruning_message(messages[[child]], transition, innovation)
        }
        if(is.null(messages[[parent]])){
            messages[[parent]] <- zero.message()
        }
        messages[[parent]]$precision <- messages[[parent]]$precision +
            incoming$precision
        messages[[parent]]$linear <- messages[[parent]]$linear + incoming$linear
        messages[[parent]]$constant <- messages[[parent]]$constant +
            incoming$constant
    }
    root <- setdiff(unique(ordered$edge[, 1L]), ordered$edge[, 2L])[[1L]]
    root.message <- messages[[root]]
    if(identical(root.model, "OUfixedRoot")){
        return(root.message$constant)
    }
    stationary <- pruning_spd_components(
        ou_stationary_trait_covariance(alpha, trait.covariance)
    )
    combined <- pruning_spd_components(
        root.message$precision + stationary$inverse
    )
    root.message$constant - 0.5 * stationary$log.determinant -
        0.5 * combined$log.determinant +
        0.5 * drop(crossprod(
            root.message$linear,
            combined$inverse %*% root.message$linear
        ))
}

pruning_multivariate_ou_fit <- function(tree, Y, design.builder, opt,
                                        fixed.alpha=NULL){
    Y <- as.matrix(Y)
    n <- nrow(Y)
    p <- ncol(Y)
    observed <- !is.na(as.vector(Y))
    bounds <- multivariate_alpha_bounds(tree, opt, fixed.alpha, p)
    initial.alpha <- bounds$start
    initial.design <- design.builder(initial.alpha)
    design.columns <- vapply(initial.design, ncol, integer(1))
    if(length(unique(design.columns)) != 1L){
        stop("trait-specific mean design blocks must have equal dimensions.")
    }
    q <- design.columns[[1L]]

    unit.covariance <- multivariate_ou_dense_covariance(
        tree, initial.alpha, diag(p), opt$root.model
    )
    marginal.scale <- vapply(seq_len(p), function(j){
        rows <- ((j - 1L) * n + 1L):(j * n)
        max(mean(diag(unit.covariance[rows, rows, drop=FALSE])),
            .Machine$double.eps)
    }, numeric(1))
    marginal.variance <- vapply(seq_len(p), function(j){
        max(stats::var(Y[, j], na.rm=TRUE) / marginal.scale[[j]], 1e-8)
    }, numeric(1))
    initial.sd <- sqrt(marginal.variance)
    initial.coefficients <- matrix(0, q, p)
    for(j in seq_len(p)){
        keep <- !is.na(Y[, j])
        coefficient <- lm.fit(initial.design[[j]][keep, , drop=FALSE],
                              Y[keep, j])$coefficients
        coefficient[!is.finite(coefficient)] <- 0
        initial.coefficients[, j] <- coefficient
    }

    shrinkage.lambda <- 0
    if(identical(opt$covariance.regularization, "shrinkage")){
        shrinkage.lambda <- opt$regularization.lambda
        if(is.na(shrinkage.lambda)){
            shrinkage.lambda <- min(0.95, p / max(p + n, 1L))
        }
    }

    decode <- function(parameters){
        index <- 1L
        if(bounds$fixed){
            alpha <- bounds$start
        } else{
            alpha <- exp(parameters[index:(index + bounds$n.alpha - 1L)])
            index <- index + bounds$n.alpha
        }
        L <- diag(exp(parameters[index:(index + p - 1L)]), p)
        index <- index + p
        if(p > 1L){
            count <- p * (p - 1L) / 2L
            L[lower.tri(L)] <- parameters[index:(index + count - 1L)]
            index <- index + count
        }
        sigma2.error <- numeric(p)
        if(isTRUE(opt$measurement_error)){
            sigma2.error <- exp(parameters[index:(index + p - 1L)])
            index <- index + p
        }
        coefficients <- matrix(parameters[index:(index + q * p - 1L)], q, p)
        list(
            alpha=alpha,
            trait.covariance=tcrossprod(L),
            sigma2.error=sigma2.error,
            coefficients=coefficients
        )
    }

    evaluate <- function(parameters){
        decoded <- decode(parameters)
        design <- design.builder(decoded$alpha)
        fitted <- matrix(NA_real_, n, p, dimnames=dimnames(Y))
        for(j in seq_len(p)){
            fitted[, j] <- design[[j]] %*% decoded$coefficients[, j]
        }
        residuals <- Y - fitted
        observation.error <- matrix(
            multivariate_observation_error_vector(
                Y, opt, decoded$sigma2.error
            ), nrow=n, ncol=p
        )
        effective.root <- if(all(decoded$alpha <= .Machine$double.eps)){
            "OUfixedRoot"
        } else opt$root.model
        log.likelihood <- tryCatch(
            pruning_multivariate_ou_loglik(
                tree, residuals, decoded$alpha, decoded$trait.covariance,
                effective.root, observation.error
            ),
            error=function(e) NA_real_
        )
        penalty <- 0
        if(shrinkage.lambda > 0){
            correlation <- stats::cov2cor(decoded$trait.covariance)
            penalty <- n * shrinkage.lambda / max(1 - shrinkage.lambda, 1e-8) *
                sum(correlation[lower.tri(correlation)]^2)
        }
        c(decoded, list(
            base.design=design,
            fitted.values=fitted,
            residuals=residuals,
            logLik=log.likelihood,
            n2llh=-2 * log.likelihood,
            regularization.penalty=penalty,
            regularization.lambda=shrinkage.lambda
        ))
    }

    parameters <- lower <- upper <- numeric(0)
    if(!bounds$fixed){
        parameters <- c(parameters, log(initial.alpha))
        lower <- c(lower, log(bounds$lower))
        upper <- c(upper, log(bounds$upper))
    }
    parameters <- c(parameters, log(initial.sd))
    lower <- c(lower, log(initial.sd) - 12)
    upper <- c(upper, log(initial.sd) + 12)
    if(p > 1L){
        count <- p * (p - 1L) / 2L
        off.bound <- 100 * max(initial.sd)
        parameters <- c(parameters, rep(0, count))
        lower <- c(lower, rep(-off.bound, count))
        upper <- c(upper, rep(off.bound, count))
    }
    if(isTRUE(opt$measurement_error)){
        error.start <- pmax(vapply(seq_len(p), function(j){
            stats::var(Y[, j], na.rm=TRUE) * 0.05
        }, numeric(1)), 1e-10)
        parameters <- c(parameters, log(error.start))
        lower <- c(lower, log(error.start) - 16)
        upper <- c(upper, log(error.start) + 8)
    }
    process.parameter.count <- length(parameters)
    parameters <- c(parameters, as.vector(initial.coefficients))
    lower <- c(lower, rep(-Inf, q * p))
    upper <- c(upper, rep(Inf, q * p))

    objective <- function(parameters){
        fit <- evaluate(parameters)
        if(!is.finite(fit$n2llh)){
            return(.Machine$double.xmax / 1000)
        }
        fit$n2llh + fit$regularization.penalty
    }
    starting.points <- list(parameters)
    if(opt$optimizer.starts > 1L){
        for(start.index in 2:opt$optimizer.starts){
            candidate <- parameters
            fraction <- (start.index - 1L) / opt$optimizer.starts
            varied <- seq_len(process.parameter.count)
            phase <- (varied * 0.61803398875 + fraction) %% 1
            candidate[varied] <- lower[varied] + phase *
                (upper[varied] - lower[varied])
            starting.points[[start.index]] <- candidate
        }
    }
    optimization.runs <- lapply(seq_along(starting.points), function(index){
        result <- tryCatch(
            optim(
                starting.points[[index]], objective, method="L-BFGS-B",
                lower=lower, upper=upper,
                control=list(maxit=2000, factr=1e7)
            ),
            error=function(e) list(
                par=starting.points[[index]], value=Inf,
                convergence=99L, message=conditionMessage(e)
            )
        )
        result$start.index <- index
        result
    })
    values <- vapply(optimization.runs, function(x) x$value, numeric(1))
    if(!any(is.finite(values))){
        stop("all optimization starts failed for the pruning OU likelihood.")
    }
    optimized <- optimization.runs[[which.min(values)]]
    fit <- evaluate(optimized$par)
    if(!is.finite(fit$logLik)){
        stop("failed to optimize the pruning OU likelihood.")
    }

    joint.design <- multivariate_joint_design(
        fit$base.design, p, observed
    )
    fit$rank <- qr(joint.design, tol=1e-9, LAPACK=FALSE)$rank
    fit$mean.rank <- fit$rank
    fit$prediction.coefficients <- fit$coefficients
    fit$alpha.df <- if(bounds$fixed) 0L else bounds$n.alpha
    fit$covariance.df <- p * (p + 1L) / 2L
    fit$error.df <- if(isTRUE(opt$measurement_error)) p else 0L
    fit$p <- fit$rank + fit$covariance.df + fit$alpha.df + fit$error.df
    fit$n <- sum(observed)
    fit$sample.size <- sum(rowSums(!is.na(Y)) > 0L)
    fit$observed <- observed
    fit$convergence <- optimized$convergence
    at.bound <- is.finite(lower) & (
        abs(optimized$par - lower) <= 1e-6 * pmax(1, abs(lower)) |
        abs(optimized$par - upper) <= 1e-6 * pmax(1, abs(upper))
    )
    hessian <- NULL
    hessian.eigenvalues <- numeric(0)
    parameter.vcov <- NULL
    hessian.skipped <- FALSE
    if(isTRUE(opt$compute.hessian) && length(optimized$par) <= 100L){
        hessian <- tryCatch(optimHess(optimized$par, objective),
                            error=function(e) NULL)
        if(!is.null(hessian) && all(is.finite(hessian))){
            hessian <- 0.5 * (hessian + t(hessian))
            hessian.eigenvalues <- eigen(
                hessian, symmetric=TRUE, only.values=TRUE
            )$values
            if(all(hessian.eigenvalues > 1e-8)){
                parameter.vcov <- tryCatch(2 * solve(hessian),
                                           error=function(e) NULL)
            }
        }
    } else if(isTRUE(opt$compute.hessian)){
        hessian.skipped <- TRUE
    }
    covariance.eigenvalues <- eigen(
        fit$trait.covariance, symmetric=TRUE, only.values=TRUE
    )$values
    fit$optimization <- list(
        selected.start=optimized$start.index,
        runs=data.frame(
            start=seq_along(optimization.runs), objective=values,
            convergence=vapply(optimization.runs, function(x) x$convergence,
                               integer(1)), stringsAsFactors=FALSE
        ),
        message=optimized$message,
        at.bound=at.bound,
        hessian=hessian,
        hessian.eigenvalues=hessian.eigenvalues,
        parameter.vcov=parameter.vcov
    )
    fit$diagnostics <- list(
        converged=identical(as.integer(optimized$convergence), 0L),
        boundary=any(at.bound),
        trait.covariance.condition=kappa(fit$trait.covariance),
        trait.covariance.min.eigenvalue=min(covariance.eigenvalues),
        hessian.positive.definite=length(hessian.eigenvalues) > 0L &&
            all(hessian.eigenvalues > 1e-8),
        hessian.skipped=hessian.skipped,
        regularization.lambda=fit$regularization.lambda,
        likelihood.engine="pruning"
    )
    if(!fit$diagnostics$converged && !isTRUE(opt$quietly)){
        warning("the selected pruning-likelihood optimization did not converge.")
    }
    fit
}
