is_full_trait_covariance <- function(opt, Y=NULL){
    mode <- opt$trait.covariance
    if(is.null(mode)){
        mode <- "diagonal"
    }
    identical(mode, "full") && (is.null(Y) || ncol(as.matrix(Y)) > 1L)
}

validate_trait_covariance_mode <- function(mode, Y, criterion){
    mode <- match.arg(mode, c("diagonal", "full"))
    Y <- as.matrix(Y)
    if(mode == "full" && ncol(Y) < 2L){
        warning("trait.covariance = \"full\" has no effect for a univariate response; using \"diagonal\".")
        mode <- "diagonal"
    }
    if(mode == "full" && !identical(criterion, "BIC")){
        stop(paste0(
            "trait.covariance = \"full\" currently supports criterion = \"BIC\". ",
            "Small-sample and phylogenetic penalties have not been derived for ",
            "this correlated multivariate OU model."
        ))
    }
    mode
}

normalize_multivariate_options <- function(opt, Y){
    Y <- as.matrix(Y)
    if(is.null(opt$alpha.structure)){
        opt$alpha.structure <- "shared"
    }
    opt$alpha.structure <- match.arg(
        opt$alpha.structure, c("shared", "diagonal")
    )
    if(is.null(opt$covariance.regularization)){
        opt$covariance.regularization <- "none"
    }
    opt$covariance.regularization <- match.arg(
        opt$covariance.regularization, c("none", "shrinkage")
    )
    if(is.null(opt$regularization.lambda)){
        opt$regularization.lambda <- NA_real_
    }
    if(length(opt$regularization.lambda) != 1L ||
       (!is.na(opt$regularization.lambda) &&
        (!is.finite(opt$regularization.lambda) ||
         opt$regularization.lambda < 0 || opt$regularization.lambda > 1))){
        stop("regularization.lambda must be NA or a number between 0 and 1.")
    }
    if(is.null(opt$likelihood.engine)){
        opt$likelihood.engine <- "auto"
    }
    opt$likelihood.engine <- match.arg(
        opt$likelihood.engine, c("auto", "dense", "pruning")
    )
    if(is.null(opt$optimizer.starts)){
        opt$optimizer.starts <- 5L
    }
    opt$optimizer.starts <- l1ou_integer_argument(
        opt$optimizer.starts, "optimizer.starts", 1L
    )
    if(is.null(opt$compute.hessian)){
        opt$compute.hessian <- TRUE
    }
    opt$compute.hessian <- l1ou_logical_argument(
        opt$compute.hessian, "compute.hessian"
    )
    if(is.null(opt$measurement_error)) opt$measurement_error <- FALSE
    opt$measurement_error <- l1ou_logical_argument(
        opt$measurement_error, "measurement_error"
    )

    if(!is_full_trait_covariance(opt, Y)){
        opt$alpha.structure <- "diagonal"
        opt$covariance.regularization <- "none"
        opt$likelihood.engine <- "dense"
    }
    opt
}

sanitize_multivariate_alpha_arguments <- function(lower, upper, starting,
                                                   n.traits){
    lower <- rep(as.numeric(lower), length.out=n.traits)
    upper <- rep(as.numeric(upper), length.out=n.traits)
    starting <- rep(as.numeric(starting), length.out=n.traits)
    if(any(!is.na(lower) & lower < 0)){
        warning("alpha.lower values must be non-negative; negative values were set to zero.")
        lower[!is.na(lower) & lower < 0] <- 0
    }
    if(any(!is.na(upper) & upper < 0)){
        stop("alpha.upper values must be non-negative.")
    }
    zero.upper <- !is.na(upper) & upper == 0
    if(any(zero.upper & (is.na(lower) | lower != 0))){
        stop("alpha.upper can be zero only when alpha is fixed at zero.")
    }
    inverted <- !is.na(lower) & !is.na(upper) & upper < lower
    if(any(inverted)){
        warning("alpha.upper values below alpha.lower were set equal to alpha.lower.")
        upper[inverted] <- lower[inverted]
    }
    finite.bounds <- !is.na(lower) & !is.na(upper)
    if(any(zero.upper) &&
       !all(finite.bounds & abs(lower - upper) <= .Machine$double.eps^0.5)){
        stop("zero alpha bounds are supported only when every alpha is fixed.")
    }
    list(lower=lower, upper=upper, starting=starting)
}

regularize_positive_definite <- function(x, relative.floor=1e-10){
    x <- 0.5 * (x + t(x))
    eig <- eigen(x, symmetric=TRUE)
    largest <- max(eig$values, .Machine$double.eps)
    eig$values <- pmax(eig$values, largest * relative.floor)
    out <- eig$vectors %*% (eig$values * t(eig$vectors))
    0.5 * (out + t(out))
}

shrink_trait_covariance <- function(sample.covariance, residuals,
                                    lambda=NA_real_){
    sample.covariance <- 0.5 * (sample.covariance + t(sample.covariance))
    target <- diag(diag(sample.covariance), nrow=nrow(sample.covariance))
    if(is.na(lambda)){
        n <- nrow(residuals)
        variability <- 0
        for(i in seq_len(n)){
            outer <- tcrossprod(residuals[i, ])
            variability <- variability + sum((outer - sample.covariance)^2)
        }
        variability <- variability / n
        separation <- sum((sample.covariance - target)^2)
        lambda <- if(separation <= .Machine$double.eps){
            1
        } else{
            min(1, max(0, variability / (n * separation)))
        }
    }
    covariance <- (1 - lambda) * sample.covariance + lambda * target
    list(
        covariance=regularize_positive_definite(covariance),
        lambda=lambda
    )
}

multivariate_ou_base_covariance <- function(tree, alpha, root.model){
    effective.root <- if(alpha <= 0) "OUfixedRoot" else root.model
    re <- sqrt_OU_covariance(
        tree,
        alpha=alpha,
        root.model=effective.root,
        sigma2=1,
        check.order=FALSE,
        check.ultrametric=FALSE
    )
    tcrossprod(re$sqrtSigma)
}

multivariate_ou_dense_covariance <- function(tree, alpha, trait.covariance,
                                              root.model){
    alpha <- as.numeric(alpha)
    p <- nrow(trait.covariance)
    if(length(alpha) == 1L){
        alpha <- rep(alpha, p)
    }
    if(length(alpha) != p || any(!is.finite(alpha)) || any(alpha < 0)){
        stop("alpha must contain one non-negative value per trait.")
    }
    if(length(unique(alpha)) == 1L){
        return(kronecker(
            trait.covariance,
            multivariate_ou_base_covariance(tree, alpha[[1L]], root.model)
        ))
    }

    shared.time <- ape::vcv.phylo(tree)
    shared.time <- shared.time[tree$tip.label, tree$tip.label, drop=FALSE]
    tip.time <- diag(shared.time)
    n <- length(tip.time)
    covariance <- matrix(0, nrow=n * p, ncol=n * p)
    random.root <- identical(root.model, "OUrandomRoot")

    for(i in seq_len(p)){
        for(j in seq_len(p)){
            rate <- alpha[[i]] + alpha[[j]]
            block <- matrix(0, nrow=n, ncol=n)
            if(rate <= .Machine$double.eps){
                if(random.root){
                    stop("OUrandomRoot requires strictly positive alpha values.")
                }
                block <- trait.covariance[i, j] * shared.time
            } else{
                row.distance <- outer(tip.time, rep(1, n)) - shared.time
                col.distance <- outer(rep(1, n), tip.time) - shared.time
                propagation <- exp(
                    -alpha[[i]] * row.distance - alpha[[j]] * col.distance
                )
                ancestral.variance <- if(random.root){
                    matrix(1, n, n)
                } else{
                    -expm1(-rate * shared.time)
                }
                block <- trait.covariance[i, j] / rate * propagation *
                    ancestral.variance
            }
            rows <- ((i - 1L) * n + 1L):(i * n)
            cols <- ((j - 1L) * n + 1L):(j * n)
            covariance[rows, cols] <- block
        }
    }
    0.5 * (covariance + t(covariance))
}

multivariate_tip_trait_covariance <- function(tree, alpha, trait.covariance,
                                               root.model){
    p <- nrow(trait.covariance)
    alpha <- rep(as.numeric(alpha), length.out=p)
    rates <- outer(alpha, alpha, "+")
    multiplier <- matrix(0, p, p)
    positive <- rates > .Machine$double.eps
    if(identical(root.model, "OUrandomRoot")){
        if(any(!positive)){
            stop("OUrandomRoot requires strictly positive alpha values.")
        }
        multiplier <- 1 / rates
    } else{
        height <- mean(ape::node.depth.edgelength(tree)[
            seq_along(tree$tip.label)
        ])
        multiplier[positive] <- -expm1(-rates[positive] * height) /
            rates[positive]
        multiplier[!positive] <- height
    }
    covariance <- trait.covariance * multiplier
    0.5 * (covariance + t(covariance))
}

multivariate_marginal_scale <- function(tree, alpha, root.model, n.traits){
    alpha <- rep(as.numeric(alpha), length.out=n.traits)
    diag(multivariate_tip_trait_covariance(
        tree, alpha, diag(n.traits), root.model
    ))
}

initial_trait_covariance <- function(Y, initial.sd){
    p <- ncol(Y)
    correlation <- suppressWarnings(stats::cor(
        Y, use="pairwise.complete.obs"
    ))
    if(p == 1L) correlation <- matrix(1, 1L, 1L)
    correlation[!is.finite(correlation)] <- 0
    diag(correlation) <- 1
    correlation <- regularize_positive_definite(correlation, 1e-6)
    covariance <- correlation * tcrossprod(initial.sd)
    regularize_positive_definite(covariance, 1e-8)
}

multivariate_observation_error_vector <- function(Y, opt, sigma2.error=NULL){
    Y <- as.matrix(Y)
    n <- nrow(Y)
    p <- ncol(Y)
    result <- numeric(n * p)

    if(!is.null(opt$input_error)){
        input.error <- as.matrix(opt$input_error)
        if(ncol(input.error) == 1L && p > 1L){
            input.error <- input.error[, rep(1L, p), drop=FALSE]
        }
        result <- result + as.vector(input.error)
    }

    if(!is.null(sigma2.error)){
        result <- result + rep(as.numeric(sigma2.error), each=n)
    }
    result
}

multivariate_ou_observed_covariance <- function(tree, Y, alpha,
                                                 trait.covariance, opt,
                                                 sigma2.error=NULL){
    Y <- as.matrix(Y)
    covariance <- multivariate_ou_dense_covariance(
        tree, alpha, trait.covariance, opt$root.model
    )
    observation.error <- multivariate_observation_error_vector(
        Y, opt, sigma2.error
    )
    diag(covariance) <- diag(covariance) + observation.error
    observed <- !is.na(as.vector(Y))
    covariance[observed, observed, drop=FALSE]
}

multivariate_joint_design <- function(base.design, n.traits, observed){
    blocks <- if(is.list(base.design)){
        base.design
    } else{
        rep(list(as.matrix(base.design)), n.traits)
    }
    if(length(blocks) != n.traits){
        stop("the multivariate mean design must contain one block per trait.")
    }
    full <- block_diag_matrix(blocks)
    full[observed, , drop=FALSE]
}

profile_multivariate_gls <- function(y, design, covariance){
    covariance <- 0.5 * (covariance + t(covariance))
    covariance.chol <- tryCatch(chol(covariance), error=function(e) NULL)
    if(is.null(covariance.chol)){
        return(NULL)
    }

    whitened.y <- tryCatch(
        drop(forwardsolve(t(covariance.chol), matrix(y, ncol=1L))),
        error=function(e) NULL
    )
    whitened.design <- tryCatch(
        forwardsolve(t(covariance.chol), design),
        error=function(e) NULL
    )
    if(is.null(whitened.y) || is.null(whitened.design)){
        return(NULL)
    }

    qr.design <- qr(whitened.design, tol=1e-9, LAPACK=FALSE)
    rank <- qr.design$rank
    if(rank < 1L){
        return(NULL)
    }
    keep <- sort(qr.design$pivot[seq_len(rank)])
    reduced <- whitened.design[, keep, drop=FALSE]
    cross <- crossprod(reduced)
    cross.chol <- tryCatch(chol(cross), error=function(e) NULL)
    if(is.null(cross.chol)){
        return(NULL)
    }
    coefficients.reduced <- drop(backsolve(
        cross.chol,
        forwardsolve(t(cross.chol), crossprod(reduced, whitened.y))
    ))
    coefficients <- rep(NA_real_, ncol(design))
    coefficients[keep] <- coefficients.reduced
    coefficients.for.prediction <- coefficients
    coefficients.for.prediction[is.na(coefficients.for.prediction)] <- 0
    fitted <- drop(design %*% coefficients.for.prediction)
    residual <- y - fitted
    whitened.residual <- drop(forwardsolve(
        t(covariance.chol), matrix(residual, ncol=1L)
    ))
    log.determinant <- 2 * sum(log(diag(covariance.chol)))

    list(
        coefficients=coefficients,
        coefficients.for.prediction=coefficients.for.prediction,
        fitted.values=fitted,
        residuals=residual,
        rank=rank,
        keep=keep,
        covariance.chol=covariance.chol,
        n2llh=length(y) * log(2 * pi) + log.determinant +
            sum(whitened.residual^2)
    )
}

multivariate_alpha_bounds <- function(tree, opt, fixed.alpha=NULL, n.traits=1L){
    n.alpha <- if(identical(opt$alpha.structure, "diagonal")) n.traits else 1L
    if(!is.null(fixed.alpha)){
        fixed.alpha <- rep(as.numeric(fixed.alpha), length.out=n.alpha)
        return(list(fixed=TRUE, lower=fixed.alpha, upper=fixed.alpha,
                    start=fixed.alpha, n.alpha=n.alpha))
    }
    tree.height <- mean(ape::node.depth.edgelength(tree)[seq_along(tree$tip.label)])
    lower <- rep(as.numeric(opt$alpha.lower.bound), length.out=n.alpha)
    lower[is.na(lower) | lower <= 0] <- 1e-7 / tree.height
    upper <- rep(as.numeric(opt$alpha.upper.bound), length.out=n.alpha)
    upper[is.na(upper)] <- alpha_upper_bound(tree)
    start <- rep(as.numeric(opt$alpha.starting.value), length.out=n.alpha)
    missing.start <- is.na(start) | start <= 0
    start[missing.start] <- 0.5 / tree.height
    start <- pmin(pmax(start, lower), upper)
    list(
        fixed=all(abs(lower - upper) <= .Machine$double.eps^0.5),
        lower=lower,
        upper=upper,
        start=start,
        n.alpha=n.alpha
    )
}

matrix_normal_ou_fit <- function(tree, Y, design.builder, opt,
                                 fixed.alpha=NULL){
    Y <- as.matrix(Y)
    n <- nrow(Y)
    p <- ncol(Y)
    bounds <- multivariate_alpha_bounds(tree, opt, fixed.alpha, p)
    if(bounds$n.alpha != 1L){
        stop("matrix-normal profiling requires alpha.structure = \"shared\".")
    }

    profile.alpha <- function(alpha){
        base.covariance <- multivariate_ou_base_covariance(
            tree, alpha, opt$root.model
        )
        covariance.chol <- tryCatch(chol(base.covariance), error=function(e) NULL)
        if(is.null(covariance.chol)){
            return(NULL)
        }
        built.design <- design.builder(alpha)
        design <- as.matrix(if(is.list(built.design)) built.design[[1L]] else built.design)
        whitened.design <- forwardsolve(t(covariance.chol), design)
        whitened.Y <- forwardsolve(t(covariance.chol), Y)
        qr.design <- qr(whitened.design, tol=1e-9, LAPACK=FALSE)
        rank <- qr.design$rank
        if(rank < 1L ||
           (n - rank < p && identical(opt$covariance.regularization, "none"))){
            return(NULL)
        }
        keep <- sort(qr.design$pivot[seq_len(rank)])
        reduced <- whitened.design[, keep, drop=FALSE]
        cross.chol <- tryCatch(chol(crossprod(reduced)), error=function(e) NULL)
        if(is.null(cross.chol)){
            return(NULL)
        }
        coefficients.reduced <- backsolve(
            cross.chol,
            forwardsolve(t(cross.chol), crossprod(reduced, whitened.Y))
        )
        coefficients <- matrix(NA_real_, nrow=ncol(design), ncol=p)
        coefficients[keep, ] <- coefficients.reduced
        prediction.coefficients <- coefficients
        prediction.coefficients[is.na(prediction.coefficients)] <- 0
        fitted <- design %*% prediction.coefficients
        residual <- Y - fitted
        whitened.residual <- forwardsolve(t(covariance.chol), residual)
        raw.covariance <- crossprod(whitened.residual) / n
        shrinkage <- list(
            covariance=regularize_positive_definite(raw.covariance),
            lambda=0
        )
        if(identical(opt$covariance.regularization, "shrinkage")){
            shrinkage <- shrink_trait_covariance(
                raw.covariance,
                whitened.residual,
                opt$regularization.lambda
            )
        }
        trait.covariance <- shrinkage$covariance
        trait.chol <- tryCatch(chol(trait.covariance), error=function(e) NULL)
        if(is.null(trait.chol)){
            return(NULL)
        }
        trait.whitened.residual <- forwardsolve(
            t(trait.chol), t(whitened.residual)
        )
        n2llh <- n * p * log(2 * pi) +
            p * 2 * sum(log(diag(covariance.chol))) +
            n * 2 * sum(log(diag(trait.chol))) +
            sum(trait.whitened.residual^2)
        list(
            n2llh=n2llh,
            alpha=alpha,
            coefficients=coefficients,
            prediction.coefficients=prediction.coefficients,
            fitted.values=fitted,
            residuals=residual,
            trait.covariance=trait.covariance,
            regularization.lambda=shrinkage$lambda,
            sigma2.error=rep(0, p),
            mean.rank=rank * p,
            base.design=design
        )
    }

    optimized <- NULL
    if(bounds$fixed){
        fit <- profile.alpha(bounds$start)
    } else{
        optimized <- optimize(
            function(log.alpha){
                fit <- profile.alpha(exp(log.alpha))
                if(is.null(fit)) .Machine$double.xmax / 1000 else fit$n2llh
            },
            interval=log(c(bounds$lower, bounds$upper)),
            tol=1e-6
        )
        fit <- profile.alpha(exp(optimized$minimum))
    }
    if(is.null(fit)){
        stop(paste0(
            "the full trait-covariance model is not estimable: require at least ",
            "nrow(Y) - rank(mean design) >= ncol(Y) complete residual contrasts."
        ))
    }
    fit$alpha.df <- ifelse(bounds$fixed, 0L, 1L)
    fit$covariance.df <- p * (p + 1L) / 2L
    fit$p <- fit$mean.rank + fit$covariance.df + fit$alpha.df
    fit$n <- n * p
    fit$sample.size <- n
    fit$logLik <- -fit$n2llh / 2
    fit$observed <- rep(TRUE, n * p)
    covariance.eigenvalues <- eigen(
        fit$trait.covariance, symmetric=TRUE, only.values=TRUE
    )$values
    at.bound <- !bounds$fixed && (
        abs(fit$alpha - bounds$lower) <= 1e-6 * max(1, abs(bounds$lower)) ||
        abs(fit$alpha - bounds$upper) <= 1e-6 * max(1, abs(bounds$upper))
    )
    alpha.hessian <- NA_real_
    if(!bounds$fixed && isTRUE(opt$compute.hessian)){
        step <- 1e-4
        center <- log(fit$alpha)
        values <- vapply(c(center - step, center, center + step), function(x){
            profiled <- profile.alpha(exp(x))
            if(is.null(profiled)) NA_real_ else profiled$n2llh
        }, numeric(1))
        if(all(is.finite(values))){
            alpha.hessian <- (values[[1L]] - 2 * values[[2L]] + values[[3L]]) /
                step^2
        }
    }
    fit$optimization <- list(
        selected.start=1L,
        runs=data.frame(start=1L, objective=fit$n2llh, convergence=0L),
        message=NULL,
        at.bound=at.bound,
        hessian=alpha.hessian,
        hessian.eigenvalues=alpha.hessian,
        parameter.vcov=if(is.finite(alpha.hessian) && alpha.hessian > 0 &&
                          identical(opt$covariance.regularization, "none"))
            matrix(2 / alpha.hessian, 1L, 1L) else NULL,
        alpha.parameter.index=if(bounds$fixed) integer() else 1L
    )
    fit$diagnostics <- list(
        converged=TRUE,
        boundary=isTRUE(at.bound),
        trait.covariance.condition=kappa(fit$trait.covariance),
        trait.covariance.min.eigenvalue=min(covariance.eigenvalues),
        hessian.positive.definite=is.finite(alpha.hessian) && alpha.hessian > 0,
        hessian.type=if(identical(opt$covariance.regularization, "none"))
            "likelihood" else "regularized-profile",
        regularization.lambda=fit$regularization.lambda,
        information.criterion.calibrated=identical(
            opt$covariance.regularization, "none"
        ),
        score.note=if(identical(opt$covariance.regularization, "none")) NULL else
            paste0(
                "BIC is evaluated at a regularized estimate and should be ",
                "treated as a sensitivity score, not a calibrated marginal-",
                "likelihood approximation."
            ),
        likelihood.engine="matrix-normal"
    )
    fit
}

general_multivariate_ou_fit <- function(tree, Y, design.builder, opt,
                                        fixed.alpha=NULL){
    Y <- as.matrix(Y)
    n <- nrow(Y)
    p <- ncol(Y)
    observed <- !is.na(as.vector(Y))
    y <- as.vector(Y)[observed]
    if(any(colSums(!is.na(Y)) < 2L)){
        stop("each trait needs at least two observed tips for full covariance estimation.")
    }
    bounds <- multivariate_alpha_bounds(tree, opt, fixed.alpha, p)
    shrinkage.lambda <- 0
    if(identical(opt$covariance.regularization, "shrinkage")){
        shrinkage.lambda <- opt$regularization.lambda
        if(is.na(shrinkage.lambda)){
            shrinkage.lambda <- min(0.95, p / max(p + n, 1L))
        }
    }
    initial.alpha <- bounds$start
    marginal.scale <- pmax(
        multivariate_marginal_scale(tree, initial.alpha, opt$root.model, p),
        .Machine$double.eps
    )
    marginal.variance <- vapply(seq_len(p), function(j){
        value <- stats::var(Y[, j], na.rm=TRUE) / marginal.scale[[j]]
        max(value, 1e-8)
    }, numeric(1))
    initial.sd <- sqrt(marginal.variance)
    initial.covariance <- initial_trait_covariance(Y, initial.sd)
    initial.L <- t(chol(initial.covariance))

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
        sigma2.error <- rep(0, p)
        if(isTRUE(opt$measurement_error)){
            sigma2.error <- exp(parameters[index:(index + p - 1L)])
        }
        list(
            alpha=alpha,
            trait.covariance=tcrossprod(L),
            sigma2.error=sigma2.error
        )
    }

    evaluate <- function(parameters){
        decoded <- decode(parameters)
        covariance <- multivariate_ou_observed_covariance(
            tree, Y, decoded$alpha, decoded$trait.covariance, opt,
            decoded$sigma2.error
        )
        base.design <- design.builder(decoded$alpha)
        design <- multivariate_joint_design(base.design, p, observed)
        profiled <- profile_multivariate_gls(y, design, covariance)
        if(is.null(profiled)){
            return(NULL)
        }
        penalty <- 0
        if(shrinkage.lambda > 0){
            correlation <- stats::cov2cor(decoded$trait.covariance)
            penalty <- n * shrinkage.lambda / max(1 - shrinkage.lambda, 1e-8) *
                sum(correlation[lower.tri(correlation)]^2)
        }
        c(profiled, decoded, list(
            base.design=base.design,
            regularization.penalty=penalty,
            regularization.lambda=shrinkage.lambda
        ))
    }

    parameters <- numeric(0)
    lower <- numeric(0)
    upper <- numeric(0)
    if(!bounds$fixed){
        parameters <- c(parameters, log(initial.alpha))
        lower <- c(lower, log(bounds$lower))
        upper <- c(upper, log(bounds$upper))
    }
    parameters <- c(parameters, log(diag(initial.L)))
    lower <- c(lower, log(initial.sd) - 12)
    upper <- c(upper, log(initial.sd) + 12)
    if(p > 1L){
        count <- p * (p - 1L) / 2L
        off.bound <- 10 * max(initial.sd)
        parameters <- c(parameters, initial.L[lower.tri(initial.L)])
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

    objective <- function(parameters){
        fit <- evaluate(parameters)
        if(is.null(fit) || !is.finite(fit$n2llh)){
            return(.Machine$double.xmax / 1000)
        }
        fit$n2llh + fit$regularization.penalty
    }
    starting.points <- list(parameters)
    if(opt$optimizer.starts > 1L){
        for(start.index in 2:opt$optimizer.starts){
            fraction <- (start.index - 1L) / opt$optimizer.starts
            candidate <- parameters
            finite.bounds <- is.finite(lower) & is.finite(upper)
            phase <- (seq_along(candidate) * 0.61803398875 + fraction) %% 1
            local.span <- pmin(upper - lower, 4 * pmax(abs(parameters), 1))
            candidate[finite.bounds] <- parameters[finite.bounds] +
                (phase[finite.bounds] - 0.5) * local.span[finite.bounds]
            candidate[!finite.bounds] <- candidate[!finite.bounds] +
                (-1)^start.index * 0.2 * start.index
            starting.points[[start.index]] <- pmin(pmax(candidate, lower), upper)
        }
    }
    optimization.runs <- lapply(seq_along(starting.points), function(index){
        result <- tryCatch(
            optim(
                starting.points[[index]],
                objective,
                method="L-BFGS-B",
                lower=lower,
                upper=upper,
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
        stop("all optimization starts failed for the full trait-covariance model.")
    }
    optimized <- optimization.runs[[which.min(values)]]
    fit <- evaluate(optimized$par)
    if(is.null(fit) || !is.finite(fit$n2llh)){
        stop("failed to optimize the full trait-covariance OU likelihood.")
    }

    design.columns <- vapply(fit$base.design, ncol, integer(1))
    if(length(unique(design.columns)) != 1L){
        stop("trait-specific mean design blocks must have the same number of columns.")
    }
    coefficients <- matrix(
        fit$coefficients,
        nrow=design.columns[[1L]], ncol=p
    )
    prediction.coefficients <- matrix(
        fit$coefficients.for.prediction,
        nrow=design.columns[[1L]], ncol=p
    )
    fitted.vector <- residual.vector <- rep(NA_real_, n * p)
    fitted.vector[observed] <- fit$fitted.values
    residual.vector[observed] <- fit$residuals
    fit$coefficients <- coefficients
    fit$prediction.coefficients <- prediction.coefficients
    fit$fitted.values <- matrix(fitted.vector, nrow=n, ncol=p)
    fit$residuals <- matrix(residual.vector, nrow=n, ncol=p)
    fit$alpha.df <- ifelse(bounds$fixed, 0L, bounds$n.alpha)
    fit$covariance.df <- p * (p + 1L) / 2L
    fit$error.df <- ifelse(isTRUE(opt$measurement_error), p, 0L)
    fit$p <- fit$rank + fit$covariance.df + fit$alpha.df + fit$error.df
    fit$n <- length(y)
    fit$sample.size <- sum(rowSums(!is.na(Y)) > 0L)
    fit$logLik <- -fit$n2llh / 2
    fit$observed <- observed
    fit$convergence <- optimized$convergence
    at.bound <- abs(optimized$par - lower) <= 1e-6 * pmax(1, abs(lower)) |
        abs(optimized$par - upper) <= 1e-6 * pmax(1, abs(upper))
    hessian <- NULL
    hessian.eigenvalues <- numeric(0)
    parameter.vcov <- NULL
    if(isTRUE(opt$compute.hessian)){
        hessian <- tryCatch(
            optimHess(optimized$par, objective),
            error=function(e) NULL
        )
        if(!is.null(hessian) && all(is.finite(hessian))){
            hessian <- 0.5 * (hessian + t(hessian))
            hessian.eigenvalues <- eigen(hessian, symmetric=TRUE,
                                         only.values=TRUE)$values
            if(all(hessian.eigenvalues > 1e-8)){
                if(identical(opt$covariance.regularization, "none")){
                    parameter.vcov <- tryCatch(
                        2 * solve(hessian), error=function(e) NULL
                    )
                }
            }
        }
    }
    covariance.eigenvalues <- eigen(
        fit$trait.covariance, symmetric=TRUE, only.values=TRUE
    )$values
    fit$optimization <- list(
        selected.start=optimized$start.index,
        runs=data.frame(
            start=seq_along(optimization.runs),
            objective=values,
            convergence=vapply(optimization.runs, function(x) x$convergence,
                               integer(1)),
            message=vapply(optimization.runs, function(x){
                if(is.null(x$message)) "" else as.character(x$message[[1L]])
            }, character(1)),
            stringsAsFactors=FALSE
        ),
        message=optimized$message,
        at.bound=at.bound,
        hessian=hessian,
        hessian.eigenvalues=hessian.eigenvalues,
        parameter.vcov=parameter.vcov,
        alpha.parameter.index=if(bounds$fixed) integer() else
            seq_len(bounds$n.alpha)
    )
    fit$diagnostics <- list(
        converged=identical(as.integer(optimized$convergence), 0L),
        boundary=any(at.bound),
        trait.covariance.condition=kappa(fit$trait.covariance),
        trait.covariance.min.eigenvalue=min(covariance.eigenvalues),
        hessian.positive.definite=length(hessian.eigenvalues) > 0L &&
            all(hessian.eigenvalues > 1e-8),
        hessian.type=if(shrinkage.lambda > 0) "penalized" else "likelihood",
        regularization.lambda=fit$regularization.lambda,
        likelihood.engine="dense"
    )
    if(!fit$diagnostics$converged && !isTRUE(opt$quietly)){
        warning("the selected multivariate optimization run did not converge.")
    }
    fit
}

fit_multivariate_ou_likelihood <- function(tree, Y, shift.configuration, opt,
                                           regimes=NULL, fixed.alpha=NULL){
    Y <- as.matrix(Y)
    if(is.null(colnames(Y))){
        colnames(Y) <- paste0("trait", seq_len(ncol(Y)))
    }
    shift.configuration <- as.integer(shift.configuration)
    design.builder <- if(is.null(regimes)){
        function(alpha){
            rep(list(cbind(1, opt$Z[, shift.configuration, drop=FALSE])),
                ncol(Y))
        }
    } else{
        normalized <- normalize_convergent_regimes(regimes, shift.configuration)
        function(alpha){
            alpha <- rep(as.numeric(alpha), length.out=ncol(Y))
            lapply(seq_len(ncol(Y)), function(j){
                generate_prediction_vec(
                    tree,
                    shift.configuration,
                    normalized,
                    alpha[[j]],
                    root.model=opt$root.model
                )
            })
        }
    }

    initial.alpha <- multivariate_alpha_bounds(
        tree, opt, fixed.alpha, ncol(Y)
    )$start
    initial.design.blocks <- design.builder(initial.alpha)
    initial.design <- as.matrix(initial.design.blocks[[1L]])
    mean.rank <- qr(initial.design, tol=1e-9, LAPACK=FALSE)$rank
    if(nrow(Y) - mean.rank < ncol(Y) &&
       identical(opt$covariance.regularization, "none")){
        stop(paste0(
            "the full trait-covariance model is not estimable: require ",
            "nrow(Y) - rank(mean design) >= ncol(Y)."
        ))
    }

    complete.separable <- identical(opt$alpha.structure, "shared") &&
        !anyNA(Y) && !isTRUE(opt$measurement_error) && is.null(opt$input_error)
    use.matrix.normal <- identical(opt$likelihood.engine, "auto") &&
        complete.separable &&
        identical(opt$covariance.regularization, "none")
    use.pruning <- identical(opt$likelihood.engine, "pruning") ||
        (identical(opt$likelihood.engine, "auto") && !use.matrix.normal &&
         length(Y) >= 250L)
    fit <- if(use.matrix.normal){
        matrix_normal_ou_fit(tree, Y, design.builder, opt, fixed.alpha)
    } else if(use.pruning){
        pruning_multivariate_ou_fit(
            tree, Y, design.builder, opt, fixed.alpha
        )
    } else{
        general_multivariate_ou_fit(tree, Y, design.builder, opt, fixed.alpha)
    }
    final.design <- design.builder(fit$alpha)
    rownames(fit$coefficients) <- colnames(final.design[[1L]])
    colnames(fit$coefficients) <- colnames(Y)
    rownames(fit$trait.covariance) <- colnames(Y)
    colnames(fit$trait.covariance) <- colnames(Y)
    fit$trait.correlation <- stats::cov2cor(fit$trait.covariance)
    fit$tip.trait.covariance <- multivariate_tip_trait_covariance(
        tree, fit$alpha, fit$trait.covariance, opt$root.model
    )
    rownames(fit$tip.trait.covariance) <- colnames(Y)
    colnames(fit$tip.trait.covariance) <- colnames(Y)
    fit$tip.trait.correlation <- stats::cov2cor(fit$tip.trait.covariance)
    fit
}

multivariate_full_information_score <- function(fit, n.shifts, criterion){
    criterion <- match.arg(criterion, "BIC")
    parameters <- fit$p + n.shifts
    -2 * fit$logLik + log(fit$sample.size) * parameters
}

fit_full_covariance_l1ou_model <- function(tree, Y, shift.configuration, opt){
    Y <- as.matrix(Y)
    shift.configuration <- as.integer(shift.configuration)
    fit <- fit_multivariate_ou_likelihood(
        tree, Y, shift.configuration, opt
    )
    score <- multivariate_full_information_score(
        fit, length(shift.configuration), opt$criterion
    )
    n <- nrow(Y)
    p <- ncol(Y)
    trait.names <- colnames(Y)
    if(is.null(trait.names)){
        trait.names <- paste0("trait", seq_len(p))
        colnames(Y) <- trait.names
    }
    mu <- fit$fitted.values
    residuals <- fit$residuals
    dimnames(mu) <- dimnames(Y)
    dimnames(residuals) <- dimnames(Y)

    intercept <- fit$coefficients[1L, ]
    names(intercept) <- trait.names
    shift.means <- numeric(0)
    shift.values <- numeric(0)
    optima <- matrix(
        rep(intercept, each=n), nrow=n, ncol=p,
        dimnames=dimnames(Y)
    )
    if(length(shift.configuration) > 0L){
        shift.means <- fit$coefficients[-1L, , drop=FALSE]
        rownames(shift.means) <- as.character(shift.configuration)
        colnames(shift.means) <- trait.names
        alpha.by.trait <- rep(as.numeric(fit$alpha), length.out=p)
        scales <- vapply(alpha.by.trait, function(a){
            edge_scaling_from_cache(
                opt$edge.age, type="orgX", alpha=a
            )[shift.configuration]
        }, numeric(length(shift.configuration)))
        shift.values <- shift.means / scales
        optima <- optima + opt$Z[, shift.configuration, drop=FALSE] %*%
            shift.values
    }

    model.opt <- opt
    model.opt$prepared.tree <- NULL
    model.opt$prepared.tree.list <- NULL
    alpha <- rep(as.numeric(fit$alpha), length.out=p)
    names(alpha) <- trait.names
    sigma2 <- diag(fit$trait.covariance)
    names(sigma2) <- trait.names
    sigma2.error <- fit$sigma2.error
    names(sigma2.error) <- trait.names
    model <- list(
        Y=Y,
        tree=tree,
        tree.scale=if(is.null(opt$tree.scale)) 1 else opt$tree.scale,
        shift.configuration=shift.configuration,
        shift.values=shift.values,
        shift.means=shift.means,
        nShifts=length(shift.configuration),
        optima=optima,
        alpha=alpha,
        sigma2=sigma2,
        sigma2_error=sigma2.error,
        trait.covariance=fit$trait.covariance,
        trait.correlation=fit$trait.correlation,
        tip.trait.covariance=fit$tip.trait.covariance,
        tip.trait.correlation=fit$tip.trait.correlation,
        covariance.regularization=opt$covariance.regularization,
        regularization.lambda=fit$regularization.lambda,
        information.criterion.calibrated=identical(
            opt$covariance.regularization, "none"
        ),
        score.note=if(identical(opt$covariance.regularization, "none")) NULL else
            paste0(
                "BIC is evaluated at a regularized estimate and should be ",
                "treated as a sensitivity score."
            ),
        optimization=fit$optimization,
        diagnostics=fit$diagnostics,
        likelihood.engine=fit$diagnostics$likelihood.engine,
        parameter.count=fit$p,
        information.parameter.count=fit$p + length(shift.configuration),
        nobs=fit$sample.size,
        observed.entries=fit$n,
        intercept=intercept,
        mu=mu,
        residuals=residuals,
        score=score,
        logLik=fit$logLik,
        joint.logLik=fit$logLik,
        l1ou.options=model.opt
    )
    class(model) <- "l1ou"
    model
}

fit_full_covariance_convergent_model <- function(tree, Y, shift.configuration,
                                                  regimes, opt, score=NULL,
                                                  base.model=NULL){
    Y <- as.matrix(Y)
    shift.configuration <- as.integer(shift.configuration)
    states <- convergent_regime_states(tree, shift.configuration, regimes)
    regimes <- states$regimes
    fixed.alpha <- NULL
    if(isTRUE(opt$fixed.alpha) && !is.null(base.model$alpha)){
        fixed.alpha <- as.numeric(base.model$alpha)
    }
    fit <- fit_multivariate_ou_likelihood(
        tree, Y, shift.configuration, opt, regimes=regimes,
        fixed.alpha=fixed.alpha
    )
    if(is.null(score)){
        score <- multivariate_full_information_score(
            fit, length(shift.configuration), opt$criterion
        )
    }
    n <- nrow(Y)
    p <- ncol(Y)
    trait.names <- colnames(Y)
    if(is.null(trait.names)){
        trait.names <- paste0("trait", seq_len(p))
        colnames(Y) <- trait.names
    }
    regime.optima <- fit$coefficients
    rownames(regime.optima) <- as.character(seq_len(nrow(regime.optima)) - 1L)
    colnames(regime.optima) <- trait.names
    optima <- regime.optima[states$tip.regime + 1L, , drop=FALSE]
    rownames(optima) <- tree$tip.label
    mu <- fit$fitted.values
    residuals <- fit$residuals
    dimnames(mu) <- dimnames(Y)
    dimnames(residuals) <- dimnames(Y)

    shift.values <- shift.means <- numeric(0)
    if(length(shift.configuration) > 0L){
        target <- regime.optima[states$shift.regime + 1L, , drop=FALSE]
        source <- regime.optima[states$shift.parent.regime + 1L, , drop=FALSE]
        shift.values <- target - source
        alpha.by.trait <- rep(as.numeric(fit$alpha), length.out=p)
        scales <- vapply(alpha.by.trait, function(a){
            edge_scaling_from_cache(
                opt$edge.age, type="orgX", alpha=a
            )[shift.configuration]
        }, numeric(length(shift.configuration)))
        shift.means <- shift.values * scales
        rownames(shift.values) <- rownames(shift.means) <-
            as.character(shift.configuration)
        colnames(shift.values) <- colnames(shift.means) <- trait.names
    }

    model <- if(is.null(base.model)) list() else base.model
    named.shifts <- shift.configuration
    names(named.shifts) <- as.character(states$shift.regime)
    alpha <- rep(as.numeric(fit$alpha), length.out=p)
    names(alpha) <- trait.names
    sigma2 <- diag(fit$trait.covariance)
    names(sigma2) <- trait.names
    sigma2.error <- fit$sigma2.error
    names(sigma2.error) <- trait.names
    model$Y <- Y
    model$tree <- tree
    model$tree.scale <- if(is.null(opt$tree.scale)) 1 else opt$tree.scale
    model$shift.configuration <- named.shifts
    model$shift.values <- shift.values
    model$shift.means <- shift.means
    model$nShifts <- length(shift.configuration)
    model$optima <- optima
    model$alpha <- alpha
    model$sigma2 <- sigma2
    model$sigma2_error <- sigma2.error
    model$trait.covariance <- fit$trait.covariance
    model$trait.correlation <- fit$trait.correlation
    model$tip.trait.covariance <- fit$tip.trait.covariance
    model$tip.trait.correlation <- fit$tip.trait.correlation
    model$covariance.regularization <- opt$covariance.regularization
    model$regularization.lambda <- fit$regularization.lambda
    model$information.criterion.calibrated <- identical(
        opt$covariance.regularization, "none"
    )
    model$score.note <- if(model$information.criterion.calibrated) NULL else
        paste0(
            "BIC is evaluated at a regularized estimate and should be treated ",
            "as a sensitivity score."
        )
    model$optimization <- fit$optimization
    model$diagnostics <- fit$diagnostics
    model$likelihood.engine <- fit$diagnostics$likelihood.engine
    model$parameter.count <- fit$p
    model$information.parameter.count <- fit$p + length(shift.configuration)
    model$nobs <- fit$sample.size
    model$observed.entries <- fit$n
    model$intercept <- regime.optima[1L, ]
    model$mu <- mu
    model$residuals <- residuals
    model$score <- score
    model$cr.score <- score
    model$logLik <- fit$logLik
    model$joint.logLik <- fit$logLik
    model$convergent.regimes <- regimes
    model$convergent <- TRUE
    model$l1ou.options$criterion <- opt$criterion
    model$l1ou.options$root.model <- opt$root.model
    model$l1ou.options$trait.covariance <- "full"
    class(model) <- "l1ou"
    model
}

full_covariance_sparse_inputs <- function(tree, Y, alpha, est.alpha, opt){
    Y <- as.matrix(Y)
    n <- nrow(Y)
    p <- ncol(Y)
    fixed.alpha <- if(est.alpha){
        if(identical(opt$alpha.structure, "diagonal")) rep(0, p) else 0
    } else{
        if(identical(opt$alpha.structure, "diagonal")){
            rep(as.numeric(alpha), length.out=p)
        } else{
            as.numeric(alpha[[1L]])
        }
    }
    null.fit <- fit_multivariate_ou_likelihood(
        tree, Y, integer(), opt, fixed.alpha=fixed.alpha
    )
    covariance <- multivariate_ou_observed_covariance(
        tree, Y, fixed.alpha, null.fit$trait.covariance, opt,
        null.fit$sigma2.error
    )
    covariance.chol <- chol(covariance)
    observed <- !is.na(as.vector(Y))
    y <- as.vector(Y)[observed]
    edge.design <- if(est.alpha){
        rep(list(cached_weighted_design_matrix(
            opt$Z, opt$edge.age, type="apprX"
        )), p)
    } else{
        alpha.by.trait <- rep(fixed.alpha, length.out=p)
        lapply(alpha.by.trait, function(a){
            cached_weighted_design_matrix(
                opt$Z, opt$edge.age, type="orgX", alpha=a
            )
        })
    }
    full.edge.design <- block_diag_matrix(edge.design)[observed, , drop=FALSE]
    intercept.design <- kronecker(diag(p), matrix(1, nrow=n, ncol=1L))[
        observed, , drop=FALSE
    ]
    whitened.y <- drop(forwardsolve(t(covariance.chol), matrix(y, ncol=1L)))
    whitened.edges <- forwardsolve(t(covariance.chol), full.edge.design)
    whitened.intercepts <- forwardsolve(t(covariance.chol), intercept.design)
    qr.intercepts <- qr(whitened.intercepts, tol=1e-9, LAPACK=FALSE)
    rank <- qr.intercepts$rank
    Q <- qr.Q(qr.intercepts, complete=TRUE)
    complement <- Q[, (rank + 1L):nrow(Q), drop=FALSE]

    list(
        y=drop(crossprod(complement, whitened.y)),
        X=crossprod(complement, whitened.edges),
        trait.covariance=null.fit$trait.covariance,
        alpha=fixed.alpha
    )
}
