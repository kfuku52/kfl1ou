check_l1ou_object <- function(object){
    if(!inherits(object, "l1ou")){
        stop("object must inherit from class \"l1ou\".")
    }
    invisible(object)
}

l1ou_named_coefficients <- function(object){
    traits <- colnames(as.matrix(object$Y))
    if(is.null(traits)) traits <- paste0("trait", seq_len(ncol(as.matrix(object$Y))))
    values <- as.numeric(object$intercept)
    names(values) <- paste0("intercept:", traits)
    if(length(object$shift.configuration) > 0L){
        shifts <- as.matrix(object$shift.values)
        if(nrow(shifts) != length(object$shift.configuration)){
            shifts <- t(shifts)
        }
        shift.names <- unlist(lapply(seq_len(ncol(shifts)), function(j){
            paste0("shift[", object$shift.configuration, "]:", traits[[j]])
        }))
        shift.values <- as.vector(shifts)
        names(shift.values) <- shift.names
        values <- c(values, shift.values)
    }
    values
}

#' Standard methods for fitted kfl1ou models
#'
#' Extract coefficients, fitted values, residuals, likelihood, sample size,
#' evolutionary covariance, or tip predictions from an object returned by
#' code{link{fit_OU}} or code{link{estimate_shift_configuration}}.
#'
#'@param object,x fitted object of class code{"l1ou"}.
#'@param type prediction type: expected tip response or adaptive optimum.
#'@param ... additional arguments (currently ignored).
#'@return The conventional extracted model component. code{vcov} returns the
#' fitted evolutionary innovation covariance among traits.
#'@name l1ou-methods
NULL

#'@rdname l1ou-methods
#'@export
coef.l1ou <- function(object, ...){
    check_l1ou_object(object)
    l1ou_named_coefficients(object)
}

#'@rdname l1ou-methods
#'@export
fitted.l1ou <- function(object, ...){
    check_l1ou_object(object)
    object$mu
}

#'@rdname l1ou-methods
#'@export
residuals.l1ou <- function(object, ...){
    check_l1ou_object(object)
    object$residuals
}

#'@rdname l1ou-methods
#'@export
logLik.l1ou <- function(object, ...){
    check_l1ou_object(object)
    value <- if(!is.null(object$joint.logLik)){
        as.numeric(object$joint.logLik)
    } else{
        sum(as.numeric(object$logLik), na.rm=TRUE)
    }
    structure(
        value,
        class="logLik",
        df=if(is.null(object$parameter.count)) NA_integer_ else
            as.numeric(object$parameter.count),
        nobs=nobs.l1ou(object)
    )
}

#'@rdname l1ou-methods
#'@export
nobs.l1ou <- function(object, ...){
    check_l1ou_object(object)
    if(!is.null(object$nobs)) return(as.integer(object$nobs))
    sum(!is.na(as.matrix(object$Y)))
}

#'@rdname l1ou-methods
#'@export
vcov.l1ou <- function(object, ...){
    check_l1ou_object(object)
    if(!is.null(object$trait.covariance)) return(object$trait.covariance)
    covariance <- diag(as.numeric(object$sigma2), nrow=length(object$sigma2))
    names <- colnames(as.matrix(object$Y))
    dimnames(covariance) <- list(names, names)
    covariance
}

#'@rdname l1ou-methods
#'@export
predict.l1ou <- function(object, type=c("response", "optimum"), ...){
    check_l1ou_object(object)
    type <- match.arg(type)
    if(type == "response") object$mu else object$optima
}

l1ou_simulation_covariance <- function(object){
    Y <- as.matrix(object$Y)
    p <- ncol(Y)
    opt <- object$l1ou.options
    if(!is.null(object$trait.covariance)){
        covariance <- multivariate_ou_dense_covariance(
            object$tree, object$alpha, object$trait.covariance,
            opt$root.model
        )
    } else{
        covariance <- block_diag_matrix(lapply(seq_len(p), function(j){
            multivariate_ou_base_covariance(
                object$tree, object$alpha[[j]], opt$root.model
            ) * object$sigma2[[j]]
        }))
    }
    error <- multivariate_observation_error_vector(
        Y, opt, object$sigma2_error
    )
    diag(covariance) <- diag(covariance) + error
    regularize_positive_definite(covariance, 1e-12)
}

#' Simulate from a fitted kfl1ou model
#'
#'@param object fitted object of class code{"l1ou"}.
#'@param nsim number of simulated data sets.
#'@param seed optional random seed.
#'@param preserve.missing logical; retain the fitted data's missingness pattern.
#'@param ... additional arguments (currently ignored).
#'@return A named list of simulated trait matrices.
#'@export
simulate.l1ou <- function(object, nsim=1L, seed=NULL,
                          preserve.missing=TRUE, ...){
    check_l1ou_object(object)
    nsim <- as.integer(nsim[[1L]])
    if(is.na(nsim) || nsim < 1L) stop("nsim must be a positive integer.")
    if(!is.null(seed)) set.seed(seed)
    covariance <- l1ou_simulation_covariance(object)
    factor <- t(chol(covariance))
    mean <- as.vector(object$mu)
    missing <- is.na(as.vector(object$Y))
    simulations <- lapply(seq_len(nsim), function(i){
        value <- mean + drop(factor %*% stats::rnorm(length(mean)))
        if(isTRUE(preserve.missing)) value[missing] <- NA_real_
        matrix(value, nrow=nrow(object$Y), ncol=ncol(object$Y),
               dimnames=dimnames(object$Y))
    })
    names(simulations) <- paste0("sim_", seq_len(nsim))
    class(simulations) <- c("simulated_l1ou", "list")
    simulations
}

#' Diagnose numerical and residual problems in a fitted model
#'
#'@param model fitted object of class code{"l1ou"}.
#'@return A list containing optimizer, covariance, standardized-residual, and
#' potential tip-outlier diagnostics.
#'@export
diagnose_l1ou <- function(model){
    check_l1ou_object(model)
    Y <- as.matrix(model$Y)
    residual <- as.vector(model$residuals)
    observed <- !is.na(residual)
    covariance <- l1ou_simulation_covariance(model)[observed, observed, drop=FALSE]
    factor <- chol(covariance)
    standardized <- drop(forwardsolve(
        t(factor), matrix(residual[observed], ncol=1L)
    ))
    normality <- if(length(standardized) >= 3L && length(standardized) <= 5000L){
        tryCatch(stats::shapiro.test(standardized), error=function(e) NULL)
    } else NULL
    trait.covariance <- vcov(model)
    eigenvalues <- eigen(trait.covariance, symmetric=TRUE,
                         only.values=TRUE)$values
    residual.matrix <- as.matrix(model$residuals)
    tip.trait.covariance <- if(is.null(model$tip.trait.covariance)){
        trait.covariance
    } else model$tip.trait.covariance
    scale <- sqrt(pmax(diag(tip.trait.covariance), .Machine$double.eps))
    tip.score <- rowSums(sweep(residual.matrix, 2L, scale, "/")^2,
                         na.rm=TRUE)
    names(tip.score) <- rownames(Y)
    list(
        convergence=if(is.null(model$diagnostics)) NA else
            model$diagnostics$converged,
        boundary=if(is.null(model$diagnostics)) NA else
            model$diagnostics$boundary,
        likelihood.engine=if(is.null(model$likelihood.engine)) "marginal" else
            model$likelihood.engine,
        trait.covariance.condition=kappa(trait.covariance),
        trait.covariance.eigenvalues=eigenvalues,
        standardized.residuals=standardized,
        standardized.residual.summary=summary(standardized),
        normality.test=normality,
        residual.trait.correlation=stats::cor(
            residual.matrix, use="pairwise.complete.obs"
        ),
        tip.outlier.score=sort(tip.score, decreasing=TRUE)
    )
}

l1ou_inference_parameters <- function(object){
    traits <- colnames(as.matrix(object$Y))
    if(is.null(traits)) traits <- paste0("trait", seq_len(ncol(as.matrix(object$Y))))
    alpha.values <- as.numeric(object$alpha)
    if(!is.null(object$trait.covariance) &&
       identical(object$l1ou.options$alpha.structure, "shared")){
        alpha.values <- alpha.values[[1L]]
        alpha.names <- "alpha:shared"
    } else{
        alpha.names <- paste0("alpha:", traits)
    }
    result <- c(l1ou_named_coefficients(object),
                setNames(alpha.values, alpha.names))
    covariance <- vcov(object)
    covariance.names <- outer(
        colnames(covariance), rownames(covariance),
        function(a, b) paste0("Omega[", a, ",", b, "]")
    )
    keep <- lower.tri(covariance, diag=TRUE)
    omega <- covariance[keep]
    names(omega) <- covariance.names[keep]
    result <- c(result, omega)
    if(!is.null(object$sigma2_error)){
        error <- as.numeric(object$sigma2_error)
        names(error) <- paste0("sigma2_error:", colnames(as.matrix(object$Y)))
        result <- c(result, error)
    }
    result
}

#' Confidence intervals for fitted kfl1ou parameters
#'
#'@param object fitted object of class code{"l1ou"}.
#'@param parm parameter names or indices; by default all reported parameters.
#'@param level confidence level.
#'@param method code{"parametric"} refits simulated data, while
#' code{"wald"} uses an available optimizer Hessian for alpha parameters.
#'@param nsim number of parametric bootstrap replicates.
#'@param seed optional random seed.
#'@param ... additional arguments (currently ignored).
#'@return A matrix with lower and upper confidence limits. Failed bootstrap
#' refits are omitted and reported in an attribute.
#'@export
confint.l1ou <- function(object, parm, level=0.95,
                         method=c("parametric", "wald"), nsim=100L,
                         seed=NULL, ...){
    check_l1ou_object(object)
    method <- match.arg(method)
    estimates <- l1ou_inference_parameters(object)
    if(missing(parm)) parm <- names(estimates)
    if(is.numeric(parm)) parm <- names(estimates)[parm]
    if(anyNA(match(parm, names(estimates)))) stop("unknown parameter in parm.")
    alpha <- (1 - level) / 2
    if(method == "wald"){
        result <- matrix(NA_real_, length(parm), 2L,
                         dimnames=list(parm, paste0(round(c(alpha, 1-alpha)*100, 1), "%")))
        vc <- object$optimization$parameter.vcov
        traits <- colnames(as.matrix(object$Y))
        if(is.null(traits)) traits <- paste0("trait", seq_len(ncol(as.matrix(object$Y))))
        n.alpha <- if(identical(object$l1ou.options$alpha.structure, "diagonal"))
            length(object$alpha) else 1L
        alpha.names <- if(n.alpha == 1L) "alpha:shared" else paste0("alpha:", traits)
        if(!is.null(vc) && nrow(vc) >= n.alpha){
            se.log <- sqrt(pmax(diag(vc)[seq_len(n.alpha)], 0))
            estimate.alpha <- if(n.alpha == 1L) object$alpha[[1L]] else object$alpha
            intervals <- cbind(
                estimate.alpha * exp(stats::qnorm(alpha) * se.log),
                estimate.alpha * exp(stats::qnorm(1-alpha) * se.log)
            )
            rownames(intervals) <- if(n.alpha == 1L) alpha.names[[1L]] else alpha.names
            overlap <- intersect(rownames(intervals), parm)
            result[overlap, ] <- intervals[overlap, , drop=FALSE]
        }
        attr(result, "note") <- "Wald intervals are available only for alpha parameters with a positive-definite optimizer Hessian."
        return(result)
    }

    nsim <- as.integer(nsim[[1L]])
    if(is.na(nsim) || nsim < 2L) stop("nsim must be at least 2.")
    simulations <- simulate(object, nsim=nsim, seed=seed)
    opt <- make_nested_l1ou_options(object$l1ou.options)
    opt$use.saved.scores <- FALSE
    fits <- lapply(simulations, function(Y){
        tryCatch(
            fit_OU(
                object$tree, Y, as.integer(object$shift.configuration),
                cr.regimes=object$convergent.regimes, l1ou.options=opt
            ),
            error=function(e) NULL
        )
    })
    successful <- !vapply(fits, is.null, logical(1))
    if(sum(successful) < 2L) stop("fewer than two bootstrap refits succeeded.")
    draws <- t(vapply(fits[successful], function(fit){
        values <- l1ou_inference_parameters(fit)
        values[names(estimates)]
    }, numeric(length(estimates))))
    result <- t(apply(draws[, parm, drop=FALSE], 2L, stats::quantile,
                      probs=c(alpha, 1-alpha), na.rm=TRUE, names=FALSE))
    colnames(result) <- paste0(round(c(alpha, 1-alpha)*100, 1), "%")
    attr(result, "successful") <- sum(successful)
    attr(result, "failed") <- sum(!successful)
    result
}
