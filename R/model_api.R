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
#' \code{\link{fit_OU}} or \code{\link{estimate_shift_configuration}}.
#'
#'@param object fitted object of class \code{"l1ou"}.
#'@param type prediction type: expected tip response or adaptive optimum.
#'@param component covariance component: conventional coefficient covariance or
#' the evolutionary innovation covariance retained for backward compatibility.
#'@param ... additional arguments (currently ignored).
#'@return The conventional extracted model component. \code{vcov} returns the
#' coefficient covariance when it is available; use
#' \code{\link{evolutionary_vcov}} for evolutionary innovation covariance.
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
        values <- as.numeric(object$logLik)
        if(length(values) == 0L || all(!is.finite(values))){
            stop("the fitted object contains no finite log-likelihood value.")
        }
        sum(values[is.finite(values)])
    }
    structure(
        value,
        class="logLik",
        df=if(!is.null(object$information.parameter.count))
            as.numeric(object$information.parameter.count) else if(
                !is.null(object$parameter.count))
            as.numeric(object$parameter.count) else NA_integer_,
        nobs=nobs.l1ou(object)
    )
}

#'@rdname l1ou-methods
#'@export
nobs.l1ou <- function(object, ...){
    check_l1ou_object(object)
    if(!is.null(object$nobs)) return(as.integer(object$nobs))
    sum(rowSums(!is.na(as.matrix(object$Y))) > 0L)
}

#' Evolutionary innovation covariance from a fitted model
#'
#'@param object fitted object of class \code{"l1ou"}.
#'@return The evolutionary innovation covariance among traits.
#'@export
evolutionary_vcov <- function(object){
    check_l1ou_object(object)
    if(!is.null(object$trait.covariance)) return(object$trait.covariance)
    covariance <- diag(as.numeric(object$sigma2), nrow=length(object$sigma2))
    names <- colnames(as.matrix(object$Y))
    dimnames(covariance) <- list(names, names)
    covariance
}

#'@rdname l1ou-methods
#'@export
vcov.l1ou <- function(object, component=c("coefficients", "evolutionary"), ...){
    check_l1ou_object(object)
    component <- match.arg(component)
    if(component == "evolutionary") return(evolutionary_vcov(object))
    if(!is.null(object$coefficient.vcov)) return(object$coefficient.vcov)
    coefficients <- coef(object)
    result <- matrix(
        NA_real_, length(coefficients), length(coefficients),
        dimnames=list(names(coefficients), names(coefficients))
    )
    warning(
        "coefficient covariance is unavailable for this fit; use ",
        "evolutionary_vcov() for the evolutionary innovation covariance.",
        call.=FALSE
    )
    result
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
    covariance <- 0.5 * (covariance + t(covariance))
    if(is.null(tryCatch(chol(covariance), error=function(e) NULL))){
        covariance <- regularize_positive_definite(covariance, 1e-12)
    }
    covariance
}

with_l1ou_seed <- function(seed, expr){
    if(is.null(seed)) return(force(expr))
    if(length(seed) != 1L || !is.numeric(seed) || !is.finite(seed) ||
       seed < 0 || seed > .Machine$integer.max){
        stop("seed must be NULL or one number between zero and .Machine$integer.max.")
    }
    old.kind <- RNGkind()
    had.seed <- exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
    if(had.seed){
        old.seed <- get(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
    }
    on.exit({
        do.call(RNGkind, as.list(old.kind))
        if(had.seed){
            assign(".Random.seed", old.seed, envir=.GlobalEnv)
        } else if(exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)){
            rm(".Random.seed", envir=.GlobalEnv)
        }
    }, add=TRUE)
    set.seed(as.integer(seed))
    force(expr)
}

draw_zero_mean_gaussian <- function(covariance){
    covariance <- 0.5 * (covariance + t(covariance))
    if(max(abs(covariance)) <= .Machine$double.eps) return(numeric(nrow(covariance)))
    factor <- tryCatch(chol(covariance), error=function(e) NULL)
    if(is.null(factor)){
        covariance <- regularize_positive_definite(covariance, 1e-12)
        factor <- chol(covariance)
    }
    drop(t(factor) %*% stats::rnorm(nrow(covariance)))
}

simulate_l1ou_tree_residual <- function(object){
    tree <- object$tree
    Y <- as.matrix(object$Y)
    n <- nrow(Y)
    p <- ncol(Y)
    alpha <- rep(as.numeric(object$alpha), length.out=p)
    trait.covariance <- if(is.null(object$trait.covariance)){
        diag(as.numeric(object$sigma2), p)
    } else object$trait.covariance
    ordered <- ape::reorder.phylo(tree, "cladewise")
    node.count <- n + ordered$Nnode
    states <- matrix(0, node.count, p)
    root <- setdiff(unique(ordered$edge[, 1L]), ordered$edge[, 2L])[[1L]]
    if(identical(object$l1ou.options$root.model, "OUrandomRoot")){
        states[root, ] <- draw_zero_mean_gaussian(
            ou_stationary_trait_covariance(alpha, trait.covariance)
        )
    }
    for(edge.index in seq_len(nrow(ordered$edge))){
        parent <- ordered$edge[edge.index, 1L]
        child <- ordered$edge[edge.index, 2L]
        branch.length <- ordered$edge.length[[edge.index]]
        transition <- exp(-alpha * branch.length)
        innovation <- ou_branch_innovation_covariance(
            alpha, trait.covariance, branch.length
        )
        states[child, ] <- transition * states[parent, ] +
            draw_zero_mean_gaussian(innovation)
    }
    result <- states[seq_len(n), , drop=FALSE]
    observation.error <- matrix(
        multivariate_observation_error_vector(
            Y, object$l1ou.options, object$sigma2_error
        ), nrow=n, ncol=p
    )
    positive <- observation.error > 0
    result[positive] <- result[positive] +
        stats::rnorm(sum(positive), sd=sqrt(observation.error[positive]))
    dimnames(result) <- dimnames(Y)
    result
}

#' Simulate from a fitted kfl1ou model
#'
#'@param object fitted object of class \code{"l1ou"}.
#'@param nsim number of simulated data sets.
#'@param seed optional random seed.
#'@param preserve.missing logical; retain the fitted data's missingness pattern.
#'@param engine simulation engine. \code{"tree"} generates OU innovations
#' directly along branches without constructing a dense covariance matrix;
#' \code{"dense"} is retained for validation and small problems.
#'@param ... additional arguments (currently ignored).
#'@return A named list of simulated trait matrices.
#'@export
simulate.l1ou <- function(object, nsim=1L, seed=NULL,
                          preserve.missing=TRUE,
                          engine=c("tree", "dense"), ...){
    check_l1ou_object(object)
    nsim <- l1ou_integer_argument(nsim, "nsim", 1L)
    if(length(preserve.missing) != 1L || !is.logical(preserve.missing) ||
       is.na(preserve.missing)){
        stop("preserve.missing must be TRUE or FALSE.")
    }
    engine <- match.arg(engine)
    simulations <- with_l1ou_seed(seed, {
        mean <- as.matrix(object$mu)
        missing <- is.na(as.matrix(object$Y))
        if(engine == "tree"){
            lapply(seq_len(nsim), function(i){
                value <- mean + simulate_l1ou_tree_residual(object)
                if(isTRUE(preserve.missing)) value[missing] <- NA_real_
                value
            })
        } else{
            covariance <- l1ou_simulation_covariance(object)
            factor <- t(chol(covariance))
            lapply(seq_len(nsim), function(i){
                value <- as.vector(mean) +
                    drop(factor %*% stats::rnorm(length(mean)))
                if(isTRUE(preserve.missing)) value[missing] <- NA_real_
                matrix(value, nrow=nrow(object$Y), ncol=ncol(object$Y),
                       dimnames=dimnames(object$Y))
            })
        }
    })
    names(simulations) <- paste0("sim_", seq_len(nsim))
    class(simulations) <- c("simulated_l1ou", "list")
    simulations
}

l1ou_diagnostic_mean_design <- function(model){
    Y <- as.matrix(model$Y)
    n <- nrow(Y)
    p <- ncol(Y)
    if(isTRUE(model$convergent) && !is.null(model$convergent.regimes)){
        states <- convergent_regime_states(
            model$tree, as.integer(model$shift.configuration),
            model$convergent.regimes
        )
        base <- stats::model.matrix(
            ~ factor(states$tip.regime) - 1L
        )
    } else{
        base <- cbind(
            intercept=1,
            generate_design_matrix(model$tree, "simpX")[
                , as.integer(model$shift.configuration), drop=FALSE
            ]
        )
    }
    block_diag_matrix(rep(list(base), p))
}


l1ou_diagnostic_projection <- function(model, max.dense.dimension=2000L){
    Y <- as.matrix(model$Y)
    residual <- as.vector(model$residuals)
    observed <- is.finite(residual)
    if(sum(observed) > max.dense.dimension){
        return(list(
            standardized=NULL,
            tip.score=NULL,
            covariance=NULL,
            observed=observed,
            note=paste0(
                "Dense residual projection skipped because ", sum(observed),
                " observed entries exceed max.dense.dimension."
            )
        ))
    }
    covariance <- l1ou_simulation_covariance(model)[
        observed, observed, drop=FALSE
    ]
    factor <- chol(covariance)
    whitened <- drop(forwardsolve(
        t(factor), matrix(residual[observed], ncol=1L)
    ))
    design <- l1ou_diagnostic_mean_design(model)[observed, , drop=FALSE]
    whitened.design <- forwardsolve(t(factor), design)
    qr.design <- qr(whitened.design, tol=1e-9, LAPACK=FALSE)
    rank <- qr.design$rank
    Q <- qr.Q(qr.design, complete=TRUE)
    standardized <- if(rank < nrow(Q)){
        drop(crossprod(Q[, (rank + 1L):nrow(Q), drop=FALSE], whitened))
    } else numeric()

    precision <- chol2inv(factor)
    gradient <- drop(precision %*% residual[observed])
    observed.indices <- which(observed)
    n <- nrow(Y)
    tip.score <- vapply(seq_len(n), function(i){
        global <- i + (seq_len(ncol(Y)) - 1L) * n
        local <- match(global, observed.indices, nomatch=0L)
        local <- local[local > 0L]
        if(!length(local)) return(NA_real_)
        precision.block <- precision[local, local, drop=FALSE]
        solved <- tryCatch(
            solve(precision.block, gradient[local]),
            error=function(e) NULL
        )
        if(is.null(solved)) return(NA_real_)
        drop(crossprod(gradient[local], solved)) / length(local)
    }, numeric(1))
    names(tip.score) <- rownames(Y)
    list(
        standardized=standardized,
        tip.score=tip.score,
        covariance=covariance,
        observed=observed,
        factor=factor,
        projection=if(rank < nrow(Q)){
            Q[, (rank + 1L):nrow(Q), drop=FALSE]
        } else matrix(0, nrow(Q), 0L),
        rank=rank,
        note=NULL
    )
}


l1ou_clade_variance_statistic <- function(tree, tip.score, min.clade.size=3L){
    keep <- is.finite(tip.score)
    if(sum(keep) < 2L * min.clade.size) return(NA_real_)
    Z <- generate_design_matrix(tree, "simpX") > 0
    values <- vapply(seq_len(ncol(Z)), function(edge){
        inside <- Z[, edge] & keep
        outside <- !Z[, edge] & keep
        if(sum(inside) < min.clade.size || sum(outside) < min.clade.size){
            return(NA_real_)
        }
        abs(log(
            (mean(tip.score[inside]) + 1e-8) /
            (mean(tip.score[outside]) + 1e-8)
        ))
    }, numeric(1))
    if(all(!is.finite(values))) NA_real_ else max(values, na.rm=TRUE)
}


l1ou_residual_statistics <- function(standardized, tree, tip.score,
                                     min.clade.size=3L){
    standardized <- standardized[is.finite(standardized)]
    if(!length(standardized)){
        return(c(max.absolute=NA, skewness=NA, excess.kurtosis=NA,
                 tail.fraction=NA, clade.variance.ratio=NA))
    }
    centered <- standardized - mean(standardized)
    scale <- sqrt(mean(centered^2))
    skewness <- if(scale > 0) mean((centered / scale)^3) else 0
    kurtosis <- if(scale > 0) mean((centered / scale)^4) - 3 else 0
    c(
        max.absolute=max(abs(standardized)),
        skewness=abs(skewness),
        excess.kurtosis=abs(kurtosis),
        tail.fraction=mean(abs(standardized) > 2),
        clade.variance.ratio=l1ou_clade_variance_statistic(
            tree, tip.score, min.clade.size
        )
    )
}


#' Diagnose numerical and residual problems in a fitted model
#'
#'@param model fitted object of class \code{"l1ou"}.
#'@param max.dense.dimension largest observed response dimension for which a
#' dense residual-whitening covariance is constructed.
#'@param nsim number of fitted-model predictive simulations used for absolute
#' adequacy and clade-rate diagnostics. Zero skips simulation calibration.
#'@param seed optional deterministic simulation seed.
#'@param min.clade.size minimum observed tips inside and outside a clade for the
#' diffusion-rate heterogeneity statistic.
#'@return A list containing optimizer, covariance, standardized-residual, and
#' potential tip-outlier diagnostics.
#'@export
diagnose_l1ou <- function(model, max.dense.dimension=2000L, nsim=0L,
                          seed=NULL, min.clade.size=3L){
    check_l1ou_object(model)
    Y <- as.matrix(model$Y)
    max.dense.dimension <- l1ou_integer_argument(
        max.dense.dimension, "max.dense.dimension", 1L
    )
    nsim <- l1ou_integer_argument(nsim, "nsim", 0L)
    min.clade.size <- l1ou_integer_argument(
        min.clade.size, "min.clade.size", 2L
    )
    projected <- l1ou_diagnostic_projection(model, max.dense.dimension)
    standardized <- projected$standardized
    whitening.note <- projected$note
    normality <- if(length(standardized) >= 3L && length(standardized) <= 5000L){
        tryCatch(stats::shapiro.test(standardized), error=function(e) NULL)
    } else NULL
    trait.covariance <- evolutionary_vcov(model)
    eigenvalues <- eigen(trait.covariance, symmetric=TRUE,
                         only.values=TRUE)$values
    residual.matrix <- as.matrix(model$residuals)
    tip.trait.covariance <- if(is.null(model$tip.trait.covariance)){
        trait.covariance
    } else model$tip.trait.covariance
    tip.score <- projected$tip.score
    observed.statistics <- l1ou_residual_statistics(
        standardized, model$tree, tip.score, min.clade.size
    )
    simulated.statistics <- NULL
    predictive.p.value <- stats::setNames(
        rep(NA_real_, length(observed.statistics)), names(observed.statistics)
    )
    if(nsim > 0L && !is.null(projected$projection)){
        covariance.full <- l1ou_simulation_covariance(model)
        observed.index <- which(projected$observed)
        precision <- chol2inv(projected$factor)
        n <- nrow(Y)
        simulated.statistics <- with_l1ou_seed(seed, t(replicate(nsim, {
            residual.full <- draw_zero_mean_gaussian(covariance.full)
            residual.obs <- residual.full[observed.index]
            whitened <- drop(forwardsolve(
                t(projected$factor), matrix(residual.obs, ncol=1L)
            ))
            simulated.standardized <- if(ncol(projected$projection)){
                drop(crossprod(projected$projection, whitened))
            } else numeric()
            gradient <- drop(precision %*% residual.obs)
            simulated.tip.score <- vapply(seq_len(n), function(i){
                global <- i + (seq_len(ncol(Y)) - 1L) * n
                local <- match(global, observed.index, nomatch=0L)
                local <- local[local > 0L]
                if(!length(local)) return(NA_real_)
                block <- precision[local, local, drop=FALSE]
                solved <- tryCatch(solve(block, gradient[local]),
                                   error=function(e) NULL)
                if(is.null(solved)) return(NA_real_)
                drop(crossprod(gradient[local], solved)) / length(local)
            }, numeric(1))
            l1ou_residual_statistics(
                simulated.standardized, model$tree, simulated.tip.score,
                min.clade.size
            )
        })))
        colnames(simulated.statistics) <- names(observed.statistics)
        predictive.p.value <- vapply(seq_along(observed.statistics), function(j){
            values <- simulated.statistics[, j]
            values <- values[is.finite(values)]
            if(!is.finite(observed.statistics[[j]]) || !length(values)){
                return(NA_real_)
            }
            (1 + sum(values >= observed.statistics[[j]])) / (1 + length(values))
        }, numeric(1))
        names(predictive.p.value) <- names(observed.statistics)
    }
    measurement.error <- if(is.null(model$sigma2_error)){
        rep(0, ncol(Y))
    } else rep(as.numeric(model$sigma2_error), length.out=ncol(Y))
    measurement.error.fraction <- measurement.error /
        pmax(measurement.error + diag(tip.trait.covariance),
             .Machine$double.eps)
    names(measurement.error.fraction) <- colnames(Y)
    result <- list(
        convergence=if(is.null(model$diagnostics)) NA else
            model$diagnostics$converged,
        boundary=if(is.null(model$diagnostics)) NA else
            model$diagnostics$boundary,
        likelihood.engine=if(is.null(model$likelihood.engine)) "marginal" else
            model$likelihood.engine,
        trait.covariance.condition=kappa(trait.covariance),
        trait.covariance.eigenvalues=eigenvalues,
        measurement.error.fraction=measurement.error.fraction,
        measurement.error.identifiability.warning=
            isTRUE(model$l1ou.options$measurement_error) &&
            any(measurement.error.fraction < 1e-4 |
                measurement.error.fraction > 0.95),
        standardized.residuals=standardized,
        whitening.note=whitening.note,
        standardized.residual.summary=summary(standardized),
        normality.test=normality,
        residual.statistics=observed.statistics,
        predictive.statistics=simulated.statistics,
        predictive.p.value=predictive.p.value,
        predictive.simulations=nsim,
        rate.heterogeneity.p.value=
            predictive.p.value[["clade.variance.ratio"]],
        rate.heterogeneity.warning=is.finite(
            predictive.p.value[["clade.variance.ratio"]]
        ) && predictive.p.value[["clade.variance.ratio"]] < 0.05,
        raw.residual.trait.correlation=stats::cor(
            residual.matrix, use="pairwise.complete.obs"
        ),
        residual.trait.correlation=stats::cor(
            residual.matrix, use="pairwise.complete.obs"
        ),
        tip.outlier.score=sort(tip.score, decreasing=TRUE)
    )
    class(result) <- "l1ou_diagnostics"
    result
}

#'@export
print.l1ou_diagnostics <- function(x, ...){
    cat("kfl1ou model diagnostics\n")
    cat("converged:", x$convergence, "\n")
    cat("boundary:", x$boundary, "\n")
    cat("likelihood engine:", x$likelihood.engine, "\n")
    cat("trait covariance condition:",
        signif(x$trait.covariance.condition, 5), "\n")
    if(isTRUE(x$measurement.error.identifiability.warning)){
        cat("measurement-error warning: variance fraction is near a boundary\n")
    }
    if(x$predictive.simulations > 0L){
        cat("predictive adequacy p-values:\n")
        print(signif(x$predictive.p.value, 4))
    }
    if(isTRUE(x$rate.heterogeneity.warning)){
        cat("rate-heterogeneity warning: clade variance is poorly calibrated\n")
    }
    if(!is.null(x$whitening.note)) cat(x$whitening.note, "\n")
    invisible(x)
}


#' Parametric diagnostic for clade-specific diffusion-rate heterogeneity
#'
#' Calibrates the largest inside-versus-outside clade contrast in conditional
#' tip residual energy under simulations from the fitted homogeneous-rate OU
#' model. A small p-value is evidence that optimum-shift inference should be
#' checked against a heterogeneous diffusion-rate model.
#'
#'@param model fitted \code{"l1ou"} object.
#'@param nsim number of fitted-model simulations.
#'@param seed optional simulation seed.
#'@param min.clade.size minimum tips on both sides of a tested split.
#'@param max.dense.dimension maximum dense observed dimension.
#'@return A list containing the observed statistic, simulated statistics,
#' p-value, and warning flag.
#'@export
check_rate_heterogeneity <- function(
        model, nsim=199L, seed=NULL, min.clade.size=3L,
        max.dense.dimension=2000L){
    diagnostics <- diagnose_l1ou(
        model,
        max.dense.dimension=max.dense.dimension,
        nsim=nsim,
        seed=seed,
        min.clade.size=min.clade.size
    )
    result <- list(
        statistic=diagnostics$residual.statistics[["clade.variance.ratio"]],
        simulated=if(is.null(diagnostics$predictive.statistics)) numeric() else
            diagnostics$predictive.statistics[, "clade.variance.ratio"],
        p.value=diagnostics$rate.heterogeneity.p.value,
        warning=diagnostics$rate.heterogeneity.warning,
        nsim=nsim,
        min.clade.size=min.clade.size
    )
    class(result) <- "l1ou_rate_heterogeneity"
    result
}


#'@export
print.l1ou_rate_heterogeneity <- function(x, ...){
    cat("Clade diffusion-rate heterogeneity diagnostic\n")
    cat("statistic:", signif(x$statistic, 6), "\n")
    cat("parametric p-value:", signif(x$p.value, 6),
        "(", x$nsim, " simulations)\n")
    invisible(x)
}

l1ou_alpha_parameter_names <- function(object){
    traits <- colnames(as.matrix(object$Y))
    if(is.null(traits)) traits <- paste0("trait", seq_len(ncol(as.matrix(object$Y))))
    if(!is.null(object$trait.covariance) &&
       identical(object$l1ou.options$alpha.structure, "shared")){
        "alpha:shared"
    } else paste0("alpha:", traits)
}

l1ou_inference_parameters <- function(object){
    alpha.values <- as.numeric(object$alpha)
    if(!is.null(object$trait.covariance) &&
       identical(object$l1ou.options$alpha.structure, "shared")){
        alpha.values <- alpha.values[[1L]]
    }
    alpha.names <- l1ou_alpha_parameter_names(object)
    result <- c(l1ou_named_coefficients(object),
                setNames(alpha.values, alpha.names))
    covariance <- evolutionary_vcov(object)
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
#'@param object fitted object of class \code{"l1ou"}.
#'@param parm parameter names or indices; by default all reported parameters.
#'@param level confidence level.
#'@param method \code{"parametric"} refits simulated data, while
#' \code{"wald"} uses an available optimizer Hessian for estimated alpha parameters.
#'@param nsim number of parametric bootstrap replicates.
#'@param seed optional random seed.
#'@param selection whether bootstrap refits condition on the selected shift
#' configuration or repeat the complete shift search. For a full search, an
#' original shift coefficient is recorded as zero when its exact edge is not
#' selected; exact-edge reselection probabilities are returned as attributes.
#'@param interval percentile or basic bootstrap interval.
#'@param ... additional arguments (currently ignored).
#'@return A matrix with lower and upper confidence limits. Failed bootstrap
#' refits are omitted and reported in an attribute.
#'@export
confint.l1ou <- function(object, parm, level=0.95,
                         method=c("parametric", "wald"), nsim=100L,
                         seed=NULL,
                         selection=c("conditional", "full"),
                         interval=c("percentile", "basic"), ...){
    check_l1ou_object(object)
    method <- match.arg(method)
    selection <- match.arg(selection)
    interval <- match.arg(interval)
    if(length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1){
        stop("level must be one number strictly between zero and one.")
    }
    estimates <- l1ou_inference_parameters(object)
    if(missing(parm)) parm <- names(estimates)
    if(is.numeric(parm)) parm <- names(estimates)[parm]
    if(anyNA(match(parm, names(estimates)))) stop("unknown parameter in parm.")
    alpha <- (1 - level) / 2
    if(method == "wald"){
        result <- matrix(NA_real_, length(parm), 2L,
                         dimnames=list(parm, paste0(round(c(alpha, 1-alpha)*100, 1), "%")))
        vc <- object$optimization$parameter.vcov
        alpha.names <- l1ou_alpha_parameter_names(object)
        n.alpha <- length(alpha.names)
        alpha.index <- object$optimization$alpha.parameter.index
        if(!is.null(vc) && length(alpha.index) == n.alpha &&
           all(alpha.index >= 1L & alpha.index <= nrow(vc))){
            se.log <- sqrt(pmax(diag(vc)[alpha.index], 0))
            estimate.alpha <- if(n.alpha == 1L) object$alpha[[1L]] else object$alpha
            intervals <- cbind(
                estimate.alpha * exp(stats::qnorm(alpha) * se.log),
                estimate.alpha * exp(stats::qnorm(1-alpha) * se.log)
            )
            rownames(intervals) <- if(n.alpha == 1L) alpha.names[[1L]] else alpha.names
            overlap <- intersect(rownames(intervals), parm)
            result[overlap, ] <- intervals[overlap, , drop=FALSE]
        }
        attr(result, "note") <- paste0(
            "Wald intervals are available only for estimated alpha parameters ",
            "with an unpenalized positive-definite optimizer Hessian."
        )
        return(result)
    }

    nsim <- l1ou_integer_argument(nsim, "nsim", 2L)
    simulations <- simulate(object, nsim=nsim, seed=seed)
    opt <- make_nested_l1ou_options(object$l1ou.options)
    opt$use.saved.scores <- FALSE
    failures <- character(length(simulations))
    fits <- lapply(seq_along(simulations), function(i){
        Y <- simulations[[i]]
        tryCatch({
            if(selection == "full"){
                estimate_shift_configuration(
                    object$tree, Y, l1ou.options=opt
                )
            } else{
                fit_OU(
                    object$tree, Y, as.integer(object$shift.configuration),
                    cr.regimes=object$convergent.regimes, l1ou.options=opt
                )
            }
        }, error=function(e){
            failures[[i]] <<- conditionMessage(e)
            NULL
        })
    })
    successful <- !vapply(fits, is.null, logical(1))
    if(sum(successful) < 2L) stop("fewer than two bootstrap refits succeeded.")
    parameter.names <- names(estimates)
    shift.parameter <- grepl("^shift\\[[0-9]+\\]:", parameter.names)
    extracted <- lapply(fits[successful], function(fit){
        values <- l1ou_inference_parameters(fit)
        aligned <- rep(NA_real_, length(estimates))
        names(aligned) <- parameter.names
        present <- parameter.names %in% names(values)
        aligned[present] <- values[parameter.names[present]]
        if(selection == "full"){
            aligned[shift.parameter & !present] <- 0
        }
        list(values=aligned, present=present)
    })
    draws <- do.call(rbind, lapply(extracted, `[[`, "values"))
    colnames(draws) <- parameter.names
    present <- do.call(rbind, lapply(extracted, `[[`, "present"))
    colnames(present) <- parameter.names
    target.partition <- shift_partition_key(
        object$tree, object$shift.configuration
    )
    partition.reselected <- vapply(fits[successful], function(fit){
        identical(
            shift_partition_key(fit$tree, fit$shift.configuration),
            target.partition
        )
    }, logical(1))
    selected.draws <- draws[, parm, drop=FALSE]
    effective <- colSums(is.finite(selected.draws))
    result <- t(vapply(seq_along(parm), function(j){
        values <- selected.draws[, j]
        values <- values[is.finite(values)]
        if(length(values) < 2L) return(c(NA_real_, NA_real_))
        limits <- stats::quantile(
            values, probs=c(alpha, 1-alpha), names=FALSE
        )
        if(interval == "basic"){
            limits <- 2 * estimates[[parm[[j]]]] - rev(limits)
        }
        limits
    }, numeric(2)))
    colnames(result) <- paste0(round(c(alpha, 1-alpha)*100, 1), "%")
    attr(result, "successful") <- sum(successful)
    attr(result, "failed") <- sum(!successful)
    attr(result, "effective") <- effective
    attr(result, "failure.messages") <- sort(
        table(failures[nzchar(failures)]), decreasing=TRUE
    )
    attr(result, "selection") <- selection
    attr(result, "interval") <- interval
    attr(result, "reselection.count") <- colSums(present[, parm, drop=FALSE])
    attr(result, "reselection.probability") <-
        colMeans(present[, parm, drop=FALSE])
    attr(result, "nonselection.value") <- if(selection == "full") 0 else NA_real_
    attr(result, "partition.reselection.count") <- sum(partition.reselected)
    attr(result, "partition.reselection.probability") <-
        mean(partition.reselected)
    result
}

#' Conditional likelihood profile for OU adaptation rates
#'
#'@param model fitted object of class \code{"l1ou"}.
#'@param alpha.grid numeric candidate values for a shared/univariate alpha, or
#' a matrix with one column per trait for a full-covariance diagonal-alpha fit.
#' Zero may be included to inspect the Brownian-motion boundary under a fixed
#' root model.
#'@param keep.fits retain successful fitted models as an attribute.
#'@return A data frame containing candidate alpha values, log likelihoods,
#' relative likelihood deviances, information scores, and failure messages.
#'@export
profile_alpha_l1ou <- function(model, alpha.grid, keep.fits=FALSE){
    check_l1ou_object(model)
    if(length(keep.fits) != 1L || !is.logical(keep.fits) || is.na(keep.fits)){
        stop("keep.fits must be TRUE or FALSE.")
    }
    Y <- as.matrix(model$Y)
    p <- ncol(Y)
    full <- !is.null(model$trait.covariance)
    diagonal.alpha <- full &&
        identical(model$l1ou.options$alpha.structure, "diagonal")
    if(diagonal.alpha){
        grid <- as.matrix(alpha.grid)
        if(ncol(grid) != p){
            stop("alpha.grid must have one column per trait.")
        }
        colnames(grid) <- colnames(Y)
    } else{
        grid <- matrix(as.numeric(alpha.grid), ncol=1L,
                       dimnames=list(NULL, "shared"))
    }
    if(nrow(grid) < 2L || any(!is.finite(grid)) || any(grid < 0)){
        stop("alpha.grid must contain at least two rows of finite non-negative values.")
    }
    fits <- vector("list", nrow(grid))
    failures <- character(nrow(grid))
    for(i in seq_len(nrow(grid))){
        opt <- make_nested_l1ou_options(model$l1ou.options)
        opt$alpha.lower.bound <- as.numeric(grid[i, ])
        opt$alpha.upper.bound <- as.numeric(grid[i, ])
        opt$alpha.starting.value <- as.numeric(grid[i, ])
        fits[[i]] <- tryCatch(
            fit_OU(
                model$tree, model$Y,
                as.integer(model$shift.configuration),
                cr.regimes=model$convergent.regimes,
                l1ou.options=opt
            ),
            error=function(e){
                failures[[i]] <<- conditionMessage(e)
                NULL
            }
        )
    }
    log.likelihood <- vapply(fits, function(fit){
        if(is.null(fit)) NA_real_ else as.numeric(logLik(fit))
    }, numeric(1))
    score <- vapply(fits, function(fit){
        if(is.null(fit)) NA_real_ else as.numeric(fit$score)
    }, numeric(1))
    best <- if(any(is.finite(log.likelihood))) max(log.likelihood, na.rm=TRUE) else
        NA_real_
    result <- data.frame(
        grid,
        logLik=log.likelihood,
        delta.deviance=if(is.finite(best)) 2 * (best - log.likelihood) else NA_real_,
        score=score,
        error=failures,
        check.names=FALSE
    )
    if(isTRUE(keep.fits)) attr(result, "fits") <- fits
    class(result) <- c("l1ou_alpha_profile", class(result))
    result
}
