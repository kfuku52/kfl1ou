branchwise_ou_covariance <- function(
        tree, alpha, rate.multiplier=rep(1, Nedge(tree)),
        root.model=c("OUfixedRoot", "OUrandomRoot")){
    root.model <- match.arg(root.model)
    alpha <- as.numeric(alpha)[[1L]]
    rate.multiplier <- as.numeric(rate.multiplier)
    if(length(rate.multiplier) != Nedge(tree) ||
       any(!is.finite(rate.multiplier)) || any(rate.multiplier <= 0)){
        stop("rate.multiplier must contain one finite positive value per edge.")
    }
    if(!is.finite(alpha) || alpha < 0){
        stop("alpha must be one finite non-negative number.")
    }
    depths <- tree_node_depths(tree)
    tip.depth <- depths[seq_len(length(tree$tip.label))]
    child.depth <- depths[tree$edge[, 2L]]
    Z <- generate_design_matrix(tree, "simpX")
    propagation <- if(alpha <= .Machine$double.eps){
        matrix(1, nrow=nrow(Z), ncol=ncol(Z))
    } else{
        exp(-alpha * outer(tip.depth, child.depth, "-"))
    }
    propagation[Z == 0] <- 0
    innovation <- if(alpha <= .Machine$double.eps){
        tree$edge.length
    } else{
        -expm1(-2 * alpha * tree$edge.length) / (2 * alpha)
    }
    weighted <- sweep(
        propagation, 2L, sqrt(rate.multiplier * innovation), `*`
    )
    covariance <- tcrossprod(weighted)
    if(identical(root.model, "OUrandomRoot")){
        if(alpha <= .Machine$double.eps){
            stop("OUrandomRoot requires a positive alpha.")
        }
        root.propagation <- exp(-alpha * tip.depth)
        covariance <- covariance + tcrossprod(root.propagation) / (2 * alpha)
    }
    0.5 * (covariance + t(covariance))
}


descendant_subtree_edges <- function(tree, stem.edge){
    selected.nodes <- tree$edge[stem.edge, 2L]
    repeat{
        children <- tree$edge[tree$edge[, 1L] %in% selected.nodes, 2L]
        children <- setdiff(children, selected.nodes)
        if(!length(children)) break
        selected.nodes <- c(selected.nodes, children)
    }
    which(tree$edge[, 2L] %in% selected.nodes)
}


profile_branch_rate_model <- function(tree, y, design, alpha, root.model,
                                      stem.edge=NULL,
                                      log.multiplier.bounds=log(c(0.05, 20))){
    n <- length(y)
    evaluate <- function(log.multiplier){
        multiplier <- rep(1, Nedge(tree))
        if(!is.null(stem.edge)){
            multiplier[descendant_subtree_edges(tree, stem.edge)] <-
                exp(log.multiplier)
        }
        covariance <- branchwise_ou_covariance(
            tree, alpha, multiplier, root.model
        )
        factor <- tryCatch(chol(covariance), error=function(e) NULL)
        if(is.null(factor)) return(NULL)
        whitened.y <- drop(forwardsolve(t(factor), matrix(y, ncol=1L)))
        whitened.X <- forwardsolve(t(factor), design)
        coefficient <- lm.fit(whitened.X, whitened.y)$coefficients
        coefficient[!is.finite(coefficient)] <- 0
        residual <- whitened.y - drop(whitened.X %*% coefficient)
        sigma2 <- max(sum(residual^2) / n, .Machine$double.eps)
        logLik <- -0.5 * (
            n * (log(2 * pi) + log(sigma2)) +
            2 * sum(log(diag(factor))) + sum(residual^2) / sigma2
        )
        list(
            logLik=logLik,
            multiplier=if(is.null(stem.edge)) 1 else exp(log.multiplier),
            sigma2=sigma2,
            coefficients=coefficient
        )
    }
    if(is.null(stem.edge)) return(evaluate(0))
    optimized <- stats::optimize(
        function(value){
            fit <- evaluate(value)
            if(is.null(fit)) .Machine$double.xmax / 1000 else -fit$logLik
        },
        interval=log.multiplier.bounds
    )
    evaluate(optimized$minimum)
}


scan_branch_rate_models <- function(model, Y=NULL, candidate.edges=NULL,
                                    min.clade.size=3L,
                                    multiplier.bounds=c(0.05, 20)){
    tree <- model$tree
    y <- if(is.null(Y)) as.numeric(model$Y[, 1L]) else as.numeric(Y[, 1L])
    Z <- generate_design_matrix(tree, "simpX")
    design <- cbind(
        intercept=1,
        Z[, as.integer(model$shift.configuration), drop=FALSE]
    )
    if(is.null(candidate.edges)){
        sizes <- colSums(Z > 0)
        candidate.edges <- which(
            sizes >= min.clade.size &
            sizes <= length(y) - min.clade.size &
            seq_len(Nedge(tree)) < Nedge(tree)
        )
    }
    candidate.edges <- sort(unique(as.integer(candidate.edges)))
    if(!length(candidate.edges)){
        stop("no rate-shift candidate edge satisfies min.clade.size.")
    }
    alpha <- as.numeric(model$alpha)[[1L]]
    null <- profile_branch_rate_model(
        tree, y, design, alpha, model$l1ou.options$root.model
    )
    alternatives <- lapply(candidate.edges, function(edge){
        profile_branch_rate_model(
            tree, y, design, alpha, model$l1ou.options$root.model,
            stem.edge=edge,
            log.multiplier.bounds=log(multiplier.bounds)
        )
    })
    logLik <- vapply(alternatives, function(fit){
        if(is.null(fit)) -Inf else fit$logLik
    }, numeric(1))
    best <- which.max(logLik)
    list(
        null=null,
        alternative=alternatives[[best]],
        edge=candidate.edges[[best]],
        candidate.edges=candidate.edges,
        logLik=logLik,
        statistic=max(0, 2 * (logLik[[best]] - null$logLik))
    )
}


#' Compare a homogeneous OU rate with one clade-specific diffusion-rate shift
#'
#' Scans eligible clades while retaining the fitted optimum-shift mean model
#' and alpha. The clade-wide diffusion multiplier and overall process variance
#' are profiled by maximum likelihood. A parametric bootstrap repeats the full
#' rate-edge scan, thereby calibrating selection of the best rate edge.
#'
#'@param model fitted univariate \code{"l1ou"} model without observation error.
#'@param candidate.edges optional rate-shift stem edges.
#'@param min.clade.size minimum tips inside and outside a candidate clade.
#'@param multiplier.bounds positive optimization bounds for the clade rate.
#'@param nsim number of homogeneous-null simulations for a scan-calibrated
#' p-value. Zero returns only the observed comparison.
#'@param seed optional simulation seed.
#'@return A rate-shift comparison object.
#'@export
compare_diffusion_rate_shift <- function(
        model, candidate.edges=NULL, min.clade.size=3L,
        multiplier.bounds=c(0.05, 20), nsim=0L, seed=NULL){
    check_l1ou_object(model)
    Y <- as.matrix(model$Y)
    if(ncol(Y) != 1L || anyNA(Y)){
        stop("compare_diffusion_rate_shift currently requires one complete trait.")
    }
    if(isTRUE(model$l1ou.options$measurement_error) ||
       !is.null(model$l1ou.options$input_error)){
        stop("rate-shift comparison currently requires no observation-error component.")
    }
    min.clade.size <- l1ou_integer_argument(
        min.clade.size, "min.clade.size", 2L
    )
    if(length(multiplier.bounds) != 2L ||
       any(!is.finite(multiplier.bounds)) || any(multiplier.bounds <= 0) ||
       multiplier.bounds[[1L]] >= multiplier.bounds[[2L]]){
        stop("multiplier.bounds must be two increasing positive numbers.")
    }
    nsim <- l1ou_integer_argument(nsim, "nsim", 0L)
    observed <- scan_branch_rate_models(
        model, candidate.edges=candidate.edges,
        min.clade.size=min.clade.size,
        multiplier.bounds=multiplier.bounds
    )
    simulated <- numeric()
    failures <- 0L
    if(nsim > 0L){
        simulations <- simulate(model, nsim=nsim, seed=seed)
        simulated <- vapply(simulations, function(Ystar){
            value <- tryCatch(
                scan_branch_rate_models(
                    model, Y=Ystar,
                    candidate.edges=observed$candidate.edges,
                    min.clade.size=min.clade.size,
                    multiplier.bounds=multiplier.bounds
                )$statistic,
                error=function(e) NA_real_
            )
            value
        }, numeric(1))
        failures <- sum(!is.finite(simulated))
        simulated <- simulated[is.finite(simulated)]
    }
    p.value <- if(nsim > 0L && length(simulated)){
        (1 + sum(simulated >= observed$statistic)) / (1 + length(simulated))
    } else NA_real_
    result <- list(
        statistic=observed$statistic,
        p.value=p.value,
        best.edge=observed$edge,
        multiplier=observed$alternative$multiplier,
        null.logLik=observed$null$logLik,
        alternative.logLik=observed$alternative$logLik,
        candidate.edges=observed$candidate.edges,
        bootstrap.statistics=simulated,
        successful=length(simulated),
        failed=failures,
        alpha.fixed=as.numeric(model$alpha)[[1L]],
        note=paste0(
            "Alpha and the optimum-shift configuration are conditional on the ",
            "input model; the bootstrap calibrates selection of the rate edge."
        )
    )
    class(result) <- "l1ou_rate_shift_comparison"
    result
}


#'@export
print.l1ou_rate_shift_comparison <- function(x, ...){
    cat("Clade diffusion-rate shift comparison\n")
    cat("best edge:", x$best.edge,
        "rate multiplier:", signif(x$multiplier, 6), "\n")
    cat("likelihood-ratio statistic:", signif(x$statistic, 6), "\n")
    if(is.finite(x$p.value)){
        cat("scan-calibrated bootstrap p-value:", signif(x$p.value, 6), "\n")
    }
    invisible(x)
}
